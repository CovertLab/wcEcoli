"""
this is the script that does the bulk of data extraction.

reads each cells simOut just once and retain the data used for downstream analysis by other scripts. This should alsways be done first 

From this data we should be able to derive def 5_CI classifiaction with a 95% CI, along with the base subgen definitions 4 and 5 

Outputs are written to plotOutDir under the FIXED basename `subgen_raw_extract`
(sc.raw_extract_prefix), NOT under plotOutFileName, so consumers can always find
them regardless of what invoked this script:
  ..._synth_per_cell.tsv      wide: seed, generation, is_successful_lineage +
                              one column per gene = completed transcripts that
                              generation. PRIMARY substrate.
  ..._max_mrna_per_cell.tsv   wide: same index + per-gene max mRNA count.
  ..._max_protein_per_cell.tsv wide: same index + per-gene max protein count.
  ..._frac_protein_zero_per_cell.tsv wide: same index + per-gene time-weighted
                              fraction of the generation's cell cycle spent at
                              zero protein copies (the protein-absence rate).
  ..._lineage_success.tsv     per seed: strict successful-lineage flags.
  ..._doubling_times.tsv      per seed: doubling time (min) for each generation,
                              -1 where the generation did not run or was
                              unreadable. Already computed for the strict filter,
                              so persisting it saves every consumer a walk over
                              every cell's Main/time.
  ..._genes.tsv               gene_id, cistron_id, monomer_id (column key).
  ..._run_metadata.json       provenance.

Run this ONCE before the downstream subgenerational plots:
  python runscripts/manual/analysisCohort.py --plot subgen_raw_extract.py <dir>
"""

import csv
import os
import pickle
from datetime import datetime

import numpy as np

from models.ecoli.analysis import cohortAnalysisPlot
from models.ecoli.analysis.cohort import subgen_common as sc
from wholecell.io.tablereader import TableReader


def _time_weighted_frac_zero(counts, time):
	"""Fraction of the cell cycle each column spends at zero copies.

	`counts` is a [n_timesteps, n_genes] per-timestep count array; `time` is the
	[n_timesteps] cumulative simulation time (s). Weight each timestep by its
	duration so the result is a true fraction of the cell's lifetime (the
	snapshot probability of observing zero copies), returning a [n_genes] vector
	in [0, 1]. Fall back to an unweighted timestep fraction if the clock is
	degenerate (near-identical at ~1 s steps).
	"""
	zero = counts == 0
	if time is not None and len(time) == len(counts) and len(time) > 1:
		dt = np.diff(np.asarray(time, dtype=float))
		total = dt.sum()
		if total > 0:
			return (dt[:, None] * zero[1:]).sum(axis=0) / total
	return zero.mean(axis=0)


class Plot(cohortAnalysisPlot.CohortAnalysisPlot):
	def do_plot(self, variantDir, plotOutDir, plotOutFileName, simDataFile,
			validationDataFile, metadata):
		analysis_run_time = datetime.now().isoformat(timespec='seconds')
		with open(simDataFile, 'rb') as f:
			sim_data = pickle.load(f)

		# check for generation number being greater than ignore_first
	
		n_generation = self.ap.n_generation
		if n_generation <= sc.IGNORE_FIRST_N_GENS:
			print('Skipping extraction - not enough generations run.')
			return

		mRNA_cistron_ids, monomer_ids, gene_ids = sc.get_mrna_gene_set(sim_data)
		n_genes = len(mRNA_cistron_ids)

		sim_metadata_path, sim_metadata = sc.load_sim_metadata(variantDir)
		total_init_sims = sim_metadata.get('total_init_sims')
		success = sc.compute_lineage_success(
			self.ap, n_generation, total_init_sims=total_init_sims)

		# Write under the fixed raw-extract basename (not plotOutFileName) so the
		# downstream script can always find these files.
		prefix = sc.raw_extract_prefix(plotOutDir)

		# Gene metadata / column key (always writable).
		self._write_genes(prefix + sc.GENES_SUFFIX,
			gene_ids, mRNA_cistron_ids, monomer_ids)
		# Lineage-success table (always writable, independent of the def-5 column).
		self._write_lineage_success(prefix + sc.LINEAGE_SUCCESS_SUFFIX,
			success, n_generation)

		cell_paths = self.ap.get_cells(
			generation=np.arange(sc.IGNORE_FIRST_N_GENS, n_generation),
			only_successful=True)
		print('Extracting from %d cells (gen >= %d)...'
			% (len(cell_paths), sc.IGNORE_FIRST_N_GENS))
		if len(cell_paths) == 0:
			print('No successful cells found. Skipping matrices.')
			return

		# Build index maps from our gene list into each listener's subcolumns.
		first_sim_out = os.path.join(cell_paths[0], 'simOut')

		rna_reader = TableReader(os.path.join(first_sim_out, 'RNACounts'))
		mRNA_cistron_ids_table = rna_reader.readAttribute('mRNA_cistron_ids')
		rna_reader.close()
		rna_id_to_index = {
			cid: i for i, cid in enumerate(mRNA_cistron_ids_table)}
		mrna_indexes = np.array([
			rna_id_to_index[c] for c in mRNA_cistron_ids])

		rnap_reader = TableReader(os.path.join(first_sim_out, 'RnapData'))
		full_cistron_ids = rnap_reader.readAttribute('cistron_ids')
		rnap_reader.close()
		full_id_to_index = {c: i for i, c in enumerate(full_cistron_ids)}
		cistron_indexes = np.array([
			full_id_to_index[c] for c in mRNA_cistron_ids])

		monomer_reader = TableReader(os.path.join(first_sim_out, 'MonomerCounts'))
		monomer_ids_table = monomer_reader.readAttribute('monomerIds')
		monomer_reader.close()
		monomer_id_to_index = {m: i for i, m in enumerate(monomer_ids_table)}
		monomer_indexes = np.array([
			monomer_id_to_index[m] for m in monomer_ids])

		# Definitions 4/5 need countRnaCistronSynthesized, present only in
		# simulations run after that listener was added. Detect up front.
		has_synth = True
		try:
			TableReader(os.path.join(first_sim_out, sc.SYNTH_TABLE)
				).readColumn(sc.SYNTH_COLUMN)
		except Exception:
			has_synth = False
			print('WARNING: %s/%s not found in this cohort. The synth matrix '
				'(and thus Definition 5 / Task 4) cannot be produced; only the '
				'max-count matrices, gene key, and lineage-success table are '
				'written. Re-run simulations after the listener change.'
				% (sc.SYNTH_TABLE, sc.SYNTH_COLUMN))

		# Single pass: reduce each cell to per-gene vectors.
		meta_rows = []
		synth_list = []
		max_mrna_list = []
		max_protein_list = []
		frac_protein_zero_list = []
		skipped = []
		included = []

		for i, cell_path in enumerate(cell_paths):
			if i % 100 == 0:
				print('  Cell %d/%d' % (i, len(cell_paths)))
			sim_out = os.path.join(cell_path, 'simOut')
			seed, gen = sc.parse_cell_id(cell_path)
			try:
				mRNA_counts = TableReader(os.path.join(sim_out, 'RNACounts')
					).readColumn('mRNA_cistron_counts')[:, mrna_indexes]
				monomer_counts = TableReader(os.path.join(sim_out, 'MonomerCounts')
					).readColumn('monomerCounts')[:, monomer_indexes]
				time = TableReader(os.path.join(sim_out, 'Main')
					).readColumn('time')
				if has_synth:
					synth_events = TableReader(
						os.path.join(sim_out, sc.SYNTH_TABLE)
						).readColumn(sc.SYNTH_COLUMN)[:, cistron_indexes]
			except Exception as e:
				print('  Warning: could not read %s: %s' % (cell_path, e))
				skipped.append(cell_path)
				continue

			is_succ = 1 if seed in success['successful_seeds'] else 0
			meta_rows.append([seed, gen, is_succ])
			max_mrna_list.append(mRNA_counts.max(axis=0))
			max_protein_list.append(monomer_counts.max(axis=0))
			frac_protein_zero_list.append(
				_time_weighted_frac_zero(monomer_counts, time))
			if has_synth:
				synth_list.append(synth_events.sum(axis=0))
			included.append(cell_path)

		if not included:
			print('No readable cells found. Skipping matrices.')
			return
		print('Extracted %d cells.' % len(included))

		meta_header = ['seed', 'generation', 'is_successful_lineage']

		# Completed-transcript counts per cell (the def-5 / Task-4 substrate).
		if has_synth:
			sc.write_per_cell_matrix(
				prefix + sc.SYNTH_PER_CELL_SUFFIX, meta_header, meta_rows,
				gene_ids, np.array(synth_list, dtype=np.int64), value_fmt=None)
			print('Wrote %s' % (prefix + sc.SYNTH_PER_CELL_SUFFIX))

		# Descriptive max copy numbers (for the def-5 expression table).
		sc.write_per_cell_matrix(
			prefix + sc.MAX_MRNA_PER_CELL_SUFFIX, meta_header, meta_rows,
			gene_ids, np.array(max_mrna_list, dtype=np.int64), value_fmt=None)
		sc.write_per_cell_matrix(
			prefix + sc.MAX_PROTEIN_PER_CELL_SUFFIX, meta_header, meta_rows,
			gene_ids, np.array(max_protein_list, dtype=np.int64), value_fmt=None)

		# Time-weighted fraction of the cell cycle at zero protein copies (the
		# protein-absence rate for the protein-memory analysis).
		sc.write_per_cell_matrix(
			prefix + sc.FRAC_PROTEIN_ZERO_PER_CELL_SUFFIX, meta_header, meta_rows,
			gene_ids, np.array(frac_protein_zero_list, dtype=np.float64),
			value_fmt='%.6g')
		print('Wrote %s' % (prefix + sc.FRAC_PROTEIN_ZERO_PER_CELL_SUFFIX))

		# compute_lineage_success already walked every cell's Main/time to build the
		# strict filter, so persisting the per-generation doubling times here is free
		# and saves every consumer that walk.
		self._write_doubling_times(
			prefix + sc.DOUBLING_TIMES_SUFFIX, success, n_generation)

		self._write_metadata(
			prefix + '_run_metadata.json', analysis_run_time,
			sim_metadata_path, sim_metadata, has_synth, n_genes,
			success, included, skipped)
		print('Done.')


	def _write_genes(self, path, gene_ids, cistron_ids, monomer_ids):
		with open(path, 'w') as f:
			w = csv.writer(f, delimiter='\t')
			w.writerow(['gene_id', 'cistron_id', 'monomer_id'])
			for g, c, m in zip(gene_ids, cistron_ids, monomer_ids):
				w.writerow([g, c, m])
		print('Wrote %s' % path)

	def _write_lineage_success(self, path, success, n_generation):
		with open(path, 'w') as f:
			w = csv.writer(f, delimiter='\t')
			w.writerow(['seed', 'n_gens_ran', 'reached_final_gen',
				'completed_all_gens', 'n_cells_at_180', 'gens_at_180',
				'max_doubling_min', 'is_successful'])
			for s in success['all_seed_ids']:
				gens = success['successful_gens'][s]
				valid = success['doubling'][s][success['doubling'][s] >= 0]
				max_dt = float(valid.max()) if valid.size else -1
				w.writerow([
					'%06d' % s, len(gens), (n_generation - 1) in gens,
					success['completed_all'][s], success['n_at_180'][s],
					','.join(str(g) for g in success['gens_at_180'][s]),
					'%.4g' % max_dt, success['in_successful'][s]])
		print('Wrote %s' % path)

	def _write_doubling_times(self, path, success, n_generation):
		"""Per-seed doubling time for each generation, in minutes (-1 if absent).

		Same layout as subgen_definition5_lineage_ci.py's _doubling_times table, so
		either file can be read interchangeably.
		"""
		with open(path, 'w') as f:
			w = csv.writer(f, delimiter='\t')
			w.writerow(['seed'] + ['gen_%d' % g for g in range(n_generation)])
			for s in success['all_seed_ids']:
				row = success['doubling'][s]
				w.writerow(['%06d' % s]
					+ ['%.4g' % v if v >= 0 else '-1' for v in row])
		print('Wrote %s' % path)

	def _write_metadata(self, path, analysis_run_time, sim_metadata_path,
			sim_metadata, has_synth, n_genes, success, included, skipped):
		# Layout comes from sc.write_run_metadata so every subgen script's
		# provenance stays comparable; the blocks below are this script's own.
		sc.write_run_metadata(
			path,
			script=os.path.basename(__file__),
			parameters={
				'ignore_first_n_gens': sc.IGNORE_FIRST_N_GENS,
				'max_doubling_min': sc.MAX_DOUBLING_MIN,
				'synth_available': has_synth,
				},
			sim_metadata_path=sim_metadata_path,
			sim_metadata=sim_metadata,
			run_time=analysis_run_time,
			extra={
				'genes': {'n_genes': n_genes},
				'lineages': {
					'n_seeds': len(success['seeds']),
					'n_successful': len(success['successful_seeds']),
					'successful_seeds': sorted(success['successful_seeds']),
					},
				'cells': {
					'n_included': len(included),
					'n_skipped': len(skipped),
					'skipped': skipped,
					},
				})


if __name__ == '__main__':
	Plot().cli()
