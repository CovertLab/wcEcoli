"""
Single-pass raw extraction for the subgenerational-expression pipeline.

Reads each cell's simOut EXACTLY ONCE and persists the granular per-cell data
that every downstream subgenerational analysis needs, so those analyses never
re-read simulation output. From this one extraction you can derive Definition 5
Form A (pooled-cell threshold), Form B (per-lineage confidence interval, the
canonical subgen label), Definition 4 (any completed transcript), and the
per-generation expression matrix (Task 4).

Definition 5 substrate is `countRnaCistronSynthesized` (completed transcripts,
insulated from tRNA attenuation) summed over each generation's timesteps, per
protein-coding mRNA cistron.

Outputs (written to plotOutDir, prefixed by plotOutFileName):
  <name>_synth_per_cell.tsv      wide: seed, generation, is_successful_lineage +
                                 one column per gene = completed transcripts
                                 that generation. PRIMARY substrate.
  <name>_max_mrna_per_cell.tsv   wide: same index + per-gene max mRNA count.
  <name>_max_protein_per_cell.tsv wide: same index + per-gene max protein count.
  <name>_lineage_success.tsv     per seed: strict successful-lineage flags.
  <name>_genes.tsv               gene_id, cistron_id, monomer_id (column key).
  <name>_run_metadata.json       provenance.

Run this ONCE before the downstream subgenerational plots:
  python runscripts/manual/analysisCohort.py --plot subgen_raw_extract.py <dir>
"""

import csv
import json
import os
import pickle
from datetime import datetime

import numpy as np

from models.ecoli.analysis import cohortAnalysisPlot
from models.ecoli.analysis.cohort import subgen_common as sc
from wholecell.io.tablereader import TableReader


class Plot(cohortAnalysisPlot.CohortAnalysisPlot):
	def do_plot(self, variantDir, plotOutDir, plotOutFileName, simDataFile,
			validationDataFile, metadata):
		analysis_run_time = datetime.now().isoformat(timespec='seconds')
		with open(simDataFile, 'rb') as f:
			sim_data = pickle.load(f)

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
		# downstream consumers can always find these files.
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

		self._write_metadata(
			prefix + '_run_metadata.json', analysis_run_time,
			sim_metadata_path, sim_metadata, has_synth, n_genes,
			success, included, skipped)
		print('Done.')

	# ------------------------------------------------------------------ writers

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

	def _write_metadata(self, path, analysis_run_time, sim_metadata_path,
			sim_metadata, has_synth, n_genes, success, included, skipped):
		repo_dir = os.path.dirname(os.path.abspath(__file__))
		meta = {
			'analysis': {
				'script': os.path.basename(__file__),
				'run_time': analysis_run_time,
				'git': sc.git_info(repo_dir),
				'parameters': {
					'ignore_first_n_gens': sc.IGNORE_FIRST_N_GENS,
					'max_doubling_min': sc.MAX_DOUBLING_MIN,
					'synth_available': has_synth,
					},
				},
			'simulation': {
				'metadata_source': sim_metadata_path,
				'git_hash': sim_metadata.get('git_hash'),
				'git_branch': sim_metadata.get('git_branch'),
				'run_time': sim_metadata.get('time'),
				'description': sim_metadata.get('description'),
				'variant': sim_metadata.get('variant'),
				'total_gens': sim_metadata.get('total_gens'),
				'total_init_sims': sim_metadata.get('total_init_sims'),
				},
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
			}
		with open(path, 'w') as f:
			json.dump(meta, f, indent=2)
		print('Wrote %s' % path)


if __name__ == '__main__':
	Plot().cli()
