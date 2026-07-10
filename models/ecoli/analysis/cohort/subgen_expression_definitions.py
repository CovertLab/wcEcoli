"""
Compute five RNA-based definitions of subgenerational gene expression as per-gene
"expression probabilities", then report how many genes each definition classifies
as subgenerational across a sweep of probability thresholds.

Definitions (each computed per cell-generation, then averaged over all cells as
`value = sum over cells / number of cells`):

	1. mRNA present this gen (count >= 1)                  -- matches Science Fig 4C
	2. any transcription (initiation) event this gen
	3. number of transcription (initiation) events this gen
	4. any successful (completed) transcription event this gen
	5. number of successful (completed) transcription events this gen

Definitions 1/2/4 are probabilities in [0, 1]; definitions 3/5 are mean events per
cell. For every definition, a value < 1 indicates subgenerational expression.

Definitions 4/5 count completed transcripts (from countRnaCistronSynthesized) and
so are insulated from tRNA attenuation of the amino-acid-biosynthesis operons,
whereas definitions 2/3 count initiation events that may be attenuated before
completing.

This is a superset of subgenerational_expression_table.py: it uses the same gene
set (mRNA cistrons with an associated protein) and reproduces that file's
p_present, max_mRNA_count, and max_protein_count columns.
"""

import pickle
import os
import csv
import json
import subprocess
from datetime import datetime

import numpy as np

from models.ecoli.analysis import cohortAnalysisPlot
from wholecell.io.tablereader import TableReader
from wholecell.utils import constants


IGNORE_FIRST_N_GENS = 8
# Drop the inherited boundary timestep from each generation when True. Default
# False to match subgenerational_expression_table.py for direct comparison.
REMOVE_FIRST_TIMESTEP = False
# Probability filters for the subgenerational classification sweep.
THRESHOLDS = [1.0, 0.99, 0.95]
# Report the absent-cell list only for near-ubiquitous genes that are absent in at
# most this many cells (i.e. expression probability close to 1).
MAX_ABSENT_CELLS_TO_REPORT = 10

# Output column label for each definition, in order.
DEFINITION_LABELS = [
	'p_present_def1',
	'p_any_init_def2',
	'mean_init_events_def3',
	'p_any_synth_def4',
	'mean_synth_events_def5',
	]


def _blank_if_nan(value):
	"""Render NaN values (defs 4/5 on pre-listener-change cohorts) as blank."""
	return '' if isinstance(value, float) and np.isnan(value) else value


def _git_info(repo_dir):
	"""Return the current git hash, branch, and dirty flag for repo_dir."""
	def run(args):
		return subprocess.check_output(
			['git', '-C', repo_dir] + args,
			stderr=subprocess.DEVNULL).decode().strip()
	try:
		return {
			'git_hash': run(['rev-parse', 'HEAD']),
			'git_branch': run(['rev-parse', '--abbrev-ref', 'HEAD']),
			'git_dirty': bool(run(['status', '--porcelain'])),
			}
	except Exception as e:
		return {'git_hash': None, 'git_branch': None, 'git_dirty': None,
			'error': str(e)}


def _load_sim_metadata(variant_dir):
	"""Load the simulation's metadata.json (git hash, run time, options).

	The sim-level metadata directory sits one level above the variant directory;
	fall back to a metadata directory inside the variant directory.
	"""
	candidates = [
		os.path.join(os.path.dirname(variant_dir),
			constants.METADATA_DIR, constants.JSON_METADATA_FILE),
		os.path.join(variant_dir,
			constants.METADATA_DIR, constants.JSON_METADATA_FILE),
		]
	for path in candidates:
		if os.path.isfile(path):
			with open(path) as f:
				return path, json.load(f)
	return None, {}


def _parse_cell_id(cell_path):
	"""Split a cell path into (seed, generation).

	Cell paths look like .../<variant>/<seed>/generation_<gen>/<daughter>.
	"""
	parts = cell_path.rstrip(os.sep).split(os.sep)
	for i, part in enumerate(parts):
		if part.startswith('generation_'):
			seed = parts[i - 1] if i > 0 else ''
			return seed, int(part.split('_')[1])
	return '', -1


class Plot(cohortAnalysisPlot.CohortAnalysisPlot):
	# When True, restrict the analysis to seeds whose lineage reached the final
	# generation (subgen_expression_definitions_complete.py sets this) so that
	# cells from lineages that died mid-run do not contribute absences.
	REQUIRE_COMPLETE_LINEAGE = False

	def do_plot(self, variantDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		analysis_run_time = datetime.now().isoformat(timespec='seconds')

		with open(simDataFile, 'rb') as f:
			sim_data = pickle.load(f)

		# Ignore data from a predefined number of generations per seed
		if self.ap.n_generation <= IGNORE_FIRST_N_GENS:
			print('Skipping analysis - not enough generations run.')
			return

		# Canonical gene set: mRNA cistrons with an associated protein/monomer,
		# in cistron order (matches subgenerational_expression_table.py).
		cistron_data = sim_data.process.transcription.cistron_data
		cistron_ids = cistron_data['id']
		cistron_id_to_protein_id = {
			protein['cistron_id']: protein['id']
			for protein in sim_data.process.translation.monomer_data
			}
		mRNA_cistron_ids = [
			cistron_id for cistron_id in cistron_ids
			if cistron_id in cistron_id_to_protein_id]
		monomer_ids = [
			cistron_id_to_protein_id[cistron_id]
			for cistron_id in mRNA_cistron_ids]
		cistron_id_to_gene_id = {
			cistron['id']: cistron['gene_id'] for cistron in cistron_data
			}
		gene_ids = [
			cistron_id_to_gene_id[cistron_id]
			for cistron_id in mRNA_cistron_ids]
		n_genes = len(mRNA_cistron_ids)

		cell_paths = self.ap.get_cells(
			generation=np.arange(IGNORE_FIRST_N_GENS, self.ap.n_generation),
			only_successful=True)

		# Optionally keep only lineages that completed every generation: a seed
		# qualifies when its last successful generation is the final one. Cells from
		# lineages that died mid-run (a major source of near-ubiquitous absences)
		# are then excluded.
		if self.REQUIRE_COMPLETE_LINEAGE:
			final_gen = self.ap.n_generation - 1
			last_gen_by_seed = {}
			for path in self.ap.get_cells(only_successful=True):
				seed, gen = _parse_cell_id(path)
				last_gen_by_seed[seed] = max(last_gen_by_seed.get(seed, -1), gen)
			complete_seeds = {
				seed for seed, last in last_gen_by_seed.items()
				if last == final_gen}
			cell_paths = np.array([
				path for path in cell_paths
				if _parse_cell_id(path)[0] in complete_seeds])
			print('Restricting to %d complete lineages (reached generation %d).'
				% (len(complete_seeds), final_gen))

		print('Analyzing %d cells...' % len(cell_paths))

		if len(cell_paths) == 0:
			print('No successful cells found. Skipping analysis.')
			return

		# Build index maps from our gene list into each listener's subcolumns.
		first_sim_out = os.path.join(cell_paths[0], 'simOut')

		# Definition 1: RNACounts/mRNA_cistron_counts (subcolumn mRNA_cistron_ids).
		# Read the id list from the listener attribute, not the process.
		rna_counts_reader = TableReader(os.path.join(first_sim_out, 'RNACounts'))
		mRNA_cistron_ids_table = rna_counts_reader.readAttribute('mRNA_cistron_ids')
		rna_counts_reader.close()
		mRNA_counts_id_to_index = {
			cistron_id: i for (i, cistron_id)
			in enumerate(mRNA_cistron_ids_table)
			}
		def1_indexes = np.array([
			mRNA_counts_id_to_index[cistron_id]
			for cistron_id in mRNA_cistron_ids])

		# Definitions 2-5 are indexed by the full cistron_ids, exposed by both
		# RnapData and TranscriptElongationListener.
		rnap_reader = TableReader(os.path.join(first_sim_out, 'RnapData'))
		full_cistron_ids = rnap_reader.readAttribute('cistron_ids')
		rnap_reader.close()
		full_cistron_id_to_index = {
			cistron_id: i for (i, cistron_id)
			in enumerate(full_cistron_ids)
			}
		full_cistron_indexes = np.array([
			full_cistron_id_to_index[cistron_id]
			for cistron_id in mRNA_cistron_ids])

		# Max protein counts: MonomerCounts/monomerCounts (subcolumn monomerIds).
		monomer_reader = TableReader(os.path.join(first_sim_out, 'MonomerCounts'))
		monomer_ids_table = monomer_reader.readAttribute('monomerIds')
		monomer_reader.close()
		monomer_id_to_index = {
			monomer_id: i for (i, monomer_id)
			in enumerate(monomer_ids_table)
			}
		monomer_indexes = np.array([
			monomer_id_to_index[monomer_id] for monomer_id in monomer_ids])

		# Definitions 4/5 require the countRnaCistronSynthesized column, which only
		# exists in simulations run after that listener was added. Detect it up
		# front; if absent, still report definitions 1-3 and mark 4/5 as NA.
		has_synth = True
		try:
			TableReader(os.path.join(first_sim_out, 'TranscriptElongationListener')
				).readColumn('countRnaCistronSynthesized')
		except Exception:
			has_synth = False
			print('WARNING: countRnaCistronSynthesized not found in this cohort; '
				'definitions 4/5 will be reported as NA. Re-run simulations after '
				'the listener change to populate them.')

		# Single memory-efficient pass: read the per-cell tables, reduce each to a
		# per-generation statistic, and accumulate O(n_genes) arrays. A cell is
		# only counted if all required tables read successfully, so every
		# definition shares an identical cell set.
		time_slice = slice(1, None) if REMOVE_FIRST_TIMESTEP else slice(None)

		sum_present = np.zeros(n_genes)
		sum_any_init = np.zeros(n_genes)
		sum_init_events = np.zeros(n_genes)
		sum_any_synth = np.zeros(n_genes)
		sum_synth_events = np.zeros(n_genes)
		max_mRNA_counts = np.zeros(n_genes, dtype=np.int64)
		max_protein_counts = np.zeros(n_genes, dtype=np.int64)
		# Time-and-population-averaged copy numbers: accumulate the sum over all
		# timesteps of all cells, then divide by the total timestep count.
		sum_mRNA_counts = np.zeros(n_genes)
		sum_protein_counts = np.zeros(n_genes)
		total_timesteps = 0
		# For near-ubiquitous genes (absent in only a handful of cells), record which
		# cells lacked the mRNA so repeated offenders can be spotted. Once a gene is
		# absent in more than the report cap it can never be "close to 1", so drop its
		# record and disqualify it to keep memory bounded.
		absent_cells_by_gene = {}
		disqualified_absent = np.zeros(n_genes, dtype=bool)
		included_cells = []
		skipped_cells = []

		for i, cell_path in enumerate(cell_paths):
			if i % 100 == 0:
				print('  Cell %d/%d' % (i, len(cell_paths)))
			sim_out = os.path.join(cell_path, 'simOut')
			try:
				mRNA_counts = TableReader(os.path.join(sim_out, 'RNACounts')
					).readColumn('mRNA_cistron_counts')[time_slice][:, def1_indexes]
				init_events = TableReader(os.path.join(sim_out, 'RnapData')
					).readColumn('rna_init_event_per_cistron')[time_slice][:, full_cistron_indexes]
				monomer_counts = TableReader(os.path.join(sim_out, 'MonomerCounts')
					).readColumn('monomerCounts')[time_slice][:, monomer_indexes]
				if has_synth:
					synth_events = TableReader(os.path.join(sim_out, 'TranscriptElongationListener')
						).readColumn('countRnaCistronSynthesized')[time_slice][:, full_cistron_indexes]
			except Exception as e:
				print('  Warning: could not read cell %s: %s' % (cell_path, e))
				skipped_cells.append(cell_path)
				continue

			# Per-generation reductions over the time axis
			init_per_cell = init_events.sum(axis=0)
			present_in_cell = mRNA_counts.sum(axis=0) > 0

			sum_present += present_in_cell                   # def 1
			sum_any_init += init_per_cell > 0                # def 2
			sum_init_events += init_per_cell                 # def 3
			if has_synth:
				synth_per_cell = synth_events.sum(axis=0)
				sum_any_synth += synth_per_cell > 0          # def 4
				sum_synth_events += synth_per_cell           # def 5

			max_mRNA_counts = np.maximum(max_mRNA_counts, mRNA_counts.max(axis=0))
			max_protein_counts = np.maximum(
				max_protein_counts, monomer_counts.max(axis=0))
			sum_mRNA_counts += mRNA_counts.sum(axis=0)
			sum_protein_counts += monomer_counts.sum(axis=0)
			total_timesteps += mRNA_counts.shape[0]

			# Track absent cells for genes still eligible to be "close to 1".
			for gene_index in np.nonzero(~present_in_cell & ~disqualified_absent)[0]:
				gene_index = int(gene_index)
				cells = absent_cells_by_gene.setdefault(gene_index, [])
				cells.append(cell_path)
				if len(cells) > MAX_ABSENT_CELLS_TO_REPORT:
					disqualified_absent[gene_index] = True
					del absent_cells_by_gene[gene_index]

			included_cells.append(cell_path)

		n_cells = len(included_cells)
		if n_cells == 0:
			print('No readable cells found. Skipping.')
			return
		print('Aggregated %d cells.' % n_cells)

		# Aggregate accumulators into the five definitions
		p_present = sum_present / n_cells
		p_any_init = sum_any_init / n_cells
		mean_init_events = sum_init_events / n_cells
		if has_synth:
			p_any_synth = sum_any_synth / n_cells
			mean_synth_events = sum_synth_events / n_cells
		else:
			p_any_synth = np.full(n_genes, np.nan)
			mean_synth_events = np.full(n_genes, np.nan)
		definition_values = [
			p_present, p_any_init, mean_init_events,
			p_any_synth, mean_synth_events]

		# Time-and-population-averaged copy numbers
		mean_mRNA_counts = sum_mRNA_counts / total_timesteps
		mean_protein_counts = sum_protein_counts / total_timesteps

		# Write the main per-gene table
		main_path = os.path.join(plotOutDir, plotOutFileName + '.tsv')
		print('Writing %s' % main_path)
		with open(main_path, 'w') as f:
			writer = csv.writer(f, delimiter='\t')
			writer.writerow(
				['gene_name', 'cistron_name', 'protein_name']
				+ DEFINITION_LABELS
				+ ['max_mRNA_count', 'mean_mRNA_count',
					'max_protein_count', 'mean_protein_count'])
			for i in range(n_genes):
				writer.writerow([
					gene_ids[i], mRNA_cistron_ids[i], monomer_ids[i][:-3],
					p_present[i], p_any_init[i], mean_init_events[i],
					_blank_if_nan(p_any_synth[i]), _blank_if_nan(mean_synth_events[i]),
					max_mRNA_counts[i], mean_mRNA_counts[i],
					max_protein_counts[i], mean_protein_counts[i],
					])

		# Definition comparison: how many genes are subgenerational under each
		# definition, swept across probability thresholds. A gene is subgen for a
		# definition at threshold t if 0 < value < t (the > 0 excludes genes that
		# are never expressed; subsets are nested as t shrinks).
		summary_path = os.path.join(
			plotOutDir, plotOutFileName + '_definition_summary.tsv')
		print('Writing %s' % summary_path)
		header = (
			['definition', 'n_genes_total', 'n_never']
			+ ['n_subgen@%g' % t for t in THRESHOLDS])
		summary_rows = []
		for label, values in zip(DEFINITION_LABELS, definition_values):
			if np.all(np.isnan(values)):
				summary_rows.append(
					[label, n_genes, 'NA'] + ['NA'] * len(THRESHOLDS))
				continue
			n_never = int(np.sum(values == 0))
			subgen_counts = [
				int(np.sum((values > 0) & (values < t))) for t in THRESHOLDS]
			summary_rows.append([label, n_genes, n_never] + subgen_counts)

		with open(summary_path, 'w') as f:
			writer = csv.writer(f, delimiter='\t')
			writer.writerow(header)
			for row in summary_rows:
				writer.writerow(row)

		# Echo the comparison to stdout for quick inspection
		print('\nSubgenerational gene counts (0 < value < threshold):')
		print('\t'.join(header))
		for row in summary_rows:
			print('\t'.join(str(x) for x in row))

		# Near-ubiquitous genes: those present in almost every cell but absent in a
		# few (1 .. MAX_ABSENT_CELLS_TO_REPORT cells). List the offending cells in
		# long format so the same cells recurring across genes -- a sign those cells
		# were unhealthy rather than the gene being genuinely subgenerational -- can
		# be spotted by sorting on absent_cell_path.
		absent_path = os.path.join(
			plotOutDir, plotOutFileName + '_near_ubiquitous_absences.tsv')
		print('Writing %s' % absent_path)
		absent_rows = []
		for gene_index, cells in absent_cells_by_gene.items():
			# A gene that is never present is the opposite of "close to 1"; it only
			# reaches here when the cohort has <= MAX_ABSENT_CELLS_TO_REPORT cells.
			if sum_present[gene_index] == 0:
				continue
			for cell_path in cells:
				absent_rows.append([
					gene_ids[gene_index], mRNA_cistron_ids[gene_index],
					monomer_ids[gene_index][:-3], p_present[gene_index],
					len(cells), cell_path,
					])
		# Sort by cell path (group recurring cells), then by gene for readability.
		absent_rows.sort(key=lambda row: (row[5], row[0]))
		with open(absent_path, 'w') as f:
			writer = csv.writer(f, delimiter='\t')
			writer.writerow([
				'gene_name', 'cistron_name', 'protein_name', 'p_present_def1',
				'n_absent_cells', 'absent_cell_path'])
			for row in absent_rows:
				writer.writerow(row)
		n_reported_genes = len({row[0] for row in absent_rows})
		print('  %d near-ubiquitous genes absent in 1-%d cells (%d gene-cell rows).'
			% (n_reported_genes, MAX_ABSENT_CELLS_TO_REPORT, len(absent_rows)))

		# Provenance metadata: which cells were included, when the analysis ran,
		# when the sims ran, and the git hash of each.
		sim_metadata_path, sim_metadata = _load_sim_metadata(variantDir)
		repo_dir = os.path.dirname(os.path.abspath(__file__))
		run_metadata = {
			'analysis': {
				'script': os.path.basename(__file__),
				'run_time': analysis_run_time,
				'git': _git_info(repo_dir),
				'parameters': {
					'ignore_first_n_gens': IGNORE_FIRST_N_GENS,
					'remove_first_timestep': REMOVE_FIRST_TIMESTEP,
					'thresholds': THRESHOLDS,
					'max_absent_cells_to_report': MAX_ABSENT_CELLS_TO_REPORT,
					'require_complete_lineage': self.REQUIRE_COMPLETE_LINEAGE,
					'countRnaCistronSynthesized_available': has_synth,
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
			'cells': {
				'n_attempted': len(cell_paths),
				'n_included': len(included_cells),
				'n_skipped': len(skipped_cells),
				'included': included_cells,
				'skipped': skipped_cells,
				},
			}
		metadata_path = os.path.join(
			plotOutDir, plotOutFileName + '_run_metadata.json')
		print('Writing %s' % metadata_path)
		with open(metadata_path, 'w') as f:
			json.dump(run_metadata, f, indent=2)


if __name__ == '__main__':
	Plot().cli()
