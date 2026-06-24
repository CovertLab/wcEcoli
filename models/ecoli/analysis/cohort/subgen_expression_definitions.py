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

import numpy as np

from models.ecoli.analysis import cohortAnalysisPlot
from wholecell.io.tablereader import TableReader


IGNORE_FIRST_N_GENS = 8
# Drop the inherited boundary timestep from each generation when True. Default
# False to match subgenerational_expression_table.py for direct comparison.
REMOVE_FIRST_TIMESTEP = False
# Probability filters for the subgenerational classification sweep.
THRESHOLDS = [1.0, 0.99, 0.95]

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


class Plot(cohortAnalysisPlot.CohortAnalysisPlot):
	def do_plot(self, variantDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
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
		print('Analyzing %d cells...' % len(cell_paths))

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
		n_cells = 0

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
				continue

			# Per-generation reductions over the time axis
			init_per_cell = init_events.sum(axis=0)

			sum_present += mRNA_counts.sum(axis=0) > 0       # def 1
			sum_any_init += init_per_cell > 0                # def 2
			sum_init_events += init_per_cell                 # def 3
			if has_synth:
				synth_per_cell = synth_events.sum(axis=0)
				sum_any_synth += synth_per_cell > 0          # def 4
				sum_synth_events += synth_per_cell           # def 5

			max_mRNA_counts = np.maximum(max_mRNA_counts, mRNA_counts.max(axis=0))
			max_protein_counts = np.maximum(
				max_protein_counts, monomer_counts.max(axis=0))
			n_cells += 1

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

		# Write the main per-gene table
		main_path = os.path.join(plotOutDir, plotOutFileName + '.tsv')
		print('Writing %s' % main_path)
		with open(main_path, 'w') as f:
			writer = csv.writer(f, delimiter='\t')
			writer.writerow(
				['gene_name', 'cistron_name', 'protein_name']
				+ DEFINITION_LABELS
				+ ['max_mRNA_count', 'max_protein_count'])
			for i in range(n_genes):
				writer.writerow([
					gene_ids[i], mRNA_cistron_ids[i], monomer_ids[i][:-3],
					p_present[i], p_any_init[i], mean_init_events[i],
					_blank_if_nan(p_any_synth[i]), _blank_if_nan(mean_synth_events[i]),
					max_mRNA_counts[i], max_protein_counts[i],
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


if __name__ == '__main__':
	Plot().cli()
