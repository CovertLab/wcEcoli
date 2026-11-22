"""
Pseudo sampling cells by taking 10k timepoints divided by the number of seeds 
avaialble in the sims. For each side takes (seed number)/10000 time points per seed 
where each snapshot takes completed mRNA counts and protein monomer counts. Also 
records seed number, generation, and time of initiation and division per cell to 
calculate cell age downstream. Ignores the first 8 generation for initiation effects 
and uses only strictly successful seeds 


  <script_name>_mrna.tsv       seed, generation, time_step, cell_birth_time,
                        cell_division_time, then one column per gene holding
                        RNACounts/full_mRNA_cistron_counts (completed
                        transcripts present at that timepoint).

  <script_name>_monomers.tsv   the same five leading columns, then one column per
                        monomer holding MonomerCounts/monomerCounts.

  <script_name>_run_metadata.json

cell_birth_time is Main/time[0] of the containing cell (its true initial
state) and cell_division_time is Main/time[-1], so
cell_division_time - cell_birth_time equals the doubling time All times are in seconds. 
Sampled timesteps never include a cell's row 0 (tables are stacked 
with remove_first=True), so cell_birth_time < time_step <= cell_division_time.
"""

import csv
import os
import pickle

import numpy as np

from models.ecoli.analysis import cohortAnalysisPlot
from models.ecoli.analysis.cohort import subgen_helper_functions as sc
from models.ecoli.analysis.cohort.subgen_subsample_mrna_protein_10k import (
	_sampled_rows, BATCH_SIZE)
from wholecell.analysis.analysis_tools import read_stacked_columns
from wholecell.io.tablereader import TableReader

IGNORE_FIRST_N_GENS = sc.IGNORE_FIRST_N_GENS
SEED_RANGE = sc.SEED_RANGE
TIMEPOINTS_TO_SAMPLE = sc.TIMEPOINTS_TO_SAMPLE
NUMPY_SEED = 0

META_COLS = [
	'seed', 'generation', 'time_step', 'cell_birth_time', 'cell_division_time']


def _write_table(path, id_cols, meta, counts):
	"""Write the five metadata columns followed by integer count columns."""
	print('Writing %s' % path)
	counts = np.rint(counts).astype(np.int64)
	with open(path, 'w', newline='') as f:
		writer = csv.writer(f, delimiter='\t', lineterminator='\n')
		writer.writerow(META_COLS + list(id_cols))
		for i in range(counts.shape[0]):
			writer.writerow(
				[row[i] for row in meta] + counts[i].tolist())


class Plot(cohortAnalysisPlot.CohortAnalysisPlot):
	def do_plot(self, variantDir, plotOutDir, plotOutFileName, simDataFile,
			validationDataFile, metadata):
		with open(simDataFile, 'rb') as f:
			sim_data = pickle.load(f)
		if self.ap.n_generation <= IGNORE_FIRST_N_GENS:
			print('Skipping analysis - not enough generations run.')
			return
		cell_paths = self.ap.get_cells(
			generation=np.arange(IGNORE_FIRST_N_GENS, self.ap.n_generation),
			seed=SEED_RANGE, only_successful=True)
		if len(cell_paths) == 0:
			print('No cells found. Skipping.')
			return

		# Strict-successful seeds
		success = sc.compute_seed_success(self.ap, self.ap.n_generation)
		seed_range_set = set(SEED_RANGE.tolist())
		seeds_to_sample = sorted(
			s for s in success['successful_seeds'] if s in seed_range_set)
		if not seeds_to_sample:
			print('No successful seeds found. Skipping.')
			return
		sample_per_seed = sc.sample_per_seed(len(seeds_to_sample))
		print('Sampling %d timepoints from each of %d successful seeds '
			'(target %d total).'
			% (sample_per_seed, len(seeds_to_sample), TIMEPOINTS_TO_SAMPLE))

		# Gene set: protein-coding mRNA cistrons, in cistron order.
		cistron_data = sim_data.process.transcription.cistron_data
		cistron_id_to_protein_id = {
			protein['cistron_id']: protein['id']
			for protein in sim_data.process.translation.monomer_data
			}
		mRNA_cistron_ids = [
			cistron_id for cistron_id in cistron_data['id']
			if cistron_id in cistron_id_to_protein_id]
		cistron_id_to_gene_id = {
			cistron['id']: cistron['gene_id'] for cistron in cistron_data}
		gene_ids = [cistron_id_to_gene_id[c] for c in mRNA_cistron_ids]
		monomer_ids = [cistron_id_to_protein_id[c] for c in mRNA_cistron_ids]
		n_genes = len(gene_ids)

		# Subcolumn indices into the RNACounts and MonomerCounts tables.
		RNA_reader = TableReader(
			os.path.join(cell_paths[0], 'simOut', 'RNACounts'))
		mRNA_ids_counts_table = RNA_reader.readAttribute('mRNA_cistron_ids')
		RNA_reader.close()
		mRNA_cistron_id_to_index = {
			cistron_id: i for i, cistron_id in enumerate(mRNA_ids_counts_table)}
		mRNA_cistron_indices = np.array([
			mRNA_cistron_id_to_index[c] for c in mRNA_cistron_ids])

		monomer_reader = TableReader(
			os.path.join(cell_paths[0], 'simOut', 'MonomerCounts'))
		monomer_ids_counts_table = monomer_reader.readAttribute('monomerIds')
		monomer_reader.close()
		monomer_id_to_index = {
			monomer_id: i for i, monomer_id in enumerate(monomer_ids_counts_table)}
		monomer_indexes = np.array([
			monomer_id_to_index[m] for m in monomer_ids])

		total_seed_ids = []
		total_generation_ids = []
		total_time_steps = []
		total_birth_times = []
		total_division_times = []
		mrna_blocks = []
		monomer_blocks = []

		# Seed the RNG once (not per seed) so each seed draws independently.
		# Keep this exactly where subgen_subsample_mrna_protein_10k.py has it,
		# with subsample_seed_timepoints as the only RNG consumer, so the
		# sampled timepoints match that script's.
		np.random.seed(NUMPY_SEED)
		for seed in SEED_RANGE:
			cell_paths_per_seed = self.ap.get_cells(
				generation=np.arange(IGNORE_FIRST_N_GENS, self.ap.n_generation),
				seed=[seed], only_successful=True)
			if seed not in success['successful_seeds']:
				continue
			if len(cell_paths_per_seed) == 0:
				continue

			sample = sc.subsample_seed_timepoints(
				cell_paths_per_seed, sample_per_seed)
			if sample is None:
				continue
			random_time_indices = sample['time_indices']
			gen_index = sample['gen_index']

			generation_ids = np.array([
				sc.parse_cell_id(cp)[1] for cp in cell_paths_per_seed
				])[gen_index]

			# True birth (row 0) and division (last row) of each cell; no
			# remove_first, since row 0 IS the birth state.
			birth_times = read_stacked_columns(
				cell_paths_per_seed, 'Main', 'time',
				fun=lambda x: x[0]).flatten()
			division_times = read_stacked_columns(
				cell_paths_per_seed, 'Main', 'time',
				fun=lambda x: x[-1]).flatten()

			# Per-cell row counts (remove_first) map global sampled indices
			# back into per-batch reads of the heavy tables.
			per_cell_rows = read_stacked_columns(
				cell_paths_per_seed, 'Main', 'time', remove_first=True,
				fun=lambda x: len(x)).flatten().astype(int)
			cell_row_offsets = np.concatenate([[0], np.cumsum(per_cell_rows)])

			full_mRNA_counts = _sampled_rows(
				cell_paths_per_seed, 'RNACounts', 'full_mRNA_cistron_counts',
				random_time_indices, cell_row_offsets, n_genes,
				mRNA_cistron_indices)
			monomer_counts = _sampled_rows(
				cell_paths_per_seed, 'MonomerCounts', 'monomerCounts',
				random_time_indices, cell_row_offsets, n_genes,
				monomer_indexes)

			n = len(random_time_indices)
			total_seed_ids.extend([int(seed)] * n)
			total_generation_ids.extend(generation_ids.tolist())
			total_time_steps.extend(sample['time_steps'].tolist())
			total_birth_times.extend(birth_times[gen_index].tolist())
			total_division_times.extend(division_times[gen_index].tolist())
			mrna_blocks.append(full_mRNA_counts)
			monomer_blocks.append(monomer_counts)

		if not mrna_blocks:
			print('No timepoints sampled. Skipping.')
			return

		meta = [total_seed_ids, total_generation_ids, total_time_steps,
			total_birth_times, total_division_times]
		_write_table(
			os.path.join(plotOutDir, plotOutFileName + '_mrna.tsv'),
			gene_ids, meta, np.vstack(mrna_blocks))
		_write_table(
			os.path.join(plotOutDir, plotOutFileName + '_monomers.tsv'),
			monomer_ids, meta, np.vstack(monomer_blocks))

		n_rows = len(total_time_steps)
		sim_metadata_path, sim_metadata = sc.load_sim_metadata(variantDir)
		sc.write_run_metadata(
			os.path.join(plotOutDir, plotOutFileName + '_run_metadata.json'),
			script=os.path.basename(__file__),
			parameters={
				'ignore_first_n_gens': IGNORE_FIRST_N_GENS,
				'max_doubling_min': sc.MAX_DOUBLING_MIN,
				'timepoints_to_sample': TIMEPOINTS_TO_SAMPLE,
				'sample_per_seed': sample_per_seed,
				'batch_size': BATCH_SIZE,
				'numpy_seed': NUMPY_SEED,
				'mrna_column': 'RNACounts/full_mRNA_cistron_counts',
				'monomer_column': 'MonomerCounts/monomerCounts',
				'cell_birth_time': 'Main/time[0]',
				'cell_division_time': 'Main/time[-1]',
				},
			sim_metadata_path=sim_metadata_path,
			sim_metadata=sim_metadata,
			extra={
				'seeds': {
					'n_successful': len(seeds_to_sample),
					'successful_seeds': seeds_to_sample,
					},
				'timepoints': {'n_sampled': n_rows},
				'genes': {'n_genes': n_genes},
				})


if __name__ == '__main__':
	Plot().cli()
