"""
Cohort analysis: export subsampled monomer counts + cell dry mass.

Randomly samples TIMEPOINTS_TO_SAMPLE total timepoints spread evenly across
seeds (SAMPLE_PER_SEED per seed), matching the pattern in
subsample_timepoints_for_scRNAseq_comp.py.

Produces two TSV files:
  <plotOutFileName>_monomer_counts_dry_mass.tsv
      One row per sampled timepoint.
      Columns: seed, time_step, generation, generation_start_time, dry_mass_fg,
               <monomer_id_0>, <monomer_id_1>, ...

  <plotOutFileName>_monomer_mws.tsv
      One row per monomer.
      Columns: monomer_id, mw_g_per_mol

dry_mass_fg is in femtograms (native listener unit).
mw_g_per_mol is from sim_data.process.translation.monomer_data['mw'] in g/mol.
Monomer column order is identical in both files (order from MonomerCounts listener).
"""

import pickle
import os

# noinspection PyUnresolvedReferences
import numpy as np

import csv

from models.ecoli.analysis import cohortAnalysisPlot
from models.ecoli.analysis.cohort import subgen_common as sc
from wholecell.analysis.analysis_tools import read_stacked_columns
from wholecell.io.tablereader import TableReader

IGNORE_FIRST_N_GENS = sc.IGNORE_FIRST_N_GENS
SEED_RANGE = sc.SEED_RANGE
TIMEPOINTS_TO_SAMPLE = 10000
SAMPLE_PER_SEED = TIMEPOINTS_TO_SAMPLE // len(SEED_RANGE)
BATCH_SIZE = 100  # max cell paths passed to TableReader at once


class Plot(cohortAnalysisPlot.CohortAnalysisPlot):
	def do_plot(self, variantDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		with open(simDataFile, 'rb') as f:
			sim_data = pickle.load(f)
			# Ignore data from predefined number of generations per seed
		if self.ap.n_generation <= IGNORE_FIRST_N_GENS:
			print('Skipping analysis - not enough generations run.')
			return
		cell_paths = self.ap.get_cells(
			generation=np.arange(IGNORE_FIRST_N_GENS, self.ap.n_generation), seed=SEED_RANGE,
			only_successful=True)

		if len(cell_paths) == 0:
			print('No valid cell paths found for this variant. Skipping analysis.')
			return

		# Strict-successful lineages (completed every generation and no cell at
		# the 180-min doubling cap).
		success = sc.compute_lineage_success(self.ap, self.ap.n_generation)

		print('Analyzing %d cells...' % len(cell_paths))

		# Read monomer IDs in the order used by the MonomerCounts listener
		monomer_reader = TableReader(
			os.path.join(cell_paths[0], 'simOut', 'MonomerCounts'))
		monomer_ids = monomer_reader.readAttribute('monomerIds')
		monomer_reader.close()

		# Align molecular weights to listener order
		monomer_data = sim_data.process.translation.monomer_data
		id_to_mw = {m['id']: m['mw'] for m in monomer_data}
		monomer_mws = np.array([id_to_mw[mid] for mid in monomer_ids])

		total_seed_ids = []
		total_random_time_steps = []
		total_gen_start_times = []
		total_generation_ids = []
		total_monomer_counts = np.empty((0, len(monomer_ids)), dtype=np.float64)
		total_dry_mass = []

		# Seed the RNG once (not per seed) so each lineage draws independently.
		np.random.seed(0)
		for seed in SEED_RANGE:
			cell_paths_per_seed = self.ap.get_cells(
				generation=np.arange(IGNORE_FIRST_N_GENS, self.ap.n_generation), seed=[seed],
				only_successful=True)

			if len(cell_paths_per_seed) == 0:
				continue

			if seed not in success['successful_seeds']:
				continue

			# Draw random timepoints and align them to generation starts (shared
			# helper reads Main/time once; clamps the draw to available timesteps).
			sample = sc.subsample_seed_timepoints(cell_paths_per_seed, SAMPLE_PER_SEED)
			if sample is None:
				continue
			random_time_indices = sample['time_indices']
			random_time_steps = sample['time_steps']
			aligned_start_times = sample['gen_start_times']
			actual_sample_size = len(random_time_indices)

			# Map each sampled time step to its generation number.
			gen_indices_per_cell = [
				int(os.path.basename(os.path.dirname(cp))[-6:])
				for cp in cell_paths_per_seed
				]
			generation_ids = np.array(gen_indices_per_cell)[sample['gen_index']]

			# Per-cell row counts (remove_first) give the row offsets used to map
			# global sampled indices back into per-batch reads of the heavy tables.
			per_cell_rows = read_stacked_columns(
				cell_paths_per_seed, 'Main', 'time', remove_first=True,
				fun=lambda x: len(x)).flatten().astype(int)
			cell_row_offsets = np.concatenate([[0], np.cumsum(per_cell_rows)])

			# Read monomer counts and dry mass in batches to bound memory; only
			# read a batch if it contains at least one sampled timestep.
			monomer_counts_sampled = np.empty((actual_sample_size, len(monomer_ids)))
			dry_mass_sampled = np.empty(actual_sample_size)

			for b in range(0, len(cell_paths_per_seed), BATCH_SIZE):
				batch = cell_paths_per_seed[b:b + BATCH_SIZE]
				row_start = cell_row_offsets[b]
				row_end = cell_row_offsets[b + len(batch)]

				mask = (random_time_indices >= row_start) & (random_time_indices < row_end)
				if not mask.any():
					continue

				local_idx = random_time_indices[mask] - row_start
				monomer_counts_sampled[mask] = read_stacked_columns(
					batch, 'MonomerCounts', 'monomerCounts',
					remove_first=True)[local_idx]
				dry_mass_sampled[mask] = read_stacked_columns(
					batch, 'Mass', 'dryMass',
					remove_first=True).flatten()[local_idx]

			total_seed_ids = total_seed_ids + [seed] * len(random_time_steps)
			total_random_time_steps = total_random_time_steps + random_time_steps.tolist()
			total_gen_start_times = total_gen_start_times + aligned_start_times.tolist()
			total_generation_ids = total_generation_ids + generation_ids.tolist()
			total_monomer_counts = np.vstack([total_monomer_counts, monomer_counts_sampled])
			total_dry_mass = total_dry_mass + dry_mass_sampled.tolist()

		# Write MW TSV
		with open(os.path.join(plotOutDir, plotOutFileName + '_monomer_mws.tsv'), 'w') as f:
			writer = csv.writer(f, delimiter='\t')
			writer.writerow(['monomer_id', 'mw_g_per_mol'])
			for mid, mw in zip(monomer_ids, monomer_mws):
				writer.writerow([mid, mw])

		# Write data to table
		table_cols = (
			['seed', 'time_step', 'generation', 'generation_start_time', 'dry_mass_fg']
			+ list(monomer_ids))
		with open(os.path.join(plotOutDir, plotOutFileName + '_monomer_counts_dry_mass.tsv'), 'w') as f:
			writer = csv.writer(f, delimiter='\t')
			writer.writerow(table_cols)
			for i in np.arange(0, len(total_random_time_steps)):
				seed = total_seed_ids[i]
				time_step = total_random_time_steps[i]
				generation = total_generation_ids[i]
				gen_start = total_gen_start_times[i]
				dry_mass = total_dry_mass[i]
				monomer_counts_row = total_monomer_counts[i].tolist()
				counts_row = [seed] + [time_step] + [generation] + [gen_start] + [dry_mass] + monomer_counts_row
				writer.writerow(counts_row)


if __name__ == '__main__':
	Plot().cli()
