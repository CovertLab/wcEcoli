"""
Report per-cell health metrics across a cohort so that outlier ("unhealthy")
cells can be identified and cross-referenced against other analyses.

For every cell-generation the script records:

	- doubling_time_min   -- generation length (Main/time span), minutes
	- mean_growth_rate    -- Mass/instantaneous_growth_rate, time-averaged
	- mean_cell_mass_fg   -- Mass/cellMass, time-averaged
	- final_cell_mass_fg  -- Mass/cellMass, last timestep
	- mean_dry_mass_fg    -- Mass/dryMass, time-averaged
	- mean_active_rnap    -- UniqueMoleculeCounts active_RNAP, time-averaged
	- mean_active_ribosome-- UniqueMoleculeCounts active_ribosome, time-averaged
	- mean_ppgpp_conc     -- GrowthLimits/ppgpp_conc, time-averaged

A healthy cell divides on schedule with a stable mass, high polymerase/ribosome
activity, and low ppGpp; unhealthy cells show a long doubling time, low growth
rate/mass, depleted RNAP/ribosome activity, and/or elevated ppGpp (stringent
response). The companion subgen_definitions.py flags near-ubiquitous
genes that are absent only in a handful of cells; joining that report to this
table on the cell path tests whether those absences concentrate in unhealthy
cells rather than reflecting genuine subgenerational regulation.

skip the first IGNORE_FIRST_N_GENS generations, and cells from successful seeds only) 
so the two tables join directly on the cell path.
"""

import os
import csv
from datetime import datetime

import numpy as np

from models.ecoli.analysis import cohortAnalysisPlot
from models.ecoli.analysis.cohort import subgen_helper_functions as sc
from wholecell.io.tablereader import TableReader


IGNORE_FIRST_N_GENS = sc.IGNORE_FIRST_N_GENS
REMOVE_FIRST_TIMESTEP = False

# Sibling report written by subgen_definitions.py: one row per
# (near-ubiquitous gene, absent cell). When present in the same plotOut
# directory, each cell is annotated with how many such genes it is missing.
NEAR_UBIQUITOUS_ABSENCE_FILE = (
	'subgen_definitions_all_near_ubiquitous_absences.tsv')

# Per-cell metric columns, in output order.
METRIC_LABELS = [
	'doubling_time_min',
	'mean_growth_rate',
	'mean_cell_mass_fg',
	'final_cell_mass_fg',
	'mean_dry_mass_fg',
	'mean_active_rnap',
	'mean_active_ribosome',
	'mean_ppgpp_conc',
	]



def _load_absence_counts(plot_out_dir):
	"""Count near-ubiquitous absences per cell from the subgen report, if present.

	subgen_definitions.py writes a long-format table with one row per
	(gene, absent_cell_path); the per-cell count is how many near-ubiquitous genes
	that cell is missing. Return (counts_by_cell_path, source_path), or ({}, None)
	when the report is absent.
	"""
	path = os.path.join(plot_out_dir, NEAR_UBIQUITOUS_ABSENCE_FILE)
	if not os.path.isfile(path):
		return {}, None
	counts = {}
	with open(path) as f:
		reader = csv.DictReader(f, delimiter='\t')
		for row in reader:
			cell_path = row['absent_cell_path']
			counts[cell_path] = counts.get(cell_path, 0) + 1
	return counts, path


class Plot(cohortAnalysisPlot.CohortAnalysisPlot):
	def do_plot(self, variantDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		analysis_run_time = datetime.now().isoformat(timespec='seconds')

		# Ignore data from a predefined number of generations per seed
		if self.ap.n_generation <= IGNORE_FIRST_N_GENS:
			print('Skipping analysis - not enough generations run.')
			return

		cell_paths = self.ap.get_cells(
			generation=np.arange(IGNORE_FIRST_N_GENS, self.ap.n_generation),
			only_successful=True)
		print('Analyzing %d cells...' % len(cell_paths))

		if len(cell_paths) == 0:
			print('No successful cells found. Skipping analysis.')
			return

		# Resolve the active_RNAP / active_ribosome indices from the first cell.
		first_sim_out = os.path.join(cell_paths[0], 'simOut')
		unique_reader = TableReader(
			os.path.join(first_sim_out, 'UniqueMoleculeCounts'))
		unique_molecule_ids = unique_reader.readAttribute('uniqueMoleculeIds')
		unique_reader.close()
		rnap_index = unique_molecule_ids.index('active_RNAP')
		ribosome_index = unique_molecule_ids.index('active_ribosome')

		time_slice = slice(1, None) if REMOVE_FIRST_TIMESTEP else slice(None)

		rows = []
		skipped_cells = []
		for i, cell_path in enumerate(cell_paths):
			if i % 100 == 0:
				print('  Cell %d/%d' % (i, len(cell_paths)))
			sim_out = os.path.join(cell_path, 'simOut')
			try:
				time = TableReader(os.path.join(sim_out, 'Main')
					).readColumn('time')[time_slice]
				cell_mass = TableReader(os.path.join(sim_out, 'Mass')
					).readColumn('cellMass')[time_slice]
				dry_mass = TableReader(os.path.join(sim_out, 'Mass')
					).readColumn('dryMass')[time_slice]
				growth_rate = TableReader(os.path.join(sim_out, 'Mass')
					).readColumn('instantaneous_growth_rate')[time_slice]
				unique_counts = TableReader(os.path.join(sim_out,
					'UniqueMoleculeCounts')).readColumn(
					'uniqueMoleculeCounts')[time_slice]
				ppgpp_conc = TableReader(os.path.join(sim_out, 'GrowthLimits')
					).readColumn('ppgpp_conc')[time_slice]
			except Exception as e:
				print('  Warning: could not read cell %s: %s' % (cell_path, e))
				skipped_cells.append(cell_path)
				continue

			seed, generation = sc.parse_cell_id(cell_path)
			# Doubling time in minutes: span of the generation's time axis.
			doubling_time_min = (time[-1] - time[0]) / 60.0
			# instantaneous_growth_rate is NaN at the first timestep (no prior step
			# to difference), so average over the finite values.
			rows.append([
				cell_path, seed, generation,
				doubling_time_min,
				float(np.nanmean(growth_rate)),
				float(cell_mass.mean()),
				float(cell_mass[-1]),
				float(dry_mass.mean()),
				float(unique_counts[:, rnap_index].mean()),
				float(unique_counts[:, ribosome_index].mean()),
				float(ppgpp_conc.mean()),
				])

		n_cells = len(rows)
		if n_cells == 0:
			print('No readable cells found. Skipping.')
			return
		print('Aggregated %d cells.' % n_cells)

		# Optional join: annotate each cell with how many near-ubiquitous genes it
		# is missing, from the sibling subgen report (if present in this plotOut).
		absence_counts, absence_path = _load_absence_counts(plotOutDir)
		if absence_path is None:
			print('Note: %s not found in plotOut; n_near_ubiquitous_absences left '
				'blank.' % NEAR_UBIQUITOUS_ABSENCE_FILE)
		else:
			print('Joining near-ubiquitous absences from %s' % absence_path)
			# Surface the suspected unhealthy cells first for easy inspection.
			rows.sort(key=lambda row: -absence_counts.get(row[0], 0))

		# Write the per-cell health table
		main_path = os.path.join(plotOutDir, plotOutFileName + '.tsv')
		print('Writing %s' % main_path)
		with open(main_path, 'w') as f:
			writer = csv.writer(f, delimiter='\t')
			writer.writerow(['cell_path', 'seed', 'generation'] + METRIC_LABELS
				+ ['n_near_ubiquitous_absences'])
			for row in rows:
				absences = absence_counts.get(row[0], 0) if absence_path else ''
				writer.writerow(row + [absences])

		# Cohort-level distribution of each metric, so outliers are obvious.
		metrics = np.array([row[3:] for row in rows], dtype=np.float64)
		summary_path = os.path.join(
			plotOutDir, plotOutFileName + '_metric_summary.tsv')
		print('Writing %s' % summary_path)
		percentiles = [1, 5, 50, 95, 99]
		header = (['metric', 'mean', 'std', 'min', 'max']
			+ ['p%d' % p for p in percentiles])
		summary_rows = []
		for j, label in enumerate(METRIC_LABELS):
			col = metrics[:, j]
			summary_rows.append([label,
				float(col.mean()), float(col.std()),
				float(col.min()), float(col.max())]
				+ [float(np.percentile(col, p)) for p in percentiles])
		with open(summary_path, 'w') as f:
			writer = csv.writer(f, delimiter='\t')
			writer.writerow(header)
			for row in summary_rows:
				writer.writerow(row)

		# Echo the distribution to stdout for a quick scan.
		print('\nPer-cell health metric distribution (%d cells):' % n_cells)
		print('\t'.join(header))
		for row in summary_rows:
			print(row[0] + '\t' + '\t'.join('%.4g' % x for x in row[1:]))

		# If the join ran, persist and print the flagged cells against their health
		# so the unhealthy-cell hypothesis can be judged (cohort medians alongside).
		if absence_path is not None:
			flagged = [row for row in rows if absence_counts.get(row[0], 0) > 0]
			medians = np.median(metrics, axis=0)

			# Durable per-flagged-cell summary with a trailing cohort-median row.
			flagged_path = os.path.join(
				plotOutDir, plotOutFileName + '_flagged_summary.tsv')
			print('Writing %s' % flagged_path)
			with open(flagged_path, 'w') as f:
				writer = csv.writer(f, delimiter='\t')
				writer.writerow(
					['cell_path', 'seed', 'generation',
						'n_near_ubiquitous_absences'] + METRIC_LABELS)
				for row in flagged:
					writer.writerow(
						[row[0], row[1], row[2], absence_counts[row[0]]] + row[3:])
				writer.writerow(
					['MEDIAN(all)', '', '', '']
					+ [float(m) for m in medians])

			# Echo the head of the comparison for a quick scan.
			print('\n%d cells missing >=1 near-ubiquitous gene. Top 15 vs cohort '
				'median:' % len(flagged))
			print('seed/gen\tn_absent\tdoubling_min\tgrowth_rate\tribosome\tppgpp')
			for row in flagged[:15]:
				print('%s/%s\t%d\t%.1f\t%.2e\t%.0f\t%.1f' % (
					row[1], row[2], absence_counts[row[0]],
					row[3], row[4], row[9], row[10]))
			print('MEDIAN(all)\t-\t%.1f\t%.2e\t%.0f\t%.1f' % (
				medians[0], medians[1], medians[6], medians[7]))

		# Full-run seed completeness: how far each seed got before it
		# died or reached the final generation. Independent of IGNORE_FIRST_N_GENS,
		# so seeds that die in the first few generations are still captured.
		sim_metadata_path, sim_metadata = sc.load_sim_metadata(variantDir)
		final_gen = self.ap.n_generation - 1
		successful_gens_by_seed = {}
		for path in self.ap.get_cells(only_successful=True):
			seed, gen = sc.parse_cell_id(path)
			successful_gens_by_seed.setdefault(seed, set()).add(gen)
		seeds_present = set()
		for path in self.ap.get_cells(only_successful=False):
			seed, _ = sc.parse_cell_id(path)
			seeds_present.add(seed)

		completeness_rows = []
		for seed in sorted(seeds_present):
			gens = successful_gens_by_seed.get(seed, set())
			if gens:
				first_gen, last_gen = min(gens), max(gens)
				missing_gens = sorted(set(range(first_gen, last_gen + 1)) - gens)
			else:
				first_gen = last_gen = -1
				missing_gens = []
			completeness_rows.append([
				seed, first_gen, last_gen, len(gens),
				last_gen == final_gen,
				';'.join(str(g) for g in missing_gens),
				])

		completeness_path = os.path.join(
			plotOutDir, plotOutFileName + '_seed_completeness.tsv')
		print('Writing %s' % completeness_path)
		with open(completeness_path, 'w') as f:
			writer = csv.writer(f, delimiter='\t')
			writer.writerow(['seed', 'first_gen', 'last_gen', 'n_gens_present',
				'reached_final', 'missing_gens'])
			for row in completeness_rows:
				writer.writerow(row)

		n_reached_final = sum(1 for row in completeness_rows if row[4])
		n_died_early = len(completeness_rows) - n_reached_final
		total_init_sims = sim_metadata.get('total_init_sims')
		n_never_ran = (total_init_sims - len(seeds_present)
			if total_init_sims is not None else None)
		print('\nSeed completeness (final generation = %d):' % final_gen)
		print('  seeds present            : %d' % len(seeds_present))
		if n_never_ran is not None:
			print('  seeds that never ran     : %d (of %d expected)'
				% (n_never_ran, total_init_sims))
		print('  reached final generation : %d' % n_reached_final)
		print('  died before final gen    : %d' % n_died_early)
		if n_died_early:
			incomplete = sorted(
				(row for row in completeness_rows if not row[4]),
				key=lambda row: row[2])
			print('  incomplete seeds (seed@last_gen): '
				+ ', '.join('%s@%d' % (row[0], row[2]) for row in incomplete))

		# Provenance metadata: which cells were included, when the analysis ran,
		# when the sims ran, and the git hash of each. `near_ubiquitous_join` is
		# specific to this script, so it rides along as an `extra` block.
		sc.write_run_metadata(
			os.path.join(plotOutDir, plotOutFileName + '_run_metadata.json'),
			os.path.basename(__file__),
			{
				'ignore_first_n_gens': IGNORE_FIRST_N_GENS,
				'remove_first_timestep': REMOVE_FIRST_TIMESTEP,
				},
			sim_metadata_path=sim_metadata_path,
			sim_metadata=sim_metadata,
			run_time=analysis_run_time,
			extra={
				'near_ubiquitous_join': {
					'source': absence_path,
					'n_cells_with_absences': sum(
						1 for row in rows if absence_counts.get(row[0], 0) > 0),
					},
				'cells': {
					'n_attempted': len(cell_paths),
					'n_included': n_cells,
					'n_skipped': len(skipped_cells),
					'skipped': skipped_cells,
					},
				'seeds': {
					'final_generation': final_gen,
					'n_seeds_present': len(seeds_present),
					'n_seeds_never_ran': n_never_ran,
					'n_reached_final': n_reached_final,
					'n_died_before_final': n_died_early,
					},
				})


if __name__ == '__main__':
	Plot().cli()
