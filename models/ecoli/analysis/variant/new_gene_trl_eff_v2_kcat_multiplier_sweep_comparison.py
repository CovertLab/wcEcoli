"""
Completion / doubling time / dry mass for the
new_gene_trl_eff_v2_kcat_multiplier_sweep, vs kcat multiplier, one line per
GFP translation-efficiency level, faceted by estimator (gap | smoothed_max).

Rows = completion rate, avg doubling time, avg dry mass; columns = estimator.
The multiplier = 1.0 line edge reproduces GFP-only; the trl-eff = 0 line is the
estimator's pure multiplier sweep (no translation burden).  Also writes
doubling_times_by_lineage.csv.  Reads only recorded listeners.
"""

import csv
import os

import numpy as np
from matplotlib import pyplot as plt

from models.ecoli.analysis import variantAnalysisPlot
from models.ecoli.sim.variants.new_gene_trl_eff_v2_kcat_multiplier_sweep import (
	ESTIMATORS, ESTIMATOR_SHORT, TRL_EFF_VALUES, MULTIPLIERS,
	is_wildtype, variant_to_cell)
from wholecell.analysis.analysis_tools import exportFigure
from wholecell.io.tablereader import TableReader


MAX_DOUBLING_TIME_MIN = 300
TRL_COLORS = plt.cm.viridis(np.linspace(0, 0.9, len(TRL_EFF_VALUES)))


def _seed_of(cell_path):
	parts = cell_path.rstrip('/').split(os.sep)
	for i, part in enumerate(parts):
		if part.startswith('generation_') and i > 0:
			return parts[i - 1]
	return parts[-3] if len(parts) >= 3 else cell_path


class Plot(variantAnalysisPlot.VariantAnalysisPlot):
	def do_plot(self, inputDir, plotOutDir, plotOutFileName, simDataFile,
				validationDataFile, metadata):
		variant_indexes = self.ap.get_variants()
		n_total_gens = self.ap.n_generation
		ignore_first_n_gens = max(n_total_gens - 8, 0)
		window = np.arange(ignore_first_n_gens, n_total_gens)

		completion, avg_dt, avg_mass = {}, {}, {}
		lineage_dts = []
		for vi in variant_indexes:
			n0 = len(self.ap.get_cells(
				variant=[vi], generation=[0], only_successful=True))
			nlast = len(self.ap.get_cells(
				variant=[vi], generation=[n_total_gens - 1],
				only_successful=True))
			completion[vi] = (nlast / n0) if n0 > 0 else 0.0

			mean_dts = np.full(n_total_gens, np.nan)
			seed_dts = {}
			for gen in range(n_total_gens):
				times = []
				for cp in self.ap.get_cells(
						variant=[vi], generation=[gen], only_successful=True):
					try:
						t = TableReader(os.path.join(cp, 'simOut', 'Main')
							).readColumn('time', squeeze=True)
					except Exception:
						continue
					dt = (t[-1] - t[0]) / 60.
					if dt <= MAX_DOUBLING_TIME_MIN:
						times.append(dt)
						seed_dts.setdefault(_seed_of(cp), {})[gen] = dt
				if times:
					mean_dts[gen] = float(np.mean(times))
			later = mean_dts[ignore_first_n_gens:]
			avg_dt[vi] = (float(np.nanmean(later))
				if np.any(np.isfinite(later)) else np.nan)
			for seed in sorted(seed_dts):
				row = {'variant': vi, 'seed': seed}
				for gen in range(n_total_gens):
					row[f'gen_{gen}_dt_min'] = seed_dts[seed].get(gen, np.nan)
				lineage_dts.append(row)

			masses = []
			for cp in self.ap.get_cells(
					variant=[vi], generation=window, only_successful=True):
				try:
					dm = TableReader(os.path.join(cp, 'simOut', 'Mass')
						).readColumn('dryMass', squeeze=True)
					masses.append(float(np.mean(dm)))
				except Exception:
					pass
			avg_mass[vi] = float(np.mean(masses)) if masses else np.nan

		if lineage_dts:
			gen_cols = [f'gen_{g}_dt_min' for g in range(n_total_gens)]
			csv_path = os.path.join(plotOutDir, 'doubling_times_by_lineage.csv')
			with open(csv_path, 'w', newline='', encoding='utf-8') as fh:
				w = csv.DictWriter(fh, fieldnames=['variant', 'seed'] + gen_cols)
				w.writeheader()
				w.writerows(lineage_dts)
			print(f'Wrote {csv_path}')

		wt = next((vi for vi in variant_indexes if is_wildtype(vi)), None)

		def _series(metric):
			s = {}
			for vi in variant_indexes:
				if is_wildtype(vi):
					continue
				e, t, m = variant_to_cell(vi)
				s.setdefault((e, t), [np.nan] * len(MULTIPLIERS))[m] = \
					metric.get(vi, np.nan)
			return s

		x = np.array(MULTIPLIERS)
		n_est = len(ESTIMATORS)
		fig, axes = plt.subplots(3, n_est, figsize=(6 * n_est, 12),
			squeeze=False)

		def _facet(row, metric, ylabel, wt_fmt):
			s = _series(metric)
			for e in range(n_est):
				ax = axes[row][e]
				for t in range(len(TRL_EFF_VALUES)):
					ax.plot(x, s.get((e, t), [np.nan] * len(MULTIPLIERS)),
						marker='o', lw=1.3, color=TRL_COLORS[t],
						label=f'trl {TRL_EFF_VALUES[t]}')
				if wt is not None and np.isfinite(metric.get(wt, np.nan)):
					ax.axhline(metric[wt], color='k', ls='--', lw=1.0,
						label=f'wildtype ({wt_fmt.format(metric[wt])})')
				ax.set_xlabel('kcat multiplier')
				ax.set_ylabel(ylabel)
				ax.set_title(f'{ESTIMATOR_SHORT[ESTIMATORS[e]]}')
				ax.invert_xaxis()

		_facet(0, completion, 'Completion rate', '{:.0%}')
		axes[0][0].set_ylim(0, 1.05)
		axes[0][0].legend(fontsize=6, ncol=2, loc='lower left')
		_facet(1, avg_dt, 'Avg doubling time (min)', '{:.1f} min')
		_facet(2, avg_mass, 'Avg dry mass (fg)', '{:.0f} fg')

		fig.suptitle(
			'new_gene_trl_eff_v2_kcat_multiplier_sweep: completion / doubling '
			'time / dry mass vs kcat multiplier (line per trl-eff)\n'
			f'(gens {ignore_first_n_gens}-{n_total_gens - 1}, all seeds)',
			fontsize=12)
		fig.tight_layout(rect=[0, 0, 1, 0.95])
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close('all')


if __name__ == '__main__':
	Plot().cli()
