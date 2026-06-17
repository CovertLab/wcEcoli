"""
Growth-limitation diagnosis for the kcat_v2_estimator_sweep (wildtype, gap /
smoothed_max / drop_top_20 estimators x a multiplier grid).

This is the wildtype-multiplier counterpart of
new_gene_trl_eff_v2_estimator_sweep_growth_limitation.py.  There, GFP ribosome
competition drove the burden; here there is no GFP -- the kcat bound is tightened
directly by the multiplier, so this is the regime where metabolism *should*
eventually become limiting.  The same decomposition reveals the signature:

  * If, as the multiplier tightens, mu falls while tRNA charging drops, ppGpp
    rises, and elongation rate falls -> the kcat cap is throttling supply
    (metabolism is the binding constraint).  Contrast the GFP sweep, where those
    supply signals stayed flat and mu tracked ribosome count (machinery loss).

Across the multiplier grid (one line per estimator), averaged over the settled
window (last 8 gens, all seeds), it plots (2x4):
  mu, active ribosomes, active RNAP, ribosome elongation rate,
  tRNA charged fraction, ppGpp, mu-vs-ribosomes scatter, growth per ribosome.

Reads only recorded listeners -- no re-simulation.
"""

import os

import numpy as np
from matplotlib import pyplot as plt

from models.ecoli.analysis import variantAnalysisPlot
from models.ecoli.sim.variants.kcat_v2_estimator_sweep import (
	ESTIMATOR_LABELS,
	ESTIMATOR_SHORT_NAMES,
	MULTIPLIERS,
	N_PER_ESTIMATOR,
	is_wildtype,
	variant_to_estimator,
)
from wholecell.analysis.analysis_tools import exportFigure, read_stacked_columns
from wholecell.io.tablereader import TableReader


ESTIMATOR_BASE_COLORS = {
	'v2_gap':          '#1f77b4',
	'v2_smoothed_max': '#2ca02c',
	'v2_drop_top_20':  '#ff7f0e',
}
WILDTYPE_COLOR = '#222222'
SEC_PER_HOUR = 3600.


def _series_by_estimator(metric, variant_indexes):
	"""Return {estimator_label: list-of-value-by-multiplier-index}."""
	by_est = {label: [np.nan] * N_PER_ESTIMATOR for label in ESTIMATOR_LABELS}
	for vi in variant_indexes:
		if is_wildtype(vi):
			continue
		est_idx, mult_idx = variant_to_estimator(vi)
		by_est[ESTIMATOR_LABELS[est_idx]][mult_idx] = metric.get(vi, np.nan)
	return by_est


def _safe_mean(cells, table, column, idx=None, scale=1.0):
	"""Stacked-column mean over the window, or NaN if unreadable/empty."""
	if len(cells) == 0:
		return np.nan
	try:
		data = read_stacked_columns(
			cells, table, column, remove_first=True, ignore_exception=True)
		if idx is not None:
			data = data[:, idx]
		return float(np.nanmean(data)) * scale
	except Exception:
		return np.nan


class Plot(variantAnalysisPlot.VariantAnalysisPlot):
	def do_plot(self, inputDir, plotOutDir, plotOutFileName, simDataFile,
				validationDataFile, metadata):
		variant_indexes = self.ap.get_variants()
		n_total_gens = self.ap.n_generation
		ignore_first_n_gens = max(n_total_gens - 8, 0)
		window = np.arange(ignore_first_n_gens, n_total_gens)

		all_cells = self.ap.get_cells(only_successful=True)
		if len(all_cells) == 0:
			print('growth_limitation: no successful cells; skipping.')
			return
		unique_ids = TableReader(
			os.path.join(all_cells[0], 'simOut', 'UniqueMoleculeCounts')
			).readAttribute('uniqueMoleculeIds')
		ribo_idx = unique_ids.index('active_ribosome')
		rnap_idx = unique_ids.index('active_RNAP')

		mu, ribo, rnap, elong, charged, ppgpp = ({} for _ in range(6))
		for vi in variant_indexes:
			cells = self.ap.get_cells(
				variant=[vi], generation=window, only_successful=True)
			mu[vi] = _safe_mean(
				cells, 'Mass', 'instantaneous_growth_rate', scale=SEC_PER_HOUR)
			ribo[vi] = _safe_mean(
				cells, 'UniqueMoleculeCounts', 'uniqueMoleculeCounts',
				idx=ribo_idx)
			rnap[vi] = _safe_mean(
				cells, 'UniqueMoleculeCounts', 'uniqueMoleculeCounts',
				idx=rnap_idx)
			elong[vi] = _safe_mean(
				cells, 'RibosomeData', 'effectiveElongationRate')
			charged[vi] = _safe_mean(
				cells, 'GrowthLimits', 'fraction_trna_charged')
			ppgpp[vi] = _safe_mean(cells, 'GrowthLimits', 'ppgpp_conc')

		mu_per_ribo = {
			vi: (mu[vi] / ribo[vi] * 1000.0
				if (np.isfinite(mu[vi]) and np.isfinite(ribo[vi])
					and ribo[vi] > 0) else np.nan)
			for vi in variant_indexes}

		wt = next((vi for vi in variant_indexes if is_wildtype(vi)), None)
		x = np.array(MULTIPLIERS)

		fig, axes = plt.subplots(2, 4, figsize=(20, 9))

		def _line(ax, metric, ylabel, title, wt_fmt=None):
			series = _series_by_estimator(metric, variant_indexes)
			for label in ESTIMATOR_LABELS:
				ax.plot(x, series[label], marker='o', lw=1.5,
					color=ESTIMATOR_BASE_COLORS[label],
					label=ESTIMATOR_SHORT_NAMES[ESTIMATOR_LABELS.index(label)])
			if wt is not None and np.isfinite(metric.get(wt, np.nan)):
				lbl = ('wildtype' if wt_fmt is None
					else f'wildtype ({wt_fmt.format(metric[wt])})')
				ax.axhline(metric[wt], color=WILDTYPE_COLOR, ls='--', lw=1.0,
					label=lbl)
			ax.set_xlabel('Multiplier')
			ax.set_ylabel(ylabel)
			ax.set_title(title)
			ax.invert_xaxis()  # tighter bound (smaller multiplier) on the right

		_line(axes[0, 0], mu, 'Growth rate (1/h)', 'Growth rate mu', '{:.3f}')
		axes[0, 0].legend(fontsize=7, loc='lower left')
		_line(axes[0, 1], ribo, 'Active ribosomes', 'Active ribosomes')
		_line(axes[0, 2], rnap, 'Active RNAP', 'Active RNAP')
		_line(axes[0, 3], elong, 'Elongation rate (aa/s)',
			'Ribosome elongation rate\n(per-ribosome output)')
		_line(axes[1, 0], charged, 'Charged tRNA fraction',
			'tRNA charged fraction\n(AA supply adequacy)')
		_line(axes[1, 1], ppgpp, 'ppGpp (uM)', 'ppGpp\n(starvation signal)')

		# Panel 7: mu vs active ribosomes scatter
		ax = axes[1, 2]
		for vi in variant_indexes:
			if is_wildtype(vi):
				ax.plot(ribo[vi], mu[vi], '*', ms=13, color=WILDTYPE_COLOR)
				continue
			est_idx, _m = variant_to_estimator(vi)
			ax.plot(ribo[vi], mu[vi], 'o', ms=5,
				color=ESTIMATOR_BASE_COLORS[ESTIMATOR_LABELS[est_idx]])
		ax.set_xlabel('Active ribosomes')
		ax.set_ylabel('Growth rate (1/h)')
		ax.set_title('Growth rate vs active ribosomes')
		handles = [plt.Line2D([], [], marker='o', ls='',
			color=ESTIMATOR_BASE_COLORS[lbl], label=ESTIMATOR_SHORT_NAMES[i])
			for i, lbl in enumerate(ESTIMATOR_LABELS)]
		handles.append(plt.Line2D([], [], marker='*', ls='', ms=11,
			color=WILDTYPE_COLOR, label='wildtype'))
		ax.legend(handles=handles, fontsize=7, loc='lower right')

		# Panel 8: growth per 1000 ribosomes
		_line(axes[1, 3], mu_per_ribo, 'Growth per 1000 ribosomes (1/h)',
			'Growth per ribosome', '{:.4f}')

		fig.suptitle(
			'kcat_v2_estimator_sweep (wildtype): does tightening the kcat bound '
			'limit growth via metabolic supply?\n'
			'(supply-limited signature = charging down / ppGpp up / elongation '
			f'down; averaged over gens {ignore_first_n_gens}-{n_total_gens - 1}, '
			'all seeds)',
			fontsize=13)
		fig.tight_layout(rect=[0, 0, 1, 0.93])
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close('all')


if __name__ == '__main__':
	Plot().cli()
