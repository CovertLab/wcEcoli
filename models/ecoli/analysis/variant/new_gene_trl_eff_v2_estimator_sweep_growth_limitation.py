"""
Growth-limitation diagnosis for the new_gene_trl_eff_v2_estimator_sweep.

Doubling time in this model is paced by bulk dry-mass accumulation (replication
initiates at a fixed mass-per-origin, division follows C+D later), so it is set
by the growth rate mu.  This analysis asks WHAT limits mu as GFP translation
efficiency rises: the amount of translation/transcription machinery
(ribosomes + RNAP), or the metabolic supply feeding them.

Across the trl-eff gradient (one line per block: gfp_only / gap / smoothed_max)
it plots, averaged over the settled window (last 8 gens, all seeds):

  Machinery & output           Supply / limitation signals
    1 growth rate mu (1/h)        5 tRNA charged fraction
    2 active ribosomes            6 ppGpp (uM)
    3 active RNAP                  7 mu vs active ribosomes (scatter)
    4 ribosome elongation (aa/s)  8 growth per 1000 ribosomes

Reading it:
  * If mu falls in lockstep with the active-ribosome count while elongation
    rate, tRNA charging, and ppGpp stay ~flat (and "growth per ribosome" is
    flat, scatter is linear through the origin) -> the burden is set by losing
    machinery (ribosome/RNAP), and metabolism is NOT limiting.
  * If elongation rate / charging fall and ppGpp rises as trl-eff increases ->
    the AA supply (metabolism) is throttling translation and contributes to
    the burden.

Uses only listeners already recorded, so no re-simulation is needed.
"""

import pickle

import numpy as np
from matplotlib import pyplot as plt

from models.ecoli.analysis import variantAnalysisPlot
from models.ecoli.sim.variants.new_gene_trl_eff_v2_estimator_sweep import (
	BLOCK_ESTIMATORS,
	BLOCK_SHORT_NAMES,
	TRL_EFF_VALUES,
	N_PER_BLOCK,
	is_wildtype,
	variant_to_block,
)
from wholecell.analysis.analysis_tools import exportFigure, read_stacked_columns
from wholecell.io.tablereader import TableReader

import os


BLOCK_COLORS = ('#7f7f7f', '#1f77b4', '#2ca02c')  # gfp_only, gap, smoothed_max
WILDTYPE_COLOR = '#222222'
SEC_PER_HOUR = 3600.


def _series_by_block(metric, variant_indexes):
	"""Return {block_idx: list-of-value-by-trl_eff-index}."""
	by_block = {b: [np.nan] * N_PER_BLOCK
		for b in range(len(BLOCK_ESTIMATORS))}
	for vi in variant_indexes:
		if is_wildtype(vi):
			continue
		block_idx, te_idx = variant_to_block(vi)
		by_block[block_idx][te_idx] = metric.get(vi, np.nan)
	return by_block


def _safe_mean(cells, table, column, idx=None, scale=1.0):
	"""Stacked-column mean over the window, or NaN if unreadable/empty."""
	if len(cells) == 0:
		return np.nan
	try:
		data = read_stacked_columns(
			cells, table, column, remove_first=True, ignore_exception=True)
		if idx is not None:
			data = data[:, idx]
		val = float(np.nanmean(data))
		return val * scale
	except Exception:
		return np.nan


class Plot(variantAnalysisPlot.VariantAnalysisPlot):
	def do_plot(self, inputDir, plotOutDir, plotOutFileName, simDataFile,
				validationDataFile, metadata):
		variant_indexes = self.ap.get_variants()
		n_total_gens = self.ap.n_generation
		ignore_first_n_gens = max(n_total_gens - 8, 0)
		window = np.arange(ignore_first_n_gens, n_total_gens)

		# Resolve active-ribosome / active-RNAP column indices from any cell.
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
				cells, 'Mass', 'instantaneous_growth_rate',
				scale=SEC_PER_HOUR)
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

		# growth per 1000 active ribosomes (flat -> ribosome-count-limited)
		mu_per_ribo = {
			vi: (mu[vi] / ribo[vi] * 1000.0
				if (np.isfinite(mu[vi]) and np.isfinite(ribo[vi])
					and ribo[vi] > 0) else np.nan)
			for vi in variant_indexes}

		wt = [vi for vi in variant_indexes if is_wildtype(vi)]
		wt = wt[0] if wt else None
		x = np.array(TRL_EFF_VALUES)

		fig, axes = plt.subplots(2, 4, figsize=(20, 9))

		def _line(ax, metric, ylabel, title, wt_fmt=None):
			series = _series_by_block(metric, variant_indexes)
			for b in range(len(BLOCK_ESTIMATORS)):
				ax.plot(x, series[b], marker='o', lw=1.5,
					color=BLOCK_COLORS[b], label=BLOCK_SHORT_NAMES[b])
			if wt is not None and np.isfinite(metric[wt]):
				lbl = ('wildtype' if wt_fmt is None
					else f'wildtype ({wt_fmt.format(metric[wt])})')
				ax.axhline(metric[wt], color=WILDTYPE_COLOR, ls='--', lw=1.0,
					label=lbl)
			ax.set_xlabel('GFP translation efficiency')
			ax.set_ylabel(ylabel)
			ax.set_title(title)

		_line(axes[0, 0], mu, 'Growth rate (1/h)',
			'Growth rate mu', '{:.3f}')
		axes[0, 0].legend(fontsize=7, loc='lower left')
		_line(axes[0, 1], ribo, 'Active ribosomes',
			'Active ribosomes')
		_line(axes[0, 2], rnap, 'Active RNAP', 'Active RNAP')
		_line(axes[0, 3], elong, 'Elongation rate (aa/s)',
			'Ribosome elongation rate\n(per-ribosome output)')
		_line(axes[1, 0], charged, 'Charged tRNA fraction',
			'tRNA charged fraction\n(AA supply adequacy)')
		_line(axes[1, 1], ppgpp, 'ppGpp (uM)',
			'ppGpp\n(starvation signal)')

		# Panel 7: mu vs active ribosomes scatter (linearity test)
		ax = axes[1, 2]
		for vi in variant_indexes:
			if is_wildtype(vi):
				ax.plot(ribo[vi], mu[vi], '*', ms=13, color=WILDTYPE_COLOR,
					label='wildtype')
				continue
			b, _te = variant_to_block(vi)
			ax.plot(ribo[vi], mu[vi], 'o', ms=5, color=BLOCK_COLORS[b])
		ax.set_xlabel('Active ribosomes')
		ax.set_ylabel('Growth rate (1/h)')
		ax.set_title('Growth rate vs active ribosomes\n'
			'(linear through origin -> ribosome-limited)')
		handles = [plt.Line2D([], [], marker='o', ls='', color=BLOCK_COLORS[b],
			label=BLOCK_SHORT_NAMES[b]) for b in range(len(BLOCK_ESTIMATORS))]
		handles.append(plt.Line2D([], [], marker='*', ls='', ms=11,
			color=WILDTYPE_COLOR, label='wildtype'))
		ax.legend(handles=handles, fontsize=7, loc='lower right')

		# Panel 8: growth per 1000 ribosomes (flat -> ribosome-count-limited)
		_line(axes[1, 3], mu_per_ribo,
			'Growth per 1000 ribosomes (1/h)',
			'Growth per ribosome\n(flat -> ribosome-count-limited)', '{:.4f}')

		fig.suptitle(
			'new_gene_trl_eff_v2_estimator_sweep: is the GFP doubling-time '
			'burden set by ribosome/RNAP loss or by metabolic supply?\n'
			f'(averaged over gens {ignore_first_n_gens}-{n_total_gens - 1}, '
			'all seeds)',
			fontsize=13)
		fig.tight_layout(rect=[0, 0, 1, 0.93])
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close('all')


if __name__ == '__main__':
	Plot().cli()
