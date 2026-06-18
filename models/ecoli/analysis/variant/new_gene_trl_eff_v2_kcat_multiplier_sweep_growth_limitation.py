"""
Growth-limitation decomposition for the new_gene_trl_eff_v2_kcat_multiplier_sweep,
vs kcat multiplier, one line per GFP translation-efficiency level, faceted by
estimator (gap | smoothed_max).

Rows = growth rate mu, active ribosomes, active RNAP, ribosome elongation rate,
tRNA charged fraction, ppGpp; columns = estimator.  As the multiplier tightens
(right), a metabolism-limited regime shows the stringent signature: ppGpp up,
elongation/charging down.  Reads only recorded listeners; no re-simulation.
"""

import os

import numpy as np
from matplotlib import pyplot as plt

from models.ecoli.analysis import variantAnalysisPlot
from models.ecoli.sim.variants.new_gene_trl_eff_v2_kcat_multiplier_sweep import (
	ESTIMATORS, ESTIMATOR_SHORT, TRL_EFF_VALUES, MULTIPLIERS,
	is_wildtype, variant_to_cell)
from models.ecoli.analysis.variant.new_gene_trl_eff_v2_estimator_sweep_growth_limitation import (
	_safe_mean)
from wholecell.analysis.analysis_tools import exportFigure
from wholecell.io.tablereader import TableReader


TRL_COLORS = plt.cm.viridis(np.linspace(0, 0.9, len(TRL_EFF_VALUES)))
SEC_PER_HOUR = 3600.


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
			mu[vi] = _safe_mean(cells, 'Mass', 'instantaneous_growth_rate',
				scale=SEC_PER_HOUR)
			ribo[vi] = _safe_mean(cells, 'UniqueMoleculeCounts',
				'uniqueMoleculeCounts', idx=ribo_idx)
			rnap[vi] = _safe_mean(cells, 'UniqueMoleculeCounts',
				'uniqueMoleculeCounts', idx=rnap_idx)
			elong[vi] = _safe_mean(cells, 'RibosomeData',
				'effectiveElongationRate')
			charged[vi] = _safe_mean(cells, 'GrowthLimits',
				'fraction_trna_charged')
			ppgpp[vi] = _safe_mean(cells, 'GrowthLimits', 'ppgpp_conc')

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
		metrics = [
			(mu, 'Growth rate (1/h)', '{:.3f}'),
			(ribo, 'Active ribosomes', None),
			(rnap, 'Active RNAP', None),
			(elong, 'Elongation rate (aa/s)', None),
			(charged, 'tRNA charged fraction', None),
			(ppgpp, 'ppGpp (uM)', None),
		]
		fig, axes = plt.subplots(len(metrics), n_est,
			figsize=(6 * n_est, 3 * len(metrics)), squeeze=False)

		for row, (metric, ylabel, wt_fmt) in enumerate(metrics):
			s = _series(metric)
			for e in range(n_est):
				ax = axes[row][e]
				for t in range(len(TRL_EFF_VALUES)):
					ax.plot(x, s.get((e, t), [np.nan] * len(MULTIPLIERS)),
						marker='o', lw=1.3, color=TRL_COLORS[t],
						label=f'trl {TRL_EFF_VALUES[t]}')
				if wt is not None and np.isfinite(metric.get(wt, np.nan)):
					lbl = ('wildtype' if wt_fmt is None
						else f'wildtype ({wt_fmt.format(metric[wt])})')
					ax.axhline(metric[wt], color='k', ls='--', lw=1.0, label=lbl)
				ax.set_xlabel('kcat multiplier')
				ax.set_ylabel(ylabel)
				ax.set_title(ESTIMATOR_SHORT[ESTIMATORS[e]])
				ax.invert_xaxis()
		axes[0][0].legend(fontsize=6, ncol=2, loc='lower left')

		fig.suptitle(
			'new_gene_trl_eff_v2_kcat_multiplier_sweep: growth-limitation '
			'decomposition vs kcat multiplier (line per trl-eff)\n'
			'tightening (right) with ppGpp up / elongation+charging down = '
			f'metabolism-limited.  gens {ignore_first_n_gens}-{n_total_gens-1}',
			fontsize=12)
		fig.tight_layout(rect=[0, 0, 1, 0.96])
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close('all')


if __name__ == '__main__':
	Plot().cli()
