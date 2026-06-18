"""
Cell-mass / enzyme evidence for the new_gene_trl_eff_v2_kcat_multiplier_sweep,
vs kcat multiplier, one line per GFP translation-efficiency level, faceted by
estimator (gap | smoothed_max).

Rows = mean metabolic catalyst conc / wildtype (= the bound ceiling), mean flux
of constrained reactions, small-molecule mass fraction, small-molecule mass
(absolute); columns = estimator.  Shows whether tightening the bound (right)
reroutes metabolism / changes the metabolite pool.  Reads recorded listeners;
no re-simulation.
"""

import os
import pickle

import numpy as np
from matplotlib import pyplot as plt

from models.ecoli.analysis import variantAnalysisPlot
from models.ecoli.sim.variants.new_gene_trl_eff_v2_kcat_multiplier_sweep import (
	ESTIMATORS, ESTIMATOR_SHORT, TRL_EFF_VALUES, MULTIPLIERS,
	is_wildtype, variant_to_cell)
from models.ecoli.analysis.variant.new_gene_trl_eff_v2_estimator_sweep_enzyme_mass import (
	_build_pair_indices, _col)
from models.ecoli.processes.metabolism import (
	COUNTS_UNITS, VOLUME_UNITS, TIME_UNITS, MASS_UNITS)
from wholecell.analysis.analysis_tools import exportFigure
from wholecell.utils import units


TRL_COLORS = plt.cm.viridis(np.linspace(0, 0.9, len(TRL_EFF_VALUES)))


class Plot(variantAnalysisPlot.VariantAnalysisPlot):
	def do_plot(self, inputDir, plotOutDir, plotOutFileName, simDataFile,
				validationDataFile, metadata):
		variant_indexes = self.ap.get_variants()
		n_total_gens = self.ap.n_generation
		ignore_first_n_gens = max(n_total_gens - 8, 0)
		window = np.arange(ignore_first_n_gens, n_total_gens)

		with open(simDataFile, 'rb') as f:
			sim_data = pickle.load(f)
		all_cells = self.ap.get_cells(only_successful=True)
		if len(all_cells) == 0:
			print('mass: no successful cells; skipping.')
			return
		rxn_indexes, cat_indexes = _build_pair_indices(sim_data, all_cells[0])
		cell_density_value = sim_data.constants.cell_density.asNumber(
			MASS_UNITS / VOLUME_UNITS)

		conc, flux, smol_frac, smol_abs = {}, {}, {}, {}
		for vi in variant_indexes:
			cells = self.ap.get_cells(
				variant=[vi], generation=window, only_successful=True)
			dcol = _col(cells, 'Mass', 'dryMass')
			scol = _col(cells, 'Mass', 'smallMoleculeMass')
			if dcol is not None and np.nanmean(dcol) > 0:
				dm = float(np.nanmean(dcol))
				smol_abs[vi] = (float(np.nanmean(scol))
					if scol is not None else np.nan)
				smol_frac[vi] = (smol_abs[vi] / dm
					if np.isfinite(smol_abs[vi]) else np.nan)
			else:
				smol_abs[vi] = np.nan
				smol_frac[vi] = np.nan

			conc[vi], flux[vi] = np.nan, np.nan
			if len(cat_indexes) > 0:
				cc = _col(cells, 'FBAResults', 'catalyst_counts')
				c2m = _col(cells, 'EnzymeKinetics', 'countsToMolar')
				if cc is not None and c2m is not None:
					c2m1 = c2m.reshape(c2m.shape[0], -1)[:, :1]
					conc[vi] = float(np.nanmean(cc[:, cat_indexes] * c2m1))
				rf = _col(cells, 'FBAResults', 'reactionFluxes')
				cm = _col(cells, 'Mass', 'cellMass')
				if rf is not None and dcol is not None and cm is not None:
					conv = (dcol.reshape(-1) / cm.reshape(-1)
						* cell_density_value)
					with np.errstate(divide='ignore', invalid='ignore'):
						fv = ((COUNTS_UNITS / MASS_UNITS / TIME_UNITS)
							* (rf[:, rxn_indexes] / conv[:, np.newaxis])
							).asNumber(units.mmol / units.g / units.h)
					fv[~np.isfinite(fv)] = np.nan
					flux[vi] = float(np.nanmean(fv))

		wt = next((vi for vi in variant_indexes if is_wildtype(vi)), None)
		conc_wt = conc.get(wt) if wt is not None else None
		conc_norm = {vi: (conc[vi] / conc_wt
			if (conc_wt and np.isfinite(conc_wt) and conc_wt > 0) else np.nan)
			for vi in variant_indexes}

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
			(conc_norm, 'Catalyst conc / wildtype\n(= bound ceiling)', None),
			(flux, 'Flux of constrained rxns (mmol/g/h)', None),
			(smol_frac, 'Small-molecule mass fraction', None),
			(smol_abs, 'Small-molecule mass (fg)', None),
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
					ax.axhline(metric[wt], color='k', ls='--', lw=1.0,
						label='wildtype')
				ax.set_xlabel('kcat multiplier')
				ax.set_ylabel(ylabel)
				ax.set_title(ESTIMATOR_SHORT[ESTIMATORS[e]])
				ax.invert_xaxis()
		axes[0][0].legend(fontsize=6, ncol=2, loc='lower left')

		fig.suptitle(
			'new_gene_trl_eff_v2_kcat_multiplier_sweep: catalyst conc / flux / '
			'small-molecule mass vs kcat multiplier (line per trl-eff)\n'
			f'(gens {ignore_first_n_gens}-{n_total_gens - 1}, all seeds)',
			fontsize=12)
		fig.tight_layout(rect=[0, 0, 1, 0.96])
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close('all')


if __name__ == '__main__':
	Plot().cli()
