"""
Mass-category change for the kcat_v2_estimator_sweep (wildtype, estimators x
multiplier grid, no GFP), vs kcat multiplier, one line per estimator.

The WT-multiplier counterpart of the GFP mass analyses.  Here tightening the
multiplier genuinely makes metabolism bind, so this is the regime to ask whether
metabolite pools actually *deplete* -- i.e. whether the small-molecule
CONCENTRATION (mass / cell volume) drops, in contrast to the GFP case where
concentration was preserved and only the absolute mass fell (a cell-size effect).

Per dry-mass category (protein, RNA, DNA, small molecule): fraction of dry mass,
absolute mass, and concentration (mass / cell volume) vs multiplier, one line
per estimator (gap / smoothed_max / drop_top_20).  Reads recorded listeners; no
re-simulation.
"""

import numpy as np
from matplotlib import pyplot as plt

from models.ecoli.analysis import variantAnalysisPlot
from models.ecoli.sim.variants.kcat_v2_estimator_sweep import (
	ESTIMATOR_LABELS, ESTIMATOR_SHORT_NAMES, MULTIPLIERS, is_wildtype)
from models.ecoli.analysis.variant.kcat_v2_estimator_sweep_growth_limitation import (
	_series_by_estimator, ESTIMATOR_BASE_COLORS, WILDTYPE_COLOR)
from models.ecoli.analysis.variant.new_gene_trl_eff_v2_estimator_sweep_enzyme_mass import (
	_col)
from wholecell.analysis.analysis_tools import exportFigure


MASS_CATEGORIES = [
	('protein', 'proteinMass'), ('RNA', 'rnaMass'),
	('DNA', 'dnaMass'), ('small mol', 'smallMoleculeMass')]


class Plot(variantAnalysisPlot.VariantAnalysisPlot):
	def do_plot(self, inputDir, plotOutDir, plotOutFileName, simDataFile,
				validationDataFile, metadata):
		variant_indexes = self.ap.get_variants()
		n_total_gens = self.ap.n_generation
		ignore_first_n_gens = max(n_total_gens - 8, 0)
		window = np.arange(ignore_first_n_gens, n_total_gens)

		absolute = {name: {} for name, _ in MASS_CATEGORIES}
		frac = {name: {} for name, _ in MASS_CATEGORIES}
		conc_density = {name: {} for name, _ in MASS_CATEGORIES}
		for vi in variant_indexes:
			cells = self.ap.get_cells(
				variant=[vi], generation=window, only_successful=True)
			dcol = _col(cells, 'Mass', 'dryMass')
			vcol = _col(cells, 'Mass', 'cellVolume')
			dm = (float(np.nanmean(dcol))
				if dcol is not None and np.nanmean(dcol) > 0 else np.nan)
			vol = (float(np.nanmean(vcol))
				if vcol is not None and np.nanmean(vcol) > 0 else np.nan)
			for name, col in MASS_CATEGORIES:
				arr = _col(cells, 'Mass', col)
				a = float(np.nanmean(arr)) if arr is not None else np.nan
				absolute[name][vi] = a
				frac[name][vi] = a / dm if np.isfinite(dm) else np.nan
				conc_density[name][vi] = a / vol if np.isfinite(vol) else np.nan

		wt = next((vi for vi in variant_indexes if is_wildtype(vi)), None)
		x = np.array(MULTIPLIERS)

		def _facet(ax, metric, ylabel, title, wt_fmt=None):
			s = _series_by_estimator(metric, variant_indexes)
			for label in ESTIMATOR_LABELS:
				ax.plot(x, s[label], marker='o', lw=1.4,
					color=ESTIMATOR_BASE_COLORS[label],
					label=ESTIMATOR_SHORT_NAMES[ESTIMATOR_LABELS.index(label)])
			if wt is not None and np.isfinite(metric.get(wt, np.nan)):
				lbl = ('wildtype' if wt_fmt is None
					else f'wildtype ({wt_fmt.format(metric[wt])})')
				ax.axhline(metric[wt], color=WILDTYPE_COLOR, ls='--', lw=1.0,
					label=lbl)
			ax.set_xlabel('kcat multiplier')
			ax.set_ylabel(ylabel)
			ax.set_title(title)
			ax.invert_xaxis()  # tighter bound (smaller multiplier) on the right

		fig, axes = plt.subplots(len(MASS_CATEGORIES), 3, figsize=(18, 16))
		for r, (name, _col_name) in enumerate(MASS_CATEGORIES):
			_facet(axes[r, 0], frac[name],
				'fraction of dry mass', f'{name} mass fraction')
			_facet(axes[r, 1], absolute[name],
				'mass (fg)', f'{name} mass (absolute)')
			_facet(axes[r, 2], conc_density[name],
				'mass / volume (fg/fL)', f'{name} concentration')
		axes[0][0].legend(fontsize=7, loc='best')

		fig.suptitle(
			'kcat_v2_estimator_sweep (wildtype): mass categories vs kcat '
			'multiplier -- fraction, absolute, concentration (line per estimator)'
			'.\nConcentration dropping as the bound tightens = genuine '
			'metabolite depletion (vs the GFP size-effect, where it held flat).',
			fontsize=12)
		fig.tight_layout(rect=[0, 0, 1, 0.96])
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close('all')


if __name__ == '__main__':
	Plot().cli()
