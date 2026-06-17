"""
Mass categories resolved by block (gfp_only / gap / smoothed_max) for the
new_gene_trl_eff_v2_estimator_sweep.

Companion to ..._enzyme_mass.py, which showed mass categories for the gfp_only
block only.  Here every mass quantity is compared across the three kcat regimes
so we can see whether imposing the kcat bounds actually changes metabolism --
in particular the small-molecule (metabolite) mass -- even where growth does not
differ.

Two PDFs (both averaged over the settled window, last 8 gens, all seeds):

  <name>.pdf            -- per category (protein/RNA/DNA/small molecule), the
                          mass FRACTION and ABSOLUTE mass vs trl-eff, one line
                          per block (+ wildtype).  The small-molecule panels
                          answer "do the constraints change metabolite mass?".

  <name>_scenarios.pdf  -- one row per scenario (gfp_only/gap/smoothed_max),
                          columns = dry-mass category fractions, RNA
                          sub-fractions (rRNA/tRNA/mRNA), and bound-vs-demand
                          (catalyst conc = the kcat ceiling, and reaction flux,
                          each normalized to wildtype).

Reuses the verified read + pair-mapping helpers from ..._enzyme_mass.py.  Reads
only recorded listeners -- no re-simulation.
"""

import pickle

import numpy as np
from matplotlib import pyplot as plt

from models.ecoli.analysis import variantAnalysisPlot
from models.ecoli.sim.variants.new_gene_trl_eff_v2_estimator_sweep import (
	BLOCK_ESTIMATORS,
	BLOCK_SHORT_NAMES,
	TRL_EFF_VALUES,
	is_wildtype,
)
from models.ecoli.analysis.variant.new_gene_trl_eff_v2_estimator_sweep_enzyme_mass import (
	_build_pair_indices,
	_col,
	_series_by_block,
	CATEGORY_COLORS,
)
from models.ecoli.processes.metabolism import (
	COUNTS_UNITS, VOLUME_UNITS, TIME_UNITS, MASS_UNITS)
from wholecell.analysis.analysis_tools import exportFigure
from wholecell.utils import units


BLOCK_COLORS = ('#7f7f7f', '#1f77b4', '#2ca02c')  # gfp_only, gap, smoothed_max
WILDTYPE_COLOR = '#222222'

MASS_CATEGORIES = [
	('protein', 'proteinMass'), ('RNA', 'rnaMass'),
	('DNA', 'dnaMass'), ('small mol', 'smallMoleculeMass')]
RNA_SUB = [('rRNA', 'rRnaMass'), ('tRNA', 'tRnaMass'), ('mRNA', 'mRnaMass')]


class Plot(variantAnalysisPlot.VariantAnalysisPlot):
	def do_plot(self, inputDir, plotOutDir, plotOutFileName, simDataFile,
				validationDataFile, metadata):
		variant_indexes = self.ap.get_variants()
		n_total_gens = self.ap.n_generation
		ignore_first_n_gens = max(n_total_gens - 8, 0)
		window = np.arange(ignore_first_n_gens, n_total_gens)
		x = np.array(TRL_EFF_VALUES)

		with open(simDataFile, 'rb') as f:
			sim_data = pickle.load(f)
		all_cells = self.ap.get_cells(only_successful=True)
		if len(all_cells) == 0:
			print('mass_by_block: no successful cells; skipping.')
			return
		rxn_indexes, cat_indexes = _build_pair_indices(sim_data, all_cells[0])
		cell_density_value = sim_data.constants.cell_density.asNumber(
			MASS_UNITS / VOLUME_UNITS)

		all_cats = MASS_CATEGORIES + RNA_SUB
		absolute = {name: {} for name, _ in all_cats}
		frac = {name: {} for name, _ in all_cats}
		dry, conc, flux = {}, {}, {}
		for vi in variant_indexes:
			cells = self.ap.get_cells(
				variant=[vi], generation=window, only_successful=True)
			dcol = _col(cells, 'Mass', 'dryMass')
			dry[vi] = (float(np.nanmean(dcol))
				if dcol is not None and np.nanmean(dcol) > 0 else np.nan)
			for name, col in all_cats:
				arr = _col(cells, 'Mass', col)
				absolute[name][vi] = (float(np.nanmean(arr))
					if arr is not None else np.nan)
				frac[name][vi] = (absolute[name][vi] / dry[vi]
					if np.isfinite(dry[vi]) else np.nan)

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

		def _norm(metric):
			base = metric.get(wt) if wt is not None else None
			if base is None or not np.isfinite(base) or base == 0:
				return {vi: np.nan for vi in metric}
			return {vi: metric[vi] / base for vi in metric}

		conc_norm = _norm(conc)
		flux_norm = _norm(flux)

		def _blocks(ax, metric, title, ylabel, wt_fmt=None):
			series = _series_by_block(metric, variant_indexes)
			for b in range(len(BLOCK_ESTIMATORS)):
				ax.plot(x, series[b], marker='o', lw=1.5,
					color=BLOCK_COLORS[b], label=BLOCK_SHORT_NAMES[b])
			if wt is not None and np.isfinite(metric.get(wt, np.nan)):
				lbl = ('wildtype' if wt_fmt is None
					else f'wildtype ({wt_fmt.format(metric[wt])})')
				ax.axhline(metric[wt], color=WILDTYPE_COLOR, ls='--', lw=1.0,
					label=lbl)
			ax.set_xlabel('GFP translation efficiency')
			ax.set_ylabel(ylabel)
			ax.set_title(title)
			ax.legend(fontsize=7)

		# ---- Figure 1: per category, fraction + absolute, lines per block ----
		fig, axes = plt.subplots(len(MASS_CATEGORIES), 2, figsize=(12, 16))
		for r, (name, _col_name) in enumerate(MASS_CATEGORIES):
			_blocks(axes[r, 0], frac[name],
				f'{name} mass fraction', 'fraction of dry mass')
			_blocks(axes[r, 1], absolute[name],
				f'{name} mass (absolute)', 'mass (fg)')
		fig.suptitle(
			'new_gene_trl_eff_v2_estimator_sweep: mass categories by block '
			'(do the kcat constraints change metabolite / small-molecule mass?)',
			fontsize=13)
		fig.tight_layout(rect=[0, 0, 1, 0.97])
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close('all')

		# ---- Figure 2: one row per scenario ----
		fig2, axes2 = plt.subplots(len(BLOCK_ESTIMATORS), 3, figsize=(16, 11))
		frac_series = {name: _series_by_block(frac[name], variant_indexes)
			for name, _ in all_cats}
		conc_n_series = _series_by_block(conc_norm, variant_indexes)
		flux_n_series = _series_by_block(flux_norm, variant_indexes)
		for b in range(len(BLOCK_ESTIMATORS)):
			# Col 0: dry-mass category fractions
			ax = axes2[b, 0]
			for name, _ in MASS_CATEGORIES:
				ax.plot(x, frac_series[name][b], marker='o', lw=1.3,
					color=CATEGORY_COLORS[name], label=name)
			ax.set_ylabel('fraction of dry mass')
			ax.set_title(f'{BLOCK_SHORT_NAMES[b]}: dry-mass category fractions')
			ax.legend(fontsize=7)

			# Col 1: RNA sub-fractions
			ax = axes2[b, 1]
			for name, _ in RNA_SUB:
				ax.plot(x, frac_series[name][b], marker='o', lw=1.3,
					color=CATEGORY_COLORS[name], label=name)
			ax.set_ylabel('fraction of dry mass')
			ax.set_title(f'{BLOCK_SHORT_NAMES[b]}: RNA sub-fractions')
			ax.legend(fontsize=7)

			# Col 2: bound (catalyst conc) vs demand (flux), normalized to WT
			ax = axes2[b, 2]
			ax.plot(x, conc_n_series[b], marker='o', lw=1.5, color='#1f77b4',
				label='catalyst conc (= bound)')
			ax.plot(x, flux_n_series[b], marker='s', lw=1.5, color='#d62728',
				label='reaction flux (= demand)')
			ax.set_ylabel('fraction of wildtype')
			ax.set_title(f'{BLOCK_SHORT_NAMES[b]}: bound vs demand')
			ax.legend(fontsize=7)

			for c in range(3):
				axes2[b, c].set_xlabel('GFP translation efficiency')

		fig2.suptitle(
			'new_gene_trl_eff_v2_estimator_sweep: mass / bound-vs-demand per '
			'scenario (rows = gfp_only / gap / smoothed_max)',
			fontsize=13)
		fig2.tight_layout(rect=[0, 0, 1, 0.96])
		exportFigure(plt, plotOutDir, plotOutFileName + '_scenarios', metadata)
		plt.close('all')


if __name__ == '__main__':
	Plot().cli()
