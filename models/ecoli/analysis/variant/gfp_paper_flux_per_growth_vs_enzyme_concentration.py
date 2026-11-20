"""
Metabolic flux per unit growth and enzyme concentration against GFP share,
relative to no GFP.

# Figure 4C

One point per design (window means over surviving lineages), hollow where
some lineages failed:
	enzyme concentration      median over the kcat-bounded enzymes (as 4B)
	                          of the design's concentration over the no-GFP
	                          control's
	flux per unit growth      median over base reactions of the design's
	                          flux over the control's (mmol per g dry weight
	                          per hour, so composition changes cancel), then
	                          divided by the design's relative growth rate:
	                          the flux spent per unit of biomass made.
	                          Reactions count only if they carry at least
	                          FLUX_FLOOR in both conditions in the same
	                          direction.
The no-GFP control sits at 1 on both. Dashed line: no change.

The supporting claim: enzyme concentration falls steadily with GFP while
metabolism per unit growth does not change.
"""

import pickle

import numpy as np
from matplotlib import pyplot as plt

from models.ecoli.analysis import variantAnalysisPlot
from wholecell.analysis import gfp_paper_constants as C
from wholecell.analysis import gfp_paper_data
from wholecell.analysis import gfp_paper_style as style

PERCENT = 100.0


class Plot(variantAnalysisPlot.VariantAnalysisPlot):
	def do_plot(self, inputDir, plotOutDir, plotOutFileName, simDataFile,
			validationDataFile, metadata):
		data = gfp_paper_data.load(inputDir, simDataFile, metadata,
			plotOutDir, self.cpus)
		with open(simDataFile, 'rb') as f:
			sim_data = pickle.load(f)
		control = gfp_paper_data.CONTROL_VARIANT
		designs = data.designs[(data.designs['promoter'] > 0)
			& (data.designs['n_lineages_survived'] > 0)
			& np.isfinite(data.designs['gfp_proteome_share'])]
		variants = [int(v) for v in designs['variant']]
		ratios, enzymes = gfp_paper_data.enzyme_concentration_ratios(data,
			sim_data, variants, plotOutDir, self.cpus)
		fluxes = gfp_paper_data.vector_profiles(data, 'base_fluxes',
			[control] + variants, plotOutDir, self.cpus,
			extra=(float(sim_data.constants.cell_density.asNumber()),))
		base = fluxes[control].mean

		enzyme_ratio, flux_per_growth, n_reactions = [], [], []
		for _, row in designs.iterrows():
			variant = int(row['variant'])
			enzyme_ratio.append(float(np.nanmedian(ratios[variant])))
			flux = fluxes[variant].mean
			ok = ((np.abs(base) >= C.FLUX_FLOOR)
				& (np.abs(flux) >= C.FLUX_FLOOR) & (np.sign(base)
				== np.sign(flux)))
			n_reactions.append(int(ok.sum()))
			flux_per_growth.append(float(np.median(flux[ok] / base[ok]))
				/ row['relative_growth_rate'] if ok.any() else np.nan)
		x = np.append(0.0, PERCENT * designs['gfp_proteome_share'].values)
		failed = np.append(False, designs['failed'].values)

		style.apply()
		fig, ax = plt.subplots(figsize=style.panel_size('single', 0.62))
		style.design_scatter(ax, x, np.append(1.0, flux_per_growth),
			style.INK_MUTED, failed=failed, label='metabolic flux per unit'
			' growth', zorder=3)
		style.design_scatter(ax, x, np.append(1.0, enzyme_ratio),
			style.COLORS['enzyme'], failed=failed, label='metabolic enzyme'
			' concentration', zorder=3)
		ax.axhline(1.0, color=style.INK_MUTED, linestyle=(0, (4, 3)),
			linewidth=0.8, zorder=1)
		ax.set_xlim(0, max(75.0, float(np.nanmax(x)) * 1.05))
		ax.set_ylim(0, max(1.2, float(np.nanmax(np.append(flux_per_growth,
			enzyme_ratio))) * 1.05))
		ax.set_xlabel('GFP share of proteome (%)')
		ax.set_ylabel('relative to no GFP')
		style.legend_with_failed(ax, bool(failed.any()), loc='lower left')

		style.save_panel(fig, plotOutDir, style.panel_name(__file__),
			script=__file__, batch=inputDir, sim_data_path=simDataFile,
			n=len(variants),
			notes='%d kcat-bounded enzymes (estimator %s); reactions per'
				' design passing the flux filter: %d-%d. Per design (variant:'
				' enzyme ratio, flux per growth): %s. Omitted (no survivors):'
				' %s.' % (len(enzymes), C.KCAT_ESTIMATOR, min(n_reactions)
				if n_reactions else 0, max(n_reactions) if n_reactions
				else 0, ', '.join('%d: %.3f, %.3f' % (v, e, f) for v, e, f
				in zip(variants, enzyme_ratio, flux_per_growth)),
				gfp_paper_data.lost_designs(data.designs[
				data.designs['promoter'] > 0])),
			**data.provenance([control] + variants))


if __name__ == '__main__':
	Plot().cli()
