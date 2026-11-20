"""
Metabolic enzyme supply capacity against translational demand for amino
acids, relative to no GFP.

# Figure 4G

One point per design against GFP share, window means over surviving
lineages, hollow where some lineages failed; the no-GFP control sits at 1:
	enzyme supply capacity   median over the kcat-bounded enzymes of
	                         concentration over the control's (as 4B, 4C);
	                         a kcat bound scales with concentration, so this
	                         is how much the ceiling on supply falls
	translational demand     amino acids built into protein per second per
	                         fL (active ribosomes x elongation rate / cell
	                         volume) over the control's
Both are concentration rates, so the gap between them compares like with
like. Lines join designs in order of GFP share to show the gap; they are a
guide, not a fitted relation. The shading is the gap.

The supporting claim: demand falls further than supply, so the enzymes the
cell keeps are more than enough for the protein it is making.
"""

import pickle

import numpy as np
from matplotlib import pyplot as plt

from models.ecoli.analysis import variantAnalysisPlot
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
			& np.isfinite(data.designs['gfp_proteome_share'])].sort_values(
			'gfp_proteome_share')
		variants = [int(v) for v in designs['variant']]
		ratios, enzymes = gfp_paper_data.enzyme_concentration_ratios(data,
			sim_data, variants, plotOutDir, self.cpus)
		base_demand = float(data.designs.loc[data.designs['variant']
			== control, 'translation_demand'].iloc[0])
		supply = np.append(1.0, [float(np.nanmedian(ratios[v]))
			for v in variants])
		demand = np.append(1.0, designs['translation_demand'].values
			/ base_demand)
		x = np.append(0.0, PERCENT * designs['gfp_proteome_share'].values)
		failed = np.append(False, designs['failed'].values)

		style.apply()
		fig, ax = plt.subplots(figsize=style.panel_size('single', 0.62))
		ax.fill_between(x, demand, supply, where=supply >= demand,
			color=style.COLORS['enzyme'], alpha=0.15, linewidth=0, zorder=1)
		for y, color, label in ((supply, style.COLORS['enzyme'],
				'enzyme supply capacity'), (demand,
				style.COLORS['ribosome'], 'translational demand')):
			ax.plot(x, y, color=color, linewidth=1.0, zorder=2)
			style.design_scatter(ax, x, y, color, failed=failed, label=label,
				zorder=3)
		ax.set_xlim(0, max(45.0, float(np.nanmax(x)) * 1.05))
		ax.set_ylim(0, 1.15)
		ax.set_xlabel('GFP share of proteome (%)')
		ax.set_ylabel('relative to no GFP')
		style.legend_with_failed(ax, bool(failed.any()), loc='lower left')

		style.save_panel(fig, plotOutDir, style.panel_name(__file__),
			script=__file__, batch=inputDir, sim_data_path=simDataFile,
			n=len(variants),
			notes='%d kcat-bounded enzymes. Per design (variant: GFP share,'
				' supply, demand): %s.' % (len(enzymes), ', '.join(
				'%d: %.3f, %.3f, %.3f' % (v, s, sp, d) for v, s, sp, d in
				zip(variants, designs['gfp_proteome_share'], supply[1:],
				demand[1:]))),
			**data.provenance([control] + variants))


if __name__ == '__main__':
	Plot().cli()
