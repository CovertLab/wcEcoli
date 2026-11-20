"""
How much each kcat-bounded enzyme's concentration changes with GFP.

# Figure 4B

Histogram over enzymes of concentration with GFP over without (showcase
design nearest SHOWCASE_PROTEOME_SHARE against the no-GFP control). An
enzyme's concentration is its window-mean catalyst count (FBAResults
catalyst_counts) over window-mean cell volume, averaged over surviving
lineages. Enzymes are the catalysts with a kcat estimate under
KCAT_ESTIMATOR, the ones a kcat bound can act on; the same set as 4C and 4G.
Bins of BIN_WIDTH from BIN_RANGE, with values outside clipped into the end
bins. Dashed line: no change; solid line: the median.

The supporting claim: GFP lowers metabolic enzyme concentrations broadly,
so a metabolic burden would be expected if enzymes were limiting.
"""

import pickle

import numpy as np
from matplotlib import pyplot as plt
from matplotlib.lines import Line2D

from models.ecoli.analysis import variantAnalysisPlot
from wholecell.analysis import gfp_paper_constants as C
from wholecell.analysis import gfp_paper_data
from wholecell.analysis import gfp_paper_style as style

PERCENT = 100.0
BIN_RANGE = (0.2, 1.3)
BIN_WIDTH = 0.05


class Plot(variantAnalysisPlot.VariantAnalysisPlot):
	def do_plot(self, inputDir, plotOutDir, plotOutFileName, simDataFile,
			validationDataFile, metadata):
		data = gfp_paper_data.load(inputDir, simDataFile, metadata,
			plotOutDir, self.cpus)
		with open(simDataFile, 'rb') as f:
			sim_data = pickle.load(f)
		showcase = gfp_paper_data.nearest_design(data.designs,
			C.SHOWCASE_PROTEOME_SHARE)
		design = int(showcase['variant'])
		ratios, enzymes = gfp_paper_data.enzyme_concentration_ratios(data,
			sim_data, [design], plotOutDir, self.cpus)
		values = ratios[design].values
		values = values[np.isfinite(values)]
		median = float(np.median(values)) if values.size else np.nan
		edges = np.arange(BIN_RANGE[0], BIN_RANGE[1] + BIN_WIDTH / 2,
			BIN_WIDTH)
		clipped = np.clip(values, edges[0] + 1e-9, edges[-1] - 1e-9)

		style.apply()
		fig, ax = plt.subplots(figsize=style.panel_size('single', 0.55))
		ax.hist(clipped, bins=edges, color=style.COLORS['enzyme'],
			edgecolor='white', linewidth=0.6, zorder=2)
		ax.axvline(1.0, color=style.INK_MUTED, linestyle=(0, (4, 3)),
			linewidth=0.8, zorder=3)
		ax.axvline(median, color=style.INK, linewidth=1.0, zorder=3)
		ax.set_xlim(*BIN_RANGE)
		ax.set_xlabel('enzyme concentration, GFP %.0f%% / no GFP'
			% (PERCENT * showcase['gfp_proteome_share']))
		ax.set_ylabel('enzymes')
		ax.legend(handles=[Line2D([], [], color=style.INK_MUTED,
			linestyle=(0, (4, 3)), label='no change'), Line2D([], [],
			color=style.INK, label='median, %.2f' % median)],
			loc='upper left', bbox_to_anchor=(1.0, 1.0))

		style.save_panel(fig, plotOutDir, style.panel_name(__file__),
			script=__file__, batch=inputDir, sim_data_path=simDataFile,
			n=int(values.size),
			notes='Showcase = variant %d. %d kcat-bounded enzymes (estimator'
				' %s) present without GFP; median ratio %.4f; %d below and %d'
				' above the plotted range, clipped into the end bins.'
				% (design, values.size, C.KCAT_ESTIMATOR, median,
				int((values < BIN_RANGE[0]).sum()),
				int((values > BIN_RANGE[1]).sum())),
			**data.provenance([gfp_paper_data.CONTROL_VARIANT, design]))


if __name__ == '__main__':
	Plot().cli()
