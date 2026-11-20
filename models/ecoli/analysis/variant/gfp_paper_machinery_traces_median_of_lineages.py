"""
RNA polymerases and ribosomes per cell through the run, median of lineages.

# Figure 2E

Each generation of each lineage is split into TRACE_POINTS equal fractions of
the cell cycle and the mean total (active plus free) taken in each, so
lineages with cycles of different length line up generation by generation.
Lines are the median across surviving lineages at each point, for the no-GFP
control (grey) and the design nearest the second entry of
FIG5_PROTEOME_SHARES (RNAP and ribosome colours). The dashed line is
induction.

This is a variant analysis rather than multigen, because a median across
lineages needs every lineage at once.

The supporting claim: once GFP is induced, RNAP and ribosome counts fall
within a few generations and stay low.
"""

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
		high = gfp_paper_data.nearest_design(data.designs,
			C.FIG5_PROTEOME_SHARES[-1])
		control = gfp_paper_data.CONTROL_VARIANT
		variants = [control, int(high['variant'])]
		traces = gfp_paper_data.cycle_traces(data, variants, plotOutDir,
			self.cpus)
		alive = data.lineages[data.lineages['survived']]
		keep = set(zip(alive['variant'], alive['seed']))
		traces = traces[[k in keep for k in zip(traces['variant'],
			traces['seed'])]]
		median = traces.groupby(['variant', 'x'])[['rnap', 'ribosomes']] \
			.median().reset_index()
		n_lineages = traces.groupby('variant')['seed'].nunique().to_dict()

		style.apply()
		fig, axes = plt.subplots(2, 1, sharex=True,
			figsize=style.panel_size('single', 0.85))
		for ax, column, entity, label in (
				(axes[0], 'rnap', 'rnap', 'RNAP counts'),
				(axes[1], 'ribosomes', 'ribosome', 'ribosome counts')):
			for variant in variants:
				line = median[median['variant'] == variant]
				color = style.COLORS['control'] if variant == control else (
					style.COLORS[entity])
				ax.plot(line['x'], line[column], color=color, linewidth=1.0,
					zorder=2 if variant == control else 3)
			ax.axvline(C.NEW_GENE_INDUCTION_GEN, color=style.INK_MUTED,
				linestyle=(0, (4, 3)), linewidth=0.8, zorder=1)
			ax.set_ylim(bottom=0)
			ax.set_ylabel(label)
			ax.yaxis.set_major_formatter(plt.FuncFormatter(
				lambda v, _: '{:,.0f}'.format(v)))
		axes[0].annotate('induction', xy=(C.NEW_GENE_INDUCTION_GEN, 1),
			xycoords=('data', 'axes fraction'), xytext=(3, -2),
			textcoords='offset points', va='top',
			fontsize=style.FONT_SIZES['annotation'], color=style.INK_MUTED)
		axes[1].set_xlim(0, C.N_GENERATIONS)
		axes[1].set_xticks(range(0, C.N_GENERATIONS + 1, 4))
		axes[1].set_xlabel('generation')
		from matplotlib.lines import Line2D
		handles = [Line2D([], [], color=style.COLORS['control'],
			label='no GFP'), Line2D([], [], color=style.INK,
			label='GFP %.0f%% of proteome, median of %d lineages' % (
			PERCENT * high['gfp_proteome_share'],
			n_lineages.get(int(high['variant']), 0)))]
		style.legend_below(fig, handles, ncol=1)

		style.save_panel(fig, plotOutDir, style.panel_name(__file__),
			script=__file__, batch=inputDir, sim_data_path=simDataFile,
			n=int(sum(n_lineages.values())),
			notes='Median over surviving lineages of %d-point within-cycle'
				' means of total RNAP and ribosomes. Lineages: %s. High GFP ='
				' variant %d.' % (gfp_paper_data.TRACE_POINTS, n_lineages,
				high['variant']),
			**data.provenance(variants))


if __name__ == '__main__':
	Plot().cli()
