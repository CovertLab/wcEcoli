"""
Ribosomes and RNA polymerases per cell against GFP's share of the proteome.

# Figure 2D

One point per design, window means over surviving lineages, hollow where some
lineages failed. Totals are active plus free (ribosomes: active plus the
smaller of free 30S and 50S; RNAP: active plus free complexes). Each axis
carries a small promoter x RBS heatmap of the same quantity, so a point can
be traced to its design.

Two versions, because counts fall partly because burdened cells are smaller:
	counts          molecules per cell (plain name)
	concentration   micromolar, count / (Avogadro x cell volume) (suffix
	                'concentration')

The supporting claim: more GFP means fewer ribosomes and RNAPs.
"""

import numpy as np
from matplotlib import pyplot as plt

from models.ecoli.analysis import variantAnalysisPlot
from wholecell.analysis import gfp_paper_data
from wholecell.analysis import gfp_paper_style as style

PERCENT = 100.0
VERSIONS = (
	(None, ('total_ribosomes', 'ribosome counts'),
		('total_rnap', 'RNAP counts'), '{:,.0f}'),
	('concentration', ('total_ribosomes_uM', 'ribosomes (µM)'),
		('total_rnap_uM', 'RNAP (µM)'), '{:,.1f}'),
)


class Plot(variantAnalysisPlot.VariantAnalysisPlot):
	def do_plot(self, inputDir, plotOutDir, plotOutFileName, simDataFile,
			validationDataFile, metadata):
		data = gfp_paper_data.load(inputDir, simDataFile, metadata,
			plotOutDir, self.cpus)
		designs = data.designs[np.isfinite(data.designs['gfp_proteome_share'])]
		x = PERCENT * designs['gfp_proteome_share'].values

		for part, ribosome, rnap, fmt in VERSIONS:
			style.apply()
			fig, axes = plt.subplots(1, 2, figsize=style.panel_size('double',
				0.38))
			notes = []
			for ax, (column, label), entity in zip(axes, (ribosome, rnap),
					('ribosome', 'rnap')):
				color = style.COLORS[entity]
				y = designs[column].values
				style.design_scatter(ax, x, y, color,
					failed=designs['failed'].values, zorder=3)
				ax.set_xlim(left=0)
				ax.set_ylim(0, float(np.nanmax(y)) * 1.1)
				ax.set_xlabel('GFP share of proteome (%)')
				ax.set_ylabel(label)
				ax.yaxis.set_major_formatter(plt.FuncFormatter(
					lambda v, _: fmt.format(v)))
				inset = ax.inset_axes([0.08, 0.08, 0.36, 0.36])
				induced = designs[designs['promoter'] > 0]
				style.design_heatmap(inset, designs, column, style.ramp_cmap(
					(style.GROWTH_RAMP[0], color)), float(np.nanmin(
					induced[column])), float(np.nanmax(induced[column])),
					annotate=False, colorbar=False)
				inset.set_xticks([])
				inset.set_yticks([])
				inset.set_xlabel('promoter →', fontsize=6, labelpad=1)
				inset.set_ylabel('RBS →', fontsize=6, labelpad=1)
				inset.set_title('per design', fontsize=6, pad=2)
				control = designs[designs['variant']
					== gfp_paper_data.CONTROL_VARIANT][column]
				notes.append('%s: control %s, lowest design %s' % (label,
					fmt.format(float(control.iloc[0])) if len(control)
					else 'n/a', fmt.format(float(np.nanmin(y)))))
			shown = designs[np.isfinite(designs[ribosome[0]])]
			if shown['failed'].any():
				fig.legend(handles=[style.failed_legend_handle()],
					loc='lower center', bbox_to_anchor=(0.5, -0.08))

			style.save_panel(fig, plotOutDir,
				style.panel_name(__file__, part), script=__file__,
				batch=inputDir, sim_data_path=simDataFile,
				n=int(len(designs)),
				notes='Window means over surviving lineages; totals are'
					' active + free. %s. Omitted (no survivors): %s.'
					% ('; '.join(notes), gfp_paper_data.lost_designs(
					data.designs)),
				**data.provenance(sorted(int(v) for v in designs['variant'])))


if __name__ == '__main__':
	Plot().cli()
