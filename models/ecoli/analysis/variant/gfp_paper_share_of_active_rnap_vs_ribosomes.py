"""
Share of the cell's active ribosomes and active RNA polymerases working on
gfp, against GFP's share of the proteome.

# Figure 3G

One pair of points per design: the share of active ribosomes making GFP
(RibosomeData n_ribosomes_per_transcript; ribosome colour) and the share of active RNAPs transcribing gfp (RNAP
colour), window means over surviving lineages. Hollow points are designs where
some lineages failed.

The short RNAP-coloured dashes mark, for each design, the RNAP share gfp would
hold if every gfp copy were packed with polymerases footprint to footprint
(gene length / RNAP footprint x promoter copies / active RNAP). It is drawn
per design because both copy number and the active RNAP pool move with
burden, so a single band would misplace the limit for the strongest designs.
An RNAP point on its dash is a design transcribing at the spacing cap.

The supporting claim: the ribosome share keeps rising with GFP while the RNAP
share is held at the spacing limit.
"""

import numpy as np
from matplotlib import pyplot as plt
from matplotlib.lines import Line2D

from models.ecoli.analysis import variantAnalysisPlot
from wholecell.analysis import gfp_paper_data
from wholecell.analysis import gfp_paper_style as style

PERCENT = 100.0


class Plot(variantAnalysisPlot.VariantAnalysisPlot):
	def do_plot(self, inputDir, plotOutDir, plotOutFileName, simDataFile,
			validationDataFile, metadata):
		data = gfp_paper_data.load(inputDir, simDataFile, metadata,
			plotOutDir, self.cpus)
		designs = data.designs[data.designs['promoter'] > 0]
		designs = designs[np.isfinite(designs['gfp_proteome_share'])]

		style.apply()
		fig, ax = plt.subplots(figsize=style.panel_size('single', 0.62))
		x = PERCENT * designs['gfp_proteome_share']
		ax.scatter(x, PERCENT * designs['gfp_rnap_share_at_spacing_limit'],
			marker='_', s=40, linewidths=1.0, color=style.COLORS['rnap'],
			alpha=0.55, zorder=1)
		style.design_scatter(ax, x, PERCENT * designs['gfp_ribosome_share'],
			style.COLORS['ribosome'], failed=designs['failed'],
			label='active ribosomes', zorder=3)
		style.design_scatter(ax, x, PERCENT * designs['gfp_rnap_share'],
			style.COLORS['rnap'], failed=designs['failed'],
			label='active RNAPs', zorder=3)

		ax.set_xlim(0, max(75, float(np.nanmax(x)) * 1.05 if len(x) else 75))
		top = np.nanmax(PERCENT * designs[['gfp_ribosome_share',
			'gfp_rnap_share']].values) if len(designs) else 75
		ax.set_ylim(0, max(75, float(top) * 1.05))
		ax.set_xlabel('GFP share of proteome (%)')
		ax.set_ylabel('share of active machines on gfp (%)')

		handles, labels = ax.get_legend_handles_labels()
		handles.append(Line2D([], [], linestyle='none', marker='_',
			markersize=7, markeredgewidth=1.0, color=style.COLORS['rnap'],
			alpha=0.55))
		labels.append('RNAP share at the spacing limit')
		if designs['failed'].any():
			handles.append(style.failed_legend_handle())
			labels.append(style.FAILED_LABEL)
		ax.legend(handles, labels, loc='upper left')

		variants = sorted(int(v) for v in designs['variant'])
		style.save_panel(fig, plotOutDir, style.panel_name(__file__),
			script=__file__, batch=inputDir, sim_data_path=simDataFile,
			n=len(variants),
			notes='Shares are window means over surviving lineages: ribosomes'
				' making GFP / active ribosomes and RNAPs on gfp / active'
				' RNAP.'
				' Proteome share is GFP mass / protein mass. Designs with no'
				' surviving lineage are omitted: %s.'
				% gfp_paper_data.lost_designs(data.designs[
				data.designs['promoter'] > 0]),
			**data.provenance(variants))


if __name__ == '__main__':
	Plot().cli()
