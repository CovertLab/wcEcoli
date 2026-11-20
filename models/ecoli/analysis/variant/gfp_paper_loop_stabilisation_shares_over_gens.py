"""
GFP's share of transcription, total transcription and GFP transcription by
generation after induction, relative to the induction generation.

# Figure 6C

For the showcase design (nearest SHOWCASE_PROTEOME_SHARE), each lineage's
per-cycle values from the induction generation on, divided by the same
lineage's value in the induction generation, then averaged over the lineages
that reached each generation (dying lineages count until they stop):
	total transcription         all RNAP initiations per second
	GFP transcription           gfp initiations per second
	GFP share of transcription  the second over the first
The dashed line is 1, no change since induction.

The supporting claim: the loop settles because GFP keeps a fixed share of a
transcription budget that shrinks, so its own transcription shrinks with it.
"""

import numpy as np
from matplotlib import pyplot as plt

from models.ecoli.analysis import variantAnalysisPlot
from wholecell.analysis import gfp_paper_constants as C
from wholecell.analysis import gfp_paper_data
from wholecell.analysis import gfp_paper_style as style

LINES = (
	('gfp_transcription_share', 'GFP share of transcription',
		style.GFP_RAMP[2], '-'),
	('rnap_initiations_per_s', 'total transcription', style.INK_MUTED, '-'),
	('gfp_transcription_per_s', 'GFP transcription', style.GFP_RAMP[3],
		(0, (4, 2))),
)


class Plot(variantAnalysisPlot.VariantAnalysisPlot):
	def do_plot(self, inputDir, plotOutDir, plotOutFileName, simDataFile,
			validationDataFile, metadata):
		data = gfp_paper_data.load(inputDir, simDataFile, metadata,
			plotOutDir, self.cpus)
		showcase = gfp_paper_data.nearest_design(data.designs,
			C.SHOWCASE_PROTEOME_SHARE)
		design = int(showcase['variant'])
		cycles = data.cycles[(data.cycles['variant'] == design)
			& data.cycles['divided']
			& (data.cycles['generation'] >= C.NEW_GENE_INDUCTION_GEN)]

		style.apply()
		fig, ax = plt.subplots(figsize=style.panel_size('single', 0.7))
		notes = []
		for column, label, color, dash in LINES:
			table = cycles.pivot_table(index='seed', columns='generation',
				values=column)
			if C.NEW_GENE_INDUCTION_GEN not in table.columns:
				continue
			base = table[C.NEW_GENE_INDUCTION_GEN]
			relative = table.div(base.where(base > 0), axis=0)
			mean = relative.mean(axis=0)
			ax.plot(mean.index, mean.values, color=color, linestyle=dash,
				linewidth=1.2, label=label, zorder=3)
			notes.append('%s: %s' % (label, ', '.join('%d:%.3f' % (g, v)
				for g, v in mean.items())))
		ax.axhline(1.0, color=style.INK_MUTED, linestyle=(0, (4, 3)),
			linewidth=0.8, zorder=1)
		ax.set_xlim(C.NEW_GENE_INDUCTION_GEN, C.ANALYSIS_LAST_GEN)
		ax.set_xticks(range(C.NEW_GENE_INDUCTION_GEN, C.ANALYSIS_LAST_GEN + 1,
			3))
		ax.set_xlabel('generation')
		ax.set_ylabel('relative to induction generation')
		ax.legend(loc='center right')

		style.save_panel(fig, plotOutDir, style.panel_name(__file__),
			script=__file__, batch=inputDir, sim_data_path=simDataFile,
			n=int(cycles['seed'].nunique()),
			notes='Showcase = variant %d (GFP %.3f of proteome). Per lineage'
				' relative to its generation-%d value, mean over lineages per'
				' generation. %s.' % (design, showcase['gfp_proteome_share'],
				C.NEW_GENE_INDUCTION_GEN, '; '.join(notes)),
			**data.provenance([design]))


if __name__ == '__main__':
	Plot().cli()
