"""
Mean spacing between RNA polymerases on gfp, against promoter strength.

# Figure 3D

Spacing is the gfp transcription unit length divided by the RNAPs on each gfp
copy: RNAPs on gfp (transcripts still being made, one polymerase each) over
gfp promoter copies, both cycle means. The dashed line is the RNAP footprint
from sim_data, the closest two polymerases can sit; spacing that reaches it
means transcription initiation is at the model's spacing cap.

One line per RBS strength. RBS should barely matter here, since ribosomes do
not change how often RNAP starts on gfp; a spread between the lines would
point to an indirect effect, such as RNAP pool size.

The supporting claim: the strongest promoters run into the RNAP spacing limit.
"""

import numpy as np
from matplotlib import pyplot as plt

from models.ecoli.analysis import variantAnalysisPlot
from wholecell.analysis import gfp_paper_constants as C
from wholecell.analysis import gfp_paper_data
from wholecell.analysis import gfp_paper_style as style

COLUMN = 'rnap_spacing_nt'


class Plot(variantAnalysisPlot.VariantAnalysisPlot):
	def do_plot(self, inputDir, plotOutDir, plotOutFileName, simDataFile,
			validationDataFile, metadata):
		data = gfp_paper_data.load(inputDir, simDataFile, metadata,
			plotOutDir, self.cpus)
		designs = data.designs[data.designs['promoter'] > 0]
		gene = data.gene

		style.apply()
		fig, ax = plt.subplots(figsize=style.panel_size('single', 0.72))
		drawn, any_failed = style.promoter_lines(ax, designs, COLUMN,
			C.RBS_LINES_RNAP_SPACING)
		style.reference_line(ax, gene.rnap_footprint_nt, 'RNAP footprint')
		style.promoter_axis(ax, designs['promoter'])
		values = np.concatenate([d[COLUMN].values for d in drawn]) if drawn \
			else None
		style.log_axis(ax, 'y', 40, 900, data=values)
		ax.set_ylabel('spacing between RNAPs on gfp (nt)')
		style.legend_with_failed(ax, any_failed, loc='upper right')

		variants = sorted(int(v) for d in drawn for v in d['variant'])
		style.save_panel(fig, plotOutDir, style.panel_name(__file__),
			script=__file__, batch=inputDir, sim_data_path=simDataFile,
			n=len(variants),
			notes='Spacing = %g nt / (RNAPs on gfp / gfp promoter copies),'
				' window means over surviving lineages. RNAP counted from %s.'
				' Footprint %g nt. Designs with no surviving lineage are'
				' omitted: %s.' % (gene.length_nt, _rnap_field(data),
				gene.rnap_footprint_nt,
				gfp_paper_data.lost_designs(designs)),
			**data.provenance(variants))


def _rnap_field(data):
	fields = data.cycles.get('rnap_field')
	return ', '.join(sorted(fields.dropna().unique())) if fields is not None \
		else 'unknown'


if __name__ == '__main__':
	Plot().cli()
