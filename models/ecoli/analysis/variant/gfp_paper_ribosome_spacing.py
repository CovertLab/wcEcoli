"""
Mean spacing between ribosomes on a gfp mRNA, against promoter strength.

# Figure 3E

Spacing is the gfp transcript length divided by ribosomes per gfp mRNA (mean
ribosomes on gfp mRNAs over mean physical gfp mRNA, per cycle). The dashed
line is the ribosome footprint from sim_data, the closest two ribosomes can
sit. One line per RBS strength.

The supporting claim: even the strongest RBS does not bring ribosomes close to
the footprint, so translation is not limited by space on the mRNA.

Batches run before the footprint units fix applied the footprint per codon,
so their effective limit was three footprints (about 72 nt), not one. The
dashed line is the corrected limit; on an older batch it sits below the limit
those simulations actually used.
"""

import numpy as np
from matplotlib import pyplot as plt

from models.ecoli.analysis import variantAnalysisPlot
from wholecell.analysis import gfp_paper_constants as C
from wholecell.analysis import gfp_paper_data
from wholecell.analysis import gfp_paper_style as style

COLUMN = 'ribosome_spacing_nt'


class Plot(variantAnalysisPlot.VariantAnalysisPlot):
	def do_plot(self, inputDir, plotOutDir, plotOutFileName, simDataFile,
			validationDataFile, metadata):
		data = gfp_paper_data.load(inputDir, simDataFile, metadata,
			plotOutDir, self.cpus)
		designs = data.designs[(data.designs['promoter'] > 0)
			& (data.designs['rbs'] > 0)]
		gene = data.gene

		style.apply()
		fig, ax = plt.subplots(figsize=style.panel_size('single', 0.72))
		drawn, any_failed = style.promoter_lines(ax, designs, COLUMN,
			C.RBS_LINES_RIBOSOME)
		style.reference_line(ax, gene.ribosome_footprint_nt,
			'ribosome footprint')
		style.promoter_axis(ax, designs['promoter'])
		values = np.concatenate([d[COLUMN].values for d in drawn]) if drawn \
			else None
		style.log_axis(ax, 'y', 15, 900, data=values)
		ax.set_ylabel('ribosome spacing on gfp mRNA (nt)')
		style.legend_with_failed(ax, any_failed, loc='upper left',
			bbox_to_anchor=(1.0, 0.85))

		variants = sorted(int(v) for d in drawn for v in d['variant'])
		style.save_panel(fig, plotOutDir, style.panel_name(__file__),
			script=__file__, batch=inputDir, sim_data_path=simDataFile,
			n=len(variants),
			notes='Spacing = %g nt / ribosomes per gfp mRNA, window means'
				' over surviving lineages. Footprint %g nt. Designs with no'
				' surviving lineage are omitted: %s.'
				% (gene.length_nt, gene.ribosome_footprint_nt,
				gfp_paper_data.lost_designs(designs)),
			**data.provenance(variants))


if __name__ == '__main__':
	Plot().cli()
