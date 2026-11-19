"""
Each gene's share of mRNA copies with GFP against without it.

# Figure 3A

One point per mRNA cistron: its share of all mRNA cistron copies in the
no-GFP control (x) and in the GFP design (y), each the mean over the analysis
window of surviving lineages. Counts are physical (every RNACounts mRNA field
summed) and mapped to cistrons with sim_data's cistron-TU matrix. The dashed
line is no change.

The GFP design is the one whose mean GFP proteome share is closest to
TRANSCRIPTOME_TARGET_PROTEOME_SHARE among designs with surviving lineages; the
share actually reached is in the axis label and the sidecar.

Ribosomal proteins and RNA polymerase subunits are highlighted, since those
are the genes the proteome panels single out. gfp has no copies in the
control, so it is drawn at the left edge of the axis; any other gene below
the axis floor (usually zero copies) in one condition is drawn at that edge
too and counted in the sidecar, and genes with zero copies in both are
omitted. The edge piles shrink as seeds are added, since most are genes
expressed less than once per cell cycle.

The supporting claim: the transcriptome is far less perturbed than the
proteome.
"""

import pickle

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
		try:
			chosen = gfp_paper_data.nearest_design(data.designs,
				C.TRANSCRIPTOME_TARGET_PROTEOME_SHARE)
		except ValueError:
			print('No induced design with surviving lineages; skipping 3A.')
			return
		design = int(chosen['variant'])
		control = gfp_paper_data.CONTROL_VARIANT

		with open(simDataFile, 'rb') as f:
			sim_data = pickle.load(f)
		genes = gfp_paper_data.cistron_profiles(sim_data, data,
			[control, design], plotOutDir, self.cpus)
		shares = genes[[control, design]] / genes[[control, design]].sum()

		style.apply()
		fig, ax = plt.subplots(figsize=style.panel_size('single', 1.0, 0.75))
		counts = style.gene_share_scatter(ax, shares[control].values,
			shares[design].values, genes['is_ribosomal_protein'].values,
			genes['is_RNAP'].values, genes['is_new_gene'].values)
		ax.set_xlabel('share of mRNA copies, no GFP')
		ax.set_ylabel('share of mRNA copies,\nGFP %.0f%% of proteome'
			% (PERCENT * chosen['gfp_proteome_share']))

		style.save_panel(fig, plotOutDir, style.panel_name(__file__),
			script=__file__, batch=inputDir, sim_data_path=simDataFile,
			n=counts['shown'],
			notes='Design: variant %d, promoter %g RBS %g, GFP %.3f of'
				' proteome, %d of %d lineages survived (target share %.2f).'
				' Shares are of all mRNA cistron copies, physical counts.'
				' %d genes below the axis floor in one condition drawn at'
				' the axis edge (including gfp); %d with zero in both'
				' omitted.'
				% (design, chosen['promoter'], chosen['rbs'],
				chosen['gfp_proteome_share'], chosen['n_lineages_survived'],
				chosen['n_lineages_total'],
				C.TRANSCRIPTOME_TARGET_PROTEOME_SHARE, counts['at_floor'],
				counts['omitted']),
			**data.provenance([control, design]))


if __name__ == '__main__':
	Plot().cli()
