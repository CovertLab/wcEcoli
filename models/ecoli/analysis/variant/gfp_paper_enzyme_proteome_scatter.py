"""
Each protein's share of the proteome with GFP against without, metabolic
enzymes highlighted.

# Figure 4A

The same scatter as 1C (mass shares of all monomers, window means over
surviving lineages, showcase design nearest SHOWCASE_PROTEOME_SHARE against
the no-GFP control), with the subunits of every FBA catalyst highlighted as
metabolic enzymes. This is the broad set; the kcat-bounded subset that 4B,
4C and 4G follow is a part of it.

The supporting claim: in the GFP-driven proteome shift, metabolic enzymes
fall with everything else.
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
		showcase = gfp_paper_data.nearest_design(data.designs,
			C.SHOWCASE_PROTEOME_SHARE)
		design = int(showcase['variant'])
		control = gfp_paper_data.CONTROL_VARIANT
		profiles = gfp_paper_data.vector_profiles(data, 'monomers',
			[control, design], plotOutDir, self.cpus)
		ids, shares = gfp_paper_data.monomer_mass_shares(sim_data, profiles)
		enzymes = set(gfp_paper_data.catalyst_monomers(sim_data,
			sim_data.process.metabolism.catalyst_ids))
		is_enzyme = np.array([i in enzymes for i in ids])
		new = np.array([i == data.gene.monomer_id for i in ids])

		style.apply()
		fig, ax = plt.subplots(figsize=style.panel_size('single', 1.0, 0.75))
		counts = style.gene_share_scatter(ax, shares[control], shares[design],
			None, None, new, highlights=[('metabolic enzymes',
			style.COLORS['enzyme'], is_enzyme & ~new)])
		ax.set_xlabel('proteome fraction, no GFP')
		ax.set_ylabel('proteome fraction,\nGFP %.0f%% of proteome'
			% (PERCENT * showcase['gfp_proteome_share']))
		enzyme_share = [float(np.nansum(shares[v][is_enzyme]))
			for v in (control, design)]

		style.save_panel(fig, plotOutDir, style.panel_name(__file__),
			script=__file__, batch=inputDir, sim_data_path=simDataFile,
			n=counts['shown'],
			notes='Showcase = variant %d. Metabolic enzymes = subunits of all'
				' %d FBA catalysts (%d monomers present). Their summed'
				' proteome share: %.3f without GFP, %.3f with. %d at the axis'
				' floor in one condition, %d omitted.' % (design,
				len(sim_data.process.metabolism.catalyst_ids),
				int(is_enzyme.sum()), enzyme_share[0], enzyme_share[1],
				counts['at_floor'], counts['omitted']),
			**data.provenance([control, design]))


if __name__ == '__main__':
	Plot().cli()
