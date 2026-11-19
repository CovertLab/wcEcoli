"""
Each protein's share of the proteome with GFP against without it.

# Figure 1C

One point per monomer: its share of total protein mass (count x molecular
weight over the sum of the same for every monomer; counts include monomers
inside complexes) in the no-GFP control (x) and the showcase design (y),
window means over surviving lineages. Ribosomal proteins and RNAP subunits
are highlighted; gfp is absent without GFP and drawn at the left edge.
Drawn with the same helper as the transcriptome scatter (3A).

The showcase design is the one nearest SHOWCASE_PROTEOME_SHARE (decision D3).
For the attenuation check that informs D3, the sidecar lists the proteins
whose share moves most (largest absolute log2 change among those above
1e-5 of the proteome in either condition), excluding gfp.

The supporting claim: GFP shifts the whole proteome, and machinery proteins
fall with the rest.
"""

import pickle

import numpy as np
from matplotlib import pyplot as plt

from models.ecoli.analysis import variantAnalysisPlot
from wholecell.analysis import gfp_paper_constants as C
from wholecell.analysis import gfp_paper_data
from wholecell.analysis import gfp_paper_style as style

PERCENT = 100.0
MOVER_FLOOR = 1e-5
N_MOVERS = 15


class Plot(variantAnalysisPlot.VariantAnalysisPlot):
	def do_plot(self, inputDir, plotOutDir, plotOutFileName, simDataFile,
			validationDataFile, metadata):
		from wholecell.utils import units
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

		monomers = sim_data.process.translation.monomer_data
		mw = dict(zip(monomers['id'], monomers['mw'].asNumber(
			units.g / units.mol)))
		cistrons = sim_data.process.transcription.cistron_data
		flags = dict(zip(cistrons['id'], zip(cistrons['is_ribosomal_protein'],
			cistrons['is_RNAP'])))
		cistron_of = dict(zip(monomers['id'], monomers['cistron_id']))
		ids = profiles[control].ids
		weights = np.array([mw.get(i, np.nan) for i in ids])
		shares = {}
		for variant in (control, design):
			mass = profiles[variant].mean * weights
			shares[variant] = mass / np.nansum(mass)
		ribosomal = np.array([flags.get(cistron_of.get(i), (False,))[0]
			for i in ids])
		rnap = np.array([flags.get(cistron_of.get(i), (False, False))[1]
			for i in ids])
		new = np.array([i == data.gene.monomer_id for i in ids])

		style.apply()
		fig, ax = plt.subplots(figsize=style.panel_size('single', 1.0, 0.75))
		counts = style.gene_share_scatter(ax, shares[control], shares[design],
			ribosomal, rnap, new)
		ax.set_xlabel('proteome fraction, no GFP')
		ax.set_ylabel('proteome fraction,\nGFP %.0f%% of proteome'
			% (PERCENT * showcase['gfp_proteome_share']))

		x, y = shares[control], shares[design]
		ok = ~new & ((x > MOVER_FLOOR) | (y > MOVER_FLOOR)) & (x > 0) & (y > 0)
		change = np.full(len(ids), np.nan)
		change[ok] = np.log2(y[ok] / x[ok])
		order = np.argsort(-np.abs(np.nan_to_num(change)))[:N_MOVERS]
		movers = ', '.join('%s %+.2f' % (ids[i], change[i]) for i in order
			if np.isfinite(change[i]))

		style.save_panel(fig, plotOutDir, style.panel_name(__file__),
			script=__file__, batch=inputDir, sim_data_path=simDataFile,
			n=counts['shown'],
			notes='Showcase = variant %d (GFP %.3f of proteome). Mass shares'
				' of all monomers. %d below the axis floor in one condition'
				' drawn at the edge (including gfp); %d at zero in both'
				' omitted. Largest log2 changes (share > %g): %s.'
				% (design, showcase['gfp_proteome_share'], counts['at_floor'],
				counts['omitted'], MOVER_FLOOR, movers),
			**data.provenance([control, design]))


if __name__ == '__main__':
	Plot().cli()
