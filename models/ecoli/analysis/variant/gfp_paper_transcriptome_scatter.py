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

import numpy as np
from matplotlib import pyplot as plt

from models.ecoli.analysis import variantAnalysisPlot
from wholecell.analysis import gfp_paper_constants as C
from wholecell.analysis import gfp_paper_data
from wholecell.analysis import gfp_paper_style as style

AXIS_FLOOR = 1e-8
BACKGROUND = '#CDD1CF'      # all other genes: neutral, not the control grey
PERCENT = 100.0


class Plot(variantAnalysisPlot.VariantAnalysisPlot):
	def do_plot(self, inputDir, plotOutDir, plotOutFileName, simDataFile,
			validationDataFile, metadata):
		data = gfp_paper_data.load(inputDir, simDataFile, metadata,
			plotOutDir, self.cpus)
		candidates = data.designs[(data.designs['promoter'] > 0)
			& (data.designs['n_lineages_survived'] > 0)
			& np.isfinite(data.designs['gfp_proteome_share'])]
		if candidates.empty:
			print('No induced design with surviving lineages; skipping 3A.')
			return
		gap = (candidates['gfp_proteome_share']
			- C.TRANSCRIPTOME_TARGET_PROTEOME_SHARE).abs()
		chosen = candidates.loc[gap.idxmin()]
		design = int(chosen['variant'])
		control = gfp_paper_data.CONTROL_VARIANT

		with open(simDataFile, 'rb') as f:
			sim_data = pickle.load(f)
		genes = gfp_paper_data.cistron_profiles(inputDir, sim_data, data,
			[control, design], plotOutDir, self.cpus)
		shares = genes[[control, design]] / genes[[control, design]].sum()
		x, y = shares[control].values, shares[design].values
		keep = (x > 0) | (y > 0)
		edge = keep & ((x < AXIS_FLOOR) | (y < AXIS_FLOOR))
		# Zero, and the rare share below the axis floor, sit on the floor.
		x = np.maximum(x, AXIS_FLOOR)
		y = np.maximum(y, AXIS_FLOOR)

		style.apply()
		fig, ax = plt.subplots(figsize=style.panel_size('single', 1.0, 0.75))
		ax.plot([AXIS_FLOOR, 1], [AXIS_FLOOR, 1], color=style.INK_MUTED,
			linestyle=(0, (4, 3)), linewidth=0.8, zorder=1)
		groups = [
			('all other genes', BACKGROUND, 6,
				~(genes['is_ribosomal_protein'] | genes['is_RNAP']
				| genes['is_new_gene']).values),
			('ribosomal proteins', style.COLORS['ribosome'], 12,
				genes['is_ribosomal_protein'].values),
			('RNA polymerase subunits', style.COLORS['rnap'], 12,
				genes['is_RNAP'].values),
			('GFP', style.COLORS['gfp'], 30, genes['is_new_gene'].values),
		]
		for label, color, size, mask in groups:
			mask = mask & keep
			ax.scatter(x[mask], y[mask], s=size, c=color, linewidths=0.3,
				edgecolors='white', label=label, zorder=2 if
				label == 'all other genes' else 3, clip_on=False)

		for setter in (ax.set_xscale, ax.set_yscale):
			setter('log')
		ticks = [1e-8, 1e-6, 1e-4, 1e-2, 1]
		ax.set_xticks(ticks)
		ax.set_yticks(ticks)
		ax.set_xticklabels(['$10^{-8}$', '$10^{-6}$', '$10^{-4}$',
			'$10^{-2}$', '1'])
		ax.set_yticklabels(['$10^{-8}$', '$10^{-6}$', '$10^{-4}$',
			'$10^{-2}$', '1'])
		ax.minorticks_off()
		ax.set_xlim(AXIS_FLOOR, 1)
		ax.set_ylim(AXIS_FLOOR, 1)
		ax.set_aspect('equal')
		ax.set_xlabel('share of mRNA copies, no GFP')
		ax.set_ylabel('share of mRNA copies,\nGFP %.0f%% of proteome'
			% (PERCENT * chosen['gfp_proteome_share']))
		ax.legend(loc='upper left', handletextpad=0.3)

		style.save_panel(fig, plotOutDir, style.panel_name(__file__),
			script=__file__, batch=inputDir, sim_data_path=simDataFile,
			n=int(keep.sum()),
			notes='Design: variant %d, promoter %g RBS %g, GFP %.3f of'
				' proteome, %d of %d lineages survived (target share %.2f).'
				' Shares are of all mRNA cistron copies, physical counts.'
				' %d genes below the axis floor in one condition drawn at'
				' the axis edge (including gfp); %d with zero in both'
				' omitted.'
				% (design, chosen['promoter'], chosen['rbs'],
				chosen['gfp_proteome_share'], chosen['n_lineages_survived'],
				chosen['n_lineages_total'],
				C.TRANSCRIPTOME_TARGET_PROTEOME_SHARE, int(edge.sum()),
				int((~keep).sum())),
			**data.provenance([control, design]))


if __name__ == '__main__':
	Plot().cli()
