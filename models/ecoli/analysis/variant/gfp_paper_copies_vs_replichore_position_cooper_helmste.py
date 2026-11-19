"""
Gene copies per cell along the replichore, simulated against the
Cooper-Helmstetter prediction, without and with high GFP.

# Figure 5B

Points: each gene's window-mean copies (RnaSynthProb gene_copy_number),
averaged over surviving lineages, then averaged (unweighted) over the genes in
each of N_BINS equal bins of replichore fraction (0 = oriC, 1 = terC, each arm
normalised by its own length). Dashed lines: Cooper-Helmstetter,
2 ** (((1 - f) C + D) / tau), from replication.get_average_copy_number, with
C and D from sim_data and tau the design's own mean doubling time.

"High GFP" is the design nearest the second entry of FIG5_PROTEOME_SHARES.
The supporting claim: the model's dosage follows Cooper-Helmstetter with and
without burden, so slower growth means fewer copies, most of all near oriC.
"""

import pickle

import numpy as np
from matplotlib import pyplot as plt

from models.ecoli.analysis import variantAnalysisPlot
from wholecell.analysis import gfp_paper_constants as C
from wholecell.analysis import gfp_paper_data
from wholecell.analysis import gfp_paper_style as style

N_BINS = 20
PERCENT = 100.0


class Plot(variantAnalysisPlot.VariantAnalysisPlot):
	def do_plot(self, inputDir, plotOutDir, plotOutFileName, simDataFile,
			validationDataFile, metadata):
		data = gfp_paper_data.load(inputDir, simDataFile, metadata,
			plotOutDir, self.cpus)
		with open(simDataFile, 'rb') as f:
			sim_data = pickle.load(f)
		high = gfp_paper_data.nearest_design(data.designs,
			C.FIG5_PROTEOME_SHARES[-1])
		control = gfp_paper_data.CONTROL_VARIANT
		variants = [control, int(high['variant'])]
		profiles = gfp_paper_data.vector_profiles(data, 'gene_copies',
			variants, plotOutDir, self.cpus)

		cistrons = sim_data.process.transcription.cistron_data
		coordinate = dict(zip(cistrons['gene_id'],
			cistrons['replication_coordinate']))
		replication = sim_data.process.replication
		right = replication.replichore_lengths[0]
		designs = data.designs.set_index('variant')

		style.apply()
		fig, ax = plt.subplots(figsize=style.panel_size('single', 0.62))
		edges = np.linspace(0, 1, N_BINS + 1)
		centres = (edges[:-1] + edges[1:]) / 2
		grid = np.linspace(0, 1, 101)
		notes = []
		for variant in variants:
			profile = profiles[variant]
			if not profile.ids:
				continue
			coords = np.array([coordinate.get(g, np.nan)
				for g in profile.ids], dtype=float)
			f = gfp_paper_data.replichore_fraction(sim_data, coords)
			copies = profile.mean
			ok = np.isfinite(f) & np.isfinite(copies)
			# right=True keeps a gene at exactly f = 1 in the last bin
			bins = np.clip(np.digitize(f[ok], edges, right=True) - 1, 0,
				N_BINS - 1)
			means = np.array([copies[ok][bins == b].mean()
				if np.any(bins == b) else np.nan for b in range(N_BINS)])
			tau = float(designs.loc[variant, 'doubling_time_min'])
			predicted = replication.get_average_copy_number(tau,
				grid * right)
			if variant == control:
				color, name = style.COLORS['control'], 'no GFP'
			else:
				color = style.COLORS['gfp']
				name = 'GFP %.0f%%' % (PERCENT * high['gfp_proteome_share'])
			ax.scatter(centres, means, s=14, color=color, zorder=3,
				label='%s, simulated' % name)
			ax.plot(grid, predicted, color=color, linestyle=(0, (4, 3)),
				linewidth=1.0, zorder=2,
				label='%s, Cooper-Helmstetter' % name)
			notes.append('%s: tau %.1f min, %d genes' % (name, tau, ok.sum()))

		ax.set_xlim(0, 1)
		ax.set_xticks([0, 0.25, 0.5, 0.75, 1])
		ax.set_xticklabels(['0', '0.25', '0.5', '0.75', '1'])
		ax.set_xlabel('replichore position (0 = oriC, 1 = terC)')
		ax.set_ylabel('gene copies per cell')
		ax.legend(loc='upper right')

		style.save_panel(fig, plotOutDir, style.panel_name(__file__),
			script=__file__, batch=inputDir, sim_data_path=simDataFile,
			n=N_BINS,
			notes='%d equal replichore bins, unweighted gene means. C = %.1f'
				' min, D = %.1f min from sim_data. %s. High GFP = variant %d.'
				% (N_BINS, replication.c_period_in_mins,
				replication.d_period_in_mins, '; '.join(notes),
				high['variant']),
			**data.provenance(variants))


if __name__ == '__main__':
	Plot().cli()
