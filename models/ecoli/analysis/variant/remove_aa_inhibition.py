"""
Compare amino acid concentrations across variants similar to data from
Sander et al. Allosteric Feedback Inhibition Enables Robust Amino Acid
Biosynthesis in E. coli by Enforcing Enzyme Overabundance. 2019. Fig 1B.

Associated variant to run:
	remove_aa_inhibition

TODO:
	- add heatmap from Fig 1B in addition to bar plots
"""

import os

from matplotlib import gridspec
from matplotlib import pyplot as plt
import numpy as np

from models.ecoli.analysis import variantAnalysisPlot
from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from models.ecoli.sim.variants.remove_aa_inhibition import AA_TO_ENZYME
from wholecell.analysis.analysis_tools import exportFigure, read_bulk_molecule_counts
from wholecell.io.tablereader import TableReader


def subplot(gs, data, xlabels, title, amino_acids):
	ax = plt.subplot(gs)
	n_variants = data.shape[0]
	x = list(range(n_variants))
	idx = np.array([list(AA_TO_ENZYME.keys()).index(aa) for aa in amino_acids])

	# Plot an amino acid concentration (or sum of multiple amino acids)
	# for each mutant (variant)
	ax.bar(x, data[:, idx].sum(axis=1))

	# Format subplot
	if len(xlabels) == n_variants:
		ax.set_xticks(x)
		ax.set_xticklabels(xlabels, fontsize=8, rotation=45)
	ax.set_title(title, fontsize=8)
	ax.spines['right'].set_visible(False)
	ax.spines['top'].set_visible(False)


class Plot(variantAnalysisPlot.VariantAnalysisPlot):
	def do_plot(self, inputDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		ap = AnalysisPaths(inputDir, variant_plot=True)
		variants = ap.get_variants()

		aa_ids = list(AA_TO_ENZYME.keys())

		aa_conc = np.zeros((len(variants), len(aa_ids)))
		for i, variant in enumerate(variants):
			variant_conc = []
			for sim_dir in ap.get_cells(variant=[variant]):
				simOutDir = os.path.join(sim_dir, "simOut")

				# Listeners used
				kinetics_reader = TableReader(os.path.join(simOutDir, 'EnzymeKinetics'))

				# Read data
				(aa_counts,) = read_bulk_molecule_counts(simOutDir, aa_ids)
				counts_to_molar = kinetics_reader.readColumn('countsToMolar').reshape(-1, 1)

				# Calculate amino acid concentration at each time step
				variant_conc.append(aa_counts[1:, :]*counts_to_molar[1:])

			# Average concentration over all time steps
			aa_conc[i, :] = np.vstack(variant_conc).mean(axis=0)

		# xtick labels
		labels = ['wt'] + list(AA_TO_ENZYME.values())

		# Create figure
		plt.figure(figsize=(10, 10))
		gs = gridspec.GridSpec(nrows=3, ncols=3)

		## Plot subplots for each amino acid
		subplot(gs[0, 0], aa_conc, labels, 'Arginine', ['ARG[c]'])
		subplot(gs[0, 1], aa_conc, labels, 'Tryptophan', ['TRP[c]'])
		subplot(gs[0, 2], aa_conc, labels, 'Histidine', ['HIS[c]'])
		subplot(gs[1, 0], aa_conc, labels, '(Iso-)leucine', ['ILE[c]', 'LEU[c]'])  # group isoforms like paper
		subplot(gs[1, 1], aa_conc, labels, 'Threonine', ['THR[c]'])
		subplot(gs[1, 2], aa_conc, labels, 'Proline', ['PRO[c]'])
		subplot(gs[2, 0], aa_conc, labels, 'Isoleucine', ['ILE[c]'])  # also plot single isoforms since model allows granularity
		subplot(gs[2, 1], aa_conc, labels, 'Leucine', ['LEU[c]'])  # also plot single isoforms since model allows granularity

		## Format and save figure
		plt.tight_layout()
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close('all')


if __name__ == "__main__":
	Plot().cli()
