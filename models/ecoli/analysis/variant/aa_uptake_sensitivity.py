"""
Compare cell cycle times and growth rates across variants.  Useful as validation
with the add_one_aa variant and can also be used to compare variants in
remove_one_aa variant.
"""

import pickle

from matplotlib import colors, pyplot as plt
import numpy as np

from models.ecoli.analysis import variantAnalysisPlot
from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from models.ecoli.sim.variants.aa_uptake_sensitivity import get_aa_index, get_factor
from wholecell.analysis.analysis_tools import exportFigure, read_stacked_columns
from wholecell.utils import units


CONTROL_LABEL = 'L-SELENOCYSTEINE'  # control because SEL is already included for uptake in minimal media


class Plot(variantAnalysisPlot.VariantAnalysisPlot):
	def do_plot(self, inputDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		ap = AnalysisPaths(inputDir, variant_plot=True)
		variants = ap.get_variants()

		with open(simDataFile, 'rb') as f:
			sim_data = pickle.load(f)
		with open(validationDataFile, 'rb') as f:
			validation_data = pickle.load(f)

		aa_ids = sim_data.molecule_groups.amino_acids
		aa_idx = {aa[:-3]: i for i, aa in enumerate(aa_ids)}
		aa_labels = ['Control' if aa == CONTROL_LABEL + '[c]' else aa[:-3] for aa in aa_ids]

		growth_rates = {}
		control_rates = []
		for variant in variants:
			cells = ap.get_cells(variant=[variant])
			growth_rate = read_stacked_columns(cells, 'Mass', 'instantaneous_growth_rate',
				remove_first=True, ignore_exception=True).mean()
			aa_id = aa_ids[get_aa_index(variant)][:-3]
			if aa_id == CONTROL_LABEL:
				control_rates.append(growth_rate)
			growth_rates[(aa_id, get_factor(variant))] = growth_rate
		control_rate = np.mean(control_rates) if control_rates else 1

		# Load validation growth rates
		val_control = validation_data.amino_acid_growth_rates['minimal']['mean']
		val_aa_idx = []
		val_normalized_growth_rates = []
		for media, rates in validation_data.amino_acid_growth_rates.items():
			aa_id = media.split('_')[-1]
			if aa_id in aa_idx:
				val_aa_idx.append(aa_idx[aa_id])
				val_normalized_growth_rates.append(units.strip_empty_units(rates['mean'] / val_control))

		# Data to plot for each amino acid condition
		x = []
		y = []
		c = []
		edges = []
		for (aa, factor), rate in growth_rates.items():
			color = 'k'
			x.append(aa_idx[aa])
			y.append(rate / control_rate)
			if factor == 0:
				color = 'gray'
				c.append(1)
			else:
				c.append(factor)
			edges.append(color)

		# Create plots
		plt.figure()

		# Plot data
		scatter = plt.scatter(x, y, c=c, cmap='RdBu', edgecolors=edges, alpha=0.8, norm=colors.LogNorm())
		plt.scatter(val_aa_idx, val_normalized_growth_rates, marker='_', c='k')
		plt.colorbar(scatter)

		# Plot formatting
		plt.xticks(range(len(aa_labels)), aa_labels, rotation=45, fontsize=6, ha='right')
		plt.ylabel('Growth rate\n(Normalized to minimal media)')

		plt.tight_layout()
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close('all')


if __name__ == "__main__":
	Plot().cli()
