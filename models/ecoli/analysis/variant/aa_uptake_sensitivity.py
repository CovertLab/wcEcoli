"""
Compare the growth rates in variants with varying levels of uptake for each
amino acid to see what level of import achieves a growth rate similar to the
validation data.  Useful with the aa_uptake_sensitivity variant.
"""

import pickle

from matplotlib import colors, pyplot as plt
import numpy as np

from models.ecoli.analysis import variantAnalysisPlot
from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from models.ecoli.sim.variants.aa_uptake_sensitivity import get_aa_index, get_factor
from wholecell.analysis.analysis_tools import exportFigure, read_stacked_columns
from wholecell.utils import units


CONTROL_AA = 'L-SELENOCYSTEINE'  # control because SEL is already included for uptake in minimal media
CONTROL_LABEL = 'Control'


class Plot(variantAnalysisPlot.VariantAnalysisPlot):
	def do_plot(self, inputDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		ap = AnalysisPaths(inputDir, variant_plot=True)
		variants = ap.get_variants()

		with open(simDataFile, 'rb') as f:
			sim_data = pickle.load(f)
		with open(validationDataFile, 'rb') as f:
			validation_data = pickle.load(f)

		aa_ids = sim_data.molecule_groups.amino_acids
		aa_labels = [CONTROL_LABEL if aa == CONTROL_AA + '[c]' else aa[:-3] for aa in aa_ids]

		# Load validation growth rates
		val_control = validation_data.amino_acid_growth_rates['minimal']['mean']
		val_aa = []
		val_normalized_growth_rates = []
		for media, rates in validation_data.amino_acid_growth_rates.items():
			aa_id = media.split('_')[-1]
			if aa_id in aa_labels:
				val_aa.append(aa_id)
				val_normalized_growth_rates.append(units.strip_empty_units(rates['mean'] / val_control))

		# Sort order of amino acids by the validation growth rate expected
		sorted_order = (
				[CONTROL_AA]  # Control first
				+ list(np.array(val_aa)[np.argsort(val_normalized_growth_rates)])  # Sorted validation growth rate
				+ [aa for aa in aa_labels if aa not in val_aa and aa != CONTROL_LABEL]  # Amino acids without data
		)
		sorted_idx = {aa: i for i, aa in enumerate(sorted_order)}
		val_aa_idx = [sorted_idx[aa] for aa in val_aa]
		sorted_order[0] = CONTROL_LABEL

		# Load simulation growth rates
		x_scatter = []
		growth_rates = []
		scale = []
		edges = []
		control_rates = []
		for variant in variants:
			# Load data
			cells = ap.get_cells(variant=[variant])
			growth_rate = read_stacked_columns(cells, 'Mass', 'instantaneous_growth_rate',
				remove_first=True, ignore_exception=True).mean()
			aa_id = aa_ids[get_aa_index(variant)][:-3]
			if aa_id == CONTROL_AA:
				control_rates.append(growth_rate)
			factor = get_factor(variant)

			# Save data to plot
			x_scatter.append(sorted_idx[aa_id])
			growth_rates.append(growth_rate)
			if factor == 0:
				scale.append(1)
				edges.append('gray')
			else:
				scale.append(factor)
				edges.append('black')

		# Normalize by the control growth rate
		control_rate = np.mean(control_rates) if control_rates else 1
		normalized_growth_rates = np.array(growth_rates) / control_rate

		# Plot data
		plt.figure()
		scatter = plt.scatter(x_scatter, normalized_growth_rates, c=scale, cmap='RdBu', edgecolors=edges, alpha=0.8, norm=colors.LogNorm())
		plt.scatter(val_aa_idx, val_normalized_growth_rates, marker='_', c='k')
		plt.colorbar(scatter)

		# Plot formatting
		self.remove_border()
		plt.xticks(range(len(sorted_order)), sorted_order, rotation=45, fontsize=6, ha='right')
		plt.yticks(fontsize=6)
		plt.ylabel('Growth rate\n(Normalized to minimal media)', fontsize=6)
		plt.title('Effect of scaling uptake rates by a factor\n'
			'Color represents uptake rate scale factor\n'
			'Horizontal bars are expected growth rates', fontsize=6)

		plt.tight_layout()
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close('all')


if __name__ == "__main__":
	Plot().cli()
