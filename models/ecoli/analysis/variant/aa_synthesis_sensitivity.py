"""
Compare the growth rates in variants with varying synthesis parameters for each
amino acid to see how sensitive growth rate is to a range of parameters.
Useful with the aa_synthesis_sensitivity variant.
"""

import pickle

from matplotlib import colors, pyplot as plt
import numpy as np

from models.ecoli.analysis import variantAnalysisPlot
from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from models.ecoli.sim.variants import aa_synthesis_sensitivity
from wholecell.analysis.analysis_tools import exportFigure, read_stacked_columns


GLT_INDEX = 5
CONTROL_INDEX = 19


class Plot(variantAnalysisPlot.VariantAnalysisPlot):
	def do_plot(self, inputDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		ap = AnalysisPaths(inputDir, variant_plot=True)
		variants = ap.get_variants()

		with open(simDataFile, 'rb') as f:
			sim_data = pickle.load(f)
		aa_ids = sim_data.molecule_groups.amino_acids

		n_aas = 1  # TODO; generalize to all and match with the sim variant file

		# Load simulation growth rates
		data = {}
		for variant in variants:
			media_index = aa_synthesis_sensitivity.get_media_index(variant, n_aas)
			aa_index = aa_synthesis_sensitivity.get_aa_index(variant, n_aas)
			param, factor = aa_synthesis_sensitivity.get_adjustment(variant)
			aa_adjusted = aa_ids[aa_index][:-3]
			param_label = f'{aa_adjusted} {param}'

			# Load data
			cells = ap.get_cells(variant=[variant])
			growth_rate = read_stacked_columns(cells, 'Mass', 'instantaneous_growth_rate',
				remove_first=True, ignore_exception=True).mean()
			elong_rate = read_stacked_columns(cells, 'RibosomeData', 'effectiveElongationRate',
				remove_first=True, ignore_exception=True).mean()

			variant_data = {'Growth rate': growth_rate * 3600, 'Elongation rate': elong_rate}
			media_data = data.get(media_index, {})
			param_data = media_data.get(param_label, {})
			param_data[factor] = variant_data
			media_data[param_label] = param_data
			data[media_index] = media_data

		# Plot data
		plt.figure(figsize=(4, 16))
		n_subplots = 4

		def subplot(label, data_index, subplot_index, normalized=True):
			plt.subplot(n_subplots, 1, subplot_index)
			glt_data = data[data_index]  # TODO: generalize to all
			param_labels = list(glt_data.keys())
			label_indices = {param: i for i, param in enumerate(param_labels)}

			x = []
			y = []
			scale = []
			for param_label, factors in glt_data.items():
				for factor, values in factors.items():
					x.append(label_indices[param_label])
					denom = data.get(CONTROL_INDEX, {}).get(param_label, {}).get(factor, {}).get(label, 1) if normalized else 1
					y.append(values[label] / denom)
					scale.append(factor)
			scatter = plt.scatter(x, y, c=scale, cmap='RdBu', alpha=0.8, norm=colors.LogNorm())
			plt.colorbar(scatter)

			## Plot formatting
			normalized_text = '\n(Normalized to minimal media)' if normalized else ''
			self.remove_border()
			plt.xticks(range(len(param_labels)), param_labels, rotation=45, fontsize=6, ha='right')
			plt.yticks(fontsize=6)
			plt.ylabel(f'{label}{normalized_text}', fontsize=6)

		subplot('Growth rate', GLT_INDEX, 1)
		plt.title('Effect of scaling synthesis parameters by a factor\n'
			'Color represents parameter scale factor', fontsize=6)
		subplot('Growth rate', GLT_INDEX, 2, normalized=False)
		subplot('Growth rate', CONTROL_INDEX, 3, normalized=False)
		subplot('Elongation rate', GLT_INDEX, 4)

		plt.tight_layout()
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close('all')


if __name__ == "__main__":
	Plot().cli()
