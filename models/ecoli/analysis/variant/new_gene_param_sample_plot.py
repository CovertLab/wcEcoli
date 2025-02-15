"""
Plotting expression factor and translation efficiency parameters used for
new gene expression from new_gene_param_sample variant in order to verify
sampling is uniform across the range of values.
"""

from matplotlib import pyplot as plt
import numpy as np

from models.ecoli.analysis import variantAnalysisPlot
from models.ecoli.sim.variants.new_gene_param_sampling_internal_shift import (
	NEW_GENE_EXPRESSION_FACTOR_CONTROL, NEW_GENE_EXPRESSION_FACTOR_MIN, NEW_GENE_EXPRESSION_FACTOR_MAX,
	NEW_GENE_TRANSLATION_EFFICIENCY_CONTROL, NEW_GENE_TRANSLATION_EFFICIENCY_MIN, NEW_GENE_TRANSLATION_EFFICIENCY_MAX
	)
from wholecell.analysis.analysis_tools import exportFigure

class Plot(variantAnalysisPlot.VariantAnalysisPlot):
	def do_plot(self, inputDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):

		variants = self.ap.get_variants()

		plt.figure()
		plt.xlabel("Expression Factor")
		plt.ylabel("Translation Efficiency")
		plt.xscale('log')
		plt.yscale('log')

		variant_name = metadata["variant"]
		params_to_use = metadata["params_to_use"]

		for variant in variants:

			condition_index = variant // 1000
			index_remainder = variant - condition_index * 1000

			if variant_name == "new_gene_param_sampling_internal_shift":
				from models.ecoli.sim.variants.new_gene_param_sampling_internal_shift import get_sampled_new_gene_expression_factor_and_translation_efficiency
				np.random.seed(index_remainder)
				expression_factor, trl_eff_value = get_sampled_new_gene_expression_factor_and_translation_efficiency(
					index_remainder)

			elif variant_name == "new_gene_param_sampling_internal_shift_narrow":
				if i == 0:
					expression_factors[i]= 0
					trl_eff_values[i] = 0
				else:
					expression_factors[i] = params_to_use[i]["expression_factor"]
					trl_eff_values[i] = params_to_use[i]["trl_eff_value"]

			else:
				print(variant_name + " is not a valid variant name for this plot")
				return

			plt.scatter(expression_factor, trl_eff_value)

		plt.tight_layout()
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close('all')

if __name__ == "__main__":
	Plot().cli()
