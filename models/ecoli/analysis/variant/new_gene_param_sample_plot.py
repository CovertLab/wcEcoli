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

		for variant in variants:
			np.random.seed(variant)
			if variant == 0:
				expression_factor = NEW_GENE_EXPRESSION_FACTOR_CONTROL
				trl_eff_value = NEW_GENE_TRANSLATION_EFFICIENCY_CONTROL
			else:
				expression_factor = np.random.uniform(NEW_GENE_EXPRESSION_FACTOR_MIN, NEW_GENE_EXPRESSION_FACTOR_MAX)
				trl_eff_value = 10 ** np.random.uniform(NEW_GENE_TRANSLATION_EFFICIENCY_MIN,NEW_GENE_TRANSLATION_EFFICIENCY_MAX)

			plt.scatter(expression_factor, trl_eff_value)

		plt.tight_layout()
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close('all')

if __name__ == "__main__":
	Plot().cli()
