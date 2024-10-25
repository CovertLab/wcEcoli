"""
plotting parameters from new_gene_param_sample variant
"""

from matplotlib import pyplot as plt
# noinspection PyUnresolvedReferences
import numpy as np

from models.ecoli.analysis import variantAnalysisPlot
from wholecell.analysis.analysis_tools import exportFigure



class Plot(variantAnalysisPlot.VariantAnalysisPlot):
	def do_plot(self, inputDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):

		variants = self.ap.get_variants()



		NEW_GENE_EXPRESSION_FACTOR_CONTROL = 0
		NEW_GENE_EXPRESSION_FACTOR_MIN = 7
		NEW_GENE_EXPRESSION_FACTOR_MAX = 10
		NEW_GENE_TRANSLATION_EFFICIENCY_CONTROL = 0
		NEW_GENE_TRANSLATION_EFFICIENCY_MIN = np.log10(0.01)
		NEW_GENE_TRANSLATION_EFFICIENCY_MAX = np.log10(10)

		plt.figure()
		plt.xlabel("Expression Factor")
		plt.ylabel("Translation Efficiency")
		plt.xscale('log')
		plt.yscale('log')

		for variant in variants:

				#if variant >1:
					#continue
				np.random.seed(variant)
				if variant == 0:
					expression_factor = NEW_GENE_EXPRESSION_FACTOR_CONTROL
					trl_eff_value = NEW_GENE_TRANSLATION_EFFICIENCY_CONTROL
				else:
					expression_factor = np.random.uniform(NEW_GENE_EXPRESSION_FACTOR_MIN, NEW_GENE_EXPRESSION_FACTOR_MAX)
					trl_eff_value = 10 ** np.random.uniform(NEW_GENE_TRANSLATION_EFFICIENCY_MIN,NEW_GENE_TRANSLATION_EFFICIENCY_MAX)

				plt.scatter(expression_factor, trl_eff_value)



		### Create Plot ###

		plt.tight_layout()
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close('all')


if __name__ == "__main__":
	Plot().cli()
from models.ecoli.analysis import analysisPlot
from runscripts.manual import analysisVariant











