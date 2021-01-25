"""
Plot to visualize basal probabilities and TF bound delta probabilities for
gene expression.  Useful to see which genes are highly expressed and which
ones have low expression or regulation.  Note that some genes will have 0 basal
probability and some TF-gene pairs will have 0 delta probability, which might
be unexpected.
"""

import pickle

from matplotlib import gridspec
from matplotlib import pyplot as plt
import numpy as np

from models.ecoli.analysis import parcaAnalysisPlot
from wholecell.analysis.analysis_tools import exportFigure


class Plot(parcaAnalysisPlot.ParcaAnalysisPlot):
	def do_plot(self, input_dir, plot_out_dir, plot_out_filename, sim_data_file, validation_data_file, metadata):
		with open(sim_data_file, 'rb') as f:
			sim_data = pickle.load(f)

		t_reg = sim_data.process.transcription_regulation
		basal_prob = t_reg.basal_prob
		reg_i = t_reg.delta_prob['deltaI']
		reg_j = t_reg.delta_prob['deltaJ']
		reg_v = t_reg.delta_prob['deltaV']

		# Sort regulation by gene
		gene_sort = np.argsort(reg_i)
		gene_sorted_v = reg_v[gene_sort]

		# Sort regulation by TF
		tf_sort = np.argsort(reg_j)
		tf_sorted_v = reg_v[tf_sort]

		plt.figure(figsize=(6, 12))
		gs = gridspec.GridSpec(nrows=3, ncols=1)

		# Plot sorted basal probabilities
		ax = plt.subplot(gs[0, :])
		ax.bar(range(len(basal_prob)), np.sort(basal_prob))
		ax.set_yscale('symlog', linthreshy=np.min(basal_prob[basal_prob > 0]))
		ax.spines['top'].set_visible(False)
		ax.spines['right'].set_visible(False)
		ax.tick_params(labelsize=6, bottom=False, labelbottom=False)
		ax.set_xlabel('Sorted genes')
		ax.set_ylabel('Basal prob')

		# Plot delta probabilities grouped by genes
		ax = plt.subplot(gs[1, :])
		ax.bar(range(len(gene_sorted_v)), gene_sorted_v)
		ax.set_yscale('symlog', linthreshy=np.min(np.abs(reg_v[np.abs(reg_v) > 0])))
		ax.spines['top'].set_visible(False)
		ax.spines['right'].set_visible(False)
		ax.tick_params(labelsize=6)
		ax.set_ylabel('Delta prob')

		# Plot delta probabilities grouped by transcription factors
		ax = plt.subplot(gs[2, :])
		ax.bar(range(len(gene_sorted_v)), tf_sorted_v)
		ax.set_yscale('symlog', linthreshy=np.min(np.abs(reg_v[np.abs(reg_v) > 0])))
		ax.spines['top'].set_visible(False)
		ax.spines['right'].set_visible(False)
		ax.tick_params(labelsize=6)
		ax.set_ylabel('Delta prob')

		plt.tight_layout()
		exportFigure(plt, plot_out_dir, plot_out_filename, metadata)
		plt.close('all')


if __name__ == "__main__":
	Plot().cli()
