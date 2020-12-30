"""
Template for parca analysis plots
"""

import pickle

from matplotlib import pyplot as plt
import numpy as np

from models.ecoli.analysis import parcaAnalysisPlot
from wholecell.analysis.analysis_tools import exportFigure
from wholecell.utils import units


class Plot(parcaAnalysisPlot.ParcaAnalysisPlot):
	def do_plot(self, input_dir, plot_out_dir, plot_out_filename, sim_data_file, validation_data_file, metadata):
		with open(sim_data_file, 'rb') as f:
			sim_data = pickle.load(f)

		# Load kcat data for allosteric reactions
		aa_synthesis_pathways = sim_data.process.metabolism.aa_synthesis_pathways
		amino_acids = sorted(aa_synthesis_pathways)
		kcat_data = np.array([
			aa_synthesis_pathways[aa]['kcat_data'].asNumber(1 / units.s)
			for aa in amino_acids
			])
		kcat_calc = np.array([
			aa_synthesis_pathways[aa]['kcat'].asNumber(1 / units.s)
			for aa in amino_acids
			])

		# Get kcat limits
		kcat_min = min(kcat_data[kcat_data > 0].min(), kcat_calc[kcat_calc > 0].min())
		kcat_max = max(kcat_data.max(), kcat_calc.max())
		kcat_range = [kcat_min, kcat_max]

		# Plot data
		plt.figure()

		plt.loglog(kcat_data, kcat_calc, 'o')
		plt.loglog(kcat_range, kcat_range, 'k--')

		plt.tight_layout()
		exportFigure(plt, plot_out_dir, plot_out_filename, metadata)
		plt.close('all')


if __name__ == "__main__":
	Plot().cli()
