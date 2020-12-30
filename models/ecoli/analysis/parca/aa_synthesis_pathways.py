"""
Template for parca analysis plots
"""

import pickle

from matplotlib import pyplot as plt
from matplotlib import gridspec
import numpy as np

from models.ecoli.analysis import parcaAnalysisPlot
from wholecell.analysis.analysis_tools import exportFigure
from wholecell.utils import units


class Plot(parcaAnalysisPlot.ParcaAnalysisPlot):
	def do_plot(self, input_dir, plot_out_dir, plot_out_filename, sim_data_file, validation_data_file, metadata):
		with open(sim_data_file, 'rb') as f:
			sim_data = pickle.load(f)

		metabolism = sim_data.process.metabolism
		aa_synthesis_pathways = metabolism.aa_synthesis_pathways
		conc = metabolism.concentration_updates.concentrations_based_on_nutrients

		# Load kcat data for allosteric reactions
		amino_acids = sorted(aa_synthesis_pathways)
		kcat_data = np.array([
			aa_synthesis_pathways[aa]['kcat_data'].asNumber(1 / units.s)
			for aa in amino_acids
			])
		kcat_calc = np.array([
			aa_synthesis_pathways[aa]['kcat'].asNumber(1 / units.s)
			for aa in amino_acids
			])

		# Get kcat limits for plot
		kcat_min = min(kcat_data[kcat_data > 0].min(), kcat_calc[kcat_calc > 0].min())
		kcat_max = max(kcat_data.max(), kcat_calc.max())
		kcat_range = [kcat_min, kcat_max]

		# Calculate expected inhibition in minimal media condition
		aa_inhibition = np.array([
			[1 / (1 + conc('minimal')[aa] / ki).asNumber() for ki in aa_synthesis_pathways[aa]['ki']]
			for aa in amino_acids
			])

		# Ranges for plot
		x_amino_acids = np.arange(len(amino_acids))

		# Plot data
		plt.figure(figsize=(5, 10))
		gs = gridspec.GridSpec(2, 1)

		## kcat comparison
		plt.subplot(gs[0, 0])
		plt.loglog(kcat_data, kcat_calc, 'o')
		plt.loglog(kcat_range, kcat_range, 'k--')

		## Inhibition comparison
		plt.subplot(gs[1, 0])
		plt.bar(x_amino_acids, aa_inhibition[:, 1], label='upper limit')
		plt.bar(x_amino_acids, aa_inhibition[:, 0], label='lower limit')
		plt.xticks(x_amino_acids, amino_acids, rotation=45, fontsize=6)
		plt.legend()

		plt.tight_layout()
		exportFigure(plt, plot_out_dir, plot_out_filename, metadata)
		plt.close('all')


if __name__ == "__main__":
	Plot().cli()
