"""
Template for parca analysis plots
"""

import pickle

from matplotlib import pyplot as plt
import numpy as np

from models.ecoli.analysis import parcaAnalysisPlot
from reconstruction.ecoli.initialization import create_bulk_container
from wholecell.analysis.analysis_tools import exportFigure


class Plot(parcaAnalysisPlot.ParcaAnalysisPlot):
	def do_plot(self, input_dir, plot_out_dir, plot_out_filename, sim_data_file, validation_data_file, metadata):
		with open(sim_data_file, 'rb') as f:
			sim_data = pickle.load(f)
		with open(validation_data_file, 'rb') as f:
			validation_data = pickle.load(f)

		# Expected bulk containers in different conditions
		rich_container = create_bulk_container(sim_data, condition='with_aa')
		basal_container = create_bulk_container(sim_data)

		# TODO: convert complexes into monomers

		# TODO: select molecule groups of interest

		# TODO: get mass fraction from validation
		rich_validation = None
		basal_validation = None
		expected_change = np.sign(rich_validation - basal_validation)


		plt.figure()

		### Create Plot ###

		plt.tight_layout()
		exportFigure(plt, plot_out_dir, plot_out_filename, metadata)
		plt.close('all')


if __name__ == "__main__":
	Plot().cli()
