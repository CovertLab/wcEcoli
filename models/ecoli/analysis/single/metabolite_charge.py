"""
Shows aggregated charge of metabolites with known concentration by cell location
"""

import pickle
import os

from matplotlib import pyplot as plt
import numpy as np

from models.ecoli.analysis import singleAnalysisPlot
from wholecell.analysis.plotting_tools import COLORS_LARGE
from wholecell.analysis.analysis_tools import exportFigure, read_bulk_molecule_counts
from wholecell.io.tablereader import TableReader


class Plot(singleAnalysisPlot.SingleAnalysisPlot):
	def do_plot(self, simOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		with open(simDataFile, 'rb') as f:
			sim_data = pickle.load(f)
		metabolite_charge = sim_data.metabolite_charge

		# Listeners used
		enzymeKineticsdata = TableReader(os.path.join(simOutDir, "EnzymeKinetics"))
		main_reader = TableReader(os.path.join(simOutDir, "Main"))

		# Metabolite data
		metaboliteNames = enzymeKineticsdata.readAttribute("metaboliteNames")

		# Load data
		initial_time = main_reader.readAttribute('initialTime')
		time = (main_reader.readColumn('time') - initial_time)/60

		(counts,) = read_bulk_molecule_counts(simOutDir, (metaboliteNames,))

		charge_vector = [metabolite_charge[m[:-3]] for m in metaboliteNames]
		charge_t = counts*charge_vector
		m_location = [m[-2] for m in metaboliteNames]

		# Separate charge data by location
		locations = np.unique(m_location)
		loc_mask = {}
		loc_charge = {}
		for loc in locations:
			loc_mask[loc] = [ml == loc for ml in m_location]
			mask = loc_mask[loc]
			loc_charge[loc] = np.sum(charge_t[:,mask], axis=1)

		# Plot total charge by location
		colors = COLORS_LARGE
		plt.figure(figsize = (8, 8))
		nrows = len(loc_charge.keys())
		ncols = 1
		loc_dict_name = {'c': 'cytoplasm', 'p': 'periplasm'}

		for i in range(nrows):
			loc = locations[i]
			ax = plt.subplot(nrows, ncols, i+1)
			box = ax.get_position()
			ax.set_position([box.x0, box.y0, box.width, box.height])
			plt.plot(time, loc_charge[loc], color = colors[i])
			ax.legend(loc = "center left", bbox_to_anchor = (1, 0.5), labels = loc)
			ax.set_title("Metabolite charge in {}".format(loc_dict_name[loc]))
			ax.set_ylabel("Total metabolite charge")
			ax.set_xlabel("Time (min)")

		plt.tight_layout()
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close('all')


if __name__ == '__main__':
	Plot().cli()
