import os
import pickle

import numpy as np
from matplotlib import pyplot as plt

from wholecell.io.tablereader import TableReader
from wholecell.utils import units
from models.ecoli.processes.metabolism import (COUNTS_UNITS, VOLUME_UNITS,
	TIME_UNITS)
from wholecell.analysis.analysis_tools import exportFigure
from models.ecoli.analysis import cohortAnalysisPlot

FLUX_UNITS = COUNTS_UNITS / VOLUME_UNITS / TIME_UNITS


class Plot(cohortAnalysisPlot.CohortAnalysisPlot):
	_suppress_numpy_warnings = True

	def do_plot(self, variantDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		plt.figure(figsize = (8.5, 11))

		# Get all cells in each seed

		validation_data = self.read_pickle_file(validationDataFile)
		sim_data = self.read_pickle_file(simDataFile)

		cellDensity = sim_data.constants.cell_density

		seed_color = {}
		line_instances = {}
		for seed_num in range(self.ap.n_seed):
			# Get all cells in this seed
			seedDir = self.ap.get_cells(seed=[seed_num])

			for simDir in seedDir:
				simOutDir = os.path.join(simDir, "simOut")

				time = TableReader(os.path.join(simOutDir, "Main")).readColumn("time")

				massListener = TableReader(os.path.join(simOutDir, "Mass"))
				cellMass = massListener.readColumn("cellMass") * units.fg
				dryMass = massListener.readColumn("dryMass") * units.fg
				massListener.close()

				fbaResults = TableReader(os.path.join(simOutDir, "FBAResults"))
				reaction_ids = np.array(
					fbaResults.readAttribute("base_reaction_ids"))
				reactionFluxes = (COUNTS_UNITS / VOLUME_UNITS / TIME_UNITS) * np.array(
					fbaResults.readColumn("base_reaction_fluxes"))
				rxn_id_to_index = {
					rxn_id: i for (i, rxn_id) in enumerate(reaction_ids)
				}

				dryMassFracAverage = np.mean(dryMass / cellMass)

				toya_reactions = validation_data.reactionFlux.toya2010fluxes["reactionID"]
				toya_fluxes = FLUX_UNITS * np.array([(dryMassFracAverage * cellDensity * x).asNumber(FLUX_UNITS) for x in validation_data.reactionFlux.toya2010fluxes["reactionFlux"]])

				netFluxes = []
				for toyaReactionID in toya_reactions:
					fluxTimeCourse = reactionFluxes[:, rxn_id_to_index[toyaReactionID]].asNumber(FLUX_UNITS).squeeze()
					netFluxes.append(fluxTimeCourse)

				trimmedReactions = FLUX_UNITS * np.array(netFluxes)

				corrCoefTimecourse = []
				for fluxes in trimmedReactions.asNumber(FLUX_UNITS).T:
					correlationCoefficient = np.corrcoef(fluxes, toya_fluxes.asNumber(FLUX_UNITS))[0,1]
					corrCoefTimecourse.append(correlationCoefficient)

				if seed_num in seed_color:
					current_line = plt.plot(time / 60., corrCoefTimecourse, color=seed_color[seed_num])
				else:
					current_line = plt.plot(time / 60., corrCoefTimecourse)
					seed_color[seed_num] = current_line[0].get_color()
					line_instances[seed_num] = current_line[0]

				plt.title("Measured vs. Simulated Central Carbon Fluxes")
				plt.xlabel("Time (min)")
				plt.ylabel("Pearson R")

		plt.legend(list(line_instances.values()), ["Seed {}".format(x) for x in line_instances.keys()], loc="best")

		plt.subplots_adjust(hspace = 0.2, wspace = 0.5)
		exportFigure(plt, plotOutDir, plotOutFileName,metadata)
		plt.close("all")


if __name__ == "__main__":
	Plot().cli()
