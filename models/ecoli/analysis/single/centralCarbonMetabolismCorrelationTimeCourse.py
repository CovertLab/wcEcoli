import os
import pickle

import numpy as np
from matplotlib import pyplot as plt

from wholecell.io.tablereader import TableReader
from wholecell.utils import units

from models.ecoli.processes.metabolism import COUNTS_UNITS, VOLUME_UNITS, TIME_UNITS
from wholecell.analysis.analysis_tools import exportFigure
from models.ecoli.analysis import singleAnalysisPlot

FLUX_UNITS = COUNTS_UNITS / VOLUME_UNITS / TIME_UNITS


class Plot(singleAnalysisPlot.SingleAnalysisPlot):
	_suppress_numpy_warnings = True

	def do_plot(self, simOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		validation_data = self.read_pickle_file(validationDataFile)
		sim_data = self.read_pickle_file(simDataFile)

		cellDensity = sim_data.constants.cell_density

		main_reader = TableReader(os.path.join(simOutDir, "Main"))
		initialTime = main_reader.readAttribute("initialTime")
		time = main_reader.readColumn("time") - initialTime

		massListener = TableReader(os.path.join(simOutDir, "Mass"))
		cellMass = massListener.readColumn("cellMass") * units.fg
		dryMass = massListener.readColumn("dryMass") * units.fg
		massListener.close()

		fbaResults = TableReader(os.path.join(simOutDir, "FBAResults"))
		reaction_ids = np.array(fbaResults.readAttribute("base_reaction_ids"))
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

		meanCorr = np.mean(np.array(corrCoefTimecourse)[~np.isnan(corrCoefTimecourse)])

		plt.figure()
		plt.plot(time / 60., corrCoefTimecourse)
		plt.axhline(y=meanCorr, color='r')
		plt.title("Measured vs. Simulated Central Carbon Fluxes")
		plt.text(.5*np.max(time / 60.),.95*meanCorr, "Mean = {:.2}".format(meanCorr), horizontalalignment="center")
		plt.xlabel("Time (min)")
		plt.ylabel("Pearson R")

		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close("all")


if __name__ == "__main__":
	Plot().cli()
