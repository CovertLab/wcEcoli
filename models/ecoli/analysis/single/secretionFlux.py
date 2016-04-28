#!/usr/bin/env python

import argparse
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt

from wholecell.io.tablereader import TableReader
import wholecell.utils.constants
from wholecell.utils import units
from reconstruction.ecoli.knowledge_base_raw import KnowledgeBaseEcoli

SECRETED = ["ACET[p]", "ETOH[p]", "D-LACTATE[p]", "CARBON-DIOXIDE[p]"]

def main(simOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata=None):
	if not os.path.isdir(simOutDir):
		raise Exception, "simOutDir does not currently exist as a directory"

	if not os.path.exists(plotOutDir):
		os.mkdir(plotOutDir)

	fig, axes = plt.subplots(1, sharex = True)

	# Plot flux of secreted molecules that are not 0
	initialTime = TableReader(os.path.join(simOutDir, "Main")).readAttribute("initialTime")
	time = TableReader(os.path.join(simOutDir, "Main")).readColumn("time") - initialTime

	fba_results = TableReader(os.path.join(simOutDir, "FBAResults"))
	exFlux = fba_results.readColumn("externalExchangeFluxes")
	exMolec = fba_results.readAttribute("externalMoleculeIDs")
	
	raw_data = KnowledgeBaseEcoli()
	for secreted in raw_data.secretions:
		nutrient = secreted.get("molecule id").encode('ascii')
		if nutrient in exMolec:
			flux = exFlux[:,exMolec.index(nutrient)]
			if np.max(np.abs(flux)) > 0.01:
				axes.plot(time / 60. / 60., -1. * flux, label='%s avg:%5.2f' % (nutrient,np.mean(flux)))
	
	plt.legend(prop={'size':4})
	
	axes.set_ylabel("Secreted nutrient\n(mmol/gDCW/hr)", fontsize=10)
	axes.set_xlabel("Time(hr)", fontsize=10)

	from wholecell.analysis.analysis_tools import exportFigure
	exportFigure(plt, plotOutDir, plotOutFileName, metadata)
	plt.close("all")



if __name__ == "__main__":
	defaultSimDataFile = os.path.join(
			wholecell.utils.constants.SERIALIZED_KB_DIR,
			wholecell.utils.constants.SERIALIZED_KB_MOST_FIT_FILENAME
			)

	parser = argparse.ArgumentParser()
	parser.add_argument("simOutDir", help = "Directory containing simulation output", type = str)
	parser.add_argument("plotOutDir", help = "Directory containing plot output (will get created if necessary)", type = str)
	parser.add_argument("plotOutFileName", help = "File name to produce", type = str)
	parser.add_argument("--simDataFile", help = "KB file name", type = str, default = defaultSimDataFile)

	args = parser.parse_args().__dict__

	main(args["simOutDir"], args["plotOutDir"], args["plotOutFileName"], args["simDataFile"])
