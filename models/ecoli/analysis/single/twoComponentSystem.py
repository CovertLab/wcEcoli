"""
Plot two component system counts

@author: Heejo Choi
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 5/20/2016
"""

import argparse
import os
import cPickle

import numpy as np
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt

from wholecell.io.tablereader import TableReader
import wholecell.utils.constants
from wholecell.utils import units

def main(simOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata = None):

	if not os.path.isdir(simOutDir):
		raise Exception, "simOutDir does not currently exist as a directory"

	if not os.path.exists(plotOutDir):
		os.mkdir(plotOutDir)

	sim_data = cPickle.load(open(simDataFile, "rb"))
	TCS_IDS = []
	moleculeTypeOrder = ["HK", "PHOSPHO-HK", "LIGAND", "HK-LIGAND", "PHOSPHO-HK-LIGAND", "RR", "PHOSPHO-RR"]
	moleculeTypeLocation = ["[i]", "[i]", "[p]", "[i]", "[i]", "[c]", "[c]"]
	moleculeTypeColor = ["b", "b", "orange", "g", "g", "r", "r"]
	for system in sim_data.moleculeGroups.twoComponentSystems:
		for idx, moleculeType in enumerate(moleculeTypeOrder):
			if moleculeType == "LIGAND" and system["molecules"][moleculeType] == "PHOSPHO-ABC-27-CPLX":
				TCS_IDS.append(str(system["molecules"][moleculeType]) + "[i]")
			else:
				TCS_IDS.append(str(system["molecules"][moleculeType]) + moleculeTypeLocation[idx])

	bulkMolecules = TableReader(os.path.join(simOutDir, "BulkMolecules"))
	moleculeIds = bulkMolecules.readAttribute("objectNames")
	moleculeIndexes = np.array([moleculeIds.index(moleculeId) for moleculeId in TCS_IDS], np.int)
	moleculeCounts = bulkMolecules.readColumn("counts")[:, moleculeIndexes]
	initialTime = TableReader(os.path.join(simOutDir, "Main")).readAttribute("initialTime")
	time = TableReader(os.path.join(simOutDir, "Main")).readColumn("time") - initialTime

	bulkMolecules.close()

	# Convert molecule counts to concentrations
	nAvogadro = sim_data.constants.nAvogadro.asNumber(1 / units.mol)
	cellDensity = sim_data.constants.cellDensity.asNumber(units.g / units.L)
	mass = TableReader(os.path.join(simOutDir, "Mass"))
	cellMass = (units.fg * mass.readColumn("cellMass")).asNumber(units.g)
	cellVolume = cellMass / cellDensity

	rows = 3
	cols = 5
	num_subentries = 7

	plt.figure(figsize = (8.5, 11))

	for idx in xrange(len(sim_data.moleculeGroups.twoComponentSystems)):
		grid_loc = idx + 1 + (cols*(num_subentries + 1))*( idx / cols)

		for subentryIdx in xrange(len(moleculeTypeOrder)):
			ax = plt.subplot(rows*(num_subentries + 2), cols, grid_loc+((cols)*(subentryIdx)))
			ax.plot(time / 60., moleculeCounts[:, (idx * num_subentries) + subentryIdx] / (cellVolume * nAvogadro), linewidth = 1, color = moleculeTypeColor[subentryIdx])
			ax.set_title(TCS_IDS[(idx * num_subentries) + subentryIdx], fontsize = 4)

			ymin = np.amin(moleculeCounts[:, (idx * num_subentries) + subentryIdx] / (cellVolume * nAvogadro))
			ymax = np.amax(moleculeCounts[:, (idx * num_subentries) + subentryIdx] / (cellVolume * nAvogadro))
			ax.set_ylim([ymin, ymax])
			ax.set_yticks([ymin, ymax])
			ax.set_yticklabels(["%0.2e" % ymin, "%0.2e" % ymax])
			ax.spines['top'].set_visible(False)
			ax.spines['bottom'].set_visible(False)
			ax.xaxis.set_ticks_position('none')
			ax.tick_params(which = 'both', direction = 'out', labelsize = 4)
			ax.set_xticks([])

		
	plt.subplots_adjust(hspace = 1, wspace = 1)

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
