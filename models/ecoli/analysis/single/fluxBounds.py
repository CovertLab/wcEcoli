"""
Plot upper and lower flux targets

@author: Sophie Landon	
@organization: University of Bristol
@date: Created 14/5/2020
"""

from __future__ import absolute_import

import os
import cPickle

import numpy as np
from matplotlib import pyplot as plt

from wholecell.io.tablereader import TableReader
from wholecell.analysis.analysis_tools import exportFigure, read_bulk_molecule_counts
from models.ecoli.analysis import singleAnalysisPlot
from models.ecoli.processes.metabolism import COUNTS_UNITS, VOLUME_UNITS, TIME_UNITS, MASS_UNITS

class Plot(singleAnalysisPlot.SingleAnalysisPlot):
	def do_plot(self, simOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		if not os.path.isdir(simOutDir):
			raise Exception(
				'simOutDir "{}" does not currently exist as a directory'.format(simOutDir))

		if not os.path.exists(plotOutDir):
			os.mkdir(plotOutDir)

		sim_data = cPickle.load(open(simDataFile))

		#aaIDs = sim_data.moleculeGroups.aaIDs
		#(aaCounts,) = read_bulk_molecule_counts(simOutDir, (aaIDs,))


		mainListener = TableReader(os.path.join(simOutDir, "Main"))
                initialTime = mainListener.readAttribute("initialTime")
                time = mainListener.readColumn("time") - initialTime
                mainListener.close()

                massListener = TableReader(os.path.join(simOutDir, "Mass"))
                cellMass = massListener.readColumn("cellMass")
                dryMass = massListener.readColumn("dryMass")
                massListener.close()

                coefficient = dryMass / cellMass * sim_data.constants.cellDensity.asNumber(MASS_UNITS / VOLUME_UNITS)
		enzymeKineticsReader = TableReader(os.path.join(simOutDir, "EnzymeKinetics"))
                targetFluxesUpper = enzymeKineticsReader.readColumn('targetFluxesUpper')
		targetFluxesLower = enzymeKineticsReader.readColumn('targetFluxesLower')
		actualFluxes = enzymeKineticsReader.readColumn('actualFluxes')

		fig = plt.figure(figsize = (34, 34))
		ax_full = fig.add_subplot(111)
		ax_full.set_xlabel('Time', fontsize = 48)
		ax_full.set_ylabel('Flux', fontsize = 48)
		ax_full.tick_params(labelcolor='w', top=False, bottom=False, left=False, right=False)

		for idx in xrange(417):

			ax = fig.add_subplot(20, 21, idx + 1)

			plt.plot(time / 60., targetFluxesUpper[:, idx], linewidth = 2, label = 'Upper Target Flux')
			plt.plot(time / 60., targetFluxesLower[:, idx], linewidth = 2, label = 'Lower Target Flux')		
			plt.plot(time / 60., actualFluxes[:, idx], linewidth = 2, label = 'Actual Flux')
			plt.tick_params(labelsize=0)
		
		handles, labels = ax.get_legend_handles_labels()
		ax_full.legend(handles, labels, bbox_to_anchor = (0.35, 1.15), fontsize = 48)
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close("all")


if __name__ == "__main__":
	Plot().cli()
