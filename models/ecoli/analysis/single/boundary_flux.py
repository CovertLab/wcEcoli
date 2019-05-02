"""
Plot boundary fluxes
@organization: Covert Lab, Department of Bioengineering, Stanford University
"""

from __future__ import absolute_import
from __future__ import division

import cPickle
import os

from matplotlib import pyplot as plt
import numpy as np

from models.ecoli.analysis import singleAnalysisPlot
from wholecell.analysis.analysis_tools import exportFigure
from wholecell.io.tablereader import TableReader
from wholecell.utils import units
from wholecell.utils import filepath

from models.ecoli.processes.metabolism import COUNTS_UNITS, VOLUME_UNITS, TIME_UNITS, MASS_UNITS

BURN_IN_STEPS = 20 # remove initialization artifacts

def set_ticks(ax, time):
	ax.spines['top'].set_visible(False)
	ax.spines['bottom'].set_visible(False)
	ax.xaxis.set_ticks_position('none')
	ax.tick_params(which='both', direction='out', labelsize=10)
	ax.set_xticks([time.min() / 60., time.max() / 60.])

class Plot(singleAnalysisPlot.SingleAnalysisPlot):
	def do_plot(self, simOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		if not os.path.isdir(simOutDir):
			raise Exception, 'simOutDir does not currently exist as a directory'

		filepath.makedirs(plotOutDir)

		with open(simDataFile, 'rb') as f:
			sim_data = cPickle.load(f)
		with open(validationDataFile, 'rb') as f:
			validation_data = cPickle.load(f)

		# Listeners used
		main_reader = TableReader(os.path.join(simOutDir, 'Main'))
		mass_reader = TableReader(os.path.join(simOutDir, "Mass"))
		enzyme_kinetics_reader = TableReader(os.path.join(simOutDir, "EnzymeKinetics"))

		# main reader
		initialTime = main_reader.readAttribute("initialTime")
		time = main_reader.readColumn("time") - initialTime
		main_reader.close()

		# mass reader
		cellMass = mass_reader.readColumn("cellMass")
		dryMass = mass_reader.readColumn("dryMass")
		mass_reader.close()

		coefficient = dryMass / cellMass * sim_data.constants.cellDensity.asNumber(MASS_UNITS / VOLUME_UNITS)

		# enzyme kinetics reader
		allTargetFluxes = (COUNTS_UNITS / MASS_UNITS / TIME_UNITS) * (enzyme_kinetics_reader.readColumn("targetFluxes").T / coefficient).T
		allActualFluxes = (COUNTS_UNITS / MASS_UNITS / TIME_UNITS) * (enzyme_kinetics_reader.readColumn("actualFluxes").T / coefficient).T
		kineticsConstrainedReactions = np.array(enzyme_kinetics_reader.readAttribute("kineticsConstrainedReactions"))
		boundaryConstrainedReactions = np.array(enzyme_kinetics_reader.readAttribute("boundaryConstrainedReactions"))
		enzyme_kinetics_reader.close()

		allTargetFluxes = allTargetFluxes.asNumber(units.mmol / units.g / units.h)
		allActualFluxes = allActualFluxes.asNumber(units.mmol / units.g / units.h)

		# boundary target fluxes
		boundaryTargetFluxes = allTargetFluxes[:, len(kineticsConstrainedReactions):]
		boundaryActualFluxes = allActualFluxes[:, len(kineticsConstrainedReactions):]

		boundary_target_flux_dict = dict(zip(boundaryConstrainedReactions, boundaryTargetFluxes.T))
		boundary_actual_flux_dict = dict(zip(boundaryConstrainedReactions, boundaryActualFluxes.T))


		## Plot
		cols = 1
		rows = len(boundaryConstrainedReactions) + 1
		n_char_of_reaction_id = 25 # number of characters of reaction_id string in title

		# initialize subplot indices
		row_idx = 1

		plt.figure(figsize=(5*cols, 2*rows))
		for reaction_id, reaction_flux in boundary_actual_flux_dict.iteritems():
			target_flux = boundary_target_flux_dict[reaction_id]

			# initialize flux column
			col = 1
			plot_index = row_idx * cols + col
			ax1 = plt.subplot(rows, cols, plot_index)

			# plot flux
			ax1.plot(time / 60., target_flux, color='Red', label='target flux')
			ax1.plot(time / 60., reaction_flux, color='Blue', label='actual flux')

			# add labels
			ax1.set_xlabel("Time (min)", fontsize = 12)
			ax1.set_ylabel("flux", fontsize = 12)
			ax1.set_title("%i. %s" % (row_idx, reaction_id[:n_char_of_reaction_id]), fontsize=14, y=1.15)
			set_ticks(ax1, time)
			ax1.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=8)

			row_idx += 1


		plt.tight_layout()
		plt.subplots_adjust(hspace=1.5, wspace=2.0)
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close('all')


if __name__ == '__main__':
	Plot().cli()
