#!/usr/bin/env python
'''
Plots environment nutrient concentrations

@organization: Covert Lab, Department of Bioengineering, Stanford University
'''

from __future__ import absolute_import

import os
import cPickle
import math

import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import itertools

from wholecell.io.tablereader import TableReader

from models.ecoli.processes.metabolism import COUNTS_UNITS, TIME_UNITS, VOLUME_UNITS
from wholecell.analysis.plotting_tools import COLORS_LARGE
from wholecell.analysis.analysis_tools import exportFigure
from models.ecoli.analysis import singleAnalysisPlot

FLUX_UNITS = COUNTS_UNITS / VOLUME_UNITS / TIME_UNITS

BURN_IN_PERIOD = 150
RANGE_THRESHOLD = 2
MOVING_AVE_WINDOW_SIZE = 200

class Plot(singleAnalysisPlot.SingleAnalysisPlot):
	def do_plot(self, simOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		if not os.path.isdir(simOutDir):
			raise Exception, 'simOutDir does not currently exist as a directory'

		if not os.path.exists(plotOutDir):
			os.mkdir(plotOutDir)

		# Load data from KB
		sim_data = cPickle.load(open(simDataFile, 'rb'))

		# Load masses
		mass = TableReader(os.path.join(simOutDir, 'Mass'))
		protein = mass.readColumn('proteinMass')
		tRna = mass.readColumn('tRnaMass')
		rRna = mass.readColumn('rRnaMass')
		mRna = mass.readColumn('mRnaMass')
		dna = mass.readColumn('dnaMass')
		smallMolecules = mass.readColumn('smallMoleculeMass')
		total_dry_mass = mass.readColumn('dryMass')

		mass_fractions = np.vstack([
			protein / protein[0],
			rRna / rRna[0],
			tRna / tRna[0],
			mRna / mRna[0],
			dna / dna[0],
			smallMolecules / smallMolecules[0],
			]).T
		massLabels = ['Protein', 'rRNA', 'tRNA', 'mRNA', 'DNA', 'Small Mol.s']

		# Load time
		initialTime = TableReader(
			os.path.join(simOutDir, 'Main')).readAttribute('initialTime')
		time = TableReader(os.path.join(simOutDir, 'Main')).readColumn(
			'time') - initialTime

		# Load environment data
		environment = TableReader(os.path.join(simOutDir, 'Environment'))
		nutrient_names = environment.readAttribute('objectNames')
		nutrient_concentrations = environment.readColumn(
			'nutrientConcentrations')
		# volume = environment.readColumn('volume')
		environment.close()

		# Load flux data
		fbaResults = TableReader(os.path.join(simOutDir, 'FBAResults'))
		reactionIDs = np.array(fbaResults.readAttribute('reactionIDs'))
		# reactionFluxes = np.array(fbaResults.readColumn('reactionFluxes'))
		externalExchangeFluxes = np.array(fbaResults.readColumn('externalExchangeFluxes'))
		fbaResults.close()


		# Build a mapping from nutrient_name to color
		idToColor = {}
		for nutrient_name, color in itertools.izip(nutrient_names, itertools.cycle(COLORS_LARGE)):
			idToColor[nutrient_name] = color

		# Build a mapping from reactionID to color
		rxnIdToColor = {}
		for reactionID, color in itertools.izip(reactionIDs, itertools.cycle(COLORS_LARGE)):
			rxnIdToColor[reactionID] = color

		fig = plt.figure(figsize=(30, 20))

		ax1_1 = fig.add_subplot(5, 2, 1)
		ax1_2 = fig.add_subplot(5, 2, 3)
		ax1_3 = fig.add_subplot(5, 2, 5)
		ax2 = fig.add_subplot(2, 2, 2)

		## cell mass
		# total mass
		ax1_1.plot(time, total_dry_mass, linewidth=2)
		ax1_1.ticklabel_format(useOffset=False)
		ax1_1.set_xlabel('Time (sec)')
		ax1_1.set_ylabel('Mass (fg)')
		ax1_1.set_title('Dry Cell Mass')

		# mass fractions
		ax1_2.plot(time, mass_fractions, linewidth=2)
		ax1_2.set_xlabel('Time (sec)')
		ax1_2.set_ylabel('Mass (normalized by t = 0)')
		ax1_2.set_title('Mass Fractions of Biomass Components')
		ax1_2.legend(massLabels, loc='best')

		# exchange fluxes
		for idx, (reactionID, externalExchangeFlux) in enumerate(zip(reactionIDs, externalExchangeFluxes.T)):
			runningMeanFlux = np.convolve(externalExchangeFlux[BURN_IN_PERIOD:], np.ones((MOVING_AVE_WINDOW_SIZE,))/MOVING_AVE_WINDOW_SIZE, mode='valid')
			meanNormFlux = runningMeanFlux / np.mean(runningMeanFlux)
			fluxRange = meanNormFlux.max() - meanNormFlux.min()

			if fluxRange > RANGE_THRESHOLD:
				ax1_3.plot(time, externalExchangeFlux, linewidth=2, label=reactionID, color=rxnIdToColor[reactionID])

		ax1_3.set_yscale('symlog')
		ax1_3.set_xlabel('Time (sec)')
		ax1_3.set_ylabel('symlog Flux {}'.format(FLUX_UNITS.strUnit()))
		ax1_3.set_title('Exchange Fluxes')

		# place legend
		ax1_3.legend(bbox_to_anchor=(0.5, -0.25), loc=9, borderaxespad=0.,
			ncol=1, prop={'size': 10})

		## environment concentrations
		for idx, nutrient_name in enumerate(nutrient_names):
			if (not math.isnan(nutrient_concentrations[0, idx]) and np.mean(
					nutrient_concentrations[:, idx]) != 0):
				ax2.plot(time, nutrient_concentrations[:, idx], linewidth=2,
					label=nutrient_name, color=idToColor[nutrient_name])

		ax2.set_yscale('symlog',linthresh=1, linscale=0.1)
		ax2.set_title('environment concentrations -- symlog')
		ax2.set_xlabel('Time (sec)')
		ax2.set_ylabel('symlog concentration (mmol/L)')

		# place legend
		ax2.legend(bbox_to_anchor=(0.5, -0.25), loc=9, borderaxespad=0.,
			ncol=3, prop={'size': 10})

		plt.subplots_adjust(wspace=0.2, hspace=0.4)

		exportFigure(plt, plotOutDir, plotOutFileName,metadata)
		plt.close('all')


if __name__ == '__main__':
	Plot().cli()