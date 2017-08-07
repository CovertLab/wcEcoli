#!/usr/bin/env python

import argparse
import os
import re

import numpy as np
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt
import matplotlib.patches as patches


from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from wholecell.io.tablereader import TableReader
import wholecell.utils.constants
from wholecell.utils import units
from wholecell.utils.sparkline import whitePadSparklineAxis

# extra imports from allReactionFluxes.py
import cPickle
import itertools

from wholecell.io.tablereader import TableReader
import wholecell.utils.constants

from models.ecoli.processes.metabolism import COUNTS_UNITS, TIME_UNITS, VOLUME_UNITS
FLUX_UNITS = COUNTS_UNITS / VOLUME_UNITS / TIME_UNITS

from wholecell.analysis.plotting_tools import COLORS_LARGE

BURN_IN_PERIOD = 150

NUMERICAL_ZERO = 1e-15

RANGE_THRESHOLD = 2
MOVING_AVE_WINDOW_SIZE = 200

PLACE_HOLDER = -1
CELL_WIDTH = 1.0
FONT_SIZE=8
trim = 0.03

variantDir = "out/20170723.031448.241997__Daily_build./wildtype_000000"
simDataFile = "out/20170723.031448.241997__Daily_build./wildtype_000000/kb/simData_Modified.cPickle"

import cPickle
sim_data = cPickle.load(open(simDataFile, "rb"))


#Get cell density
cellDensity = sim_data.constants.cellDensity.asNumber(units.fg/units.L)

# Get all cells in each seed
ap = AnalysisPaths(variantDir, cohort_plot = True)
all_cells = ap.get_cells(generation=[0,1])

heightData = []
reactionFluxes = np.array([])
reactionIDs = []
time = []
indices = []
for idx, simDir in enumerate(all_cells):

	simOutDir = os.path.join(simDir, "simOut")

	timeThisCell = TableReader(os.path.join(simOutDir, "Main")).readColumn("time")
	time.append(timeThisCell)




	## Cell mass
	mass = TableReader(os.path.join(simOutDir, "Mass"))
	cellMass = mass.readColumn("cellMass")
	
	cellVolume = cellMass/cellDensity*10**15

	#micrometer
	cellHeight = cellVolume/(np.pi*(CELL_WIDTH/2)**2)
	#import ipdb; ipdb.set_trace()
	heightData.append(cellHeight)
	#import ipdb; ipdb.set_trace()

	#reaction IDs and fluxes 
	fbaResults = TableReader(os.path.join(simOutDir, "FBAResults"))
	#get IDs on the first one 
	if idx == 0:		
		reactionIDsThisCell = np.array(fbaResults.readAttribute("reactionIDs"))
		reactionIDs.append(reactionIDsThisCell) 

	#get reactionFluxes for each cell, append horizontally to reactionFluxes ndarray
	reactionFluxesThisCell = np.array(fbaResults.readColumn("reactionFluxes"))
	reactionFluxesThisCell[np.abs(reactionFluxesThisCell) < NUMERICAL_ZERO] = 0

	if idx == 0:
		reactionFluxes = reactionFluxesThisCell
	else:
		reactionFluxes = np.vstack((reactionFluxes, reactionFluxesThisCell))

	
	fbaResults.close()

	#check if, within this generation, any flux changes dramatically. Add rxn index to indices if so. 
	for idx, (reactionID, reactionFlux) in enumerate(zip(reactionIDsThisCell, reactionFluxesThisCell.T)):
		runningMeanFlux = np.convolve(reactionFlux[BURN_IN_PERIOD:], np.ones((MOVING_AVE_WINDOW_SIZE,))/MOVING_AVE_WINDOW_SIZE, mode='valid')

		meanNormFlux = runningMeanFlux / np.mean(runningMeanFlux)
		fluxRange = meanNormFlux.max() - meanNormFlux.min()

		if fluxRange > RANGE_THRESHOLD:
			if idx not in indices:
				indices.append(idx)

	

heightData = np.array(heightData)
reactionIDs = np.array(reactionIDs)
reactionFluxes = np.array(reactionFluxes)
indices = np.array(indices)
time = np.array(time)
#this is saved in whatever directory this script is run in, so wcEcoli in this case
import ipdb; ipdb.set_trace()
cPickle.dump(heightData, open("heightData.cPickle", "wb"))
cPickle.dump(reactionIDs[:,indices], open("reactionIDData.cPickle", "wb"))
cPickle.dump(reactionFluxes[:,indices], open("reactionFluxData.cPickle", "wb"))
cPickle.dump(time, open("time.cPickle", "wb"))
