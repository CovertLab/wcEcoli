#!/usr/bin/env python
"""
Plot reaction max rate over course of the simulation.

@date: Created 7/02/2015
@author: Morgan Paull
@organization: Covert Lab, Department of Bioengineering, Stanford University
"""

from __future__ import division

import argparse
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")

from matplotlib import pyplot as plt
import matplotlib.animation as animation

from wholecell.io.tablereader import TableReader
import wholecell.utils.constants

DISABLED = True

def main(simOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata = None):
	
	if DISABLED:
		print "Currently disabled."
		return

	if not os.path.isdir(simOutDir):
		raise Exception, "simOutDir does not currently exist as a directory"

	if not os.path.exists(plotOutDir):
		os.mkdir(plotOutDir)

	enzymeKineticsdata = TableReader(os.path.join(simOutDir, "EnzymeKinetics"))
	
	enzymeKineticsArray = enzymeKineticsdata.readColumn("reactionConstraints")
	overconstraintMultiples = enzymeKineticsdata.readColumn("overconstraintMultiples")

	reactionIDs = enzymeKineticsdata.readAttribute("reactionIDs")
	
	initialTime = TableReader(os.path.join(simOutDir, "Main")).readAttribute("initialTime")
	time = TableReader(os.path.join(simOutDir, "Main")).readColumn("time") - initialTime
	
	enzymeKineticsdata.close()


	fbaData = TableReader(os.path.join(simOutDir, "FBAResults"))

	reactionRatesUnconstrained = fbaData.readColumn('reactionFluxes')

	rateEstimatesArray = enzymeKineticsArray.copy()

	unconstrainedFluxes = fbaData.readColumn('reactionFluxes')
	fluxNames = fbaData.readAttribute('reactionIDs')
	simulationSteps = fbaData.readColumn('simulationStep')
	fbaData.close()



	# Only look at fluxes which had a noninfinite estimate at at least one point
	rateEstimates =  rateEstimatesArray[1:]
	rateEstimates += 1.0
	rateEstimates[np.where(rateEstimates == np.inf)] = 0

	fluxesWithEstimates = unconstrainedFluxes[:,np.where(np.sum(rateEstimates, axis=0) > 0)[0]]
	fluxNamesEstimates = np.array([fluxNames[x] for x in np.where(np.sum(rateEstimates, axis=0) > 0)[0]])

	rateEstimates = rateEstimates[:,np.where(np.sum(rateEstimates, axis=0) > 0)[0]]
	rateEstimates -= 1.0
	rateEstimates = np.vstack((np.zeros(rateEstimates.shape[1]),rateEstimates))



	
	relativeRates = fluxesWithEstimates / rateEstimates
	# Remove NaNs from the arrays
	# relativeRates[np.where(np.isnan(relativeRates))] = 0
	# Only plot rates which are overconstrained
	relativeRates[np.where(relativeRates < 1)] = 0
	relativeRates[np.isnan(relativeRates)] = 0

	num_x = relativeRates.shape[1]

	# First set up the figure, the axis, and the plot element we want to animate
	fig, ax = plt.subplots(1)

	# import ipdb; ipdb.set_trace()

	ax.set_xlim(0, num_x)
	ax.set_ylim(0, 100)
	plt.title("Kinetic Rates Divided By Reaction Fluxes")
	plt.xlabel("Reaction")
	plt.ylabel("Fold Difference")


	rects = ax.bar(range(1, num_x+1), relativeRates[0,:],  align='center', color='purple', alpha=.7, label="Fold Difference")
	plt.legend(framealpha=.5)
	plt.xticks(range(1, num_x+1), fluxNamesEstimates, rotation='vertical', fontsize=5)

	# initialization function: plot the background of each frame
	def init():
		for idx, patch in enumerate(rects.patches):
			patch.set_height(relativeRates[0,idx])
		return rects,

	# animation function.  This is called sequentially
	def animate(i):
		if i % 10 == 0:
			print "Frame %i" % (i)
		for idx, patch in enumerate(rects.patches):
			# print "relativeRates[%i, %i] is: " % (idx, i)
			# print relativeRates[idx, i]
			patch.set_height(relativeRates[i, idx])
		return rects,

	# call the animator.  blit=True means only re-draw the parts that have changed.
	anim = animation.FuncAnimation(fig, animate, init_func=init,
		frames=relativeRates.shape[0], interval=20, blit=True)

	# save the animation as an mp4.  This requires ffmpeg or mencoder to be
	# installed.  The extra_args ensure that the x264 codec is used, so that
	# the video can be embedded in html5.  You may need to adjust this for
	# your system: for more information, see
	# http://matplotlib.sourceforge.net/api/animation_api.html
	anim.save(os.path.join(plotOutDir, plotOutFileName) + '.mp4', fps=15, extra_args=['-vcodec', 'libx264'])




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
