import numpy as np
#import matplotlib as mpl
#mpl.use('TKAgg')
from matplotlib import pyplot as plt
from matplotlib import animation
from multiprocessing import Pool
import time
start_time = time.time()
#imports from allReactionFluxes.py
import argparse
import os
import cPickle
import matplotlib
matplotlib.use("Agg")
import itertools

from wholecell.io.tablereader import TableReader
import wholecell.utils.constants

from models.ecoli.processes.metabolism import COUNTS_UNITS, TIME_UNITS, VOLUME_UNITS
FLUX_UNITS = COUNTS_UNITS / VOLUME_UNITS / TIME_UNITS

from wholecell.analysis.plotting_tools import COLORS_LARGE
from wholecell.analysis.analysis_tools import exportFigure


def fluxBar(arg):
	timestep, fluxes, idData = arg
	#Build a mapping from reaction to color 
	idToColor = {}
	for reactionID, color in itertools.izip(idData.T, itertools.cycle(COLORS_LARGE)):
	  idToColor[reactionID[0]] = color

	N = fluxes.shape[0]
	ind = np.arange(N)
	fig, ax = plt.subplots()
	width = 1
	rectsText = ax.bar(ind, fluxes, width, color = idToColor.values(), edgecolor = 'none')
	ax.set_xticks(ind + width/2)
  	ax.tick_params(length = 0)
 	ax.set_xticklabels(idToColor.keys(), fontsize = 4, rotation = 45, ha = 'right')
  	ax.set_ylim(0,.8)
  	ax.set_xlabel('Reaction')
  	ax.set_ylabel('Flux [UNITS?]')
  	ax.set_title('Fluxes at timestep ' + str(timestep))
	plt.savefig('/scratch/users/mhinks/fluxesTS' + str(timestep) + '.png', bbox_inches = 'tight')
	plt.close()
	print "Saved time %f" % timestep


def main():
	
	#import data
	reactionFluxData = cPickle.load(open("reactionFluxData.cPickle","rb"))
	reactionIDData = cPickle.load(open("reactionIDData.cPickle","rb"))
	timeData = cPickle.load(open("time.cPickle","rb"))
	time = np.hstack((timeData[0], timeData[1]))
	args = []
	for i in xrange(len(time)):
		args.append((time[i],reactionFluxData[i,:],reactionIDData))

	print "Initializing worker pool"
	pool = Pool(processes = 16)
	pool.map(fluxBar,args)

	pool.close()
	pool.join()



#fluxBar takes in a timestep and the array of rxn fluxes for that time step, 
#creates a bar chart and saves it as a .png 



#fluxBar((tsTest, fluxesTest))

if __name__ == "__main__":
	main()