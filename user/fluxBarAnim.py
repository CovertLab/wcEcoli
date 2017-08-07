import numpy as np
#import matplotlib as mpl
#mpl.use('TKAgg')
from matplotlib import pyplot as plt
from matplotlib import animation

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


#get the flux data from cPickle
reactionFluxData = cPickle.load(open("reactionFluxData.cPickle","rb"))
reactionIDData = cPickle.load(open("reactionIDData.cPickle","rb"))
timeData = cPickle.load(open("time.cPickle","rb"))
# time = timeData[0]

#Build a mapping from reaction to color 
idToColor = {}
for reactionID, color in itertools.izip(reactionIDData.T, itertools.cycle(COLORS_LARGE)):
  idToColor[reactionID[0]] = color


#Initialize the figure and axes
rfSubset = reactionFluxData[0,:]
N = rfSubset.shape[0]
ind = np.arange(N)
width = 1
fig, ax = plt.subplots()
rectsText = ax.bar(ind, rfSubset, width, color = idToColor.values(), edgecolor = 'none')

time = 0

def init():
  #set parameters of the axes that won't change within each call of animate
  ax.set_xticks(ind + width/2)
  ax.tick_params(length = 0)
  ax.set_xticklabels(idToColor.keys(), fontsize = 4, rotation = 45, ha = 'right')
  ax.set_ylim(0,.8)
  return rectsText


def animate(i):
  global time 
  #for this time step, plot all reaction fluxes
  print time 
  if time > reactionFluxData.shape[0]:
    time = 0

  currRF = reactionFluxData[time,:]
  rectsText = ax.bar(ind, currRF, width, color = idToColor.values(), edgecolor = 'none')


  time += 50
  return rectsText,
   

anim = animation.FuncAnimation(fig, animate, 
                               init_func = init,  
                               frames=100, 
                               interval=100,
                               blit=True)

anim.save("fluxBarAnim" + ".mp4")

#anim.save("fluxBarAnim", writer="mencoder", fps=2)


#plt.show()




#make function that given a timestep makes png image of the plot (put timestep in title) if name = main for callable by command line

































