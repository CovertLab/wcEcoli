import numpy as np
#import matplotlib as mpl
#mpl.use('TKAgg')
from matplotlib import pyplot as plt
from matplotlib.patches import Ellipse
from matplotlib.patches import BoxStyle
import matplotlib.patches as mpatches
from matplotlib.collections import PatchCollection 
from matplotlib import animation
import cPickle

#get the height data from cPickle
heightData = cPickle.load(open("heightData.cPickle","rb"))

#Create new Figure and an Axes which fills it. 
fig = plt.figure()
fig.set_dpi(100)
fig.set_size_inches(7, 6.5)
ax = plt.axes(xlim=(0, 2.8), ylim=(0, 2.8))
ylim = ax.get_ylim()

ax2 = ax

#time series of heights of mother cell - compute this from the volumes so don't have to mess with downstream code
heightTS = heightData #np.array([[np.linspace(1, np.random.uniform(1,2.4), 30)] for i in range(0,10)])
currTime = 0
currCell = 0

currHeightTS = heightTS[currCell]

#Create cell data (maybe in future have cell class which is two arcs and a rectangle)
initHeight = currHeightTS[0]
initialY =initHeight/2
initialX = 1.4

motherCell = Ellipse((initialX,initialY), width=1, height=initHeight, angle=0.0) 
daughterCell = Ellipse((initialX,10), width=1, height=initHeight, angle=0.0)
#import ipdb; ipdb.set_trace()
ax.add_patch(motherCell)
ax.add_patch(daughterCell)
patches = [motherCell, daughterCell]
divided = False
ax2.add_patch(motherCell)

#Ecoli = mpatches.FancyBboxPatch((initialX,initialY), width=.25, height=initHeight, boxstyle = 
#BoxStyle("Round",pad=0.3, rounding_size=0.45))
#ax.add_patch(Ecoli)

#TODO
#calculate height properly so this looks normal

def init():
  return motherCell, daughterCell


def animate(i):
  global currTime
  global currCell

  
  currHeightTS = heightTS[currCell]

  #if we aren't at the end of the heightTS
  if currTime < currHeightTS.size:
    #get old height and center
    oldHeight = motherCell.height
    x, y = motherCell.center
    #update the height 
    motherCell.height = currHeightTS[currTime]
    #
    #print motherCell.height
    #update the center
    y += (motherCell.height-oldHeight)/2
    motherCell.center = (x,y)
    #increase the index by 1
    currTime += 100 
    if divided == True:
      daughterCell.height=motherCell.height
      x2, y2 = daughterCell.center
      y2 += 3*(motherCell.height-oldHeight)/2
      daughterCell.center = (x2, y2)

  else:
    #reset height index
    currTime = 0  
    #either move to next cell, or start over 
    if currCell < len(heightTS)-1:
      currCell += 1
    else:
      currCell = 0

    #import ipdb; ipdb.set_trace()
    #reset motherCell height and position based on on next TS
    motherCell.height = heightTS[currCell][currTime]
    motherCell.center = (initialX,motherCell.height/2)

    global divided
    divided = True
    x, y = motherCell.center
    daughterCell.height = motherCell.height
    daughterCell.center = (x, y+motherCell.height)
    #import ipdb; ipdb.set_trace()
    
  

  return motherCell, daughterCell, divided
   

anim = animation.FuncAnimation(fig, animate, 
                               init_func = init,  
                               frames=360, 
                               interval=100,
                               blit=False)

#plt.show()
anim.save("Mothercell" + ".mp4")










































