import numpy as np
#import matplotlib as mpl
#mpl.use('TKAgg')
from matplotlib import pyplot as plt
from matplotlib.patches import Ellipse
from matplotlib.patches import Rectangle
from matplotlib.patches import BoxStyle
from matplotlib.patches import Arc
from matplotlib.patches import Wedge
import matplotlib.patches as mpatches
from matplotlib.collections import PatchCollection 
from matplotlib import animation
import cPickle

import time
start_time = time.time()

from multiprocessing import Pool

def Ecoli(x, y, width, height):
	"This will return a PatchCollection that looks like an Ecoli cell"

	rectHeight = height-width
	rectX, rectY = x-width/2, y-rectHeight/2
	radius = width/2

	topWedge = Wedge((x, y+rectHeight/2), r=radius, theta1=0.0, theta2=180.0)
	bottomWedge = Wedge((x,rectY), r=radius, theta1=180.0, theta2=0.0)
	rect = Rectangle((rectX,rectY), width=width, height=rectHeight)
	collection = PatchCollection([topWedge,bottomWedge,rect], color='g')

	return collection



def makeFrame(arg):
	parentHeight, daughterHeight, timestep = arg

	#Create new Figure and an Axes which fills it. 
	fig = plt.figure()
	fig.set_dpi(100)
	fig.set_size_inches(7, 6.5)
	ax = plt.axes(xlim=(0, 2.8), ylim=(0, 2.8))
	ylim = ax.get_ylim()

	#Make params of the new cells 
	x = 1.4
	pY = parentHeight/2.0
	dY = pY + daughterHeight - parentHeight 
	oldDY = daughterHeight + pY
	w = 1.0

	#Make new cells
	motherCell = Ecoli(x=x,y=pY, width=w, height=parentHeight) 
	daughterCell = Ecoli(x=x,y=dY, width=w, height=parentHeight)
	oldDaughterCell = Ecoli(x=x,y=oldDY, width=w, height=parentHeight)

	#Add new cells to axes
	ax.add_collection(motherCell)
	ax.add_collection(daughterCell)
	ax.add_collection(oldDaughterCell)

	plt.savefig('/scratch/users/mhinks/MotherCell/MoCellFrame' + str(timestep) + '.png', bbox_inches = 'tight')



def main():
	#import data
	heightData = cPickle.load(open("heightData.cPickle","rb"))
	timeData = cPickle.load(open("time.cPickle","rb"))
	#time = np.hstack((timeData[0], timeData[1]))
	args = []

	#for each cell, store its initial size 
	for i in range(0,2):
		currCell = heightData[0]
		ogHeight = currCell[0]
		time = timeData[i]
		for k in xrange(len(currCell)/100):
			args.append((ogHeight,currCell[k*100],time[k*100]))


	print "Initializing worker pool"
	pool = Pool(processes = 16)
	pool.map(makeFrame,args)

	pool.close()
	pool.join()
	



if __name__ == "__main__":
	main()















