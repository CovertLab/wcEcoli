'''
Creates series of svg files that can be stitched together using avconv to create an animation of Ecoli cells dividing like a "mother machine"
Input: growthAnimation.tsv file in /scratch/PI/mcovert/thorst/allenTalk with 3 columns (current time, mass of cell, final timestep in a cell cycle)
Output: time series of svg images to specified folder
'''

import numpy as np
import os
import cPickle
import csv 

def boilerplate_start(h):
	# h.write("<html>\n")
	# h.write("<body>\n")
	h.write("<svg width=\"120\" height=\"400\"  viewBox=\"0 0 120 500\" xmlns=\"http://www.w3.org/2000/svg\" version=\"1.1\" xmlns:xlink=\"http://www.w3.org/1999/xlink\">\n")
	h.write(" <g>\n")
	h.write("  <rect x=\"0\" y=\"0\" width=\"120\" height=\"500\" fill = \"white\" fill-opacity=\"1.0\"/> \n")
	h.write(" </g>")

def boilerplate_end(h):
	h.write("</svg>\n")
	# h.write("</body>\n")
	# h.write("</html>")

def svgEcoliParams(height):
	#compute parameters needed to write an svg file containing a cell of the appropriate size
	hW = 50
	height *= 100
	innerRect = height - 2*hW
	tR = (hW,innerRect+hW)
	tL = (-hW,innerRect+hW)
	bL = (-hW,hW)
	bR = (hW,hW)
	svgParams = [tR, tL, bL, bR]
	return svgParams

def writeEcoliSVG(h,params):
	#insert all the h.write lines that would make this cell, but don't actually display it yet
	tR, tL, bL, bR = params
	h.write(" <defs>\n")
	h.write("  <symbol id=\"Ecoli\">\n")
	h.write("   <g transform=\"translate(60,500) scale(1,-1)\">\n")
	h.write("    <g style=\"stroke: black; stroke-width: 0\">\n")
	h.write("     <path d=\"M %f %f\n" % (tR[0], tR[1]))
	h.write("        A 50 50, 0, 0, 1, %f %f\n" % (tL[0], tL[1]))							
	h.write("        L %f %f\n" % (bL[0], bL[1]))							
	h.write("        A 50 50, 0, 0, 1, %f %f\n" % (bR[0], bR[1]))							
	h.write("        Z")							
	h.write("     \" fill-opacity=\"1.0\" />\n")							
	h.write("    </g>\n")
	h.write("   </g>\n")
	h.write("  </symbol>\n")
	h.write(" </defs>\n")

def writeSVG(dirname, currHeight, finalCellsHeight, pinching, timestep):
	#put starting boilerplate (1 line)
	if not os.path.exists(dirname):
		os.makedirs(dirname)

	h = open(os.path.join(dirname, "moCell_%05d.svg" % timestep), "w")
	boilerplate_start(h)

	#figure out if you're going to simulate pinching or not, 
	#figure out positions and heights of cells according to that check
	if pinching == False:
		#figure out params of and write three identically sized cells that will grow and move up
		svgParams = svgEcoliParams(currHeight)
		writeEcoliSVG(h,svgParams)
		h.write("<use xlink:href=\"#Ecoli\" x=\"0\" y=\"0\" />\n")
		h.write("<use xlink:href=\"#Ecoli\" x=\"0\" y=\"-%f\" />\n" % (currHeight*100))
		h.write("<use xlink:href=\"#Ecoli\" x=\"0\" y=\"-%f\" />\n" % (2*currHeight*100))
		
	else:
		#figure out params of and write two identically sized cells that will overlap and pull apart to simulate pinching (two cells)
		svgParams = svgEcoliParams(finalCellsHeight)
		writeEcoliSVG(h,svgParams)
		h.write("<use xlink:href=\"#Ecoli\" x=\"0\" y=\"0\" />\n")
		h.write("<use xlink:href=\"#Ecoli\" x=\"0\" y=\"-%f\" />\n" % ((currHeight-finalCellsHeight)*100))
		#figure out params of and write the previous daughter cell that will also do this pinching (two cells)
		svgParams = svgEcoliParams(currHeight)
		writeEcoliSVG(h,svgParams)
		h.write("<use xlink:href=\"#Ecoli\" x=\"0\" y=\"-%f\" />\n" % (currHeight*100))
		h.write("<use xlink:href=\"#Ecoli\" x=\"0\" y=\"-%f\" />\n" % ((currHeight+(currHeight-finalCellsHeight))*100))

	#end boilerplate
	boilerplate_end(h)
	h.close()


def main():
	#go through a loop of the data, and call writeSVG on each set of relevant params. 
	directorySomehow = "/scratch/PI/mcovert/mhinks/MoCellRevamp"
	if not os.path.exists(directorySomehow):
		os.makedirs(directorySomehow)

	#get data from tsv file
	folder = "/scratch/PI/mcovert/thorst/allenTalk"
	reader = csv.reader(open(os.path.join(folder, "growthAnimation.tsv"), "rU"), delimiter="\t")
	reader.next()
	data = np.array(list(reader), dtype = np.float)

	#create a list of numpy arrays. Each numpy array is the time series of one cell, and contains
	#tuples of the time and mass of the cell at each time point. 
	allData = []
	currCell = []	
	for line in data:
		if line[0] <= line[2]:
			currCell.append((line[0], line[1]))
			if line[0] == line[2]:
				allData.append(currCell)
				currCell = []	

	#for each cell
	for cell in xrange(len(allData)):
		cellDensity = 1.0999999999999999e+18

		currCellData = allData[cell]
		finalTime = currCellData[len(currCellData)-1][0]
		finalMass = currCellData[len(currCellData)-1][1]
		finalVol = finalMass/cellDensity*10**15 #micrometers cubed 
		#at the last time step, we want to have two cells whose combined mass is the mass
		#specified in the last time step
		eachCellVol = finalVol/2.0 
		r = 0.5 #micrometers
		a = eachCellVol/(np.pi*r**2)-4./3*r
		eachCellHeight = a + 2.0*r
		finalHeight = eachCellHeight*2.
		

		for k in xrange(0,len(currCellData)):
			currCellMass = currCellData[k][1]			
			currCellVolume = currCellMass/cellDensity*10**15 #in micrometers cubed

			#removedVol is the volume of the intersecting spherical caps of the one cell pinching apart
			#into two cells. This volume is necessary to determine the offset of the daughter cell from
			#the parent as pinching occurs
			removedVol = finalVol - currCellVolume
			sphereVol = (4./3)*np.pi*r**3 #max overlap possible happens when it is a sphere

			#if the cell is not pinching, compute the height of the cell with a regular rod shape
			if removedVol > sphereVol:
				a = currCellVolume/(np.pi*r**2)-4./3*r
				currHeight = a + 2.0*r
				writeSVG(directorySomehow, currHeight, finalHeight, False, currCellData[k][0])
				
			#if the cell is pinching, compute height of cell accounting for the wonkiness in volume that
			#occurs during pinching 
			else:
				
				capVol = removedVol/2.0
				h = np.linspace(0.,.5,10000)
				#after a page of algebra, we know that we need a root of this cubic equation, 
				#which is easiest to find numerically like this
				cubic = np.absolute(np.pi/3.*h**3 - r*np.pi*h**2 + capVol)
				hVal = h[np.argmin(cubic)]
				currHeight = finalHeight-2.*hVal
				writeSVG(directorySomehow, currHeight, eachCellHeight, True, currCellData[k][0])
				

if __name__ == "__main__":
	main()









