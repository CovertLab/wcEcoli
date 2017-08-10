import numpy as np
import os
import cPickle

def boilerplate_start(h):
	h.write("<html>\n")
	h.write("<body>\n")
	h.write("<svg width=\"120\" height=\"275\"  viewBox=\"0 0 120 275\">\n")

def boilerplate_end(h):
	h.write("</svg>\n")
	h.write("</body>\n")
	h.write("</html>")

def svgEcoliParams(height):
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
	h.write("   <g transform=\"translate(60,275) scale(1,-1)\">\n")
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
	#h.write("<use xlink:href=\"#Ecoli\" x=\"0\" y=\"-%f\" />\n" % yval)


def writeSVG(dirname, currHeight, finalHeight, threeCells, timestep):
	#put starting boilerplate (1 line)
	if not os.path.exists(dirname):
		os.makedirs(dirname)

	h = open(os.path.join(dirname, "moCell_%d.html" % timestep), "w")
	boilerplate_start(h)

	#figure out if you're going to make two or three cells, 
	#figure out positions and heights of cells according to that check
	if threeCells == False:
		svgParams = svgEcoliParams(currHeight)
		writeEcoliSVG(h,svgParams)
		h.write("<use xlink:href=\"#Ecoli\" x=\"0\" y=\"0\" />\n")
		h.write("<use xlink:href=\"#Ecoli\" x=\"0\" y=\"-%f\" />\n" % (currHeight*100))
		
	else:
		svgParams = svgEcoliParams(finalHeight/2.0)
		writeEcoliSVG(h,svgParams)
		h.write("<use xlink:href=\"#Ecoli\" x=\"0\" y=\"0\" />\n")
		h.write("<use xlink:href=\"#Ecoli\" x=\"0\" y=\"-%f\" />\n" % ((currHeight-finalHeight/2.0)*100))
		h.write("<use xlink:href=\"#Ecoli\" x=\"0\" y=\"-%f\" />\n" % (currHeight*100))
	
	#end boilerplate
	boilerplate_end(h)
	h.close()


def main():
	#go through a loop of the data, and call writeSVG on each set of relevant params. 
	directorySomehow = "svg_plots/MoCell2"
	if not os.path.exists(directorySomehow):
		os.makedirs(directorySomehow)
	#import data
	heightData = cPickle.load(open("heightData.cPickle","rb"))
	timeData = cPickle.load(open("time.cPickle","rb"))

	args = []

	#for each cell, store its initial size 
	for cell in range(0,2):
		currCellData = heightData[cell]
		finalHeight = currCellData[len(currCellData)-1]
		time = timeData[cell]
		finalTime = time[len(time)-1]

		for k in xrange(0,len(currCellData)):

			if time[k] < finalTime - 1800:
				writeSVG(directorySomehow, currCellData[k], finalHeight, False, time[k])
			else:
				writeSVG(directorySomehow, currCellData[k], finalHeight, True, time[k])
				print "I just switched to three cells at time %f" % time[k] 

	


def test():
	dirname = "svg_plots"
	if not os.path.exists(dirname):
		os.makedirs(dirname)
	writeSVG(dirname, 1.8, 2.5, False, 1)


main()

#







