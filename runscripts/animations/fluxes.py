'''
Creates series of svg files that can be stitched together using avconv to create an animation of fluxes
Input: fluxes.tsv file in same file location with time column and column for each flux
Output: time series of svg images to folder svg
'''

import numpy as np
import os
import csv
from matplotlib import pyplot as plt
from multiprocessing import Pool

def write_svg((dirname, idx, fluxes)):
	fluxes = [float(f) for f in fluxes]
	plt.bar(range(len(fluxes)), fluxes)
	plt.axis([0, 20, -5, 5])
	plt.savefig(os.path.join(dirname, "flux_%05d.svg" % idx))
	plt.close("all")

def main():
	folder = os.path.dirname(__file__)
	dirname = os.path.join(folder, "out")
	if not os.path.exists(dirname):
		os.makedirs(dirname)

	reader = csv.reader(open(os.path.join(folder, "fluxes.tsv"), "r"), delimiter="\t")
	reader.next()

	data = []
	for line in reader:
		data.append([dirname, float(line[0]), line[1:]])

	pool = Pool(processes = 16)
	pool.map(write_svg, data)
	pool.close()

if __name__ == "__main__":
	main()
