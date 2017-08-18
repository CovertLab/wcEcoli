'''
Creates series of svg files that can be stitched together using avconv to create an animation of protein concentrations
Input: protein.tsv file in same file location with time column and column for each protein conc
Output: time series of svg images to folder svg
'''

import numpy as np
import os
import csv
from matplotlib import pyplot as plt
from multiprocessing import Pool

def write_svg((dirname, idx, proteins)):
	proteins = [float(p) for p in proteins]
	plt.bar(range(len(proteins)), proteins)
	plt.axis([0, 20, -5, 5])
	plt.savefig(os.path.join(dirname, "protein_%05d.svg" % idx))
	plt.close("all")

def main():
	folder = os.path.dirname(__file__)
	dirname = os.path.join(folder, "svg")
	if not os.path.exists(dirname):
		os.makedirs(dirname)

	reader = csv.reader(open(os.path.join(folder, "protein.tsv"), "r"), delimiter="\t")
	reader.next()

	data = []
	for line in reader:
		data.append([dirname, float(line[0]), line[1:]])

	pool = Pool(processes = 16)
	pool.map(write_svg, data)
	pool.close()

if __name__ == "__main__":
	main()
