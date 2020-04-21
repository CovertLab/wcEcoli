from __future__ import absolute_import


import os
import csv
import numpy as np
from matplotlib import pyplot as plt
import itertools

from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from wholecell.io.tablereader import TableReader

from wholecell.analysis.plotting_tools import COLORS_LARGE
from wholecell.analysis.analysis_tools import exportFigure
from models.ecoli.analysis import variantAnalysisPlot


class Plot(variantAnalysisPlot.VariantAnalysisPlot):
	def do_plot(self, inputDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		# Look through functional_genes.tsv to match each KO sim to corresponding type.
		# First, pull list of geneIDs from functional_genes.tsv,
		# then at each file, read metadata and check with vectors of geneID (match indices)

		functional_genes_path = '/Users/petersuzuki/repos/wcEcoli/runscripts/reflect/functional_genes.tsv' #ToDo: make relative path
		#functional_genes_path = '/home/users/psuzuki/functional_genes.tsv'

		geneIDs = []
		process_type = []
		# Read functional_genes.tsv and get geneIDs as list of string
		with open(functional_genes_path) as tsvfile:
			tsvreader = csv.reader(tsvfile, delimiter="\t")
			tsvreader.next()
			tsvreader.next() # skip header
			for line in tsvreader:
				geneID = line[1].strip('_RNA')
				geneIDs.append(geneID)
				process_type.append(line[-1])


		massNames = [
					"dryMass",
					"proteinMass",
					#"tRnaMass",
					"rRnaMass",
					'mRnaMass',
					"dnaMass"
					]

		cleanNames = [
					"Dry\nmass",
					"Protein\nmass",
					#"tRNA\nmass",
					"rRNA\nmass",
					"mRNA\nmass",
					"DNA\nmass"
					]

		if not os.path.isdir(inputDir):
			raise Exception, "inputDir does not currently exist as a directory"

		ap = AnalysisPaths(inputDir, variant_plot=True)
		all_cells = ap.get_cells()

		# Build a mapping from variant id to color
		idToColor = {}
		for idx, (cell_id, color) in enumerate(itertools.izip(all_cells, itertools.cycle(COLORS_LARGE))):
			idToColor[idx] = color

		if not os.path.exists(plotOutDir):
			os.mkdir(plotOutDir)

		mass_fold_change = np.zeros(len(massNames))
		variant_info = []

		#fig, axesList = plt.subplots(len(massNames), sharex=True,figsize=[10,10])

		currentMaxTime = 0
		for cellIdx, simDir in enumerate(all_cells):
			with open(os.path.join(simDir[:-32],'metadata','short_name')) as f:
				variant_name = [line for line in f][0]

			# Check variant_name against functional_genes.tsv listing to get process
			variant_type = []
			for index in range(len(process_type)):
				if geneIDs[index] in variant_name:
					variant_type = process_type[index]
			if 'wildtype' in variant_name:
				variant_type = "Wildtype"
			if variant_type == []:
				variant_type = "Unknown"


			simOutDir = os.path.join(simDir, "simOut")
			if len(os.listdir(simOutDir)) == 0:
				value = [variant_name, variant_type]
				variant_info.append(value)
				continue

			time = TableReader(os.path.join(simOutDir, "Main")).readColumn("time")
			if time.size > 1:
				cellCycleTime = ((time[-1] - time[0]) / 60. / 60.)

				mass = TableReader(os.path.join(simOutDir, "Mass"))

				for massIdx, massType in enumerate(massNames):
					massToPlot = mass.readColumn(massType)
					mass_fold_change[massIdx] = massToPlot[-1] / massToPlot[0]

				value = [variant_name, variant_type] + mass_fold_change.tolist()
				value.append(cellCycleTime)
				variant_info.append(value)
			else:
				value = [variant_name, variant_type]
				variant_info.append(value)


		titles = ['Variant Name','Variant Type'] + massNames + ['Cell Cycle Time']

		filename = plotOutDir + plotOutFileName + '.csv'

		# writing to csv file
		with open(filename, 'w') as csvfile:
			# creating a csv writer object
			csvwriter = csv.writer(csvfile)

			csvwriter.writerow(titles)
			csvwriter.writerows(variant_info)

		#exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		#plt.close("all")


if __name__ == "__main__":
	Plot().cli()