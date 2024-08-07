"""
Plot dynamic traces of genes with low expression (on average less than 1 mRNA)

G7355_RNA	0.23	ypjD	Predicted inner membrane protein
EG11783_RNA	0.23	intA	CP4-57 prophage; integrase
G7742_RNA	0.23	yrfG	Purine nucleotidase
G6253_RNA	0.23	ylaC	Predicted inner membrane protein
EG10632_RNA	0.23	nagA	N-acetylglucosamine-6-phosphate deacetylase
EG11484_RNA	0.23	yigZ	Predicted elongation factor
G7889_RNA	0.23	lptG	LptG (part of LPS transport system)
EG10997_RNA	0.24	mnmE	GTPase, involved in modification of U34 in tRNA
EG10780_RNA	0.24	pspE	Thiosulfate sulfurtransferase
EG11060_RNA	0.24	ushA	UDP-sugar hydrolase / 5'-ribonucleotidase / 5'-deoxyribonucleotidase

(sorted from sim_data.transcription.process.rnaData, mrna only, elements in range -1505:-1495)
"""

import os

import numpy as np
from matplotlib import pyplot as plt

from wholecell.io.tablereader import TableReader
from wholecell.analysis.analysis_tools import exportFigure
from models.ecoli.analysis import singleAnalysisPlot


class Plot(singleAnalysisPlot.SingleAnalysisPlot):
	def do_plot(self, simOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		RNA_counts_reader = TableReader(os.path.join(simOutDir, 'RNACounts'))
		mRNA_cistron_counts = RNA_counts_reader.readColumn('mRNA_cistron_counts')
		all_mRNA_cistron_idx = {cistron: i for i, cistron in enumerate(RNA_counts_reader.readAttribute('mRNA_cistron_ids'))}

		cistron_ids = [
			"G7355_RNA", "EG11783_RNA", "G7742_RNA", "G6253_RNA", "EG10632_RNA",
			"EG11484_RNA", "G7889_RNA", "EG10997_RNA", "EG10780_RNA", "EG11060_RNA",
			]
		names = [
			"ypjD - Predicted inner membrane protein",
			"intA - CP4-57 prophage; integrase",
			"yrfG - Purine nucleotidase",
			"ylaC - Predicted inner membrane protein",
			"nagA - N-acetylglucosamine-6-phosphate deacetylase",
			"yigZ - Predicted elongation factor",
			"lptG - LptG (part of LPS transport system)",
			"mnmE - GTPase, involved in modification of U34 in tRNA",
			"pspE - Thiosulfate sulfurtransferase",
			"ushA - UDP-sugar hydrolase / 5'-ribonucleotidase / 5'-deoxyribonucleotidase",
		]

		cistron_indexes = np.array([all_mRNA_cistron_idx[x] for x in cistron_ids], int)
		cistron_counts = mRNA_cistron_counts[:, cistron_indexes]

		main_reader = TableReader(os.path.join(simOutDir, "Main"))
		initialTime = main_reader.readAttribute("initialTime")
		time = main_reader.readColumn("time") - initialTime

		plt.figure(figsize = (8.5, 11))

		for subplotIdx in range(1, 10):

			plt.subplot(3, 3, subplotIdx)

			plt.plot(time / 60., cistron_counts[:, subplotIdx])
			plt.xlabel("Time (min)")
			plt.ylabel("mRNA counts")
			plt.title(names[subplotIdx].split(" - ")[0])

		plt.subplots_adjust(hspace = 0.5, top = 0.95, bottom = 0.05)
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close("all")


if __name__ == "__main__":
	Plot().cli()
