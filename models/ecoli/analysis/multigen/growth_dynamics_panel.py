from __future__ import absolute_import
from __future__ import division

import os

import numpy as np
from matplotlib import pyplot as plt
import matplotlib.patches as patches
import cPickle

from wholecell.io.tablereader import TableReader
from wholecell.utils import units
from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from wholecell.analysis.analysis_tools import exportFigure
from models.ecoli.analysis import multigenAnalysisPlot

class Plot(multigenAnalysisPlot.MultigenAnalysisPlot):
	def do_plot(self, seedOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		if not os.path.isdir(seedOutDir):
			raise Exception, "seedOutDir does not currently exist as a directory"

		if not os.path.exists(plotOutDir):
			os.mkdir(plotOutDir)

		ap = AnalysisPaths(seedOutDir, multi_gen_plot = True)

		allDirs = ap.get_cells()

		# Load data from KB
		sim_data = cPickle.load(open(simDataFile, "rb"))
		nAvogadro = sim_data.constants.nAvogadro
		cellDensity = sim_data.constants.cellDensity
		oriC = sim_data.constants.oriCCenter.asNumber()
		terC = sim_data.constants.terCCenter.asNumber()
		genomeLength = len(sim_data.process.replication.genome_sequence)

		recruitmentColNames = sim_data.process.transcription_regulation.recruitmentColNames
		tfs = sorted(set([x.split("__")[-1] for x in recruitmentColNames if x.split("__")[-1] != "alpha"]))
		trpRIndex = [i for i, tf in enumerate(tfs) if tf == "CPLX-125"][0]

		tfBoundIds = [target + "__CPLX-125" for target in sim_data.tfToFC["CPLX-125"].keys()]
		synthProbIds = [target + "[c]" for target in sim_data.tfToFC["CPLX-125"].keys()]

		plt.figure(figsize = (10, 10))
		nRows = 3
		nCols = 3

		for simDir in allDirs:
			simOutDir = os.path.join(simDir, "simOut")
			# Load time
			initialTime = 0#TableReader(os.path.join(simOutDir, "Main")).readAttribute("initialTime")
			time = TableReader(os.path.join(simOutDir, "Main")).readColumn("time") - initialTime

			# Load mass data
			# Total cell mass is needed to compute concentrations (since we have cell density)
			# Protein mass is needed to compute the mass fraction of the proteome that is trpA
			massReader = TableReader(os.path.join(simOutDir, "Mass"))
			cellMass = units.fg * massReader.readColumn("cellMass")
			cellMass_no_conv = massReader.readColumn("cellMass")
			proteinMass = units.fg * massReader.readColumn("proteinMass")

			# Get instantanous growth rate
			growth_rate = massReader.readColumn("instantaniousGrowthRate")

			massReader.close()

			# Load data from ribosome data listener
			#ribosomeDataReader = TableReader(os.path.join(simOutDir, "RibosomeData"))
			#nTrpATranslated = ribosomeDataReader.readColumn("numTrpATerminated")
			#ribosomeDataReader.close()

			# Load data from bulk molecules
			bulkMoleculesReader = TableReader(os.path.join(simOutDir, "BulkMolecules"))
			bulkMoleculeIds = bulkMoleculesReader.readAttribute("objectNames")

			# Get the concentration of intracellular trp
			trpId = ["TRP[c]"]
			trpIndex = np.array([bulkMoleculeIds.index(x) for x in trpId])
			trpCounts = bulkMoleculesReader.readColumn("counts")[:, trpIndex].reshape(-1)
			trpMols = 1. / nAvogadro * trpCounts
			volume = cellMass / cellDensity
			trpConcentration = trpMols * 1. / volume

			# Get the amount of active trpR (that isn't promoter bound)
			#trpRActiveId = ["CPLX-125[c]"]
			#trpRActiveIndex = np.array([bulkMoleculeIds.index(x) for x in trpRActiveId])
			#trpRActiveCounts = bulkMoleculesReader.readColumn("counts")[:, trpRActiveIndex].reshape(-1)

			# Get the amount of inactive trpR
			#trpRInactiveId = ["PC00007[c]"]
			#trpRInactiveIndex = np.array([bulkMoleculeIds.index(x) for x in trpRInactiveId])
			#trpRInactiveCounts = bulkMoleculesReader.readColumn("counts")[:, trpRInactiveIndex].reshape(-1)

			# Get the amount of monomeric trpR
			#trpRMonomerId = ["PD00423[c]"]
			#trpRMonomerIndex = np.array([bulkMoleculeIds.index(x) for x in trpRMonomerId])
			#trpRMonomerCounts = bulkMoleculesReader.readColumn("counts")[:, trpRMonomerIndex].reshape(-1)

			# Get the promoter-bound status for all regulated genes
			tfBoundIndex = np.array([bulkMoleculeIds.index(x) for x in tfBoundIds])
			tfBoundCounts = bulkMoleculesReader.readColumn("counts")[:, tfBoundIndex]

			# Get the amount of monomeric trpA
			trpAProteinId = ["TRYPSYN-APROTEIN[c]"]
			trpAProteinIndex = np.array([bulkMoleculeIds.index(x) for x in trpAProteinId])
			trpAProteinCounts = bulkMoleculesReader.readColumn("counts")[:, trpAProteinIndex].reshape(-1)

			# Get the amount of complexed trpA
			trpABComplexId = ["TRYPSYN[c]"]
			trpABComplexIndex = np.array([bulkMoleculeIds.index(x) for x in trpABComplexId])
			trpABComplexCounts = bulkMoleculesReader.readColumn("counts")[:, trpABComplexIndex].reshape(-1)

			# Get the amount of trpA mRNA
			trpARnaId = ["EG11024_RNA[c]"]
			trpARnaIndex = np.array([bulkMoleculeIds.index(x) for x in trpARnaId])
			trpARnaCounts = bulkMoleculesReader.readColumn("counts")[:, trpARnaIndex].reshape(-1)

			#Get mass per Origin

			mass_per_oriC = TableReader(os.path.join(simOutDir, "ReplicationData")).readColumn("criticalMassPerOriC")

			#RNA over protein:
			#Get active ribosome counts
			uniqueMoleculeCounts = TableReader(os.path.join(simOutDir, "UniqueMoleculeCounts"))
			ribosomeIndex = uniqueMoleculeCounts.readAttribute("uniqueMoleculeIds").index("activeRibosome")
			ribosomeCounts = uniqueMoleculeCounts.readColumn("uniqueMoleculeCounts")[:, ribosomeIndex]
			uniqueMoleculeCounts.close()

			#Find the sequence index and length (to find the fork position later):
			sequenceIdx = TableReader(os.path.join(simOutDir, "ReplicationData")).readColumn("sequenceIdx")
			sequenceLength = TableReader(os.path.join(simOutDir, "ReplicationData")).readColumn("sequenceLength")
			sequenceLength[sequenceLength == -1] = np.nan

			bulkMoleculesReader.close()

			# Computations:

			# Compute total counts and concentration of trpA in monomeric and complexed form
			# (we know the stoichiometry)
			trpAProteinTotalCounts = trpAProteinCounts + 2 * trpABComplexCounts
			#trpAProteinTotalMols = 1. / nAvogadro * trpAProteinTotalCounts
			#trpAProteinTotalConcentration = trpAProteinTotalMols * 1. / volume

			# Compute concentration of trpA mRNA
			#trpARnaMols = 1. / nAvogadro * trpARnaCounts
			#trpARnaConcentration = trpARnaMols * 1. / volume

			# Compute the trpA mass in the cell
			#trpAMw = sim_data.getter.getMass(trpAProteinId)
			#trpAMass = 1. / nAvogadro * trpAProteinTotalCounts * trpAMw

			# Compute the proteome mass fraction
			#proteomeMassFraction = trpAMass.asNumber(units.fg) / proteinMass.asNumber(units.fg)

			# Get the synthesis probability for all regulated genes
			#rnaSynthProbReader = TableReader(os.path.join(simOutDir, "RnaSynthProb"))

			#rnaIds = rnaSynthProbReader.readAttribute("rnaIds")
			#synthProbIndex = np.array([rnaIds.index(x) for x in synthProbIds])
			#synthProbs = rnaSynthProbReader.readColumn("rnaSynthProb")[:, synthProbIndex]

			#trpRBound = rnaSynthProbReader.readColumn("nActualBound")[:,trpRIndex]

			#rnaSynthProbReader.close()

			# Calculate total trpR - active, inactive, bound and monomeric
			#trpRTotalCounts = 2 * (trpRActiveCounts + trpRInactiveCounts + trpRBound) + trpRMonomerCounts

			# Compute moving averages
			width = 100

			tfBoundCountsMA = np.array([np.convolve(tfBoundCounts[:,i], np.ones(width) / width, mode = "same")
					for i in range(tfBoundCounts.shape[1])]).T
			#synthProbsMA = np.array([np.convolve(synthProbs[:,i], np.ones(width) / width, mode = "same")
			#		for i in range(synthProbs.shape[1])]).T

			# Find the index of initialization:
			idxInit = np.where(mass_per_oriC >= 1)[0]

			# Calculate the growth rate:
			growth_rate = (1 / units.s) * growth_rate
			growth_rate = growth_rate.asNumber(1 / units.min)
			growth_rate[abs(growth_rate - np.median(growth_rate)) > 1.25 * np.nanstd(growth_rate)] = np.nan

			# Calculate Ribosome Concentration:
			ribosomeConcentration = ((1 / sim_data.constants.nAvogadro) * ribosomeCounts) / ((1.0 / sim_data.constants.cellDensity) * (cellMass))
			ribosomeConcentration = ribosomeConcentration.asNumber(units.umol / units.L)

			# Fork Position:
			reverseIdx = 1
			reverseCompIdx = 3
			reverseSequences = np.logical_or(sequenceIdx == reverseIdx, sequenceIdx == reverseCompIdx)
			sequenceLength[reverseSequences] = -1 * sequenceLength[reverseSequences]

			# Down sample dna polymerase position, every position is only plotted once here
			unique, index, value = np.unique(sequenceLength, return_index=True, return_inverse=True)
			m = np.zeros_like(value, dtype=bool)
			m[index] = True
			m = m.reshape(sequenceLength.shape)
			sequenceLength[~m] = np.nan

			# Pairs of forks
			pairsOfForks = (sequenceIdx != -1).sum(axis = 1) / 4

			##############################################################
			ax9 = plt.subplot(nRows, nCols, 9)
			ax9.plot(time, trpConcentration.asNumber(units.umol / units.L), color = '#0d71b9')
			plt.ylabel("Internal TRP Conc. (uM)", fontsize = 6)

			#ymin, ymax = ax.get_ylim()
			ax9.set_ylim([0, 400])
			ax9.set_yticks([0, 400])
			ax9.set_yticklabels([0, 400])
			ax9.spines['top'].set_visible(False)
			#ax.spines['bottom'].set_visible(False)
			ax9.spines['right'].set_visible(False)
			ax9.xaxis.set_ticks_position('bottom')
			ax9.tick_params(which = 'both', direction = 'out', labelsize = 6)

			#ax.set_xticklabels([0., 360.])
			ax9.set_xticks([0, time.max()])
			ax9.set_xticklabels([0., np.round(time.max() / 60., decimals = 0)])
		
			##############################################################
			'''
			##############################################################
			ax = plt.subplot(nRows, 1, 2)
			ax.plot(time, trpRActiveCounts)
			ax.plot(time, trpRInactiveCounts)
			ax.plot(time, trpRTotalCounts)
			plt.ylabel("TrpR Counts", fontsize = 6)
			plt.legend(["Active (dimer)", "Inactive (dimer)", "Total (monomeric)"], fontsize = 6)

			ymin, ymax = ax.get_ylim()
			ax.set_yticks([ymin, ymax])
			ax.set_yticklabels(["%0.0f" % ymin, "%0.0f" % ymax])
			ax.spines['top'].set_visible(False)
			ax.spines['bottom'].set_visible(False)
			ax.xaxis.set_ticks_position('none')
			ax.tick_params(which = 'both', direction = 'out', labelsize=6)
			ax.set_xticks([])
			##############################################################
			'''
			##############################################################
			ax = plt.subplot(nRows, nCols, 8)
			ax.plot(time, tfBoundCountsMA, color = '#0d71b9')
			#comment out color in order to see colors per generation
			plt.ylabel("TrpR Bound To Promoters\n(Moving Average)", fontsize = 6)

			#ymin, ymax = ax.get_ylim()
			ax.set_ylim([0, 1])
			ax.set_yticks([0, 1])
			ax.set_yticklabels([0, 1])
			ax.spines['top'].set_visible(False)
			#ax.spines['bottom'].set_visible(False)
			ax.spines['right'].set_visible(False)
			ax.xaxis.set_ticks_position('bottom')
			ax.tick_params(which = 'both', direction = 'out', labelsize = 6)
			
			#xmin, xmax = ax.get_xlim()
			ax.set_xticks([0, time.max()])
			ax.set_xticklabels([0., np.round(time.max() / 60., decimals = 0)])	
			##############################################################
			'''
			##############################################################
			ax = plt.subplot(nRows, 1, 4)
			ax.plot(time, synthProbsMA)
			plt.ylabel("Regulated Gene Synthesis Prob.\n(Moving Average)", fontsize = 6)

			ymin, ymax = ax.get_ylim()
			ax.set_yticks([ymin, ymax])
			ax.set_yticklabels(["%0.2e" % ymin, "%0.2e" % ymax])
			ax.spines['top'].set_visible(False)
			ax.spines['bottom'].set_visible(False)
			ax.xaxis.set_ticks_position('none')
			ax.tick_params(which = 'both', direction = 'out', labelsize = 6)
			ax.set_xticks([])
			##############################################################
			'''
			##############################################################
			ax = plt.subplot(nRows, nCols, 6)
			ax.plot(time, trpAProteinTotalCounts, color = '#0d71b9')
			plt.ylabel("TrpA Counts", fontsize = 6)

			#ymin, ymax = ax.get_ylim()
			#ax.set_yticks([ymin, ymax])
			#ax.set_yticklabels(["%0.0f" % ymin, "%0.0f" % ymax])
			ax.set_ylim([800, 6000])
			ax.set_yticks([800,6000])
			ax.set_yticklabels([800, 6000])
			ax.spines['top'].set_visible(False)
			ax.spines['bottom'].set_visible(False)
			ax.spines['right'].set_visible(False)
			ax.xaxis.set_ticks_position('none')
			ax.tick_params(which = 'both', direction = 'out', labelsize = 6)
			ax.set_xticks([])
			##############################################################
			'''
			##############################################################
			ax = plt.subplot(nRows, 1, 6)
			ax.plot(time, trpAProteinTotalConcentration.asNumber(units.umol / units.L), color = '#0d71b9')
			plt.ylabel("TrpA Concentration", fontsize = 6)

			ymin, ymax = ax.get_ylim()
			ax.set_yticks([ymin, ymax])
			ax.set_yticklabels(["%0.2f" % ymin, "%0.2f" % ymax])
			ax.spines['top'].set_visible(False)
			ax.spines['bottom'].set_visible(False)
			ax.xaxis.set_ticks_position('none')
			ax.tick_params(which = 'both', direction = 'out', labelsize = 6)
			ax.set_xticks([])
			##############################################################
			'''
			##############################################################
			ax = plt.subplot(nRows, nCols, 3)
			ax.plot(time, trpARnaCounts, color = '#0d71b9')
			plt.ylabel("TrpA mRNA Counts", fontsize = 6)

			ymin, ymax = ax.get_ylim()
			ax.set_ylim([0, ymax])
			ax.set_yticks([0, ymax])
			ax.set_yticklabels([0, "%0.0f" % ymax])
			ax.spines['top'].set_visible(False)
			ax.spines['bottom'].set_visible(False)
			ax.spines['right'].set_visible(False)
			ax.xaxis.set_ticks_position('none')
			ax.tick_params(which = 'both', direction = 'out', labelsize = 6)
			ax.set_xticks([])
			##############################################################
			'''
			##############################################################
			ax = plt.subplot(nRows, 1, 8)
			ax.plot(time, trpARnaConcentration.asNumber(units.umol / units.L), color = '#0d71b9')
			plt.ylabel("TrpA mRNA Concentration", fontsize = 6)

			ymin, ymax = ax.get_ylim()
			ax.set_yticks([ymin, ymax])
			ax.set_yticklabels(["%0.2e" % ymin, "%0.2e" % ymax])
			ax.spines['top'].set_visible(False)
			ax.spines['bottom'].set_visible(False)
			ax.xaxis.set_ticks_position('none')
			ax.tick_params(which = 'both', direction = 'out', labelsize = 6)
			ax.set_xticks([])
			##############################################################

			##############################################################
			ax = plt.subplot(nRows, 1, 9)
			ax.plot(time / 3600., proteomeMassFraction, color = '#0d71b9')
			plt.ylabel("TrpA MAss FRaction of Proteome", fontsize = 6)

			ymin, ymax = ax.get_ylim()
			ax.set_yticks([ymin, ymax])
			ax.set_yticklabels(["%0.2e" % ymin, "%0.2e" % ymax])
			ax.spines['top'].set_visible(False)
			ax.spines['bottom'].set_visible(False)
			# ax.xaxis.set_ticks_position('none')
			ax.tick_params(which = 'both', direction = 'out', labelsize = 6)
			ax.set_xticks(ax.get_xlim())
			##############################################################

			##############################################################
			ax = plt.subplot(nRows, 1, 10)
			ax.plot(time, nTrpATranslated, color = '#0d71b9')
			plt.ylabel("Number of TrpA translation events", fontsize = 6)

			ymin, ymax = ax.get_ylim()
			ax.set_yticks([ymin, ymax])
			ax.set_yticklabels(["%0.0f" % ymin, "%0.0f" % ymax])
			ax.spines['top'].set_visible(False)
			ax.spines['bottom'].set_visible(False)
			ax.xaxis.set_ticks_position('none')
			ax.tick_params(which = 'both', direction = 'out', labelsize = 6)
			ax.set_xticks([])
			##############################################################
			'''
			##############################################################
			ax = plt.subplot(nRows, nCols, 1)
			ax.plot(time, cellMass_no_conv, color = '#0d71b9')
			ax.plot(time[idxInit], cellMass_no_conv[idxInit],  markersize=6, linewidth=0, marker="o", color = "#ed2224", markeredgewidth=0)
			plt.ylabel("Cell Mass (fg)", fontsize = 6)

			#ymin, ymax = ax.get_ylim()
			ax.set_ylim([1000, 4600])
			ax.set_yticks([1000, 4600])
			#ax.set_yticklabels(["%0.0f" % ymin, "%0.0f" % ymax])
			ax.spines['top'].set_visible(False)
			ax.spines['bottom'].set_visible(False)
			ax.spines['right'].set_visible(False)
			ax.xaxis.set_ticks_position('none')
			ax.tick_params(which = 'both', direction = 'out', labelsize = 6)
			ax.set_xticks([])
			##############################################################

			##############################################################
			ax = plt.subplot(nRows, nCols, 4)
			ax.plot(time, growth_rate, color = '#0d71b9')
			plt.ylabel("Instantaneouse growth Rate", fontsize = 6)

			ax.set_ylabel(r"$\mu$ $(\frac{gDCW}{gDCW \cdot \, min})$")
			
			ymin, ymax = ax.get_ylim()
			ax.set_ylim([0, 0.032])
			ax.set_yticks([0, 0.032])
			#ax.set_yticklabels(["%0.0f" % ymin, "%0.0f" % ymax])
			ax.spines['top'].set_visible(False)
			ax.spines['bottom'].set_visible(False)
			ax.spines['right'].set_visible(False)
			ax.xaxis.set_ticks_position('none')
			ax.tick_params(which = 'both', direction = 'out', labelsize = 6)
			ax.set_xticks([])
			##############################################################

			##############################################################
			ax = plt.subplot(nRows, nCols, 7)
			ax.plot(time, ribosomeConcentration, color = '#0d71b9')
			ax.set_ylim([15., 25.])
			plt.ylabel("Active Ribosome (umol/L)", fontsize = 6)

			ymin, ymax = ax.get_ylim()
			ax.set_yticks([ymin, ymax])
			ax.set_yticklabels(["%0.0f" % ymin, "%0.0f" % ymax])
			ax.spines['top'].set_visible(False)
			#ax.spines['bottom'].set_visible(False)
			ax.spines['right'].set_visible(False)
			ax.xaxis.set_ticks_position('bottom')
			ax.tick_params(which = 'both', direction = 'out', labelsize = 6)
			ax.set_xticks([0, time.max()])			
			ax.set_xticklabels([0., np.round(time.max() / 60., decimals = 0)])

			##############################################################

			##############################################################
			ax = plt.subplot(nRows, nCols, 2)
			ax.plot(time, sequenceLength, marker='o', markersize=2, linewidth=0, color = '#0d71b9')
			plt.ylabel("DNA polymerase position", fontsize = 6)
			ymin, ymax = ax.get_ylim()
			ax.set_yticks([-1 * genomeLength / 2, 0, genomeLength / 2])
			ax.set_yticklabels(['-terC', 'oriC', '+terC'])
			ax.spines['top'].set_visible(False)
			ax.spines['bottom'].set_visible(False)
			ax.spines['right'].set_visible(False)
			ax.xaxis.set_ticks_position('none')
			ax.tick_params(which = 'both', direction = 'out', labelsize = 6)
			ax.set_xticks([])
			##############################################################

			##############################################################
			ax = plt.subplot(nRows, nCols, 5)
			ax.plot(time, pairsOfForks, linewidth=2, color = '#0d71b9')
			ax.set_yticks(np.arange(0,7))
			ax.set_ylim([0, 6.])
			ax.set_yticks([0, 6.])
			ax.set_yticklabels(["0","6"])
			plt.ylabel("Relative rate of dNTP polymerization", fontsize = 6)

			ax.spines['top'].set_visible(False)
			ax.spines['bottom'].set_visible(False)
			ax.spines['right'].set_visible(False)
			ax.xaxis.set_ticks_position('none')
			ax.tick_params(which = 'both', direction = 'out', labelsize = 6)
			ax.set_xticks([])
			##############################################################


		plt.subplots_adjust(hspace = 1)

		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close("all")


if __name__ == "__main__":
	Plot().cli()