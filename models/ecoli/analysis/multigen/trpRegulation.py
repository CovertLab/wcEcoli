"""
Plot trp regulation
"""

import os
import pickle

import numpy as np
from matplotlib import pyplot as plt

from wholecell.io.tablereader import TableReader
from wholecell.utils import units
from wholecell.analysis.analysis_tools import exportFigure
from models.ecoli.analysis import multigenAnalysisPlot


class Plot(multigenAnalysisPlot.MultigenAnalysisPlot):
	def do_plot(self, seedOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):

		allDirs = self.ap.get_cells()

		# Load data from KB
		sim_data = self.read_pickle_file(simDataFile)
		nAvogadro = sim_data.constants.n_avogadro
		cellDensity = sim_data.constants.cell_density

		# Get list of TF and transcription unit IDs from first simOut directory
		simOutDir = os.path.join(allDirs[0], "simOut")
		rna_synth_prob_reader = TableReader(os.path.join(simOutDir, "RnaSynthProb"))
		tf_ids = rna_synth_prob_reader.readAttribute("tf_ids")
		rna_ids = rna_synth_prob_reader.readAttribute("rnaIds")

		trpRIndex = tf_ids.index("CPLX-125")
		target_ids = sim_data.relation.tf_id_to_target_RNAs["CPLX-125"]
		target_idx = np.array([rna_ids.index(target_id) for target_id in target_ids])

		plt.figure(figsize = (10, 15))
		nRows = 10

		for simDir in allDirs:
			simOutDir = os.path.join(simDir, "simOut")
			# Load time
			initialTime = 0  # TableReader(os.path.join(simOutDir, "Main")).readAttribute("initialTime")
			time = TableReader(os.path.join(simOutDir, "Main")).readColumn("time") - initialTime

			# Load mass data
			# Total cell mass is needed to compute concentrations (since we have cell density)
			# Protein mass is needed to compute the mass fraction of the proteome that is trpA
			massReader = TableReader(os.path.join(simOutDir, "Mass"))
			cellMass = units.fg * massReader.readColumn("cellMass")
			proteinMass = units.fg * massReader.readColumn("proteinMass")

			# Load data from ribosome data listener
			ribosomeDataReader = TableReader(os.path.join(simOutDir, "RibosomeData"))
			nTrpATranslated = ribosomeDataReader.readColumn("numTrpATerminated")

			# Load data from bulk molecules
			bulkMoleculesReader = TableReader(os.path.join(simOutDir, "BulkMolecules"))
			bulkMoleculeIds = bulkMoleculesReader.readAttribute("objectNames")
			bulkMoleculeCounts = bulkMoleculesReader.readColumn("counts")

			# Load data from mRNA counts listener
			RNA_counts_reader = TableReader(
				os.path.join(simOutDir, 'RNACounts'))
			all_mRNA_cistron_ids = RNA_counts_reader.readAttribute('mRNA_cistron_ids')
			mRNA_cistron_counts = RNA_counts_reader.readColumn('mRNA_cistron_counts')

			# Load data from RnaSynthProb listener
			rna_synth_prob_reader = TableReader(
				os.path.join(simOutDir, "RnaSynthProb"))
			n_bound_TF_per_TU = rna_synth_prob_reader.readColumn(
				"n_bound_TF_per_TU").reshape((-1, len(rna_ids), len(tf_ids)))

			# Get the concentration of intracellular trp
			trpId = ["TRP[c]"]
			trpIndex = np.array([bulkMoleculeIds.index(x) for x in trpId])
			trpCounts = bulkMoleculeCounts[:, trpIndex].reshape(-1)
			trpMols = 1. / nAvogadro * trpCounts
			volume = cellMass / cellDensity
			trpConcentration = trpMols * 1. / volume

			# Get the amount of active trpR (that isn't promoter bound)
			trpRActiveId = ["CPLX-125[c]"]
			trpRActiveIndex = np.array([bulkMoleculeIds.index(x) for x in trpRActiveId])
			trpRActiveCounts = bulkMoleculeCounts[:, trpRActiveIndex].reshape(-1)

			# Get the amount of inactive trpR
			trpRInactiveId = ["PC00007[c]"]
			trpRInactiveIndex = np.array([bulkMoleculeIds.index(x) for x in trpRInactiveId])
			trpRInactiveCounts = bulkMoleculeCounts[:, trpRInactiveIndex].reshape(-1)

			# Get the amount of monomeric trpR
			trpRMonomerId = ["PD00423[c]"]
			trpRMonomerIndex = np.array([bulkMoleculeIds.index(x) for x in trpRMonomerId])
			trpRMonomerCounts = bulkMoleculeCounts[:, trpRMonomerIndex].reshape(-1)

			# Get the promoter-bound status for all regulated genes
			tfBoundCounts = n_bound_TF_per_TU[:, target_idx, trpRIndex]

			# Get the amount of monomeric trpA
			trpAProteinId = ["TRYPSYN-APROTEIN[c]"]
			trpAProteinIndex = np.array([bulkMoleculeIds.index(x) for x in trpAProteinId])
			trpAProteinCounts = bulkMoleculeCounts[:, trpAProteinIndex].reshape(-1)

			# Get the amount of complexed trpA
			trpABComplexId = ["TRYPSYN[c]"]
			trpABComplexIndex = np.array([bulkMoleculeIds.index(x) for x in trpABComplexId])
			trpABComplexCounts = bulkMoleculeCounts[:, trpABComplexIndex].reshape(-1)

			# Get the amount of trpA mRNA cistrons
			trpA_cistron_id = ["EG11024_RNA"]
			trpA_cistron_index = np.array([all_mRNA_cistron_ids.index(x) for x in trpA_cistron_id])
			trpA_cistron_counts = mRNA_cistron_counts[:, trpA_cistron_index].reshape(-1)

			# Compute total counts and concentration of trpA in monomeric and complexed form
			# (we know the stoichiometry)
			trpAProteinTotalCounts = trpAProteinCounts + 2 * trpABComplexCounts
			trpAProteinTotalMols = 1. / nAvogadro * trpAProteinTotalCounts
			trpAProteinTotalConcentration = trpAProteinTotalMols * 1. / volume

			# Compute concentration of trpA mRNA
			trpA_cistron_mols = 1. / nAvogadro * trpA_cistron_counts
			trpA_cistron_concentration = trpA_cistron_mols * 1. / volume

			# Compute the trpA mass in the cell
			trpAMw = sim_data.getter.get_masses(trpAProteinId)
			trpAMass = 1. / nAvogadro * trpAProteinTotalCounts * trpAMw

			# Compute the proteome mass fraction
			proteomeMassFraction = trpAMass.asNumber(units.fg) / proteinMass.asNumber(units.fg)

			# Get the synthesis probability for all regulated genes
			synthProbs = rna_synth_prob_reader.readColumn("actual_rna_synth_prob")[:, target_idx]
			trpRBound = rna_synth_prob_reader.readColumn("nActualBound")[:, trpRIndex]

			rna_synth_prob_reader.close()

			# Calculate total trpR - active, inactive, bound and monomeric
			trpRTotalCounts = 2 * (trpRActiveCounts + trpRInactiveCounts + trpRBound) + trpRMonomerCounts

			# Compute moving averages
			width = 100

			tfBoundCountsMA = np.array([np.convolve(tfBoundCounts[:,i], np.ones(width) / width, mode = "same")
					for i in range(tfBoundCounts.shape[1])]).T
			synthProbsMA = np.array([np.convolve(synthProbs[:,i], np.ones(width) / width, mode = "same")
					for i in range(synthProbs.shape[1])]).T


			##############################################################
			ax = self.subplot(nRows, 1, 1)
			ax.plot(time, trpConcentration.asNumber(units.umol / units.L), color = "b")
			plt.ylabel("Internal TRP Conc. [uM]", fontsize = 6)

			ymin, ymax = ax.get_ylim()
			ax.set_yticks([ymin, ymax])
			ax.set_yticklabels(["%0.0f" % ymin, "%0.0f" % ymax])
			ax.spines['top'].set_visible(False)
			ax.spines['bottom'].set_visible(False)
			ax.xaxis.set_ticks_position('none')
			ax.tick_params(which = 'both', direction = 'out', labelsize = 6)
			ax.set_xticks([])
			##############################################################

			##############################################################
			ax = self.subplot(nRows, 1, 2)
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

			##############################################################
			ax = self.subplot(nRows, 1, 3)
			ax.plot(time, tfBoundCountsMA)
			plt.ylabel("TrpR Bound To Promoters\n(Moving Average)", fontsize = 6)

			ymin, ymax = ax.get_ylim()
			ax.set_yticks([ymin, ymax])
			ax.set_yticklabels(["%0.0f" % ymin, "%0.0f" % ymax])
			ax.spines['top'].set_visible(False)
			ax.spines['bottom'].set_visible(False)
			ax.xaxis.set_ticks_position('none')
			ax.tick_params(which = 'both', direction = 'out', labelsize = 6)
			ax.set_xticks([])
			##############################################################

			##############################################################
			ax = self.subplot(nRows, 1, 4)
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

			##############################################################
			ax = self.subplot(nRows, 1, 5)
			ax.plot(time, trpAProteinTotalCounts, color = "b")
			plt.ylabel("TrpA Counts", fontsize = 6)

			ymin, ymax = ax.get_ylim()
			ax.set_yticks([ymin, ymax])
			ax.set_yticklabels(["%0.0f" % ymin, "%0.0f" % ymax])
			ax.spines['top'].set_visible(False)
			ax.spines['bottom'].set_visible(False)
			ax.xaxis.set_ticks_position('none')
			ax.tick_params(which = 'both', direction = 'out', labelsize = 6)
			ax.set_xticks([])
			##############################################################

			##############################################################
			ax = self.subplot(nRows, 1, 6)
			ax.plot(time, trpAProteinTotalConcentration.asNumber(units.umol / units.L), color = "b")
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

			##############################################################
			ax = self.subplot(nRows, 1, 7)
			ax.plot(time, trpA_cistron_counts, color = "b")
			plt.ylabel("TrpA mRNA cistron counts", fontsize = 6)

			ymin, ymax = ax.get_ylim()
			ax.set_yticks([ymin, ymax])
			ax.set_yticklabels(["%0.0f" % ymin, "%0.0f" % ymax])
			ax.spines['top'].set_visible(False)
			ax.spines['bottom'].set_visible(False)
			ax.xaxis.set_ticks_position('none')
			ax.tick_params(which = 'both', direction = 'out', labelsize = 6)
			ax.set_xticks([])
			##############################################################

			##############################################################
			ax = self.subplot(nRows, 1, 8)
			ax.plot(time, trpA_cistron_concentration.asNumber(units.umol / units.L), color = "b")
			plt.ylabel("TrpA mRNA cistron\nconcentration", fontsize = 6)

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
			ax = self.subplot(nRows, 1, 9)
			ax.plot(time / 3600., proteomeMassFraction, color = "b")
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
			ax = self.subplot(nRows, 1, 10)
			ax.plot(time, nTrpATranslated, color = "b")
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

		plt.subplots_adjust(hspace = 1)

		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close("all")


if __name__ == "__main__":
	Plot().cli()
