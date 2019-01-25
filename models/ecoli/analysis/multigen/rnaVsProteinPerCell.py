"""
Plots average RNA counts per cell vs average protein counts per cell.

@author: Heejo Choi
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 11/1/2017
"""

from __future__ import absolute_import

import os
import cPickle

import numpy as np
import matplotlib.pyplot as plt

from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from wholecell.io.tablereader import TableReader
from wholecell.containers.bulk_objects_container import BulkObjectsContainer
import matplotlib.lines as mlines
from wholecell.analysis.analysis_tools import exportFigure
from models.ecoli.analysis import multigenAnalysisPlot

complexToMonomer = {
	"CPLX0-7620[c]": "PD00260[c]", # CPLX0-7620's monomer is EG10359-MONOMER, which is ID'ed as PD00260 (proteins.tsv)
	"CPLX0-8801[c]": "G6420-MONOMER[c]",
	"CPLX0-7677[c]": "EG11969-MONOMER[c]",
	"CPLX0-7702[c]": "G6263-MONOMER[c]",
	"CPLX0-7701[c]": "G6420-MONOMER[c]",
	}

monomerToTranslationMonomer = {
	"MONOMER0-1781[c]": "EG11519-MONOMER[c]", # MONOMER0-1781 is a complex, EG11519 is its monomer
	"EG10359-MONOMER[c]": "PD00260[c]", # EG10359-MONOMER is not the ID of fur monomer, it's PD00260 (proteins.tsv)
	}

PLOT_ZEROS_ON_LINE = 2.5e-6


class Plot(multigenAnalysisPlot.MultigenAnalysisPlot):
	def do_plot(self, seedOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		if not os.path.isdir(seedOutDir):
			raise Exception, "seedOutDir does not currently exist as a directory"

		if not os.path.exists(plotOutDir):
			os.mkdir(plotOutDir)

		# Get all cells
		ap = AnalysisPaths(seedOutDir, multi_gen_plot = True)
		allDir = ap.get_cells()

		# Load sim data
		sim_data = cPickle.load(open(simDataFile, "rb"))

		# Load rna IDs, ordered by monomer IDs
		rnaIds = sim_data.process.transcription.rnaData["id"][sim_data.relation.rnaIndexToMonomerMapping]

		# Identify sub-generationally transcribed genes
		nonzeroSumRnaCounts_allGens = []
		for i, simDir in enumerate(allDir):
			simOutDir = os.path.join(simDir, "simOut")

			# Read counts of transcripts
			bulkMolecules = TableReader(os.path.join(simOutDir, "BulkMolecules"))
			if i == 0:
				moleculeIds = bulkMolecules.readAttribute("objectNames")
				rnaIndices = np.array([moleculeIds.index(x) for x in rnaIds])
			rnaCounts = bulkMolecules.readColumn("counts")[:, rnaIndices]
			bulkMolecules.close()

			# Sum counts over timesteps
			sumRnaCounts = rnaCounts.sum(axis=0)

			# Flag where the sum is nonzero (True if nonzero, False if zero)
			nonzeroSumRnaCounts = sumRnaCounts != 0
			nonzeroSumRnaCounts_allGens.append(nonzeroSumRnaCounts)

		# Average (mean) over generations
		nonzeroSumRnaCounts_allGens = np.array(nonzeroSumRnaCounts_allGens)
		avgRnaCounts = nonzeroSumRnaCounts_allGens.mean(axis=0)

		# Identify subgenerationally transcribed genes
		subgenRnaIndices = np.where(np.logical_and(avgRnaCounts != 0., avgRnaCounts != 1.))[0]
		subgenRnaIds = rnaIds[subgenRnaIndices]

		# Generate 'colors' for use in plotting (y = 0, b = subgen, r = 1)
		colors = ["y" if avgRnaCount == 0 else "r" if avgRnaCount == 1 else "b" for avgRnaCount in avgRnaCounts]
		mrnaIds = rnaIds.tolist()

		# Make views for monomers
		ids_complexation = sim_data.process.complexation.moleculeNames
		ids_complexation_complexes = sim_data.process.complexation.ids_complexes
		ids_equilibrium = sim_data.process.equilibrium.moleculeNames
		ids_translation = sim_data.process.translation.monomerData["id"].tolist()
		ids_protein = sorted(set(ids_complexation + ids_equilibrium + ids_translation))
		bulkContainer = BulkObjectsContainer(ids_protein, dtype = np.float64)
		view_complexation_complexes = bulkContainer.countsView(ids_complexation_complexes)

		# Identify monomers that are subunits for multiple complexes
		monomersInvolvedInManyComplexes = []
		monomersInvolvedInComplexes = []
		for complexId in ids_complexation_complexes:
			subunitIds = sim_data.process.complexation.getMonomers(complexId)["subunitIds"]
			for subunitId in subunitIds:
				if subunitId in monomersInvolvedInComplexes:
					monomersInvolvedInManyComplexes.append(subunitId)
				monomersInvolvedInComplexes.append(subunitId)
		monomersInvolvedInManyComplexes_id = list(set(monomersInvolvedInManyComplexes))
		monomersInvolvedInManyComplexes_dict = {}
		for x in monomersInvolvedInManyComplexes_id:
			monomersInvolvedInManyComplexes_dict[x] = {}

		# Get average (over timesteps) counts for all generations
		avgRnaCounts_forAllCells = np.zeros(rnaIds.shape[0], np.float64)
		avgProteinCounts_forAllCells = np.zeros(rnaIds.shape[0], np.float64)
		for i, simDir in enumerate(allDir):
			simOutDir = os.path.join(simDir, "simOut")

			# Account for bulk molecules
			bulkMolecules = TableReader(os.path.join(simOutDir, "BulkMolecules"))
			bulkMoleculeCounts = bulkMolecules.readColumn("counts")
			moleculeIds = bulkMolecules.readAttribute("objectNames")
			rnaIndexes = np.array([moleculeIds.index(moleculeId) for moleculeId in rnaIds], np.int)
			avgRnaCounts = bulkMoleculeCounts[:, rnaIndexes].mean(axis = 0)
			bulkMolecules.close()

			# Average counts of monomers
			monomerCountsReader = TableReader(os.path.join(simOutDir, "MonomerCounts"))
			monomerCounts = monomerCountsReader.readColumn("monomerCounts")

			if i == 0:
				# Skip first few time steps since complexes have not formed yet
				avgMonomerCounts = np.average(monomerCounts[5:, :], axis = 0)
			else:
				avgMonomerCounts = np.average(monomerCounts, axis = 0)

			# Get counts of "functional units" (ie. complexed forms)
			avgProteinCounts = avgMonomerCounts[:]
			avgComplexCounts = view_complexation_complexes.counts()

			for j, complexId in enumerate(ids_complexation_complexes):
				# Map all subunits to the average counts of the complex (ignores counts of monomers)
				# Some subunits are involved in multiple complexes - these cases are kept track
				subunitIds = sim_data.process.complexation.getMonomers(complexId)["subunitIds"]

				for subunitId in subunitIds:
					if subunitId not in ids_translation:
						if subunitId in monomerToTranslationMonomer:
							# couple monomers have different ID in ids_translation
							subunitId = monomerToTranslationMonomer[subunitId]
						elif "CPLX" in subunitId:
							# few transcription factors are complexed with ions
							subunitId = complexToMonomer[subunitId]
						elif "RNA" in subunitId:
							continue

					if subunitId not in monomersInvolvedInManyComplexes_id:
						avgProteinCounts[ids_translation.index(subunitId)] = avgComplexCounts[j]
					else:
						if complexId not in monomersInvolvedInManyComplexes_dict[subunitId]:
							monomersInvolvedInManyComplexes_dict[subunitId][complexId] = 0.
						monomersInvolvedInManyComplexes_dict[subunitId][complexId] += avgComplexCounts[j]

			# Store
			avgRnaCounts_forAllCells += avgRnaCounts
			avgProteinCounts_forAllCells += avgProteinCounts


		# Per cell
		avgRnaCounts_perCell = avgRnaCounts_forAllCells / float(len(allDir))
		avgProteinCounts_perCell = avgProteinCounts_forAllCells / float(len(allDir))

		# Plot
		fig, ax = plt.subplots(1, 1, figsize = (10, 10))

		for monomer in monomersInvolvedInManyComplexes_id:
			index = ids_translation.index(monomer)
			color_index = mrnaIds.index(rnaIds[index])
			color = colors[color_index]

			for complexId in monomersInvolvedInManyComplexes_dict[monomer]:
				avgComplexCount = monomersInvolvedInManyComplexes_dict[monomer][complexId] / float(len(allDir))

				if avgComplexCount == 0:
					ax.loglog(avgRnaCounts_perCell[index], 2.5e-6, alpha = 0.5, marker = ".", lw = 0., color = color)

				else:
					if avgRnaCounts_perCell[index] == 0:
						ax.loglog(PLOT_ZEROS_ON_LINE, avgComplexCount, alpha = 0.5, marker = ".", lw = 0., color = color)
					else:
						ax.loglog(avgRnaCounts_perCell[index], avgComplexCount, alpha = 0.5, marker = ".", lw = 0., color = color)

		# plot monomers that are not involved in complexes or involved in only 1 complex
		monomersInvolvedInManyComplexes_index = [ids_translation.index(x) for x in monomersInvolvedInManyComplexes_id]
		A = [x for x in xrange(len(ids_translation)) if x not in monomersInvolvedInManyComplexes_index]
		for i in A:
			color = colors[mrnaIds.index(rnaIds[i])]
			ax.loglog(avgRnaCounts_perCell[i], avgProteinCounts_perCell[i], alpha = 0.5, marker = ".", lw = 0., color = color)

		# Plot genes with zero transcripts on an arbitrary line
		noTranscripts_indices = [x for x in np.where(avgRnaCounts_perCell == 0)[0] if x not in monomersInvolvedInManyComplexes_index]
		for i in noTranscripts_indices:
			color = colors[mrnaIds.index(rnaIds[i])]
			ax.loglog(PLOT_ZEROS_ON_LINE, avgProteinCounts_perCell[i], alpha = 0.5, marker = ".", lw = 0., color = color)

		ax.set_title("Each (translatable) gene's functional unit is represented as a point\n(ie. x points per gene "
					 "where x is number of complexes the monomer is involved in)\n(avg across {0} "
					 "generations)".format(len(allDir)))
		ax.set_xlabel("<RNA> per cell")
		ax.set_ylabel("<Functional units (protein)> per cell")
		ax.tick_params(which = "both", direction = "out")

		plt.subplots_adjust(hspace = 0.5, wspace = 0.5, left = 0.1, bottom = 0.1, top = 0.9, right = 0.95)

		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close("all")


if __name__ == "__main__":
	Plot().cli()
