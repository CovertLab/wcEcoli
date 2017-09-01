#!/usr/bin/env python
"""
Exports data for Allen Talk presentation animations to tsvs
@author: Travis Horst
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 8/7/2017
"""

import argparse
import os
import csv
import cPickle

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches

from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from wholecell.containers.bulk_objects_container import BulkObjectsContainer
from wholecell.io.tablereader import TableReader
import wholecell.utils.constants
from wholecell.utils import units

def getProteinMonomersDissociated(simOutDir, sim_data):

	tcsComplexToMonomers = sim_data.process.two_component_system.complexToMonomer
	ids_complexation = sim_data.process.complexation.moleculeNames
	ids_complexation_complexes = [ids_complexation[i] for i in np.where((sim_data.process.complexation.stoichMatrix() == 1).sum(axis = 1))[0]]
	ids_equilibrium = sim_data.process.equilibrium.moleculeNames
	ids_equilibrium_complexes = [ids_equilibrium[i] for i in np.where((sim_data.process.equilibrium.stoichMatrix() == 1).sum(axis = 1))[0]]
	ids_twoComponent = sim_data.process.two_component_system.moleculeNames.tolist()
	ids_twoComponent_complexes = sim_data.process.two_component_system.complexToMonomer.keys()
	ids_translation = sim_data.process.translation.monomerData["id"].tolist()
	ids_protein = sorted(set(ids_complexation + ids_equilibrium + ids_twoComponent + ids_translation))

	bulkContainer = BulkObjectsContainer(ids_protein, dtype = np.float64)
	view_complexation = bulkContainer.countsView(ids_complexation)
	view_complexation_complexes = bulkContainer.countsView(ids_complexation_complexes)
	view_equilibrium = bulkContainer.countsView(ids_equilibrium)
	view_equilibrium_complexes = bulkContainer.countsView(ids_equilibrium_complexes)
	view_twoComponent = bulkContainer.countsView(ids_twoComponent)
	view_twoComponent_complexes = bulkContainer.countsView(ids_twoComponent_complexes)
	view_translation = bulkContainer.countsView(ids_translation)

	bulkMolecules = TableReader(os.path.join(simOutDir, "BulkMolecules"))
	moleculeIds = bulkMolecules.readAttribute("objectNames")
	proteinIndexes = np.array([moleculeIds.index(moleculeId) for moleculeId in ids_protein], np.int)
	proteinCountsBulk = bulkMolecules.readColumn("counts")[:, proteinIndexes]
	bulkMolecules.close()

	# Account for unique molecules
	uniqueMoleculeCounts = TableReader(os.path.join(simOutDir, "UniqueMoleculeCounts"))
	ribosomeIndex = uniqueMoleculeCounts.readAttribute("uniqueMoleculeIds").index("activeRibosome")
	rnaPolyIndex = uniqueMoleculeCounts.readAttribute("uniqueMoleculeIds").index("activeRnaPoly")
	nActiveRibosome = uniqueMoleculeCounts.readColumn("uniqueMoleculeCounts")[:, ribosomeIndex]
	nActiveRnaPoly = uniqueMoleculeCounts.readColumn("uniqueMoleculeCounts")[:, rnaPolyIndex]
	uniqueMoleculeCounts.close()

	tfs = sim_data.conditionActiveTfs['with_aa']
	targets = []
	targetsPos = []
	targetsNeg = []
	targetIndexes = []

	for target in sim_data.process.transcription_regulation.targetTf:
		for tf in sim_data.process.transcription_regulation.targetTf[target]:
			if tf in tfs:
				if sim_data.tfToFC[tf][target] < 0.3:
					targetsNeg += [target]
				elif sim_data.tfToFC[tf][target] > 2:
					targetsPos += [target]
				targetIndexes += np.where(sim_data.process.transcription.rnaData["id"] == target + '[c]')[0].tolist()
				break

	targets = targetsPos + targetsNeg

	rnaIdToMonomerId = dict(sim_data.process.translation.monomerData.struct_array[["rnaId", "id"]].tolist())
	targetMonomers = [rnaIdToMonomerId[x + "[c]"] for x in targets]

	P = proteinCountsBulk.T.copy()

	PRowNames = bulkContainer._objectNames
	PRowNameToIdx = bulkContainer._objectIndex
	namesToIdxs = bulkContainer._namesToIndexes

	rowIdxs = namesToIdxs(sim_data.moleculeGroups.s30_fullComplex + sim_data.moleculeGroups.s50_fullComplex)
	P[rowIdxs, :] += nActiveRibosome

	rowIdxs = namesToIdxs(sim_data.moleculeGroups.rnapFull)
	P[rowIdxs, :] += nActiveRnaPoly

	rowIdxsInput = view_twoComponent_complexes._indexes
	rowIdxsOutput = view_twoComponent._indexes
	S = sim_data.process.two_component_system.stoichMatrixMonomers()
	P[rowIdxsOutput, :] += np.dot(S, P[rowIdxsInput, :] * -1).astype(np.int64)

	rowIdxsInput = view_equilibrium_complexes._indexes
	rowIdxsOutput = view_equilibrium._indexes
	S = sim_data.process.equilibrium.stoichMatrixMonomers()
	P[rowIdxsOutput, :] += np.dot(S, P[rowIdxsInput, :] * -1).astype(np.int64)

	rowIdxsInput = view_complexation_complexes._indexes
	rowIdxsOutput = view_complexation._indexes
	S = sim_data.process.complexation.stoichMatrixMonomers()
	P[rowIdxsOutput, :] += np.dot(S, P[rowIdxsInput, :] * -1).astype(np.int64)

	rowIdxs = namesToIdxs(targetMonomers)
	return P[rowIdxs, :].T, targetMonomers


def main(seedOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata = None):
	if not os.path.isdir(seedOutDir):
		raise Exception, "seedOutDir does not currently exist as a directory"

	if not os.path.exists(plotOutDir):
		os.mkdir(plotOutDir)

	tsvOutDir = os.path.join(plotOutDir, "tsv_data")
	if not os.path.exists(tsvOutDir):
		os.mkdir(tsvOutDir)

	reactants = [
		'GLC[p]',
		'FRUCTOSE-6P[c]',
		'2-PG[c]',
		'ACETYL-COA[c]',
		'SUC[c]',
		'D-6-P-GLUCONO-DELTA-LACTONE[c]',
		'RIBULOSE-5P[c]',
		'WATER[c]',
		'ACETYLSERINE[c]',
		'GLT[c]',
		'HISTIDINAL[c]',
		'GLT[c]',
		'L-ASPARTATE[c]',
		'O-PHOSPHO-L-HOMOSERINE[c]',
		'DADP[c]',
		'TDP[c]',
		'GDP[c]',
		'ADP[c]',
		'UDP[c]',
		'ADP[c]',
		]
	products = [
		'GLC-6-P[c]',
		'FRUCTOSE-16-DIPHOSPHATE[c]',
		'PHOSPHO-ENOL-PYRUVATE[c]',
		'CIT[c]',
		'FUM[c]',
		'CPD-2961[c]',
		'RIBOSE-5P[c]',
		'SER[c]',
		'CYS[c]',
		'PHE[c]',
		'HIS[c]',
		'LEU[c]',
		'ASN[c]',
		'THR[c]',
		'DATP[c]',
		'TTP[c]',
		'GTP[c]',
		'ATP[c]',
		'UTP[c]',
		'DADP[c]',
		]

	# Get all cells
	ap = AnalysisPaths(seedOutDir, multi_gen_plot = True)
	allDir = ap.get_cells()

	sim_data = cPickle.load(open(simDataFile, "rb"))
	genomeLength = sim_data.process.replication.genome_length
	rxnStoich = sim_data.process.metabolism.reactionStoich

	# get data from each generation
	time = []
	transcriptionEvents = []
	proteinConc = []
	cellMass = []
	nOric = []
	dnaIdx = []
	dnaLength = []
	chromosomes = []
	ribosomeConc = []
	aaConc = []
	fluxes = []
	for gen, simDir in enumerate(allDir):
		print gen
		simOutDir = os.path.join(simDir, "simOut")
		bulkReader = TableReader(os.path.join(simOutDir, 'BulkMolecules'))
		bulkMoleculeNames = bulkReader.readAttribute('objectNames')
		bulkCounts = bulkReader.readColumn("counts")
		bulkReader.close()

		# sim time with first gen handled as cumsum because not saved as a float
		mainReader = TableReader(os.path.join(simOutDir, "Main"))
		if gen == 0:
			timeStep = mainReader.readColumn("timeStepSec")
			time += [0] + timeStep[:-1].cumsum().tolist()
		else:
			time += mainReader.readColumn("time").tolist()
		mainReader.close()

		regulatedProteins, regulatedProteinNames = getProteinMonomersDissociated(simOutDir, sim_data)

		# cell mass
		massReader = TableReader(os.path.join(simOutDir, "Mass"))
		mass = massReader.readColumn("cellMass")
		cellMass += mass.tolist()
		countsToMolar = sim_data.constants.cellDensity.asNumber(units.fg / units.L) / sim_data.constants.nAvogadro.asNumber(1 / units.mol) / mass
		massReader.close()

		# DNA replication
		chromIdx = bulkMoleculeNames.index('CHROM_FULL[c]')
		chromosomes += bulkCounts[:, chromIdx].tolist()
		replicationReader = TableReader(os.path.join(simOutDir, 'ReplicationData'))
		nOric += replicationReader.readColumn('numberOfOric').tolist()

		# active ribosomes
		uniqueMoleculeCounts = TableReader(os.path.join(simOutDir, "UniqueMoleculeCounts"))
		ribosomeIdx = uniqueMoleculeCounts.readAttribute("uniqueMoleculeIds").index("activeRibosome")
		activeRibosomes = uniqueMoleculeCounts.readColumn("uniqueMoleculeCounts")[:, ribosomeIdx]
		ribosomeConc += (activeRibosomes * countsToMolar).tolist()

		# amino acid concentrations
		aaIdx = [bulkMoleculeNames.index(aa) for aa in sim_data.moleculeGroups.aaIDs]
		aaCounts = bulkCounts[:, aaIdx]

		# metabolism fluxes
		fbaReader = TableReader(os.path.join(simOutDir, "FBAResults"))
		reactionIds = fbaReader.readAttribute("reactionIDs")
		flux = fbaReader.readColumn("reactionFluxes")
		totalFlux = np.zeros_like(flux[:, :len(reactants)])
		for idx, (reactant, product) in enumerate(zip(reactants, products)):
			for rxn in rxnStoich:
				if reactant in rxnStoich[rxn] and product in rxnStoich[rxn]:
					if rxnStoich[rxn][reactant] < 0 and rxnStoich[rxn][product] > 0:
						direction = 1
					elif rxnStoich[rxn][reactant] > 0 and rxnStoich[rxn][product] < 0:
						direction = -1
					else:
						continue

					totalFlux[:, idx] += flux[:, reactionIds.index(rxn)] * direction
		fbaReader.close()

		if gen == 0:
			dnaLength = replicationReader.readColumn('sequenceLength')
			dnaIdx = replicationReader.readColumn('sequenceIdx')
			aaConc = (aaCounts.T * countsToMolar).T
			proteinConc = (regulatedProteins.T * countsToMolar).T
			fluxes = totalFlux
		else:
			dnaLength = np.vstack((dnaLength, replicationReader.readColumn('sequenceLength')))
			dnaIdx = np.vstack((dnaIdx, replicationReader.readColumn('sequenceIdx')))
			aaConc = np.vstack((aaConc, (aaCounts.T * countsToMolar).T))
			proteinConc = np.vstack((proteinConc, (regulatedProteins.T * countsToMolar).T))
			fluxes = np.vstack((fluxes, totalFlux))

	# write data for cell mass
	cellMassWriter = csv.writer(open(os.path.join(tsvOutDir, 'cellMass.tsv'), 'w'), delimiter='\t')
	cellMassWriter.writerow(['Time', 'Cell Mass (fg)'])
	for t, m in zip(time, cellMass):
		cellMassWriter.writerow([t, m])

	# write all data for replication
	replicationWriter = csv.writer(open(os.path.join(tsvOutDir, 'replication.tsv'), 'w'), delimiter='\t')
	replicationWriter.writerow(['Time', 'Full Chromosomes', 'OriC'])
	for t, chrom, oric, lengths in zip(time, chromosomes, nOric, dnaLength):
		replicationWriter.writerow([t, chrom, oric])

	# process replication data for chromosome animation
	fraction1 = []
	fraction2 = []
	for lengths in dnaLength:
		if lengths[0] == -1:
			fraction1 += [0]
			fraction2 += [0]
			continue
		temp1 = lengths[0]
		temp2 = 0
		for length in lengths[1:]:
			if temp1 / length < 0.9 or temp1 / length > 1.1:
				temp2 = length
				break

		if temp1 > temp2:
			if temp1 * 2 / genomeLength > 1:
				fraction1 += [1-1e-5]
			else:
				fraction1 += [temp1 * 2 / genomeLength]
			fraction2 += [temp2 * 2 / genomeLength]
		else:
			if temp2 * 2 / genomeLength > 1:
				fraction1 += [1-1e-5]
			else:
				fraction1 += [temp2 * 2 / genomeLength]
			fraction2 += [temp1 * 2 / genomeLength]

	# adjust time steps where chromosome replication has completed but chromosome has not formed
	for idx in np.where(np.array(chromosomes[1:]) - np.array(chromosomes[:-1]) == 1)[0]:
		chromosomes[idx] = 2

	repAnimationWriter = csv.writer(open(os.path.join(tsvOutDir, 'repAnimation.tsv'), 'w'), delimiter='\t')
	repAnimationWriter.writerow(['Time', 'Full Chromosome', 'Fraction 1', 'Fraction 2'])
	for t, c, f1, f2 in zip(time, chromosomes, fraction1, fraction2):
		repAnimationWriter.writerow([t, c, f1, f2])

	# write active ribosome data
	ribosomeWriter = csv.writer(open(os.path.join(tsvOutDir, 'ribosomes.tsv'), 'w'), delimiter='\t')
	ribosomeWriter.writerow(['Time', 'Active Ribosomes (M)'])
	for t, ribosomes in zip(time, ribosomeConc):
		ribosomeWriter.writerow([t, ribosomes])

	# write data for amino acid concentrations
	aaWriter = csv.writer(open(os.path.join(tsvOutDir, 'aas.tsv'), 'w'), delimiter='\t')
	aaWriter.writerow(['Time', 'Total AA Conc [M]'] + sim_data.moleculeGroups.aaIDs)
	for t, C in zip(time, aaConc):
		aaWriter.writerow([t, np.sum(C)] + C.tolist())

	# write data for cell growth animation
	nextDivision = np.array(time)
	prev = 0
	for idx in np.where(np.array(cellMass[1:]) < np.array(cellMass[:-1]) * 0.75)[0]:
		nextDivision[prev:idx] = time[idx]
		prev = idx + 1
	nextDivision[prev:] = time[-1]

	growthAnimationWriter = csv.writer(open(os.path.join(tsvOutDir, 'growthAnimation.tsv'), 'w'), delimiter='\t')
	growthAnimationWriter.writerow(['Time', 'Cell Mass (fg)', 'Time of Next Division'])
	for t, m, td in zip(time, cellMass, nextDivision):
		growthAnimationWriter.writerow([t, m, td])

	# write data for fluxes
	aveFlux = np.mean(fluxes[10000:11000, :], axis=0)
	adjustedFluxes = np.log2(fluxes / aveFlux)
	adjustedFluxes[np.isinf(adjustedFluxes)] = adjustedFluxes[~np.isinf(adjustedFluxes)].min() - 1

	fluxWriter = csv.writer(open(os.path.join(tsvOutDir, 'fluxes.tsv'), 'w'), delimiter='\t')
	fluxWriter.writerow(['Time'] + ['%s to %s' % (r, p) for (r, p) in zip(reactants, products)])
	for t, flux in zip(time, adjustedFluxes):
		fluxWriter.writerow([t] + flux.tolist())

	# write data for protein conc
	aveProteinConc = np.mean(proteinConc[10000:11000, :], axis=0)
	adjustedConc = np.log2(proteinConc / aveProteinConc)
	adjustedConc[np.isinf(adjustedConc)] = adjustedConc[~np.isinf(adjustedConc)].min() - 1

	proteinWriter = csv.writer(open(os.path.join(tsvOutDir, 'protein.tsv'), 'w'), delimiter='\t')
	proteinWriter.writerow(['Time'] + regulatedProteinNames)
	for t, conc in zip(time, adjustedConc):
		proteinWriter.writerow([t] + conc.tolist())

	import ipdb; ipdb.set_trace()

if __name__ == "__main__":
	defaultSimDataFile = os.path.join(
			wholecell.utils.constants.SERIALIZED_KB_DIR,
			wholecell.utils.constants.SERIALIZED_KB_MOST_FIT_FILENAME
			)

	parser = argparse.ArgumentParser()
	parser.add_argument("simOutDir", help = "Directory containing simulation output", type = str)
	parser.add_argument("plotOutDir", help = "Directory containing plot output (will get created if necessary)", type = str)
	parser.add_argument("plotOutFileName", help = "File name to produce", type = str)
	parser.add_argument("--simDataFile", help = "KB file name", type = str, default = defaultSimDataFile)
	parser.add_argument("--validationDataFile")

	args = parser.parse_args().__dict__

	main(args["simOutDir"], args["plotOutDir"], args["plotOutFileName"], args["simDataFile"], args["validationDataFile"])
