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
from wholecell.io.tablereader import TableReader
import wholecell.utils.constants
from wholecell.utils import units

def main(seedOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata = None):
	if not os.path.isdir(seedOutDir):
		raise Exception, "seedOutDir does not currently exist as a directory"

	if not os.path.exists(plotOutDir):
		os.mkdir(plotOutDir)

	tsvOutDir = os.path.join(plotOutDir, "tsv_data")
	if not os.path.exists(tsvOutDir):
		os.mkdir(tsvOutDir)

	# Get all cells
	ap = AnalysisPaths(seedOutDir, multi_gen_plot = True)
	allDir = ap.get_cells()

	sim_data = cPickle.load(open(simDataFile, "rb"))
	genomeLength = sim_data.process.replication.genome_length

	# get names of regulated genes
	tfs = sim_data.conditionActiveTfs['with_aa']
	targets = []
	targetIndexes = []

	for target in sim_data.process.transcription_regulation.targetTf:
		for tf in sim_data.process.transcription_regulation.targetTf[target]:
			if tf in tfs:
				targets += [target]
				targetIndexes += np.where(sim_data.process.transcription.rnaData["id"] == target + '[c]')[0].tolist()
				break

	# get data from each generation
	time = []
	transcriptionEvents = []
	cellMass = []
	nOric = []
	dnaIdx = []
	dnaLength = []
	chromosomes = []
	ribosomeConc = []
	aaConc = []
	for gen, simDir in enumerate(allDir):
		simOutDir = os.path.join(simDir, "simOut")
		bulkReader = TableReader(os.path.join(simOutDir, 'BulkMolecules'))
		bulkMoleculeNames = bulkReader.readAttribute('objectNames')

		# sim time with first gen handled as cumsum because not saved as a float
		mainReader = TableReader(os.path.join(simOutDir, "Main"))
		if gen == 0:
			timeStep = mainReader.readColumn("timeStepSec")
			time += [0] + timeStep[:-1].cumsum().tolist()
		else:
			time += mainReader.readColumn("time").tolist()
		mainReader.close()

		# transcription events
		rnapDataReader = TableReader(os.path.join(simOutDir, "RnapData"))
		rnaInitEvent = rnapDataReader.readColumn("rnaInitEvent")[:, targetIndexes]
		rnapDataReader.close()

		# cell mass
		massReader = TableReader(os.path.join(simOutDir, "Mass"))
		mass = massReader.readColumn("cellMass")
		cellMass += mass.tolist()

		# DNA replication
		chromIdx = bulkMoleculeNames.index('CHROM_FULL[c]')
		chromosomes += bulkReader.readColumn('counts')[:, chromIdx].tolist()
		replicationReader = TableReader(os.path.join(simOutDir, 'ReplicationData'))
		nOric += replicationReader.readColumn('numberOfOric').tolist()

		# active ribosomes
		uniqueMoleculeCounts = TableReader(os.path.join(simOutDir, "UniqueMoleculeCounts"))
		ribosomeIdx = uniqueMoleculeCounts.readAttribute("uniqueMoleculeIds").index("activeRibosome")
		activeRibosomes = uniqueMoleculeCounts.readColumn("uniqueMoleculeCounts")[:, ribosomeIdx]
		molRibosomes = activeRibosomes / sim_data.constants.nAvogadro.asNumber(1 / units.mol)
		cellVolume = mass / sim_data.constants.cellDensity.asNumber(units.fg / units.L)
		ribosomeConc += (molRibosomes / cellVolume).tolist()

		# amino acid concentrations
		aaIdx = [bulkMoleculeNames.index(aa) for aa in sim_data.moleculeGroups.aaIDs]
		aaCounts = bulkReader.readColumn('counts')[:, aaIdx]
		aaMol = aaCounts / sim_data.constants.nAvogadro.asNumber(1 / units.mol)

		if gen == 0:
			transcriptionEvents = (rnaInitEvent != 0)
			dnaLength = replicationReader.readColumn('sequenceLength')
			dnaIdx = replicationReader.readColumn('sequenceIdx')
			aaConc = (aaMol.T / cellVolume).T
		else:
			transcriptionEvents = np.vstack((transcriptionEvents, (rnaInitEvent != 0)))
			dnaLength = np.vstack((dnaLength, replicationReader.readColumn('sequenceLength')))
			dnaIdx = np.vstack((dnaIdx, replicationReader.readColumn('sequenceIdx')))
			aaConc = np.vstack((aaConc, (aaMol.T / cellVolume).T))


	# write data for transcription events
	transcriptionWriter = csv.writer(open(os.path.join(tsvOutDir, 'transcriptionEvents.tsv'), 'w'), delimiter='\t')
	transcriptionWriter.writerow(['Time'] + targets)
	for t, row in zip(time, transcriptionEvents):
		transcriptionWriter.writerow([t] + row.tolist())

	# write data for cell mass
	cellMassWriter = csv.writer(open(os.path.join(tsvOutDir, 'cellMass.tsv'), 'w'), delimiter='\t')
	cellMassWriter.writerow(['Time', 'Cell Mass (fg)'])
	for t, m in zip(time, cellMass):
		cellMassWriter.writerow([t, m])

	# write all data for replication
	replicationWriter = csv.writer(open(os.path.join(tsvOutDir, 'replication.tsv'), 'w'), delimiter='\t')
	replicationWriter.writerow(['Time', 'Full Chromosomes', 'OriC'])
	for t, chrom, oric, lengths in zip(time, chromosomes, nOric, dnaLength):
		replicationWriter.writerow([t, chrom, oric] + lengths.tolist())

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
