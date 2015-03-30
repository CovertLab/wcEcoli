#!/usr/bin/env python

"""
TranscriptInitiation

Transcription initiation sub-model.

TODO:
- use transcription units instead of single genes
- match sigma factors to promoters
- implement transcriptional regulation
- modulate initiation probabilities as a function of gene copy number
- match measured levels of active RNA polymerase instead of initiating to completion

@author: John Mason
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 4/26/14
"""

from __future__ import division

import numpy as np

import wholecell.processes.process
from wholecell.utils import units

import itertools

class TranscriptInitiation(wholecell.processes.process.Process):
	""" TranscriptInitiation """

	_name = "TranscriptInitiation"

	# Constructor
	def __init__(self):
		# Parameters
		self.rnaSynthProb = None

		# Views
		self.activeRnaPolys = None
		self.inactiveRnaPolys = None

		super(TranscriptInitiation, self).__init__()


	# Construct object graph
	def initialize(self, sim, kb):
		super(TranscriptInitiation, self).initialize(sim, kb)

		# Quick hack to check expression before and after fitting
		#f=open('data_initEqualSynthProbs', 'w')
		#for i in range(0,len(kb.unfitSynthProb)): 
		#	f.write(str(kb.unfitSynthProb[i])+' ')
		#	f.write(str(kb.process.transcription.rnaData['synthProb'][i])+ '\n')
		#f.close()

		# Load parameters
		#rnaSynthProbPopAvg = kb.process.transcription.rnaData['synthProb']
		self.rnaSynthProb = kb.process.transcription.rnaData['synthProb']
		self.rnaIds = kb.process.transcription.rnaData['geneId']
		self.bulkChromosome = sim.states['BulkChromosome']
		self.geneIds = kb.state.bulkChromosome.bulkData['id']
		self.geneView = self.bulkChromosome.container.countsView(self.geneIds)
		self.mapGeneRna = [np.where(self.geneIds==x)[0][0] for x in self.rnaIds]

		# self.activationProb = kb.transcriptionActivationRate.asNumber(1/units.s) * self.timeStepSec # TODO: consider the validity of this math

		rnaLengths = kb.process.transcription.rnaData["length"]

		expectedTranscriptionTime = 1./kb.constants.rnaPolymeraseElongationRate * rnaLengths

		expectedTranscriptionTimesteps = np.ceil(
			(1/(self.timeStepSec * units.s) * expectedTranscriptionTime).asNumber()
			)

		averageTranscriptionTimesteps = np.dot(kb.process.transcription.rnaData["synthProbTimeAvg"], expectedTranscriptionTimesteps)
		expectedTerminationRate = 1./averageTranscriptionTimesteps

		expectedFractionTimeInactive = np.dot(
			1 - (1/(self.timeStepSec * units.s) * expectedTranscriptionTime).asNumber() / expectedTranscriptionTimesteps,
			kb.process.transcription.rnaData["synthProbTimeAvg"]
			)

		effectiveFractionActive = kb.fracActiveRnap * 1 / (1 - expectedFractionTimeInactive)

		self.activationProb = effectiveFractionActive * expectedTerminationRate / (1 - effectiveFractionActive)

		# Views

		self.activeRnaPolys = self.uniqueMoleculesView('activeRnaPoly')

		self.inactiveRnaPolys = self.bulkMoleculeView("APORNAP-CPLX[c]")
		self.highlyRegulated = kb.process.transcription.rnaData['isHighlyRegulated']
		
		# Do this in the knowledge base prior to fitting
		#### Calculate rna synth probabilities
		##geneEndCoordinates = kb.geneData['endCoordinate']

		##minDistFromOriC = np.minimum(np.abs(kb.oriCCenter.asNumber()-geneEndCoordinates-kb.genomeLength),
		##				np.abs(geneEndCoordinates-kb.oriCCenter.asNumber()))

		##ageReplicated = minDistFromOriC / kb.dnaPolymeraseElongationRate.asNumber()

		##self.rnaSynthProb = rnaSynthProbPopAvg / (2 * np.exp(-np.log(2)*ageReplicated/kb.cellCycleLen.asNumber()))

		##self.rnaSynthProb /= self.rnaSynthProb.sum()

		### Instead try account for exponential increase in number of RNAP molecules (assuming concentration increases exponentially too)
		##geneEndCoordinates = kb.geneData['endCoordinate']

		##minDistFromOriC = np.minimum(np.abs(kb.oriCCenter.asNumber()-geneEndCoordinates-kb.genomeLength),
		##				np.abs(geneEndCoordinates-kb.oriCCenter.asNumber()))

		##ageReplicated = minDistFromOriC / kb.dnaPolymeraseElongationRate.asNumber()
		##RNAP0=1800.

		##self.rnaSynthProb = rnaSynthProbPopAvg / (2*RNAP0*2*np.log(2)-ageReplicated*RNAP0*2*np.log(2)/kb.cellCycleLen.asNumber())


	def calculateRequest(self):
		self.inactiveRnaPolys.requestAll()


	# Calculate temporal evolution
	def evolveState(self):
		#Recalculate the synthesis probabilities, taking gene copy number into account
		#But only for the genes that aren't highly regulated: do this by always setting copy number for highly regulated genes to 1
		#p=copy number * synthprob
		#then renormalize (only for not highly regulated)
		copyNumber = self.geneView.counts()[self.mapGeneRna]
		copyNumber[self.highlyRegulated] = np.ones(len(copyNumber[self.highlyRegulated]))
		#copyNumber = [self.geneView.counts()[i] for i in [np.where(self.geneIds==x)[0][0] for x in self.rnaIds]]

		synthProbWeightedByCopyNumber = self.rnaSynthProb * copyNumber
		synthProbWeightedByCopyNumber[~self.highlyRegulated] = (1-np.sum(synthProbWeightedByCopyNumber[self.highlyRegulated])) * synthProbWeightedByCopyNumber[~self.highlyRegulated]/np.sum(synthProbWeightedByCopyNumber[~self.highlyRegulated])

		#synthProbRenormalized = synthProbWeightedByCopyNumber/np.sum(synthProbWeightedByCopyNumber)

		# Sample a multinomial distribution of synthesis probabilities to 
		# determine what molecules are initialized

		inactiveRnaPolyCount = self.inactiveRnaPolys.count()

		rnaPolyToActivate = np.int64(self.activationProb * inactiveRnaPolyCount)

		if rnaPolyToActivate == 0:
			return

		nNewRnas = self.randomState.multinomial(rnaPolyToActivate,
			synthProbWeightedByCopyNumber)

		nonzeroCount = (nNewRnas > 0)

		assert nNewRnas.sum() == rnaPolyToActivate

		# Build list of RNA indexes

		rnaIndexes = np.empty(rnaPolyToActivate, np.int64)

		startIndex = 0
		for rnaIndex, counts in itertools.izip(
				np.arange(nNewRnas.size)[nonzeroCount],
				nNewRnas[nonzeroCount]
				):

			rnaIndexes[startIndex:startIndex+counts] = rnaIndex

			startIndex += counts

		# Create the active RNA polymerases

		activeRnaPolys = self.activeRnaPolys.moleculesNew(
			"activeRnaPoly",
			rnaPolyToActivate
			)

		activeRnaPolys.attrIs(
			rnaIndex = rnaIndexes
			)

		self.inactiveRnaPolys.countDec(nNewRnas.sum())

		self.writeToListener("rnaCounts", "rnaSynthProb", synthProbWeightedByCopyNumber)
