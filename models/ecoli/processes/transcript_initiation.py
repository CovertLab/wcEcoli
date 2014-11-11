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

		# Load parameters

		rnaSynthProbPopAvg = kb.rnaData['synthProb']
		self.rnaIds = kb.rnaData['geneId']
		self.bulkChromosome = sim.states['BulkChromosome']
		self.geneIds = kb.geneData['name']
		self.geneView = self.bulkChromosome.container.countsView(self.geneIds)
		self.mapGeneRna = [np.where(self.geneIds==x)[0][0] for x in self.rnaIds]

		# self.activationProb = kb.transcriptionActivationRate.asNumber(1/units.s) * self.timeStepSec # TODO: consider the validity of this math

		rnaLengths = kb.rnaData["length"]

		expectedTranscriptionTime = 1./kb.rnaPolymeraseElongationRate * rnaLengths

		expectedTranscriptionTimesteps = np.ceil(
			(1/(self.timeStepSec * units.s) * expectedTranscriptionTime).asNumber()
			)

		averageTranscriptionTimesteps = np.dot(kb.rnaData["synthProb"], expectedTranscriptionTimesteps)

		expectedTerminationRate = 1./averageTranscriptionTimesteps

		expectedFractionTimeInactive = np.dot(
			1 - (1/(self.timeStepSec * units.s) * expectedTranscriptionTime).asNumber() / expectedTranscriptionTimesteps,
			kb.rnaData["synthProb"]
			)

		effectiveFractionActive = kb.fracActiveRnap * 1 / (1 - expectedFractionTimeInactive)

		self.activationProb = effectiveFractionActive * expectedTerminationRate / (1 - effectiveFractionActive)

		# Views

		self.activeRnaPolys = self.uniqueMoleculesView('activeRnaPoly')

		self.inactiveRnaPolys = self.bulkMoleculeView("APORNAP-CPLX[c]")
		
		# Calculate rna synth probabilities
		geneEndCoordinates = kb.geneData['endCoordinate']

		minDistFromOriC = np.minimum(np.abs(kb.oriCCenter.asNumber()-geneEndCoordinates-kb.genomeLength),
						np.abs(geneEndCoordinates-kb.oriCCenter.asNumber()))

		ageReplicated = minDistFromOriC / kb.dnaPolymeraseElongationRate.asNumber()

		self.rnaSynthProb = rnaSynthProbPopAvg / (2 * np.exp(-np.log(2)*ageReplicated/kb.cellCycleLen.asNumber()))

		self.rnaSynthProb /= self.rnaSynthProb.sum()


	def calculateRequest(self):
		self.inactiveRnaPolys.requestAll()


	# Calculate temporal evolution
	def evolveState(self):
		#Recalculate the synthesis probabilities, taking gene copy number into account
		#p=copy number * synthprob
		#then renormalize
		copyNumber = self.geneView.counts()[self.mapGeneRna]
		#copyNumber = [self.geneView.counts()[i] for i in [np.where(self.geneIds==x)[0][0] for x in self.rnaIds]]

		synthProbWeightedByCopyNumber = self.rnaSynthProb * copyNumber

		synthProbRenormalized = synthProbWeightedByCopyNumber/np.sum(synthProbWeightedByCopyNumber)

		# Sample a multinomial distribution of synthesis probabilities to 
		# determine what molecules are initialized

		inactiveRnaPolyCount = self.inactiveRnaPolys.count()

		rnaPolyToActivate = np.int64(self.activationProb * inactiveRnaPolyCount)

		if rnaPolyToActivate == 0:
			return

		nNewRnas = self.randomState.multinomial(rnaPolyToActivate,
			synthProbRenormalized)

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

		self.writeToListener("rnaCounts", "rnaSynthProb", synthProbRenormalized)
