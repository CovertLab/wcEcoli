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

#PPGPP_POWER = 1.5
PPGPP_POWER = 3.5

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
		self.rnaSynthProb = kb.process.transcription.rnaData["synthProb"]

		self.rRNAIdx = kb.process.transcription.rnaData['isRRna']
		self.tRNAIdx = kb.process.transcription.rnaData['isTRna']
		self.mRNAIdx = kb.process.transcription.rnaData['isMRna']
		self.rnaLengths = kb.process.transcription.rnaData['length'].asNumber(units.nt)
		self.nAvogadro = kb.constants.nAvogadro
		self.cellDensity = kb.constants.cellDensity

		# self.activationProb = kb.transcriptionActivationRate.asNumber(1/units.s) * self.timeStepSec # TODO: consider the validity of this math

		rnaLengths = kb.process.transcription.rnaData["length"]

		expectedTranscriptionTime = 1./kb.growthRateParameters.rnaPolymeraseElongationRate * rnaLengths

		expectedTranscriptionTimesteps = np.ceil(
			(1/(self.timeStepSec * units.s) * expectedTranscriptionTime).asNumber()
			)

		averageTranscriptionTimesteps = np.dot(kb.process.transcription.rnaData["synthProb"], expectedTranscriptionTimesteps)

		expectedTerminationRate = 1./averageTranscriptionTimesteps

		expectedFractionTimeInactive = np.dot(
			1 - (1/(self.timeStepSec * units.s) * expectedTranscriptionTime).asNumber() / expectedTranscriptionTimesteps,
			kb.process.transcription.rnaData["synthProb"]
			)

		effectiveFractionActive = kb.fracActiveRnap * 1 / (1 - expectedFractionTimeInactive)

		self.activationProb = effectiveFractionActive * expectedTerminationRate / (1 - effectiveFractionActive)

		# Views

		self.activeRnaPolys = self.uniqueMoleculesView('activeRnaPoly')

		self.inactiveRnaPolys = self.bulkMoleculeView("APORNAP-CPLX[c]")

		self.ppGpp = self.bulkMoleculeView("GUANOSINE-5DP-3DP[c]")
		self.ppGpp_base_conc = kb.process.metabolism.metabolitePoolConcentrations[kb.process.metabolism.metabolitePoolIDs.index("GUANOSINE-5DP-3DP[c]")]
		self.ppGpp_scaling_factor = 1 / (self.ppGpp_base_conc ** kb.ppGpp_power)
		
		###### VARIANT CODE #######
		# self.ppGppFeedback = kb.ppGppFeedback
		self.ppGppFeedback = True
		self.ppGpp_power = kb.ppGpp_power
		###### VARIANT CODE #######

	def calculateRequest(self):
		self.inactiveRnaPolys.requestAll()


	# Calculate temporal evolution
	def evolveState(self):
		# Sample a multinomial distribution of synthesis probabilities to 
		# determine what molecules are initialized

		inactiveRnaPolyCount = self.inactiveRnaPolys.count()

		rnaPolyToActivate = np.int64(self.activationProb * inactiveRnaPolyCount)

		if rnaPolyToActivate == 0:
			return

		# Scale synthesis probabilities of stable and unstable RNA by ppGpp concentration
		cellMass = (self.readFromListener("Mass", "cellMass") * units.fg)
		cellVolume = cellMass / self.cellDensity
		ppGpp_conc = (1 / self.nAvogadro) * (1 / cellVolume) * self.ppGpp.total()[0]

		###### VARIANT CODE #######
		if self.ppGppFeedback:
			stable_rna_scale = (1/(self.ppGpp_scaling_factor * (ppGpp_conc ** self.ppGpp_power))).normalize()
			stable_rna_scale.checkNoUnit()
			stable_rna_scale = np.fmin(1, stable_rna_scale.asNumber())
		else:
			stable_rna_scale = 1.
		###### VARIANT CODE #######

		scaledRnaSynthProb = self.rnaSynthProb.copy()
		scaledRnaSynthProb[self.tRNAIdx] = scaledRnaSynthProb[self.tRNAIdx] * stable_rna_scale
		scaledRnaSynthProb[self.rRNAIdx] = scaledRnaSynthProb[self.rRNAIdx] * stable_rna_scale
		scaledRnaSynthProb = scaledRnaSynthProb / scaledRnaSynthProb.sum()

		nNewRnas = self.randomState.multinomial(rnaPolyToActivate,
			scaledRnaSynthProb)

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


		activeRnaPolys = self.activeRnaPolys.moleculesNew(
			"activeRnaPoly",
			rnaPolyToActivate
			)

		activeRnaPolys.attrIs(
			rnaIndex = rnaIndexes
			)

		self.inactiveRnaPolys.countDec(nNewRnas.sum())
		
		mRnaInitalized = (self.mRNAIdx * nNewRnas * self.rnaLengths).sum()
		tRnaInitalized = (self.tRNAIdx * nNewRnas * self.rnaLengths).sum()
		rRnaInitalized = (self.rRNAIdx * nNewRnas * self.rnaLengths).sum()

		# Write initiation information to listner
		self.writeToListener("InitiatedTranscripts", "mRnaInitalized", mRnaInitalized)
		self.writeToListener("InitiatedTranscripts", "tRnaInitalized", tRnaInitalized)
		self.writeToListener("InitiatedTranscripts", "rRnaInitalized", rRnaInitalized)

		self.writeToListener("RnapData", "didInitialize", nNewRnas.sum())
