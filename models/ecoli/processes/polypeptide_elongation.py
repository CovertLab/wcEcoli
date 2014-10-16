#!/usr/bin/env python

"""
PolypeptideElongation

Translation elongation sub-model.

TODO:
- see the initiation process for more TODOs

@author: Derek Macklin
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 4/30/14
"""

from __future__ import division

from itertools import izip

import numpy as np

import wholecell.processes.process
from wholecell.utils.polymerize import buildSequences, polymerize, computeMassIncrease, PAD_VALUE
from wholecell.utils.random import stochasticRound
from wholecell.utils import units

class PolypeptideElongation(wholecell.processes.process.Process):
	""" PolypeptideElongation """

	_name = "PolypeptideElongation"

	# Constructor
	def __init__(self):
		# Parameters
		self.elngRate = None
		self.proteinLengths = None
		self.proteinSequences = None
		self.proteinSequencesNascent = None
		self.h2oWeight = None
		self.aaWeightsIncorporated = None
		self.gtpPerElongation = None
		self.synthetase_turnover = None

		# Views
		self.activeRibosomes = None
		self.bulkMonomers = None
		self.bulkMonomersNascent = None
		self.aas = None
		self.h2o = None
		self.trna_groups = None
		self.synthetase_groups = None
		self.ribosome30S = None
		self.ribosome50S = None

		super(PolypeptideElongation, self).__init__()


	# Construct object graph
	def initialize(self, sim, kb):
		super(PolypeptideElongation, self).initialize(sim, kb)

		# Load parameters

		self.elngRate = float(kb.ribosomeElongationRate.asNumber(units.aa / units.s)) * self.timeStepSec

		self.aa_trna_groups = kb.aa_trna_groups
		self.aa_synthetase_groups = kb.aa_synthetase_groups
		self.synthetase_turnover = kb.trna_synthetase_rates.asNumber(units.aa/units.s)

		enzIds = ["RRLA-RRNA[c]", "RRSA-RRNA[c]", "RRFA-RRNA[c]"]

		proteinIds = kb.proteinData['id']
		proteinNascentIds = kb.proteinNascentData['id']

		self.proteinLengths = kb.proteinData["length"].asNumber()
		self.proteinLengthsNascent = kb.proteinNascentData["length"].asNumber()

		#self.proteinSequences = kb.translationSequences
		self.proteinSequencesNascent = kb.translationSequencesNascent

		self.aaWeightsIncorporated = kb.translationMonomerWeights

		self.endWeight = kb.translationEndWeight

		self.gtpPerElongation = kb.gtpPerTranslation

		# Views

		self.activeRibosomes = self.uniqueMoleculesView('activeRibosome')
		self.bulkMonomers = self.bulkMoleculesView(proteinIds)
		self.bulkMonomersNascent = self.bulkMoleculesView(proteinNascentIds)

		self.aas = self.bulkMoleculesView(kb.aaIDs)
		self.trna_groups = [self.bulkMoleculesView(x) for x in self.aa_trna_groups.itervalues()]
		self.synthetase_groups = [self.bulkMoleculesView(x) for x in self.aa_synthetase_groups.itervalues()]
		self.h2o = self.bulkMoleculeView('H2O[c]')

		self.gtp = self.bulkMoleculeView("GTP[c]")
		self.gdp = self.bulkMoleculeView("GDP[c]")
		self.pi = self.bulkMoleculeView("PI[c]")
		self.h   = self.bulkMoleculeView("H[c]")

		self.ribosome30S = self.bulkMoleculeView(kb.s30_fullComplex)
		self.ribosome50S = self.bulkMoleculeView(kb.s50_fullComplex)

		self.mtf = self.bulkMoleculeView("EG11268-MONOMER[c]")
		self.mtfKcat = kb.mtfKcat
		self.fthf = self.bulkMoleculeView("10FTHF[c]")
		self.numberAAs = len(kb.aaIDs)
		self.formate = self.bulkMoleculeView("FOR[c]")

	def calculateRequest(self):
		self.activeRibosomes.requestAll()

		activeRibosomes = self.activeRibosomes.allMolecules()

		proteinIndexes, peptideLengths = activeRibosomes.attrs(
			'proteinIndex', 'peptideLength'
			)

		sequences = buildSequences(
			self.proteinSequencesNascent,
			proteinIndexes,
			peptideLengths,
			self.elngRate
			)

		sequenceHasAA = (sequences != PAD_VALUE)

		aasRequested = np.bincount(sequences[sequenceHasAA],minlength=self.numberAAs)

		self.aas.requestIs(
			aasRequested
			)

		# Should essentially request all tRNAs
		# and all synthetases
		trnasRequested = aasRequested
		for i,group in enumerate(self.trna_groups):
			group.requestIs(trnasRequested[i])
		synthetaseRequested = aasRequested
		for i,group in enumerate(self.synthetase_groups):
			group.requestIs(synthetaseRequested[i])

		gtpsHydrolyzed = np.int64(np.ceil(
			self.gtpPerElongation * np.fmin(
				sequenceHasAA.sum(),
				self.aas.total().sum()
				)
			))

		self.gtp.requestIs(gtpsHydrolyzed)

		self.h2o.requestIs(gtpsHydrolyzed) # note: this is roughly a 2x overestimate

		self.mtf.requestAll()
		self.fthf.requestAll()

	# Calculate temporal evolution
	def evolveState(self):
		activeRibosomes = self.activeRibosomes.molecules()

		if len(activeRibosomes) == 0:
			return

		proteinIndexes, peptideLengths, massDiffProtein = activeRibosomes.attrs(
			'proteinIndex', 'peptideLength', 'massDiff_protein'
			)

		# Build sequence array

		sequences = buildSequences(
			self.proteinSequencesNascent,
			proteinIndexes,
			peptideLengths,
			self.elngRate
			)

		# Calculate elongation resource capacity
		aaCountInSequence = np.bincount(sequences[(sequences != PAD_VALUE)],minlength=self.numberAAs)

		metIndex=12
		fmetIndex=21
		Nmet=aaCountInSequence[metIndex]
		Nfmet=aaCountInSequence[fmetIndex]
		aaCounts = self.aas.counts()
		metCount = aaCounts[metIndex]
		metCountAllocated = metCount * (Nmet/(Nmet+Nfmet))
		fmetCountAllocated = metCount * (Nfmet/(Nmet+Nfmet))
		aaCountsUpdated = aaCounts
		aaCountsUpdated[metIndex]=metCountAllocated
		aaCountsUpdated[fmetIndex]=fmetCountAllocated
		synthetaseCounts = np.array([x.counts().sum() for x in self.synthetase_groups],dtype = np.int64)
		metSynthetaseCount = synthetaseCounts[metIndex]
		metSynthetaseCountAllocated = metSynthetaseCount * (Nmet/(Nmet+Nfmet))
		fmetSynthetaseCountAllocated = metSynthetaseCount * (Nfmet/(Nmet+Nfmet))
		synthetaseCountsUpdated = synthetaseCounts
		synthetaseCountsUpdated[metIndex] = metSynthetaseCountAllocated
		synthetaseCountsUpdated[fmetIndex] = fmetSynthetaseCountAllocated
		mtfCapacity = self.mtfKcat * self.mtf.count()
		trnasCapacity = self.synthetase_turnover * np.array([x.counts().sum() for x in self.trna_groups],dtype = np.int64)
		#synthetaseCapacity = self.synthetase_turnover * np.array([x.counts().sum() for x in self.synthetase_groups],dtype = np.int64)
		synthetaseCapacity = self.synthetase_turnover * synthetaseCountsUpdated
		#elongationResourceCapacity = np.minimum(aaCounts, synthetaseCapacity, trnasCapacity) #np.minimum does not take the min of 3 arrays!
		elongationResourceCapacity = np.min([aaCountsUpdated, synthetaseCapacity, trnasCapacity],axis=0)
		elongationResourceCapacity[fmetIndex] = np.min([aaCountsUpdated[fmetIndex], synthetaseCapacity[fmetIndex], 
							trnasCapacity[fmetIndex], mtfCapacity, self.fthf.count()])

		# Calculate update

		reactionLimit = self.gtp.count() // self.gtpPerElongation

		sequenceElongations, aasUsed, nElongations = polymerize(
			sequences,
#			aaCounts, # elongationResourceCapacity,
			elongationResourceCapacity,
			reactionLimit,
			self.randomState
			)

		#FMet aasUsed are actually Met aasUsed
		aasUsed[metIndex]=aasUsed[metIndex]+aasUsed[fmetIndex]
		aasUsed[fmetIndex]=0

		massIncreaseProtein = computeMassIncrease(
			sequences,
			sequenceElongations,
			self.aaWeightsIncorporated
			)

		updatedLengths = peptideLengths + sequenceElongations

		didInitialize = (
			(sequenceElongations > 1) &
			(peptideLengths == 0)
			)

		updatedMass = massDiffProtein + massIncreaseProtein

		updatedMass[didInitialize] += self.endWeight

		# Update active ribosomes, terminating if neccessary

		activeRibosomes.attrIs(
			peptideLength = updatedLengths,
			massDiff_protein = updatedMass
			)

		terminalLengths = self.proteinLengthsNascent[proteinIndexes]

		didTerminate = (updatedLengths == terminalLengths)

		terminatedProteins = np.bincount(
			proteinIndexes[didTerminate],
			minlength = self.proteinSequencesNascent.shape[0]
			)

		activeRibosomes.delByIndexes(np.where(didTerminate)[0])

		nTerminated = didTerminate.sum()
		nInitialized = didInitialize.sum()

		# Update bulk molecules

		self.aas.countsDec(aasUsed)

		self.bulkMonomersNascent.countsInc(terminatedProteins)

		self.ribosome30S.countInc(nTerminated)
		self.ribosome50S.countInc(nTerminated)

		self.h2o.countInc(nElongations - nInitialized)

		gtpUsed = np.int64(stochasticRound(
			self.randomState, 
			nElongations * self.gtpPerElongation
			))

		self.gtp.countDec(gtpUsed)
		self.gdp.countInc(gtpUsed)
		self.pi.countInc(gtpUsed)
		self.h.countInc(gtpUsed)
		self.fthf.countDec(nInitialized)
		self.formate.countInc(nInitialized)

		self.h2o.countDec(gtpUsed)

		# Report stalling information

		expectedElongations = np.fmin(
			self.elngRate,
			terminalLengths - peptideLengths
			)

		ribosomeStalls = expectedElongations - sequenceElongations

		self.writeToListener("RibosomeStalling", "ribosomeStalls", ribosomeStalls)
		self.writeToListener("RibosomeStalling", "aaCountInSequence", aaCountInSequence)
		self.writeToListener("RibosomeStalling", "aaCounts", aaCounts)
		self.writeToListener("RibosomeStalling", "trnasCapacity", trnasCapacity)
		self.writeToListener("RibosomeStalling", "synthetaseCapacity", synthetaseCapacity)
