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

SYNTHETASE_KM_SCALE = 0.1

class PolypeptideElongation(wholecell.processes.process.Process):
	""" PolypeptideElongation """

	_name = "PolypeptideElongation"

	# Constructor
	def __init__(self):
		super(PolypeptideElongation, self).__init__()


	# Construct object graph
	def initialize(self, sim, kb):
		super(PolypeptideElongation, self).initialize(sim, kb)

		# Load parameters
		self.elngRate = float(kb.constants.ribosomeElongationRate.asNumber(units.aa / units.s)) * self.timeStepSec
		self.nAvogadro = kb.constants.nAvogadro
		self.cellDensity = kb.constants.cellDensity
		self.aa_trna_groups = kb.process.translation.AA_TRNA_GROUPS
		self.aa_synthetase_groups = kb.process.translation.AA_SYNTHETASE_GROUPS
		self.synthetase_turnover = kb.trna_synthetase_rates.asNumber(units.aa/units.s)

		proteinIds = kb.process.translation.monomerData['id']
		self.proteinLengths = kb.process.translation.monomerData["length"].asNumber()
		self.proteinSequences = kb.process.translation.translationSequences
		self.aaWeightsIncorporated = kb.process.translation.translationMonomerWeights
		self.endWeight = kb.process.translation.translationEndWeight
		self.gtpPerElongation = kb.constants.gtpPerTranslation

		##########
		# Setting synthetase Km's to be fraction of steady state amino acid concentrations
		aaIdxs = [kb.process.metabolism.metabolitePoolIDs.index(aaID) for aaID in kb.moleculeGroups.aaIDs]
		aaConcentrations = kb.process.metabolism.metabolitePoolConcentrations[aaIdxs]
		self.synthetase_km = SYNTHETASE_KM_SCALE * aaConcentrations
		##########

		# Views
		self.activeRibosomes = self.uniqueMoleculesView('activeRibosome')
		self.bulkMonomers = self.bulkMoleculesView(proteinIds)

		self.aas = self.bulkMoleculesView(kb.moleculeGroups.aaIDs)
		self.trna_groups = [self.bulkMoleculesView(x) for x in self.aa_trna_groups.itervalues()]
		self.synthetase_groups = [self.bulkMoleculesView(x) for x in self.aa_synthetase_groups.itervalues()]

		self.h2o = self.bulkMoleculeView('H2O[c]')

		self.gtp = self.bulkMoleculeView("GTP[c]")
		self.gdp = self.bulkMoleculeView("GDP[c]")
		self.pi = self.bulkMoleculeView("PI[c]")
		self.h   = self.bulkMoleculeView("H[c]")

		self.ppgpp = self.bulkMoleculeView("PPGPP[c]")
		self.atp = self.bulkMoleculeView("ATP[c]")
		self.amp = self.bulkMoleculeView("AMP[c]")

		self.ribosome30S = self.bulkMoleculeView(kb.moleculeGroups.s30_fullComplex[0])
		self.ribosome50S = self.bulkMoleculeView(kb.moleculeGroups.s50_fullComplex[0])

	def calculateRequest(self):
		self.activeRibosomes.requestAll()

		activeRibosomes = self.activeRibosomes.allMolecules()

		if len(activeRibosomes) == 0:
			return

		proteinIndexes, peptideLengths = activeRibosomes.attrs(
			'proteinIndex', 'peptideLength'
			)

		sequences = buildSequences(
			self.proteinSequences,
			proteinIndexes,
			peptideLengths,
			self.elngRate
			)

		sequenceHasAA = (sequences != PAD_VALUE)

		aasInSequences = np.bincount(sequences[sequenceHasAA])

		# Should essentially request all tRNAs
		# and all synthetases
		trnasRequested = aasInSequences
		for i,group in enumerate(self.trna_groups):
			group.requestIs(trnasRequested[i])
		synthetaseRequested = aasInSequences
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

		# Calculate potential stalls with heuristic
		cellMass = (self.readFromListener("Mass", "cellMass") * units.fg)
		cellVolume = cellMass / self.cellDensity
		aaTotalConc = (1 / self.nAvogadro) * (1 / cellVolume) * self.aas.total()
		synthetaseCapacity = self.synthetase_turnover * np.array([x.total().sum() for x in self.synthetase_groups],dtype = np.int64)
		synthetaseSaturation = (aaTotalConc / (self.synthetase_km + aaTotalConc)).normalize()
		synthetaseSaturation.checkNoUnit()
		stallsPerAA = np.fmax(aasInSequences - synthetaseCapacity * synthetaseSaturation.asNumber(),0)
		totalStalls = np.ceil(stallsPerAA.sum())

		self.atp.requestIs(totalStalls)
		self.gdp.requestIs(totalStalls)

		#########	
		# Limiting request by synthetase capacity * saturation	
		self.aas.requestIs(
			np.floor(np.fmin(aasInSequences, synthetaseCapacity * synthetaseSaturation.asNumber()))
			)
		########

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
			self.proteinSequences,
			proteinIndexes,
			peptideLengths,
			self.elngRate
			)

		# Calculate elongation resource capacity

		aaCountInSequence = np.bincount(sequences[(sequences != PAD_VALUE)])
		aaCounts = self.aas.counts()
		trnasCapacity = self.synthetase_turnover * np.array([x.counts().sum() for x in self.trna_groups],dtype = np.int64)
		synthetaseCapacity = self.synthetase_turnover * np.array([x.counts().sum() for x in self.synthetase_groups],dtype = np.int64)
		elongationResourceCapacity = np.minimum(aaCounts, synthetaseCapacity, trnasCapacity)

		# Calculate expected stalls huristic
		cellMass = (self.readFromListener("Mass", "cellMass") * units.fg)
		cellVolume = cellMass / self.cellDensity
		aaTotalConc = (1 / self.nAvogadro) * (1 / cellVolume) * self.aas.total()
		synthetaseSaturation = (aaTotalConc / (self.synthetase_km + aaTotalConc)).normalize()
		synthetaseSaturation.checkNoUnit()
		synthetaseSaturation = synthetaseSaturation.asNumber()
		stallsPerAA = np.fmax(aaCountInSequence - synthetaseCapacity * synthetaseSaturation,0)
		totalStalls = np.ceil(stallsPerAA.sum())

		# Calculate update
		reactionLimit = self.gtp.count() // self.gtpPerElongation

		sequenceElongations, aasUsed, nElongations = polymerize(
			sequences,
			elongationResourceCapacity,
			reactionLimit,
			self.randomState
			)

		massIncreaseProtein = computeMassIncrease(
			sequences,
			sequenceElongations,
			self.aaWeightsIncorporated
			)

		updatedLengths = peptideLengths + sequenceElongations

		didInitialize = (
			(sequenceElongations > 0) &
			(peptideLengths == 0)
			)

		updatedMass = massDiffProtein + massIncreaseProtein

		updatedMass[didInitialize] += self.endWeight

		# Update active ribosomes, terminating if neccessary

		activeRibosomes.attrIs(
			peptideLength = updatedLengths,
			massDiff_protein = updatedMass
			)

		terminalLengths = self.proteinLengths[proteinIndexes]

		didTerminate = (updatedLengths == terminalLengths)

		terminatedProteins = np.bincount(
			proteinIndexes[didTerminate],
			minlength = self.proteinSequences.shape[0]
			)

		activeRibosomes.delByIndexes(np.where(didTerminate)[0])

		nTerminated = didTerminate.sum()
		nInitialized = didInitialize.sum()

		# Update bulk molecules

		self.aas.countsDec(aasUsed)

		self.bulkMonomers.countsInc(terminatedProteins)

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

		self.h2o.countDec(gtpUsed)

		usedCounts = np.array([self.atp.count(), self.gtp.count()])
		if (totalStalls > usedCounts).any():
			print 'Actual stalls greater than computed during request! Off by {}.'.format((usedCounts - totalStalls)[(usedCounts - totalStalls) < 0])
			totalStalls = usedCounts.min()
		
		self.atp.countDec(totalStalls)
		self.gdp.countDec(totalStalls)
		self.amp.countInc(totalStalls)
		self.ppgpp.countInc(totalStalls)
		self.h.countInc(totalStalls)

		# Report stalling information
		expectedElongations = np.fmin(
			self.elngRate,
			terminalLengths - peptideLengths
			)

		ribosomeStalls = expectedElongations - sequenceElongations

		self.writeToListener("RibosomeData", "ribosomeStalls", ribosomeStalls)
		self.writeToListener("RibosomeData", "aaCountInSequence", aaCountInSequence)
		self.writeToListener("RibosomeData", "aaCounts", aaCounts)
		# self.writeToListener("RibosomeData", "trnasCapacity", trnasCapacity)
		# self.writeToListener("RibosomeData", "synthetaseCapacity", synthetaseCapacity)

		self.writeToListener("GrowthRateControl", "totalStalls", totalStalls)
		self.writeToListener("GrowthRateControl", "synthetaseSaturation", synthetaseSaturation)

		self.writeToListener("RibosomeData", "expectedElongations", expectedElongations.sum())
		self.writeToListener("RibosomeData", "actualElongations", sequenceElongations.sum())
