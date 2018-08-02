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
import copy
import re

import wholecell.processes.process
from wholecell.utils.polymerize import buildSequences, polymerize, computeMassIncrease
from wholecell.utils.random import stochasticRound
from wholecell.utils import units

uM = units.umol / units.L


class PolypeptideElongation(wholecell.processes.process.Process):
	""" PolypeptideElongation """

	_name = "PolypeptideElongation"

	def __init__(self):
		super(PolypeptideElongation, self).__init__()

	def initialize(self, sim, sim_data):
		super(PolypeptideElongation, self).initialize(sim, sim_data)

		# Load parameters
		self.nAvogadro = sim_data.constants.nAvogadro
		self.cellDensity = sim_data.constants.cellDensity
		self.aaNames = sim_data.moleculeGroups.aaIDs
		proteinIds = sim_data.process.translation.monomerData['id']
		self.proteinLengths = sim_data.process.translation.monomerData["length"].asNumber()
		self.proteinSequences = sim_data.process.translation.translationSequences
		self.aaWeightsIncorporated = sim_data.process.translation.translationMonomerWeights
		self.endWeight = sim_data.process.translation.translationEndWeight
		self.gtpPerElongation = sim_data.constants.gtpPerTranslation
		# self.ribosomeElongationRate = float(sim_data.growthRateParameters.ribosomeElongationRate.asNumber(units.aa / units.s))

		self.maxRibosomeElongationRate = float(sim_data.constants.ribosomeElongationRateMax.asNumber(units.aa / units.s))

		self.ribosomeElongationRateDict = sim_data.process.translation.ribosomeElongationRateDict

		self.translation_aa_supply = sim_data.translationSupplyRate

		# Used for figure in publication
		self.trpAIndex = np.where(proteinIds == "TRYPSYN-APROTEIN[c]")[0][0]

		# Create view onto activly elongating 70S ribosomes
		self.activeRibosomes = self.uniqueMoleculesView('activeRibosome')

		# Create views onto 30S and 70S ribosomal subunits for termination
		self.ribosome30S = self.bulkMoleculeView(sim_data.moleculeIds.s30_fullComplex)
		self.ribosome50S = self.bulkMoleculeView(sim_data.moleculeIds.s50_fullComplex)

		# Create view onto all proteins
		self.bulkMonomers = self.bulkMoleculesView(proteinIds)

		# Create views onto all polymerization reaction small molecules
		self.aas = self.bulkMoleculesView(self.aaNames)
		self.h2o = self.bulkMoleculeView('WATER[c]')
		self.gtp = self.bulkMoleculeView("GTP[c]")
		self.gdp = self.bulkMoleculeView("GDP[c]")
		self.pi = self.bulkMoleculeView("PI[c]")
		self.h   = self.bulkMoleculeView("PROTON[c]")

		# Set for timestep calculation
		self.gtpUsed = 0
		self.gtpAvailable = 0

		self.ribosome30S = self.bulkMoleculeView(sim_data.moleculeIds.s30_fullComplex)
		self.ribosome50S = self.bulkMoleculeView(sim_data.moleculeIds.s50_fullComplex)

		self.translationSupply = sim._translationSupply

		self.elngRateFactor = 1.


		self.synthetaseNames = []
		synthetaseRxns = [rxn for rxn in sim_data.process.metabolism.reactionCatalysts if re.findall('Charged\-(.*?)\-tRNAs', rxn)]
		for rxn in synthetaseRxns:
			for syn in sim_data.process.metabolism.reactionCatalysts[rxn]:
				if syn not in self.synthetaseNames:
					self.synthetaseNames.append(syn)

		# Build matrix to map AA and synthetases
		self.aaFromSynthetase = np.zeros((len(self.aaNames), len(self.synthetaseNames)))
		for rxn in synthetaseRxns:
			aa = re.findall('Charged\-(.*?)\-tRNAs', rxn)[0]
			if aa == "ALA":
				aa = "L-ALPHA-ALANINE"
			elif aa == "ASP":
				aa = "L-ASPARTATE"

			for syn in sim_data.process.metabolism.reactionCatalysts[rxn]:
				aaIdx = self.aaNames.index(aa + "[c]")
				synIdx = self.synthetaseNames.index(syn)
				self.aaFromSynthetase[aaIdx, synIdx] = 1

		self.chargingStoichMatrix = sim_data.process.transcription.chargingStoichMatrix
		self.unchargedTrnaNames = sim_data.process.transcription.rnaData['id'][sim_data.process.transcription.rnaData['isTRna']]
		self.chargedTrnaNames = sim_data.process.transcription.chargedTrnaNames
		self.chargingMoleculeNames = sim_data.process.transcription.chargingMolecules
		self.uncharged_trna = self.bulkMoleculesView(self.unchargedTrnaNames)
		self.charged_trna = self.bulkMoleculesView(self.chargedTrnaNames)
		self.chargingMolecules = self.bulkMoleculesView(self.chargingMoleculeNames)
		self.aaFromTrna = sim_data.process.transcription.aaFromTrna
		self.synthetases = self.bulkMoleculesView(self.synthetaseNames)

	def calculateRequest(self):
		# Set ribosome elongation rate based on simulation medium environment and elongation rate factor
		# which is used to create single-cell variability in growth rate
		# The maximum number of amino acids that can be elongated in a single timestep is set to 22 intentionally as the minimum number of padding values
		# on the protein sequence matrix is set to 22. If timesteps longer than 1.0s are used, this feature will lead to errors in the effective ribosome
		# elongation rate.

		current_nutrients = self._external_states['Environment'].nutrients

		if self.translationSupply:
			self.ribosomeElongationRate = np.min([self.maxRibosomeElongationRate, int(stochasticRound(self.randomState,
				self.maxRibosomeElongationRate * self.timeStepSec()))]) # Will be set to maxRibosomeElongationRate if timeStepSec > 1.0s
		else:
			self.ribosomeElongationRate = np.min([22, int(stochasticRound(self.randomState,
				self.elngRateFactor * self.ribosomeElongationRateDict[current_nutrients].asNumber(units.aa / units.s) * self.timeStepSec()))])

		# Request all active ribosomes
		self.activeRibosomes.requestAll()

		activeRibosomes = self.activeRibosomes.allMolecules()

		if len(activeRibosomes) == 0:
			return

		# Build sequences to request appropriate amount of amino acids to
		# polymerize for next timestep
		proteinIndexes, peptideLengths = activeRibosomes.attrs(
					'proteinIndex', 'peptideLength'
					)

		sequences = buildSequences(
			self.proteinSequences,
			proteinIndexes,
			peptideLengths,
			self.ribosomeElongationRate
			)

		sequenceHasAA = (sequences != polymerize.PAD_VALUE)
		aasInSequences = np.bincount(sequences[sequenceHasAA], minlength=21)

		# ppGpp parameters
		# TODO - put in flat file
		kS = 100.  # / units.s  # synthetase charging rate
		KMtf = 1.  # * uM  # Micahelis constant for synthetases and uncharged tRNAs
		KMaa = 100.  # * uM # Michaelis constant for synthetases and amino acids
		krib = 22.  # / units.s  # ribosome elongation rate
		krta = 1.  # * uM  # dissociation constant of charged tRNA-ribosome
		krtf = 500.  # * uM  # dissociation constant of uncharged tRNA-ribosome

		# conversion from counts to molarity
		cell_mass = self.readFromListener("Mass", "cellMass") * units.fg
		cell_volume = cell_mass / self.cellDensity
		counts_to_molar = 1 / (self.nAvogadro * cell_volume)

		# get counts - convert synthetase and tRNA to a per AA basis
		synthetase_counts = np.dot(self.aaFromSynthetase, self.synthetases.total())
		aa_counts = self.aas.total()
		uncharged_trna_counts = np.dot(self.aaFromTrna, self.uncharged_trna.total())
		charged_trna_counts = np.dot(self.aaFromTrna, self.charged_trna.total()) + 1
		total_trna_counts = uncharged_trna_counts + charged_trna_counts
		ribosome_counts = len(self.activeRibosomes.allMolecules())

		# get concentration
		mask = np.ones(21, dtype=bool)
		mask[-2] = False
		f = aasInSequences[mask] / np.sum(aasInSequences[mask])
		synthetase_conc = (counts_to_molar * synthetase_counts).asNumber(uM)[mask]
		aa_conc = (counts_to_molar * aa_counts).asNumber(uM)[mask]
		uncharged_trna_conc = (counts_to_molar * uncharged_trna_counts).asNumber(uM)[mask]
		charged_trna_conc = (counts_to_molar * charged_trna_counts).asNumber(uM)[mask]
		total_trna_conc = (counts_to_molar * total_trna_counts).asNumber(uM)[mask]
		ribosome_conc = (counts_to_molar * ribosome_counts).asNumber(uM)

		## solve to steady state with short time steps
		# TODO - functionalize and add to initial conditions
		dt = 0.001
		diff = 1
		i = 0
		while diff > 1e-3:
			v_charging = (kS * synthetase_conc * uncharged_trna_conc * aa_conc
				/ (KMaa * KMtf *
				(1 + uncharged_trna_conc/KMtf + aa_conc/KMaa + uncharged_trna_conc*aa_conc/KMtf/KMaa))
				)
			numerator_ribosome = 1 + np.sum(f * krta / charged_trna_conc + uncharged_trna_conc / charged_trna_conc * krta / krtf)
			v_rib = krib*ribosome_conc / numerator_ribosome

			delta_conc = (v_charging - v_rib*f) * dt
			uncharged_trna_conc -= delta_conc
			charged_trna_conc += delta_conc
			diff = np.linalg.norm(delta_conc)
			if np.any(charged_trna_conc < 0) or np.any(uncharged_trna_conc < 0):
				import ipdb; ipdb.set_trace()
			i += 1
			# print diff
		print 'iterations: {}'.format(i)

		if self.translationSupply:
			translationSupplyRate = self.translation_aa_supply[current_nutrients] * self.elngRateFactor

			self.writeToListener("RibosomeData", "translationSupply", translationSupplyRate.asNumber())


			dryMass = (self.readFromListener("Mass", "dryMass") * units.fg)

			molAasRequested = translationSupplyRate * dryMass * self.timeStepSec() * units.s

			countAasRequested = units.convertNoUnitToNumber(molAasRequested * self.nAvogadro)

			countAasRequested = np.fmin(countAasRequested, aasInSequences) # Check if this is required. It is a better request but there may be fewer elongations.
		else:
			countAasRequested = aasInSequences

		# self.aas.requestIs(
		# 	countAasRequested
		# 	)

		charging_aa_request = v_rib * f * self._sim.timeStepSec() / counts_to_molar.asNumber(uM)
		supply_aa_request = countAasRequested[mask]
		fraction = charging_aa_request / supply_aa_request
		total_trna_conc = uncharged_trna_conc + charged_trna_conc
		fraction_charged = np.zeros(len(self.aaNames))
		fraction_charged[mask] = charged_trna_conc / total_trna_conc
		# print fraction
		print min(fraction), max(fraction)
		# print fraction_charged

		total_trna = self.charged_trna.total() + self.uncharged_trna.total()
		final_charged_trna = np.dot(fraction_charged, self.aaFromTrna * total_trna)

		charged_trna_request = self.charged_trna.total() - final_charged_trna
		charged_trna_request[charged_trna_request < 0] = 0
		uncharged_trna_request = final_charged_trna - self.charged_trna.total()
		uncharged_trna_request[uncharged_trna_request < 0] = 0

		countAasRequested[mask] = charging_aa_request
		self.aa_counts_for_translation = np.array(countAasRequested)
		countAasRequested += np.dot(self.aaFromTrna, uncharged_trna_request)

		self.aas.requestIs(countAasRequested)
		self.charged_trna.requestIs(charged_trna_request)
		self.uncharged_trna.requestIs(uncharged_trna_request)

		self.writeToListener("GrowthLimits", "aaPoolSize", self.aas.total())
		self.writeToListener("GrowthLimits", "aaRequestSize", countAasRequested)

		# Request GTP for polymerization based on sequences
		gtpsHydrolyzed = np.int64(np.ceil(self.gtpPerElongation * countAasRequested.sum()))

		self.writeToListener("GrowthLimits", "gtpPoolSize", self.gtp.total()[0])
		self.writeToListener("GrowthLimits", "gtpRequestSize", gtpsHydrolyzed)

		# GTP hydrolysis is carried out in Metabolism process for growth associated maintenence
		# THis is set here for metabolism to use
		self.gtpRequest = gtpsHydrolyzed

	def evolveState(self):
		# Write allocation data to listener
		self.writeToListener("GrowthLimits", "gtpAllocated", self.gtp.count())
		self.writeToListener("GrowthLimits", "aaAllocated", self.aas.counts())

		# Get number of active ribosomes
		activeRibosomes = self.activeRibosomes.molecules()

		self.writeToListener("GrowthLimits", "activeRibosomeAllocated", len(activeRibosomes))

		if len(activeRibosomes) == 0:
			return

		# Build amino acids sequences for each ribosome to polymerize
		proteinIndexes, peptideLengths, massDiffProtein = activeRibosomes.attrs(
			'proteinIndex', 'peptideLength', 'massDiff_protein'
			)

		sequences = buildSequences(
			self.proteinSequences,
			proteinIndexes,
			peptideLengths,
			self.ribosomeElongationRate
			)

		if sequences.size == 0:
			return

		# Calculate elongation resource capacity
		aaCountInSequence = np.bincount(sequences[(sequences != polymerize.PAD_VALUE)])
		total_aa_counts = self.aas.counts()
		aa_counts_for_translation = np.fmin(total_aa_counts, self.aa_counts_for_translation)

		# Using polymerization algorithm elongate each ribosome up to the limits
		# of amino acids, sequence, and GTP
		result = polymerize(
			sequences,
			aa_counts_for_translation,
			10000000, # Set to a large number, the limit is now taken care of in metabolism
			self.randomState
			)

		sequenceElongations = result.sequenceElongation
		aasUsed = result.monomerUsages
		nElongations = result.nReactions

		# Update masses of ribosomes attached to polymerizing polypeptides
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

		# number that can be charged
		# TODO - use charging reactions
		aa_for_charging = total_aa_counts - aasUsed
		n_aa_charged = np.fmin(aa_for_charging, np.dot(self.aaFromTrna, self.uncharged_trna.counts()))
		fraction_uncharged_trna_per_aa = self.uncharged_trna.counts() / np.dot(
			np.dot(self.aaFromTrna, self.uncharged_trna.counts()),
			self.aaFromTrna
			)
		n_trna_charged = np.dot(n_aa_charged, self.aaFromTrna) * fraction_uncharged_trna_per_aa
		net_charged = n_trna_charged - self.charged_trna.counts()
		net_charged[~np.isfinite(net_charged)] = 0

		self.charged_trna.countsInc(net_charged)
		self.uncharged_trna.countsDec(net_charged)
		self.aas.countsDec(np.dot(self.aaFromTrna, net_charged))

		# Write current average elongation to listener
		currElongRate = (sequenceElongations.sum() / len(activeRibosomes)) / self.timeStepSec()
		self.writeToListener("RibosomeData", "effectiveElongationRate", currElongRate)

		# Update active ribosomes, terminating if neccessary
		activeRibosomes.attrIs(
			peptideLength = updatedLengths,
			massDiff_protein = updatedMass
			)

		# Ribosomes that reach the end of their sequences are terminated and
		# dissociated into 30S and 50S subunits. The polypeptide that they are polymerizing
		# is converted into a protein in BulkMolecules
		terminalLengths = self.proteinLengths[proteinIndexes]

		didTerminate = (updatedLengths == terminalLengths)

		terminatedProteins = np.bincount(
			proteinIndexes[didTerminate],
			minlength = self.proteinSequences.shape[0]
			)

		activeRibosomes.delByIndexes(np.where(didTerminate)[0])
		self.bulkMonomers.countsInc(terminatedProteins)

		nTerminated = didTerminate.sum()
		nInitialized = didInitialize.sum()

		self.ribosome30S.countInc(nTerminated)
		self.ribosome50S.countInc(nTerminated)

		# Update counts of amino acids and water to reflect polymerization reactions
		self.aas.countsDec(aasUsed)
		self.h2o.countInc(nElongations - nInitialized)

		# Report stalling information
		expectedElongations = np.fmin(
			self.ribosomeElongationRate,
			terminalLengths - peptideLengths
			)

		ribosomeStalls = expectedElongations - sequenceElongations

		# Write data to listeners
		self.writeToListener("GrowthLimits", "aasUsed", aasUsed)
		self.writeToListener("GrowthLimits", "gtpUsed", self.gtpUsed)

		self.writeToListener("RibosomeData", "ribosomeStalls", ribosomeStalls)
		self.writeToListener("RibosomeData", "aaCountInSequence", aaCountInSequence)
		self.writeToListener("RibosomeData", "aaCounts", aa_counts_for_translation)

		self.writeToListener("RibosomeData", "expectedElongations", expectedElongations.sum())
		self.writeToListener("RibosomeData", "actualElongations", sequenceElongations.sum())
		self.writeToListener("RibosomeData", "actualElongationHist", np.histogram(sequenceElongations, bins = np.arange(0,23))[0])
		self.writeToListener("RibosomeData", "elongationsNonTerminatingHist", np.histogram(sequenceElongations[~didTerminate], bins=np.arange(0,23))[0])

		self.writeToListener("RibosomeData", "didTerminate", didTerminate.sum())
		self.writeToListener("RibosomeData", "terminationLoss", (terminalLengths - peptideLengths)[didTerminate].sum())
		self.writeToListener("RibosomeData", "numTrpATerminated", terminatedProteins[self.trpAIndex])

		self.writeToListener("RibosomeData", "processElongationRate", self.ribosomeElongationRate / self.timeStepSec())

	def isTimeStepShortEnough(self, inputTimeStep, timeStepSafetyFraction):
		"""
		Assumes GTP is the readout for failed translation with respect to the timestep.
		"""

		# Until more padding values are added to the protein sequence matrix, limit the maximum timestep length to 1 second
		# Since the current upper limit on a.a's elongated by ribosomes during a single timestep is set to 22, timesteps
		# longer than 1.0s do not lead to errors, but does slow down the ribosome elongation rate of the resulting simulation.
		# Must be modified if timesteps longer than 1.0s are desired.
		if inputTimeStep > 1.0:
			return False

		activeRibosomes = float(self.activeRibosomes.total()[0])
		self.gtpAvailable = float(self.gtp.total()[0])

		# Without an estimate on ribosome counts, require a short timestep until estimates available
		if activeRibosomes == 0:
			if inputTimeStep <= .2:
				return True
			else:
				return False

		dt = inputTimeStep * timeStepSafetyFraction
		gtpExpectedUsage = activeRibosomes * self.ribosomeElongationRate * self.gtpPerElongation * dt

		if gtpExpectedUsage < self.gtpAvailable:
			return True
		else:
			return False

	def wasTimeStepShortEnough(self):
		"""
		If translation used more than 90 percent of gtp, timeStep was too short.
		"""

		# If gtpAvailable is 0 and the timeStep is short, use the gtp produced this timeStep as the estimate
		if self.gtpAvailable == 0 and self.timeStepSec() <= .2:
			self.gtpAvailable = self.gtp.total()[0]

		if (self.gtpAvailable * .9) < self.gtpUsed:
			return False
		else:
			return True
