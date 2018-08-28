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

import copy
from itertools import izip
import re

import numpy as np

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
		self.water = self.bulkMoleculeView('WATER[c]')
		self.gtp = self.bulkMoleculeView("GTP[c]")
		self.proton = self.bulkMoleculeView("PROTON[c]")

		# Set for timestep calculation
		self.gtpUsed = 0
		self.gtpAvailable = 0

		self.ribosome30S = self.bulkMoleculeView(sim_data.moleculeIds.s30_fullComplex)
		self.ribosome50S = self.bulkMoleculeView(sim_data.moleculeIds.s50_fullComplex)

		self.translationSupply = sim._translationSupply

		self.elngRateFactor = 1.

		self.synthetase_names = []
		synthetase_rxns = [rxn for rxn in sim_data.process.metabolism.reactionCatalysts if re.findall('Charged\-(.*?)\-tRNAs', rxn)]
		for rxn in synthetase_rxns:
			for synthetase in sim_data.process.metabolism.reactionCatalysts[rxn]:
				if synthetase not in self.synthetase_names:
					self.synthetase_names.append(synthetase)

		# Build matrix to map AA and synthetases
		self.aa_from_synthetase = np.zeros((len(self.aaNames), len(self.synthetase_names)))
		for rxn in synthetase_rxns:
			aa = re.findall('Charged\-(.*?)\-tRNAs', rxn)[0]
			if aa == "ALA":
				aa = "L-ALPHA-ALANINE"
			elif aa == "ASP":
				aa = "L-ASPARTATE"

			for syn in sim_data.process.metabolism.reactionCatalysts[rxn]:
				aa_idx = self.aaNames.index(aa + "[c]")
				syn_idx = self.synthetase_names.index(syn)
				self.aa_from_synthetase[aa_idx, syn_idx] = 1

		self.charging_stoich_matrix = sim_data.process.transcription.charging_stoich_matrix()
		self.uncharged_trna_names = sim_data.process.transcription.rnaData['id'][sim_data.process.transcription.rnaData['isTRna']]
		self.charged_trna_names = sim_data.process.transcription.charged_trna_names
		self.charging_molecule_names = sim_data.process.transcription.charging_molecules
		self.uncharged_trna = self.bulkMoleculesView(self.uncharged_trna_names)
		self.charged_trna = self.bulkMoleculesView(self.charged_trna_names)
		self.charging_molecules = self.bulkMoleculesView(self.charging_molecule_names)
		self.aa_from_trna = sim_data.process.transcription.aa_from_trna
		self.synthetases = self.bulkMoleculesView(self.synthetase_names)

		# ppGpp parameters
		# TODO - put in flat file
		self.kS = 100.  # / units.s  # synthetase charging rate
		self.KMtf = 1.  # * uM  # Micahelis constant for synthetases and uncharged tRNAs
		self.KMaa = 100.  # * uM # Michaelis constant for synthetases and amino acids
		self.krib = 22.  # / units.s  # ribosome elongation rate
		self.krta = 1.  # * uM  # dissociation constant of charged tRNA-ribosome
		self.krtf = 500.  # * uM  # dissociation constant of uncharged tRNA-ribosome

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
			self.ribosomeElongationRate,
			)

		sequenceHasAA = (sequences != polymerize.PAD_VALUE)
		aasInSequences = np.bincount(sequences[sequenceHasAA], minlength=21)

		# Conversion from counts to molarity
		cell_mass = self.readFromListener("Mass", "cellMass") * units.fg
		cell_volume = cell_mass / self.cellDensity
		counts_to_molar = 1 / (self.nAvogadro * cell_volume)

		# Get counts and convert synthetase and tRNA to a per AA basis
		synthetase_counts = np.dot(self.aa_from_synthetase, self.synthetases.total())
		aa_counts = self.aas.total()
		uncharged_trna_counts = np.dot(self.aa_from_trna, self.uncharged_trna.total())
		charged_trna_counts = np.dot(self.aa_from_trna, self.charged_trna.total())
		ribosome_counts = len(self.activeRibosomes.allMolecules())

		# Get concentration
		mask = np.ones(21, dtype=bool)  # TODO - handle selenocysteine
		mask[-2] = False
		f = aasInSequences[mask] / np.sum(aasInSequences[mask])
		synthetase_conc = (counts_to_molar * synthetase_counts)[mask]
		aa_conc = (counts_to_molar * aa_counts)[mask]
		uncharged_trna_conc = (counts_to_molar * uncharged_trna_counts)[mask]
		charged_trna_conc = (counts_to_molar * charged_trna_counts)[mask]
		ribosome_conc = (counts_to_molar * ribosome_counts)

		# Calculate steady state tRNA levels and resulting elongation rate
		fraction_charged, v_rib = self.calculate_trna_charging(
			synthetase_conc, uncharged_trna_conc, charged_trna_conc, aa_conc, ribosome_conc, f
			)

		if self.translationSupply:
			translationSupplyRate = self.translation_aa_supply[current_nutrients] * self.elngRateFactor

			self.writeToListener("RibosomeData", "translationSupply", translationSupplyRate.asNumber())


			dryMass = (self.readFromListener("Mass", "dryMass") * units.fg)

			molAasRequested = translationSupplyRate * dryMass * self.timeStepSec() * units.s

			count_aas_requested = units.convertNoUnitToNumber(molAasRequested * self.nAvogadro)

			count_aas_requested = np.fmin(count_aas_requested, aasInSequences) # Check if this is required. It is a better request but there may be fewer elongations.
		else:
			count_aas_requested = aasInSequences

		charging_aa_request = v_rib * f * self._sim.timeStepSec() / counts_to_molar.asNumber(uM)

		updated_fraction_charged = np.zeros(len(self.aaNames))
		updated_fraction_charged[mask] = fraction_charged

		total_trna = self.charged_trna.total() + self.uncharged_trna.total()
		final_charged_trna = np.dot(updated_fraction_charged, self.aa_from_trna * total_trna)

		charged_trna_request = self.charged_trna.total() - final_charged_trna
		charged_trna_request[charged_trna_request < 0] = 0
		uncharged_trna_request = final_charged_trna - self.charged_trna.total()
		uncharged_trna_request[uncharged_trna_request < 0] = 0

		count_aas_requested[mask] = charging_aa_request
		self.aa_counts_for_translation = np.array(count_aas_requested)

		frac_trna_of_aa = total_trna / np.dot(np.dot(self.aa_from_trna, total_trna), self.aa_from_trna)
		total_charging_reactions = np.dot(count_aas_requested, self.aa_from_trna) * frac_trna_of_aa + uncharged_trna_request

		# Only request molecules that will be consumed in the charging reactions
		requested_molecules = -np.dot(self.charging_stoich_matrix, total_charging_reactions)
		requested_molecules[requested_molecules < 0] = 0
		self.charging_molecules.requestIs(requested_molecules)

		# Request charged tRNA that will become uncharged
		self.charged_trna.requestIs(charged_trna_request)

		# Request water for transfer of AA from tRNA for initial polypeptide
		# This is severe overestimate being worst case of every elongation is
		# to initialize a polypeptide but excess of water shouldn't matter
		self.water.requestIs(np.sum(charging_aa_request))

		self.writeToListener("GrowthLimits", "fraction_trna_charged", np.dot(updated_fraction_charged, self.aa_from_trna))
		self.writeToListener("GrowthLimits", "aaPoolSize", self.aas.total())
		self.writeToListener("GrowthLimits", "aaRequestSize", count_aas_requested)

		# Request GTP for polymerization based on sequences
		gtpsHydrolyzed = np.int64(np.ceil(self.gtpPerElongation * count_aas_requested.sum()))

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
			self.ribosomeElongationRate,
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
		aas_used = result.monomerUsages
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

		# get tRNA counts
		uncharged_trna = self.uncharged_trna.counts()
		charged_trna = self.charged_trna.counts()
		total_trna = uncharged_trna + charged_trna
		fraction_trna_per_aa = total_trna / np.dot(np.dot(self.aa_from_trna, total_trna), self.aa_from_trna)
		fraction_trna_per_aa[~np.isfinite(fraction_trna_per_aa)] = 0

		# number that can be charged
		aa_for_charging = total_aa_counts - aas_used
		n_aa_charged = np.fmin(aa_for_charging, np.dot(self.aa_from_trna, uncharged_trna))
		n_trna_charged = self.distribution_from_aa(n_aa_charged, fraction_trna_per_aa)
		net_charged = n_trna_charged - charged_trna

		charged_and_elongated = self.distribution_from_aa(aas_used, fraction_trna_per_aa)
		total_charging_reactions = charged_and_elongated + net_charged
		self.charging_molecules.countsInc(np.dot(self.charging_stoich_matrix, total_charging_reactions))

		# account for uncharging of tRNA during elongation
		self.charged_trna.countsDec(charged_and_elongated)
		self.uncharged_trna.countsInc(charged_and_elongated)

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

		# Update counts protons to reflect polymerization reactions and transfer of AA from tRNA
		# Peptide bond formation releases a water but transferring AA from tRNA consumes a OH-
		# Net production of H+ for each elongation, consume extra water for each initialization
		# since a peptide bond doesn't form
		self.proton.countInc(nElongations)
		self.water.countDec(nInitialized)

		# Report stalling information
		expectedElongations = np.fmin(
			self.ribosomeElongationRate,
			terminalLengths - peptideLengths
			)

		ribosomeStalls = expectedElongations - sequenceElongations

		# Write data to listeners
		self.writeToListener("GrowthLimits", "net_charged", net_charged)
		self.writeToListener("GrowthLimits", "aasUsed", aas_used)
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

	def calculate_trna_charging(self, synthetase_conc, uncharged_trna_conc, charged_trna_conc, aa_conc, ribosome_conc, f):
		'''
		Calculates the steady state value of tRNA based on charging and incorporation through polypeptide elongation.
		The fraction of charged/uncharged is also used to determine how quickly the ribosome is elongating.

		Inputs:
			synthetase_conc (array of floats with concentration units) - concentration of synthetases associated
				with each amino acid
			uncharged_trna_conc (array of floats with concentration units) - concentration of uncharged tRNA associated
				with each amino acid
			charged_trna_conc (array of floats with concentration units) - concentration of charged tRNA associated
				with each amino acid
			aa_conc (array of floats with concentration units) - concentration of each amino acid
			ribosome_conc (float with concentration units) - concentration of active ribosomes
			f (array of floats) - fraction of each amino acid to be incorporated to total amino acids incorporated

		Returns:
			fraction_charged (array of floats) - fraction of total tRNA that is charged for each tRNA species
			v_rib (float) - ribosomal elongation rate in units of uM/s
		'''

		synthetase_conc = synthetase_conc.asNumber(uM)
		uncharged_trna_conc = uncharged_trna_conc.asNumber(uM)
		charged_trna_conc = charged_trna_conc.asNumber(uM)
		aa_conc = aa_conc.asNumber(uM)
		ribosome_conc = ribosome_conc.asNumber(uM)

		## solve to steady state with short time steps
		dt = 0.001
		diff = 1
		while diff > 1e-3:
			v_charging = (self.kS * synthetase_conc * uncharged_trna_conc * aa_conc
				/ (self.KMaa * self.KMtf *
				(1 + uncharged_trna_conc/self.KMtf + aa_conc/self.KMaa + uncharged_trna_conc*aa_conc/self.KMtf/self.KMaa))
				)
			numerator_ribosome = 1 + np.sum(f * (self.krta / charged_trna_conc + uncharged_trna_conc / charged_trna_conc * self.krta / self.krtf))
			v_rib = self.krib*ribosome_conc / numerator_ribosome

			delta_conc = (v_charging - v_rib*f) * dt
			uncharged_trna_conc -= delta_conc
			charged_trna_conc += delta_conc
			diff = np.linalg.norm(delta_conc)

		fraction_charged = charged_trna_conc / (uncharged_trna_conc + charged_trna_conc)
		return fraction_charged, v_rib

	def distribution_from_aa(self, n_aa, f_trna):
		'''
		Distributes counts of amino acids to tRNAs that are associated with each amino acid
		Uses self.aa_from_trna mapping to distribute from amino acids to tRNA

		Inputs:
			n_aa (array of ints) - counts per amino acid to distribute to each tRNA
			f_trna (array of floats) - fraction of each tRNA to determine the distribution

		Returns:
			array of ints - distributed counts for each tRNA

		TODO:
			include limits on what can be charged/uncharged - might be necessary for some edge cases
		'''

		trna_counts = np.zeros(f_trna.shape, int)
		for i, row in enumerate(self.aa_from_trna):
			count = n_aa[i]
			idx = row == 1
			frac = f_trna[idx]

			counts = np.round(frac * count)
			diff = int(count - np.sum(counts))

			# Adjust if rounding doesn't lead to correct distribution
			# Only remove from ones that have been assigned if more assigned than available
			if diff < 0:
				for _ in range(-diff):
					new_frac = counts / np.sum(counts)
					adjustment = self.randomState.multinomial(1, new_frac)
					counts -= adjustment
			# Add additional counts to get up to counts to distribute
			elif diff > 0:
				adjustment = self.randomState.multinomial(diff, frac)
				counts += adjustment

			trna_counts[idx] = counts
			if np.any(counts < 0):
				import ipdb; ipdb.set_trace()

		return trna_counts

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
