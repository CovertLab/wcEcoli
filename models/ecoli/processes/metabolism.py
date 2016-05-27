#!/usr/bin/env python

"""
Metabolism

Metabolism sub-model. Encodes molecular simulation of microbial metabolism using flux-balance analysis.

TODO:
- enzyme-limited reactions (& fit enzyme expression)
- option to call a reduced form of metabolism (assume optimal)

@author: Derek Macklin
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 4/2/2013
"""

from __future__ import division

from itertools import izip

import numpy as np

import wholecell.processes.process
from wholecell.utils import units

from wholecell.utils.random import stochasticRound
from wholecell.utils.constants import REQUEST_PRIORITY_METABOLISM

from wholecell.utils.modular_fba import FluxBalanceAnalysis
from wholecell.utils.enzymeKinetics import EnzymeKinetics
from wholecell.utils.fitting import massesAndCountsToAddForPools

COUNTS_UNITS = units.dmol
VOLUME_UNITS = units.L
MASS_UNITS = units.g
TIME_UNITS = units.s

NONZERO_ENZYMES = True

SECRETION_PENALTY_COEFF = 1e-5

USE_RATELIMITS = False # Enable/disable kinetic rate limits in the model

USE_MANUAL_FLUX_COEFF = True # enable to overrid flux coefficients in the knowledgebase and use these local values instead
MAX_FLUX_COEFF = 1 # Multiple of predicted rate at which to set the max fluxes
MIN_FLUX_COEFF = 0 # Multiple of predicted rate at which to set the min fluxes


class Metabolism(wholecell.processes.process.Process):
	""" Metabolism """

	_name = "Metabolism"

	# Constructor
	def __init__(self):

		super(Metabolism, self).__init__()

	# Construct object graph
	def initialize(self, sim, sim_data):
		super(Metabolism, self).initialize(sim, sim_data)

		# Load constants
		self.nAvogadro = sim_data.constants.nAvogadro
		self.cellDensity = sim_data.constants.cellDensity

		self.metabolitePoolIDs = sorted(sim_data.process.metabolism.concDict)

		self.environment = sim_data.environment
		self.exchangeConstraints = sim_data.process.metabolism.exchangeConstraints

		self.doublingTime = sim_data.doubling_time

		# Load enzyme kinetic rate information
		self.reactionRateInfo = sim_data.process.metabolism.reactionRateInfo
		self.enzymeNames = sim_data.process.metabolism.enzymeNames
		self.constraintIDs = sim_data.process.metabolism.constraintIDs
		self.constraintMultiplesDict = sim_data.process.metabolism.constraintMultiplesDict
		self.constraintToReactionDict = sim_data.process.metabolism.constraintToReactionDict

		if USE_MANUAL_FLUX_COEFF:
			self.max_flux_coefficient = MAX_FLUX_COEFF
			self.min_flux_coefficient = MIN_FLUX_COEFF
		else:
			self.max_flux_coefficient = sim_data.constants.kineticRateLimitFactorUpper
			self.min_flux_coefficient = sim_data.constants.kineticRateLimitFactorLower

		self.objective = dict(
			(key, sim_data.process.metabolism.concDict[key].asNumber(COUNTS_UNITS / VOLUME_UNITS)) for key in sim_data.process.metabolism.concDict
			)

		# TODO: make sim_data method?
		extIDs = sim_data.externalExchangeMolecules[sim_data.environment]
		self.extMoleculeMasses = sim_data.getter.getMass(extIDs).asNumber(MASS_UNITS/COUNTS_UNITS) # TODO: delete this line?

		self.getMass = sim_data.getter.getMass
		self.massReconstruction = sim_data.mass
		self.avgCellToInitialCellConvFactor = sim_data.mass.avgCellToInitialCellConvFactor

		self.moleculeMasses = dict(zip(
			extIDs,
			self.getMass(extIDs).asNumber(MASS_UNITS/COUNTS_UNITS)
			))

		self.ngam = sim_data.constants.nonGrowthAssociatedMaintenance

		initWaterMass = sim_data.mass.avgCellWaterMassInit
		initDryMass = sim_data.mass.avgCellDryMassInit

		initCellMass = (
			initWaterMass
			+ initDryMass
			)

		self.energyCostPerWetMass = sim_data.constants.darkATP * initDryMass / initCellMass

		self.reactionStoich = sim_data.process.metabolism.reactionStoich
		self.externalExchangeMolecules = sim_data.externalExchangeMolecules[sim_data.environment]
		self.reversibleReactions = sim_data.process.metabolism.reversibleReactions

		# Set up FBA solver
		self.fba_object_options = {
			"reactionStoich" : self.reactionStoich.copy(), # TODO: copy in class
			"externalExchangedMolecules" : self.externalExchangeMolecules,
			"objective" : self.objective,
			"objectiveType" : "pools",
			"reversibleReactions" : self.reversibleReactions,
			"moleculeMasses" : self.moleculeMasses,
			"secretionPenaltyCoeff" : 0., # The "inconvenient constant"--limit secretion (e.g., of CO2); a value of 1e-5 seems to work
			"solver" : "glpk",
			"maintenanceCostGAM" : energyCostPerWetMass.asNumber(COUNTS_UNITS / MASS_UNITS),
			"maintenanceReaction" : {
				"ATP[c]": -1, "WATER[c]": -1, "ADP[c]": +1, "Pi[c]": +1, "PROTON[c]": +1,
				} # TODO: move to KB TODO: check reaction stoich
		}

		self.fba = FluxBalanceAnalysis(**self.fba_object_options)

		# Set up enzyme kinetics object
		self.enzymeKinetics = EnzymeKinetics(
			reactionRateInfo = sim_data.process.metabolism.reactionRateInfo,
			noCustoms=True
			)

		# Build dictionary from reaction name to index in the reaction array
		self.reactionNameToReactionIndexDict = {}
		for idx, reactionID in enumerate(self.fba.reactionIDs()):
			self.reactionNameToReactionIndexDict[reactionID] = idx

		# Determine which kinetic limits to use
		self.reactionsWithKineticLimits = [True]*len(self.fba.reactionIDs())
		self.constraintMultiplesDict = [self.constraintMultiplesDict[x] for x in self.constraintIDs]

		# Set constraints
		## External molecules
		self.externalMoleculeIDs = self.fba.externalMoleculeIDs()

		## Set enzymes unlimited
		self.fba.enzymeLevelsIs(np.inf)

		# Views
		self.metaboliteNames = self.fba.outputMoleculeIDs()
		self.metabolites = self.bulkMoleculesView(self.metaboliteNames)
		self.poolMetabolites = self.bulkMoleculesView(self.metabolitePoolIDs)
		self.enzymes = self.bulkMoleculesView(self.enzymeNames)

			
		outputMoleculeIDs = self.fba.outputMoleculeIDs()

		assert outputMoleculeIDs == self.fba.internalMoleculeIDs()

		# Set the priority to a low value
		self.bulkMoleculesRequestPriorityIs(REQUEST_PRIORITY_METABOLISM)


	def calculateRequest(self):
		self.metabolites.requestAll()
		self.enzymes.requestAll()

	# Calculate temporal evolution
	def evolveState(self):

		# Solve for metabolic fluxes
		metaboliteCountsInit = self.metabolites.counts()

		cellMass = (self.readFromListener("Mass", "cellMass") * units.fg)
		dryMass = (self.readFromListener("Mass", "dryMass") * units.fg)

		cellVolume = cellMass / self.cellDensity

		countsToMolar = 1 / (self.nAvogadro * cellVolume)

		polypeptideElongationEnergy = countsToMolar * 0
		if hasattr(self._sim.processes["PolypeptideElongation"], "gtpRequest"):
			polypeptideElongationEnergy = countsToMolar * self._sim.processes["PolypeptideElongation"].gtpRequest

		# Set external molecule levels
		coefficient = dryMass / cellMass * self.cellDensity * (self.timeStepSec() * units.s)


		externalMoleculeLevels, newObjective = self.exchangeConstraints(
			self.externalMoleculeIDs,
			coefficient,
			COUNTS_UNITS / VOLUME_UNITS,
			self.environment,
			self.time()
			)

		if newObjective != None and newObjective != self.objective:
			# Build new fba instance with new objective
			self.objective = newObjective
			self.fba_object_options["objective"] = self.objective
			self.fba_object_options = {key:val for key, val in self.fba_object_options.iteritems() if key !="maintenanceCostGAM" and key !="maintenanceReaction"}
			self.fba = FluxBalanceAnalysis(**self.fba_object_options)

			massComposition = self.massReconstruction.getFractionMass(self.doublingTime)
			massInitial = (massComposition["proteinMass"] + massComposition["rnaMass"] + massComposition["dnaMass"]) / self.avgCellToInitialCellConvFactor
			objIds = sorted(self.objective)
			objConc = (COUNTS_UNITS / VOLUME_UNITS) * np.array([self.objective[x] for x in objIds])
			mws = self.getMass(objIds)
			massesToAdd, _ = massesAndCountsToAddForPools(massInitial, objIds, objConc, mws, self.cellDensity, self.nAvogadro)
			smallMoleculePoolsDryMass = units.hstack((massesToAdd[:objIds.index('WATER[c]')], massesToAdd[objIds.index('WATER[c]') + 1:]))
			totalDryMass = units.sum(smallMoleculePoolsDryMass) + massInitial
			self.writeToListener("CellDivision", "expectedDryMassIncrease", totalDryMass)

		# Set external molecule levels
		self.fba.externalMoleculeLevelsIs(externalMoleculeLevels)

		# self.fba.maxReactionFluxIs(self.fba._reactionID_NGAM, (self.ngam * coefficient).asNumber(COUNTS_UNITS / VOLUME_UNITS))
		# self.fba.minReactionFluxIs(self.fba._reactionID_NGAM, (self.ngam * coefficient).asNumber(COUNTS_UNITS / VOLUME_UNITS))

		# self.fba.maxReactionFluxIs(self.fba._reactionID_polypeptideElongationEnergy, polypeptideElongationEnergy.asNumber(COUNTS_UNITS / VOLUME_UNITS))
		# self.fba.minReactionFluxIs(self.fba._reactionID_polypeptideElongationEnergy, polypeptideElongationEnergy.asNumber(COUNTS_UNITS / VOLUME_UNITS))


		#  Find metabolite concentrations from metabolite counts
		metaboliteConcentrations =  countsToMolar * metaboliteCountsInit

		# Make a dictionary of metabolite names to metabolite concentrations
		metaboliteConcentrationsDict = dict(zip(self.metaboliteNames, metaboliteConcentrations.asNumber(COUNTS_UNITS/VOLUME_UNITS)))

		self.fba.internalMoleculeLevelsIs(
			metaboliteConcentrations.asNumber(COUNTS_UNITS / VOLUME_UNITS)
			)

		#  Find enzyme concentrations from enzyme counts
		enzymeCountsInit = self.enzymes.counts()

		enzymeConcentrations = countsToMolar * enzymeCountsInit

		if NONZERO_ENZYMES:
			enzymeConcentrations = countsToMolar * (enzymeCountsInit + 1)

		# Make a dictionary of enzyme names to enzyme concentrations
		enzymeConcentrationsDict = dict(zip(self.enzymeNames, enzymeConcentrations.asNumber(COUNTS_UNITS/VOLUME_UNITS)))

		defaultRate = self.enzymeKinetics.defaultRate

		# Remove any enzyme kinetics paramters for which the needed enzyme and substrate information is not available
		if not self.enzymeKinetics.inputsChecked:
			knownConstraints, unusableConstraints, unknownVals = self.enzymeKinetics.checkKnownSubstratesAndEnzymes(metaboliteConcentrationsDict, enzymeConcentrationsDict, removeUnknowns=True)

		constraintsDict = self.enzymeKinetics.allConstraintsDict(metaboliteConcentrationsDict, enzymeConcentrationsDict)
		reactionsDict = self.enzymeKinetics.allReactionsDict(metaboliteConcentrationsDict, enzymeConcentrationsDict)

		# self.allConstraintsLimits = np.ones(len(self.fba.reactionIDs())) * defaultRate
		# for idx, constraintID in enumerate(self.constraintIDs):
		# 	if constraintID in constraintsDict:
		# 		self.allConstraintsLimits[idx] = constraintsDict[constraintID] * self.timeStepSec()
		# 	else:
		# 		self.allConstraintsLimits[idx] == defaultRate

		self.allConstraintsLimits = np.ones(len(self.fba.reactionIDs())) * defaultRate
		for idx, constraintID in enumerate(self.constraintIDs):
			reactionID = self.constraintToReactionDict[constraintID]
			if reactionID in reactionsDict:
				self.allConstraintsLimits[idx] = np.amax(reactionsDict[reactionID].values()) * self.timeStepSec()
			else:
				self.allConstraintsLimits[idx] == defaultRate

		self.reactionConstraints = np.ones(len(self.fba.reactionIDs())) * np.inf

		currentRateLimits = {}
		# Set reaction fluxes to be between  self.max_flux_coefficient and self.min_flux_coefficient of the predicted rate
		for index, constraintID in enumerate(self.constraintIDs[:(len(self.constraintIDs) // 1)]):
			reactionID = self.constraintToReactionDict[constraintID]

			# Only use this kinetic limit if it's enabled
			if self.constraintMultiplesDict[index]:

				# Skip any reactions not known to the network
				if reactionID in self.reactionNameToReactionIndexDict:
					reactionIndex = self.reactionNameToReactionIndexDict[reactionID]
				else:
					continue

				constraintEstimate = self.allConstraintsLimits[index]
				
				# Determine the min and max fluxes to set for this reaction
				maxFlux = constraintEstimate*self.max_flux_coefficient*self.constraintMultiplesDict[index]
				if self.min_flux_coefficient:
					minFlux = constraintEstimate*self.min_flux_coefficient
				else:
					minFlux = self.min_flux_coefficient
				
				# Make sure to never set negative maximum rates
				assert (constraintEstimate >= 0 and constraintEstimate != np.nan)
				
				# Ensure that this reaction hasn't already been constrained more than this yet
				if reactionID in currentRateLimits and currentRateLimits[reactionID] > maxFlux and currentRateLimits[reactionID] < np.inf:
					# This rate has already been constrained more than this constraint, so skip it
					continue

				# Set the rate limits only if the option flag is enabled
				if USE_RATELIMITS:
					# Set the max reaction rate for this reaction
					self.fba.maxReactionFluxIs(reactionID, maxFlux, raiseForReversible = False)
					# Set the minimum reaction rate for this reaction
					# self.fba.minReactionFluxIs(reactionID, minFlux, raiseForReversible = False)

				# Record what constraint was just applied to this reaction
				currentRateLimits[reactionID] = maxFlux
				self.reactionConstraints[reactionIndex] = maxFlux
			else:

				# Set the reaction max to the default rate (usually infinity)
				self.fba.maxReactionFluxIs(reactionID, defaultRate, raiseForReversible = False)
				# self.fba.minReactionFluxIs(reactionID, 0, raiseForReversible = False)

				# Record that this reaction is at default
				self.reactionConstraints[reactionIndex] = defaultRate

		deltaMetabolites = (1 / countsToMolar) * (COUNTS_UNITS / VOLUME_UNITS * self.fba.outputMoleculeLevelsChange())

		metaboliteCountsFinal = np.fmax(stochasticRound(
			self.randomState,
			metaboliteCountsInit + deltaMetabolites.asNumber()
			), 0).astype(np.int64)

		self.metabolites.countsIs(metaboliteCountsFinal)

		self.overconstraintMultiples = self.fba.reactionFluxes() / self.reactionConstraints

		exFluxes = ((COUNTS_UNITS / VOLUME_UNITS) * self.fba.externalExchangeFluxes() / coefficient).asNumber(units.mmol / units.g / units.h)
		outFluxes = ((COUNTS_UNITS / VOLUME_UNITS) * self.fba.outputMoleculeLevelsChange() / coefficient).asNumber(units.mmol / units.g / units.h)


		# TODO: report as reactions (#) per second & store volume elsewhere
		self.writeToListener("FBAResults", "reactionFluxes",
			self.fba.reactionFluxes() / self.timeStepSec())
		self.writeToListener("FBAResults", "externalExchangeFluxes",
			exFluxes)
		# self.writeToListener("FBAResults", "objectiveValue", # TODO
		# 	self.fba.objectiveValue() / deltaMetabolites.size) # divide to normalize by number of metabolites
		self.writeToListener("FBAResults", "outputFluxes",
			outFluxes)

		self.writeToListener("EnzymeKinetics", "reactionConstraints",
			self.reactionConstraints)

		self.writeToListener("EnzymeKinetics", "allConstraintsLimits",
			self.allConstraintsLimits)

		self.writeToListener("EnzymeKinetics", "overconstraintMultiples",
			self.overconstraintMultiples)

		self.writeToListener("EnzymeKinetics", "metaboliteCountsInit",
			metaboliteCountsInit)

		self.writeToListener("EnzymeKinetics", "metaboliteCountsFinal",
			metaboliteCountsFinal)

		self.writeToListener("EnzymeKinetics", "enzymeCountsInit",
			enzymeCountsInit)

		self.writeToListener("EnzymeKinetics", "metaboliteConcentrations",
			metaboliteConcentrations.asNumber(COUNTS_UNITS / VOLUME_UNITS))

		self.writeToListener("EnzymeKinetics", "countsToMolar",
			countsToMolar.asNumber(COUNTS_UNITS / VOLUME_UNITS))

		self.writeToListener("EnzymeKinetics", "counts_units",
			str(COUNTS_UNITS))

		self.writeToListener("EnzymeKinetics", "mass_units",
			str(MASS_UNITS))

		self.writeToListener("EnzymeKinetics", "volume_units",
			str(VOLUME_UNITS))




		# TODO
		# NOTE: the calculation for the objective components doesn't yet have
		# an interface, since it will vary in calculation and shape for every
		# objective type

		# objectiveComponents_raw = (np.array(self.fba._f).flatten() * self.fba._solutionFluxes)[self.fba._objIndexes]
		# objectiveComponents = objectiveComponents_raw[::2] + objectiveComponents_raw[1::2]

		# self.writeToListener("FBAResults", "objectiveComponents",
		# 	objectiveComponents
		# 	)

		# TODO:
		# - which media exchanges/reactions are limiting, if any
		# - objective details (value, component values)


		self.writeToListener("FBAResults", "outputFluxes",
			outFluxes)
