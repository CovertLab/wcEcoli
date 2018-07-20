#!/usr/bin/env python

"""
Metabolism

Metabolism sub-model. Encodes molecular simulation of microbial metabolism using flux-balance analysis.

TODO:
- option to call a reduced form of metabolism (assume optimal)
- handle oneSidedReaction constraints

@author: Derek Macklin
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 4/2/2013
"""

from __future__ import division

import numpy as np
from scipy.sparse import csr_matrix

import wholecell.processes.process
from wholecell.utils import units

from wholecell.utils.random import stochasticRound
from wholecell.utils.constants import REQUEST_PRIORITY_METABOLISM

from wholecell.utils.modular_fba import FluxBalanceAnalysis

COUNTS_UNITS = units.mmol
VOLUME_UNITS = units.L
MASS_UNITS = units.g
TIME_UNITS = units.s

NONZERO_ENZYMES = False

USE_KINETICS = True
KINETICS_BURN_IN_PERIOD = 0

class Metabolism(wholecell.processes.process.Process):
	""" Metabolism """

	_name = "Metabolism"

	# Constructor
	def __init__(self):

		super(Metabolism, self).__init__()

	# Construct object graph
	def initialize(self, sim, sim_data):
		super(Metabolism, self).initialize(sim, sim_data)

		# initialize exchange_data based on initial nutrient condition
		self.exchange_data = self._initExchangeData(sim_data)

		#TODO (Eran) this can be remove once transport is in place
		self.exchange_data_dict = sim_data.exchange_data_dict.copy()

		## lists for updateImportConstraint
		# TODO (Eran) remove these once transport kinetics are working

		# import exchange molecules can be both constrained/unconstrained
		self.importExchangeMolecules_noGLC = self.exchange_data['importExchangeMolecules'][:]
		self.importExchangeMolecules_noGLC.remove('GLC[p]')

		# dictionary of with conditions specifying sets of molecules that determine
		# glc's upper bound for FBA import constraint.
		self.glc_vmax_conditions = {
			'glc_vmax_condition_1': ['GLC[p]'],
			'glc_vmax_condition_2': ['CA+2[p]', 'MG+2[p]'],
			'glc_vmax_condition_3': ['CPD-183[p]', 'INDOLE[p]', 'NITRATE[p]', 'NITRITE[p]', 'CPD-520[p]', 'TUNGSTATE[p]'],
			'glc_vmax_condition_4': ['OXYGEN-MOLECULE[p]'],
			}

		# Load constants
		self.nAvogadro = sim_data.constants.nAvogadro
		self.cellDensity = sim_data.constants.cellDensity
		self.ngam = sim_data.constants.nonGrowthAssociatedMaintenance

		self.exchangeConstraints = sim_data.process.metabolism.exchangeConstraints

		self.getBiomassAsConcentrations = sim_data.mass.getBiomassAsConcentrations
		self.nutrientToDoublingTime = sim_data.nutrientToDoublingTime

		# Create objective for homeostatic constraints
		nutrients_time_series_label = sim_data.external_state.environment.nutrients_time_series_label

		concDict = sim_data.process.metabolism.concentrationUpdates.concentrationsBasedOnNutrients(
			self.exchange_data
			)
		self.concModificationsBasedOnCondition = self.getBiomassAsConcentrations(
			sim_data.conditionToDoublingTime[sim_data.condition]
			)
		concDict.update(self.concModificationsBasedOnCondition)
		self.homeostaticObjective = dict((key, concDict[key].asNumber(COUNTS_UNITS / VOLUME_UNITS)) for key in concDict)

		# Load initial mass
		initWaterMass = sim_data.mass.avgCellWaterMassInit
		initDryMass = sim_data.mass.avgCellDryMassInit
		initCellMass = initWaterMass + initDryMass

		energyCostPerWetMass = sim_data.constants.darkATP * initDryMass / initCellMass

		# Identify all molecules in external environment that can be exchanged for the given time series
		# TODO (eran) initialize externalExchangeMolecules without exchange_data_dict
		externalExchangedMolecules = sim_data.exchange_data_dict["secretionExchangeMolecules"]
		self.metaboliteNamesFromNutrients = set()

		for time, nutrient_label, volume in sim_data.external_state.environment.nutrients_time_series[
				nutrients_time_series_label]:
			# get exchange data for each nutrient condition in time series
			nutrient_label_exchange_data = sim_data.process.metabolism._getExchangeData(nutrient_label)
			externalExchangedMolecules += nutrient_label_exchange_data["importExchangeMolecules"]
			self.metaboliteNamesFromNutrients.update(
				sim_data.process.metabolism.concentrationUpdates.concentrationsBasedOnNutrients(
					nutrient_label_exchange_data, sim_data.process.metabolism.nutrientsToInternalConc
					)
				)
		externalExchangedMolecules = sorted(set(externalExchangedMolecules))

		# save nutrient names for environment view
		self.environment_nutrients_names = externalExchangedMolecules

		self.metaboliteNamesFromNutrients = sorted(self.metaboliteNamesFromNutrients)

		moleculeMasses = dict(zip(externalExchangedMolecules, sim_data.getter.getMass(externalExchangedMolecules).asNumber(MASS_UNITS / COUNTS_UNITS)))

		# Data structures to compute reaction bounds based on enzyme presence/absence
		self.catalystsList = sim_data.process.metabolism.catalystsList
		self.reactionsWithCatalystsList = sim_data.process.metabolism.reactionCatalystsList
		self.reactionCatalystsDict = sim_data.process.metabolism.reactionCatalysts

		catalysisMatrixI = sim_data.process.metabolism.catalysisMatrixI
		catalysisMatrixJ = sim_data.process.metabolism.catalysisMatrixJ
		catalysisMatrixV = sim_data.process.metabolism.catalysisMatrixV

		shape = (catalysisMatrixI.max() + 1, catalysisMatrixJ.max() + 1)
		self.catalysisMatrix = csr_matrix((catalysisMatrixV, (catalysisMatrixI, catalysisMatrixJ)), shape = shape)

		self.catalyzedReactionBoundsPrev = np.inf * np.ones(len(self.reactionsWithCatalystsList))

		# Function to compute reaction targets based on kinetic parameters and molecule concentrations
		self.getKineticConstraints = sim_data.process.metabolism.getKineticConstraints

		self.useAllConstraints = sim_data.process.metabolism.useAllConstraints
		self.constraintsToDisable = sim_data.process.metabolism.constraintsToDisable
		self.kineticsConstrainedReactions = sim_data.process.metabolism.constrainedReactionList
		if hasattr(sim_data.process.metabolism, "kineticTargetShuffleRxns") and sim_data.process.metabolism.kineticTargetShuffleRxns != None:
			self.kineticsConstrainedReactions = sim_data.process.metabolism.kineticTargetShuffleRxns
			self.useAllConstraints = True
		self.kineticsEnzymesList = sim_data.process.metabolism.enzymeIdList
		self.kineticsSubstratesList = sim_data.process.metabolism.kineticsSubstratesList

		constraintToReactionMatrixI = sim_data.process.metabolism.constraintToReactionMatrixI
		constraintToReactionMatrixJ = sim_data.process.metabolism.constraintToReactionMatrixJ
		constraintToReactionMatrixV = sim_data.process.metabolism.constraintToReactionMatrixV
		shape = (constraintToReactionMatrixI.max() + 1, constraintToReactionMatrixJ.max() + 1)
		self.constraintToReactionMatrix = np.zeros(shape, np.float64)
		self.constraintToReactionMatrix[constraintToReactionMatrixI, constraintToReactionMatrixJ] = constraintToReactionMatrixV
		self.constraintIsKcatOnly = sim_data.process.metabolism.constraintIsKcatOnly

		# Select solver and associated kinetic objective weight (lambda)
		solver = "glpk-linear"
		if "linear" in solver:
			kineticObjectiveWeight = sim_data.constants.metabolismKineticObjectiveWeightLinear
		else:
			kineticObjectiveWeight = sim_data.constants.metabolismKineticObjectiveWeightQuadratic

		# Set up FBA solver
		# reactionRateTargets value is just for initialization, it gets reset each timestep during evolveState
		self.fbaObjectOptions = {
			"reactionStoich" : sim_data.process.metabolism.reactionStoich,
			"externalExchangedMolecules" : externalExchangedMolecules,
			"objective" : self.homeostaticObjective,
			"objectiveType" : "homeostatic_kinetics_mixed",
			"objectiveParameters" : {
					"kineticObjectiveWeight" : kineticObjectiveWeight,
					"reactionRateTargets" : {reaction : 1 for reaction in self.kineticsConstrainedReactions},
					"oneSidedReactionTargets" : [],
					},
			"moleculeMasses" : moleculeMasses,
			"secretionPenaltyCoeff" : sim_data.constants.secretion_penalty_coeff, # The "inconvenient constant"--limit secretion (e.g., of CO2)
			"solver" : solver,
			"maintenanceCostGAM" : energyCostPerWetMass.asNumber(COUNTS_UNITS / MASS_UNITS),
			"maintenanceReaction" : sim_data.process.metabolism.maintenanceReaction,
		}
		if USE_KINETICS == False:
			self.fbaObjectOptions["objectiveType"] = "homeostatic"
		self.fba = FluxBalanceAnalysis(**self.fbaObjectOptions)

		self.internalExchangeIdxs = np.array([self.metaboliteNamesFromNutrients.index(x) for x in self.fba.getOutputMoleculeIDs()])

		# Disable all rates during burn-in
		if USE_KINETICS:
			if KINETICS_BURN_IN_PERIOD > 0:
				self.fba.disableKineticTargets()
				self.burnInComplete = False
			else:
				self.burnInComplete = True
				if not self.useAllConstraints:
					for rxn in self.constraintsToDisable:
						self.fba.disableKineticTargets(rxn)

		# Values will get updated at each time point
		self.currentNgam = 1 * (COUNTS_UNITS / VOLUME_UNITS)
		self.currentPolypeptideElongationEnergy = 1 * (COUNTS_UNITS / VOLUME_UNITS)

		# External molecules
		self.externalMoleculeIDs = self.fba.getExternalMoleculeIDs()

		## Views
		# views of environment
		self.environment_nutrients = self.environmentView(self.environment_nutrients_names)

		# views of environment, specific for determining FBA import constraints
		self.glc_vmax_condition_1 = self.environmentView(self.glc_vmax_conditions['glc_vmax_condition_1'])
		self.glc_vmax_condition_2 = self.environmentView(self.glc_vmax_conditions['glc_vmax_condition_2'])
		self.glc_vmax_condition_3 = self.environmentView(self.glc_vmax_conditions['glc_vmax_condition_3'])
		self.glc_vmax_condition_4 = self.environmentView(self.glc_vmax_conditions['glc_vmax_condition_4'])

		# views for metabolism
		self.metaboliteNames = self.fba.getOutputMoleculeIDs()
		self.metabolites = self.bulkMoleculesView(self.metaboliteNamesFromNutrients)
		self.catalysts = self.bulkMoleculesView(self.catalystsList)
		self.kineticsEnzymes = self.bulkMoleculesView(self.kineticsEnzymesList)
		self.kineticsSubstrates = self.bulkMoleculesView(self.kineticsSubstratesList)

		outputMoleculeIDs = self.fba.getOutputMoleculeIDs()

		assert outputMoleculeIDs == self.fba.getInternalMoleculeIDs()

		# Set the priority to a low value
		self.bulkMoleculesRequestPriorityIs(REQUEST_PRIORITY_METABOLISM)

		self.AAs = [x[:-3] for x in sorted(sim_data.amino_acid_1_to_3_ordered.values())]

		self.shuffleIdxs = None
		if hasattr(sim_data.process.metabolism, "kineticTargetShuffleIdxs") and sim_data.process.metabolism.kineticTargetShuffleIdxs != None:
			self.shuffleIdxs = sim_data.process.metabolism.kineticTargetShuffleIdxs

		self.shuffleCatalyzedIdxs = None
		if hasattr(sim_data.process.metabolism, "catalystShuffleIdxs") and sim_data.process.metabolism.catalystShuffleIdxs != None:
			self.shuffleCatalyzedIdxs = sim_data.process.metabolism.catalystShuffleIdxs

	def calculateRequest(self):
		self.metabolites.requestAll()
		self.catalysts.requestAll()
		self.kineticsEnzymes.requestAll()
		self.kineticsSubstrates.requestAll()

	def evolveState(self):
		metaboliteCountsInit = self.metabolites.counts()

		cellMass = (self.readFromListener("Mass", "cellMass") * units.fg)
		dryMass = (self.readFromListener("Mass", "dryMass") * units.fg)

		cellVolume = cellMass / self.cellDensity
		countsToMolar = 1 / (self.nAvogadro * cellVolume)

		current_nutrients = self._external_states['Environment'].nutrients
		self.concModificationsBasedOnCondition = self.getBiomassAsConcentrations(
			self.nutrientToDoublingTime.get(current_nutrients, self.nutrientToDoublingTime["minimal"])
			)

		# Coefficient to convert between flux (mol/g DCW/hr) basis and concentration (M) basis
		coefficient = dryMass / cellMass * self.cellDensity * (self.timeStepSec() * units.s)


		# TODO (Eran) set nutrients with 0 concentration to 0 vmax

		# Update FBA import constraint variables based on current nutrient concentrations
		self._updateImportConstraint()

		# Set external molecule levels
		externalMoleculeLevels, newObjective = self.exchangeConstraints(
			self.externalMoleculeIDs,
			coefficient,
			COUNTS_UNITS / VOLUME_UNITS,
			self.exchange_data,
			self.concModificationsBasedOnCondition,
			)

		updatedObjective = False
		if newObjective != None and newObjective != self.homeostaticObjective:
			# Build new fba instance with new objective
			self.fbaObjectOptions["objective"] = newObjective
			self.fba = FluxBalanceAnalysis(**self.fbaObjectOptions)
			self.internalExchangeIdxs = np.array([self.metaboliteNamesFromNutrients.index(x) for x in self.fba.getOutputMoleculeIDs()])
			self.homeostaticObjective = newObjective
			updatedObjective = True

		# After completing the burn-in, enable kinetic rates
		if (USE_KINETICS) and (not self.burnInComplete) and (self._sim.time() > KINETICS_BURN_IN_PERIOD):
			self.burnInComplete = True
			self.fba.enableKineticTargets()

		# Allow flexibility for solver in first time step after an environment shift
		if updatedObjective:
			self.fba.disableKineticTargets()
			self.burnInComplete = False

		#  Find metabolite concentrations from metabolite counts
		metaboliteConcentrations =  countsToMolar * metaboliteCountsInit[self.internalExchangeIdxs]

		# Make a dictionary of metabolite names to metabolite concentrations
		metaboliteConcentrationsDict = dict(zip(self.metaboliteNames, metaboliteConcentrations))
		self.fba.setInternalMoleculeLevels(metaboliteConcentrations.asNumber(COUNTS_UNITS / VOLUME_UNITS))

		# Set external molecule levels
		self._setExternalMoleculeLevels(externalMoleculeLevels, metaboliteConcentrations)

		# Change the ngam and polypeptide elongation energy penalty only if they are noticably different from the current value
		ADJUSTMENT_RATIO = .01

		# Calculate new NGAM and update if necessary
		self.newNgam = self.ngam * coefficient
		ngam_diff = np.abs(self.currentNgam.asNumber() - self.newNgam.asNumber()) / (self.currentNgam.asNumber() + 1e-20)
		if ngam_diff > ADJUSTMENT_RATIO:
			self.currentNgam = self.newNgam
			flux = (self.ngam * coefficient).asNumber(COUNTS_UNITS / VOLUME_UNITS)
			self.fba.setReactionFluxBounds(self.fba._reactionID_NGAM, lowerBounds=flux, upperBounds=flux)

		# Calculate GTP usage based on how much was needed in polypeptide elongation in previous step and update if necessary
		newPolypeptideElongationEnergy = countsToMolar * 0
		if hasattr(self._sim.processes["PolypeptideElongation"], "gtpRequest"):
			newPolypeptideElongationEnergy = countsToMolar * self._sim.processes["PolypeptideElongation"].gtpRequest
		poly_diff = np.abs((self.currentPolypeptideElongationEnergy.asNumber() - newPolypeptideElongationEnergy.asNumber())) / (self.currentPolypeptideElongationEnergy.asNumber() + 1e-20)
		if poly_diff > ADJUSTMENT_RATIO:
			self.currentPolypeptideElongationEnergy = newPolypeptideElongationEnergy
			flux = self.currentPolypeptideElongationEnergy.asNumber(COUNTS_UNITS / VOLUME_UNITS)
			self.fba.setReactionFluxBounds(self.fba._reactionID_polypeptideElongationEnergy, lowerBounds=flux, upperBounds=flux)

		# Read counts for catalysts and enzymes (catalysts with kinetics constraints)
		catalystsCountsInit = self.catalysts.counts()

		kineticsEnzymesCountsInit = self.kineticsEnzymes.counts()
		kineticsEnzymesConcentrations = countsToMolar * kineticsEnzymesCountsInit

		kineticsSubstratesCountsInit = self.kineticsSubstrates.counts()
		kineticsSubstratesConcentrations = countsToMolar * kineticsSubstratesCountsInit

		# Add one of every enzyme to ensure none are zero
		if NONZERO_ENZYMES:
			catalystsCountsInit += 1
			kineticsEnzymesConcentrations = countsToMolar * (kineticsEnzymesCountsInit + 1)

		# Calculate and set kinetic targets if kinetics is enabled
		if USE_KINETICS and self.burnInComplete:
			# Set hard upper bounds constraints based on enzyme presence (infinite upper bound) or absence (upper bound of zero)
			catalyzedReactionBounds = np.inf * np.ones(len(self.reactionsWithCatalystsList))
			rxnPresence = self.catalysisMatrix.dot(catalystsCountsInit)
			catalyzedReactionBounds[rxnPresence == 0] = 0
			if self.shuffleCatalyzedIdxs is not None:
				catalyzedReactionBounds = catalyzedReactionBounds[self.shuffleCatalyzedIdxs]


			# Only update reaction limits that are different from previous time step
			updateIdxs = np.where(catalyzedReactionBounds != self.catalyzedReactionBoundsPrev)[0]
			updateRxns = [self.reactionsWithCatalystsList[idx] for idx in updateIdxs]
			updateVals = catalyzedReactionBounds[updateIdxs]
			self.fba.setReactionFluxBounds(updateRxns, upperBounds=updateVals, raiseForReversible=False)
			self.catalyzedReactionBoundsPrev = catalyzedReactionBounds

			# Set target fluxes for reactions based on their most relaxed constraint
			constraintValues = self.getKineticConstraints(
				kineticsEnzymesConcentrations.asNumber(units.umol / units.L),
				kineticsSubstratesConcentrations.asNumber(units.umol / units.L),
				)
			reactionTargets = (units.umol / units.L / units.s) * np.max(self.constraintToReactionMatrix * constraintValues, axis = 1)

			# Shuffle parameters (only performed in very specific cases)
			if self.shuffleIdxs is not None:
				reactionTargets = (units.umol / units.L / units.s) * reactionTargets.asNumber()[self.shuffleIdxs]

			# record which constraint was used, add constraintToReactionMatrix to ensure the index is one of the constraints if multiplication is 0
			reactionConstraint = np.argmax(self.constraintToReactionMatrix * constraintValues + self.constraintToReactionMatrix, axis = 1)

			# set reaction flux targets
			targets = (TIME_UNITS * self.timeStepSec() * reactionTargets).asNumber(COUNTS_UNITS / VOLUME_UNITS)
			self.fba.setKineticTarget(self.kineticsConstrainedReactions, targets, raiseForReversible = False)

			if not self.useAllConstraints:
				self.fba.disableKineticTargets(self.constraintsToDisable)

		# Solve FBA problem and update metabolite counts
		deltaMetabolites = (1 / countsToMolar) * (COUNTS_UNITS / VOLUME_UNITS * self.fba.getOutputMoleculeLevelsChange())

		metaboliteCountsFinal = np.zeros_like(metaboliteCountsInit)
		metaboliteCountsFinal[self.internalExchangeIdxs] = np.fmax(stochasticRound(
			self.randomState,
			metaboliteCountsInit[self.internalExchangeIdxs] + deltaMetabolites.asNumber()
			), 0).astype(np.int64)

		self.metabolites.countsIs(metaboliteCountsFinal)

		exFluxes = ((COUNTS_UNITS / VOLUME_UNITS) * self.fba.getExternalExchangeFluxes() / coefficient).asNumber(units.mmol / units.g / units.h)

		# change in nutrient counts, used in non-infinite environments
		delta_nutrients = ((1 / countsToMolar) * (COUNTS_UNITS / VOLUME_UNITS) * self.fba.getExternalExchangeFluxes()).asNumber().astype(int)
		self.environment_nutrients.changeCounts(delta_nutrients)


		# Write outputs to listeners
		self.writeToListener("FBAResults", "deltaMetabolites", metaboliteCountsFinal - metaboliteCountsInit)
		self.writeToListener("FBAResults", "reactionFluxes", self.fba.getReactionFluxes() / self.timeStepSec())
		self.writeToListener("FBAResults", "externalExchangeFluxes", exFluxes)
		self.writeToListener("FBAResults", "objectiveValue", self.fba.getObjectiveValue())
		self.writeToListener("FBAResults", "shadowPrices", self.fba.getShadowPrices(self.metaboliteNames))
		self.writeToListener("FBAResults", "reducedCosts", self.fba.getReducedCosts(self.fba.getReactionIDs()))
		self.writeToListener("FBAResults", "targetConcentrations", [self.homeostaticObjective[mol] for mol in self.fba.getHomeostaticTargetMolecules()])
		self.writeToListener("FBAResults", "homeostaticObjectiveValues", self.fba.getHomeostaticObjectiveValues())

		self.writeToListener("EnzymeKinetics", "metaboliteCountsInit", metaboliteCountsInit)
		self.writeToListener("EnzymeKinetics", "metaboliteCountsFinal", metaboliteCountsFinal)
		self.writeToListener("EnzymeKinetics", "enzymeCountsInit", kineticsEnzymesCountsInit)
		self.writeToListener("EnzymeKinetics", "metaboliteConcentrations", metaboliteConcentrations.asNumber(COUNTS_UNITS / VOLUME_UNITS))
		self.writeToListener("EnzymeKinetics", "countsToMolar", countsToMolar.asNumber(COUNTS_UNITS / VOLUME_UNITS))
		self.writeToListener("EnzymeKinetics", "actualFluxes", self.fba.getReactionFluxes(self.kineticsConstrainedReactions) / self.timeStepSec())

		if USE_KINETICS and self.burnInComplete:
			self.writeToListener("EnzymeKinetics", "targetFluxes", targets / self.timeStepSec())
			self.writeToListener("EnzymeKinetics", "reactionConstraint", reactionConstraint)

	# limit amino acid uptake to what is needed to meet concentration objective to prevent use as carbon source
	def _setExternalMoleculeLevels(self, externalMoleculeLevels, metaboliteConcentrations):
		for aa in self.AAs:
			if aa + "[p]" in self.fba.getExternalMoleculeIDs():
				idx = self.externalMoleculeIDs.index(aa + "[p]")
			elif aa + "[c]" in self.fba.getExternalMoleculeIDs():
				idx = self.externalMoleculeIDs.index(aa + "[c]")
			else:
				continue

			concDiff = self.homeostaticObjective[aa + "[c]"] - metaboliteConcentrations[self.metaboliteNames.index(aa + "[c]")].asNumber(COUNTS_UNITS / VOLUME_UNITS)
			if concDiff < 0:
				concDiff = 0

			if externalMoleculeLevels[idx] > concDiff:
				externalMoleculeLevels[idx] =  concDiff

		self.fba.setExternalMoleculeLevels(externalMoleculeLevels)


	def _initExchangeData(self, sim_data):
		'''
		Returns a dictionary with the five categories of exchange data used by FBA.

		The categories of molecules include:
			- externalExchangeMolecules: All exchange molecules, both import and
				secretion exchanged molecules.
			- importExchangeMolecules: molecules that can be imported from the
				environment into the cell.
			- importConstrainedExchangeMolecules: exchange molecules that have
				an upper bound on their flux.
			- importUnconstrainedExchangeMolecules: exchange molecules that do
				not have an upper bound on their flux.
			- secretionExchangeMolecules: molecules that can be secreted by the
				cell into the environment.

		'''

		exchange_data_dict = sim_data.exchange_data_dict.copy()
		nutrient_label = sim_data.external_state.environment.nutrients_time_series[
			sim_data.external_state.environment.nutrients_time_series_label
			][0][1]

		# all nutrients from nutrients_def
		externalExchangeMolecules = exchange_data_dict['externalExchangeMolecules'][nutrient_label]
		importExchangeMolecules = exchange_data_dict['importExchangeMolecules'][nutrient_label]
		importConstrainedExchangeMolecules = exchange_data_dict['importConstrainedExchangeMolecules'][nutrient_label]
		importUnconstrainedExchangeMolecules = exchange_data_dict['importUnconstrainedExchangeMolecules'][nutrient_label]
		secretionExchangeMolecules = exchange_data_dict['secretionExchangeMolecules']

		return {
			"externalExchangeMolecules": externalExchangeMolecules,
			"importExchangeMolecules": importExchangeMolecules,
			"importConstrainedExchangeMolecules": importConstrainedExchangeMolecules,
			"importUnconstrainedExchangeMolecules": importUnconstrainedExchangeMolecules,
			"secretionExchangeMolecules": secretionExchangeMolecules,
			}

	def _updateImportConstraint(self):
		'''
		Update importExchangeMolecules for FBA based on current nutrient concentrations.

		This provides a simple type of transport to accomodate changing nutrient
		concentrations in the environment. Transport is modeled as a binary switch:
		When there is a high concentrations of environment nutrients, transporters
		are unconstrained and transport nutrients as needed. When concentrations
		fall below a threshold, k_m, the transport is constrained to 0 and don't
		let any nutrients through.



		TODO (Eran) Glucose is treated differently.

		Notes
		-----
		- TODO (ERAN) eliminate this when kinetic transport process is operational
		- TODO (Eran) with importConstrained changing, importUnconstrained also needs to change.
		- TODO (Eran) importConstrained + importUnconstrained = importExchange

		'''

		k_m = 1 # (units.mmol / units.L). John suggests 10 micromolar

		# currently constrained molecules
		constrained_ids = self.exchange_data['importConstrainedExchangeMolecules'].keys()

		## identify nutrients that crossed the concentration threshold (k_m)
		below_thresh_ids = []
		above_thresh_ids = []

		# go through all environmental nutrient concentrations
		for idx, conc in enumerate(self.environment_nutrients.totalConcentrations()):
			nutrient_name = self.environment_nutrients_names[idx]

			# Separate nutrients that are above and below threshold, k_m
			# Only use nutrients in importExchangeMolecules_noGLC (GLC always be constrained)
			if nutrient_name in self.importExchangeMolecules_noGLC:
				if (conc <= k_m and not np.isnan(conc)):
					below_thresh_ids.append(nutrient_name)
				elif (conc >= k_m and not np.isnan(conc)):
					above_thresh_ids.append(nutrient_name)

		new_below_thresh_ids = np.setdiff1d(below_thresh_ids, constrained_ids)
		new_above_thresh_ids = np.intersect1d(above_thresh_ids,	constrained_ids)

		# if newly below threshold, add molecule to import constraint and set max flux to 0
		for id in new_below_thresh_ids:
			self.exchange_data['importConstrainedExchangeMolecules'][id] = 0 * (units.mmol / units.g / units.h)
			self.exchange_data['importUnconstrainedExchangeMolecules'].remove(id)

		# remove molecule from import constraint if newly above threshold
		for id in new_above_thresh_ids:
			self.exchange_data['importUnconstrainedExchangeMolecules'].append(id)
			self.exchange_data['importConstrainedExchangeMolecules'].pop(id, None)

		## Glucose always has an import constraint, with a flux upper bound depending on environment
		# if any molecules in condition_1 have a concentration of 0, set glc flux upper bound  to 0
		if any(self.glc_vmax_condition_1.totalConcentrations() == 0):
			self.exchange_data['importConstrainedExchangeMolecules']['GLC[p]']._value = 0

		# if any molecules in condition_2 have a concentration of 0, set glc flux upper bound  to 10
		elif any(self.glc_vmax_condition_2.totalConcentrations() == 0):
			self.exchange_data['importConstrainedExchangeMolecules']['GLC[p]']._value = 10

		# if any molecules in condition_3 do not have a concentration of 0, set glc flux upper bound  to 10
		elif any(self.glc_vmax_condition_3.totalConcentrations() != 0):
			self.exchange_data['importConstrainedExchangeMolecules']['GLC[p]']._value = 10

		# if any molecules in condition_4 have a concentration of 0, set glc flux upper bound  to 100
		elif any(self.glc_vmax_condition_4.totalConcentrations() == 0):
			self.exchange_data['importConstrainedExchangeMolecules']['GLC[p]']._value = 100

		# if normal condition, set glc vmax to 20
		else:
			self.exchange_data['importConstrainedExchangeMolecules']['GLC[p]']._value = 20
