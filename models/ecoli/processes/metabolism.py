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
CONC_UNITS = COUNTS_UNITS / VOLUME_UNITS

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

		# Local sim_data references
		metabolism = sim_data.process.metabolism
		constants = sim_data.constants
		mass = sim_data.mass

		# Load constants
		self.nAvogadro = constants.nAvogadro
		self.cellDensity = constants.cellDensity
		self.ngam = constants.nonGrowthAssociatedMaintenance
		energyCostPerWetMass = constants.darkATP * mass.cellDryMassFraction

		self.exchangeConstraints = metabolism.exchangeConstraints

		self._biomass_concentrations = {}
		self._getBiomassAsConcentrations = mass.getBiomassAsConcentrations
		self.nutrientToDoublingTime = sim_data.nutrientToDoublingTime

		self.use_trna_charging = sim._trna_charging

		# Include ppGpp concentration target in objective if not handled kinetically in other processes
		self.include_ppgpp = not sim._ppgpp_regulation or not self.use_trna_charging
		self.ppgpp_id = sim_data.moleculeIds.ppGpp
		self.getppGppConc = sim_data.growthRateParameters.getppGppConc

		# Use information from the environment
		environment = self._external_states['Environment']

		# go through all media in the timeline and add to metaboliteNames
		self.metaboliteNamesFromNutrients = set()
		exchange_molecules = set()
		if self.include_ppgpp:
			self.metaboliteNamesFromNutrients.add(self.ppgpp_id)
		for time, media_id in environment.current_timeline:
			self.metaboliteNamesFromNutrients.update(
				metabolism.concentrationUpdates.concentrationsBasedOnNutrients(media_id)
				)
			exchanges = sim_data.external_state.exchange_data_from_media(media_id)
			exchange_molecules.update(exchanges['externalExchangeMolecules'])
		self.metaboliteNamesFromNutrients = list(sorted(self.metaboliteNamesFromNutrients))
		exchange_molecules = list(sorted(exchange_molecules))
		moleculeMasses = dict(zip(exchange_molecules,
			sim_data.getter.getMass(exchange_molecules).asNumber(MASS_UNITS / COUNTS_UNITS)))

		concDict = metabolism.concentrationUpdates.concentrationsBasedOnNutrients(
			environment.current_media_id
			)
		doubling_time = sim_data.conditionToDoublingTime[sim_data.condition]
		concModificationsBasedOnCondition = self.getBiomassAsConcentrations(
			doubling_time
			)
		concDict.update(concModificationsBasedOnCondition)

		if self.include_ppgpp:
			concDict[self.ppgpp_id] = self.getppGppConc(doubling_time)

		self.homeostaticObjective = dict((key, concDict[key].asNumber(CONC_UNITS)) for key in concDict)


		# Data structures to compute reaction bounds based on enzyme presence/absence
		self.catalyst_ids = metabolism.catalyst_ids
		self.reactions_with_catalyst = metabolism.reactions_with_catalyst

		catalysisMatrixI = metabolism.catalysisMatrixI
		catalysisMatrixJ = metabolism.catalysisMatrixJ
		catalysisMatrixV = metabolism.catalysisMatrixV

		shape = (catalysisMatrixI.max() + 1, catalysisMatrixJ.max() + 1)
		self.catalysisMatrix = csr_matrix((catalysisMatrixV, (catalysisMatrixI, catalysisMatrixJ)), shape = shape)

		# Function to compute reaction targets based on kinetic parameters and molecule concentrations
		self.getKineticConstraints = metabolism.getKineticConstraints

		# Remove disabled reactions so they don't get included in the FBA problem setup
		if hasattr(metabolism, "kineticTargetShuffleRxns") and metabolism.kineticTargetShuffleRxns is not None:
			self.kinetics_constrained_reactions = metabolism.kineticTargetShuffleRxns
			self.active_constraints_mask = np.ones(len(self.kinetics_constrained_reactions), dtype=bool)
		else:
			kinetic_constraint_reactions = metabolism.kinetic_constraint_reactions
			constraintsToDisable = metabolism.constraintsToDisable
			self.active_constraints_mask = np.array([(rxn not in constraintsToDisable) for rxn in kinetic_constraint_reactions])
			self.kinetics_constrained_reactions = list(np.array(kinetic_constraint_reactions)[self.active_constraints_mask])

		# Add kinetic reaction targets from boundary
		self.boundary_constrained_reactions = environment.transport_fluxes.keys()
		self.all_constrained_reactions = self.kinetics_constrained_reactions + self.boundary_constrained_reactions

		self.kinetic_constraint_enzymes = metabolism.kinetic_constraint_enzymes
		self.kinetic_constraint_substrates = metabolism.kinetic_constraint_substrates

		# Set solver and kinetic objective weight (lambda)
		solver = metabolism.solver
		kinetic_objective_weight = metabolism.kinetic_objective_weight
		kinetic_objective_weight_in_range = metabolism.kinetic_objective_weight_in_range

		# Disable kinetics completely if weight is 0 or specified in file above
		if not USE_KINETICS or kinetic_objective_weight == 0:
			objective_type = 'homeostatic'
			self.use_kinetics = False
			kinetic_objective_weight = 0
		else:
			objective_type = 'homeostatic_kinetics_mixed'
			self.use_kinetics = True

		# Set up FBA solver
		# reactionRateTargets value is just for initialization, it gets reset each timestep during evolveState
		self.fbaObjectOptions = {
			"reactionStoich": metabolism.reactionStoich,
			"externalExchangedMolecules": exchange_molecules,
			"objective": self.homeostaticObjective,
			"objectiveType": objective_type,
			"objectiveParameters": {
					"kineticObjectiveWeight": kinetic_objective_weight,
					'kinetic_objective_weight_in_range': kinetic_objective_weight_in_range,
					"reactionRateTargets": {reaction: 1 for reaction in self.all_constrained_reactions},
					"oneSidedReactionTargets": [],
					},
			"moleculeMasses": moleculeMasses,
			"secretionPenaltyCoeff": metabolism.secretion_penalty_coeff, # The "inconvenient constant"--limit secretion (e.g., of CO2)
			"solver": solver,
			"maintenanceCostGAM": energyCostPerWetMass.asNumber(COUNTS_UNITS / MASS_UNITS),
			"maintenanceReaction": metabolism.maintenanceReaction,
		}
		self.fba = FluxBalanceAnalysis(**self.fbaObjectOptions)

		self.internalExchangeIdxs = np.array([self.metaboliteNamesFromNutrients.index(x) for x in self.fba.getOutputMoleculeIDs()])

		# Disable all rates during burn-in
		if self.use_kinetics:
			if KINETICS_BURN_IN_PERIOD > 0:
				self.fba.disableKineticTargets()
				self.burnInComplete = False
			else:
				self.burnInComplete = True

		## Construct views
		# views on metabolism bulk molecules
		self.metaboliteNames = self.fba.getOutputMoleculeIDs()
		self.metabolites = self.bulkMoleculesView(self.metaboliteNamesFromNutrients)
		self.catalysts = self.bulkMoleculesView(self.catalyst_ids)
		self.kineticsEnzymes = self.bulkMoleculesView(self.kinetic_constraint_enzymes)
		self.kineticsSubstrates = self.bulkMoleculesView(self.kinetic_constraint_substrates)

		assert self.metaboliteNames == self.fba.getInternalMoleculeIDs()

		# Set the priority to a low value
		self.bulkMoleculesRequestPriorityIs(REQUEST_PRIORITY_METABOLISM)

		self.aa_names_no_location = [x[:-3] for x in sorted(sim_data.amino_acid_1_to_3_ordered.values())]

		self.shuffleIdxs = None
		if hasattr(metabolism, "kineticTargetShuffleIdxs") and metabolism.kineticTargetShuffleIdxs is not None:
			self.shuffleIdxs = metabolism.kineticTargetShuffleIdxs

		self.shuffleCatalyzedIdxs = None
		if hasattr(metabolism, "catalystShuffleIdxs") and metabolism.catalystShuffleIdxs is not None:
			self.shuffleCatalyzedIdxs = metabolism.catalystShuffleIdxs

		# Track updated AA concentration targets with tRNA charging
		self.aa_targets = {}
		self.aa_targets_not_updated = set(['L-SELENOCYSTEINE[c]'])
		self.aa_names = sim_data.moleculeGroups.aaIDs
		self.aas = self.bulkMoleculesView(self.aa_names)

	def calculateRequest(self):
		self.metabolites.requestAll()
		self.catalysts.requestAll()
		self.kineticsEnzymes.requestAll()
		self.kineticsSubstrates.requestAll()

	def evolveState(self):
		# After completing the burn-in, enable kinetic rates
		if self.use_kinetics and (not self.burnInComplete) and (self._sim.time() > KINETICS_BURN_IN_PERIOD):
			self.burnInComplete = True
			self.fba.enableKineticTargets()

		time_step = self.timeStepSec() * units.s
		metaboliteCountsInit = self.metabolites.counts()
		cellMass = (self.readFromListener("Mass", "cellMass") * units.fg)
		dryMass = (self.readFromListener("Mass", "dryMass") * units.fg)

		cellVolume = cellMass / self.cellDensity
		countsToMolar = 1 / (self.nAvogadro * cellVolume)

		# Coefficient to convert between flux (mol/g DCW/hr) basis and concentration (M) basis
		coefficient = dryMass / cellMass * self.cellDensity * time_step

		# Get environment updates
		environment = self._external_states['Environment']
		current_media_id = environment.current_media_id
		exchange_data = environment.get_exchange_data()

		# make sure there are no new flux targets from the boundary
		assert set(environment.transport_fluxes.keys()).issubset(self.all_constrained_reactions)

		# Determine updates to concentrations depending on the current state
		doubling_time = self.nutrientToDoublingTime.get(current_media_id, self.nutrientToDoublingTime["minimal"])
		conc_updates = self.getBiomassAsConcentrations(doubling_time)
		if self.use_trna_charging:
			conc_updates.update(self.update_amino_acid_targets(countsToMolar))
		if self.include_ppgpp:
			conc_updates[self.ppgpp_id] = self.getppGppConc(doubling_time)

		# Set molecule availability (internal and external)
		set_molecule_levels(self.fba, self.exchangeConstraints, coefficient,
			current_media_id, exchange_data, conc_updates, countsToMolar,
			metaboliteCountsInit, self.internalExchangeIdxs,
			self.aa_names_no_location, self.metaboliteNames)

		# Set reaction limits for maintenance and catalysts present
		translation_gtp = self._sim.processes["PolypeptideElongation"].gtp_to_hydrolyze
		catalyst_counts = self.catalysts.counts()
		set_reaction_bounds(self.fba, self.ngam, coefficient, countsToMolar, translation_gtp,
			self.reactions_with_catalyst, catalyst_counts, self.catalysisMatrix)

		# Constrain reactions based on targets
		kinetic_enzyme_counts = self.kineticsEnzymes.counts()
		kinetic_substrate_counts = self.kineticsSubstrates.counts()
		transport_targets = environment.transport_fluxes.values()
		if self.use_kinetics and self.burnInComplete:
			targets = set_reaction_targets(self.fba, kinetic_enzyme_counts,
				kinetic_substrate_counts, countsToMolar, self.getKineticConstraints,
				transport_targets, time_step, self.active_constraints_mask,
				self.all_constrained_reactions)
		else:
			targets = np.zeros(len(self.all_constrained_reactions))

		# Solve FBA problem and update metabolite counts
		deltaMetabolites = (1 / countsToMolar) * (CONC_UNITS * self.fba.getOutputMoleculeLevelsChange())

		metaboliteCountsFinal = np.zeros_like(metaboliteCountsInit)
		metaboliteCountsFinal[self.internalExchangeIdxs] = np.fmax(stochasticRound(
			self.randomState,
			metaboliteCountsInit[self.internalExchangeIdxs] + deltaMetabolites.asNumber()
			), 0).astype(np.int64)

		self.metabolites.countsIs(metaboliteCountsFinal)

		exchange_fluxes = CONC_UNITS * self.fba.getExternalExchangeFluxes()
		converted_exchange_fluxes = (exchange_fluxes / coefficient).asNumber(units.mmol / units.g / units.h)

		# update environmental nutrient counts
		delta_nutrients = ((1 / countsToMolar) * exchange_fluxes).asNumber().astype(int)
		environment.molecule_exchange(self.fba.getExternalMoleculeIDs(), delta_nutrients)

		import_exchange, import_constraint = environment.get_import_constraints(exchange_data)

		# Write outputs to listeners
		time_step_unitless = time_step.asNumber(TIME_UNITS)
		self.writeToListener("FBAResults", "import_exchange", import_exchange)
		self.writeToListener("FBAResults", "import_constraint", import_constraint)
		self.writeToListener("FBAResults", "deltaMetabolites", metaboliteCountsFinal - metaboliteCountsInit)
		self.writeToListener("FBAResults", "reactionFluxes", self.fba.getReactionFluxes() / time_step_unitless)
		self.writeToListener("FBAResults", "externalExchangeFluxes", converted_exchange_fluxes)
		self.writeToListener("FBAResults", "objectiveValue", self.fba.getObjectiveValue())
		self.writeToListener("FBAResults", "shadowPrices", self.fba.getShadowPrices(self.metaboliteNames))
		self.writeToListener("FBAResults", "reducedCosts", self.fba.getReducedCosts(self.fba.getReactionIDs()))
		self.writeToListener("FBAResults", "targetConcentrations", [self.homeostaticObjective[mol] for mol in self.fba.getHomeostaticTargetMolecules()])
		self.writeToListener("FBAResults", "homeostaticObjectiveValues", self.fba.getHomeostaticObjectiveValues())
		self.writeToListener("FBAResults", "kineticObjectiveValues", self.fba.getKineticObjectiveValues())
		self.writeToListener("EnzymeKinetics", "metaboliteCountsInit", metaboliteCountsInit)
		self.writeToListener("EnzymeKinetics", "metaboliteCountsFinal", metaboliteCountsFinal)
		self.writeToListener("EnzymeKinetics", "enzymeCountsInit", kinetic_enzyme_counts)
		self.writeToListener("EnzymeKinetics", "countsToMolar", countsToMolar.asNumber(CONC_UNITS))
		self.writeToListener("EnzymeKinetics", "actualFluxes", self.fba.getReactionFluxes(self.all_constrained_reactions) / time_step_unitless)
		self.writeToListener("EnzymeKinetics", "targetFluxes", targets / time_step_unitless)
		# TODO: add lower and upper targets

	def getBiomassAsConcentrations(self, doubling_time):
		'''
		Caches the result of the sim_data function to improve performance since
		function requires computation but won't change for a given doubling_time.

		Args:
			doubling_time (float with time units): doubling time of the cell to
				get the metabolite concentrations for

		Returns:
			dict {str : float with concentration units}: dictionary with metabolite
				IDs as keys and concentrations as values
		'''

		if doubling_time not in self._biomass_concentrations:
			self._biomass_concentrations[doubling_time] = self._getBiomassAsConcentrations(doubling_time)

		return self._biomass_concentrations[doubling_time]

	def update_amino_acid_targets(self, counts_to_molar):
		'''
		Finds new amino acid concentration targets based on difference in supply
		and number of amino acids used in polypeptide_elongation

		Args:
			counts_to_molar (float with mol/volume units): conversion from counts
				to molar for the current state of the cell

		Returns:
			dict {AA name (str): AA conc (float with mol/volume units)}:
				new concentration targets for each amino acid

		Skips updates to certain molecules defined in self.aa_targets_not_updated:
		- L-SELENOCYSTEINE: rare amino acid that led to high variability when updated

		TODO:
		- remove access to PolypeptideElongation class attribute (aa_count_diff)
		'''

		count_diff = self._sim.processes['PolypeptideElongation'].aa_count_diff

		if len(self.aa_targets):
			for aa, diff in count_diff.items():
				if aa in self.aa_targets_not_updated:
					continue
				self.aa_targets[aa] += diff
		# First time step of a simulation so set target to current counts to prevent
		# concentration jumps between generations
		else:
			for aa, counts in zip(self.aa_names, self.aas.total_counts()):
				if aa in self.aa_targets_not_updated:
					continue
				self.aa_targets[aa] = counts

		return {aa: counts * counts_to_molar for aa, counts in self.aa_targets.items()}

def update_external_molecule_levels(external_exchange_molecule_ids,
		aa_names_no_location, objective, metabolite_names,
		metabolite_concentrations, external_molecule_levels):
	"""
	Limit amino acid uptake to what is needed to meet concentration objective
	to prevent use as carbon source.

	Args:
		external_exchange_molecule_ids (Tuple[str]): IDs of molecules that
			are exchanged with the environment (with internal location tag)
		aa_names_no_location (List[str]): IDs of all amino acids without a
			location tag
		objective (Dict[str, Unum]): homeostatic objective for internal
			molecules (molecule ID: concentration in counts/volume units)
		metabolite_names (List[str]): IDs of all metabolites with a concentration
			target
		metabolite_concentrations (Unum[float]): concentration for each molecule
			in metabolite_names
		external_molecule_levels (np.ndarray[float]): current limits on external
			molecule availability

	Returns:
		external_molecule_levels (np.ndarray[float]): updated limits on external
			molcule availability
	"""

	for aa in aa_names_no_location:
		if aa + "[p]" in external_exchange_molecule_ids:
			idx = external_exchange_molecule_ids.index(aa + "[p]")
		elif aa + "[c]" in external_exchange_molecule_ids:
			idx = external_exchange_molecule_ids.index(aa + "[c]")
		else:
			continue

		concDiff = objective[aa + "[c]"] - metabolite_concentrations[metabolite_names.index(aa + "[c]")].asNumber(CONC_UNITS)
		if concDiff < 0:
			concDiff = 0

		if external_molecule_levels[idx] > concDiff:
			external_molecule_levels[idx] =  concDiff

	return external_molecule_levels

def set_molecule_levels(fba, exchange_constraints, coefficient, current_media_id,
		exchange_data, conc_updates, counts_to_molar, metabolite_counts,
		exchange_idxs, aa_names_no_location, metabolite_names):
	"""
	Set internal and external molecule levels available to the FBA solver.

	Args:
		fba (FluxBalanceAnalysis): FBA solver
		exchange_constraints (Callable): get external molecules and objective
			based on environment and cell state
		coefficient (Unum): coefficient to convert from mmol/g DCW/hr to mM basis
			(mass.time/volume units)
		current_media_id (str): ID of current media
		exchange_data (Dict[str, Any]): exchange data for the current environment
		conc_updates (Dict[str, Unum]): updates to concentrations targets for
			molecules (molecule ID: concentration in counts/volume units)
		counts_to_molar (Unum): conversion from counts to molar (counts/volume units)
		metabolite_counts (np.ndarray[int]): counts for each metabolite with a
			concentration target
		exchange_idxs (np.ndarray[int]): mapping of fba object metabolites to
			simulation metabolite counts
		aa_names_no_location (List[str]): IDs of all amino acids without a
			location tag
		metabolite_names (List[str]): IDs of all metabolites with a concentration
			target
	"""

	# Update objective from media exchanges
	# TODO: check if this is all exchange or only import?
	# - should it be one or the other for places used
	# - update arg doc string for function above
	external_exchange_molecule_ids = fba.getExternalMoleculeIDs()
	external_molecule_levels, objective = exchange_constraints(
		external_exchange_molecule_ids,
		coefficient,
		CONC_UNITS,
		current_media_id,
		exchange_data,
		conc_updates,
		)
	fba.update_homeostatic_targets(objective)

	# Internal concentrations
	metabolite_conc = counts_to_molar * metabolite_counts[exchange_idxs]
	fba.setInternalMoleculeLevels(metabolite_conc.asNumber(CONC_UNITS))

	# External concentrations
	external_molecule_levels = update_external_molecule_levels(
		external_exchange_molecule_ids, aa_names_no_location, objective,
		metabolite_names, metabolite_conc, external_molecule_levels)
	fba.setExternalMoleculeLevels(external_molecule_levels)

def set_reaction_bounds(fba, ngam, coefficient, counts_to_molar, gtp_to_hydrolyze,
		catalyzed_reactions, catalyst_counts, catalyst_matrix):
	"""
	Set reaction bounds for constrained reactions in the FBA object.

	Args:
		fba (FluxBalanceAnalysis): FBA solver
		ngam (Unum): non-growth associated maintenance (counts/mass/time units)
		coefficient (Unum): coefficient to convert from mmol/g DCW/hr to mM basis
			(mass.time/volume units)
		counts_to_molar (Unum): conversion from counts to molar (counts/volume units)
		gtp_to_hydrolyze (float): number of GTP molecules to hydrolyze to
			account for consumption in translation
		catalyzed_reactions (List[str]): reaction IDs for reactions that have
			an enzyme catalyst
		catalyst_counts (np.ndarray[int]): counts of enzyme catalysts
		catalyst_matrix (csr_matrix[int]): mapping of enzymes to reactions
			they catalyze (n reactions, m enzymes)
	"""

	# Maintenance reactions
	## Calculate new NGAM
	flux = (ngam * coefficient).asNumber(CONC_UNITS)
	fba.setReactionFluxBounds(
		fba._reactionID_NGAM,
		lowerBounds=flux, upperBounds=flux,
		)

	## Calculate GTP usage based on how much was needed in polypeptide
	## elongation in previous step.
	flux = (counts_to_molar * gtp_to_hydrolyze).asNumber(CONC_UNITS)
	fba.setReactionFluxBounds(
		fba._reactionID_polypeptideElongationEnergy,
		lowerBounds=flux, upperBounds=flux,
		)

	# Set hard upper bounds constraints based on enzyme presence
	# (infinite upper bound) or absence (upper bound of zero)
	reaction_bounds = np.inf * np.ones(len(catalyzed_reactions))
	no_rxn_mask = catalyst_matrix.dot(catalyst_counts) == 0
	reaction_bounds[no_rxn_mask] = 0
	fba.setReactionFluxBounds(catalyzed_reactions,
		upperBounds=reaction_bounds, raiseForReversible=False)

	# TODO: remove this variant and other attributes in class
	# if self.shuffleCatalyzedIdxs is not None:
	# 	catalyzedReactionBounds = catalyzedReactionBounds[self.shuffleCatalyzedIdxs]

def set_reaction_targets(fba, kinetic_enzyme_counts, kinetic_substrate_counts,
		counts_to_molar, get_kinetic_constraints, transport_targets, time_step,
		active_constraints_mask, constrained_reactions):
	"""
	Set reaction targets for constrained reactions in the FBA object.

	Args:
		fba (FluxBalanceAnalysis): FBA solver
		kinetic_enzyme_counts (np.ndarray[int]): counts of enzymes used in
			kinetic constraints
		kinetic_substrate_counts (np.ndarray[int]): counts of substrates used
			in kinetic constraints
		counts_to_molar (Unum): conversion from counts to molar (counts/volume units)
		get_kinetic_constraints (Callable[np.ndarray[float], np.ndarray[float]]):
			returns kinetic constraint targets for each reaction
		transport_targets (List[float]): targets for transport reactions
		time_step (Unum): current time step (time units)
		active_constraints_mask (np.ndarray[bool]): mask for all kinetic target
			reactions (True if active constraint, False if not)
		constrained_reactions (List[str]): reaction IDs for reactions with
			a target, includes both kinetics and transport

	Returns:
		mean_targets (np.ndarray[float]): mean target for each reaction in
			constrained_reactions
	"""

	# Unit basis for kinetic constraints
	# TODO: handle unit conversion in get_kinetic_constraints function?
	constraint_conc_units = units.umol / units.L
	constraint_time_units = units.s

	enzyme_conc = counts_to_molar * kinetic_enzyme_counts
	substrate_conc = counts_to_molar * kinetic_substrate_counts

	## Set target fluxes for reactions based on their most relaxed constraint
	reaction_targets = (constraint_conc_units / constraint_time_units) * get_kinetic_constraints(
		enzyme_conc.asNumber(constraint_conc_units),
		substrate_conc.asNumber(constraint_conc_units),
		)

	# TODO: remove this variant and other attributes in class
	# Shuffle parameters (only performed in very specific cases)
	# if self.shuffleIdxs is not None:
	# 	reaction_targets = (units.umol / units.L / units.s) * reaction_targets.asNumber()[self.shuffleIdxs, :]

	## Calculate reaction flux target for current time step
	targets = (time_step * reaction_targets).asNumber(CONC_UNITS)[active_constraints_mask, :]

	# add boundary targets
	lower_targets = np.concatenate((targets[:, 0], transport_targets), axis=0)
	mean_targets = np.concatenate((targets[:, 1], transport_targets), axis=0)
	upper_targets = np.concatenate((targets[:, 2], transport_targets), axis=0)

	## Set kinetic targets only if kinetics is enabled
	fba.set_scaled_kinetic_objective(time_step.asNumber(constraint_time_units))
	fba.setKineticTarget(
		constrained_reactions, mean_targets,
		lower_targets=lower_targets, upper_targets=upper_targets)

	return mean_targets
