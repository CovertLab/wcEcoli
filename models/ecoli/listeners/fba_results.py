"""
FBAResults

Records dynamics of FBA output.
"""

import numpy as np
from scipy.sparse import csr_matrix

import wholecell.listeners.listener
from reconstruction.ecoli.dataclasses.process.metabolism import REVERSE_TAG

class FBAResults(wholecell.listeners.listener.Listener):
	""" FBAResults """

	_name = "FBAResults"

	# Constructor
	def __init__(self, *args, **kwargs):
		super(FBAResults, self).__init__(*args, **kwargs)


	# Construct object graph
	def initialize(self, sim, sim_data):
		super(FBAResults, self).initialize(sim, sim_data)

		# Read required values
		self.metabolism = sim.processes["Metabolism"].model
		self.media_id = sim.external_states['Environment'].current_media_id
		fba = self.metabolism.fba
		self.fba_reaction_ids = fba.getReactionIDs()
		self.base_reaction_ids = sim_data.process.metabolism.base_reaction_ids
		fba_reaction_ids_to_base_reaction_ids = sim_data.process.metabolism.reaction_id_to_base_reaction_id
		self.externalMoleculeIDs = fba.getExternalMoleculeIDs()
		self.outputMoleculeIDs = fba.getOutputMoleculeIDs()
		self.kineticTargetFluxNames = fba.getKineticTargetFluxNames()
		self.homeostaticTargetMolecules = fba.getHomeostaticTargetMolecules()

		self.objectiveValue = 0.0

		# exchange with environment
		self.all_external_exchange_molecules = sim_data.external_state.all_external_exchange_molecules
		self.conc_update_molecules = sim.processes["Metabolism"].conc_update_molecules

		# Get conversion matrix to compile individual fluxes in the FBA
		# solution to the fluxes of base reactions
		fba_reaction_id_to_index = {
			rxn_id: i for (i, rxn_id) in enumerate(self.fba_reaction_ids)
			}
		base_reaction_id_to_index = {
			rxn_id: i for (i, rxn_id) in enumerate(self.base_reaction_ids)
			}
		base_rxn_indexes = []
		fba_rxn_indexes = []
		v = []

		for fba_rxn_id in self.fba_reaction_ids:
			base_rxn_id = fba_reaction_ids_to_base_reaction_ids[fba_rxn_id]
			base_rxn_indexes.append(base_reaction_id_to_index[base_rxn_id])
			fba_rxn_indexes.append(fba_reaction_id_to_index[fba_rxn_id])
			if fba_rxn_id.endswith(REVERSE_TAG):
				v.append(-1)
			else:
				v.append(1)

		base_rxn_indexes = np.array(base_rxn_indexes)
		fba_rxn_indexes = np.array(fba_rxn_indexes)
		v = np.array(v)
		shape = (len(self.base_reaction_ids), len(self.fba_reaction_ids))

		self.reaction_mapping_matrix = csr_matrix(
			(v, (base_rxn_indexes, fba_rxn_indexes)),
			shape=shape)

	# Allocate memory
	def allocate(self):
		super(FBAResults, self).allocate()

		self.reactionFluxes = np.zeros(len(self.fba_reaction_ids), np.float64)
		self.base_reaction_fluxes = np.zeros(len(self.base_reaction_ids), np.float64)
		self.externalExchangeFluxes = np.zeros(len(self.externalMoleculeIDs), np.float64)
		self.shadowPrices = np.zeros(len(self.outputMoleculeIDs), np.float64)
		self.reducedCosts = np.zeros(len(self.fba_reaction_ids), np.float64)
		self.homeostaticObjectiveValues = np.zeros(len(self.homeostaticTargetMolecules))
		self.kineticObjectiveValues = np.zeros(len(self.kineticTargetFluxNames))
		self.deltaMetabolites = np.zeros(len(self.metabolism.metaboliteNamesFromNutrients), np.float64)
		self.targetConcentrations = np.zeros(len(self.homeostaticTargetMolecules))

		# Args for metabolism functions
		self.conc_updates = np.zeros(len(self.conc_update_molecules))
		self.catalyst_counts = np.zeros(len(self.metabolism.catalyst_ids))
		self.translation_gtp = 0.
		self.coefficient = 0.
		self.unconstrained_molecules = [False] * len(self.all_external_exchange_molecules)
		self.constrained_molecules = [False] * len(self.all_external_exchange_molecules)
		self.uptake_constraints = [np.nan] * len(self.all_external_exchange_molecules)

	def update(self):
		# Compile reaction fluxes to those of base reactions
		self.base_reaction_fluxes = self.reaction_mapping_matrix.dot(
			self.reactionFluxes)

	def tableCreate(self, tableWriter):
		subcolumns = {
			'reactionFluxes': 'reactionIDs',
			'base_reaction_fluxes': 'base_reaction_ids',
			'externalExchangeFluxes': 'externalMoleculeIDs',
			'shadowPrices': 'outputMoleculeIDs',
			'reducedCosts': 'reactionIDs',
			'homeostaticObjectiveValues': 'homeostaticTargetMolecules',
			'kineticObjectiveValues': 'kineticTargetFluxNames',
			'deltaMetabolites': 'metaboliteNames',
			'targetConcentrations': 'homeostaticTargetMolecules',
			'importConstraint': 'all_external_exchange_molecules',
			'importExchange': 'all_external_exchange_molecules',
			'conc_updates': 'conc_update_molecules',
			'catalyst_counts': 'catalyst_ids',
			}

		tableWriter.writeAttributes(
			reactionIDs=list(self.fba_reaction_ids),
			base_reaction_ids=self.base_reaction_ids,
			externalMoleculeIDs=self.externalMoleculeIDs,
			outputMoleculeIDs=self.outputMoleculeIDs,
			homeostaticTargetMolecules=self.homeostaticTargetMolecules,
			kineticTargetFluxNames=self.kineticTargetFluxNames,
			metaboliteNames=self.metabolism.metaboliteNamesFromNutrients,
			all_external_exchange_molecules=self.all_external_exchange_molecules,
			conc_update_molecules=self.conc_update_molecules,
			catalyst_ids=self.metabolism.catalyst_ids,
			subcolumns=subcolumns,
			)


	def tableAppend(self, tableWriter):
		tableWriter.append(
			time=self.time(),
			simulationStep=self.simulationStep(),
			reactionFluxes=self.reactionFluxes,
			base_reaction_fluxes=self.base_reaction_fluxes,
			externalExchangeFluxes=self.externalExchangeFluxes,
			shadowPrices=self.shadowPrices,
			reducedCosts=self.reducedCosts,
			objectiveValue=self.objectiveValue,
			homeostaticObjectiveValues=self.homeostaticObjectiveValues,
			kineticObjectiveValues=self.kineticObjectiveValues,
			deltaMetabolites=self.deltaMetabolites,
			targetConcentrations=self.targetConcentrations,
			media_id=self.media_id,
			conc_updates=self.conc_updates,
			catalyst_counts=self.catalyst_counts,
			translation_gtp=self.translation_gtp,
			coefficient=self.coefficient,
			unconstrained_molecules=self.unconstrained_molecules,
			constrained_molecules=self.constrained_molecules,
			uptake_constraints=self.uptake_constraints,
			)
	

	def get_dict(self):
		return {
			'fba_results': {
				'reaction_fluxes': self.reactionFluxes,
				'base_reaction_fluxes': self.base_reaction_fluxes,
				'external_exchange_fluxes': self.externalExchangeFluxes,
				'shadow_prices': self.shadowPrices,
				'reduced_costs': self.reducedCosts,
				'objective_value': self.objectiveValue,
				'homeostatic_objective_values': self.homeostaticObjectiveValues,
				'kinetic_objective_values': self.kineticObjectiveValues,
				'delta_metabolites': self.deltaMetabolites,
				'target_concentrations': self.targetConcentrations,
				'media_id': self.media_id,
				'conc_updates': self.conc_updates,
				'catalyst_counts': self.catalyst_counts,
				'translation_gtp': self.translation_gtp,
				'coefficient': self.coefficient,
				'unconstrained_molecules': self.unconstrained_molecules,
				'constrained_molecules': self.constrained_molecules,
				'uptake_constraints': self.uptake_constraints,
			}
		}
