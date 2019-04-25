from __future__ import absolute_import, division, print_function

import time
import numpy as np
from scipy import constants

from agent.inner import CellSimulation

# Raw data class
from reconstruction.ecoli.knowledge_base_raw import KnowledgeBaseEcoli

TUMBLE_JITTER = 2.0 # (radians)
DEFAULT_COLOR = [color/255 for color in [255, 69, 0]]

variants = [
	'minimal',
	'minimal_minus_oxygen',
	'minimal_plus_amino_acids']

amino_acids = [
	'L-ALPHA-ALANINE',
	'ARG',
	'ASN',
	'L-ASPARTATE',
	'CYS',
	'GLT',
	'GLN',
	'GLY',
	'HIS',
	'ILE',
	'LEU',
	'LYS',
	'MET',
	'PHE',
	'PRO',
	'SER',
	'THR',
	'TRP',
	'TYR',
	'L-SELENOCYSTEINE',
	'VAL'
]

class TransportProcess(CellSimulation):
	''''''

	def __init__(self):
		self.initial_time = 0.0
		self.local_time = 0.0
		# self.timestep = 1.0
		self.environment_change = {}
		self.volume = 1.0  # (fL)
		self.division_time = 100
		self.nAvogadro = constants.N_A

		# initial state
		self.external_concentrations = {}
		self.internal_concentrations = {}
		self.motile_force = [0.01, 0.01] # initial magnitude and relative orientation
		self.division = []

		raw_data = KnowledgeBaseEcoli()
		self.all_transport_reactions = {}
		for row in raw_data.transport_reactions:
			reaction_id = row["reaction id"]
			stoichiometry = row["stoichiometry"]
			reversible = row["is reversible"]
			catalyzed = row["catalyzed by"]
			self.all_transport_reactions[reaction_id] = {
				"stoichiometry": stoichiometry,
				"is reversible": reversible,
				"catalyzed by": catalyzed,
			}

		# make a dictionary with saved average fluxes for all transport reactions, in the three conditions
		# fluxes are in mmol/L
		self.flux_lookup = {variant: {} for variant in variants}
		for row in raw_data.transport_data.transport_avg_minimal:
			reaction_id = row["reaction id"]
			flux = row["average flux"]
			self.flux_lookup['minimal'][reaction_id] = flux
		for row in raw_data.transport_data.transport_avg_minimal_minus_oxygen:
			reaction_id = row["reaction id"]
			flux = row["average flux"]
			self.flux_lookup['minimal_minus_oxygen'][reaction_id] = flux
		for row in raw_data.transport_data.transport_avg_minimal_plus_amino_acids:
			reaction_id = row["reaction id"]
			flux = row["average flux"]
			self.flux_lookup['minimal_plus_amino_acids'][reaction_id] = flux

		# exchange_ids declares which molecules' exchange will be controlled by transport
		aa_p_ids = [aa_id + "[p]" for aa_id in amino_acids]
		exchange_molecules = ["OXYGEN-MOLECULE[p]", "GLC[p]"]
		exchange_ids = exchange_molecules + aa_p_ids
		self.transport_reactions_ids = self.reactions_from_exchange(exchange_ids)

		# TODO (Eran) -- use synchronize to get media_id upon initialization
		self.media_id = 'minimal'
		self.current_flux_lookup = self.flux_lookup[self.media_id]
		self.transport_fluxes = self.get_fluxes(self.current_flux_lookup, self.transport_reactions_ids)


	def update_state(self):
		# nAvogadro is in 1/mol --> convert to 1/mmol. volume is in fL --> convert to L
		self.molar_to_counts = (self.nAvogadro * 1e-3) * (self.volume * 1e-15)
		self.transport_fluxes = self.get_fluxes(self.current_flux_lookup, self.transport_reactions_ids)
		delta_counts = self.flux_to_counts(self.transport_fluxes)

		for molecule in self.external_concentrations.keys():
			# TODO -- use external exchange map rather than (molecule + '[p]')
			if (molecule + '[p]') in delta_counts:
				self.environment_change[molecule] = delta_counts[molecule + '[p]']

	def check_division(self):
		# update division state based on time since initialization
		if self.local_time >= self.initial_time + self.division_time:
			self.division = [{'time': self.local_time}, {'time': self.local_time}]
		return self.division

	def time(self):
		return self.local_time

	def apply_outer_update(self, update):
		self.external_concentrations = update['concentrations']
		self.media_id = update['media_id']

		# update lookup table
		self.current_flux_lookup = self.flux_lookup[self.media_id]

		# reset environment change
		self.environment_change = {}
		for molecule in self.external_concentrations.iterkeys():
			self.environment_change[molecule] = 0

	def run_incremental(self, run_until):
		# TODO -- implement time steps! Will require an accumulate_deltas() as in local_environment.py
		# update state once per message exchange
		self.update_state()
		# self.check_division()
		self.local_time = run_until

		time.sleep(1.0)  # pause for better coordination with Lens visualization. TODO: remove this

	def generate_inner_update(self):
		return {
			'volume': self.volume,
			'motile_force': self.motile_force,
			'environment_change': self.environment_change,
			'division': self.division,
			'color': DEFAULT_COLOR,
			'transport_fluxes': self.transport_fluxes,
			}

	def synchronize_state(self, state):
		if 'time' in state:
			self.initial_time = state['time']


	## Flux-related functions
	def get_fluxes(self, flux_lookup, transport_reactions_ids):
		transport_fluxes = {transport_id: flux_lookup[transport_id] for transport_id in transport_reactions_ids}
		return transport_fluxes

	def flux_to_counts(self, fluxes):
		rxn_counts = {reaction_id: int(self.molar_to_counts * flux) for reaction_id, flux in fluxes.iteritems()}
		delta_counts = {}
		for reaction_id, rxn_count in rxn_counts.iteritems():
			stoichiometry = self.all_transport_reactions[reaction_id]['stoichiometry']
			substrate_counts = {substrate_id: coeff * rxn_count for substrate_id, coeff in stoichiometry.iteritems()}
			# add to delta_counts
			for substrate, delta in substrate_counts.iteritems():
				if substrate in delta_counts:
					delta_counts[substrate] += delta
				else:
					delta_counts[substrate] = delta
		return delta_counts

	def reactions_from_exchange(self, include_exchanges):
		include_reactions = []
		for reaction_id, specs in self.all_transport_reactions.iteritems():
			reaction_molecules = specs['stoichiometry'].keys()
			for exchange in include_exchanges:
				if exchange in reaction_molecules:
					include_reactions.append(reaction_id)
		return include_reactions
