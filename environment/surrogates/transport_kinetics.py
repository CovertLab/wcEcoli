from __future__ import absolute_import, division, print_function

import os
import csv
import time
from scipy import constants
import json

from reconstruction.spreadsheets import JsonReader
from itertools import ifilter

from agent.inner import CellSimulation
from environment.kinetic_rate_laws.kinetic_rate_laws import KineticFluxModel


TUMBLE_JITTER = 2.0 # (radians)
DEFAULT_COLOR = [color/255 for color in [255, 51, 51]]

CSV_DIALECT = csv.excel_tab
LOOKUP_DIR = os.path.join('environment', 'condition', 'look_up_tables')
TRANSPORT_REACTIONS_FILE = os.path.join(LOOKUP_DIR, 'transport_reactions.tsv')
CONC_LOOKUP_MINIMAL = os.path.join(LOOKUP_DIR, 'avg_concentrations', 'minimal.tsv')
CONC_LOOKUP_ANAEROBIC = os.path.join(LOOKUP_DIR, 'avg_concentrations', 'minimal_minus_oxygen.tsv')
CONC_LOOKUP_AA = os.path.join(LOOKUP_DIR, 'avg_concentrations', 'minimal_plus_amino_acids.tsv')
KINETIC_PARAMETERS_FILE = os.path.join('environment', 'kinetic_rate_laws', 'parameters', 'glt.json')
EXTERNAL_MOLECULES_FILE = os.path.join('environment', 'condition', 'environment_molecules.tsv')

class TransportKinetics(CellSimulation):
	'''
	A surrogate that uses kinetic rate laws to determine transport flux

	# TODO (Eran) -- bring back units management
	'''

	def __init__(self, state):
		self.initial_time = state.get('time', 0.0)
		self.local_time = state.get('time', 0.0)
		self.media_id = state.get('media_id', 'minimal')
		self.timestep = 1.0
		self.environment_change = {}
		self.volume = 1.0  # (fL) TODO (Eran) volume needs to change for transport fluxes to translate to increasing delta counts
		self.division_time = 100
		self.nAvogadro = constants.N_A

		# Initial state
		self.external_concentrations = {}
		self.internal_concentrations = {}
		self.motile_force = [0.01, 0.01] # initial magnitude and relative orientation
		self.division = []

		# Make dict of transport reactions
		self.all_transport_reactions = {}
		with open(TRANSPORT_REACTIONS_FILE, 'rU') as csvfile:
			reader = JsonReader(
				ifilter(lambda x: x.lstrip()[0] != "#", csvfile), # Strip comments
				dialect = CSV_DIALECT)
			for row in reader:
				reaction_id = row['reaction id']
				stoichiometry = row['stoichiometry']
				reversible = row['is reversible']
				catalyzed = row['catalyzed by']
				self.all_transport_reactions[reaction_id] = {
					'stoichiometry': stoichiometry,
					'is reversible': reversible,
					'catalyzed by': catalyzed,
				}

		# Make map of external molecule_ids with a location tag (as used in reaction stoichiometry) to molecule_ids in the environment
		self.molecule_to_external_map = {}
		self.external_to_molecule_map = {}
		with open(EXTERNAL_MOLECULES_FILE, 'rU') as csvfile:
			reader = JsonReader(
				ifilter(lambda x: x.lstrip()[0] != "#", csvfile), # Strip comments
				dialect = CSV_DIALECT)
			for row in reader:
				molecule_id = row['molecule id']
				location = row['exchange molecule location']
				self.molecule_to_external_map[molecule_id + location] = molecule_id
				self.external_to_molecule_map[molecule_id] = molecule_id + location

		# Load saved wcEcoli concentrations of molecules from minimal condition
		wcm_concs = {}
		with open(CONC_LOOKUP_AA, 'rU') as csvfile:
			reader = JsonReader(
				ifilter(lambda x: x.lstrip()[0] != "#", csvfile),  # Strip comments
				dialect=csv.excel_tab)
			for row in reader:
				molecule_id = row['molecule id']
				avg_conc = row['average concentration mmol/L']
				wcm_concs[molecule_id] = avg_conc

		# Load dict of saved parameters
		with open(KINETIC_PARAMETERS_FILE, 'r') as fp:
			kinetic_parameters = json.load(fp)

		# List of reactions to construct
		# This is used by transport_composite to set boundary_fluxes in wcEcoli
		self.transport_reactions_ids = kinetic_parameters.keys()

		# Make dict for all reactions in transport_reactions_ids
		kinetic_reactions = {
			reaction_id: specs
			for reaction_id, specs in self.all_transport_reactions.iteritems()
			if reaction_id in self.transport_reactions_ids}

		# Make the kinetic model
		self.kinetic_rate_laws = KineticFluxModel(kinetic_reactions, kinetic_parameters)

		# Get list of molecule_ids used by kinetic rate laws
		# This is used by transport_composite to set boundary_views in wcEcoli
		self.molecule_ids = self.kinetic_rate_laws.molecule_ids

		# Get concentrations of all molecule_ids from wcm
		self.concentrations = {mol_id: wcm_concs[mol_id] for mol_id in self.molecule_ids}
		self.concentrations['GLT[p]'] = 0.6  # mmol/L

		# Set initial fluxes
		self.transport_fluxes = self.kinetic_rate_laws.get_fluxes(self.concentrations)


	def update_state(self):
		# nAvogadro is in 1/mol --> convert to 1/mmol. volume is in fL --> convert to L
		self.millimolar_to_counts = (self.nAvogadro * 1e-3) * (self.volume * 1e-15)

		# Get transport fluxes, convert to change in counts
		self.transport_fluxes = self.kinetic_rate_laws.get_fluxes(self.concentrations)
		delta_counts = self.flux_to_counts(self.transport_fluxes)

		# Get the deltas for environmental molecules
		environment_deltas = {}
		for molecule_id in delta_counts.keys():
			if molecule_id in self.molecule_to_external_map:
				external_molecule_id = self.molecule_to_external_map[molecule_id]
				environment_deltas[external_molecule_id] = delta_counts[molecule_id]

		# Accumulate in environment_change
		self.accumulate_deltas(environment_deltas)

	def accumulate_deltas(self, environment_deltas):
		for molecule_id, count in environment_deltas.iteritems():
			self.environment_change[molecule_id] += count

	def check_division(self):
		# Update division state based on time since initialization
		if self.local_time >= self.initial_time + self.division_time:
			self.division = [{'time': self.local_time}, {'time': self.local_time}]
		return self.division

	def time(self):
		return self.local_time

	def apply_outer_update(self, update):
		self.media_id = update['media_id']
		self.external_concentrations = update.get('concentrations', {})
		boundary_concentrations = update.get('boundary_view', {})  # in mmol/L

		# Map from external_id to concentration key
		new_concentrations = {
			self.external_to_molecule_map[mol_id]: conc
			for mol_id, conc in self.external_concentrations.iteritems()}

		# Update concentrations dict
		self.concentrations.update(boundary_concentrations)
		self.concentrations.update(new_concentrations)

		# Reset environment change
		self.environment_change = {}
		for molecule in self.external_concentrations.iterkeys():
			self.environment_change[molecule] = 0

	def run_incremental(self, run_until):
		'''run until run_until'''
		while self.time() < run_until:
			self.local_time += self.timestep
			self.update_state()
			# self.check_division()

		time.sleep(1.0)  # pause for better coordination with Lens visualization. TODO: remove this

	def generate_inner_update(self):
		# Round off changes in counts
		self.environment_change = {mol_id: int(counts) for mol_id, counts in self.environment_change.iteritems()}
		return {
			'volume': self.volume,
			'motile_force': self.motile_force,
			'environment_change': self.environment_change,
			'division': self.division,
			'color': DEFAULT_COLOR,
			'transport_fluxes': self.transport_fluxes,
			}

	# TODO (eran) -- move this function to rate_law_utilities
	## Flux-related functions
	def flux_to_counts(self, fluxes):
		# rxn_counts are not rounded off here, need to be rounded off before generate_inner_update is sent
		rxn_counts = {reaction_id: self.millimolar_to_counts * flux for reaction_id, flux in fluxes.iteritems()}
		delta_counts = {}
		for reaction_id, rxn_count in rxn_counts.iteritems():
			stoichiometry = self.all_transport_reactions[reaction_id]['stoichiometry']
			substrate_counts = {substrate_id: coeff * rxn_count for substrate_id, coeff in stoichiometry.iteritems()}
			# Add to delta_counts
			for substrate, delta in substrate_counts.iteritems():
				if substrate in delta_counts:
					delta_counts[substrate] += delta
				else:
					delta_counts[substrate] = delta
		return delta_counts
