from __future__ import absolute_import, division, print_function

import os
import csv
import time
import ast
from scipy import constants
import numpy as np

from reconstruction.spreadsheets import JsonReader
from itertools import ifilter
from wholecell.utils import units

from agent.inner import CellSimulation
from environment.make_rate_laws import KineticFluxModel


COUNTS_UNITS = units.mol
VOLUME_UNITS = units.L
MASS_UNITS = units.g
TIME_UNITS = units.s
CONC_UNITS = COUNTS_UNITS / VOLUME_UNITS
FLUX_UNITS = COUNTS_UNITS / VOLUME_UNITS / TIME_UNITS


TUMBLE_JITTER = 2.0 # (radians)
DEFAULT_COLOR = [color/255 for color in [255, 51, 51]]

CSV_DIALECT = csv.excel_tab
TRANSPORT_REACTIONS_FILE = os.path.join('environment', 'condition', 'look_up_tables', 'transport_reactions.tsv')
KINETIC_PARAMETERS_FILE = os.path.join('environment', 'condition', 'parameters', 'glt_family.tsv')
EXTERNAL_MOLECULES_FILE = os.path.join('environment', 'condition', 'environment_molecules.tsv')
# WCM_SIMDATA_FILE = os.path.join('environment', 'condition', 'look_up_tables', 'wcm_sim_data.json')

mM_to_M = 1E-3 # convert mmol/L to mol/L

class TransportKinetics(CellSimulation):
	'''	'''

	def __init__(self, state):
		self.initial_time = state.get('time', 0.0)
		self.local_time = state.get('time', 0.0)
		self.media_id = state.get('media_id', 'minimal')
		self.timestep = 1.0
		self.environment_change = {}
		self.volume = 1.0  # (fL)
		self.division_time = 100
		self.nAvogadro = constants.N_A * 1e-3  # convert 1/mol to 1/mmol.

		# initial state
		self.external_concentrations = {}
		self.internal_concentrations = {}
		self.motile_force = [0.01, 0.01] # initial magnitude and relative orientation
		self.division = []

		# make dict of transport reactions
		self.all_transport_reactions = {}
		with open(TRANSPORT_REACTIONS_FILE, 'rU') as csvfile:
			reader = JsonReader(csvfile, dialect=CSV_DIALECT)
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

		# Make kinetic_parameters in a nested format: {reaction_id: {transporter_id : {param_id: param_value}}}
		# TODO -- don't use literal_eval
		# TODO -- strip # before JsonReader
		self.kinetic_parameters = {}
		with open(KINETIC_PARAMETERS_FILE, 'rU') as csvfile:
			reader = JsonReader(csvfile, dialect=CSV_DIALECT)
			for row in reader:
				reaction_id = row['reaction id']
				transporter = row['transporter']
				k_avg = ast.literal_eval(row['k_avg'])
				max_conc = ast.literal_eval(row['max_conc'])

				# Combine kinetics into dictionary
				k_param = {'k_avg': k_avg}
				k_param.update(max_conc)
				transporter_kinetics = {transporter: k_param}

				if reaction_id in self.kinetic_parameters:
					self.kinetic_parameters[reaction_id].update(transporter_kinetics)
				else:
					self.kinetic_parameters[reaction_id] = transporter_kinetics

		# Get list of external molecules
		# TODO -- bring back # in EXTERNAL_MOLECULES_FILE, strip them before JsonReader
		self.external_molecule_ids = []
		with open(EXTERNAL_MOLECULES_FILE, 'rU') as csvfile:
			# reader = JsonReader(csvfile, dialect=CSV_DIALECT)
			reader = JsonReader(
				ifilter(lambda x: x.lstrip()[0] != "#", csvfile), # Strip comments
				dialect = CSV_DIALECT)
			for row in reader:
				self.external_molecule_ids.append(row['molecule id'])

		# with open(WCM_SIMDATA_FILE, 'r') as f:
		# 	wcm_sim_out = json.loads(f.read())
		#

		# make the kinetic model
		self.kinetic_rate_laws = KineticFluxModel(self.all_transport_reactions, self.kinetic_parameters)

		# Get list of molecule_ids used by kinetic rate laws
		self.molecule_ids = self.kinetic_rate_laws.molecule_ids

		# Get internal molecule ids by removing external molecule ids
		self.internal_molecule_ids = [mol_id for mol_id in self.molecule_ids if mol_id not in self.external_molecule_ids]
		self.all_molecule_ids = self.internal_molecule_ids + self.external_molecule_ids

		# TODO -- generate concentration dict from WCM data
		self.concentrations = {molecule_id: 1.0 for molecule_id in self.all_molecule_ids}

		# get initial fluxes
		self.transport_fluxes = self.kinetic_rate_laws.get_fluxes(self.concentrations)





		# TODO -- merge external and internal concentrations.









	def update_state(self):
		# nAvogadro is in 1/mol --> convert to 1/mmol. volume is in fL --> convert to L
		self.molar_to_counts = (self.nAvogadro) * (self.volume * 1e-15)

		self.transport_fluxes = self.kinetic_rate_laws.get_fluxes(self.concentrations)

		delta_counts = self.flux_to_counts(self.transport_fluxes)

		environment_deltas = {}
		for molecule in self.external_concentrations.keys():
			# TODO -- use external exchange map rather than (molecule + '[p]')
			molecule_p = molecule + '[p]'
			if molecule_p in delta_counts:
				environment_deltas[molecule] = delta_counts[molecule_p]

		# accumulate in environment_change
		self.accumulate_deltas(environment_deltas)

	def accumulate_deltas(self, environment_deltas):
		for molecule_id, count in environment_deltas.iteritems():
			self.environment_change[molecule_id] += count

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

		# reset environment change
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
		return {
			'volume': self.volume,
			'motile_force': self.motile_force,
			'environment_change': self.environment_change,
			'division': self.division,
			'color': DEFAULT_COLOR,
			'transport_fluxes': self.transport_fluxes,
			}


	## Flux-related functions
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
