from __future__ import absolute_import, division, print_function

import os
import csv
import time
import ast
from scipy import constants
import numpy as np

from agent.inner import CellSimulation
from reconstruction.spreadsheets import JsonReader


from wholecell.utils import units

COUNTS_UNITS = units.mol
VOLUME_UNITS = units.L
MASS_UNITS = units.g
TIME_UNITS = units.s
CONC_UNITS = COUNTS_UNITS / VOLUME_UNITS
FLUX_UNITS = COUNTS_UNITS / VOLUME_UNITS / TIME_UNITS


TUMBLE_JITTER = 2.0 # (radians)
DEFAULT_COLOR = [color/255 for color in [255, 51, 51]]

CSV_DIALECT = csv.excel_tab
TRANSPORT_REACTIONS_FILE = os.path.join("environment", "condition", "look_up_tables", "transport_reactions.tsv")
KINETIC_PARAMETERS_FILE = os.path.join("environment", "condition", "parameters", "glt_family.tsv")
EXTERNAL_MOLECULES_FILE = os.path.join("environment", "condition", "environment_molecules.tsv")

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

mM_to_M = 1E-3 # convert mmol/L to mol/L

class TransportKinetics(CellSimulation):
	''''''

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
				reaction_id = row["reaction id"]
				stoichiometry = row["stoichiometry"]
				reversible = row["is reversible"]
				catalyzed = row["catalyzed by"]
				self.all_transport_reactions[reaction_id] = {
					"stoichiometry": stoichiometry,
					"is reversible": reversible,
					"catalyzed by": catalyzed,
				}

		# Make kinetic_parameters in a nested format: {reaction_id: {transporter_id : {param_id: param_value}}}
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
		self.external_molecule_ids = []
		with open(EXTERNAL_MOLECULES_FILE, 'rU') as csvfile:
			reader = JsonReader(csvfile, dialect=CSV_DIALECT)
			for row in reader:
				self.external_molecule_ids.append(row["molecule id"])


		# Get list of reaction_ids and molecule_ids
		self.transport_reactions_ids = self.kinetic_parameters.keys()  # use all kinetic parameters
		self.molecule_ids = self._get_molecules(self.transport_reactions_ids)

		# Get internal molecule ids by removing external molecule ids
		self.internal_molecule_ids = [mol_id for mol_id in self.molecule_ids if mol_id not in self.external_molecule_ids]
		self.all_molecule_ids = self.internal_molecule_ids + self.external_molecule_ids

		# # Get molecule IDs
		# bulk_molecule_ids = sim_data.internal_state.bulkMolecules.bulkData['id']
		# self.internal_molecule_ids = [mol_id for mol_id in internal_molecule_ids if mol_id in bulk_molecule_ids]
		# self.all_molecule_ids = self.internal_molecule_ids + self.external_molecule_ids

		# Get indices, for reading out arrays in calculateRequest and evolveState
		self.internal_molecule_indices = [index for index, mol_id in enumerate(self.all_molecule_ids) if mol_id in self.internal_molecule_ids]
		self.external_molecule_indices = [index for index, mol_id in enumerate(self.all_molecule_ids) if mol_id in self.external_molecule_ids]


		# Build rate laws
		self.rate_laws = self._make_rate_laws()


		# get initial fluxes
		# self.transport_fluxes = self.get_fluxes(self.current_flux_lookup, self.transport_reactions_ids)


	def update_state(self):
		# nAvogadro is in 1/mol --> convert to 1/mmol. volume is in fL --> convert to L
		self.molar_to_counts = (self.nAvogadro) * (self.volume * 1e-15)

		self.transport_fluxes = []

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

	def get_fluxes(self, concentrations_dict):
		'''
		Use rate law functions to calculate flux
		Args:
			concentrations_dict (dict) {molecule_id: concentration} - a dictionary of all relevant
		molecules and their concentrations, in mmol/L.
		Returns:
			reaction_fluxes_array - an array with fluxes for all reactions in reaction_id
			exchange_fluxes_array - an array with exchange fluxes for all molecules in molecule_id
		'''

		# Initialize reaction_fluxes and exchange_fluxes dictionaries
		reaction_fluxes = {reaction["reaction id"]: 0.0 for reaction in self.transport_reactions}
		exchange_fluxes = {mol_id: 0.0 for mol_id in self.all_molecule_ids}

		# Save reaction flux and add to exchange flux.
		for reaction in self.transport_reactions:
			reaction_id = reaction["reaction id"]
			stoichiometry = reaction["stoichiometry"]

			# Get transporter ids with compartments, and find the flux through each
			transporters = self.kinetic_parameters[reaction_id].keys()
			for transporter in transporters:

				flux = self.rate_laws[reaction_id][transporter](concentrations_dict)
				reaction_fluxes[reaction_id] += flux

				for molecule, coefficient in stoichiometry.iteritems():
					# Only use molecules in molecule_ids
					if molecule in self.all_molecule_ids:
						if coefficient < 0:
							exchange_fluxes[molecule] -= flux
						elif coefficient > 0:
							exchange_fluxes[molecule] += flux

		reaction_fluxes_array = np.asarray([reaction_fluxes[reaction_id] for reaction_id in self.reaction_ids])
		exchange_fluxes_array = np.asarray([exchange_fluxes[molecule_id] for molecule_id in self.all_molecule_ids])

		return reaction_fluxes_array, exchange_fluxes_array




	## Configure reactions
	def _get_molecules(self, reaction_ids):
		'''
		Inputs:
			reaction_ids - a list of all reaction ids that will be used by transport
		Returns:
			self.molecule_ids - a list of all molecules used by these reactions
		'''

		molecule_ids = []
		for reaction_id, specs in self.all_transport_reactions.iteritems():
			if reaction_id in reaction_ids:
				stoichiometry = specs["stoichiometry"]
				substrates = stoichiometry.keys()
				transporters = self.kinetic_parameters[reaction_id].keys()

				# Add all relevant molecules_ids
				molecule_ids.extend(substrates)
				molecule_ids.extend(transporters)

		# Save unique molecule ids
		return list(set(molecule_ids))

	def _make_rate_laws(self):
		'''
		Make rate laws for all reactions
		Returns:
			rate_laws (dict) {'reaction id': {transporter: rate_law}} - A dictionary with
				rate laws saved as functions, in nested.
		'''

		# Make rate laws for each reaction
		rate_laws = {reaction["reaction id"]: {} for reaction in self.transport_reactions}
		for reaction in self.transport_reactions:

			reaction_id = reaction["reaction id"]
			transporters = self.kinetic_parameters[reaction_id]

			# Rate law for each transporter
			for transporter, params in transporters.iteritems():

				if transporter == 'None':
					del params['k_avg']
					# If there are remaining parameters (max_conc), build rate law
					if params:
						rate_law = self._construct_rate_law_no_transporter(params)

						# Save the rate law in this reaction under transporter 'None'
						rate_laws[reaction_id][transporter] = rate_law

						print('{}, has no transporter'.format(reaction_id))
					else:
						rate_law = self._construct_empty_rate_law()
						# Save the rate law in this reaction under transporter 'None'
						rate_laws[reaction_id][transporter] = rate_law

						print('{}, has no transporter, no parameters'.format(reaction_id))

				else:
					rate_law = self._construct_rate_law_piecewise_linear(
						transporter,
						params
					)

					# Save the rate law for each transporter in this reaction
					rate_laws[reaction_id][transporter] = rate_law

		return rate_laws

	def _construct_empty_rate_law(self):
		def rate_law(concentrations):
			return 0.0

		return rate_law

	def _construct_rate_law_piecewise_linear(self, transporter, parameters):
		'''
		Makes a single piecewise linear rate law.
		Inputs:
			transporter (str)
			parameters (dict) {parameter: value} - requires a 'k_avg' entry, and max_conc entries for all
				relevant substrates.
		Returns:
			rate_law - a function for the kinetic rate law of the given transporter, which takes a
				dictionary of concentrations as input, and which returns flux as output
		'''

		k_avg = parameters.get('k_avg')
		max_conc = dict(parameters)
		del max_conc['k_avg']

		# add units. max_conc is in mmol/L --> convert to mol/L
		max_conc = {mol_id: value * mM_to_M for mol_id, value in max_conc.iteritems()}

		# Construct rate law for this transporter
		def rate_law(concentrations):
			transporter_concentration = concentrations[transporter]

			substrates_normalized = 1.0
			for substrate, substrate_max_conc in max_conc.iteritems():
				# If concentration is above the max, assume saturated
				if concentrations[substrate] > substrate_max_conc:
					continue
				else:
					substrates_normalized *= concentrations[substrate] / substrate_max_conc

			flux = k_avg * transporter_concentration * substrates_normalized
			return flux

		return rate_law

	def _construct_rate_law_no_transporter(self, parameters):
		'''
		Makes a single piecewise linear rate law, based entirely on substrate concentrations.
		Inputs:
			parameters (dict) {parameter: value} - requires a 'k_avg' entry, and max_conc
				entries for all relevant substrates.
		Returns:
			rate_law - a function for the kinetic rate law, which takes a
				dictionary of concentrations as input, and which returns flux as output
		'''

		max_conc = dict(parameters)

		# add units. max_conc is in mmol/L --> convert to mol/L
		max_conc = {mol_id: value * mM_to_M for mol_id, value in max_conc.iteritems()}

		# Construct rate law for this transporter
		def rate_law(concentrations):
			substrates_normalized = 1.0
			for substrate, substrate_max_conc in max_conc.iteritems():
				# If concentration is above the max, assume saturated
				if concentrations[substrate] > substrate_max_conc:
					continue
				else:
					substrates_normalized *= concentrations[substrate] / substrate_max_conc

			flux = substrates_normalized

			# TODO -- we don't want a flux of 1... need a different rate law
			# return flux
			return 0.0

		return rate_law



