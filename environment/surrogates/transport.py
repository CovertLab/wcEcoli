from __future__ import absolute_import, division, print_function

import os
import time

import csv
import numpy as np
from scipy import constants

from agent.inner import CellSimulation
import wholecell.utils.filepath as fp


ROOT_PATH = os.path.dirname(os.path.dirname(os.path.abspath('environment/boot.py')))

TUMBLE_JITTER = 0.4 # (radians)
N_AVOGADRO = constants.N_A # Avogadro's number

# def select_transport_lookup_table(media):
# 	# Select a transport lookup table based on the surrogate's media condition
# 	transport_lookup_table_path = \
# 		ROOT_PATH + '/environment/condition/tables/aa_transport_lookup_' + media + '.tsv'
# 	transport_lookup_table = pd.read_csv(transport_lookup_table_path, sep='\t', header = 0).to_dict()
# 	return transport_lookup_table

'''The following functions are called sequentially to parse TSV files that contain the lookup tables from the WCM
transport experiments. In order: 
1. determine_substrate_set() determines the set (i.e. unique values) of molecular substrates being modified during this 
experiment based on the lattice media condition.
2. determine_stoichiometry() determines the relationship between the substrates and their flux value for each
transport reaction.
3. determine_reaction_flux() samples a flux distribution from the lookup tables and picks a single flux value for each unique 
substrate for the duration of this surrogate's generation or until the media changes.
4. determine_substrate_flux() combines the output of the previous 3 functions. First, it determines the affect of 
each transport reaction on every substrate by multiplying each substrate's flux per reaction by its stoichiometry. 
Then, it combines all non-unique substrate values into one dictionary {substrate: total flux over all reactions}'''

def determine_substrate_set(media):
	# Determine the set of substrates that undergo flux in the selected lookup table
	# substrate_set = dict of form {substrate : 0}, with all keys being unique
	substrate_set = fp.read_json_file(
		ROOT_PATH + '/environment/condition/tables/aa_transport_lookup_' + media + '/substrate_set.json')
	return substrate_set

def determine_stoichiometry(media):
	# Determine the stoichiometry of each transport reaction
	# stoichiometry = dict of form {index : {substrate : values}} where index corresponds to reaction id
	stoichiometry = fp.read_json_file(
		ROOT_PATH + '/environment/condition/tables/aa_transport_lookup_' + media + '/stoichiometry.json')
	return stoichiometry

def determine_reaction_flux(media):
	# Select single flux value for each transport reaction based on sampling the selected lookup table
	# distribution = dict of form {index : [a, b, c, ...]} where index corresponds to reaction id & a/b/c... = fluxes
	# flux = single sampled value of distribution at each index
	distribution = fp.read_json_file(
		ROOT_PATH + '/environment/condition/tables/aa_transport_lookup_' + media + '/flux.json')
	flux = []
	for index in range(len(distribution)):
		flux.append(np.random.choice(distribution[str(index)]))
	return flux

def determine_substrate_flux(substrate_set, stoichiometry, flux):
	# Use randomly-selected fluxes for each reaction * stoichiometry of each reaction to determine total flux
	# substrate_set = dict of form {substrate : n}, with all keys being unique and all n's as sums of reaction_flux
	for index in range(len(stoichiometry)):
		for substrate in stoichiometry[str(index)]:
			value = stoichiometry[str(index)][substrate] * flux[index]
			if substrate in substrate_set:
				substrate_set[substrate] += value
				print("working on determine_substrate_flux()...")
	print("finished working on determine_substrate_flux()")
	return substrate_set

class Transport(CellSimulation):
	'''
	Surrogate that uses simple look-up tables to determine flux of various substrates into and out of itself.
	'''

	def __init__(self, config, agent_id):
		# give self a unique ID
		# self.random_id = np.random.randint(1, 1000000)
		self.listener_id = None

		self.initial_time = 0.0
		self.local_time = 0.0
		self.environment_change = {}
		self.volume = config.get('volume', 1.0)  # fL
		self.division_volume = 2.0

		self.timeline = config.get('timeline')
		self.media = config.get('media',
								ROOT_PATH + '/environment/condition/tables/aa_transport_lookup_minimal.tsv')

		# self.media_name = truncated version of the media name without filepath:
		beginning = self.media.find('lookup_')
		end = self.media.find('.tsv')
		self.media_name = self.media[beginning+7:end]
		print(self.media_name)

		self.substrate_set = {}
		self.stoichiometry = {}
		self.reaction_flux = {}
		self.substrate_flux = {}

		self.substrate_set = determine_substrate_set(self.media_name)
		self.stoichiometry = determine_stoichiometry(self.media_name)
		self.reaction_flux = determine_reaction_flux(self.media_name)
		self.substrate_flux = determine_substrate_flux(self.substrate_set, self.stoichiometry, self.reaction_flux)

		# get seed for random number generator, default 0
		self.seed = config.get('seed', 0)

		self.motile_force = [0.0, 0.0]  # initial magnitude and relative orientation
		self.division = []

		# TODO (Lee): get self.id to return agent_id, not a random string
		self.agent_id = agent_id

	def update_volume(self):
		'''
		Add volume based on total substrates present in the surrogate and write output to transport_experiment.JSON.

		'''
		# TODO (Lee): exchange this linear addition for a true growth function
		delta_volume = 0

		# TODO (Lee): check against 'environment_molecules.tsv' for actual list of internal and external molecules
		for substrate in self.substrate_set:
			delta_volume += float(self.substrate_set[substrate])
		delta_volume = delta_volume / self.volume
		self.volume += delta_volume

	def divide(self):
		return self.division

	def check_division(self):
		# if self.local_time >= self.initial_time + self.division_time:
		if self.volume >= self.division_volume:
			# halve internal substrate counts and pass them along
			daughter_substrate_counts = {}
			for substrate in self.substrate_set:
				daughter_substrate_counts[substrate] = self.substrate_set[substrate] * 0.5

			# define a common dictionary to pass to daughters
			common = dict(
				time=self.local_time,
				substrates = daughter_substrate_counts,
				volume=self.volume * 0.5,
				media=self.media)
			self.division = [
				dict(common,
					 seed=37 * self.seed + 47 * 1 + 997),
				dict(common,
					 seed=37 * self.seed + 47 * 2 + 997)]

		return self.division


	def time(self):
		return self.local_time

	def transport_listener(self):
		listener_topics = [self.agent_id, self.local_time, self.volume, self.media_name]
		with open(ROOT_PATH + '/out/manual/transport_listener_' + str(self.listener_id) + '.csv', mode="a") as f:
			writer = csv.writer(f)
			writer.writerow(listener_topics)

	def apply_outer_update(self, update):
		# Update media conditions based on function select_media in lattice.py:
		if self.media != update['media']:
			self.media = update['media']
			self.substrate_set = determine_substrate_set(self.media_name)
			self.stoichiometry = determine_stoichiometry(self.media_name)
			self.reaction_flux = determine_reaction_flux(self.media_name)
			self.substrate_flux = determine_substrate_flux(self.substrate_set, self.stoichiometry, self.reaction_flux)
		self.external_concentrations = update['concentrations']
		self.environment_change = self.substrate_set
		self.listener_id = update['listener_id']
		self.transport_listener()

	def run_incremental(self, run_until):
		'''
		Checks the argument run_until, and then updates internal state until run_until is reached
		Args:
			run_until: final timepoint of the experiment
		'''
		self.update_volume()
		self.local_time = run_until
		self.check_division()

	def generate_inner_update(self):
		return {
			'agent id' : self.agent_id,
			'volume': self.volume,
			'motile_force': self.motile_force,
			'environment_change': self.environment_change,
			'division': self.division,
			'media': self.media,
			}

	def synchronize_state(self, state):
		if 'time' in state:
			self.initial_time = state['time']
