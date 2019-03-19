from __future__ import absolute_import, division, print_function

import time
import numpy as np
from scipy import constants
import csv

from agent.inner import CellSimulation


TUMBLE_JITTER = 0.4 # (radians)

def sample_fluxes(flux_distributions, media):
	'''
	Randomly sample a flux from the distribution of fluxes available for all reactions.
	Args:
		flux_distributions: Nested dictionary of the form {media: {reaction: [fluxes]}} from self.flux_distributions
		media: One of three WCM environments (minimal, minimal_plus_amino_acids, minimal_minus_oxygen)

	Returns: dictionary of the form {reaction: flux}

	'''
	return {reaction: np.random.choice(distribution)
			for reaction, distribution in flux_distributions[media].iteritems()}

# TODO: only change fluxes when media changes

N_AVOGADRO = constants.N_A

class Transport(CellSimulation):
	'''
	Surrogate that uses simple look-up tables to determine flux of various substrates into and out of itself.
	'''

	def __init__(self, config):

		# give self a unique ID
		self.random_id = np.random.randint(1,10000)

		self.initial_time = 0.0
		self.local_time = 0.0
		self.environment_change = {}
		self.volume = 1.0 # fL
		self.division_time = 100

		# Initial state
		# self.internal_counts = {
		# 	'GLC': 1000,
		# 	'GLT': 1000,
		# 	'CYS': 1000
		# }


		self.delta_nutrients = 0

		# Initialize media environment and flux distributions/associated stoichiometry
		# self.media = 'minimal'
		self.media = config.get('media')
		self.flux_distributions = config.get('flux')
		self.stoichiometry = config.get('stoichiometry')

		# Random substrate counts
		# self.internal_substrate_counts = config.get('substrate_counts')
		self.internal_substrate_counts = {'GLC[c]': np.random.randint(1000,100000), 'GLT[c]': np.random.randint(1000,100000), 'CYS[c]': np.random.randint(1000,100000)}
		self.external_substrate_counts = {'GLC[p]': 0, 'GLT[p]': 0, 'CYS[p]': 0}


		self.motile_force = [0.0, 0.0] # initial magnitude and relative orientation
		self.division = []

	def update_state(self):
		'''
		Update the internal state of the surrogate, including transport fluxes and internal substrate counts.

		TODO: modify commented code below to convert environmental substrate concentrations to counts
		TODO: prepare a message about how many molecules to change in the environment to pass in self.environment_change
		Returns: Nothing, except a printout of fluxes

		'''
		fluxes = sample_fluxes(self.flux_distributions, self.media) # TODO: make this once per cell cycle instead of once per time step, move to init
		print("fluxes = ")
		print(fluxes)

		# Set up delta_nutrients and modify internal_substrate_counts
		self.delta_nutrients = {}
		delta_volume = 0
		print(N_AVOGADRO)
		for reaction in fluxes:
			for substrate in self.stoichiometry[reaction]:
				self.delta_nutrients[substrate] = int((self.stoichiometry[reaction][substrate] * fluxes[reaction]) * (N_AVOGADRO * self.volume))
				delta_volume += self.stoichiometry[reaction][substrate] * 0.1

		for substrate in self.internal_substrate_counts:
			self.internal_substrate_counts[substrate] += self.delta_nutrients[substrate]
			print("internal substrates: " + substrate)
			print(self.internal_substrate_counts[substrate])

		for substrate in self.external_substrate_counts:
			self.external_substrate_counts[substrate] += self.delta_nutrients[substrate]
			print("external substrates: " + substrate)
			print(self.external_substrate_counts[substrate])

		# TODO: partition delta_nutrients into internal and external molecules: internal passes to self.internal_counts and external passes to self.environment_change

		# add arbitrary volume
		self.volume += delta_volume # TODO: relate volume change to delta_nutrients and self.internal_counts

		# write output to 'test_data.csv'
		data_row = [self.random_id, self.local_time, self.volume, self.media,
			self.internal_substrate_counts['GLC[c]']/(self.volume* N_AVOGADRO), self.internal_substrate_counts['CYS[c]']/(self.volume* N_AVOGADRO), self.internal_substrate_counts['CYS[c]']/(self.volume* N_AVOGADRO)]

		# with open('test_data.csv', 'ab') as fd:
		# 	fd.write(data_row)
		# 	fd.close()

		with open('test_data.csv', 'a') as f:
			writer = csv.writer(f)
			writer.writerow(data_row)
			f.close()


	def check_division(self):
		# update division state based on time since initialization

		if self.local_time >= self.initial_time + self.division_time:
			self.division = [{'time': self.local_time}, {'time': self.local_time}]

		return self.division

	def time(self):
		return self.local_time

	def apply_outer_update(self, update):
		# Update media conditions based on function select_media in lattice.py:
		self.media = update['media']
		self.external_concentrations = update['concentrations']
		self.environment_change = {}
		# quick and dirty environmental modification, this will change ASAP when molecules are renamed and lookup tables are implemented. For now, this is for demonstration purposes only:
		self.environment_change['GLC[p]'] = self.external_substrate_counts['GLC[p]']
		self.environment_change['CYS[p]'] = self.external_substrate_counts['CYS[p]']
		self.environment_change['GLT[p]'] = self.external_substrate_counts['GLT[p]']
		print('environment_change dict: ')
		print(self.environment_change)

	def run_incremental(self, run_until):
		'''
		Checks the argument run_until, and then updates internal state until run_until is reached
		Args:
			run_until: final timepoint of the experiment
		'''
		self.update_state()
		self.local_time = run_until

		time.sleep(1.0)  # pause for better coordination with Lens visualization. TODO: remove this

	def generate_inner_update(self):
		return {
			'volume': self.volume,
			'motile_force': self.motile_force,
			'environment_change': self.environment_change,
			'division': self.division,
			'media': self.media,
			}

	def synchronize_state(self, state):
		if 'time' in state:
			self.initial_time = state['time']
