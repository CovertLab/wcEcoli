from __future__ import absolute_import, division, print_function

import csv
import time

import numpy as np
from scipy import constants

from agent.inner import CellSimulation

TUMBLE_JITTER = 0.4 # (radians)
N_AVOGADRO = constants.N_A # Avogadro's number

def sample_fluxes(flux_distributions, media):
	'''
	Randomly sample a flux from the distribution of fluxes available for all reactions. If sample_fluxes is called once
	per timestep, a timeseries of internal_substrate_counts will appear erratic and random.
	Args:
		flux_distributions: Nested dictionary of the form {media: {reaction: [fluxes]}} from self.flux_distributions
		media: One of the WCM environments (e.g. minimal, minimal_plus_amino_acids, minimal_minus_oxygen)

	Returns: dictionary of the form {reaction: flux}

	'''
	fluxes = {}
	for reaction, distribution in flux_distributions[media].iteritems():
		fluxes[reaction] = np.random.choice(distribution)
	return fluxes


class Transport(CellSimulation):
	'''
	Surrogate that uses simple look-up tables to determine flux of various substrates into and out of itself.
	'''

	def __init__(self, config):
		# give self a unique ID
		self.random_id = np.random.randint(1, 10000)

		self.initial_time = 0.0
		self.local_time = 0.0
		self.environment_change = {}
		self.volume = config.get('volume', 1.0)  # fL
		self.division_time = 30

		# Get media condition from select_media function in lattice.py
		self.media = config.get('media', None)
		# get flux_distributions from make_flux_distributions() in control.py
		# note: self.flux_distributions is an entire list of choices, while self.fluxes is a single flux choice
		self.flux_distributions = config.get('flux')
		# get stoichiometry from make_stoichiometry() in control.py
		self.stoichiometry = config.get('stoichiometry')
		# get seed for random number generator, default 1
		self.seed = config.get('seed', 0)

		# Initialize internal and external substrate counts at t = 0
		self.internal_substrate_counts = config.get('internal_substrate_counts')
		self.external_substrate_counts = config.get('external_substrate_counts')

		self.motile_force = [0.0, 0.0]  # initial magnitude and relative orientation
		self.division = []

	def update_state(self):
		'''
		Update the internal state of the surrogate, including transport fluxes and internal/external substrate counts as delta_nutrients. Add volume based on total substrates present in the surrogate and write output to test_data.csv

		'''
		# fluxes = sample_fluxes(self.flux_distributions, self.media)

		# Set up delta_nutrients and modify internal_substrate_counts and external_substrate_counts
		delta_nutrients = {}
		delta_volume = 0
		for reaction in self.fluxes:
			for substrate in self.stoichiometry[reaction]:
				delta_nutrients[substrate] = int(
					(self.stoichiometry[reaction][substrate] * self.fluxes[reaction]) * (N_AVOGADRO * self.volume))
				delta_volume += self.stoichiometry[reaction][substrate] * 0.01

		# Printout for debugging
		print('*************** START update_state')
		print('media = ')
		print(self.media)
		print('fluxes = ')
		print(self.fluxes.items())
		print('stoichiometry = ')
		print(self.stoichiometry.items())
		print('delta_volume = ')
		print('volume = ')
		print(self.volume)
		print(delta_volume)
		print('delta_nutrients = ')
		print(delta_nutrients.items())
		print('internal_substrate_counts = ')
		print(self.internal_substrate_counts.items())
		print('internal_substrate_counts = ')
		print(self.external_substrate_counts.items())
		print('*************** END update_state')

		for substrate in self.internal_substrate_counts:
			self.internal_substrate_counts[substrate] += delta_nutrients[substrate]

		for substrate in self.external_substrate_counts:
			self.external_substrate_counts[substrate] += delta_nutrients[substrate]

		# TODO: exchange this arbitrary addition for a growth function
		self.volume += delta_volume

		# write output to 'test_data.csv'
		data_row = [self.random_id, self.local_time, self.volume, self.media,
					self.internal_substrate_counts['GLC[c]'] / (self.volume * N_AVOGADRO),
					self.internal_substrate_counts['CYS[c]'] / (self.volume * N_AVOGADRO),
					self.internal_substrate_counts['CYS[c]'] / (self.volume * N_AVOGADRO)]
		with open('test_data.csv', 'a') as f:
			writer = csv.writer(f)
			writer.writerow(data_row)
			f.close()

	def divide(self):
		return self.division

	def check_division(self):
		# update division state based on time since initialization

		if self.local_time >= self.initial_time + self.division_time:

			# halve internal substrate counts and pass them along
			daughter_substrate_counts = {}
			for substrate in self.internal_substrate_counts:
				daughter_substrate_counts[substrate] = self.internal_substrate_counts[substrate] * 0.5

			# define a common dictionary to pass to daughters
			common = dict(
				time=self.local_time,
				flux=self.flux_distributions,
				stoichiometry = self.stoichiometry,
				internal_substrate_counts=daughter_substrate_counts,
				external_substrate_counts=self.external_substrate_counts,
				volume=self.volume * 0.5)
			self.division = [
				dict(common,
					 seed=37 * self.seed + 47 * 1 + 997),
				dict(common,
					 seed=37 * self.seed + 47 * 2 + 997)]

		return self.division

	def time(self):
		return self.local_time

	def apply_outer_update(self, update):
		# Update media conditions based on function select_media in lattice.py:
		if self.media != update['media']:
				self.media = update['media']
				self.fluxes = sample_fluxes(self.flux_distributions, self.media)
		self.external_concentrations = update['concentrations']
		self.environment_change = {}

		# TODO: generalize this ASAP:
		# quick and dirty environmental modification, this will change ASAP when molecules are renamed and lookup tables are implemented. For now, this is for demonstration purposes only:
		self.environment_change['GLC[p]'] = self.external_substrate_counts['GLC[p]']
		self.environment_change['CYS[p]'] = self.external_substrate_counts['CYS[p]']
		self.environment_change['GLT[p]'] = self.external_substrate_counts['GLT[p]']

	def run_incremental(self, run_until):
		'''
		Checks the argument run_until, and then updates internal state until run_until is reached
		Args:
			run_until: final timepoint of the experiment
		'''
		self.update_state()
		self.local_time = run_until
		self.check_division()

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
