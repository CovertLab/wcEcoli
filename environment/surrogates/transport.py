from __future__ import absolute_import, division, print_function

import ast
import time
import csv

import pandas as pd
import numpy as np
from scipy import constants

from agent.inner import CellSimulation
import os
import environment

ROOT_PATH = os.path.dirname(os.path.dirname(os.path.abspath(environment.__file__)))

TUMBLE_JITTER = 0.4 # (radians)
N_AVOGADRO = constants.N_A # Avogadro's number

def select_transport_table(media):
	transport_table = pd.read_csv(media, sep='\t', header=0).to_dict()
	return transport_table

# Generate a SUBSTRATES dictionary with substrate keys padded with zeros
def make_substrates(transport_table):
	SUBSTRATES = []
	for index in transport_table['stoichiometry']:
		keys = ast.literal_eval(transport_table['stoichiometry'][index]).keys()
		i = 0
		while i < len(keys):
			if keys[i] not in SUBSTRATES:
				SUBSTRATES.append(keys[i])
			i += 1
	zero_pad = [0] * len(SUBSTRATES)
	SUBSTRATES = dict(zip(SUBSTRATES, zero_pad))
	return SUBSTRATES

# Generate a FLUX list with single flux values for each reaction sampled from a distribution
def sample_fluxes(transport_table):
	i = 0
	FLUX = []
	while i < len(transport_table['reaction id']):
		FLUX.append(np.random.choice(ast.literal_eval(transport_table['distribution'][i])))
		i += 1
	return FLUX

# Modify SUBSTRATES dictionary
def transport(transport_table, SUBSTRATES, FLUX):
	len(transport_table['reaction id'])
	i = 0
	while i < len(transport_table['reaction id']):
		stoich = ast.literal_eval(transport_table['stoichiometry'][i])
		for substrate in stoich:
			value = stoich[substrate] * FLUX[i]
			if substrate in SUBSTRATES:
				SUBSTRATES[substrate] += value
		i += 1
	return SUBSTRATES


class Transport(CellSimulation):
	'''
	Surrogate that uses simple look-up tables to determine flux of various substrates into and out of itself.
	'''

	def __init__(self, config):
		# give self a unique ID
		# self.random_id = np.random.randint(1, 1000000)

		self.initial_time = 0.0
		self.local_time = 0.0
		self.environment_change = {}
		self.volume = config.get('volume', 1.0)  # fL
		self.division_volume = 2.0

		# Get media condition and select associated look-up table
		# TODO: Eliminate dependency on new media/timeline files
		self.timeline = config.get('timeline')
		self.media = config.get('media', ROOT_PATH + '/environment/condition/tables/aa_transport_lookup_amino_acids.tsv')
		self.transport_table = select_transport_table(self.media)

		# Build substrate list from look-up table
		self.substrates = config.get('substrates', make_substrates(self.transport_table))
		# Sample flux distributions from a distribution
		self.flux = sample_fluxes(self.transport_table)
		# get seed for random number generator, default 0
		self.seed = config.get('seed', 0)

		self.motile_force = [0.0, 0.0]  # initial magnitude and relative orientation
		self.division = []

	def update_state(self):
		'''
		Update the internal and external substrates of the surrogate. Add volume based on total substrates present in
		the surrogate and write output to test_data.csv

		'''
		self.substrates = transport(self.transport_table, self.substrates, self.flux)

		# TODO: exchange this arbitrary addition for a growth function
		delta_volume = 0

		# TODO: check against 'environment_molecules.tsv' for actual list of internal and external molecules
		for substrate in self.substrates:
			delta_volume += float(self.substrates[substrate])
		delta_volume = delta_volume / self.volume
		self.volume += delta_volume

		# write output to 'test_data.csv'
		data_row = [self.seed, self.local_time, self.volume, self.media,
					float(self.substrates['GLT[c]']) / (self.volume * N_AVOGADRO)]
		with open(ROOT_PATH + '/environment/analysis/test_data.csv', 'a') as f:
			writer = csv.writer(f)
			writer.writerow(data_row)
			f.close()

	def divide(self):
		return self.division

	def check_division(self):
		# if self.local_time >= self.initial_time + self.division_time:
		if self.volume >= self.division_volume:
			# halve internal substrate counts and pass them along
			daughter_substrate_counts = {}
			for substrate in self.substrates:
				daughter_substrate_counts[substrate] = self.substrates[substrate] * 0.5

			# define a common dictionary to pass to daughters
			common = dict(
				time=self.local_time,
				substrates = daughter_substrate_counts,
				volume=self.volume * 0.5,
				media = self.media)
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
		self.external_concentrations = update['concentrations']
		self.environment_change = self.substrates

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
