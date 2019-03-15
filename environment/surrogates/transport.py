from __future__ import absolute_import, division, print_function

import time
import numpy as np
import random

from agent.inner import CellSimulation


TUMBLE_JITTER = 0.4 # (radians)
# TODO: convert sample_fluxes to changes in molecules, round molecules to whole integers
# TODO:
def sample_fluxes(flux_distributions, media):
	'''
	Randomly sample a flux from the distribution of fluxes available for all reactions.
	Args:
		flux_distributions: Nested dictionary of the form {media: {reaction: [fluxes]}}
		media: One of three WCM environments (minimal, minimal_plus_amino_acids, minimal_minus_oxygen)

	Returns: dictionary of the form {reaction: flux}

	'''
	return {reaction: np.random.choice(distribution)
			for reaction, distribution in flux_distributions[media].iteritems()}

# TODO: only change fluxes when media changes


class Transport(CellSimulation):
	'''
	Simple transport surrogate that inherits chemotaxis behavior from environment.surrogates.chemotaxis.py
	'''

	def __init__(self, config):
		self.initial_time = 0.0
		self.local_time = 0.0
		# self.timestep = 1.0
		self.environment_change = {}
		self.volume = 1.0 # fL
		self.division_time = 100

		# initial state
		self.state = ['tumble']
		# self.external_concentrations = {
		# 	'GLC[p]': 0.0
		# }
		self.internal_counts = {
			'GLC': 1000,
			'GLT': 1000,
			'CYS': 1000
		} # TODO: initialize this with something from the config
		self.substrate_counts = config.get('substrate_counts')

		# load in fluxes
		self.flux_distributions = config.get('flux')
		self.stoichiometry = config.get('stoichiometry')

		# save these internal_concentrations so that they can be output later through a listener
		self.motile_force = [0.0, 0.0] # initial magnitude and relative orientation
		self.division = []

		# Initialize media environment and fluxes
		self.media = 'minimal'
		# self.fluxes = config.get('flux_distributions', {'minimal':{'transporter':[0,0,0,0,1]}})
		# self._molecule_ids = update['concentrations'].keys()


	def update_state(self):
		fluxes = sample_fluxes(self.flux_distributions, self.media) # TODO: make this once per cell cycle
		print(fluxes)
		# self.internal_counts
		# countsToMolar = 1 / (self.nAvogadro * self.volume)
		# delta_nutrients = ((1 / countsToMolar) * exchange_fluxes).astype(int)

	def update_behavior(self):
		# update behavior based on the current state of the system

		if self.state is 'run':
			force = 0.02
			torque = 0.0
			self.motile_force = [force, torque]
		elif self.state is 'tumble':
			force = 0.005
			torque = np.random.normal(scale=TUMBLE_JITTER)
			self.motile_force = [force, torque]

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

	def run_incremental(self, run_until):
		# update state once per message exchange
		self.update_state()
		self.update_behavior()
		# self.check_division()
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
