from __future__ import absolute_import, division, print_function

import time
import numpy as np

from agent.inner import CellSimulation


TUMBLE_JITTER = 2.0 # (radians)
DEFAULT_COLOR = [76, 0 , 153]
DEFAULT_COLOR = [color/255 for color in DEFAULT_COLOR]

class Endocrine(CellSimulation):
	''' Endocrine Surrogate '''

	def __init__(self):
		self.initial_time = 0.0
		self.local_time = 0.0
		# self.timestep = 1.0
		self.environment_change = {}
		self.volume = 1.0
		self.division_time = 100

		# initial state
		# self.state = ['tumble']
		self.external_concentrations = {
			'GLC': 0.0
		}
		self.internal_concentrations = {
			'A': 5.0,
		}

		self.motile_force = [0.0, 0.0] # initial magnitude and relative orientation
		self.division = []


	def update_state(self):
		# update state based on internal and external concentrations
		self.activation = self.external_concentrations['GLC']

		# update intracellular concentrations
		self.internal_concentrations['CheY-P'] = self.external_concentrations['GLC']


	def update_behavior(self):
		pass

	def check_division(self):
		# update division state based on time since initialization
		if self.local_time >= self.initial_time + self.division_time:
			self.division = [{'time': self.local_time}, {'time': self.local_time}]

		return self.division

	def time(self):
		return self.local_time

	def apply_outer_update(self, update):
		self.external_concentrations = update['concentrations']

		self.environment_change = {}
		for molecule in self.external_concentrations.iterkeys():
			self.environment_change[molecule] = 0

	def run_incremental(self, run_until):
		# update state once per message exchange
		self.update_state()
		# self.update_behavior()
		# self.check_division()
		self.local_time = run_until

		time.sleep(0.2)  # pause for better coordination with Lens visualization. TODO: remove this

	def generate_inner_update(self):
		return {
			'volume': self.volume,
			'motile_force': self.motile_force,
			'environment_change': self.environment_change,
			'division': self.division,
			'color': DEFAULT_COLOR,
			}

	def synchronize_state(self, state):
		if 'time' in state:
			self.initial_time = state['time']
