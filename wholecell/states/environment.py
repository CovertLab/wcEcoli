#!/usr/bin/env python

"""
Environment.py

External State which represents the class of environmental molecules.

@organization: Covert Lab, Department of Bioengineering, Stanford University
"""

import numpy as np

import wholecell.states.external_state


class Environment(wholecell.states.external_state.ExternalState):
	_name = 'Environment'

	def __init__(self, *args, **kwargs):

		super(Environment, self).__init__(*args, **kwargs)


	def initialize(self, sim, sim_data):
		super(Environment, self).initialize(sim, sim_data)

		self.nutrients_time_series = sim_data.nutrientsTimeSeries
		self.nutrient_data = sim_data.external_state.environment.nutrient_data

		# current condition and timeseries
		self.nutrients_time_series_label = sim_data.external_state.environment.nutrients_time_series_label
		self.condition = self.nutrients_time_series[self.nutrients_time_series_label][0][1]


	def update(self):

		while len(self.nutrients_time_series[self.nutrients_time_series_label]) and \
			self.time() > self.nutrients_time_series[self.nutrients_time_series_label][0][0]:
				_, self.condition = self.nutrients_time_series[self.nutrients_time_series_label].popleft()
