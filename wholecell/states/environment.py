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

	def initialize(self, sim, sim_data):
		super(Environment, self).initialize(sim, sim_data)

		self.nutrients_time_series = sim_data.nutrientsTimeSeries
		self.nutrient_data = sim_data.external_state.environment.nutrient_data

		self.nutrients_time_series_label = sim_data.external_state.environment.nutrients_time_series_label
		self.condition = self.nutrients_time_series[self.nutrients_time_series_label][0][1]

		self.time_series = [t[0] for t in self.nutrients_time_series[self.nutrients_time_series_label]]


	def update(self):
		current_index = [i for i, t in enumerate(self.time_series) if self.time()>=t][-1]
		self.condition = self.nutrients_time_series[self.nutrients_time_series_label][current_index][1]
