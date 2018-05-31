#!/usr/bin/env python

"""
External state that represents environmental molecules and conditions.

- nutrient_data a dictionary containing:
		- externalExchangeMolecules
		- importExchangeMolecules
		- importConstrainedExchangeMolecules
		- importUnconstrainedExchangeMolecules
		- secretionExchangeMolecules

- condition: a string specifying the current condition.

- nutrient_time_series: a list of tuples that include time and condition in which shifts occur.

- time_series: a list of all times at which the condition changes.

@organization: Covert Lab, Department of Bioengineering, Stanford University
"""

import numpy as np

import wholecell.states.external_state


class Environment(wholecell.states.external_state.ExternalState):
	_name = 'Environment'

	def initialize(self, sim, sim_data):
		super(Environment, self).initialize(sim, sim_data)

		self.nutrient_data = sim_data.external_state.environment.nutrient_data
		self.nutrient_time_series = sim_data.external_state.environment.nutrients_time_series[
			sim_data.external_state.environment.nutrients_time_series_label
			]
		self.condition = self.nutrient_time_series[0][1]
		self.time_series = [t[0] for t in self.nutrient_time_series]


	def update(self):
		current_index = [i for i, t in enumerate(self.time_series) if self.time()>=t][-1]
		self.condition = self.nutrient_time_series[current_index][1]