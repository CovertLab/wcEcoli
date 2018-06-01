#!/usr/bin/env python

"""
External state that represents environmental molecules and conditions.

- nutrient_data a dictionary containing keys:
     - externalExchangeMolecules
     - importExchangeMolecules
     - importConstrainedExchangeMolecules
     - importUnconstrainedExchangeMolecules
     - secretionExchangeMolecules

- nutrients_time_series: a list of tuples that include time and nutrients in which shifts occur.

- nutrients: a string specifying the current nutrients.

- times: a list of all times at which the nutrients shift.

@organization: Covert Lab, Department of Bioengineering, Stanford University
"""

import numpy as np

import wholecell.states.external_state


class Environment(wholecell.states.external_state.ExternalState):
	_name = 'Environment'

	def initialize(self, sim, sim_data):
		super(Environment, self).initialize(sim, sim_data)

		self.nutrient_data = sim_data.external_state.environment.nutrient_data
		self.nutrients_time_series = sim_data.external_state.environment.nutrients_time_series[
			sim_data.external_state.environment.nutrients_time_series_label
			]
		self.nutrients = self.nutrients_time_series[0][1]
		self.times = [t[0] for t in self.nutrients_time_series]


	def update(self):
		current_index = [i for i, t in enumerate(self.times) if self.time()>=t][-1]
		self.nutrients = self.nutrients_time_series[current_index][1]
