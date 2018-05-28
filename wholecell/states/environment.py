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

		self.nutrientsTimeSeries = sim_data.nutrientsTimeSeries
		self.nutrientData = sim_data.external_state.environment.nutrientData

		# current condition and timeseries
		self.nutrientsTimeSeriesLabel = sim_data.external_state.environment.nutrientsTimeSeriesLabel
		self.condition = self.nutrientsTimeSeries[self.nutrientsTimeSeriesLabel][0][1]


	def update(self):

		while len(self.nutrientsTimeSeries[self.nutrientsTimeSeriesLabel]) and \
			self.time() > self.nutrientsTimeSeries[self.nutrientsTimeSeriesLabel][0][0]:
				_, self.condition = self.nutrientsTimeSeries[self.nutrientsTimeSeriesLabel].popleft()
