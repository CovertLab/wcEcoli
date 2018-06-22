"""
Simulation data for external state

This base class includes all data associated with states external to the cells.

	- environment.nutrients_time_series: a dictionary of all time series.

	- environment.nutrients_time_series_label: a string specifying the time series
		used for the current simulation.


@organization: Covert Lab, Department of Bioengineering, Stanford University
"""

from __future__ import division

from wholecell.utils import units
from wholecell.utils.unit_struct_array import UnitStructArray

from reconstruction.ecoli.dataclasses.state.environment import Environment

from reconstruction.ecoli.dataclasses.state import stateFunctions as sf

import re
import numpy as np

class ExternalState(object):
	""" External State """

	def __init__(self, raw_data, sim_data):

		self.environment = Environment(raw_data, sim_data)

		self.environment.nutrients_time_series_label = "000000_basal"

		self.environment.nutrients_time_series = {}
		for label in dir(raw_data.condition.timeseries):
			if label.startswith("__"):
				continue

			self.environment.nutrients_time_series[label] = []
			timeseries = getattr(raw_data.condition.timeseries, label)
			for row in timeseries:
				self.environment.nutrients_time_series[label].append((
					row["time"].asNumber(units.s),
					row["nutrients"].encode("utf-8")
					))

