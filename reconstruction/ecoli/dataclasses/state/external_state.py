"""
Simulation data for external state

This base class includes all data associated with states external to the cells.
Initializes the environment using conditions and time series from raw_data.

	- environment.nutrients_time_series: a dictionary of all time series.
	- environment.nutrients_time_series_label: a string specifying the time series
		used for the current simulation.
	- environment.nutrients: a dictionary of environmental nutrients (keys) and
		their concentrations (values).
	- environment.environment_dict: a dictionary of all environments, each one
		itself a dictionary nutrients (keys) and their concentrations (values).
	- environment.nutrients_defs: a list of dicts,  with keys 'id' and 'import
		constrained'. These are the nutrient's molecule ID and a boolean indicating
		whether the nutrient is import constrained.
	- environment.secretions: a list of dicts, with keys 'molecule id', 'lower
		bound' and 'upper bound' for the flux constraints.

Notes:
-----
	- secretions includes non-growth associate maintenance nutrients, which
		are not used in metabolism.

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

		#TODO (Eran) is Environment(raw_data, sim_data) needed?
		self.environment = Environment(raw_data, sim_data)

		# default parameters
		# self.environment.volume = "infinite"
		self.environment.nutrients_time_series_label = "000000_basal"

		# create a dictionary of all nutrient time series
		self.environment.nutrients_time_series = {}
		for label in dir(raw_data.condition.timeseries):
			if label.startswith("__"):
				continue
			self.environment.nutrients_time_series[label] = []
			timeseries = getattr(raw_data.condition.timeseries, label)
			for row in timeseries:
				self.environment.nutrients_time_series[label].append((
					row["time"].asNumber(units.s),
					row["nutrients"].encode("utf-8"),
					row["volume_liters"].encode("utf-8")
					))

		# create a dictionary with all saved environments
		self.environment.environment_dict = {}
		for label in dir(raw_data.condition.environment):
			if label.startswith("__"):
				continue
			self.environment.environment_dict[label] = {}
			environments = getattr(raw_data.condition.environment, label)
			for row in environments:
				self.environment.environment_dict[label].update({row["moleculeid"]: row["concentration"]})

		# initial state
		initial_environment = "minimal"
		self.environment.nutrients = self.environment.environment_dict[initial_environment]
