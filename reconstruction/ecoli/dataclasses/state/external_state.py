"""
SimulationData state associated data

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

		self._buildEnvironment(raw_data, sim_data)


	def _buildEnvironment(self, raw_data, sim_data):

		self.environment.nutrient_data = self.environment.getNutrientData(raw_data)
		self.environment.condition = "basal"
		self.environment.nutrients_time_series_label = "000000_basal"

		# TODO (ERAN) remove gene perturbations from here, keep in sim_data
		self.environment.conditions = {}
		for row in raw_data.condition.condition_defs:
			condition = row["condition"].encode("utf-8")
			self.environment.conditions[condition] = {}
			self.environment.conditions[condition]["nutrients"] = row["nutrients"].encode("utf-8")
			self.environment.conditions[condition]["perturbations"] = row["genotype perturbations"]

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