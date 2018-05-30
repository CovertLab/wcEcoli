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

		self.conditions = {}
		for row in raw_data.condition.condition_defs:
			condition = row["condition"].encode("utf-8")
			self.conditions[condition] = {}
			self.conditions[condition]["nutrients"] = row["nutrients"].encode("utf-8")
			self.conditions[condition]["perturbations"] = row["genotype perturbations"]
