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

		self.environment.nutrientData = self.environment.getNutrientData(raw_data)
		self.environment.condition = "basal"
		self.environment.nutrientsTimeSeriesLabel = "000000_basal"
