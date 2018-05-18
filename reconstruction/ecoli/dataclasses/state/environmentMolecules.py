"""
SimulationData for environment molecules state

@author: Eran Agmon
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 05/17/2018
"""

from __future__ import division

import numpy as np

from wholecell.utils import units
from wholecell.utils.unit_struct_array import UnitStructArray

from reconstruction.ecoli.dataclasses.state.stateFunctions import addToStateEnvironment

class EnvironmentMolecules(object):
	""" EnvironmentMolecules """

	def __init__(self, raw_data, sim_data):
		environmentData = np.zeros(
			0,
			dtype = [
				("id", "a50"),
				# TODO (Eran) -- mass might not be needed here
				("mass", "{}f8".format(len(sim_data.molecular_weight_order))),
				]
			)

		# Add units to values
		field_units = {
			"id"		:	None,
			"mass"				:	units.g / units.mol,
			}

		self.environmentData = UnitStructArray(environmentData, field_units)

	def addToEnvironmentState(self, ids, masses):
		self.environmentData = addToStateEnvironment(self.environmentData, ids, masses)