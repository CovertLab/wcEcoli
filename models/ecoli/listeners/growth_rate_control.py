#!/usr/bin/env python

"""
GrowthRateControl

@author: Nick Ruggero
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 1/7/14
"""

from __future__ import division

import numpy as np

import wholecell.listeners.listener

VERBOSE = False

class GrowthRateControl(wholecell.listeners.listener.Listener):
	""" GrowthRateControl """

	_name = 'GrowthRateControl'

	# Constructor
	def __init__(self, *args, **kwargs):
		super(GrowthRateControl, self).__init__(*args, **kwargs)
		self.countUnits = "counts"
		self.saturationUnits = "fraction saturated"

	# Construct object graph
	def initialize(self, sim, kb):
		super(GrowthRateControl, self).initialize(sim, kb)

	# Allocate memory
	def allocate(self):
		super(GrowthRateControl, self).allocate()

		# Computed, saved attributes

		# Attributes broadcast by the PolypeptideElongation process
		self.totalStalls = 0
		self.synthetaseSaturation = np.zeros(21, np.float64)
		self.spoT_saturation = 0.

	def tableCreate(self, tableWriter):
		# Store units as metadata
		tableWriter.writeAttributes(
			totalStalls = self.countUnits,
			synthetaseSaturation = self.saturationUnits,
			spoT_saturation = self.saturationUnits
			)

	def tableAppend(self, tableWriter):
		tableWriter.append(
			time = self.time(),
			timeStep = self.timeStep(),
			totalStalls = self.totalStalls,
			synthetaseSaturation = self.synthetaseSaturation,
			spoT_saturation = self.spoT_saturation
			)
