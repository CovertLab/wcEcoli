#!/usr/bin/env python

"""
GrowthRateControl

@author: Nick Ruggero
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 1/7/14
"""

from __future__ import division

import numpy as np
import tables

import wholecell.listeners.listener

VERBOSE = False

class GrowthRateControl(wholecell.listeners.listener.Listener):
	""" GrowthRateControl """

	_name = 'GrowthRateControl'

	# Constructor
	def __init__(self, *args, **kwargs):
		super(GrowthRateControl, self).__init__(*args, **kwargs)


	# Construct object graph
	def initialize(self, sim, kb):
		super(GrowthRateControl, self).initialize(sim, kb)

		# Computed, saved attributes
		self.ratioStableToToalInitalized = None

		# Attributes broadcast by the TranscriptInitation process
		self.totalStalls = None
		self.synthetaseSaturation = None
		#self.rRnaInitalized = None

	# Allocate memory
	def allocate(self):
		super(GrowthRateControl, self).allocate()

		# Computed, saved attributes
		#self.ratioStableToToalInitalized = np.nan

		# Attributes broadcast by the PolypeptideElongation process
		self.totalStalls = np.nan
		self.synthetaseSaturation = np.zeros(21, np.int64)
		#self.rRnaInitalized = np.nan

	def update(self):
		# totalRna = self.tRnaInitalized + self.rRnaInitalized + self.mRnaInitalized
		# stableRna = self.tRnaInitalized + self.rRnaInitalized
		# if totalRna > 0:
		# 	self.ratioStableToToalInitalized = stableRna / totalRna
		# else:
		# 	self.ratioStableToToalInitalized = np.nan
		pass

	def pytablesCreate(self, h5file, expectedRows):
		dtype = {
			"time": tables.Float64Col(),
			"timeStep": tables.Int64Col(),
			"totalStalls": tables.Float64Col(),
			"synthetaseSaturation": tables.Float64Col(self.synthetaseSaturation.size)
			}

		table = h5file.create_table(
			h5file.root,
			self._name,
			dtype,
			title = self._name,
			filters = tables.Filters(complevel = 9, complib="zlib"),
			)


	def pytablesAppend(self, h5file):
		table = h5file.get_node("/", self._name)

		entry = table.row

		entry["time"] = self.time()
		entry["timeStep"] = self.timeStep()
		entry["totalStalls"] = self.ratioStableToToalInitalized
		entry["synthetaseSaturation"] = self.synthetaseSaturation

		entry.append()

		table.flush()

	# TODO: load method
