#!/usr/bin/env python

"""
InitiatedTranscripts

@author: Nick Ruggero
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 10/28/14
"""

from __future__ import division

import numpy as np
import tables

import wholecell.listeners.listener

VERBOSE = False

class InitiatedTranscripts(wholecell.listeners.listener.Listener):
	""" InitiatedTranscripts """

	_name = 'InitiatedTranscripts'

	# Constructor
	def __init__(self, *args, **kwargs):
		super(InitiatedTranscripts, self).__init__(*args, **kwargs)


	# Construct object graph
	def initialize(self, sim, kb):
		super(InitiatedTranscripts, self).initialize(sim, kb)

		# Computed, saved attributes
		self.ratioStableToToalInitalized = None

		# Attributes broadcast by the TranscriptInitation process
		self.mRnaInitalized = None
		self.tRnaInitalized = None
		self.rRnaInitalized = None

	# Allocate memory
	def allocate(self):
		super(InitiatedTranscripts, self).allocate()

		# Computed, saved attributes
		self.ratioStableToToalInitalized = np.nan

		# Attributes broadcast by the PolypeptideElongation process
		self.mRnaInitalized = np.nan
		self.tRnaInitalized = np.nan
		self.rRnaInitalized = np.nan

	def update(self):
		totalRna = self.tRnaInitalized + self.rRnaInitalized + self.mRnaInitalized
		stableRna = self.tRnaInitalized + self.rRnaInitalized
		self.ratioStableToToalInitalized = stableRna / totalRna

	def pytablesCreate(self, h5file, expectedRows):
		dtype = {
			"time": tables.Float64Col(),
			"timeStep": tables.Int64Col(),
			"ratioStableToToalInitalized": tables.Float64Col(),
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
		entry["ratioStableToToalInitalized"] = self.ratioStableToToalInitalized

		entry.append()

		table.flush()

	# TODO: load method
