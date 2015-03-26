#!/usr/bin/env python

"""
InitiatedTranscripts

@author: Nick Ruggero
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 10/28/14
"""

from __future__ import division

import numpy as np

import wholecell.listeners.listener

VERBOSE = False

class InitiatedTranscripts(wholecell.listeners.listener.Listener):
	""" InitiatedTranscripts """

	_name = 'InitiatedTranscripts'

	# Constructor
	def __init__(self, *args, **kwargs):
		super(InitiatedTranscripts, self).__init__(*args, **kwargs)

		self.ratioUnits = "stable/total"

	# Construct object graph
	def initialize(self, sim, kb):
		super(InitiatedTranscripts, self).initialize(sim, kb)

	# Allocate memory
	def allocate(self):
		super(InitiatedTranscripts, self).allocate()

		# Computed, saved attributes
		self.ratioStableToToalInitalized = 0.

		# Attributes broadcast by the PolypeptideElongation process
		self.mRnaInitalized = 0
		self.tRnaInitalized = 0
		self.rRnaInitalized = 0

	def update(self):
		totalRna = self.tRnaInitalized + self.rRnaInitalized + self.mRnaInitalized
		stableRna = self.tRnaInitalized + self.rRnaInitalized
		if totalRna > 0:
			self.ratioStableToToalInitalized = stableRna / totalRna
		else:
			self.ratioStableToToalInitalized = 0.

	def tableCreate(self, tableWriter):
		# Store units as metadata
		tableWriter.writeAttributes(
			ratioStableToToalInitalized = self.ratioUnits,
			)

	def tableAppend(self, tableWriter):
		tableWriter.append(
			time = self.time(),
			timeStep = self.timeStep(),
			ratioStableToToalInitalized = self.ratioStableToToalInitalized
			)

