#!/usr/bin/env python

"""
ControlLoop

@author: Nick Ruggero
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 2/17/17
"""

from __future__ import division

import numpy as np

import wholecell.listeners.listener

# from numpy.lib.recfunctions import merge_arrays

VERBOSE = False

class ControlLoop(wholecell.listeners.listener.Listener):
	""" ControlLoop """

	_name = 'ControlLoop'

	# Constructor
	def __init__(self, *args, **kwargs):
		super(ControlLoop, self).__init__(*args, **kwargs)


	# Construct object graph
	def initialize(self, sim, sim_data):
		super(ControlLoop, self).initialize(sim, sim_data)

		self.errorInElongationRate = None
		self.proportionalTerm = None
		self.integralTerm = None

		self.bias = None

		self.rRnaSynthRate_expected = None
		self.rRnaSynthRate_updated = None

	# Allocate memory
	def allocate(self):
		super(ControlLoop, self).allocate()

		self.errorInElongationRate = 0.
		self.proportionalTerm = 0.
		self.integralTerm = 0.
		self.rRnaSynthRate_updated = 0.
		self.rRnaSynthRate_expected = 0.

		self.bias = 0.

	def update(self):
		pass

	def tableCreate(self, tableWriter):
		pass


	def tableAppend(self, tableWriter):
		tableWriter.append(
			time = self.time(),
			simulationStep = self.simulationStep(),
			errorInElongationRate = self.errorInElongationRate,
			proportionalTerm = self.proportionalTerm,
			integralTerm = self.integralTerm,
			rRnaSynthRate_expected = self.rRnaSynthRate_expected,
			rRnaSynthRate_updated = self.rRnaSynthRate_updated,
			bias = self.bias,
			)
