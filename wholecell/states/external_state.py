#!/usr/bin/env python

"""
External State

State variable base class. Defines the interface states expose to the simulation.

@organization: Covert Lab, Department of Bioengineering, Stanford University
"""

from __future__ import division

import numpy as np

class ExternalState(object):
	""" External State """

	_name = None

	# Constructor
	def __init__(self):
		# References to views
		self._views = []

		# Random number stream
		self.randomState = None

		self.seed = None


	# Construct state-process graph, calculate constants
	def initialize(self, sim_data, process_keys):
		self._process_keys = process_keys


	def update(self, time):
		pass


	# Allocate memory
	def allocate(self):
		pass


	# Views
	def viewAdd(self, view):
		self._views.append(view)


	# Saving

	def tableCreate(self, tableWriter):
		pass


	def tableAppend(self, tableWriter):
		pass


	# Basic accessors

	@classmethod
	def name(cls):
		return cls._name

