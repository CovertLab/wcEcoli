#!/usr/bin/env python

"""
Complexation

Macromolecular complexation sub-model. Encodes molecular simulation of macromolecular complexation

TODO:
- allow for shuffling when appropriate (maybe in another process)
- handle protein complex dissociation

@author: Derek Macklin
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 4/4/2013
"""

from __future__ import division

import numpy as np
from arrow import StochasticSystem

import wholecell.processes.process

# Maximum unsigned int value + 1 for randint() to seed srand from C stdlib
RAND_MAX = 2**31

class Complexation(wholecell.processes.process.Process):
	""" Complexation """

	_name = "Complexation"

	# Constructor
	def __init__(self):

		super(Complexation, self).__init__()

	# Construct object graph
	def initialize(self, sim, sim_data):
		super(Complexation, self).initialize(sim, sim_data)

		self.gillespie_time_step = 0.001 # instead of self._sim.timeStepSec()

		# Create matrices and vectors that describe reaction stoichiometries 
		complexation_forward_matrix = sim_data.process.complexation.stoichMatrix().astype(np.int64)
		complexation_reverse_matrix = -sim_data.process.complexation.stoichMatrix().astype(np.int64)
		self.stoichMatrix = np.append(complexation_forward_matrix, complexation_reverse_matrix, 1)

		# semi-quantitative rate constants
		forward_rates = sim_data.process.complexation.rates
		reverse_rates = np.repeat(10, len(forward_rates))
		self.rates = np.append(forward_rates, reverse_rates)
		# import ipdb; ipdb.set_trace()

		# build stochastic system simulation
		seed = self.randomState.randint(RAND_MAX)
		self.system = StochasticSystem(self.stoichMatrix.T, self.rates, random_seed=seed)

		# Build views
		moleculeNames = sim_data.process.complexation.moleculeNames
		self.molecules = self.bulkMoleculesView(moleculeNames)
		# import ipdb; ipdb.set_trace()


	def calculateRequest(self):
		moleculeCounts = self.molecules.total_counts()

		result = self.system.evolve(self.gillespie_time_step, moleculeCounts)
		updatedMoleculeCounts = result['outcome']

		self.molecules.requestIs(np.fmax(moleculeCounts - updatedMoleculeCounts, 0))


	def evolveState(self):
		moleculeCounts = self.molecules.counts()

		result = self.system.evolve(self.gillespie_time_step, moleculeCounts)
		updatedMoleculeCounts = result['outcome']
		events = result['occurrences']
		if len(result['time']) > 0:
			gillespie_time = max(result['time'])
		else:
			gillespie_time = 0
		# import ipdb; ipdb.set_trace()

		self.molecules.countsIs(updatedMoleculeCounts)

		# Write outputs to listeners
		self.writeToListener("ComplexationListener", "complexationEvents", events)
		self.writeToListener("ComplexationListener", "gillespieTime", gillespie_time)
