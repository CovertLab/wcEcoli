#!/usr/bin/env python

"""
Simplification (Dissociation)

Dissociate complexes

@author: Derek Macklin
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 6/25/2015
"""

from __future__ import division

import numpy as np

import wholecell.processes.process
from wholecell.utils.mc_complexation import mccBuildMatrices, mccFormComplexesWithPrebuiltMatrices

from wholecell.utils.constants import REQUEST_PRIORITY_SIMPLIFICATION

class Simplification(wholecell.processes.process.Process):
	""" Complexation """

	_name = "Simplification"

	# Constructor
	def __init__(self):

		super(Simplification, self).__init__()


	# Construct object graph
	def initialize(self, sim, kb):
		super(Simplification, self).initialize(sim, kb)

		# Create stoichiometric matrix
		# (it's the opposite of complexation)

		self.stoichMatrix = -1 * kb.process.complexation.stoichMatrix().astype(np.int64, order = "F")

		# Get probabilities of dissociation
		self.pRevs = kb.process.simplification.pRevs.astype(np.double)

		# Build views

		moleculeNames = kb.process.complexation.moleculeNames
		complexNames = [moleculeNames[x] for x in np.where(self.stoichMatrix < 0)[0]]
		# Can't just use kb.process.complexation.complexNames because we want ordering to match the stoichiometric matrix

		self.molecules = self.bulkMoleculesView(moleculeNames)
		self.complexes = self.bulkMoleculesView(complexNames)

		# Set priority to a low value

		self.bulkMoleculesRequestPriorityIs(REQUEST_PRIORITY_SIMPLIFICATION)

		import ipdb; ipdb.set_trace()

	def calculateRequest(self):

		self.complexes.requestAll()


	def evolveState(self):
		print self.molecules.total()[-3], self.molecules.total()[-1]
		if self.time() > 3000:
			import ipdb; ipdb.set_trace()

		complexesToDissociate = self.randomState.binomial(
			self.complexes.counts(),
			self.pRevs
			)
		self.molecules.countsInc(
			np.dot(self.stoichMatrix, complexesToDissociate)
			)