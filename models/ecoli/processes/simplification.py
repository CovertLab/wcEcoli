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

		self.stoichMatrix = kb.process.simplification.stoichMatrix

		# Get probabilities of dissociation
		self.pDissoc = kb.process.simplification.pDissoc.astype(np.double)

		# Build views

		moleculeNames = kb.process.simplification.moleculeIds
		complexNames = kb.process.simplification.complexIds

		self.molecules = self.bulkMoleculesView(moleculeNames)
		self.complexes = self.bulkMoleculesView(complexNames)

		# Set priority to a low value

		self.bulkMoleculesRequestPriorityIs(REQUEST_PRIORITY_SIMPLIFICATION)


	def calculateRequest(self):

		self.complexes.requestAll()


	def evolveState(self):

		complexesToDissociate = self.randomState.binomial(
			self.complexes.counts(),
			self.pDissoc
			)

		self.molecules.countsInc(
			np.dot(self.stoichMatrix, complexesToDissociate)
			)