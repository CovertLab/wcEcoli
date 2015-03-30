#!/usr/bin/env python

"""
Rna Counts

Rna counts listener. Tracks rna counts.

@author: Kalli Kappel
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 11/10/2014
"""

from __future__ import division

import numpy as np

import wholecell.listeners.listener
from wholecell.utils.fitting import normalize
from wholecell.utils import units

class rnaCounts(wholecell.listeners.listener.Listener):
	""" rnaCounts """

	_name = 'rnaCounts'

	# Constructor
	def __init__(self, *args, **kwargs):
		super(rnaCounts, self).__init__(*args, **kwargs)

		self.countUnits = "counts"
		self.probUnits = "dimensionless"

	# Construct object graph
	def initialize(self, sim, kb):
		super(rnaCounts, self).initialize(sim, kb)

		self.bulkMolecules = sim.states['BulkMolecules']

		self.states = sim.states

		self.sim = sim

		self.rnaIds = kb.process.transcription.rnaData['id']
		self.rnaIdxs = [ # TODO: use a bulk container view?
			np.where(
				x == self.bulkMolecules._moleculeIDs
				)[0][0] for x in self.rnaIds
			]
		self.rnaView = self.bulkMolecules.container.countsView(self.rnaIds)
		self.rnaSynthProb = None

	# Allocate memory
	def allocate(self):
		super(rnaCounts, self).allocate()

		self.rnaCounts = np.zeros(len(self.rnaIds), dtype = np.int64)
		self.rnaSynthProb = np.zeros(len(self.rnaIds), dtype = np.int64)

	def update(self):
		self.rnaCounts = self.rnaView.counts()

	def tableCreate(self, tableWriter):
		# Store units as metadata
		tableWriter.writeAttributes( # TODO: reconsider attribute names
			rnaCounts = self.countUnits,
			rnaSynthProb = self.probUnits
			)

	def tableAppend(self, tableWriter):
		tableWriter.append(
			time = self.time(),
			timeStep = self.timeStep(),
			rnaCounts = self.rnaCounts,
			rnaSynthProb = self.rnaSynthProb,
			)
