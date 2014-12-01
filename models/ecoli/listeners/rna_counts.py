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
import tables

import wholecell.listeners.listener
from wholecell.utils.fitting import normalize
from wholecell.utils import units

class rnaCounts(wholecell.listeners.listener.Listener):
	""" rnaCounts """

	_name = 'rnaCounts'

	# Constructor
	def __init__(self, *args, **kwargs):
		super(rnaCounts, self).__init__(*args, **kwargs)

		self.usageUnits = "counts"

	# Construct object graph
	def initialize(self, sim, kb):
		super(rnaCounts, self).initialize(sim, kb)

		self.bulkMolecules = sim.states['BulkMolecules']

		self.states = sim.states

		self.sim = sim

		self.rnaIds = kb.rnaData['id']
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

	def pytablesCreate(self, h5file, expectedRows):

		shape = self.rnaCounts.shape
		# Columns
		d = {
			"time": tables.Float64Col(),
			"timeStep": tables.Int64Col(),
			"rnaCounts": tables.UInt64Col(shape),
			"rnaSynthProb": tables.Float64Col(shape),
			}

		# Create table
		# TODO: Add compression options (using filters)
		t = h5file.create_table(
			h5file.root,
			self._name,
			d,
			title = self._name,
			filters = tables.Filters(complevel = 9, complib="zlib"),
			expectedrows = expectedRows
			)

		# Store units as metadata
		t.attrs.rnaCounts_units = self.usageUnits
		#t.attrs.rnaIds = self.rnaIds


	def pytablesAppend(self, h5file):

		t = h5file.get_node("/", self._name)
		entry = t.row

		entry["time"] = self.time()
		entry["timeStep"] = self.timeStep()
		entry["rnaCounts"] = self.rnaCounts
		entry["rnaSynthProb"] = self.rnaSynthProb

		entry.append()

		t.flush()
