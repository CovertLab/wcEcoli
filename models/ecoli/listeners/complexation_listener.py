#!/usr/bin/env python

"""
ComplexationListener

Records dynamics of complexation output.

"""

import numpy as np

import wholecell.listeners.listener

class ComplexationListener(wholecell.listeners.listener.Listener):
	""" ComplexationListener """

	_name = "ComplexationListener"

	# Constructor
	def __init__(self, *args, **kwargs):
		super(ComplexationListener, self).__init__(*args, **kwargs)


	# Construct object graph
	def initialize(self, sim, sim_data):
		super(ComplexationListener, self).initialize(sim, sim_data)

		self.monomerIDs = sim_data.process.translation.monomer_data["id"].tolist()
		self.complexIDs = sim_data.process.complexation.ids_complexes
		self.reactionIDs = sim_data.process.complexation.ids_reactions


	# Allocate memory
	def allocate(self):
		super(ComplexationListener, self).allocate()

		self.complexationEvents = np.zeros(len(self.reactionIDs), np.int64)

		self.complexCounts = np.zeros(len(self.complexIDs), np.int64)

		self.monomersComplexed = np.zeros(len(self.monomerIDs), np.int64)

		self.complexedMonomerCounts = np.zeros(len(self.monomerIDs), np.int64)

	def tableCreate(self, tableWriter):
		# TODO: add subcolumns for monomer ids in each complex, as well as counts?
		subcolumns = {
			'complexationEvents': 'reactionIDs',
			'complexCounts': 'complexIDs',
			'monomersComplexed': 'monomerIDs',
			'complexedMonomerCounts': 'monomerIDs'}


		tableWriter.writeAttributes(
			monomerIDs = self.monomerIDs,
			complexIDs = self.complexIDs,
			reactionIDs = self.reactionIDs,
			subcolumns = subcolumns)


	def tableAppend(self, tableWriter):
		tableWriter.append(
			time = self.time(),
			simulationStep = self.simulationStep(),
			complexationEvents = self.complexationEvents,
			complexCounts = self.complexCounts,
			monomersComplexed = self.monomersComplexed,
			complexedMonomerCounts = self.complexedMonomerCounts
			)
