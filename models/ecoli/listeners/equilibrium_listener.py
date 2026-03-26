#!/usr/bin/env python

"""
EquilibriumListener

Records dynamics of equilibrium output.

"""

import numpy as np

import wholecell.listeners.listener

class EquilibriumListener(wholecell.listeners.listener.Listener):
	""" EquilibriumListener """

	_name = "EquilibriumListener"

	# Constructor
	def __init__(self, *args, **kwargs):
		super(EquilibriumListener, self).__init__(*args, **kwargs)


	# Construct object graph
	def initialize(self, sim, sim_data):
		super(EquilibriumListener, self).initialize(sim, sim_data)

		# Get the IDs of the equilibrium reactions:
		self.reactionIDs = sim_data.process.equilibrium.rxn_ids

	# Allocate memory
	def allocate(self):
		super(EquilibriumListener, self).allocate()

		self.reactionRates = np.zeros(len(self.reactionIDs), np.float64)

		self.complexationEvents = np.zeros(len(self.reactionIDs), np.int64)

	def tableCreate(self, tableWriter):
		subcolumns = {
			'reactionRates': 'reactionIDs',
			'complexationEvents': 'reactionIDs'
			}

		tableWriter.writeAttributes(
			reactionIDs = self.reactionIDs,
			subcolumns = subcolumns)


	def tableAppend(self, tableWriter):
		tableWriter.append(
			time = self.time(),
			simulationStep = self.simulationStep(),
			reactionRates = self.reactionRates,
			complexationEvents = self.complexationEvents
		)
