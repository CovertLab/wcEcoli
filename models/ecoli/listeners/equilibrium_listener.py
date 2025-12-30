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

		self.monomerIDs = sim_data.process.translation.monomer_data["id"].tolist()
		self.complexIDs = sim_data.process.equilibrium.ids_complexes
		self.reactionIDs = sim_data.process.equilibrium.rxn_ids


	# Allocate memory
	def allocate(self):
		super(EquilibriumListener, self).allocate()

		self.reactionRates = np.zeros(len(self.reactionIDs), np.float64)

		self.complexCounts = np.zeros(len(self.complexIDs), np.int64)

		self.monomersComplexed = np.zeros(len(self.monomerIDs), np.int64)

		self.complexedMonomerCounts = np.zeros(len(self.monomerIDs), np.int64)


	def tableCreate(self, tableWriter):
		subcolumns = {
			'reactionRates': 'reactionIDs',
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
			reactionRates = self.reactionRates,
			complexCounts = self.complexCounts,
			monomersComplexed = self.monomersComplexed,
			complexedMonomerCounts = self.complexedMonomerCounts
			)
