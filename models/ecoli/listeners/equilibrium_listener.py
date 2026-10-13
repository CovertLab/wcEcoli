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

		# Get the IDs of the equilibrium reactions and complexes:
		self.reactionIDs = sim_data.process.equilibrium.rxn_ids
		self.complexIDs = sim_data.process.equilibrium.ids_complexes

	# Allocate memory
	def allocate(self):
		super(EquilibriumListener, self).allocate()

		self.reactionRates = np.zeros(len(self.reactionIDs), np.float64)

		self.complexationEvents = np.zeros(len(self.reactionIDs), np.int64)

		# Per-complex dissociation counts written by ProteinDegradation when it
		# breaks equilibrium complexes apart to degrade their subunits:
		self.complexesDissociated = np.zeros(len(self.complexIDs), np.int64)

		# Total unique equilibrium complexes dissociated this timestep:
		self.dissociationIterations = 0

	def tableCreate(self, tableWriter):
		subcolumns = {
			'reactionRates': 'reactionIDs',
			'complexationEvents': 'reactionIDs',
			'complexesDissociated': 'complexIDs',
			}

		tableWriter.writeAttributes(
			reactionIDs = self.reactionIDs,
			complexIDs = self.complexIDs,
			subcolumns = subcolumns)


	def tableAppend(self, tableWriter):
		tableWriter.append(
			time = self.time(),
			simulationStep = self.simulationStep(),
			reactionRates = self.reactionRates,
			complexationEvents = self.complexationEvents,
			complexesDissociated = self.complexesDissociated,
			dissociationIterations = self.dissociationIterations,
		)
