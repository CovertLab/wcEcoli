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

		# Get the IDs of the complexation reactions and complexes:
		self.reactionIDs = sim_data.process.complexation.ids_reactions
		self.complexIDs = sim_data.process.complexation.ids_complexes

	# Allocate memory
	def allocate(self):
		super(ComplexationListener, self).allocate()

		self.complexationEvents = np.zeros(len(self.reactionIDs), np.int64)

		# Per-complex dissociation counts written by ProteinDegradation when it
		# breaks complexes apart to degrade their subunits:
		self.complexesDissociated = np.zeros(len(self.complexIDs), np.int64)

		# Total unique complexes dissociated this timestep (overlap accounting):
		self.dissociationIterations = 0

	def tableCreate(self, tableWriter):
		subcolumns = {
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
			complexationEvents = self.complexationEvents,
			complexesDissociated = self.complexesDissociated,
			dissociationIterations = self.dissociationIterations,
			)
