#!/usr/bin/env python

"""
FBAResults

Records dynamics of FBA output.

@author: John Mason
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 6/24/2014
"""

from __future__ import division

import numpy as np

import wholecell.listeners.listener

class FBAResults(wholecell.listeners.listener.Listener):
	""" FBAResults """

	_name = "FBAResults"

	# Constructor
	def __init__(self, *args, **kwargs):
		super(FBAResults, self).__init__(*args, **kwargs)

		self.countUnits = "counts"

	# Construct object graph
	def initialize(self, sim, sim_data):
		super(FBAResults, self).initialize(sim, sim_data)

		self.metabolism = sim.processes["Metabolism"]

		self.objectiveValue = 0.0


	# Allocate memory
	def allocate(self):
		super(FBAResults, self).allocate()

		fba = self.metabolism.fba

		self.reactionIDs = fba.reactionIDs()
		self.externalMoleculeIDs = fba.externalMoleculeIDs()
		self.outputMoleculeIDs = fba.outputMoleculeIDs()

		self.reactionFluxes = np.zeros(len(self.reactionIDs), np.float64)
		self.externalExchangeFluxes = np.zeros(len(self.externalMoleculeIDs), np.float64)
		self.outputFluxes = np.zeros(len(self.outputMoleculeIDs), np.float64)
		self.objectiveComponents = np.zeros_like(self.outputFluxes)
		self.deltaMetabolites = np.zeros(len(self.outputMoleculeIDs), np.float64)
		self.PCC = 0.0
		self.PCClog = 0.0
		self.slope = 0.0
		self.slopelog = 0.0


	def tableCreate(self, tableWriter):
		tableWriter.writeAttributes(
			reactionIDs = list(self.reactionIDs),
			externalMoleculeIDs = self.externalMoleculeIDs,
			outputMoleculeIDs = self.outputMoleculeIDs,
			PCC = self.countUnits,
			PCClog = self.countUnits,
			slope = self.countUnits,
			slopelog = self.countUnits,
			)


	def tableAppend(self, tableWriter):
		tableWriter.append(
			time = self.time(),
			simulationStep = self.simulationStep(),
			reactionFluxes = self.reactionFluxes,
			externalExchangeFluxes = self.externalExchangeFluxes,
			outputFluxes = self.outputFluxes,
			objectiveValue = self.objectiveValue,
			objectiveComponents = self.objectiveComponents,
			deltaMetabolites = self.deltaMetabolites,
			PCC = self.PCC,
			PCClog = self.PCClog,
			slope = self.slope,
			slopelog = self.slopelog,
			)
