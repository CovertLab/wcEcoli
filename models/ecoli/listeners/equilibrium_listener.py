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
		self.stoichMatrixMonomers = (
			sim_data.process.equilibrium.stoich_matrix_monomers().astype(np.int64))

		# Extract relevant monomer information within the internal states:
		self.bulkMolecules = sim.internal_states["BulkMolecules"]
		bulk_molecule_IDs = self.bulkMolecules.container.objectNames()

		# Get IDs of molecules involved in complexation reactions:
		equilibrium_molecule_IDs = sim_data.process.equilibrium.molecule_names
		equilibrium_complex_IDs = sim_data.process.equilibrium.ids_complexes

		# Extract all monomer IDs so that matches within moleculeNames can be tracked:
		self.monomer_IDs = sim_data.process.translation.monomer_data["id"].tolist()

		# Find where monomers IDs are within molecules:
		matching_monomers_mask = np.isin(equilibrium_molecule_IDs, self.monomer_IDs)

		# Get the indices of the matching monomers (i.e. where it is nonzero):
		matching_indices = np.where(matching_monomers_mask)[0]

		# Obtain the indices of the monomer IDs within bulkIDs for each matching index:
		monomer_indices = []
		for i in matching_indices:
			molecule_ID = equilibrium_molecule_IDs[i]
			if molecule_ID in self.monomer_IDs:
				monomer_index = self.monomer_IDs.index(molecule_ID)
				monomer_indices.append(monomer_index)

		self.matching_monomer_indices = monomer_indices
		self.matching_molecule_indices = matching_indices

		# Construct dictionary to quickly find bulk molecule indexes from IDs:
		molecule_dict = {mol: i for i, mol in enumerate(bulk_molecule_IDs)}

		# Get indexes of all relevant bulk molecules:
		self.equilibrium_complex_idx = np.array(
			[molecule_dict[x] for x in equilibrium_complex_IDs])


	# Allocate memory
	def allocate(self):
		super(EquilibriumListener, self).allocate()

		self.reactionRates = np.zeros(len(self.reactionIDs), np.float64)

		self.complexationEvents = np.zeros(len(self.reactionIDs), np.int64)

		self.complexCounts = np.zeros(len(self.complexIDs), np.int64)

		self.freeMonomersComplexed = np.zeros(len(self.monomerIDs), np.int64)

		self.freeMonomersReleased = np.zeros(len(self.monomerIDs), np.int64)

		self.complexedMonomerCounts = np.zeros(len(self.monomerIDs), np.int64)

	def update(self):
		# Get current counts of all bulk molecules:
		bulkMoleculeCounts = self.bulkMolecules.container.counts()

		# Determine the number of complexes present in the cell currently by
		# indexing into bulkMoleculeCounts where complexed molecules are located:
		self.complexCounts = bulkMoleculeCounts[self.equilibrium_complex_idx]

		# Determine the # of monomers that are "currently" in complexes:
		monomers_in_complexes = np.negative(np.dot(self.stoichMatrixMonomers,
												   self.complexCounts))
		complexed_monomers = np.zeros(len(self.monomer_IDs), np.int64)
		complexed_monomers[self.matching_monomer_indices] = monomers_in_complexes[
			self.matching_molecule_indices]
		self.complexedMonomerCounts = complexed_monomers

	def tableCreate(self, tableWriter):
		subcolumns = {
			'reactionRates': 'reactionIDs',
			'complexationEvents': 'reactionIDs',
			'complexCounts': 'complexIDs',
			'freeMonomersComplexed': 'monomerIDs',
			'freeMonomersReleased': 'monomerIDs',
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
			complexationEvents = self.complexationEvents,
			complexCounts = self.complexCounts,
			freeMonomersComplexed = self.freeMonomersComplexed,
			freeMonomersReleased = self.freeMonomersReleased,
			complexedMonomerCounts = self.complexedMonomerCounts
			)
