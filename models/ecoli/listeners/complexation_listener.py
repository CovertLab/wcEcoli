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
		self._stoichMatrix = sim_data.process.complexation.stoich_matrix_monomers().astype(
			np.int64)

		# Get IDs of all bulk molecules
		self.bulkMolecules = sim.internal_states["BulkMolecules"]
		bulk_molecule_IDs = self.bulkMolecules.container.objectNames()

		# Get IDs of molecules involved in complexation reactions:
		complexation_molecule_IDs = sim_data.process.complexation.molecule_names
		complexation_complex_IDs = sim_data.process.complexation.ids_complexes

		# Extract all monomer IDs so that matches within moleculeNames can be tracked:
		self.monomer_IDs = sim_data.process.translation.monomer_data["id"].tolist()

		# Find where monomers IDs are within molecules:
		matching_monomers_mask = np.isin(complexation_molecule_IDs, self.monomer_IDs)

		# Get the indices of the matching monomers (i.e. where it is nonzero):
		matching_indices = np.where(matching_monomers_mask)[0]

		# Obtain the indices of the monomer IDs within bulkIDs for each matching index:
		monomer_indices = []
		for i in matching_indices:
			molecule_ID = complexation_molecule_IDs[i]
			if molecule_ID in self.monomer_IDs:
				monomer_index = self.monomer_IDs.index(molecule_ID)
				monomer_indices.append(monomer_index)

		self.matching_monomer_indices = monomer_indices
		self.matching_molecule_indices = matching_indices

		# Construct dictionary to quickly find bulk molecule indexes from IDs:
		molecule_dict = {mol: i for i, mol in enumerate(bulk_molecule_IDs)}

		# Get indexes of all relevant bulk molecules:
		self.complexation_complex_idx = np.array(
			[molecule_dict[x] for x in complexation_complex_IDs])

	# Allocate memory
	def allocate(self):
		super(ComplexationListener, self).allocate()

		self.complexationEvents = np.zeros(len(self.reactionIDs), np.int64)

		self.complexCounts = np.zeros(len(self.complexIDs), np.int64)

		self.monomersComplexed = np.zeros(len(self.monomerIDs), np.int64)

		self.complexedMonomerCounts = np.zeros(len(self.monomerIDs), np.int64)

		self.complexesDissociated = np.zeros(len(self.complexIDs), np.int64)

	def initialUpdate(self):
		# Get current counts of bulk and unique molecules
		bulkMoleculeCounts = self.bulkMolecules.container.counts()

		# Update the complex counts at the start of the cell:
		self.complexCounts = bulkMoleculeCounts[self.complexation_complex_idx]

		# Update the counts of monomers currently complexed at the start of the cell:
		monomers_in_complexes = np.negative(np.dot(self._stoichMatrix,
												   self.complexCounts))  # np.negative makes the monomers within it turn positive
		complexed_monomers = np.zeros(len(self.monomer_IDs), np.int64)
		complexed_monomers[self.matching_monomer_indices] = monomers_in_complexes[
			self.matching_molecule_indices]
		self.complexedMonomerCounts = complexed_monomers

		# TODO: decide if it is ok to delete large matrices here to save memory

	def tableCreate(self, tableWriter):
		# TODO: add subcolumns for monomer ids in each complex, as well as counts?
		subcolumns = {
			'complexationEvents': 'reactionIDs',
			'complexCounts': 'complexIDs',
			'monomersComplexed': 'monomerIDs',
			'complexedMonomerCounts': 'monomerIDs',
			'complexesDissociated': 'complexIDs'}


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
			complexedMonomerCounts = self.complexedMonomerCounts,
			complexesDissociated = self.complexesDissociated
			)
