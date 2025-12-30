"""
Complexation

Macromolecular complexation sub-model. Encodes molecular simulation of macromolecular complexation

TODO:
- allow for shuffling when appropriate (maybe in another process)
- handle protein complex dissociation
"""

import numpy as np
from stochastic_arrow import StochasticSystem

import wholecell.processes.process

# Maximum unsigned int value + 1 for randint() to seed srand from C stdlib
RAND_MAX = 2**31

class Complexation(wholecell.processes.process.Process):
	""" Complexation """

	_name = "Complexation"

	# Constructor
	def __init__(self):

		super(Complexation, self).__init__()

	# Construct object graph
	def initialize(self, sim, sim_data):
		super(Complexation, self).initialize(sim, sim_data)

		# Create matrices and vectors that describe reaction stoichiometries
		self.stoichMatrix = sim_data.process.complexation.stoich_matrix().astype(np.int64)
		self._stoichMatrix = sim_data.process.complexation.stoich_matrix_monomers().astype(np.int64)

		# semi-quantitative rate constants
		self.rates = sim_data.process.complexation.rates

		# build stochastic system simulation
		seed = self.randomState.randint(RAND_MAX)
		self.system = StochasticSystem(self.stoichMatrix.T, random_seed=seed)

		# Build views
		moleculeNames = sim_data.process.complexation.molecule_names
		self.molecules = self.bulkMoleculesView(moleculeNames)

		# Extract relevant monomer information within the internal states:
		self.bulkMolecules = sim.internal_states["BulkMolecules"]
		self.bulk_molecule_IDs = self.bulkMolecules.container.objectNames()

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
		molecule_dict = {mol: i for i, mol in enumerate(self.bulk_molecule_IDs)}

		# Get indexes of all relevant bulk molecules:
		self.complexation_complex_idx = np.array([molecule_dict[x] for x in complexation_complex_IDs])


	def calculateRequest(self):
		moleculeCounts = self.molecules.total_counts()

		result = self.system.evolve(
			self._sim.timeStepSec(), moleculeCounts, self.rates)
		updatedMoleculeCounts = result['outcome']

		self.molecules.requestIs(np.fmax(moleculeCounts - updatedMoleculeCounts, 0))


	def evolveState(self):
		moleculeCounts = self.molecules.counts()

		result = self.system.evolve(
			self._sim.timeStepSec(), moleculeCounts, self.rates)
		updatedMoleculeCounts = result['outcome']
		events = result['occurrences']

		self.molecules.countsIs(updatedMoleculeCounts)

		# Write outputs to listeners
		self.writeToListener("ComplexationListener", "complexationEvents", events)

		# Determine the total counts of each complex:
		bulkMoleculeCounts = self.bulkMolecules.container.counts()
		complex_counts = bulkMoleculeCounts[self.complexation_complex_idx]
		self.writeToListener("ComplexationListener", "complexCounts", complex_counts)

		# Determine how many free monomers were made into complexes:
		downstream_molecule_changes = np.negative(np.dot(self._stoichMatrix, events)) # np.negative() makes the counts positive
		monomer_changes = np.zeros(len(self.monomer_IDs), np.int64)
		monomer_changes[self.matching_monomer_indices] = downstream_molecule_changes[self.matching_molecule_indices]
		self.writeToListener("ComplexationListener", "monomersComplexed", monomer_changes)

		# Determine how the number of monomers in complexes exist for each monomer:
		monomers_in_complexes = np.negative(np.dot(self._stoichMatrix, complex_counts)) # np.negative makes the monomers within it turn positive
		complexed_monomers = np.zeros(len(self.monomer_IDs), np.int64)
		complexed_monomers[self.matching_monomer_indices] = monomers_in_complexes[self.matching_molecule_indices]
		self.writeToListener("ComplexationListener", "complexedMonomerCounts", complexed_monomers)

		# TODO: add complexes generated?
		hi = 5










