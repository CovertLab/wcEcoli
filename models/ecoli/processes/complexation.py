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

		# semi-quantitative rate constants
		self.rates = sim_data.process.complexation.rates

		# build stochastic system simulation
		seed = self.randomState.randint(RAND_MAX)
		self.system = StochasticSystem(self.stoichMatrix.T, random_seed=seed)

		# Build views
		moleculeNames = sim_data.process.complexation.molecule_names
		self.molecules = self.bulkMoleculesView(moleculeNames)

		# Extract relevant monomer information:
		self.bulkMolecules = sim.internal_states["BulkMolecules"]
		bulk_molecule_ids = self.bulkMolecules.container.objectNames()

		# Get IDs of molecules involved in complexation reactions:
		complexation_molecule_ids = sim_data.process.complexation.molecule_names
		complexation_complex_ids = sim_data.process.complexation.ids_complexes

		# Find where the complex IDs are within the molecule IDs:
		self.complex_IDs_within_molecule_IDs = []
		for complex in complexation_complex_ids:
			matching_index = np.where(complexation_molecule_ids == complex)
			self.complex_IDs_within_molecule_IDs.append(matching_index)

		# extract all monomer IDs so that matches within moleculeNames can be tracked:
		self.monomer_ids = sim_data.process.translation.monomer_data["id"].tolist()

		# Find where monomers are within molecules so that its easy to calculate the number of proteins in the complexes per time step
		# todo: implement this (or consider just doing this in an analysis script?) need to double check that the stoich matrix is saved first

		# Construct dictionary to quickly find bulk molecule indexes from IDs
		molecule_dict = {mol: i for i, mol in enumerate(bulk_molecule_ids)}

		def get_molecule_indexes(keys):
			return np.array([molecule_dict[x] for x in keys])

		# Get indexes of all relevant bulk molecules
		self.monomer_idx = get_molecule_indexes(self.monomer_ids)
		self.complexation_molecule_idx = get_molecule_indexes(complexation_molecule_ids)
		self.complexation_complex_idx = get_molecule_indexes(complexation_complex_ids)



	def calculateRequest(self):
		moleculeCounts = self.molecules.total_counts()

		hi = 6

		result = self.system.evolve(
			self._sim.timeStepSec(), moleculeCounts, self.rates)
		updatedMoleculeCounts = result['outcome']

		hi = 8

		self.molecules.requestIs(np.fmax(moleculeCounts - updatedMoleculeCounts, 0))

		hi = 7


	def evolveState(self):
		moleculeCounts = self.molecules.counts()

		hi = 5

		result = self.system.evolve(
			self._sim.timeStepSec(), moleculeCounts, self.rates)
		updatedMoleculeCounts = result['outcome']
		events = result['occurrences']

		hi = 5

		self.molecules.countsIs(updatedMoleculeCounts)

		hi = 6

		# Write outputs to listeners
		self.writeToListener("ComplexationListener", "complexationEvents", events)
		hi = 5
		# Determine the total counts of the complexes:
		bulkMoleculeCounts = self.bulkMolecules.container.counts()
		complex_counts = bulkMoleculeCounts[self.complexation_complex_idx]
		self.writeToListener("ComplexationListener", "complexCounts", complex_counts)

		# Determine how the counts of each molecule involvd in complexation changed this timestep:
		molecule_changes = np.dot(self.stoichMatrix, events)

		# Determine how many complexes were generated this time step:
		# find where complex indexes match molecule indexes:
		complex_indexes = set(self.complexation_complex_idx, self.complexation_molecule_idx)
		# Next, use this to figure out where the complexes must be within the events matrix using these overlapps?

		# Determine the number of each monomer within the total complexes:
		# note, to do this, would need the other matrix ...


		# TODO: write the number of molecules used to generate the monomers here, not to monomer_counts.










