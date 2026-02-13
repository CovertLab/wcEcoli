"""
Equilibrium

Equilibrium binding sub-model
"""

import numpy as np

from wholecell.utils import units
import wholecell.processes.process


class Equilibrium(wholecell.processes.process.Process):
	""" Equilibrium """

	_name = "Equilibrium"

	# Constructor
	def __init__(self):

		super(Equilibrium, self).__init__()


	# Construct object graph
	def initialize(self, sim, sim_data):
		super(Equilibrium, self).initialize(sim, sim_data)

		# Simulation options
		self.jit = sim._jit

		# Get constants
		self.nAvogadro = sim_data.constants.n_avogadro.asNumber(1 / units.mol)
		self.cellDensity = sim_data.constants.cell_density.asNumber(units.g / units.L)

		# Create matrix and method
		self.stoichMatrix = sim_data.process.equilibrium.stoich_matrix().astype(np.int64)
		self._stoichMatrix = sim_data.process.equilibrium.stoich_matrix_monomers().astype(np.int64)
		self.fluxesAndMoleculesToSS = sim_data.process.equilibrium.fluxes_and_molecules_to_SS
		self.product_indices = [idx for idx in np.where(np.any(self.stoichMatrix > 0, axis=1))[0]]

		# Build views
		self.moleculeNames = sim_data.process.equilibrium.molecule_names
		self.molecules = self.bulkMoleculesView(self.moleculeNames)

		# Extract relevant monomer information within the internal states:
		self.bulkMolecules = sim.internal_states["BulkMolecules"]
		self.bulk_molecule_IDs = self.bulkMolecules.container.objectNames()

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
		molecule_dict = {mol: i for i, mol in enumerate(self.bulk_molecule_IDs)}

		# Get indexes of all relevant bulk molecules:
		self.equilibrium_complex_idx = np.array(
			[molecule_dict[x] for x in equilibrium_complex_IDs])


	def calculateRequest(self):
		# Get molecule counts
		moleculeCounts = self.molecules.total_counts()

		# Get cell mass and volume
		cellMass = (self.readFromListener("Mass", "cellMass") * units.fg).asNumber(units.g)
		cellVolume = cellMass / self.cellDensity

		# Solve ODEs to steady state
		self.rxnFluxes, self.req = self.fluxesAndMoleculesToSS(
			moleculeCounts, cellVolume, self.nAvogadro, self.randomState,
			jit=self.jit,
			)

		# Request counts of molecules needed
		self.molecules.requestIs(self.req)


	def evolveState(self):
		# Get counts of molecules allocated to this process
		moleculeCounts = self.molecules.counts()
		rxnFluxes = self.rxnFluxes.copy()

		# If we didn't get allocated all the molecules we need, make do with
		# what we have (decrease reaction fluxes so that they make use of what
		# we have, but not more). Reduces at least one reaction every iteration
		# so the max number of iterations is the number of reactions that were
		# originally expected to occur + 1 to reach the break statement.
		max_iterations = int(np.abs(rxnFluxes).sum()) + 1
		for it in range(max_iterations):
			# Check if any metabolites will have negative counts with current reactions
			negative_metabolite_idxs = np.where(np.dot(self.stoichMatrix, rxnFluxes) + moleculeCounts < 0)[0]
			if len(negative_metabolite_idxs) == 0:
				break

			# Reduce reactions that consume metabolites with negative counts
			limited_rxn_stoich = self.stoichMatrix[negative_metabolite_idxs, :]
			fwd_rxn_idxs = np.where(np.logical_and(limited_rxn_stoich < 0, rxnFluxes > 0))[1]
			rev_rxn_idxs = np.where(np.logical_and(limited_rxn_stoich > 0, rxnFluxes < 0))[1]
			rxnFluxes[fwd_rxn_idxs] -= 1
			rxnFluxes[rev_rxn_idxs] += 1
			rxnFluxes[fwd_rxn_idxs] = np.fmax(0, rxnFluxes[fwd_rxn_idxs])
			rxnFluxes[rev_rxn_idxs] = np.fmin(0, rxnFluxes[rev_rxn_idxs])
		else:
			raise ValueError('Could not get positive counts in equilibrium with'
				' allocated molecules.')

		# Increment changes in molecule counts
		deltaMolecules = np.dot(self.stoichMatrix, rxnFluxes)
		self.molecules.countsInc(deltaMolecules)

		# Write outputs to listeners
		self.writeToListener("EquilibriumListener", "complexationEvents", rxnFluxes)
		self.writeToListener("EquilibriumListener", "reactionRates", (
			deltaMolecules[self.product_indices] / self.timeStepSec()
			))

		# Determine the total counts of each complex:
		bulkMoleculeCounts = self.bulkMolecules.container.counts()
		complex_counts = bulkMoleculeCounts[self.equilibrium_complex_idx]
		self.writeToListener("EquilibriumListener", "complexCounts", complex_counts)

		# TODO: check the np.negative sign convention here
		# Determine how many free monomers were used to generate complexes:
		opposite_deltaMolecules = np.negative(
			deltaMolecules)  # monomers that were used to form complexes will be positive here to stay consisent with the nomenclature used for monomersDegraded (in monomer_counts.py), and monomers that were generated from complex disassociation will be positive.
		free_monomers_complexed = np.zeros(len(self.monomer_IDs), np.int64)
		free_monomers_complexed[self.matching_monomer_indices] = opposite_deltaMolecules[
			self.matching_molecule_indices]
		self.writeToListener("EquilibriumListener", "freeMonomersComplexed",
							 free_monomers_complexed)

		# Determine how the number of monomers in complexes:
		monomers_in_complexes = np.negative(np.dot(self._stoichMatrix,
												   complex_counts))
		complexed_monomers = np.zeros(len(self.monomer_IDs), np.int64)
		complexed_monomers[self.matching_monomer_indices] = monomers_in_complexes[
			self.matching_molecule_indices]
		self.writeToListener("EquilibriumListener", "complexedMonomerCounts", complexed_monomers)

