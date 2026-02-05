"""
ProteinDegradation

Protein degradation sub-model. Encodes molecular simulation of protein degradation as a Poisson process

TODO:
- protein complexes
- add protease functionality
"""

import numpy as np

import wholecell.processes.process
from wholecell.utils.constants import REQUEST_PRIORITY_DEGRADATION
from wholecell.utils import units

class ProteinDegradation(wholecell.processes.process.Process):
	""" ProteinDegradation """

	_name = "ProteinDegradation"

	# Constructor
	def __init__(self):

		# Parameters
		self.proteinLengths = None		# Protein lengths
		self.proteinDegSMatrix = None	# Protein degradation stoichiometry matrix [metabolite x rna]

		# Views
		self.metabolites = None
		self.h2o = None
		self.proteins = None
		self.complexes = None

		super(ProteinDegradation, self).__init__()

	# Construct object graph
	def initialize(self, sim, sim_data):
		super(ProteinDegradation, self).initialize(sim, sim_data)

		# Load protein degradation rates (based on N-end rule)
		self.rawDegRate = sim_data.process.translation.monomer_data['deg_rate'].asNumber(1 / units.s)

		# Build metabolite IDs for S matrix
		h2oId = [sim_data.molecule_ids.water]
		metaboliteIds = sim_data.molecule_groups.amino_acids + h2oId
		aaIdxs = np.arange(0, len(sim_data.molecule_groups.amino_acids))
		h2oIdx = metaboliteIds.index(sim_data.molecule_ids.water)

		# Build protein IDs for S matrix
		proteinIds = sim_data.process.translation.monomer_data['id']

		# Load protein length
		self.proteinLengths = sim_data.process.translation.monomer_data['length']

		# Build S matrix
		self.proteinDegSMatrix = np.zeros((len(metaboliteIds), len(proteinIds)), np.int64)
		self.proteinDegSMatrix[aaIdxs, :] = np.transpose(sim_data.process.translation.monomer_data['aa_counts'].asNumber())
		self.proteinDegSMatrix[h2oIdx, :]  = -(np.sum(self.proteinDegSMatrix[aaIdxs, :], axis = 0) - 1)

		# Build Views
		self.metabolites = self.bulkMoleculesView(metaboliteIds)
		self.h2o = self.bulkMoleculeView(sim_data.molecule_ids.water)
		self.proteins = self.bulkMoleculesView(proteinIds)

		# Complex view for future use
		# TODO: integrate complexes into degradation
		complexIds = np.array(sim_data.process.complexation.ids_complexes)
		self.complexes = self.bulkMoleculesView(complexIds)

		# Add the IDs for proteins and complexes to self:
		self.protein_IDs = proteinIds
		self.complex_IDs = complexIds

		# Obtain Avogadro's number for conversion of counts to molar concentration:
		self.n_avogadro = sim_data.constants.n_avogadro.asNumber(units.mol ** -1)
		self.bulkMoleculesRequestPriorityIs(REQUEST_PRIORITY_DEGRADATION)

		# proteins of interest:
		the_proteins = ["EG12179-MONOMER[c]", "CFA-MONOMER[c]", "G7596-MONOMER[c]", "EG11111-MONOMER[c]", "ADHP-MONOMER[c]", "EG12177-MONOMER[c]"]
		protein_ID_to_index = {p: i for i, p in enumerate(self.protein_IDs)}
		complex_ID_to_index = {c: i for i, c in enumerate(self.complex_IDs)}
		self.protein_indices = [protein_ID_to_index[p] for p in the_proteins]
		self.protein_IDs_to_index = {p: protein_ID_to_index[p] for p in the_proteins}

		# Figure out the complexes and associated stoich:
		molecules_to_all_downstream_compelexes = sim_data.process.complexation.molecules_to_all_downstream_complexes_dict
		self.proteins_to_complexes = {}
		self.proteins_to_complex_stoich = {}
		self.complex_IDs_to_index = {}
		for p in the_proteins:
			complexes = molecules_to_all_downstream_compelexes[p]
			complex_name = next(iter(complexes.keys()))
			self.proteins_to_complexes[p] = complex_name
			# Get the stoich for the complex:
			stoich = complexes[complex_name][1]['stoichiometry']
			self.proteins_to_complex_stoich[p] = stoich
			# Get the index of the complex:
			self.complex_IDs_to_index[complex_name] = complex_ID_to_index[complex_name]

		self.complex_indices = [complex_ID_to_index[c] for c in self.complex_IDs_to_index.keys()]

		# Make a mini matrix? # TODO: determine if this is needed and how to scale up
		self.matrix_monomers_i = np.array(self.protein_indices) # locations of the proteins of interest in self.protein_IDs
		self.matrix_complexes_j = np.array(self.complex_indices) # locations of the complexes of interest in self.complex_IDs
		self.matrix_stoich_v = np.array(self.proteins_to_complex_stoich.items()) # stoich values for the proteins of interest
		i = len(self.matrix_monomers_i)
		j = len(self.matrix_complexes_j)
		mini_matrix = np.zeros((j, i))

		#mini_matrix[j, i] = self.matrix_stoich_v
		#self.mini_matrix = mini_matrix

		for protein, stoich in self.proteins_to_complex_stoich.items():
			protein_index = protein_ID_to_index[protein]
			matrix_p_idx = np.where(self.matrix_monomers_i == protein_index)  # Find protein index
			complex_name = self.proteins_to_complexes[protein]
			complex_index = complex_ID_to_index[complex_name]
			matrix_c_idx = np.where(self.matrix_complexes_j == complex_index) # Find complex index

			if protein_index > 0 and complex_index > 0:
				mini_matrix[matrix_c_idx, matrix_p_idx] = stoich  # Set the stoichiometry

		self.mini_matrix = mini_matrix  # Store the result
		hi = 5
		# Determine the dissociation rates for each complex of interest:
		self.complex_dissociation_rates = np.zeros(len(self.complex_IDs))
		# Find the dissociation rate for each complex of interest:
		for protein in self.proteins_to_complexes.keys():
			protein_deg_rate = self._proteinDegRates()[self.protein_IDs_to_index[protein]]
			complex_name = self.proteins_to_complexes[protein]
			complex_index = self.complex_IDs_to_index[complex_name]
			self.complex_dissociation_rates[complex_index] = protein_deg_rate

	def calculateRequest(self):
		# TODO: need to add proteins dissociated listener and have that be something that is added to the protein counts officially!
		# TODO: make that a sepearte listener from degraded proteins AND add a complexes dissociated listener


		# First, determine how many complexes dissociate into proteins:
		self.total_complex_counts_before = self.complexes.total_counts() # this has numbers!
		self.complex_counts_before = self.complexes.counts() # all zeros

		# Determine how many complexes to dissociate based on the dissociation rates and counts of each complex:
		nComplexesToDissociate = np.fmin(
			self.randomState.poisson(self.complex_dissociation_rates * self.complexes.total_counts()),
			self.complexes.total_counts()
			)

		# TODO: need to dissociate the complexes here! try Dec first then try requestIs if that does not work
		# TODO: convention seems to have countsDec in the evolveState, so maybe have this there if I decide to do countsDec and not evolve state??
		self.complexes.requestIs(nComplexesToDissociate)
		hi = 5


		# Add dissociated proteins from complexes to the protein counts:
		update_proteins = np.zeros(len(self.protein_IDs))
		# Use the matrix to calculate the number of proteins dissociated by multiplying the nComplexesToDissociate by the mini matrix:
		proteins_dissociated = np.dot(self.mini_matrix, nComplexesToDissociate[self.matrix_complexes_j])
		for idx, protein_index in enumerate(self.matrix_monomers_i):
			# todo: would it also be ok to just do idx and then find protein_index from self.matrix_monomers_i[idx]?
			update_proteins[protein_index] += proteins_dissociated[idx]

		# TODO: have something that is saved here and written to a listener once it is in evolveState()
		self.proteinsProduced = update_proteins
		self.complexesDissociated = nComplexesToDissociate
		#self.complexes.countsDec(nComplexesToDissociate)  # dissociate selected complexes
		# todo: add requestIs statement for complexes dissociated listener here

		# figure out how to actually update total_counts() as a result:
		self.proteinstotalcountsbefore = self.proteins.total_counts()
		self.protein_counts = self.proteins.counts()


		# Determine the updated protein counts after accounting for dissociated complexes:
		updated_counts = update_proteins + self.proteins.total_counts()


		# TODO: FIGURE out what .countsInc() updates (does not seem to update protiens.total_counts(), but maybe it does update proteins.counts()?
		# TODO: need to actually update the protein counts, cannot simply just temporarily add them here
		# ^ could it be as simple as adding then subtracting them?

		# Determine how many proteins to degrade based on the degradation rates and counts of each protein
		nProteinsToDegrade = np.fmin(
			self.randomState.poisson(self._proteinDegRates() * updated_counts),
			updated_counts
			)

		self.proteinsDegraded = nProteinsToDegrade

		# TODO: also update hydrolysis reactions (or some other type of molecule that complexation uses?) to account for complexes dissociating
		# Determine the number of hydrolysis reactions
		nReactions = np.dot(self.proteinLengths.asNumber(), nProteinsToDegrade)

		# Determine the amount of water required to degrade the selected proteins
		# Assuming one N-1 H2O is required per peptide chain length N
		self.h2o.requestIs(nReactions - np.sum(nProteinsToDegrade))


		#self.proteins.countsDec(nProteinsToDegrade)  # degrade selected proteins
		self.intermediate_counts = self.proteins.counts()
		self.intermediate_total_counts = self.proteins.total_counts()


		# Figure out the total change in the proteins, but remember that the proteins degraded are positive:
		self.netActivity = nProteinsToDegrade - update_proteins

		# TODO: can i put the counts increased before this at least????
		self.proteins.requestIs(nProteinsToDegrade) # there is no need to do a countsDec here since it will be accepted? I think not matter what?




	def evolveState(self):

		# Degrade selected proteins, release amino acids from those proteins back into the cell, 
		# and consume H_2O that is required for the degradation process
		hi = 5
		self.metabolites.countsInc(np.dot(
			self.proteinDegSMatrix,
			self.proteins.counts()
			))

		hi = 5


		# Update the counts of the complexes, the counts of the proteins (due to dissociation)
		#self.complexes.countsDec(self.complexesDissociated) # todo: double check if the indexing is correct becuase I have no clue why this did not work. should not have been negative.
		self.writeToListener("ComplexationListener", "complexesDissociated", self.complexesDissociated)
		# Update the protein counts due to dissociation:
		self.proteins.countsInc(self.proteinsProduced)
		self.writeToListener("MonomerCounts", "monomersProducedViaDissociation", self.proteinsProduced) # todo: this is not great practice given the monomers were used in the previous step, but not sure how to update that there accordingly yet.
		# Record how many monomers were calculated to degrade:
		counts_degraded = self.proteins.counts()
		self.writeToListener("MonomerCounts", "monomersDegraded", counts_degraded)

		hi = 4

		# Reset the degraded protein counts:
		self.proteins.countsIs(0)

		#

		# TODO: since complexation happened already, update the complexedMonomers
	#  count listener here, as well as the complexCounts listener and read the
	#  complexes dissociated listener here (and maybe add monomersProducedViaDissociation listener?)



	def _proteinDegRates(self):
		return self.rawDegRate * self.timeStepSec()
