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

	def calculateRequest(self):

		# First, determine how many complexes dissociate into proteins:
		complex_dissociation_rates = np.zeros(len(self.complex_IDs))
		# For now, only have one complex dissociating: CPLX0-7705 -> PD00196 + PD00196
		complex_dissociation_rates[498] = self._proteinDegRates()[3798]
		total_complex_counts = self.complexes.total_counts()


		hi = 5

		nComplexesToDissociate = np.fmin(
			self.randomState.poisson(complex_dissociation_rates * total_complex_counts),
			total_complex_counts
			)

		# TODO: add a dissociation listener in evolve state, just have it get passed through and read there

		hi = 5
		# Add dissociated proteins from complexes to the protein counts:
		update_proteins = np.zeros(len(self.protein_IDs))
		update_proteins[3798] = nComplexesToDissociate[498] * 2	# CPLX0-7705 dissociates into 2 PD00196 proteins
		updated_counts = update_proteins + self.proteins.total_counts()
		hi = 5
		# TODO: FIGURE out what .countsInc() updates (does not seem to update protiens.total_counts(), but maybe it does update proteins.counts()?
		# TODO: need to actually update the protein counts, cannot simply just temporarily add them here
		# ^ could it be as simple as adding then subtracting them?

		# Determine how many proteins to degrade based on the degradation rates and counts of each protein
		nProteinsToDegrade = np.fmin(
			self.randomState.poisson(self._proteinDegRates() * updated_counts),
			updated_counts
			)

		# TODO: also update hydrolysis reactions to account for complexes dissociating

		# Determine the number of hydrolysis reactions
		nReactions = np.dot(self.proteinLengths.asNumber(), nProteinsToDegrade)

		# Determine the amount of water required to degrade the selected proteins
		# Assuming one N-1 H2O is required per peptide chain length N
		self.h2o.requestIs(nReactions - np.sum(nProteinsToDegrade))
		self.proteins.requestIs(nProteinsToDegrade)

		# carry the dissociatated complexes to the evolve state:
		self.complexes.requestIs(nComplexesToDissociate) # does this need to be negative?
		self.proteinsProduced = update_proteins


	def evolveState(self):

		# Degrade selected proteins, release amino acids from those proteins back into the cell, 
		# and consume H_2O that is required for the degradation process
		self.metabolites.countsInc(np.dot(
			self.proteinDegSMatrix,
			self.proteins.counts()
			))

		# Record how many monomers were calculated to degrade:
		counts_degraded = self.proteins.counts()
		self.writeToListener("MonomerCounts", "monomersDegraded", counts_degraded)

		# Reset the degraded protein counts:
		self.proteins.countsIs(0)

		#

		# TODO: since complexation happened already, update the complexedMonomers
	#  count listener here, as well as the complexCounts listener and read the
	#  complexes dissociated listener here (and maybe add monomersProducedViaDissociation listener?)



	def _proteinDegRates(self):
		return self.rawDegRate * self.timeStepSec()
