"""
ProteinDegradation

Protein degradation sub-model. Encodes molecular simulation of protein degradation as a Poisson process

TODO:
- protein complexes
- add protease functionality
"""

from __future__ import absolute_import, division, print_function

import numpy as np

import wholecell.processes.process
from wholecell.utils.constants import REQUEST_PRIORITY_DEGRADATION
from wholecell.utils import units
from wholecell.utils.migration.write_json import write_json

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

		super(ProteinDegradation, self).__init__()

	# Construct object graph
	def initialize(self, sim, sim_data):
		super(ProteinDegradation, self).initialize(sim, sim_data)

		# Load protein degradation rates (based on N-end rule)
		self.rawDegRate = sim_data.process.translation.monomer_data['deg_rate'].asNumber(1 / units.s)

		shuffleIdxs = None
		if hasattr(sim_data.process.translation, "monomerDegRateShuffleIdxs") and sim_data.process.translation.monomerDegRateShuffleIdxs is not None:
			shuffleIdxs = sim_data.process.translation.monomerDegRateShuffleIdxs
			self.rawDegRate = self.rawDegRate[shuffleIdxs]

		# Build metabolite IDs for S matrix
		h2oId = [sim_data.molecule_ids.water]
		metaboliteIds = sim_data.molecule_groups.amino_acids + h2oId
		aaIdxs = np.arange(0, len(sim_data.molecule_groups.amino_acids))
		h2oIdx = metaboliteIds.index(sim_data.molecule_ids.water)

		# Build protein IDs for S matrix
		self.proteinIds = sim_data.process.translation.monomer_data['id']

		# Load protein length
		self.proteinLengths = sim_data.process.translation.monomer_data['length']

		# Build S matrix
		self.proteinDegSMatrix = np.zeros((len(metaboliteIds), len(self.proteinIds)), np.int64)
		self.proteinDegSMatrix[aaIdxs, :] = np.transpose(sim_data.process.translation.monomer_data['aa_counts'].asNumber())
		self.proteinDegSMatrix[h2oIdx, :]  = -(np.sum(self.proteinDegSMatrix[aaIdxs, :], axis = 0) - 1)

		# Build Views
		self.metabolites = self.bulkMoleculesView(metaboliteIds)
		self.h2o = self.bulkMoleculeView(sim_data.molecule_ids.water)
		self.proteins = self.bulkMoleculesView(self.proteinIds)

		self.bulkMoleculesRequestPriorityIs(REQUEST_PRIORITY_DEGRADATION)

		# saving updates
		self.save_time = 1
		self.update_to_save = {}
		self.saved = False

		self.update_to_save["protein_ids"] = self.proteinIds
		self.update_to_save["metabolite_ids"] = metaboliteIds


	def calculateRequest(self):

		# Determine how many proteins to degrade based on the degradation rates and counts of each protein
		nProteinsToDegrade = np.fmin(
			self.randomState.poisson(self._proteinDegRates() * self.proteins.total_counts()),
			self.proteins.total_counts()
			)

		# Determine the number of hydrolysis reactions
		nReactions = np.dot(self.proteinLengths.asNumber(), nProteinsToDegrade)

		# Determine the amount of water required to degrade the selected proteins
		# Assuming one N-1 H2O is required per peptide chain length N
		self.h2o.requestIs(nReactions - np.sum(nProteinsToDegrade))
		self.proteins.requestIs(nProteinsToDegrade)

		# update for migration
		self.update_to_save["proteins_to_degrade"] = nProteinsToDegrade
		# if self._sim.time() > self.save_time:
		# 	import ipdb; ipdb.set_trace()

	def evolveState(self):

		# Degrade selected proteins, release amino acids from those proteins back into the cell, 
		# and consume H_2O that is required for the degradation process
		metabolite_update = np.dot(
			self.proteinDegSMatrix,
			self.proteins.counts()
		)
		self.metabolites.countsInc(metabolite_update)
		self.proteins.countsIs(0)

		# update for migration
		self.update_to_save["metabolite_update"] = metabolite_update

		# save the update
		if not self.saved and "proteins_to_degrade" in self.update_to_save.keys() and self._sim.time() > self.save_time:
			write_json(f'out/migration/prot_deg_update_t{int(self._sim.time())}.json', self.update_to_save)
			self.saved = True
			#import ipdb; ipdb.set_trace()

	def _proteinDegRates(self):
		return self.rawDegRate * self.timeStepSec()
