"""
ProteinDegradation

Protein degradation sub-model. Encodes molecular simulation of protein degradation as a Poisson process

TODO:
- protein complexes
- add protease functionality
"""
# todo: show Nora the above comments

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

		# add complexes:
		complexIds = np.array(sim_data.process.complexation.ids_complexes)
		self.complexes = self.bulkMoleculesView(complexIds)

		# find lon complex:
		lon_complex_idx = np.where(complexIds == 'CPLX0-2881[c]')[0][0]
		#import ipdb ; ipdb.set_trace()
		print("Lon Complex index:",lon_complex_idx)
		print("Lon Complex id:", complexIds[lon_complex_idx])

		# TODO: ask nora how she was able to access the protein complexes per timestep. she has a file where she averages over them, try to figure out where that listener is that is recording itand then how to convert to concenration once you ahve that . add it here as a self.bulkMOlecules View thing here to be able know the complex ids at the start and be able to type self.compelxes in the evolve state and

		self.bulkMoleculesRequestPriorityIs(REQUEST_PRIORITY_DEGRADATION)
		hi = 5
		#import ipdb; ipdb.set_trace()

	def calculateRequest(self):
		print("protein degradation calcREQUEST started")
		# Determine how many proteins to degrade based on the degradation rates and counts of each protein
		nProteinsToDegrade = np.fmin(
			self.randomState.poisson(self._proteinDegRates() * self.proteins.total_counts()),
			self.proteins.total_counts()
			)
		# todo: determine if fmin is finding the minimum for each protein (between either the proteins total counts or the degrdation rate times the proteins total counts)
		# todo: determine if total_counts() adds up the counts of all proteins or still separates by protein
		hi = 5
		#print(self.proteins._totalCount[3863], self.proteins.total()[3863], self.proteins._counts()[3863], self.proteins._requestedCount[3863], self.proteins.total_counts()[3863],self.proteins.counts()[3863], nProteinsToDegrade[3863])
		print(self.proteins._totalCount[403], self.proteins.counts()[403], nProteinsToDegrade[403])
		self.writeToListener('MonomerCounts', 'protein_deg_CR1__totalCount', self.proteins._totalCount.copy())

		#import ipdb; ipdb.set_trace() # CR1


		# Determine the number of hydrolysis reactions
		nReactions = np.dot(self.proteinLengths.asNumber(), nProteinsToDegrade)

		# Determine the amount of water required to degrade the selected proteins
		# Assuming one N-1 H2O is required per peptide chain length N
		self.h2o.requestIs(nReactions - np.sum(nProteinsToDegrade))
		self.proteins.requestIs(nProteinsToDegrade)
		#print(self.proteins._totalCount[3863], self.proteins.total()[3863], self.proteins._counts()[3863], self.proteins._requestedCount[3863], self.proteins.total_counts()[3863],self.proteins.counts()[3863], nProteinsToDegrade[3863])
		print(self.proteins._totalCount[403], self.proteins.counts()[403], nProteinsToDegrade[403])

		self.writeToListener('MonomerCounts', 'protein_deg_CR2__totalCount', self.proteins._totalCount.copy())
		counts_for_CR2 = self.proteins.counts()
		self.writeToListener('MonomerCounts', 'protein_deg_CR2_counts', counts_for_CR2) # todo: ask riley why .copy() doesnt work here
		#import ipdb; ipdb.set_trace()

		#import ipdb; ipdb.set_trace() # CR2
		print("protein degradation calcREQUEST ended")




	def evolveState(self):
		print("protein degradation evolveState started")
		# todo: determine if this is where the proteins are actually degraded
		# Degrade selected proteins, release amino acids from those proteins back into the cell, 
		# and consume H_2O that is required for the degradation process
		self.metabolites.countsInc(np.dot(
			self.proteinDegSMatrix,
			self.proteins.counts()
			))
		#print(self.proteins._totalCount[3863], self.proteins.total()[3863], self.proteins._counts()[3863], self.proteins._requestedCount[3863], self.proteins.total_counts()[3863],self.proteins.counts()[3863])
		print(self.proteins._totalCount[403], self.proteins.counts()[403])
		counts_for_ES1 = self.proteins.counts()
		self.writeToListener('MonomerCounts', 'protein_deg_ES1_counts', counts_for_ES1)

		#import ipdb; ipdb.set_trace() # ES1

		self.proteins.countsIs(0) # does this set the counts to zero?
		#print(self.proteins._totalCount[3863], self.proteins.total()[3863], self.proteins._counts()[3863], self.proteins._requestedCount[3863], self.proteins.total_counts()[3863],self.proteins.counts()[3863])

		# import ipdb; ipdb.set_trace() # ES2
		print("protein degradation evolveState ended")




	def _proteinDegRates(self):
		return self.rawDegRate * self.timeStepSec()
