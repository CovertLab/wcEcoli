"""
ProteinDegradation

Protein degradation sub-model. Encodes molecular simulation of protein degradation as a Poisson process

TODO:
- protein complexes
- add protease functionality
"""

""" USER INPUTS """
USE_LON_DEGRADATION = False # whether to use lon degradation rates or not
METHOD = 2
""" END USER INPUTS """
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

		# get the names of the proteins and complexes:
		self.protein_IDs = proteinIds
		self.complex_IDs = complexIds

		# get the molecular weights of the proteins and complexes:
		self.protein_mws = sim_data.process.translation.monomer_data['mw']
		self.n_avogadro = sim_data.constants.n_avogadro.asNumber(units.mol**-1)



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
		print(self.proteins._totalCount[2342], self.proteins.counts()[2342], nProteinsToDegrade[2342])
		self.writeToListener('MonomerCounts', 'protein_deg_CR1__totalCount', self.proteins._totalCount.copy())
		#import ipdb; ipdb.set_trace()
		#import ipdb; ipdb.set_trace() # CR1
		if USE_LON_DEGRADATION == True:
			# todo: (MIA) should I be using the whole cell's mass or the proteinMass? (which also can be extracted from the listener I believe)
			cell_mass = self.readFromListener("Mass", "cellMass") * units.fg # if confused on how this works, look at 05.05.2025 notes
			# convert cell_mass to g:
			cell_mass_g = cell_mass.asNumber(units.g) # now in units of g
			cell_mass_g = cell_mass_g # todo: determine if i should put this back in: * units.g # physically puts a "[g]" unit on it

			# load in the protein mass:
			protein_mass = self.readFromListener("Mass", "proteinMass") * units.fg # note this is much smaller than the cell mass

			# load in the cell volume:
			# Todo: (MIA) determine if this is the correct way to get the cell volume, becuase it is not working but it appears to be in the listener script?
			#cell_volume = self.readFromListener("Mass", "cellVolume") * units.L # this matches the simulation output units

			# read in the cell density to manually calculate the cell volume:
			cell_density = self.readFromListener("Mass", "cellDensity") # todo: figure out if I should add this bc it messes up the math later since the mol doesnt stay:* units.g / units.L # this is apparently constant? i.e. sim_data.constants.cell_density is a variable
			cell_volume = cell_mass_g / cell_density # L
			counts_to_molar = 1 / (self.n_avogadro * cell_volume)  #  mol/L ? # todo: are these the right units?


			print("nemA free monomers to be degraded pre-active degradation:",
				  self.proteins.counts()[2342])

			if METHOD == 1:
				# Degrade selected proteins
				# lon complex index: 297

				# todo: pick up with just simply getting the conc of the lon complex, then the conc of the proteins.
				lon_complex_counts = self.complexes.total_counts()[297] # todo: go back and figure out which lon complex ids are most important (total_counts, _totalCounts, etc.)


				# just do one protein for now 'G6890-MONOMER[c]'
				interest_protein_counts = self.proteins.total_counts()[2342]

				# based off the first matlab calculation, degrade using the fsolve answer (calculated with 6 proteins present):
				# k = [P]kcat/(km + [S])
				kcat = 0.071 # 1/s, https://jbioleng.biomedcentral.com/articles/10.1186/1754-1611-6-9#Sec29
				km = 0.0575 # calculated with fsolve in matlab based on kcat (and taking into account other 6 proteins)
				# todo: since this is per second, might need to convert to per timestep (based on how many are in a second)

				k_active = lon_complex_counts * kcat / (km + interest_protein_counts)
				proteins_degraded = k_active * interest_protein_counts
				print("proteins degraded pre-rounding:", proteins_degraded)
				proteins_degraded = round(proteins_degraded)
				print("nemA monomers to be degraded post-active degradation:",proteins_degraded)

				# reassign the counts of the proteins to be degraded:
				nProteinsToDegrade[2342] = proteins_degraded
				print("nemA nProteinsToDegrade[2342] post-active degradation:",nProteinsToDegrade[2342])

			if METHOD == 2:
				print("METHOD 2 selected")
				# todo: pick up with just simply getting the conc of the lon complex, then the conc of the proteins.


				# list of the interesting proteins:
				protease_FM = 'EG10542-MONOMER[c]'
				protease = 'CPLX0-2881[c]'
				protease_idx = np.where(self.complex_IDs == protease)[0][0]
				lon_complex_counts = self.complexes.total_counts()[protease_idx]  # todo: go back and figure out which lon complex ids are most important (total_counts, _totalCounts, etc.)
				protease_concentration = lon_complex_counts * counts_to_molar

				# determine the lon complex mass:
				# todo: find a way to do this with the matricies. for now, just use the lon complex counts

				print("lon complex concentration:", protease_concentration)

				# see what the total protein counts are:
				protein_counts_now = self.readFromListener("MonomerCounts", "monomerCounts")[527]
				print("total Lon protein counts: ", protein_counts_now)
				print("lon complex counts:", lon_complex_counts)
				print("lon free monomer counts: ", self.proteins.total_counts()[527])
				print("lon FMs + 6*lon complex counts: ", (self.proteins.total_counts()[527] + (lon_complex_counts*6)))




				substrates = ['G6890-MONOMER[c]',
					   'PD03938[c]',
					   'G6737-MONOMER[c]',
					   'RPOD-MONOMER[c]',
					   'PD02936[c]',
					   'RED-THIOREDOXIN2-MONOMER[c]']

				# substrate counts:
				s1 = self.proteins.total_counts()[2342]
				s2 = self.proteins.total_counts()[3863]
				s3 = self.proteins.total_counts()[2229]
				s4 = self.proteins.total_counts()[4002]
				s5 = self.proteins.total_counts()[3854]
				s6 = self.proteins.total_counts()[3977] # todo: double check this is correct

				def get_substrate_counts(substrates):
					counts = []
					for substrate in substrates:
						substrate_idx = np.where(self.protein_IDs == substrate)[0][0]
						substrate_count = self.proteins.total_counts()[substrate_idx]
						counts.append(substrate_count)
					return counts
				def get_substrate_mws(substrates):
					mws = []
					for substrate in substrates:
						substrate_idx = np.where(self.protein_IDs == substrate)[0][0]
						substrate_mw = self.protein_mws[substrate_idx]
						mws.append(substrate_mw)
					return mws

				def get_substrate_concentrations(substrates):
					concentrations = []
					for substrate in substrates:
						substrate_idx =  np.where(self.protein_IDs == substrate)[0][0]
						substrate_concentration = self.proteins.total_counts()[substrate_idx] * counts_to_molar
						concentrations.append(substrate_concentration)
					return concentrations

				# dont actually need these:
				#substrate_mws = get_substrate_mws(substrates) # g/mol

				# get the substrate concentrations:
				substrate_concentrations = get_substrate_concentrations(substrates)  # mol/L [=] M


				kcat = 0.071  # 1/s, https://jbioleng.biomedcentral.com/articles/10.1186/1754-1611-6-9#Sec29
				# todo: the guess Km value was 3.7x10-6 M
				Kms1 = [0.0575, 0.0444, 0.0599, 0.0630, 0.0627, 0.0614] # M --> source: original matlab file

				# Kms from new solving with b = 1.24

				# original order: PD03938, G6890-MONOMER, G6737-MONOMER, RPOD-MONOMER, PD02936, RED-THIOREDOXIN2-MONOMER --> SWITCHED PROPERLY TO SUBSTRATES ORDER
				Kms2 = [ 3.1190099e-4, 5.29487648e-6, 5.2999865e-4,   4.1087516e-4, 1.23777405e-4, 9.36398294e-3] # M --> source: original matlab file

				Kms = Kms2





				# calculate the k_active value:
				def calculate_k_active(protease, substrates, kcat, Kms):
					# calulate the complex concentration:
					protease_idx = np.where(self.complex_IDs == protease)[0][0]
					lon_complex_counts = self.complexes.total_counts()[protease_idx]  # todo: go back and figure out which lon complex ids are most important (total_counts, _totalCounts, etc.)
					protease_concentration = lon_complex_counts * counts_to_molar

					substrate_concentrations = []
					for substrate in substrates:
						substrate_idx =  np.where(self.protein_IDs == substrate)[0][0]
						substrate_concentration = self.proteins.total_counts()[substrate_idx] * counts_to_molar
						substrate_concentrations.append(substrate_concentration)

					# calculate the beta value:
					beta = 1
					for i in range(len(substrate_concentrations)):
						term = substrate_concentrations[i] / (Kms[i])
						beta += term
						print("term[{}]:".format(i), term)
					print("beta value:", beta)

					# calculate the k_active value:
					k_actives = []
					for i in range(len(substrate_concentrations)):
						k_active = (protease_concentration * kcat * substrate_concentrations[i]) / (Kms[i]*beta)
						k_actives.append(k_active)
					print("k_active values:", k_actives)
					import ipdb; ipdb.set_trace()
					return k_actives



				def degrade_proteins(substrates, substrate_concentrations, k_actives, nProteinsToDegrade):
					proteins_degraded = []
					for i in range(len(substrate_concentrations)):
						substrate_idx = np.where(self.protein_IDs == substrates[i])[0][0]
						print("substrate nProteinsToDegrade[{}] pre-active degradation:".format(
							substrates[i]), nProteinsToDegrade[substrate_idx])
						proteins_degraded = k_actives[i] * substrate_concentrations[i]
						# convert to counts and round:
						proteins_degraded = proteins_degraded / counts_to_molar
						#import ipdb; ipdb.set_trace()
						print("proteins degraded pre-rounding:", proteins_degraded)
						proteins_degraded = round(proteins_degraded)
						# degrade the proteins
						nProteinsToDegrade[substrate_idx] = proteins_degraded
						print("substrate nProteinsToDegrade[{}] post-active degradation:".format(substrates[i]), nProteinsToDegrade[substrate_idx])

					return nProteinsToDegrade

				# call the function:
				k_actives = calculate_k_active(protease, substrates, kcat, Kms)
				nProteinsToDegrade = degrade_proteins(substrates, substrate_concentrations, k_actives, nProteinsToDegrade)

				# confirm it works:
				print("post-active degradation counts for nemA: ", nProteinsToDegrade[2342])

			else:
				print("no method selected, proceeding with normal degradation")

		# check if the if worked:
		print("nemA nProteinsToDegrade[2342] post-if statement:", nProteinsToDegrade[2342])

		# Determine the number of hydrolysis reactions
		nReactions = np.dot(self.proteinLengths.asNumber(), nProteinsToDegrade)

		# Determine the amount of water required to degrade the selected proteins
		# Assuming one N-1 H2O is required per peptide chain length N
		self.h2o.requestIs(nReactions - np.sum(nProteinsToDegrade))
		self.proteins.requestIs(nProteinsToDegrade)
		print("self.proteins.counts()[2342] counts after if statement:", self.proteins.counts()[2342])
		print("self.proteins._totalCount[2342] counts after if statement:", self.proteins._totalCount[2342])
		print("nProteinsToDegrade[2342] counts after if statement:", nProteinsToDegrade[2342])

		#print(self.proteins._totalCount[3863], self.proteins.total()[3863], self.proteins._counts()[3863], self.proteins._requestedCount[3863], self.proteins.total_counts()[3863],self.proteins.counts()[3863])
		#print(self.proteins._totalCount[3863], self.proteins.total()[3863], self.proteins._counts()[3863], self.proteins._requestedCount[3863], self.proteins.total_counts()[3863],self.proteins.counts()[3863], nProteinsToDegrade[3863])
		#print(self.proteins._totalCount[2342], self.proteins.counts()[2342], nProteinsToDegrade[2342])

		self.writeToListener('MonomerCounts', 'protein_deg_CR2__totalCount', self.proteins._totalCount.copy())
		counts_for_CR2 = self.proteins.counts()
		self.writeToListener('MonomerCounts', 'protein_deg_CR2_counts', counts_for_CR2) # todo: ask riley why .copy() doesnt work here
		#import ipdb; ipdb.set_trace()

		#import ipdb; ipdb.set_trace() # CR2
		print("protein degradation calcREQUEST ended")




	def evolveState(self):
		print("protein degradation evolveState started")
		# todo: determine if this is where the proteins are actually degraded
		# todo: test implementation of the degradation rates:
		# only degrade proteins if True:
		# if USE_LON_DEGRADATION == True:
		# 	print("nemA free monomers degraded pre-active degradation:", self.proteins.counts()[2342])
		# 	# Degrade selected proteins
		# 	# lon complex index: 297
		# 	lon_complex_counts = self.complexes.total_counts()[297] # todo: go back and figure out which lon complex ids are most important (total_counts, _totalCounts, etc.)
		#
		# 	# just do one protein for now 'G6890-MONOMER[c]'
		# 	interest_protein_counts = self.proteins.total_counts()[2342]
		#
		# 	# based off the first matlab calculation, degrade using the fsolve answer (calculated with 6 proteins present):
		# 	# k = [P]kcat/(km + [S])
		# 	kcat = 0.071 # 1/s, https://jbioleng.biomedcentral.com/articles/10.1186/1754-1611-6-9#Sec29
		# 	km = 0.0575 # calculated with fsolve in matlab based on kcat (and taking into account other 6 proteins)
		#
		# 	k_active = lon_complex_counts * kcat / (km + interest_protein_counts)
		# 	proteins_degraded = k_active * interest_protein_counts
		# 	print("proteins degraded pre-round:", proteins_degraded)
		# 	proteins_degraded = round(proteins_degraded)
		# 	print("nemA free monomers degraded post-active degradation:",proteins_degraded)
		#
		# 	# reassign the counts of the proteins to be degraded:
		# 	self.proteins.requestIs('G6890-MONOMER[c]') = proteins_degraded



		# Degrade selected proteins, release amino acids from those proteins back into the cell, 
		# and consume H_2O that is required for the degradation process
		print("protein counts for 2342 after if statement:", self.proteins.counts()[2342])

		self.metabolites.countsInc(np.dot(
			self.proteinDegSMatrix,
			self.proteins.counts()
			))
		#print(self.proteins._totalCount[3863], self.proteins.total()[3863], self.proteins._counts()[3863], self.proteins._requestedCount[3863], self.proteins.total_counts()[3863],self.proteins.counts()[3863])
		print(self.proteins._totalCount[2342], self.proteins.counts()[2342])
		counts_for_ES1 = self.proteins.counts()
		self.writeToListener('MonomerCounts', 'protein_deg_ES1_counts', counts_for_ES1)





		#import ipdb; ipdb.set_trace() # ES1

		self.proteins.countsIs(0) # does this set the counts to zero?
		#print(self.proteins._totalCount[3863], self.proteins.total()[3863], self.proteins._counts()[3863], self.proteins._requestedCount[3863], self.proteins.total_counts()[3863],self.proteins.counts()[3863])

		# import ipdb; ipdb.set_trace() # ES2
		print("protein degradation evolveState ended")




	def _proteinDegRates(self):
		return self.rawDegRate * self.timeStepSec()
