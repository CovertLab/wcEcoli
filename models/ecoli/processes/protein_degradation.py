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

DYNAMIC_PD = False # todo: make this a simulation or ParCa parameter option
DYNAMIC_PD_TYPE = 1 # 1 is first order, 2 is optimized degradation with beta

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
		complexIds = np.array(sim_data.process.complexation.ids_complexes)
		self.complexes = self.bulkMoleculesView(complexIds)

		# Add the IDs for proteins and complexes to self:
		self.protein_IDs = proteinIds
		self.complex_IDs = complexIds

		# Obtain Avogadro's number for conversion of counts to molar concentration:
		self.n_avogadro = sim_data.constants.n_avogadro.asNumber(units.mol ** -1)

		self.bulkMoleculesRequestPriorityIs(REQUEST_PRIORITY_DEGRADATION)

	def calculateRequest(self):

		if DYNAMIC_PD == True:
			# Designate dynamic protein degradation model type
			self.dynamic_PD_usage = True

			# Determine how many proteins to degrade based on the degradation rates and counts of each protein
			nProteinsToDegrade = np.fmin(
				self.randomState.poisson(self._proteinDegRates() * self.proteins.total_counts()),
				self.proteins.total_counts()
			) # NOTE: this will be overwritten for some proteins below

			# Extract the counts_to_molar conversion factor:
			def get_counts_to_molar():
				cell_mass = self.readFromListener("Mass",
												  "cellMass") * units.fg  # if confused on how this works, look at 05.05.2025 notes
				# convert cell_mass to g:
				cell_mass_g = cell_mass.asNumber(units.g)  # now in units of g
				cell_mass_g = cell_mass_g  # todo: determine if this should be put back in: * units.g # physically puts a "[g]" unit on it
				cell_density = self.readFromListener("Mass",
													 "cellDensity")  # todo: figure out if needed (or if the constant value can be assumed here)
				cell_volume = cell_mass_g / cell_density  # L
				counts_to_molar = 1 / (
						self.n_avogadro * cell_volume)  # mol/L ? # todo: confirm units

				return counts_to_molar


			# Function for extracting the protease concentration:
			def get_protease_concentration(protease):
				counts_to_molar = get_counts_to_molar()
				protease_idx = np.where(self.complex_IDs == protease)[0][0]
				lon_complex_counts = self.complexes.total_counts()[
					protease_idx]  # todo: go back and figure out which lon complex ids are most important (total_counts, _totalCounts, etc.)
				protease_concentration = lon_complex_counts * counts_to_molar
				return protease_concentration

			# Function for extracting the substrate concentrations:
			def get_substrate_concentrations(substrates):
				counts_to_molar = get_counts_to_molar()
				concentrations = []
				for substrate in substrates:
					substrate_idx = np.where(self.protein_IDs == substrate)[0][0]
					substrate_concentration = self.proteins.total_counts()[
												  substrate_idx] * counts_to_molar
					concentrations.append(substrate_concentration) # mol/L [=] M
				return concentrations

			# Function for degrading proteins based on k_active values:
			def degrade_proteins(substrates, k_actives, nProteinsToDegrade):
				counts_to_molar = get_counts_to_molar()
				# Obtain the substrate concentrations:
				substrate_concentrations = get_substrate_concentrations(substrates)
				for i in range(len(substrate_concentrations)):
					# Find the protein index in the bulk view:
					substrate_idx = np.where(self.protein_IDs == substrates[i])[0][0]

					# Calculate the number of proteins degraded:
					proteins_degraded = k_actives[i] * substrate_concentrations[i]

					# Convert to counts and round:
					proteins_degraded = proteins_degraded / counts_to_molar
					proteins_degraded = round(proteins_degraded)

					# Degrade the proteins
					nProteinsToDegrade[
						substrate_idx] = proteins_degraded  # todo: check to see that this is working

				return nProteinsToDegrade

			# List the protease:
			protease = 'CPLX0-2881[c]'  # this is Lon

			# List the substrates (NOTE: PAY ATTENTION TO SUBSTRATE ORDER!):
			substrates = ['G6890-MONOMER[c]',
						  'PD03938[c]',
						  'G6737-MONOMER[c]',
						  'RPOD-MONOMER[c]',
						  'PD02936[c]',
						  'RED-THIOREDOXIN2-MONOMER[c]']

			if DYNAMIC_PD_TYPE == 1:
				# First order dynamic protein degradation model
				self.dynamic_PD_type = "first_order"
				self.dynamic_PD_substrates = substrates

				# Define parameters for dynamic degradation model:
				kcat = 0.071  # 1/s, https://jbioleng.biomedcentral.com/articles/10.1186/1754-1611-6-9#Sec29
				# these are from a first order approximation using the half life and average concentrations:
				Kms = [0.1867, 0.0692, 0.2517, 0.4397, 0.4095, 0.3202] # todo: add units
				# NOTE: pay attention to substrate order when inputting here too!
				self.dynamic_PD_Kms = Kms # todo: add units

				# Calculate the active degradation rate constant for each substrate:
				def calculate_first_order_k_active(protease, substrates, kcat, Kms):
					# Calulate the complex concentration:
					protease_concentration = get_protease_concentration(protease)

					# Obtain the substrate concentrations:
					substrate_concentrations = get_substrate_concentrations(substrates)

					# Calculate the k_active values:
					k_actives = []
					for i in range(len(substrate_concentrations)):
						k_active = (protease_concentration * kcat) / (Kms[i] + substrate_concentrations[i])
						k_actives.append(k_active)

					return k_actives

				# Call the functions:
				k_actives = calculate_first_order_k_active(protease, substrates, kcat, Kms)
				# reassign selected proteins in nProteinsToDegrade based on active degradation:
				nProteinsToDegrade = degrade_proteins(substrates, k_actives, nProteinsToDegrade)

			if DYNAMIC_PD_TYPE == 2:
				# Optimized dynamic protein degradation model
				self.dynamic_PD_type = "optimized_degradation_with_beta"
				self.dynamic_PD_substrates = substrates

				# Define parameters for dynamic degradation model:
				kcat = 0.071 # 1/s, https://jbioleng.biomedcentral.com/articles/10.1186/1754-1611-6-9#Sec29
				# from H.S.:
				Kms = [0.0025348 , 0.0001845, 0.00414169, 0.0062004 , 0.00176283, 0.00508851]
				# NOTE: pay attention to substrate order when inputting here too!
				self.dynamic_PD_Kms = Kms  # todo: add units

				# Calculate the active degradation rate constant for each substrate:
				def calculate_k_active(protease, substrates, kcat, Kms):
					# Calulate the complex concentration:
					protease_concentration = get_protease_concentration(protease)

					# Obtain the substrate concentrations:
					substrate_concentrations = get_substrate_concentrations(substrates)

					# Calculate the beta value:
					beta = 1
					for i in range(len(substrate_concentrations)):
						term = substrate_concentrations[i] / (Kms[i])
						beta += term

					# Calculate the k_active value:
					k_actives = []
					for i in range(len(substrate_concentrations)):
						k_active = (protease_concentration * kcat) / (Kms[i]*beta)
						k_actives.append(k_active)
					return k_actives

				# Call the functions:
				k_actives = calculate_k_active(protease, substrates, kcat, Kms)
				nProteinsToDegrade = degrade_proteins(substrates, k_actives, nProteinsToDegrade)

		else:
			# Use default model protein degradtion:
			self.dynamic_PD_usage = False
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


	def _proteinDegRates(self):
		return self.rawDegRate * self.timeStepSec()
