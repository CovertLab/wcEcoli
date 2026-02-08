"""
SimulationData for translation process
"""

import numpy as np
from wholecell.sim.simulation import MAX_TIME_STEP
from wholecell.utils import data, units
from wholecell.utils.unit_struct_array import UnitStructArray
from wholecell.utils.polymerize import polymerize
from wholecell.utils.random import make_elongation_rates

PROCESS_MAX_TIME_STEP = 2.

class Translation(object):
	""" Translation """

	def __init__(self, raw_data, sim_data):
		self.max_time_step = min(MAX_TIME_STEP, PROCESS_MAX_TIME_STEP)
		self.next_aa_pad = 1  # Need an extra amino acid in sequences lengths to find next one

		self._build_monomer_data(raw_data, sim_data)
		self._build_translation(raw_data, sim_data)
		self._build_translation_efficiency(raw_data, sim_data)
		self._build_elongation_rates(raw_data, sim_data)

	def __getstate__(self):
		"""Return the state to pickle with translationSequences removed and
		only storing data from translationSequences with pad values stripped.
		"""

		state = data.dissoc_strict(self.__dict__, ('translation_sequences',))
		state['sequences'] = np.array([
			seq[seq != polymerize.PAD_VALUE]
			for seq in self.translation_sequences], dtype=object)
		state['sequence_shape'] = self.translation_sequences.shape
		return state

	def __setstate__(self, state):
		"""Restore translationSequences and remove processed versions of the data."""
		sequences = state.pop('sequences')
		sequence_shape = state.pop('sequence_shape')
		self.__dict__.update(state)

		self.translation_sequences = np.full(sequence_shape, polymerize.PAD_VALUE, dtype=np.int8)
		for i, seq in enumerate(sequences):
			self.translation_sequences[i, :len(seq)] = seq

	def _build_monomer_data(self, raw_data, sim_data):
		# Get set of all cistrons IDs with an associated gene and right and left
		# end positions
		rna_id_to_gene_id = {
			gene['rna_ids'][0]: gene['id'] for gene in raw_data.genes}
		gene_id_to_left_end_pos = {
			gene['id']: gene['left_end_pos'] for gene in raw_data.genes
			}
		gene_id_to_right_end_pos = {
			gene['id']: gene['right_end_pos'] for gene in raw_data.genes
			}

		all_mRNA_cistrons = {
			rna['id'] for rna in raw_data.rnas
			if rna['id'] in rna_id_to_gene_id
			    and gene_id_to_left_end_pos[rna_id_to_gene_id[rna['id']]] is not None
			    and gene_id_to_right_end_pos[rna_id_to_gene_id[rna['id']]] is not None
				and rna['type'] == 'mRNA'
			}

		# Get mappings from monomer IDs to cistron IDs
		# TODO (rjuene): Handle cistrons with two or more associated proteins
		monomer_id_to_cistron_id = {
			rna['monomer_ids'][0]: rna['id']
			for rna in raw_data.rnas
			if len(rna['monomer_ids']) > 0}

		# Select proteins with valid sequences and mappings to valid mRNAs
		all_proteins = []
		for protein in raw_data.proteins:
			if sim_data.getter.is_valid_molecule(protein['id']):
				try:
					rna_id = monomer_id_to_cistron_id[protein['id']]
				except KeyError:
					continue
				if rna_id in all_mRNA_cistrons:
					all_proteins.append(protein)

		# Get protein IDs with compartments
		protein_ids = [protein['id'] for protein in all_proteins]
		protein_compartments = sim_data.getter.get_compartments(protein_ids)
		assert all([len(loc) == 1 for loc in protein_compartments])
		protein_ids_with_compartments = [
			f'{protein_id}[{loc[0]}]' for (protein_id, loc)
			in zip(protein_ids, protein_compartments)
			]
		n_proteins = len(protein_ids)
		
		# Get cistron IDs associated to each monomer
		cistron_ids = [
			monomer_id_to_cistron_id[protein['id']] for protein in all_proteins]

		# Get lengths and amino acids counts of each protein
		protein_seqs = sim_data.getter.get_sequences(protein_ids)
		lengths = [len(seq) for seq in protein_seqs]
		aa_counts = [
			[seq.count(aa) for aa in sim_data.amino_acid_code_to_id_ordered.keys()]
			for seq in protein_seqs]
		n_amino_acids = len(sim_data.amino_acid_code_to_id_ordered)

		# Get molecular weights
		mws = sim_data.getter.get_masses(protein_ids).asNumber(units.g / units.mol)

		# Calculate degradation rates based on the N-rule (Tobias et al., 1991)
		deg_rate_units = 1 / units.s
		n_end_rule_deg_rates = {
			row['aa_code']: (np.log(2)/(row['half life'])).asNumber(deg_rate_units)
			for row in raw_data.protein_half_lives_n_end_rule}

		# Get degradation rates from measured data (Macklin et al., 2020)
		measured_deg_rates = {
			p['id']: (np.log(2) / p['half life']).asNumber(deg_rate_units)
			for p in raw_data.protein_half_lives_measured
			}

		# Get degradation rates from pulsed-SILAC data (Nagar et al., 2021)
		pulsed_silac_deg_rates = {
			p['id']: (np.log(2) / p['half_life']).asNumber(deg_rate_units)
			for p in raw_data.protein_half_lives_pulsed_silac
		}

		# Get degradation rates from carbon-limited data (Gupta et al., 2024)
		clim_deg_rates = {
			p['id']: (np.log(2) / p['half_life']).asNumber(deg_rate_units)
			for p in
			raw_data.protein_half_lives_Clim4_STD_ratio_threshold_2_keep_NaNs
		}

		# Generate the mapping from monomer IDs to common names:
		self.monomer_id_to_common_name_dict = (
			self.generate_monomer_ID_to_common_name_dict(raw_data))

		# Extract the protease degradation classification type and
		# protease degradation contributions to degradation for each protein
		# estimated from Gupta et al. 2024 data (note: this data is not used
		# functionally in the model):
		self.protease_dict = {
			p['id']: {'protease_assignment': p['protease_assignment'],
					  'ClpP_fraction': p['ClpP'],
					  'Lon_fraction': p['Lon'],
					  'HslV_fraction': p['HslV'],
					  'Unexplained_fraction': p['Unexplained']
					  }
			for p in raw_data.protease_assignments
		}

		# Initialize protein information arrays to be saved with monomer data:
		common_name = np.full(len(all_proteins), None)
		deg_rate = np.zeros(len(all_proteins))
		half_life_source_ID = np.full(len(all_proteins), None)
		self.protease_assignment = np.full(len(all_proteins), None)
		self.clpp_contribution = np.full(len(all_proteins), None)
		self.lon_contribution = np.full(len(all_proteins), None)
		self.hslv_contribution = np.full(len(all_proteins), None)
		self.unexplained_contribution = np.full(len(all_proteins), None)


		# Obtain the selected protein degradation rate combination from raw_data:
		selected_PDR_combination = raw_data.protein_degradation_combo_option
		# NOTE: the default option is listed as
		# DEFAULT_PROTEIN_DEGRADATION_COMBO in wholecell/utils/constants.py

		# Assign proteins to degradation rates based on the selected combination:
		if selected_PDR_combination == "PDR_combo_2020":
			# Uses measured rates from Macklin et al., 2020 first, followed by
			# the N-end rule from Tobias et al., 1991
			for i, protein in enumerate(all_proteins):
				common_name[i] = self.get_common_name(protein['id'])
				self.determine_protease_involvement(protein['id'], i)
				# Use measured degradation rates if available
				if protein['id'] in measured_deg_rates:
					deg_rate[i] = measured_deg_rates[protein['id']]
					half_life_source_ID[i] = 'CL_measured_deg_rates_2020'
				# If measured rates are unavailable, use N-end rule
				else:
					seq = protein['seq']
					assert seq[0] == 'M'  # All protein sequences should start with methionine
					# Set N-end residue as second amino acid if initial methionine
					# is cleaved
					n_end_residue = seq[protein['cleavage_of_initial_methionine']]
					deg_rate[i] = n_end_rule_deg_rates[n_end_residue]
					half_life_source_ID[i] = 'N_end_rule'

		if selected_PDR_combination == "PDR_combo_2022":
			# Uses measured rates from Macklin et al., 2020 first, followed by
			# pulsed silac rates from Nagar et al., 2021, and finally N-end rule
			# from Tobias et al., 1991
			for i, protein in enumerate(all_proteins):
				common_name[i] = self.get_common_name(protein['id'])
				self.determine_protease_involvement(protein['id'], i)
				# Use measured degradation rates if available
				if protein['id'] in measured_deg_rates:
					deg_rate[i] = measured_deg_rates[protein['id']]
					half_life_source_ID[i] = 'CL_measured_deg_rates_2020'
				elif protein['id'] in pulsed_silac_deg_rates:
					deg_rate[i] = pulsed_silac_deg_rates[protein['id']]
					half_life_source_ID[i] = 'Nagar_et_al_ML_2021'
				# If measured rates are unavailable, use N-end rule
				else:
					seq = protein['seq']
					assert seq[0] == 'M'  # All protein sequences should start
					# with methionine
					# Set N-end residue as second amino acid if initial methionine
					# is cleaved
					n_end_residue = seq[protein['cleavage_of_initial_methionine']]
					deg_rate[i] = n_end_rule_deg_rates[n_end_residue]
					half_life_source_ID[i] = 'N_end_rule'

		if selected_PDR_combination == "PDR_combo_2025":
			# Uses measured rates from Macklin et al., 2020 first, followed by
			# Carbon limited rates from Gupta et al., 2024, and finally N-end
			# rule from Tobias et al., 1991
			for i, protein in enumerate(all_proteins):
				common_name[i] = self.get_common_name(protein['id'])
				self.determine_protease_involvement(protein['id'], i)
				# Use measured degradation rates if available
				if protein['id'] in measured_deg_rates:
					deg_rate[i] = measured_deg_rates[protein['id']]
					half_life_source_ID[i] = 'CL_measured_deg_rates_2020'
				elif protein['id'] in clim_deg_rates:
					deg_rate[i] = clim_deg_rates[protein['id']]
					half_life_source_ID[i] = 'Gupta_et_al_MS_2024'
				# If measured rates are unavailable, use N-end rule
				else:
					seq = protein['seq']
					assert seq[0] == 'M'  # All protein sequences should start
					# with methionine
					# Set N-end residue as second amino acid if initial methionine
					# is cleaved
					n_end_residue = seq[protein['cleavage_of_initial_methionine']]
					deg_rate[i] = n_end_rule_deg_rates[n_end_residue]
					half_life_source_ID[i] = 'N_end_rule'


		max_protein_id_length = max(
			len(protein_id) for protein_id in protein_ids_with_compartments)
		max_cistron_id_length = max(
			len(cistron_id) for cistron_id in cistron_ids)
		max_HL_source_ID_length = max(
			len(source_ID) for source_ID in half_life_source_ID)
		max_protease_ID_length = max(
			len(protease_ID) for protease_ID in self.protease_assignment
			if protease_ID is not None)
		max_common_name_length = max(
			len(name) for name in common_name if name is not None)

		monomer_data = np.zeros(
			n_proteins,
			dtype = [
				('id', 'U{}'.format(max_protein_id_length)),
				('cistron_id', 'U{}'.format(max_cistron_id_length)),
				('common_name', 'U{}'.format(max_common_name_length)),
				('degradation_rate', 'f8'),
				('half_life_source', 'U{}'.format(max_HL_source_ID_length)),
				('protease_assignment', 'U{}'.format(max_protease_ID_length)),
				('ClpP_contribution_fraction', 'f8'),
				('Lon_contribution_fraction', 'f8'),
				('HslV_contribution_fraction', 'f8'),
				('Unexplained_contribution_fraction', 'f8'),
				('length', 'i8'),
				('aa_counts', '{}i8'.format(n_amino_acids)),
				('mw', 'f8'),
				]
			)

		monomer_data['id'] = protein_ids_with_compartments
		monomer_data['cistron_id'] = cistron_ids
		monomer_data['common_name'] = common_name
		monomer_data['degradation_rate'] = deg_rate
		monomer_data['half_life_source'] = half_life_source_ID
		monomer_data['protease_assignment'] = self.protease_assignment
		monomer_data['ClpP_contribution_fraction'] = self.clpp_contribution
		monomer_data['Lon_contribution_fraction'] = self.lon_contribution
		monomer_data['HslV_contribution_fraction'] = self.hslv_contribution
		monomer_data['Unexplained_contribution_fraction'] = self.unexplained_contribution
		monomer_data['length'] = lengths
		monomer_data['aa_counts'] = aa_counts
		monomer_data['mw'] = mws

		field_units = {
			'id': None,
			'cistron_id': None,
			'common_name': None,
			'degradation_rate': deg_rate_units,
			'half_life_source': None,
			'protease_assignment': None,
			'ClpP_contribution_fraction': None,
			'Lon_contribution_fraction': None,
			'HslV_contribution_fraction': None,
			'Unexplained_contribution_fraction': None,
			'length': units.aa,
			'aa_counts': units.aa,
			'mw': units.g / units.mol,
			}

		self.monomer_data = UnitStructArray(monomer_data, field_units)
		self.n_monomers = len(self.monomer_data)

	def _build_translation(self, raw_data, sim_data):
		sequences = sim_data.getter.get_sequences(
			[protein_id[:-3] for protein_id in self.monomer_data['id']])

		max_len = np.int64(
			self.monomer_data["length"].asNumber().max()
			+ self.max_time_step * sim_data.constants.ribosome_elongation_rate_max.asNumber(units.aa / units.s)
			) + self.next_aa_pad

		self.translation_sequences = np.full((len(sequences), max_len), polymerize.PAD_VALUE, dtype=np.int8)
		aa_ids_single_letter = sim_data.amino_acid_code_to_id_ordered.keys()
		aaMapping = {aa: i for i, aa in enumerate(aa_ids_single_letter)}
		for i, sequence in enumerate(sequences):
			for j, letter in enumerate(sequence):
				self.translation_sequences[i, j] = aaMapping[letter]

		aaIDs = list(sim_data.amino_acid_code_to_id_ordered.values())

		self.translation_monomer_weights = (
			(
				sim_data.getter.get_masses(aaIDs)
				- sim_data.getter.get_masses([sim_data.molecule_ids.water])
				)
			/ sim_data.constants.n_avogadro
			).asNumber(units.fg)
		self.translation_end_weight = (sim_data.getter.get_masses([sim_data.molecule_ids.water]) / sim_data.constants.n_avogadro).asNumber(units.fg)

		# Load active ribosome footprint on RNA
		molecule_id_to_footprint_sizes = {
			row['molecule_id']: row['footprint_size']
			for row in raw_data.footprint_sizes}
		try:
			self.active_ribosome_footprint_size = \
				molecule_id_to_footprint_sizes['active_ribosome']
		except KeyError:
			raise ValueError(
				'RNA footprint size for ribosomes not found.')

	def _build_translation_efficiency(self, raw_data, sim_data):
		monomer_ids = [
			protein_id[:-3] for protein_id in self.monomer_data['id']]

		# Get mappings from monomer IDs to gene IDs
		monomer_id_to_rna_id = {
			rna['monomer_ids'][0]: rna['id']
			for rna in raw_data.rnas
			if len(rna['monomer_ids']) > 0}
		rna_id_to_gene_id = {
			gene['rna_ids'][0]: gene['id'] for gene in raw_data.genes}
		monomer_id_to_gene_id = {
			monomer_id: rna_id_to_gene_id[monomer_id_to_rna_id[monomer_id]]
			for monomer_id in monomer_ids}

		# Get mappings from gene IDs to translation efficiencies
		gene_id_to_trl_eff = {
			x["geneId"]: x["translationEfficiency"]
			for x in raw_data.translation_efficiency
			if type(x["translationEfficiency"]) == float}

		trl_effs = []
		for monomer_id in monomer_ids:
			gene_id = monomer_id_to_gene_id[monomer_id]

			if gene_id in gene_id_to_trl_eff:
				trl_effs.append(gene_id_to_trl_eff[gene_id])
			else:
				trl_effs.append(np.nan)

		# If efficiency is unavailable, the average of existing effciencies
		# is used
		self.translation_efficiencies_by_monomer = np.array(trl_effs)
		self.translation_efficiencies_by_monomer[
			np.isnan(self.translation_efficiencies_by_monomer)
			] = np.nanmean(self.translation_efficiencies_by_monomer)

	def _build_elongation_rates(self, raw_data, sim_data):
		protein_ids = self.monomer_data['id']
		ribosomal_protein_ids = sim_data.molecule_groups.ribosomal_proteins

		protein_indexes = {
			protein: index
			for index, protein in enumerate(protein_ids)}

		ribosomal_proteins = {
			rprotein: protein_indexes.get(rprotein, -1)
			for rprotein in ribosomal_protein_ids}

		self.ribosomal_protein_indexes = np.array([
			index
			for index in sorted(ribosomal_proteins.values())
			if index >= 0], dtype=np.int64)

		self.basal_elongation_rate = sim_data.constants.ribosome_elongation_rate_basal.asNumber(units.aa / units.s)
		self.max_elongation_rate = sim_data.constants.ribosome_elongation_rate_max.asNumber(units.aa / units.s)
		self.elongation_rates = np.full(
			self.n_monomers,
			self.basal_elongation_rate,
			dtype=np.int64)

		self.elongation_rates[self.ribosomal_protein_indexes] = self.max_elongation_rate

	# TODO (mia): when common names are added ot the proteins.tsv file, replace rnas.tsv here with that.
	# Function that generates a mapping from monomer IDs to common names:
	def generate_monomer_ID_to_common_name_dict(self, raw_data):
		"""
        Extracts the common name listed in rnas.tsv for each monomer ID.
        Returns:
        	a dictionary mapping monomer IDs to common names.
        """
		# TODO (mia): when proteins.tsv has common names added to it, replace
		#  raw_data.rnas with raw_data.proteins here so each protein maps
		#  to one default common name.
		rnas_info = raw_data.rnas
		monomer_IDs_to_common_names = {}
		for i in range(len(rnas_info)):
			rna_info = rnas_info[i]
			common_name = rna_info['common_name']
			monomer_IDs = rna_info['monomer_ids']
			if len(monomer_IDs) > 0:
				# Some rnas have multiple associated monomer IDs, so assign the
				# same common name to each of those monomer IDs for now:
				for monomer_ID in monomer_IDs:
					monomer_IDs_to_common_names[monomer_ID] = common_name

		return monomer_IDs_to_common_names

	def get_common_name(self, protein_ID):
		"""
        Retreive the assigned common name/gene symbol listed in rnas.tsv for a protein.
        Args:
            protein_ID: the ID of the protein
        Returns:
            common_name: The common name of the protein as specified in rnas.tsv.
        """
		# First remove the compartment tag if the input protein has one:
		if '[' in protein_ID:
			protein = protein_ID[:-3]  # subtract the compartment tag
		else:
			protein = protein_ID

		# Next, determine the common name (if there is one) for the protein:
		if protein not in self.monomer_id_to_common_name_dict.keys():
			# If the protein ID is not found in the mapping,
			# assign the common name to  None:
			common_name = None
		else:
			# If the protein ID is found in the mapping, return the common name
			common_name = self.monomer_id_to_common_name_dict[protein]
		return common_name

	# Function that will map proteins to the protease assignment if there is one:
	def determine_protease_involvement(self, protein_ID, index):
		"""
        Maps a protein to its protease assignment and estimated fractional
        contributions to degradation done by each protease.
        Args:
            protein_ID: ID of the protein to be mapped
            index: index of the protein being evaluated in the all_proteins list

        Returns: An update to the protease_assignment,
        ClpP_contribution, Lon_contribution, HslV_contribution,
        and Unexplained_contribution arrays in place within protease_dict.
        """
		if protein_ID in self.protease_dict.keys():
			self.protease_assignment[index] = self.protease_dict[protein_ID]['protease_assignment']
			self.clpp_contribution[index] = self.protease_dict[protein_ID]['ClpP_fraction']
			self.lon_contribution[index] = self.protease_dict[protein_ID]['Lon_fraction']
			self.hslv_contribution[index] = self.protease_dict[protein_ID]['HslV_fraction']
			self.unexplained_contribution[index] = self.protease_dict[protein_ID][
				'Unexplained_fraction']
		else:
			# If the protein does not have a degradation classification,
			# retain the "None" values from initialization:
			pass

	def make_elongation_rates(
			self,
			random,
			base,
			time_step,
			variable_elongation=False):

		return make_elongation_rates(
			random,
			self.n_monomers,
			base,
			self.ribosomal_protein_indexes,
			self.max_elongation_rate,
			time_step,
			variable_elongation)
