"""
PolypeptideInitiation

Polypeptide initiation sub-model.
"""

import numpy as np

import wholecell.processes.process
from wholecell.utils import units
from wholecell.utils.fitting import normalize


class PolypeptideInitiation(wholecell.processes.process.Process):
	""" PolypeptideInitiation """

	_name = "PolypeptideInitiation"

	def __init__(self):
		super(PolypeptideInitiation, self).__init__()

	def initialize(self, sim, sim_data):
		super(PolypeptideInitiation, self).initialize(sim, sim_data)

		# Load parameters
		self.proteinLengths = sim_data.process.translation.monomer_data["length"].asNumber()
		self.translationEfficiencies = normalize(
			sim_data.process.translation.translation_efficiencies_by_monomer)
		self.fracActiveRibosomeDict = sim_data.process.translation.ribosomeFractionActiveDict
		self.ribosomeElongationRateDict = sim_data.process.translation.ribosomeElongationRateDict
		self.variable_elongation = sim._variable_elongation_translation
		self.make_elongation_rates = sim_data.process.translation.make_elongation_rates
		self.rna_id_to_cistron_indexes = sim_data.process.transcription.rna_id_to_cistron_indexes
		self.cistron_start_end_pos_in_tu = sim_data.process.transcription.cistron_start_end_pos_in_tu
		self.tu_ids = sim_data.process.transcription.rna_data['id']
		self.n_TUs = len(self.tu_ids)
		self.active_ribosome_footprint_size = \
			sim_data.process.translation.active_ribosome_footprint_size

		# Get mapping from cistrons to protein monomers and TUs
		self.cistron_to_monomer_mapping = sim_data.relation.cistron_to_monomer_mapping
		self.cistron_tu_mapping_matrix = sim_data.process.transcription.cistron_tu_mapping_matrix
		self.monomer_index_to_cistron_index = {
			i: sim_data.process.transcription._cistron_id_to_index[monomer['cistron_id']]
			for (i, monomer) in enumerate(sim_data.process.translation.monomer_data)
			}
		self.monomer_index_to_tu_indexes = sim_data.relation.monomer_index_to_tu_indexes
		self.monomerIds = sim_data.process.translation.monomer_data[
			'id']

		# Create view on to active 70S ribosomes
		self.active_ribosomes = self.uniqueMoleculesView('active_ribosome')

		# Create views onto bulk 30S and 50S ribosomal subunits
		self.ribosome30S = self.bulkMoleculeView(sim_data.molecule_ids.s30_full_complex)
		self.ribosome50S = self.bulkMoleculeView(sim_data.molecule_ids.s50_full_complex)

		# Create view onto RNAs
		self.RNAs = self.uniqueMoleculesView('RNA')

	def calculateRequest(self):
		current_media_id = self._external_states['Environment'].current_media_id

		self.ribosome30S.requestAll()
		self.ribosome50S.requestAll()

		self.fracActiveRibosome = self.fracActiveRibosomeDict[current_media_id]

		# Read ribosome elongation rate from last timestep
		self.ribosomeElongationRate = self.readFromListener(
			"RibosomeData", "effectiveElongationRate")
		# If the ribosome elongation rate is zero (which is always the case for
		# the first timestep), set ribosome elongation rate to the one in
		# dictionary
		if self.ribosomeElongationRate == 0:
			self.ribosomeElongationRate = self.ribosomeElongationRateDict[
				current_media_id].asNumber(units.aa / units.s)
		self.elongation_rates = self.make_elongation_rates(
			self.randomState,
			self.ribosomeElongationRate,
			1,  # want elongation rate, not lengths adjusted for time step
			self.variable_elongation)

		# Ensure rates are never zero
		self.elongation_rates = np.fmax(self.elongation_rates, 1)

	def evolveState(self):
		# Calculate number of ribosomes that could potentially be initialized
		# based on counts of free 30S and 50S subunits
		inactiveRibosomeCount = np.min([
			self.ribosome30S.count().sum(),
			self.ribosome50S.count().sum(),
			])

		# Calculate actual number of ribosomes that should be activated based on
		# probabilities
		# Get attributes of active (translatable) mRNAs
		TU_index_RNAs, transcript_lengths, can_translate, is_full_transcript, unique_index_RNAs = self.RNAs.attrs(
			'TU_index', 'transcript_length', 'can_translate', 'is_full_transcript', 'unique_index')
		TU_index_mRNAs = TU_index_RNAs[can_translate]
		length_mRNAs = transcript_lengths[can_translate]
		unique_index_mRNAs = unique_index_RNAs[can_translate]
		is_full_transcript_mRNAs = is_full_transcript[can_translate]
		is_incomplete_transcript_mRNAs = np.logical_not(is_full_transcript_mRNAs)

		# Calculate counts of each mRNA cistron from fully transcribed
		# transcription units
		TU_index_full_mRNAs = TU_index_mRNAs[is_full_transcript_mRNAs]
		TU_counts_full_mRNAs = np.bincount(
			TU_index_full_mRNAs, minlength=self.n_TUs)
		cistron_counts = self.cistron_tu_mapping_matrix.dot(
			TU_counts_full_mRNAs)

		# Calculate counts of each mRNA cistron from partially transcribed
		# transcription units
		TU_index_incomplete_mRNAs = TU_index_mRNAs[is_incomplete_transcript_mRNAs]
		length_incomplete_mRNAs = length_mRNAs[is_incomplete_transcript_mRNAs]

		for (TU_index, length) in zip(TU_index_incomplete_mRNAs, length_incomplete_mRNAs):
			cistron_indexes = self.rna_id_to_cistron_indexes(self.tu_ids[TU_index])
			cistron_start_positions = np.array([
				self.cistron_start_end_pos_in_tu[(cistron_index, TU_index)][0]
				for cistron_index in cistron_indexes
				])

			cistron_counts[cistron_indexes] += length > cistron_start_positions

		# Calculate initiation probabilities for ribosomes based on mRNA counts
		# and associated mRNA translational efficiencies
		protein_init_prob = normalize(
			cistron_counts[self.cistron_to_monomer_mapping] * self.translationEfficiencies
		)
		target_protein_init_prob = protein_init_prob.copy()
		self.writeToListener(
			"RibosomeData", "target_prob_translation_per_transcript",
			target_protein_init_prob)

		# Calculate actual number of ribosomes that should be activated based
		# on probabilities
		self.activationProb = self._calculateActivationProb(
			self.fracActiveRibosome,
			self.proteinLengths,
			self.elongation_rates,
			target_protein_init_prob,
			self.timeStepSec())

		n_ribosomes_to_activate = np.int64(self.activationProb * inactiveRibosomeCount)

		if n_ribosomes_to_activate == 0:
			return

		# Cap the initiation probabilities at the maximum level physically
		# allowed from the known ribosome footprint sizes based on the
		# number of mRNAs
		max_p = (self.ribosomeElongationRate / self.active_ribosome_footprint_size
			* (units.s) * self.timeStepSec() / 
			n_ribosomes_to_activate).asNumber()
		max_p_per_protein = max_p*cistron_counts[self.cistron_to_monomer_mapping]
		self.writeToListener("RibosomeData", "max_p", max_p)
		self.writeToListener(
			"RibosomeData", "max_p_per_protein", max_p_per_protein)
		is_overcrowded = (protein_init_prob > max_p_per_protein)

		# Initalize flag to record if the number of ribosomes activated at this
		# time step needed to be reduced to prevent overcrowding
		is_n_ribosomes_to_activate_reduced = False

		# If needed, resolve overcrowding
		while np.any(protein_init_prob > max_p_per_protein):
			if protein_init_prob[~is_overcrowded].sum() != 0:
				# Resolve overcrowding through rescaling (preferred)
				protein_init_prob[is_overcrowded] = max_p_per_protein[
					is_overcrowded]
				scale_the_rest_by = (
					(1. - protein_init_prob[is_overcrowded].sum())
					/ protein_init_prob[~is_overcrowded].sum())
				protein_init_prob[~is_overcrowded] *= scale_the_rest_by
				is_overcrowded |= (protein_init_prob > max_p_per_protein)
			else:
				# If we cannot resolve the overcrowding through rescaling,
				# we need to activate fewer ribosomes. Set the number of
				# ribosomes to activate so that there will be no overcrowding.
				is_n_ribosomes_to_activate_reduced = True
				max_index = np.argmax(
					protein_init_prob[is_overcrowded] / max_p_per_protein[is_overcrowded])
				max_init_prob = protein_init_prob[is_overcrowded][max_index]
				associated_cistron_counts = cistron_counts[
					self.cistron_to_monomer_mapping][is_overcrowded][max_index]
				n_ribosomes_to_activate = np.int64((
					self.ribosomeElongationRate
					/ self.active_ribosome_footprint_size
					* (units.s) * self.timeStepSec() / max_init_prob
					* associated_cistron_counts).asNumber())

				# Update maximum probabilities based on new number of activated
				# ribosomes.
				max_p = (
					self.ribosomeElongationRate
					/ self.active_ribosome_footprint_size
					* (units.s) * self.timeStepSec()
					/ n_ribosomes_to_activate).asNumber()
				max_p_per_protein = max_p * cistron_counts[self.cistron_to_monomer_mapping]
				is_overcrowded = (protein_init_prob > max_p_per_protein)
				assert is_overcrowded.sum() == 0 # We expect no overcrowding

		# Compute actual transcription probabilities of each transcript
		actual_protein_init_prob = protein_init_prob.copy()
		self.writeToListener(
			"RibosomeData", "actual_prob_translation_per_transcript",
			actual_protein_init_prob)
		self.writeToListener(
			"RibosomeData", "mRNA_is_overcrowded", is_overcrowded)
		self.writeToListener(
			"RibosomeData", "is_n_ribosomes_to_activate_reduced",
			is_n_ribosomes_to_activate_reduced)

		# Sample multinomial distribution to determine which mRNAs have full
		# 70S ribosomes initialized on them
		n_new_proteins = self.randomState.multinomial(
			n_ribosomes_to_activate,
			protein_init_prob
		)

		protein_indexes = np.empty(n_ribosomes_to_activate, np.int64)
		mRNA_indexes = np.empty(n_ribosomes_to_activate, np.int64)
		positions_on_mRNA = np.empty(n_ribosomes_to_activate, np.int64)
		nonzeroCount = (n_new_proteins > 0)
		start_index = 0

		for protein_index, counts in zip(
				np.arange(n_new_proteins.size)[nonzeroCount],
				n_new_proteins[nonzeroCount]):
			# Set protein index
			protein_indexes[start_index:start_index + counts] = protein_index

			cistron_index = self.monomer_index_to_cistron_index[protein_index]

			attribute_indexes = []
			cistron_start_positions = []

			for TU_index in self.monomer_index_to_tu_indexes[protein_index]:
				attribute_indexes_this_TU = np.where(TU_index_mRNAs == TU_index)[0]
				cistron_start_position = self.cistron_start_end_pos_in_tu[
					(cistron_index, TU_index)][0]

				is_transcript_long_enough = length_mRNAs[attribute_indexes_this_TU] >= cistron_start_position

				attribute_indexes.extend(attribute_indexes_this_TU[is_transcript_long_enough])
				cistron_start_positions.extend(
					[cistron_start_position] * len(attribute_indexes_this_TU[is_transcript_long_enough]))

			n_mRNAs = len(attribute_indexes)

			# Distribute ribosomes among these mRNAs
			n_ribosomes_per_RNA = self.randomState.multinomial(
				counts, np.full(n_mRNAs, 1. / n_mRNAs))

			# Get unique indexes of each mRNA
			mRNA_indexes[start_index:start_index + counts] = np.repeat(
				unique_index_mRNAs[attribute_indexes], n_ribosomes_per_RNA)

			positions_on_mRNA[start_index:start_index + counts] = np.repeat(
				cistron_start_positions, n_ribosomes_per_RNA
				)

			start_index += counts

		# Create active 70S ribosomes and assign their attributes
		self.active_ribosomes.moleculesNew(
			n_ribosomes_to_activate,
			protein_index=protein_indexes,
			peptide_length=np.zeros(n_ribosomes_to_activate, dtype=np.int64),
			mRNA_index=mRNA_indexes,
			pos_on_mRNA=positions_on_mRNA,
		)

		# Decrement free 30S and 50S ribosomal subunit counts
		self.ribosome30S.countDec(n_new_proteins.sum())
		self.ribosome50S.countDec(n_new_proteins.sum())

		# Write number of initialized ribosomes to listener
		self.writeToListener("RibosomeData", "didInitialize", n_new_proteins.sum())
		self.writeToListener("RibosomeData", "ribosome_init_event_per_monomer",
							 n_new_proteins)

	def _calculateActivationProb(
			self, fracActiveRibosome, proteinLengths, ribosomeElongationRates,
			protein_init_prob, timeStepSec):
		"""
		Calculates the expected ribosome termination rate based on the ribosome
		elongation rate
		Params:
			- allTranslationTimes: Vector of times required to translate each
			protein
			- allTranslationTimestepCounts: Vector of numbers of timesteps
			required to translate each protein
			- averageTranslationTimeStepCounts: Average number of timesteps
			required to translate a protein, weighted by initiation
			probabilities
			- expectedTerminationRate: Average number of terminations in one
			timestep for one protein
		"""
		allTranslationTimes = 1. / ribosomeElongationRates * proteinLengths
		allTranslationTimestepCounts = np.ceil(allTranslationTimes / timeStepSec)
		averageTranslationTimestepCounts = np.dot(allTranslationTimestepCounts, protein_init_prob)
		expectedTerminationRate = 1.0 / averageTranslationTimestepCounts

		# Modify given fraction of active ribosomes to take into account early
		# terminations in between timesteps
		# allFractionTimeInactive: Vector of probabilities an "active" ribosome
		# 	will in effect be "inactive" because it has terminated during a
		# 	timestep
		# averageFractionTimeInactive: Average probability of an "active"
		# 	ribosome being in effect "inactive", weighted by initiation
		#	probabilities
		# effectiveFracActiveRnap: New higher "goal" for fraction of active
		# 	ribosomes, considering that the "effective" fraction is lower than
		# 	what the listener sees
		allFractionTimeInactive = 1 - allTranslationTimes / timeStepSec / allTranslationTimestepCounts
		averageFractionTimeInactive = np.dot(allFractionTimeInactive, protein_init_prob)
		effectiveFracActiveRibosome = fracActiveRibosome * 1 / (1 - averageFractionTimeInactive)

		# Return activation probability that will balance out the expected
		# termination rate
		activationProb = effectiveFracActiveRibosome * expectedTerminationRate / (1 - effectiveFracActiveRibosome)

		# The upper bound for the activation probability is temporarily set to
		# 1.0 to prevent negative molecule counts. This will lower the fraction
		# of active ribosomes for timesteps longer than roughly 1.8s.
		if activationProb >= 1.0:
			activationProb = 1

		return activationProb
