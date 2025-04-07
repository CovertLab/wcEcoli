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
		self.monomer_id_to_index = {
			monomer_id: i for i, monomer_id in enumerate(self.monomerIds)}

		# Create view on to active 70S ribosomes
		self.active_ribosomes = self.uniqueMoleculesView('active_ribosome')

		# Create views onto bulk 30S and 50S ribosomal subunits
		self.ribosome30S = self.bulkMoleculeView(sim_data.molecule_ids.s30_full_complex)
		self.ribosome50S = self.bulkMoleculeView(sim_data.molecule_ids.s50_full_complex)

		# Create view onto RNAs
		self.RNAs = self.uniqueMoleculesView('RNA')

		# TODO: clean up implementation if this proves useful
		self.new_genes_new_renomalization_method = sim_data.new_genes_new_renomalization_method
		self.new_gene_renormalization = sim_data.new_gene_renormalization

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

		if not self.new_genes_new_renomalization_method:
			# Calculate initiation probabilities for ribosomes based on mRNA counts
			# and associated mRNA translational efficiencies
			protein_init_prob = normalize(
				cistron_counts[self.cistron_to_monomer_mapping] * self.translationEfficiencies
			)
		# TODO: clean up implementation if this proves useful
		else:
			# Fixed proability for ribosome and RNAP proteins
			protein_init_prob_production_machinery_monomer_ids = [
				"EG10893-MONOMER[c]", "RPOC-MONOMER[c]", "RPOB-MONOMER[c]", "EG10864-MONOMER[c]",
				"EG10865-MONOMER[c]", "EG10866-MONOMER[c]", "EG10867-MONOMER[c]",
				"EG10868-MONOMER[c]", "EG10869-MONOMER[c]", "EG10870-MONOMER[c]",
				"EG10871-MONOMER[c]", "EG10872-MONOMER[c]", "EG10873-MONOMER[c]",
				"EG10874-MONOMER[c]", "EG10875-MONOMER[c]", "EG10876-MONOMER[c]",
				"EG10877-MONOMER[c]", "EG10878-MONOMER[c]", "EG10879-MONOMER[c]",
				"EG10880-MONOMER[c]", "EG10881-MONOMER[c]", "EG10882-MONOMER[c]",
				"EG10883-MONOMER[c]", "EG10884-MONOMER[c]", "EG10885-MONOMER[c]",
				"EG10886-MONOMER[c]", "EG10887-MONOMER[c]", "EG10888-MONOMER[c]",
				"EG10889-MONOMER[c]", "EG10890-MONOMER[c]", "EG10891-MONOMER[c]",
				"EG10892-MONOMER[c]", "EG10900-MONOMER[c]", "EG10901-MONOMER[c]",
				"EG10902-MONOMER[c]", "EG10903-MONOMER[c]", "EG10904-MONOMER[c]",
				"EG10905-MONOMER[c]", "EG10906-MONOMER[c]", "EG10907-MONOMER[c]",
				"EG10908-MONOMER[c]", "EG10909-MONOMER[c]", "EG10910-MONOMER[c]",
				"EG10911-MONOMER[c]", "EG10912-MONOMER[c]", "EG10913-MONOMER[c]",
				"EG10914-MONOMER[c]", "EG10915-MONOMER[c]", "EG10916-MONOMER[c]",
				"EG10917-MONOMER[c]", "EG10918-MONOMER[c]", "EG10919-MONOMER[c]",
				"EG10920-MONOMER[c]", "EG11231-MONOMER[c]", "EG11232-MONOMER[c]",
				"EG50001-MONOMER[c]", "EG50002-MONOMER[c]"
			]
			protein_init_prob_production_machinery_wt_averages = [ # from Sherlock sims Var 0 20250126 (192 cells)
				0.001568891222931201,0.0006442700785720917,0.0006821817327190883,
				0.00540970102753629,0.00475859642347963,0.005034707531523795,
				0.004920232468467931,0.005350714785314307,0.005140925789781231,
				0.007660048068656697,0.004636646477916895,0.0057325854993190775,
				0.01763480093354898,0.004320080554452435,0.005479543999452174,
				0.004917244347880968,0.004399804088259584,0.004488859197106797,
				0.005060966573133763,0.0046912427457965105,0.005435711110032523,
				0.004570079546941582,0.00481457550561655,0.005409008419909916,
				0.005494705370257568,0.007657301217291726,0.004337003245307404,
				0.004940956947911407,0.005207787503348046,0.009239694631968902,
				0.005609539275637401,0.005527568848369066,0.005110946831165395,
				0.005434301337342072,0.0045139089962265636,0.004871321917916269,
				0.005008789416411867,0.008853681017458502,0.00650624715296144,
				0.005202775094430449,0.005329672381466186,0.005097581750432442,
				0.004967106914824536,0.005605083106438874,0.005061679312008517,
				0.005257429158061467,0.009294907769563793,0.0051826962725286165,
				0.004307513801147991,0.010736563374561539,0.004618135225982999,
				0.004506030836185574,0.01242081688111554,0.004659388090561654,
				0.004924109782446105,0.005288877847889438,0.005461438868113937
			]

			sum_protein_production_machinery_init_probs = np.sum(np.array(protein_init_prob_production_machinery_wt_averages))
			remaining_protein_init_prob = 1 - sum_protein_production_machinery_init_probs

			# Find index for all ids in protein_init_prob_production_machinery_monomer_ids
			indices_with_fixed_protein_init_prob = []
			for protein_init_prob_production_machinery_monomer_id in protein_init_prob_production_machinery_monomer_ids:
				indices_with_fixed_protein_init_prob.append(
					self.monomer_id_to_index[protein_init_prob_production_machinery_monomer_id])
			if not self.new_gene_renormalization: # then new gene should also have fixed protein init prob
				# exit the sim
				# TODO: think about later
				raise ValueError("New gene fixed protein init prob would need to be provided")

			indices_with_fixed_protein_init_prob_mask = np.zeros(len(self.monomerIds), dtype=bool)
			indices_with_fixed_protein_init_prob_mask[indices_with_fixed_protein_init_prob] = True

			# Calculate initiation probabilities for ribosomes based on mRNA counts
			# and associated mRNA translational efficiencies
			protein_init_prob = normalize(
				cistron_counts[self.cistron_to_monomer_mapping] * self.translationEfficiencies
			)

			# Set fixed init probs for production machinery
			protein_init_prob[indices_with_fixed_protein_init_prob] = np.array(protein_init_prob_production_machinery_wt_averages)

			# Now renormalize the rest of the protein_init_probs
			indices_to_adjust_init_prob = np.sum(protein_init_prob[~indices_with_fixed_protein_init_prob_mask])
			protein_init_prob[~indices_with_fixed_protein_init_prob_mask] *= remaining_protein_init_prob / indices_to_adjust_init_prob

			assert np.isclose(np.sum(protein_init_prob), 1, atol=1e-5), f"Sum of protein init probs is not 1: {np.sum(protein_init_prob)}"

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

		# TODO: think about this renormalization with machinery upregulation

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
