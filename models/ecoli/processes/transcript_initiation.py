"""
TranscriptInitiation

Transcription initiation sub-model.

TODO:
- use transcription units instead of single genes
- match sigma factors to promoters
"""

from typing import cast

import numpy as np
import scipy.sparse

import wholecell.processes.process
from wholecell.utils import units

class TranscriptInitiation(wholecell.processes.process.Process):
	""" TranscriptInitiation """

	_name = "TranscriptInitiation"


	# Constructor
	def __init__(self):
		super(TranscriptInitiation, self).__init__()

	def initialize(self, sim, sim_data):
		super(TranscriptInitiation, self).initialize(sim, sim_data)

		# Sim options
		self.trna_attenuation = sim._trna_attenuation
		self.ppgpp_regulation = sim._ppgpp_regulation

		# Load parameters
		self.fracActiveRnapDict = sim_data.process.transcription.rnapFractionActiveDict
		self.rnaLengths = sim_data.process.transcription.rna_data["length"]
		self.rnaPolymeraseElongationRateDict = sim_data.process.transcription.rnaPolymeraseElongationRateDict
		self.variable_elongation = sim._variable_elongation_transcription
		self.make_elongation_rates = sim_data.process.transcription.make_elongation_rates
		self.active_rnap_footprint_size = sim_data.process.transcription.active_rnap_footprint_size

		# Initialize matrices used to calculate synthesis probabilities
		self.basal_prob = sim_data.process.transcription_regulation.basal_prob.copy()
		if self.trna_attenuation:
			self.attenuated_rna_indices = sim_data.process.transcription.attenuated_rna_indices
			self.attenuation_adjustments = sim_data.process.transcription.attenuation_basal_prob_adjustments
			self.basal_prob[self.attenuated_rna_indices] += self.attenuation_adjustments
		self.n_TUs = len(self.basal_prob)
		self.delta_prob_matrix = sim_data.process.transcription_regulation.get_delta_prob_matrix(
			dense=True, ppgpp=self.ppgpp_regulation)

		# Determine changes from genetic perturbations
		self.genetic_perturbations = {}
		perturbations = getattr(sim_data, "genetic_perturbations", {})

		if len(perturbations) > 0:
			probability_indexes = [
				(index, sim_data.genetic_perturbations[rna_data['id']])
					for index, rna_data in enumerate(sim_data.process.transcription.rna_data)
					if rna_data['id'] in sim_data.genetic_perturbations]

			self.genetic_perturbations = {
				'fixedRnaIdxs': [pair[0] for pair in probability_indexes],
				'fixedSynthProbs': [pair[1] for pair in probability_indexes]
				}

		# Views
		self.active_RNAPs = self.uniqueMoleculesView('active_RNAP')
		self.inactive_RNAPs = self.bulkMoleculeView("APORNAP-CPLX[c]")
		self.full_chromosomes = self.uniqueMoleculesView('full_chromosome')
		self.promoters = self.uniqueMoleculesView('promoter')
		self.RNAs = self.uniqueMoleculesView('RNA')

		# ID Groups
		self.idx_rRNA = np.where(sim_data.process.transcription.rna_data['is_rRNA'])[0]
		self.idx_mRNA = np.where(sim_data.process.transcription.rna_data['is_mRNA'])[0]
		self.idx_tRNA = np.where(sim_data.process.transcription.rna_data['is_tRNA'])[0]
		self.idx_rprotein = np.where(sim_data.process.transcription.rna_data['includes_ribosomal_protein'])[0]
		self.idx_rnap = np.where(sim_data.process.transcription.rna_data['includes_RNAP'])[0]

		# Synthesis probabilities for different categories of genes
		self.rnaSynthProbFractions = sim_data.process.transcription.rnaSynthProbFraction
		self.rnaSynthProbRProtein = sim_data.process.transcription.rnaSynthProbRProtein
		self.rnaSynthProbRnaPolymerase = sim_data.process.transcription.rnaSynthProbRnaPolymerase

		# Coordinates and transcription directions of transcription units
		self.replication_coordinate = sim_data.process.transcription.rna_data[
			"replication_coordinate"]
		self.transcription_direction = sim_data.process.transcription.rna_data[
			'is_forward']

		# ppGpp control related
		self.n_avogadro = sim_data.constants.n_avogadro
		self.cell_density = sim_data.constants.cell_density
		self.ppgpp = self.bulkMoleculeView(sim_data.molecule_ids.ppGpp)
		self.synth_prob = sim_data.process.transcription.synth_prob_from_ppgpp
		self.copy_number = sim_data.process.replication.get_average_copy_number
		self.get_rnap_active_fraction_from_ppGpp = sim_data.process.transcription.get_rnap_active_fraction_from_ppGpp

	def calculateRequest(self):
		# Get all inactive RNA polymerases
		self.inactive_RNAPs.requestAll()

		# Read current environment
		current_media_id = self._external_states['Environment'].current_media_id

		if self.full_chromosomes.total_count() > 0:
			# Get attributes of promoters
			TU_index, bound_TF = self.promoters.attrs("TU_index", "bound_TF")

			if self.ppgpp_regulation:
				cell_mass = self.readFromListener("Mass", "cellMass") * units.fg
				cell_volume = cell_mass / self.cell_density
				counts_to_molar = 1 / (self.n_avogadro * cell_volume)
				ppgpp_conc = self.ppgpp.total_count() * counts_to_molar
				basal_prob, _ = self.synth_prob(ppgpp_conc, self.copy_number)
				if self.trna_attenuation:
					basal_prob[self.attenuated_rna_indices] += self.attenuation_adjustments
				self.fracActiveRnap = self.get_rnap_active_fraction_from_ppGpp(ppgpp_conc)
				ppgpp_scale = basal_prob[TU_index]
				ppgpp_scale[ppgpp_scale == 0] = 1  # Use original delta prob if no ppGpp basal prob
			else:
				basal_prob = self.basal_prob
				self.fracActiveRnap = self.fracActiveRnapDict[current_media_id]
				ppgpp_scale = 1

			# Calculate probabilities of the RNAP binding to each promoter
			self.promoter_init_probs = (basal_prob[TU_index] + ppgpp_scale *
				np.multiply(self.delta_prob_matrix[TU_index, :], bound_TF).sum(axis=1))

			if len(self.genetic_perturbations) > 0:
				self._rescale_initiation_probs(
					self.genetic_perturbations["fixedRnaIdxs"],
					self.genetic_perturbations["fixedSynthProbs"],
					TU_index)

			# Adjust probabilities to not be negative
			self.promoter_init_probs[self.promoter_init_probs < 0] = 0.0
			self.promoter_init_probs /= self.promoter_init_probs.sum()

			if not self.ppgpp_regulation:
				# Adjust synthesis probabilities depending on environment
				synthProbFractions = self.rnaSynthProbFractions[current_media_id]

				# Create masks for different types of RNAs
				is_mrna = np.isin(TU_index, self.idx_mRNA)
				is_trna = np.isin(TU_index, self.idx_tRNA)
				is_rrna = np.isin(TU_index, self.idx_rRNA)
				is_rprotein = np.isin(TU_index, self.idx_rprotein)
				is_rnap = np.isin(TU_index, self.idx_rnap)
				is_fixed = is_trna | is_rrna | is_rprotein | is_rnap

				# Rescale initiation probabilities based on type of RNA
				self.promoter_init_probs[is_mrna] *= synthProbFractions["mRna"] / self.promoter_init_probs[is_mrna].sum()
				self.promoter_init_probs[is_trna] *= synthProbFractions["tRna"] / self.promoter_init_probs[is_trna].sum()
				self.promoter_init_probs[is_rrna] *= synthProbFractions["rRna"] / self.promoter_init_probs[is_rrna].sum()

				# Set fixed synthesis probabilities for RProteins and RNAPs
				self._rescale_initiation_probs(
					np.concatenate((self.idx_rprotein, self.idx_rnap)),
					np.concatenate((
						self.rnaSynthProbRProtein[current_media_id],
						self.rnaSynthProbRnaPolymerase[current_media_id]
						)),
					TU_index)

				assert self.promoter_init_probs[is_fixed].sum() < 1.0

				# Scale remaining synthesis probabilities accordingly
				scaleTheRestBy = (1. - self.promoter_init_probs[is_fixed].sum()) / self.promoter_init_probs[~is_fixed].sum()
				self.promoter_init_probs[~is_fixed] *= scaleTheRestBy

		# If there are no chromosomes in the cell, set all probs to zero
		else:
			self.promoter_init_probs = np.zeros(self.promoters.total_count())

		self.rnaPolymeraseElongationRate = self.rnaPolymeraseElongationRateDict[current_media_id]
		self.elongation_rates = self.make_elongation_rates(
			self.randomState,
			self.rnaPolymeraseElongationRate.asNumber(units.nt / units.s),
			1,  # want elongation rate, not lengths adjusted for time step
			self.variable_elongation)


	def evolveState(self):
		# no synthesis if no chromosome
		if self.full_chromosomes.total_count() == 0:
			self.writeToListener(
				"RnaSynthProb", "target_rna_synth_prob", np.zeros(self.n_TUs))
			self.writeToListener(
				"RnaSynthProb", "actual_rna_synth_prob", np.zeros(self.n_TUs))
			self.writeToListener(
				"RnaSynthProb", "tu_is_overcrowded", np.zeros(self.n_TUs))
			return

		# Get attributes of promoters
		TU_index, coordinates_promoters, domain_index_promoters, bound_TF = self.promoters.attrs(
			"TU_index", "coordinates", "domain_index", "bound_TF")

		# Construct matrix that maps promoters to transcription units
		n_promoters = self.promoters.total_count()
		TU_to_promoter = scipy.sparse.csr_matrix(
			(np.ones(n_promoters), (TU_index, np.arange(n_promoters))),
			shape = (self.n_TUs, n_promoters))

		# Compute target synthesis probabilities of each transcription unit
		target_TU_synth_probs = TU_to_promoter.dot(self.promoter_init_probs)
		self.writeToListener(
			"RnaSynthProb", "target_rna_synth_prob", target_TU_synth_probs)

		# Calculate RNA polymerases to activate based on probabilities
		# Note: ideally we should be using the actual TU synthesis probabilities
		# here, but the calculation of actual probabilities requires the number
		# of RNAPs to activate. The difference should be very small.
		self.activationProb = self._calculateActivationProb(
			self.fracActiveRnap,
			self.rnaLengths,
			(units.nt / units.s) * self.elongation_rates,
			target_TU_synth_probs)
		n_RNAPs_to_activate = np.int64(
			self.activationProb * self.inactive_RNAPs.count())

		if n_RNAPs_to_activate == 0:
			return

		# Cap the initiation probabilities at the maximum level physically
		# allowed from the known RNAP footprint sizes
		max_p = (self.rnaPolymeraseElongationRate / self.active_rnap_footprint_size
			* (units.s) * self.timeStepSec() / n_RNAPs_to_activate).asNumber()
		self.writeToListener("RnaSynthProb", "max_p", max_p)
		is_overcrowded = (self.promoter_init_probs > max_p)

		while np.any(self.promoter_init_probs > max_p):
			self.promoter_init_probs[is_overcrowded] = max_p
			scale_the_rest_by = (
				(1. - self.promoter_init_probs[is_overcrowded].sum())
				/ self.promoter_init_probs[~is_overcrowded].sum()
				)
			self.promoter_init_probs[~is_overcrowded] *= scale_the_rest_by
			is_overcrowded |= (self.promoter_init_probs > max_p)

		# Compute actual synthesis probabilities of each transcription unit
		actual_TU_synth_probs = TU_to_promoter.dot(self.promoter_init_probs)
		tu_is_overcrowded = TU_to_promoter.dot(is_overcrowded).astype(bool)
		self.writeToListener(
			"RnaSynthProb", "actual_rna_synth_prob", actual_TU_synth_probs)
		self.writeToListener(
			"RnaSynthProb", "tu_is_overcrowded", tu_is_overcrowded)

		# Sample a multinomial distribution of initiation probabilities to
		# determine what promoters are initialized
		n_initiations = self.randomState.multinomial(
			n_RNAPs_to_activate, self.promoter_init_probs)

		# Build array of transcription unit indexes for partially transcribed
		# RNAs and domain indexes for RNAPs
		TU_index_partial_RNAs = np.repeat(TU_index, n_initiations)
		domain_index_rnap = np.repeat(domain_index_promoters, n_initiations)

		# Build arrays of starting coordinates and transcription directions
		coordinates = self.replication_coordinate[TU_index_partial_RNAs]
		is_forward = self.transcription_direction[TU_index_partial_RNAs]

		# Create the active RNA polymerases and get their unique indexes
		RNAP_indexes = self.active_RNAPs.moleculesNew(
			n_RNAPs_to_activate,
			domain_index = domain_index_rnap,
			coordinates = coordinates,
			is_forward = is_forward)

		# Decrement counts of inactive RNAPs
		self.inactive_RNAPs.countDec(n_initiations.sum())

		# Add partially transcribed RNAs
		is_mRNA = np.isin(TU_index_partial_RNAs, self.idx_mRNA)
		self.RNAs.moleculesNew(
			n_RNAPs_to_activate,
			TU_index=TU_index_partial_RNAs,
			transcript_length=np.zeros(cast(int, n_RNAPs_to_activate)),
			is_mRNA=is_mRNA,
			is_full_transcript=np.zeros(cast(int, n_RNAPs_to_activate), dtype=bool),
			can_translate=is_mRNA,
			RNAP_index=RNAP_indexes)

		# Create masks for ribosomal RNAs
		rna_init_event = TU_to_promoter.dot(n_initiations)
		rRNA_initiations = rna_init_event[self.idx_rRNA]

		# Write outputs to listeners
		self.writeToListener(
			"RibosomeData", "rRNA_initiated_TU", rRNA_initiations)
		self.writeToListener(
			"RibosomeData", "rRNA_init_prob_TU", rRNA_initiations / float(n_RNAPs_to_activate))
		self.writeToListener("RibosomeData", "total_rna_init", n_RNAPs_to_activate)
		self.writeToListener("RnapData", "didInitialize", n_RNAPs_to_activate)
		self.writeToListener("RnaSynthProb", "total_rna_init", n_RNAPs_to_activate)
		self.writeToListener("RnapData", "rnaInitEvent", rna_init_event)


	def _calculateActivationProb(self, fracActiveRnap, rnaLengths, rnaPolymeraseElongationRates, synthProb):
		"""
		Calculate expected RNAP termination rate based on RNAP elongation rate
		- allTranscriptionTimes: Vector of times required to transcribe each
		transcript
		- allTranscriptionTimestepCounts: Vector of numbers of timesteps
		required to transcribe each transcript
		- averageTranscriptionTimeStepCounts: Average number of timesteps
		required to transcribe a transcript, weighted by synthesis
		probabilities of each transcript
		- expectedTerminationRate: Average number of terminations in one
		timestep for one transcript
		"""
		allTranscriptionTimes = 1. / rnaPolymeraseElongationRates * rnaLengths
		timesteps = (1. / (self.timeStepSec() * units.s) * allTranscriptionTimes).asNumber()
		allTranscriptionTimestepCounts = np.ceil(timesteps)
		averageTranscriptionTimestepCounts = np.dot(
			synthProb, allTranscriptionTimestepCounts)
		expectedTerminationRate = 1. / averageTranscriptionTimestepCounts

		"""
		Modify given fraction of active RNAPs to take into account early
		terminations in between timesteps
		- allFractionTimeInactive: Vector of probabilities an "active" RNAP
		will in effect be "inactive" because it has terminated during a
		timestep
		- averageFractionTimeInactive: Average probability of an "active" RNAP
		being in effect "inactive", weighted by synthesis probabilities
		- effectiveFracActiveRnap: New higher "goal" for fraction of active
		RNAP, considering that the "effective" fraction is lower than what the
		listener sees
		"""
		allFractionTimeInactive = 1 - (
			1. / (self.timeStepSec() * units.s) * allTranscriptionTimes).asNumber() / allTranscriptionTimestepCounts
		averageFractionTimeInactive = np.dot(allFractionTimeInactive, synthProb)
		effectiveFracActiveRnap = fracActiveRnap / (1 - averageFractionTimeInactive)

		# Return activation probability that will balance out the expected termination rate
		activation_prob = effectiveFracActiveRnap * expectedTerminationRate / (1 - effectiveFracActiveRnap)

		if activation_prob > 1:
			activation_prob = 1.

		return activation_prob


	def _rescale_initiation_probs(
			self, fixed_indexes, fixed_synth_probs, TU_index):
		"""
		Rescales the initiation probabilities of each promoter such that the
		total synthesis probabilities of certain types of RNAs are fixed to
		a predetermined value. For instance, if there are two copies of
		promoters for RNA A, whose synthesis probability should be fixed to
		0.1, each promoter is given an initiation probability of 0.05.
		"""
		for idx, synth_prob in zip(fixed_indexes, fixed_synth_probs):
			fixed_mask = (TU_index == idx)
			self.promoter_init_probs[fixed_mask] = synth_prob / fixed_mask.sum()
