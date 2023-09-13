"""
RibosomeData
"""

import numpy as np

import wholecell.listeners.listener


VERBOSE = False

class RibosomeData(wholecell.listeners.listener.Listener):
	""" RibosomeData """

	_name = 'RibosomeData'

	# Constructor
	def __init__(self, *args, **kwargs):
		super(RibosomeData, self).__init__(*args, **kwargs)


	# Construct object graph
	def initialize(self, sim, sim_data):
		super(RibosomeData, self).initialize(sim, sim_data)

		self.uniqueMolecules = sim.internal_states['UniqueMolecules']

		self.monomerIds = sim_data.process.translation.monomer_data['id'].tolist()
		self.nMonomers = len(self.monomerIds)
		self.rRNA_cistron_tu_mapping_matrix = sim_data.process.transcription.rRNA_cistron_tu_mapping_matrix
		self.cistron_tu_mapping_matrix = sim_data.process.transcription.cistron_tu_mapping_matrix
		self.rRNA_is_5S = sim_data.process.transcription.cistron_data['is_5S_rRNA'][sim_data.process.transcription.cistron_data['is_rRNA']]
		self.rRNA_is_16S = sim_data.process.transcription.cistron_data['is_16S_rRNA'][sim_data.process.transcription.cistron_data['is_rRNA']]
		self.rRNA_is_23S = sim_data.process.transcription.cistron_data['is_23S_rRNA'][sim_data.process.transcription.cistron_data['is_rRNA']]
		self.n_rRNAs = sim_data.process.transcription.rna_data['is_rRNA'].sum()
		self.n_RNAs = len(sim_data.process.transcription.rna_data)

	# Allocate memory
	def allocate(self):
		super(RibosomeData, self).allocate()

		# Attributes broadcast by the PolypeptideElongation process
		self.aaCountInSequence = np.zeros(21, np.int64)
		self.aaCounts = np.zeros(21, np.int64)
		self.actualElongations = 0
		self.actualElongationHist = np.zeros(22, np.int64)
		self.elongationsNonTerminatingHist = np.zeros(22, np.int64)
		self.didTerminate = 0
		self.didInitialize = 0
		self.terminationLoss = 0
		self.effectiveElongationRate = 0.
		self.rRNA_initiated_TU = np.zeros(self.n_rRNAs, np.int64)
		self.rRNA_init_prob_TU = np.zeros(self.n_rRNAs, np.float64)
		self.total_rRNA_initiated = 0
		self.total_rRNA_init_prob = 0.
		self.rRNA16S_initiated = 0
		self.rRNA23S_initiated = 0
		self.rRNA5S_initiated = 0
		self.rRNA16S_init_prob = 0.
		self.rRNA23S_init_prob = 0.
		self.rRNA5S_init_prob = 0.
		self.total_rna_init = 0
		self.processElongationRate = 0.
		self.translationSupply = np.zeros(21, np.float64)
		self.numTrpATerminated = 0.
		self.target_prob_translation_per_transcript = np.zeros(self.nMonomers,
													 np.float64)
		self.actual_prob_translation_per_transcript = np.zeros(self.nMonomers,
														np.float64)
		self.mRNA_is_overcrowded = np.zeros(self.nMonomers, np.float64)
		self.ribosome_init_event_per_monomer = np.zeros(self.nMonomers,
													   np.int64)

		# Attributes computed by the listener
		self.n_ribosomes_per_transcript = np.zeros(self.nMonomers, np.int64)
		self.n_ribosomes_on_partial_mRNA_per_transcript = np.zeros(self.nMonomers, np.int64)
		self.n_ribosomes_on_each_mRNA = np.zeros([], np.int64)

	def update(self):
		# Get attributes of RNAs and ribosomes
		RNAs = self.uniqueMolecules.container.objectsInCollection('RNA')
		ribosomes = self.uniqueMolecules.container.objectsInCollection(
			'active_ribosome')
		is_full_transcript_RNA, unique_index_RNA, can_translate, TU_index = RNAs.attrs(
			'is_full_transcript', 'unique_index', 'can_translate', 'TU_index')
		protein_index_ribosomes, mRNA_index_ribosomes, massDiff_protein_ribosomes = ribosomes.attrs(
			'protein_index', 'mRNA_index', 'massDiff_protein')

		# Get mask for ribosomes that are translating proteins on partially
		# transcribed mRNAs
		ribosomes_on_nascent_mRNA_mask = np.isin(
			mRNA_index_ribosomes,
			unique_index_RNA[np.logical_not(is_full_transcript_RNA)])

		# Get counts of ribosomes for each type
		self.n_ribosomes_per_transcript = np.bincount(
			protein_index_ribosomes, minlength=self.nMonomers)
		self.n_ribosomes_on_partial_mRNA_per_transcript = np.bincount(
			protein_index_ribosomes[ribosomes_on_nascent_mRNA_mask],
			minlength=self.nMonomers)

		rRNA_cistrons_produced = self.rRNA_cistron_tu_mapping_matrix.dot(self.rRNA_initiated_TU)
		rRNA_cistrons_init_prob = self.rRNA_cistron_tu_mapping_matrix.dot(self.rRNA_init_prob_TU)
		self.total_rRNA_initiated = np.sum(self.rRNA_initiated_TU)
		self.total_rRNA_init_prob = np.sum(self.rRNA_init_prob_TU)
		self.rRNA5S_initiated = np.sum(rRNA_cistrons_produced[self.rRNA_is_5S])
		self.rRNA16S_initiated = np.sum(rRNA_cistrons_produced[self.rRNA_is_16S])
		self.rRNA23S_initiated = np.sum(rRNA_cistrons_produced[self.rRNA_is_23S])
		self.rRNA5S_init_prob = np.sum(rRNA_cistrons_init_prob[self.rRNA_is_5S])
		self.rRNA16S_init_prob = np.sum(rRNA_cistrons_init_prob[self.rRNA_is_16S])
		self.rRNA23S_init_prob = np.sum(rRNA_cistrons_init_prob[self.rRNA_is_23S])

		# Get fully transcribed translatable mRNA index
		is_full_mRNA = can_translate & is_full_transcript_RNA
		self.mRNA_unique_index = unique_index_RNA[is_full_mRNA]
		self.mRNA_TU_index = TU_index[is_full_mRNA]


		# Get counts of ribosomes attached to the same mRNA
		bincount_minlength = max(self.mRNA_unique_index) + 1
		bincount_ribosome_on_mRNA = np.bincount(mRNA_index_ribosomes, minlength=bincount_minlength)
		self.n_ribosomes_on_each_mRNA = bincount_ribosome_on_mRNA[self.mRNA_unique_index]

		# Initialize array representing protein mass for each mRNA
		self.protein_mass_on_polysomes = np.zeros(
			len(self.n_ribosomes_on_each_mRNA))

		# Get protein mass on each polysome
		self.protein_mass_on_polysomes = np.bincount(
			mRNA_index_ribosomes, weights=massDiff_protein_ribosomes,
			minlength=bincount_minlength)[self.mRNA_unique_index]

	def tableCreate(self, tableWriter):
		subcolumns = {
			'target_prob_translation_per_transcript': 'monomerIds',
			'actual_prob_translation_per_transcript': 'monomerIds',
			'mRNA_is_overcrowded': 'monomerIds',
			'n_ribosomes_per_transcript': 'monomerIds',
			'n_ribosomes_on_partial_mRNA_per_transcript': 'monomerIds',
			}

		tableWriter.writeAttributes(
			monomerIds = self.monomerIds,
			subcolumns = subcolumns)

		tableWriter.set_variable_length_columns(
			'n_ribosomes_on_each_mRNA',
			'mRNA_TU_index',
			'protein_mass_on_polysomes',
			)

	def tableAppend(self, tableWriter):
		tableWriter.append(
			time = self.time(),
			simulationStep = self.simulationStep(),
			aaCountInSequence = self.aaCountInSequence,
			aaCounts = self.aaCounts,
			actualElongations = self.actualElongations,
			actualElongationHist = self.actualElongationHist,
			elongationsNonTerminatingHist = self.elongationsNonTerminatingHist,
			didTerminate = self.didTerminate,
			didInitialize = self.didInitialize,
			terminationLoss = self.terminationLoss,
			effectiveElongationRate = self.effectiveElongationRate,
			total_rRNA_initiated = self.total_rRNA_initiated,
			total_rRNA_init_prob = self.total_rRNA_init_prob,
			rRNA16S_initiated = self.rRNA16S_initiated,
			rRNA23S_initiated = self.rRNA23S_initiated,
			rRNA5S_initiated = self.rRNA5S_initiated,
			rRNA16S_init_prob = self.rRNA16S_init_prob,
			rRNA23S_init_prob = self.rRNA23S_init_prob,
			rRNA5S_init_prob = self.rRNA5S_init_prob,
			total_rna_init = self.total_rna_init,
			processElongationRate = self.processElongationRate,
			translationSupply = self.translationSupply,
			numTrpATerminated = self.numTrpATerminated,
			target_prob_translation_per_transcript=self
			.target_prob_translation_per_transcript,
			actual_prob_translation_per_transcript=self
			.actual_prob_translation_per_transcript,
			mRNA_is_overcrowded = self.mRNA_is_overcrowded,
			ribosome_init_event_per_monomer = self.ribosome_init_event_per_monomer,
			n_ribosomes_per_transcript = self.n_ribosomes_per_transcript,
			n_ribosomes_on_partial_mRNA_per_transcript = self.n_ribosomes_on_partial_mRNA_per_transcript,
			n_ribosomes_on_each_mRNA = self.n_ribosomes_on_each_mRNA,
			mRNA_TU_index = self.mRNA_TU_index,
			protein_mass_on_polysomes = self.protein_mass_on_polysomes,
			)
