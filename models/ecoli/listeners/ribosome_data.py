"""
RibosomeData
"""

from __future__ import absolute_import, division, print_function

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
		self.rrn16S_produced = 0
		self.rrn23S_produced = 0
		self.rrn5S_produced = 0
		self.rrn16S_init_prob = 0.
		self.rrn23S_init_prob = 0.
		self.rrn5S_init_prob = 0.
		self.total_rna_init = 0
		self.processElongationRate = 0.
		self.translationSupply = np.zeros(21, np.float64)
		self.numTrpATerminated = 0.
		self.probTranslationPerTranscript = np.zeros(self.nMonomers, np.float64)

		# Attributes computed by the listener
		self.n_ribosomes_per_transcript = np.zeros(self.nMonomers, np.int64)
		self.n_ribosomes_on_partial_mRNA_per_transcript = np.zeros(self.nMonomers, np.int64)
		self.mRNA_unique_index = np.zeros(30, np.int64)
		self.num_ribosome_on_mRNA = np.zeros(30, np.int64)

	def update(self):
		# Get attributes of RNAs and ribosomes
		RNAs = self.uniqueMolecules.container.objectsInCollection('RNA')
		ribosomes = self.uniqueMolecules.container.objectsInCollection(
			'active_ribosome')
		is_full_transcript_RNA, unique_index_RNA, is_mRNA = RNAs.attrs(
			'is_full_transcript', 'unique_index', 'is_mRNA')
		protein_index_ribosomes, mRNA_index_ribosomes = ribosomes.attrs(
			'protein_index', 'mRNA_index')

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


		# GET MRNA UNIQUE INDEX
		self.mRNA_unique_index = unique_index_RNA[is_mRNA]

		# print(self.mRNA_unique_index[:8])
		# is_ribosome_attached = np.isin(mRNA_unique_index,mRNA_index_ribosomes) # TASK: FIGURE OUT ISIN AND USE THIS TO ADJUST BINCOUNT LENGTH AND EXTRACT INFO

		# Get counts of ribosomes attached to the same mRNA
		bincount_minlength = max(self.mRNA_unique_index)+1
		bincount_ribosome_on_mRNA = np.bincount(mRNA_index_ribosomes, minlength=bincount_minlength)
		self.num_ribosome_on_mRNA = bincount_ribosome_on_mRNA[self.mRNA_unique_index]

		# import ipdb;
		# ipdb.set_trace()
	def tableCreate(self, tableWriter):
		subcolumns = {
			'probTranslationPerTranscript': 'monomerIds',
			'n_ribosomes_per_transcript': 'monomerIds',
			'n_ribosomes_on_partial_mRNA_per_transcript': 'monomerIds',
			}

		tableWriter.writeAttributes(
			monomerIds = self.monomerIds,
			subcolumns = subcolumns)

		tableWriter.set_variable_length_columns(
			'num_ribosome_on_mRNA',
			'mRNA_unique_index'
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
			rrn16S_produced = self.rrn16S_produced,
			rrn23S_produced = self.rrn23S_produced,
			rrn5S_produced = self.rrn5S_produced,
			rrn16S_init_prob = self.rrn16S_init_prob,
			rrn23S_init_prob = self.rrn23S_init_prob,
			rrn5S_init_prob = self.rrn5S_init_prob,
			total_rna_init = self.total_rna_init,
			processElongationRate = self.processElongationRate,
			translationSupply = self.translationSupply,
			numTrpATerminated = self.numTrpATerminated,
			probTranslationPerTranscript = self.probTranslationPerTranscript,
			n_ribosomes_per_transcript = self.n_ribosomes_per_transcript,
			n_ribosomes_on_partial_mRNA_per_transcript = self.n_ribosomes_on_partial_mRNA_per_transcript,
			num_ribosome_on_mRNA = self.num_ribosome_on_mRNA,
			mRNA_unique_index = self.mRNA_unique_index
			)
