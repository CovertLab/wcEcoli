"""
RnaSynthProb

Records RNA synthesis probabilities
"""

import numpy as np

import wholecell.listeners.listener


class RnaSynthProb(wholecell.listeners.listener.Listener):
	""" RnaSynthProb """

	_name = "RnaSynthProb"

	# Constructor
	def __init__(self, *args, **kwargs):
		super(RnaSynthProb, self).__init__(*args, **kwargs)


	# Construct object graph
	def initialize(self, sim, sim_data):
		super(RnaSynthProb, self).initialize(sim, sim_data)

		self.uniqueMolecules = sim.internal_states['UniqueMolecules']

		self.transcriptInitiation = sim.processes["TranscriptInitiation"]
		self.rna_ids = sim_data.process.transcription.rna_data["id"]
		self.n_TU = len(self.rna_ids)
		self.cistron_ids = sim_data.process.transcription.cistron_data["id"]
		self.gene_ids = sim_data.process.transcription.cistron_data['gene_id']
		self.n_cistron = len(self.cistron_ids)
		self.cistron_tu_mapping_matrix = sim_data.process.transcription.cistron_tu_mapping_matrix

		self.tf_ids = sim_data.process.transcription_regulation.tf_ids
		self.n_TF = len(self.tf_ids)


	# Allocate memory
	def allocate(self):
		super(RnaSynthProb, self).allocate()

		self.target_rna_synth_prob = np.zeros(self.n_TU, np.float64)
		self.actual_rna_synth_prob = np.zeros(self.n_TU, np.float64)
		self.tu_is_overcrowded = np.zeros(self.n_TU, np.float64)
		self.promoter_copy_number = np.zeros(self.n_TU, np.int16)
		self.rna_synth_prob_per_cistron = np.zeros(self.n_cistron, np.float64)
		self.total_rna_init = 0
		self.expected_rna_init_per_cistron = np.zeros(self.n_cistron, np.float64)

		self.pPromoterBound = np.zeros(self.n_TF, np.float64)
		self.nPromoterBound = np.zeros(self.n_TF, np.float64)
		self.nActualBound = np.zeros(self.n_TF, np.float64)
		self.n_available_promoters = np.zeros(self.n_TF, np.float64)

		# These arrays gets flattened at tableAppend(). Resulting array should
		# be reshaped before use.
		self.n_bound_TF_per_TU = np.zeros((self.n_TU, self.n_TF), np.int16)
		self.n_bound_TF_per_cistron = np.zeros((self.n_TU, self.n_cistron), np.int16)

		# Properties of bound TFs
		self.bound_TF_indexes = np.array([], np.int64)
		self.bound_TF_coordinates = np.array([], np.int64)
		self.bound_TF_domains = np.array([], np.int64)


	def update(self):
		promoters = self.uniqueMolecules.container.objectsInCollection('promoter')
		TU_indexes, all_coordinates, all_domains, bound_TFs = promoters.attrs(
			"TU_index", "coordinates", "domain_index", "bound_TF"
			)
		genes = self.uniqueMolecules.container.objectsInCollection('gene')
		cistron_indexes = genes.attr("cistron_index")

		self.promoter_copy_number = np.bincount(TU_indexes, minlength=self.n_TU)
		self.gene_copy_number = np.bincount(cistron_indexes, minlength=self.n_cistron)

		bound_promoter_indexes, TF_indexes = np.where(bound_TFs)

		self.bound_TF_indexes = TF_indexes
		self.bound_TF_coordinates = all_coordinates[bound_promoter_indexes]
		self.bound_TF_domains = all_domains[bound_promoter_indexes]

		actual_rna_synth_prob_per_cistron = self.cistron_tu_mapping_matrix.dot(
			self.actual_rna_synth_prob)
		# The expected value of rna initiations per cistron. Realized values
		# during simulation will be different, because they will be integers
		# drawn from a multinomial distribution
		self.expected_rna_init_per_cistron = actual_rna_synth_prob_per_cistron * self.total_rna_init

		if actual_rna_synth_prob_per_cistron.sum() != 0:
			self.actual_rna_synth_prob_per_cistron = actual_rna_synth_prob_per_cistron / actual_rna_synth_prob_per_cistron.sum()
		else:
			self.actual_rna_synth_prob_per_cistron = actual_rna_synth_prob_per_cistron
		target_rna_synth_prob_per_cistron = self.cistron_tu_mapping_matrix.dot(
			self.target_rna_synth_prob)
		if target_rna_synth_prob_per_cistron.sum() != 0:
			self.target_rna_synth_prob_per_cistron = target_rna_synth_prob_per_cistron / target_rna_synth_prob_per_cistron.sum()
		else:
			self.target_rna_synth_prob_per_cistron = target_rna_synth_prob_per_cistron

		self.n_bound_TF_per_cistron = self.cistron_tu_mapping_matrix.dot(
			self.n_bound_TF_per_TU).astype(np.int16).T


	def tableCreate(self, tableWriter):
		subcolumns = {
			'promoter_copy_number': 'rnaIds',
			'gene_copy_number': 'gene_ids',
			'target_rna_synth_prob': 'rnaIds',
			'actual_rna_synth_prob': 'rnaIds',
			'tu_is_overcrowded': 'rnaIds',
			'actual_rna_synth_prob_per_cistron': 'cistron_ids',
			'target_rna_synth_prob_per_cistron': 'cistron_ids',
			'expected_rna_init_per_cistron': 'cistron_ids',
			'pPromoterBound': 'tf_ids',
			'nPromoterBound': 'tf_ids',
			'nActualBound': 'tf_ids',
			'n_available_promoters': 'tf_ids',
			'n_bound_TF_per_TU': 'rnaIds',
			'n_bound_TF_per_cistron': 'cistron_ids'
			}

		tableWriter.writeAttributes(
			rnaIds = list(self.rna_ids),
			cistron_ids = list(self.cistron_ids),
			tf_ids = list(self.tf_ids),
			gene_ids = list(self.gene_ids),
			subcolumns = subcolumns)

		tableWriter.set_variable_length_columns(
			'bound_TF_indexes',
			'bound_TF_coordinates',
			'bound_TF_domains',
			)


	def tableAppend(self, tableWriter):
		tableWriter.append(
			time = self.time(),
			simulationStep = self.simulationStep(),
			target_rna_synth_prob = self.target_rna_synth_prob,
			actual_rna_synth_prob =	self.actual_rna_synth_prob,
			tu_is_overcrowded = self.tu_is_overcrowded,
			actual_rna_synth_prob_per_cistron = self.actual_rna_synth_prob_per_cistron,
			target_rna_synth_prob_per_cistron = self.target_rna_synth_prob_per_cistron,
			expected_rna_init_per_cistron = self.expected_rna_init_per_cistron,
			promoter_copy_number = self.promoter_copy_number,
			gene_copy_number = self.gene_copy_number,
			pPromoterBound = self.pPromoterBound,
			nPromoterBound = self.nPromoterBound,
			nActualBound = self.nActualBound,
			n_available_promoters = self.n_available_promoters,
			n_bound_TF_per_TU = self.n_bound_TF_per_TU,
			n_bound_TF_per_cistron = self.n_bound_TF_per_cistron,
			bound_TF_indexes = self.bound_TF_indexes,
			bound_TF_coordinates = self.bound_TF_coordinates,
			bound_TF_domains = self.bound_TF_domains,
			)
