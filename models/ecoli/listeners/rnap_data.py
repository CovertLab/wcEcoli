"""
RnapData
"""

import numpy as np

import wholecell.listeners.listener
from models.ecoli.processes.transcript_elongation import get_mapping_arrays

VERBOSE = False

class RnapData(wholecell.listeners.listener.Listener):
	""" RnapData """

	_name = 'RnapData'

	# Constructor
	def __init__(self, *args, **kwargs):
		super(RnapData, self).__init__(*args, **kwargs)


	# Construct object graph
	def initialize(self, sim, sim_data):
		super(RnapData, self).initialize(sim, sim_data)

		self.rnaIds = sim_data.process.transcription.rna_data['id']
		self.stable_RNA_indexes = np.where(
			np.logical_or(sim_data.process.transcription.rna_data['is_rRNA'],
			sim_data.process.transcription.rna_data['is_tRNA']))[0]
		self.nRnaSpecies = self.rnaIds.size
		self.cistron_ids = sim_data.process.transcription.cistron_data['id']
		self.n_cistrons = self.cistron_ids.size
		self.cistron_tu_mapping_matrix = sim_data.process.transcription.cistron_tu_mapping_matrix
		self.uniqueMolecules = sim.internal_states['UniqueMolecules']


	# Allocate memory
	def allocate(self):
		super(RnapData, self).allocate()

		# Attributes broadcast by the PolypeptideElongation process
		self.actualElongations = 0
		self.didTerminate = 0
		self.didInitialize = 0
		self.terminationLoss = 0
		self.rnaInitEvent = np.zeros(self.nRnaSpecies, np.int64)
		self.incomplete_transcription_event = np.zeros(self.nRnaSpecies, np.int64)
		self.rna_init_event_per_cistron = np.zeros(self.n_cistrons, np.int64)
		self.didStall = 0

		# Collisions with replisomes
		self.n_total_collisions = 0
		self.n_headon_collisions = 0
		self.n_codirectional_collisions = 0
		self.n_removed_ribosomes = 0

		# Entries with variable lengths
		self.active_rnap_coordinates = np.array([], np.int64)
		self.active_rnap_domain_indexes = np.array([], np.int32)
		self.active_rnap_unique_indexes = np.array([], np.int64)
		self.active_rnap_n_bound_ribosomes = np.array([], np.int64)
		self.headon_collision_coordinates = np.array([], np.int64)
		self.codirectional_collision_coordinates = np.array([], np.int64)


	def update(self):
		active_rnaps = self.uniqueMolecules.container.objectsInCollection(
			'active_RNAP')
		RNAs = self.uniqueMolecules.container.objectsInCollection('RNA')
		active_ribosomes = self.uniqueMolecules.container.objectsInCollection('active_ribosome')

		# Read coordinates of all active RNAPs
		coordinates, domain_indexes, RNAP_unique_indexes = active_rnaps.attrs(
			"coordinates", "domain_index", "unique_index")
		self.active_rnap_coordinates = coordinates
		self.active_rnap_domain_indexes = domain_indexes
		self.active_rnap_unique_indexes = RNAP_unique_indexes

		RNA_RNAP_index, is_full_transcript, RNA_unique_indexes, TU_indexes = RNAs.attrs(
			'RNAP_index', 'is_full_transcript', 'unique_index', 'TU_index')
		is_partial_transcript = np.logical_not(is_full_transcript)
		is_stable_RNA = np.isin(TU_indexes, self.stable_RNA_indexes)
		partial_RNA_RNAP_indexes = RNA_RNAP_index[is_partial_transcript]
		partial_RNA_unique_indexes = RNA_unique_indexes[is_partial_transcript]
		self.partial_stable_RNA_RNAP_indexes = RNA_RNAP_index[np.logical_and(is_stable_RNA, is_partial_transcript)]

		ribosome_RNA_index = active_ribosomes.attr('mRNA_index')

		RNA_index_counts = dict(
			zip(*np.unique(ribosome_RNA_index, return_counts=True)))

		partial_RNA_to_RNAP_mapping, _ = get_mapping_arrays(
			partial_RNA_RNAP_indexes, RNAP_unique_indexes)

		self.active_rnap_n_bound_ribosomes = np.array(
			[RNA_index_counts.get(partial_RNA_unique_indexes[i], 0)
				for i in partial_RNA_to_RNAP_mapping])

		# Calculate hypothetical RNA initiation events per cistron
		self.rna_init_event_per_cistron = self.cistron_tu_mapping_matrix.dot(
			self.rnaInitEvent)


	def tableCreate(self, tableWriter):
		subcolumns = {
			'rnaInitEvent': 'rnaIds',
			'rna_init_event_per_cistron': 'cistron_ids',
			'incomplete_transcription_event': 'rnaIds',
			}

		tableWriter.writeAttributes(
			rnaIds = list(self.rnaIds),
			cistron_ids = list(self.cistron_ids),
			subcolumns = subcolumns)

		tableWriter.set_variable_length_columns(
			'active_rnap_coordinates',
			'active_rnap_domain_indexes',
			'active_rnap_unique_indexes',
			'active_rnap_on_stable_RNA_indexes',
			'active_rnap_n_bound_ribosomes',
			'headon_collision_coordinates',
			'codirectional_collision_coordinates',
			)


	def tableAppend(self, tableWriter):
		tableWriter.append(
			time = self.time(),
			simulationStep = self.simulationStep(),
			active_rnap_coordinates=self.active_rnap_coordinates,
			active_rnap_domain_indexes=self.active_rnap_domain_indexes,
			active_rnap_unique_indexes=self.active_rnap_unique_indexes,
			active_rnap_on_stable_RNA_indexes=self.partial_stable_RNA_RNAP_indexes,
			active_rnap_n_bound_ribosomes=self.active_rnap_n_bound_ribosomes,
			actualElongations = self.actualElongations,
			codirectional_collision_coordinates=self.codirectional_collision_coordinates,
			didTerminate = self.didTerminate,
			didInitialize = self.didInitialize,
			didStall = self.didStall,
			headon_collision_coordinates=self.headon_collision_coordinates,
			incomplete_transcription_event = self.incomplete_transcription_event,
			n_total_collisions=self.n_total_collisions,
			n_headon_collisions=self.n_headon_collisions,
			n_codirectional_collisions=self.n_codirectional_collisions,
			n_removed_ribosomes=self.n_removed_ribosomes,
			terminationLoss = self.terminationLoss,
			rnaInitEvent = self.rnaInitEvent,
			rna_init_event_per_cistron = self.rna_init_event_per_cistron,
			)
