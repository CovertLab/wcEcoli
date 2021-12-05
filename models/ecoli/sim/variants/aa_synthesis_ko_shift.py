"""
Knockout expression of mechanistic amino acid synthesis genes and then shift from rich to minimal.

Modifies:
	attributes from condition variant
	attributes from gene_knockout variant

Expected variant indices (depends on SYNTHESIS_GENES):
	0: control
	1-7: enzyme to knockout
"""

from .gene_knockout import gene_knockout


SYNTHESIS_GENES = [
	'alaC', 'alaA',  # Ala synthesis
	'asnA', 'asnB',  # Asn synthesis
	'gltB', 'gltD', 'gdhA',  # Glt synthesis
	]
SHIFT_TIME = 2 * 3600  # 2 hrs


def aa_synthesis_ko_shift(sim_data, index):
	timeline_id = 'remove_aa'
	sim_data.condition = 'with_aa'
	sim_data.external_state.current_timeline_id = timeline_id
	sim_data.external_state.saved_timelines[timeline_id] = [
		(0, 'minimal_plus_amino_acids'), (SHIFT_TIME, 'minimal')
	]

	# Map genes to RNA for a knockout
	replication = sim_data.process.replication
	transcription = sim_data.process.transcription
	symbol_to_cistron = {gene['symbol']: gene['cistron_id'] for gene in replication.gene_data}
	cistron_to_index = {cistron['id']: i for i, cistron in enumerate(transcription.cistron_data)}

	n_variants = len(SYNTHESIS_GENES) + 1
	if index > n_variants:
		raise ValueError(f'Variant index {index} is not supported. Choose between 0 and {n_variants}')

	rna_index = cistron_to_index[symbol_to_cistron[SYNTHESIS_GENES[index]]] + 1 if index > 0 else 0

	return gene_knockout(sim_data, rna_index)
