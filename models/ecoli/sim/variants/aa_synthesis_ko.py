"""
Knockout expression of mechanistic amino acid synthesis genes in rich conditions.

Modifies:
	attributes from condition variant
	attributes from gene_knockout variant

Expected variant indices (depends on the enzymes in sim_data.process.metabolism.synthesis_enzymes):
	0: control
	1-58: enzyme to knockout
"""

from .condition import condition
from .gene_knockout import gene_knockout


def aa_synthesis_ko(sim_data, index):
	_, sim_data = condition(sim_data, 2)

	metabolism = sim_data.process.metabolism
	complexation = sim_data.process.complexation
	translation = sim_data.process.translation
	transcription = sim_data.process.transcription

	# Enzymes involved in mechanistic amino acid synthesis
	synthesis_enzymes = metabolism.aa_enzymes
	synthesis_monomers = [
		subunit
		for enzyme in synthesis_enzymes
		for subunit in complexation.get_monomers(enzyme)['subunitIds']
		]

	# Map monomers to RNA for a knockout
	monomer_to_cistron = {monomer['id']: monomer['cistron_id'] for monomer in translation.monomer_data}
	cistron_to_index = {cistron['id']: i for i, cistron in enumerate(transcription.cistron_data)}
	cistrons = [
		monomer_to_cistron[monomer]
		for monomer in synthesis_monomers
		if monomer in monomer_to_cistron
		]

	n_variants = len(cistrons) + 1
	if index > n_variants:
		raise ValueError(f'Variant index {index} is not supported. Choose between 0 and {n_variants}')

	rna_index = cistron_to_index[cistrons[index - 1]] + 1 if index > 0 else 0

	return gene_knockout(sim_data, rna_index)
