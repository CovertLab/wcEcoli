"""
Run an explicit list of new gene designs at many seeds.

`new_gene_internal_shift` lays its designs out as a strict factorial of
expression factors by translation efficiencies, which is what a parameter scan
wants. The paper's single-design panels want the opposite shape: a handful of
chosen designs, each with enough seeds to say something about lineage-to-
lineage variation. Expressed as a factorial that would mean running a whole
grid to use four of its cells.

This variant takes the designs as a list instead, so variant index i is simply
NEW_GENE_DESIGNS[i]. Everything else -- the induction schedule, how expression
and translation efficiency are applied, the media condition -- is reused from
new_gene_internal_shift so the two cannot drift apart.

Index 0 must be the no-GFP control, and it must be run with the same seed
numbers as the other designs: the paper's growth comparisons are all made
against a seed-matched control, and the Figure 1B divergence test is paired
seed by seed.

Modifies:
	sim_data.condition
	sim_data.external_state.current_timeline_id

Modifies (after shift):
	as new_gene_internal_shift.induce_new_genes

Expected variant indices (int, non-negative):
	0: control, new genes knocked out for the whole simulation
	i > 0: NEW_GENE_DESIGNS[i], applied from NEW_GENE_INDUCTION_GEN onward

	A media condition may be selected as in new_gene_internal_shift, by
	adding 1000 * condition_index to the variant index.
"""

from models.ecoli.sim.variants.new_gene_internal_shift import (
	CONTROL_OUTPUT,
	NEW_GENE_INDUCTION_GEN,
	condition,
	determine_new_gene_ids_and_indices,
)

# Each entry is (expression factor, translation efficiency), using the same
# meanings as new_gene_internal_shift: expression is multiplied by
# 10 ** (factor - 1), and translation efficiency is set to the given value.
# Index 0 is the control and must stay (0, 0).
#
# The three induced designs are placeholders until the grid scan has been
# analysed (decision D3). Candidates differ in whether the gfp transcription
# unit is flagged overcrowded, which is what produces the attenuation
# artifacts in the proteome scatter.
# A ladder in RBS strength at a fixed promoter, plus one stronger promoter.
# RBS is the axis that sets ribosomes per mRNA, which is the only quantity
# directional mRNA decay can change for a monocistronic gene like gfp (its
# single cistron starts at position 0, so the degradation front clears it in
# the first timestep and only ribosome protection remains). Varying RBS
# therefore varies the size of the effect under test, while index 1 holds
# transcription at the same level with no translation at all and so should
# show no effect: a built-in negative control.
#
# Completion, growth and GFP share below are from the 20240502 grid and are
# indicative only; they move with the current code.
NEW_GENE_DESIGNS = [
	(0, 0),         # 0: no GFP control
	(8, 0),         # 1: transcribed, not translated   done 1.00  GFP  0%
	(8, 1),         # 2: low burden                    done 0.88  GFP  6%
	(8, 5),         # 3: moderate, fully stable        done 1.00  GFP 25%
	(8, 10),        # 4: high, some lineages lost      done 0.88  GFP 35%
	(9, 5),         # 5: highest, most lineages lost   done 0.75  GFP 44%
]

assert NEW_GENE_DESIGNS[0] == (0, 0), (
	'variant index 0 must be the no-GFP control')
assert all(expression > 0 for expression, _ in NEW_GENE_DESIGNS[1:]), (
	'only index 0 may have an expression factor of 0; a second zero-expression'
	' entry would silently duplicate the control')


def get_new_gene_design(index):
	"""
	Map a variant index to its (expression factor, translation efficiency).

	Returns:
		(design_index, expression_factor, trl_eff_value), where
		expression_factor is already raised to 10 ** (factor - 1) for induced
		designs and 0 for the control.
	"""
	if index < 0 or index >= len(NEW_GENE_DESIGNS):
		raise ValueError(
			'variant index %s is outside NEW_GENE_DESIGNS, which has %d'
			' entries' % (index, len(NEW_GENE_DESIGNS)))
	expression, trl_eff = NEW_GENE_DESIGNS[index]
	expression_factor = 0 if index == 0 else 10 ** (expression - 1)
	return index, expression_factor, trl_eff


def induce_new_gene_design(sim_data, index):
	"""
	Apply one design's expression factor and translation efficiency.

	Mirrors new_gene_internal_shift.induce_new_genes, differing only in where
	the two numbers come from.
	"""
	_, expression_factor, trl_eff_value = get_new_gene_design(index)

	(new_gene_mRNA_ids, new_gene_indices, new_gene_monomer_ids,
		new_gene_monomer_indices) = determine_new_gene_ids_and_indices(sim_data)

	for i in range(len(new_gene_indices)):
		gene_index = new_gene_indices[i]
		monomer_index = new_gene_monomer_indices[i]
		sim_data.adjust_new_gene_final_expression(
			[gene_index], [expression_factor])
		sim_data.process.translation.translation_efficiencies_by_monomer[
			monomer_index] = trl_eff_value


def get_new_gene_design_list_total_indices(sim_data):
	"""Number of variant indices this variant defines."""
	return len(NEW_GENE_DESIGNS)


def new_gene_design_list_print_variant(index):
	"""One-line description of a variant index, for logs and metadata."""
	if index == 0:
		return 'control (no GFP)'
	expression, trl_eff = NEW_GENE_DESIGNS[index]
	return 'expression factor %s, translation efficiency %s' % (
		expression, trl_eff)


def new_gene_design_list(sim_data, index):
	"""
	Apply the variant: hold the design knocked out through the wildtype
	generations, then induce it at NEW_GENE_INDUCTION_GEN and leave it on.

	Unlike new_gene_internal_shift there is no knockout shift; these runs
	exist to characterise a sustained steady state.
	"""
	condition_index = index // 1000
	condition(sim_data, condition_index)
	design_index = index - condition_index * 1000

	# Validate before scheduling, so a bad index fails at variant creation
	# rather than part-way through a generation.
	_, expression_factor, trl_eff_value = get_new_gene_design(design_index)

	setattr(sim_data, 'internal_shift_dict', {})
	if design_index != 0 and NEW_GENE_INDUCTION_GEN != -1:
		# One entry per post-induction generation, because sim_data is
		# reloaded from file between generations.
		sim_data.internal_shift_dict[NEW_GENE_INDUCTION_GEN] = [
			(induce_new_gene_design, design_index)]

	if design_index == 0:
		return CONTROL_OUTPUT, sim_data

	expression, trl_eff = NEW_GENE_DESIGNS[design_index]
	return dict(
		shortName='ng_design_%d' % design_index,
		desc='New gene design %d: expression factor %s, translation'
			' efficiency %s, induced at generation %s'
			% (design_index, expression, trl_eff, NEW_GENE_INDUCTION_GEN),
		), sim_data
