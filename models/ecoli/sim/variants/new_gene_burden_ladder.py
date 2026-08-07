"""
Sweep new-gene (GFP) translation efficiency at a single fixed expression
factor to produce one clean burden axis.

This variant is deliberately minimal: it is reused unchanged across every
insertion-position and copy-number batch so that results from different
genomes are directly comparable. The only thing that varies between those
batches is NEW_GENES (which flat-file directory supplies the construct), not
the variant definition.

Layout (8 indices total):

  0:     Knockout control. Expression factor 0 -- no construct transcription
         and no construct protein.
  1:     Transcription-only control. Full expression, translation efficiency 0
         -- the construct is transcribed (~7% of mRNA mass at exp 8) but makes
         no protein, so it imposes transcriptional load without ribosome load.
         This is also the correct baseline for the product decomposition: it
         has non-zero output, so unlike the knockout it does not make the
         dosage share undefined.
  2-7:   exp = EXPRESSION_FACTOR, trl_eff = TRL_EFF_VALUES[index - 1]

Modifies:
	sim_data.condition
	sim_data.external_state.current_timeline_id

Modifies (after shift):
	sim_data.process.transcription.rna_synth_prob
	sim_data.process.transcription.rna_expression
	sim_data.process.transcription.exp_free
	sim_data.process.transcription.exp_ppgpp
	sim_data.process.transcription.attenuation_basal_prob_adjustments
	sim_data.process.transcription_regulation.basal_prob
	sim_data.process.transcription_regulation.delta_prob
	sim_data.process.translation.translation_efficiencies_by_monomer
"""

from models.ecoli.sim.variants.new_gene_internal_shift import (
	condition,
	determine_new_gene_ids_and_indices,
	NEW_GENE_INDUCTION_GEN,
	NEW_GENE_KNOCKOUT_GEN,
	CONTROL_OUTPUT,
)

# Ascending burden. Index 0 of this list is variant index 1.
#
# Extended to 10.0 on 2026-08-07, alongside the drop to expression factor 8.
# Lowering the expression factor alone would have compressed the burden range
# badly -- over the old 0.1-5.0 ladder, exp 8 gives tau 54 -> 73 min against
# 52.5 -> 103 at exp 8.5. Pushing translation efficiency to 10 recovers most of
# it (tau 54 -> 87, RNAP -50%) and costs nothing on the crowding axis, because
# the cap depends on promoter strength rather than on translation. Dashboard
# figures, exp 8: overcrowded fraction 0.0 at every translation efficiency
# including 10, and 88% of sims still reach generation 24.
#
# The leading 0.0 is the transcription-only control described above.
TRL_EFF_VALUES = [0.0, 0.5, 1.0, 2.5, 5.0, 7.5, 10.0]

# Per the convention in new_gene_internal_shift, an expression factor of x
# multiplies new gene expression by 10^(x - 1).
#
# Lowered from 8.5 to 8 on 2026-08-07. At 8.5 the construct's promoter is
# overcrowded in 94.6% of timesteps -- pinned at the RNAP-footprint cap in
# transcript_initiation.py:269 -- which makes its per-copy initiation rate
# insensitive to burden by construction and drives the measured dosage share
# toward 100% artefactually. The existing exp/trl_eff dashboard shows the
# transition is a cliff between 8 and 9: overcrowded fraction 0.0 at 8 and 1.0
# at 9, with initiation rate 0.45-0.73 against 0.98-0.99. 8 is the highest
# value with headroom at every translation efficiency, and it also clears the
# generation-24 completion problem in the exp 9/10 corner.
#
# NOTE: Batches 1 and 2 were run at 8.5. Anything compared across that boundary
# is not comparable.
EXPRESSION_FACTOR = 8

N_TRL_EFF = len(TRL_EFF_VALUES)  # 7
N_VARIANTS = N_TRL_EFF + 1  # 8 (knockout + trl_eff 0 control + 6 burden levels)


def is_control(index):
	"""Return True if this variant is the control (GFP knockout)."""
	return index == 0


def _induce(sim_data, index):
	"""
	Induce new genes at this variant's expression factor and translation
	efficiency.

	For index == 0 (control), expression is set to 0 (knockout). For
	index >= 1, expression is set to 10^(EXPRESSION_FACTOR - 1) and
	translation efficiency is set to TRL_EFF_VALUES[index - 1]. Index 1 is
	therefore full expression at zero translation efficiency.

	Every new gene index returned by determine_new_gene_ids_and_indices is
	given the same settings, so tandem copies at one locus are all induced
	together with identical per-copy expression.
	"""
	new_gene_mRNA_ids, new_gene_indices, new_gene_monomer_ids, \
		new_gene_monomer_indices = determine_new_gene_ids_and_indices(
		sim_data)

	if index == 0:
		# Control: knock out new gene expression
		for i in range(len(new_gene_indices)):
			gene_index = new_gene_indices[i]
			sim_data.adjust_new_gene_final_expression([gene_index], [0])
	else:
		expression_factor = 10 ** (EXPRESSION_FACTOR - 1)
		trl_eff_value = TRL_EFF_VALUES[index - 1]

		for i in range(len(new_gene_indices)):
			gene_index = new_gene_indices[i]
			monomer_index = new_gene_monomer_indices[i]

			sim_data.adjust_new_gene_final_expression(
				[gene_index], [expression_factor])
			sim_data.process.translation.translation_efficiencies_by_monomer[
				monomer_index] = trl_eff_value


def new_gene_burden_ladder(sim_data, index):
	"""
	Apply variant. Sets minimal media, then schedules new gene induction at
	this variant's expression factor and translation efficiency.
	"""
	assert 0 <= index < N_VARIANTS, (
		f"new_gene_burden_ladder index must be in [0, {N_VARIANTS - 1}],"
		f" got {index}")

	# Always minimal media
	condition(sim_data, 0)

	# Initialize internal shift dictionary
	setattr(sim_data, 'internal_shift_dict', {})

	# Add the new gene induction to the internal_shift instructions
	if NEW_GENE_INDUCTION_GEN != -1:
		sim_data.internal_shift_dict[NEW_GENE_INDUCTION_GEN] = [
			(_induce, index)]
	if NEW_GENE_KNOCKOUT_GEN != -1:
		sim_data.internal_shift_dict[NEW_GENE_KNOCKOUT_GEN] = [
			(_induce, 0)]  # knockout = set expression to 0

	# Build description
	if is_control(index):
		return dict(CONTROL_OUTPUT), sim_data

	trl_eff_value = TRL_EFF_VALUES[index - 1]
	shortName = f'trl_eff_{trl_eff_value}_exp_{EXPRESSION_FACTOR}'
	desc = (
		f'Expression factor 10^{EXPRESSION_FACTOR - 1}, '
		f'translation efficiency {trl_eff_value}.')
	return dict(shortName=shortName, desc=desc), sim_data
