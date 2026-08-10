"""
The new_gene_burden_ladder translation-efficiency sweep, run in rich media
(minimal + amino acids) instead of glucose minimal.

Everything except the medium is inherited from new_gene_burden_ladder: the same
eight variant indices, the same TRL_EFF_VALUES, the same EXPRESSION_FACTOR, the
same wildtype-coordinate pin. Keeping them shared is the point -- the rich and
minimal batches are meant to differ in exactly one input, so any difference in
the result is attributable to the medium.

Why rich media is a different experiment and not just a faster one
------------------------------------------------------------------
The fitted doubling time is 25.0 min against minimal's 44.0
(condition_defs.tsv), so replication rounds overlap far more heavily. Cooper-
Helmstetter, n = 2^(((1-f)C + D)/tau) with C = 40 and D = 20, then puts the
origin-to-terminus copy-number gradient at about 2.4x in rich media against
1.65x in minimal. Every gene in the transcription and translation machinery is
origin-proximal -- rRNA operons at a median 12% of the replichore, RNAP
subunits between 4% and 31%, and 41 of 54 ribosomal protein genes inside 25% --
so rich media hands that machinery a much larger dosage advantage, and
therefore a much larger absolute amount to lose when burden slows the cell.

Whether it loses proportionally more is what this batch measures.

A warning about the crowding cap
--------------------------------
transcript_initiation.py:269 caps per-promoter initiation at

    max_p = (elongation_rate / footprint) * timestep / n_RNAPs_to_activate

which is inversely proportional to the polymerase pool. Rich media has a
substantially larger pool, so max_p is SMALLER and the construct is MORE likely
to be pinned at the ceiling at the same expression factor -- not less. The
crowding gate that cleared EXPRESSION_FACTOR = 8 in minimal
(n_pinned_tus = 0 at every variant) does NOT transfer, and has to be re-run
before any rich-media result is interpreted. See the gate in
POSITION_SWEEP_LAUNCH.md, applied to this batch.

Modifies:
	sim_data.condition
	sim_data.external_state.current_timeline_id
	sim_data.process.transcription.rna_data['wt_replication_coordinate']

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

from models.ecoli.sim.variants.new_gene_burden_ladder import (
	EXPRESSION_FACTOR,
	is_control,
	N_TRL_EFF,
	N_VARIANTS,
	PIN_WT_COORDINATE,
	REFERENCE_WT_COORDINATE,
	TRL_EFF_VALUES,
	_induce,
	_pin_wt_coordinate,
)
from models.ecoli.sim.variants.new_gene_internal_shift import (
	condition,
	NEW_GENE_INDUCTION_GEN,
	NEW_GENE_KNOCKOUT_GEN,
	CONTROL_OUTPUT,
)

# Index into sim_data.ordered_conditions. 0 is basal (glucose minimal), 1 is
# with_aa (minimal + amino acids), the model's rich medium. Same convention as
# models/ecoli/sim/variants/condition.py.
CONDITION_INDEX = 1


def new_gene_burden_ladder_rich(sim_data, index):
	"""
	Apply variant. Sets rich media, then schedules new gene induction at this
	variant's expression factor and translation efficiency.

	Identical to new_gene_burden_ladder apart from the medium, so that the two
	batches differ in exactly one input.
	"""
	assert 0 <= index < N_VARIANTS, (
		f"new_gene_burden_ladder_rich index must be in [0, {N_VARIANTS - 1}],"
		f" got {index}")

	# Must happen before the simulation starts; the generation-8 internal shift
	# only rewrites expression arrays, so the pin persists through induction.
	if PIN_WT_COORDINATE:
		_pin_wt_coordinate(sim_data)

	# Rich media, rather than the minimal media the parent variant sets
	condition(sim_data, CONDITION_INDEX)

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
		out = dict(CONTROL_OUTPUT)
		out['shortName'] = 'control_rich'
		out['desc'] = 'Control simulation, rich media.'
		return out, sim_data

	trl_eff_value = TRL_EFF_VALUES[index - 1]
	shortName = f'rich_trl_eff_{trl_eff_value}_exp_{EXPRESSION_FACTOR}'
	desc = (
		f'Rich media. Expression factor 10^{EXPRESSION_FACTOR - 1}, '
		f'translation efficiency {trl_eff_value}.')
	return dict(shortName=shortName, desc=desc), sim_data
