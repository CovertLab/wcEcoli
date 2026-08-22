"""
The new_gene_burden_ladder translation-efficiency sweep, run with the expected
gene copy number frozen at its unburdened value.

Everything except that one number is inherited from new_gene_burden_ladder: the
same eight variant indices, the same TRL_EFF_VALUES, the same EXPRESSION_FACTOR,
the same wildtype-coordinate pin, the same minimal medium. Keeping them shared is
the point -- this batch and the P4 minimal batch are meant to differ in exactly
one input, so any difference in the result is attributable to the freeze.

What is being frozen and why
----------------------------
transcription.synth_prob_from_ppgpp turns fitted per-cell expression into a
per-copy initiation probability as

	prob = normalize(expression * loss / n_avg_copy)

with `loss = growth + deg_rate` and `n_avg_copy` the Cooper-Helmstetter
expectation at doubling time `tau`. Because every promoter of a gene receives the
full per-copy probability, the share a gene actually draws goes as

	n_actual * prob  ~  expression * loss * (n_actual / n_avg_copy)

so the availability term a = n_actual / n_avg_copy is what decides whether a loss
of real gene dosage reaches transcription at all. If the expectation tracked the
cell, a would stay at 1 and dosage loss would cancel exactly.

It does not track the cell. `tau` is inferred from ppGpp through
interpolate_linearized_fit against a curve fitted to five wild-type growth
conditions (transcription.py:180), and under GFP burden that sensor sees only
part of the slowdown: about half of it in minimal media, and none of it in rich,
where ppGpp falls 26.6 -> 25.9 uM while the realized doubling time rises
27.7 -> 42.4 min. So a fraction of the dosage loss propagates and the rest is
cancelled, and the fraction is set by a fitting device rather than by biology.

Freezing `n_avg_copy` at a constant removes the cancellation entirely: a then
tracks n_actual, and the full measured dosage loss reaches transcription. The gap
between this ladder and the unfrozen one is the amount of the machinery's
loss that arrives via gene copy number -- the loop gain. That number is not
currently measured, and "there is a copy-number feedback" and "the copy-number
feedback matters" are different claims.

This is a counterfactual, not a bug fix. Real E. coli does not divide promoter
strength by expected gene dosage; the division is a device for making fitted
expression growth-rate-consistent. Freezing it is not more realistic than
leaving it, it just makes the device's contribution measurable by removing it.

Only the per-timestep response is frozen
----------------------------------------
The freeze is read by transcript_initiation and passed as `tau_override`, so
the Parca fit and reconstruction.ecoli.initialization both keep using the
doubling time they infer from the nominal condition ppGpp. Initial conditions are
therefore unchanged and only the response to burden differs.

Modifies:
	sim_data.process.transcription.frozen_expectation_tau
	everything new_gene_burden_ladder modifies
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

# Doubling time in minutes at which the expected copy number is held, for every
# variant including the knockout.
#
# Measured, not nominal: this is the mean ppGpp-inferred tau of the knockout
# control in the P4 minimal batch (new_gene_dosage_compensation_p4_minimal.csv,
# variant 0, tau_inferred = 58.52), which is the value the model itself was
# already using before GFP was induced. Anchoring here rather than at the
# condition's nominal 44 min keeps the low rungs of this ladder close to the
# unfrozen batch, so the two can be read against each other rung by rung.
#
# Only the constancy matters, not the particular value. Once the divisor stops
# moving, a = n_actual / n_avg_copy tracks n_actual alone, so the change in a
# across the ladder -- the quantity this batch measures -- is identical for any
# choice. A different value shifts the level of a and the baseline dosage ratio
# between gene classes, which is why it is worth anchoring on the unburdened
# measurement, but it cannot change the slope.
FROZEN_EXPECTATION_TAU = 58.5


def new_gene_burden_ladder_frozen(sim_data, index):
	"""
	Apply variant. Freezes the expected gene copy number, sets minimal media,
	then schedules new gene induction at this variant's expression factor and
	translation efficiency.

	Identical to new_gene_burden_ladder apart from the freeze, so that the two
	batches differ in exactly one input.
	"""
	assert 0 <= index < N_VARIANTS, (
		f"new_gene_burden_ladder_frozen index must be in"
		f" [0, {N_VARIANTS - 1}], got {index}")

	# Read by transcript_initiation when the process is constructed, and passed
	# to synth_prob_from_ppgpp as tau_override on every timestep. Applied to
	# every index including the knockout, which is what makes the control a
	# control: the freeze is a property of the batch, not of the burden level.
	sim_data.process.transcription.frozen_expectation_tau = (
		FROZEN_EXPECTATION_TAU)

	# Must happen before the simulation starts; the generation-8 internal shift
	# only rewrites expression arrays, so the pin persists through induction.
	if PIN_WT_COORDINATE:
		_pin_wt_coordinate(sim_data)

	# Always minimal media, as in the parent variant
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
		out = dict(CONTROL_OUTPUT)
		out['shortName'] = 'control_frozen'
		out['desc'] = (
			f'Control simulation, expected copy number frozen at'
			f' {FROZEN_EXPECTATION_TAU} min.')
		return out, sim_data

	trl_eff_value = TRL_EFF_VALUES[index - 1]
	shortName = f'frozen_trl_eff_{trl_eff_value}_exp_{EXPRESSION_FACTOR}'
	desc = (
		f'Expected copy number frozen at {FROZEN_EXPECTATION_TAU} min.'
		f' Expression factor 10^{EXPRESSION_FACTOR - 1},'
		f' translation efficiency {trl_eff_value}.')
	return dict(shortName=shortName, desc=desc), sim_data
