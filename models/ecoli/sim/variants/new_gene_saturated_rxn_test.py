"""
Test whether the kcat-saturated reactions identified in
`saturated_rxns_var22_vs_var23.csv` are causally responsible for the
doubling-time slowdown observed in the 1.0x max kcat regime of
`new_gene_trl_eff_sweep`.

Constrains kcat for **only** a hand-picked subset of reactions (4 or 16) and
leaves the rest unconstrained, alongside controls that constrain everything
("full") or nothing. By comparing matched gfp expression at four constraint
scopes (none / full / 16-rxn / 4-rxn) we can isolate how much of the dt
phenotype the saturated subset alone reproduces.

Layout (17 indices total):

  0:        Control: no gfp, no kcat constraint
  1-3:      gfp expression, no kcat constraint, trl_eff = 4.0/4.5/5.0
  4-6:      gfp expression, full kcat constraint @ 1.0x, trl_eff = 4.0/4.5/5.0
  7-9:      gfp expression, 4-rxn kcat constraint @ 1.0x, trl_eff = 4.0/4.5/5.0
  10-12:    gfp expression, 16-rxn kcat constraint @ 1.0x, trl_eff = 4.0/4.5/5.0
  13:       no gfp, 4-rxn kcat constraint @ 1.0x
  14:       no gfp, 16-rxn kcat constraint @ 1.0x
  15:       no gfp, 4-rxn kcat constraint @ 0.4x
  16:       no gfp, 16-rxn kcat constraint @ 0.4x

Modifies:
	sim_data.condition
	sim_data.external_state.current_timeline_id

Modifies (after shift, gfp variants only):
	sim_data.process.transcription.rna_synth_prob
	sim_data.process.transcription.rna_expression
	sim_data.process.transcription.exp_free
	sim_data.process.transcription.exp_ppgpp
	sim_data.process.transcription.attenuation_basal_prob_adjustments
	sim_data.process.transcription_regulation.basal_prob
	sim_data.process.transcription_regulation.delta_prob
	sim_data.process.translation.translation_efficiencies_by_monomer

Modifies (constrained variants only):
	sim_data.process.metabolism.kcat_estimate_quantile
	sim_data.process.metabolism.kcat_estimate_multiplier
	sim_data.process.metabolism.selected_kcat_estimates
"""

from models.ecoli.sim.variants.new_gene_internal_shift import (
	condition,
	determine_new_gene_ids_and_indices,
	NEW_GENE_INDUCTION_GEN,
	NEW_GENE_KNOCKOUT_GEN,
	CONTROL_OUTPUT,
)

TRL_EFF_LEVELS = [4.0, 4.5, 5.0]
EXPRESSION_FACTOR = 8.5

# 4-reaction subset: high-saturation gainers from var22 -> var23 in the
# `saturated_rxns_var22_vs_var23.csv` analysis.
RXNS_4 = [
	('HOLO-ACP-SYNTH-RXN__HOLO-ACP-SYNTH-CPLX', 'HOLO-ACP-SYNTH-CPLX[c]'),
	('L-GLN-FRUCT-6-P-AMINOTRANS-RXN',          'L-GLN-FRUCT-6-P-AMINOTRANS-CPLX[c]'),
	('GLUTAMIN-RXN__CARBPSYN-CPLX',             'CARBPSYN-CPLX[c]'),
	('DXPREDISOM-RXN',                          'DXPREDISOM-CPLX[i]'),
]

# 16-reaction full set: every reaction flagged by the saturation analysis.
RXNS_16 = [
	('KETOGLUTREDUCT-RXN (reverse)',                       'PGLYCDEHYDROG-CPLX[c]'),
	('3-DEHYDROQUINATE-SYNTHASE-RXN',                      'AROB-MONOMER[c]'),
	('KETOGLUTREDUCT-RXN',                                 'PGLYCDEHYDROG-CPLX[c]'),
	('METBALT-RXN',                                        'O-SUCCHOMOSERLYASE-CPLX[c]'),
	('GLUTAMIN-RXN__CTPSYN-CPLX',                          'CTPSYN-CPLX[c]'),
	('HOLO-ACP-SYNTH-RXN__HOLO-ACP-SYNTH-CPLX',            'HOLO-ACP-SYNTH-CPLX[c]'),
	('L-GLN-FRUCT-6-P-AMINOTRANS-RXN',                     'L-GLN-FRUCT-6-P-AMINOTRANS-CPLX[c]'),
	('GLUTAMIN-RXN__CARBPSYN-CPLX',                        'CARBPSYN-CPLX[c]'),
	('RXN-9535__FABB-CPLX',                                'FABB-CPLX[c]'),
	('GLUTRACE-RXN (reverse)',                             'GLUTRACE-MONOMER[c]'),
	('DXPREDISOM-RXN',                                     'DXPREDISOM-CPLX[i]'),
	('NACGLCTRANS-RXN (reverse)',                          'NACGLCTRANS-MONOMER[m]'),
	('CHORISMATEMUT-RXN__CHORISMUTPREPHENDEHYDROG-CPLX',   'CHORISMUTPREPHENDEHYDROG-CPLX[c]'),
	('RXN-16701 (reverse)',                                'PGLYCDEHYDROG-CPLX[c]'),
	('RXN-16701',                                          'PGLYCDEHYDROG-CPLX[c]'),
	('CTPSYN-RXN',                                         'CTPSYN-CPLX[c]'),
]

assert set(RXNS_4).issubset(set(RXNS_16)), \
	"RXNS_4 must be a strict subset of RXNS_16"

# (gfp, rxn_set, multiplier, trl_eff)
# rxn_set is one of: None, 'full', 'rxns4', 'rxns16'
VARIANTS = [
	(False, None,     None, None),                            # 0
	*[(True,  None,     None, t) for t in TRL_EFF_LEVELS],    # 1-3
	*[(True,  'full',   1.0,  t) for t in TRL_EFF_LEVELS],    # 4-6
	*[(True,  'rxns4',  1.0,  t) for t in TRL_EFF_LEVELS],    # 7-9
	*[(True,  'rxns16', 1.0,  t) for t in TRL_EFF_LEVELS],    # 10-12
	(False, 'rxns4',  1.0, None),                             # 13
	(False, 'rxns16', 1.0, None),                             # 14
	(False, 'rxns4',  0.4, None),                             # 15
	(False, 'rxns16', 0.4, None),                             # 16
]
N_VARIANTS = len(VARIANTS)  # 17


def _rxn_set_keys(rxn_set):
	"""Return the list of (rxn_id, catalyst_id) tuples for a given rxn_set tag."""
	if rxn_set == 'rxns4':
		return RXNS_4
	if rxn_set == 'rxns16':
		return RXNS_16
	raise ValueError(f"Unknown rxn_set tag: {rxn_set!r}")


def _condition_label(index):
	"""Return a short, descriptive label for a variant index."""
	gfp, rxn_set, multiplier, trl_eff = VARIANTS[index]
	if rxn_set is None:
		constraint = 'unconstrained'
	else:
		mult_str = f'{multiplier:.1f}x'
		if rxn_set == 'full':
			constraint = f'fullkcat_{mult_str}'
		else:
			constraint = f'{rxn_set}_{mult_str}'
	if gfp:
		gfp_str = f'gfp_trl{trl_eff}'
	else:
		gfp_str = 'nogfp'
	return f'{gfp_str}_{constraint}'


def _induce(sim_data, variant_index):
	"""
	Apply the new-gene induction at NEW_GENE_INDUCTION_GEN for this variant.

	For non-gfp variants, set expression to 0 (knockout). For gfp variants, set
	expression to 10^(EXPRESSION_FACTOR-1) and translation efficiency to the
	variant's trl_eff value.
	"""
	gfp, _rxn_set, _multiplier, trl_eff = VARIANTS[variant_index]
	new_gene_mRNA_ids, new_gene_indices, new_gene_monomer_ids, \
		new_gene_monomer_indices = determine_new_gene_ids_and_indices(sim_data)

	if not gfp:
		for i in range(len(new_gene_indices)):
			gene_index = new_gene_indices[i]
			sim_data.adjust_new_gene_final_expression([gene_index], [0])
		return

	expression_factor = 10 ** (EXPRESSION_FACTOR - 1)
	for i in range(len(new_gene_indices)):
		gene_index = new_gene_indices[i]
		monomer_index = new_gene_monomer_indices[i]
		sim_data.adjust_new_gene_final_expression(
			[gene_index], [expression_factor])
		sim_data.process.translation.translation_efficiencies_by_monomer[
			monomer_index] = trl_eff


def new_gene_saturated_rxn_test(sim_data, index):
	"""
	Apply variant. Looks up the (gfp, rxn_set, multiplier, trl_eff) row from
	VARIANTS, applies the kcat constraint (if any), and schedules the new gene
	induction at NEW_GENE_INDUCTION_GEN.
	"""
	assert 0 <= index < N_VARIANTS, \
		f"variant index {index} out of range [0, {N_VARIANTS})"
	gfp, rxn_set, multiplier, _trl_eff = VARIANTS[index]

	# Apply kcat constraint (if any) before media setup so the saved sim_data
	# pickle reflects the constrained metabolism.
	if rxn_set is not None:
		assert 'max' in sim_data.process.metabolism.kcat_estimates, \
			"sim_data.process.metabolism.kcat_estimates['max'] missing" \
			" -- Parca must populate this"
		metabolism = sim_data.process.metabolism
		max_kcats = metabolism.kcat_estimates['max']
		if rxn_set == 'full':
			selected = {
				key: value * multiplier for key, value in max_kcats.items()
			}
		else:
			keys = _rxn_set_keys(rxn_set)
			missing = [k for k in keys if k not in max_kcats]
			assert not missing, (
				f"keys missing from kcat_estimates['max']: {missing}")
			selected = {key: max_kcats[key] * multiplier for key in keys}
		metabolism.kcat_estimate_quantile = 'max'
		metabolism.kcat_estimate_multiplier = multiplier
		metabolism.selected_kcat_estimates = selected

	# Always minimal media (condition index 0).
	condition(sim_data, 0)

	# Schedule the new-gene induction shift.
	setattr(sim_data, 'internal_shift_dict', {})
	if NEW_GENE_INDUCTION_GEN != -1:
		sim_data.internal_shift_dict[NEW_GENE_INDUCTION_GEN] = [
			(_induce, index)]
	if NEW_GENE_KNOCKOUT_GEN != -1:
		# _induce(sim_data, index) for a non-gfp variant performs a knockout,
		# so we reuse it for the knockout shift if one is scheduled.
		sim_data.internal_shift_dict[NEW_GENE_KNOCKOUT_GEN] = [
			(_induce, 0)]  # variant 0 is non-gfp -> sets expression to 0

	# Build human-readable metadata.
	if index == 0:
		# Wildtype-equivalent control.
		return CONTROL_OUTPUT, sim_data

	label = _condition_label(index)
	if gfp:
		desc = (
			f"Expression factor 10^{EXPRESSION_FACTOR - 1},"
			f" translation efficiency {_trl_eff}, {label}.")
	else:
		desc = f"New gene knockout, {label}."
	return dict(shortName=label, desc=desc), sim_data
