"""
relA-null: remove RelA's ppGpp synthesis, leaving SpoT as the only source.

RelA's only role in the model is ppGpp synthesis on ribosomes stalled by
uncharged tRNA, so setting its rate constant to zero is equivalent to a relA
deletion for ppGpp, without the polar effects of removing relA's
transcription units (which also carry mazE/mazF/mazG).

The condition variant is applied first, so the medium can be selected as well.

Modifies:
	attributes modified by the condition variant
	sim_data.constants.k_RelA_ppGpp_synthesis

Expected variant indices:
	0: basal control
	1: basal relA-null (k_RelA = 0)
	2: rich (with_aa) control
	3: rich (with_aa) relA-null (k_RelA = 0)
"""

from .condition import condition


# Condition index (row of condition_defs.tsv) for each medium
CONDITIONS = [0, 1]  # basal, with_aa
CONDITION_NAMES = ['basal', 'with_aa']


def rela_null(sim_data, index):
	condition_index = CONDITIONS[index // 2]
	is_null = index % 2 == 1

	_, sim_data = condition(sim_data, condition_index)
	medium = CONDITION_NAMES[index // 2]

	if not is_null:
		if index == 0:
			return dict(
				shortName="control",
				desc="Control simulation"
				), sim_data
		return dict(
			shortName=f"control_{medium}",
			desc=f"Control simulation in {medium}."
			), sim_data

	sim_data.constants.k_RelA_ppGpp_synthesis *= 0

	if index == 1:
		return dict(
			shortName="relA_null",
			desc="RelA ppGpp synthesis rate set to 0 (relA-null)."
			), sim_data
	return dict(
		shortName=f"relA_null_{medium}",
		desc=f"RelA ppGpp synthesis rate set to 0 (relA-null) in {medium}."
		), sim_data
