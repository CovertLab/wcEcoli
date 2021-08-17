"""
TODO

Modifies:
	TODO

Expected variant indices (dependent on length of FACTORS and sim_data.molecule_groups.amino_acids):
	0-6: range of values for first parameter, first amino acid, first media condition
	0-34: range of parameters and values for first amino acid, first media condition
	0-734: range of parameters and values over all amino acids for first media condition
	0-1469: all changes
"""

from .add_one_aa import add_one_aa

import numpy as np


FACTORS = [0, 0.1, 0.5, 1.5, 2, 5, 10]  # TODO: run factor of 1 once for each media condition
PARAMETERS = ['aa_kcats', 'aa_kis', 'aa_upstream_kms', 'aa_reverse_kms', 'aa_degradation_kms']
N_PARAM_VALUES = len(FACTORS) * len(PARAMETERS)
MEDIA_IDS = [5, 19]  # Glt and control for now - TODO: run this for all AA additions?


def get_media_index(index, n_aas):
	return MEDIA_IDS[index // (N_PARAM_VALUES * n_aas)]

def get_aa_index(index, n_aas):
	sub_index = index % (N_PARAM_VALUES * n_aas)
	return sub_index // N_PARAM_VALUES

def get_adjustment(index):
	sub_index = index % N_PARAM_VALUES
	param_index = sub_index // len(FACTORS)
	factor_index = sub_index % len(FACTORS)
	return PARAMETERS[param_index], FACTORS[factor_index]

def aa_synthesis_sensitivity(sim_data, index):
	n_aas = len(sim_data.molecule_groups.amino_acids)

	# Use add_one_aa variant to add a specific amino acid to the media
	media_idx = get_media_index(index, n_aas)
	_, sim_data = add_one_aa(sim_data, media_idx)

	# Change the uptake rate for that amino acid to check the sensitivity
	aa_idx = get_aa_index(index, n_aas)
	param, factor = get_adjustment(index)
	values = getattr(sim_data.process.metabolism, param)
	if np.all(values[aa_idx] == values[aa_idx] * factor):
		raise ValueError('No change to params - not running variant sims.')
	values[aa_idx] *= factor

	aa_id = sim_data.molecule_groups.amino_acids[aa_idx]

	return dict(
		shortName=f'TODO',
		desc=f'TODO'
		), sim_data
