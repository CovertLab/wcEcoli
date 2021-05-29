"""
Add one amino acid to the minimal media condition.

Modifies:
	sim_data.external_state.saved_media

Expected variant indices (dependent on order of sim_data.moleculeGroups.aaIDs):
	0-20: adding one amino acid to media
	19: control (adding L-selenocysteine which is already required in media)

TODO:
	Create new media ID for new mixtures?
"""

import numpy as np


def aa_uptake(sim_data, index):
	aa = sim_data.molecule_groups.amino_acids[index][:-3]
	sim_data.external_state.saved_media['minimal'][aa] = np.inf

	return dict(
		shortName = "{}_added".format(aa),
		desc = "Add {} to minimal media.".format(aa)
		), sim_data
