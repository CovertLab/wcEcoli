"""
Remove one amino acid from the minimal media plus amino acids condition.

Modifies:
	sim_data.external_state.saved_media

Expected variant indices (dependent on order of sim_data.moleculeGroups.aaIDs):
	0-20: adding one amino acid to media
	19: control (L-selenocysteine must be in media)

TODO:
	Create new media ID for new mixtures?
"""

def remove_one_aa(sim_data, index):
	condition_label = sim_data.ordered_conditions[2]
	sim_data.condition = condition_label
	sim_data.external_state.current_timeline_id = condition_label
	sim_data.external_state.saved_timelines[condition_label] = [
		(0, sim_data.conditions[condition_label]["nutrients"])
	]
	aa = sim_data.molecule_groups.amino_acids[index][:-3]
	if index != 19:
		sim_data.external_state.saved_media['minimal_plus_amino_acids'][aa] = 0

	return dict(
		shortName = "{}_removed".format(aa),
		desc = "Remove {} from rich media.".format(aa)
		), sim_data
