"""
REFERENCE: docs/misc/internal_shift.md

UPDATE: describe what variant should do and any analysis plots to use with it
UPDATE: change Modifies and Expected variant indices sections below (examples shown can be removed)
UPDATE: describe the nature of the shifts, how many there are, etc.

Modifies (after shift):
	sim_data.something

Expected variant indices (dependent on sorted order of sim_data.conditionActiveTfs):
	0: control
	1: fill in
	2: fill in
	etc.
"""

def internal_shift_function(sim_data, index):
	# UPDATE: modify sim_data attributes based on the variant index
	sim_data.something = index


# UPDATE: give a descriptive function name that matches the file name
def template(sim_data, index):
	# Initialize internal shift dictionary
	sim_data.internal_shifts.internal_shift_dict = {}

	# Add desired shifts to the dictionary
	# Multiple shift functions per generation are allowed, just add tuples to
	# the list, they will be applied in list order
	sim_data.internal_shifts.internal_shift_dict[shift_gen_index] = [
		(internal_shift_function, index)]

	# UPDATE: set strings to give a description of the variant/what has changed
	return dict(
		shortName=f'{index}',
		desc=f'something is set at {index}.'
		), sim_data
