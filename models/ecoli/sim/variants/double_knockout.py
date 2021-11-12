"""
Knockout expression of two genes

Modifies:

Expected variant indices:
"""

from .gene_knockout import gene_knockout


KNOCKOUTS = [
	(1, 2),
	]


def double_knockout(sim_data, index):
	for ko in KNOCKOUTS[index]:
		_, sim_data = gene_knockout(sim_data, ko+1)

	return dict(
		shortName = '',
		desc = '',
		), sim_data
