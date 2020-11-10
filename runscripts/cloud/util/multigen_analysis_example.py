"""A small example for CustomWorkflow to run within a Docker Firetask.
Given a seed dir, run example multigen analysis on its sim generations.
"""

import os
import sys

from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from wholecell.io.tablereader import TableReader


def print_time_range(sim_out_dir):
	# type: (str) -> None
	"""Print the time range from the given cell simOut."""
	times = TableReader(os.path.join(sim_out_dir, "Main")).readColumn("time")
	print(f'    Time range: [{times[0]} .. {times[-1]}]')


def print_bulk_count_range(sim_out_dir):
	# type: (str) -> None
	"""Print the start and end BulkMolecules counts from the given cell simOut."""
	counts = TableReader(os.path.join(sim_out_dir, "BulkMolecules")).readColumn("counts")
	print(f'    BulkMolecules counts start / end: {counts[0]} / {counts[-1]}')


def analyze(seed_dir):
	# type: (str) -> None
	"""Examine each simOut dir within the seed_dir."""
	paths = AnalysisPaths(seed_dir, multi_gen_plot=True)
	for cell_dir in paths.get_cells():
		print(cell_dir)
		sim_out_dir = os.path.join(cell_dir, 'simOut', '')
		print_time_range(sim_out_dir)
		print_bulk_count_range(sim_out_dir)

	# TODO(jerry): Write an output file into seed_dir/../count_out/<seed_key>/complex/


if __name__ == '__main__':
	seed_dir_arg = sys.argv[1]  # e.g. /wcEcoli/out/counts/wildtype_000000/000000
	analyze(seed_dir_arg)
