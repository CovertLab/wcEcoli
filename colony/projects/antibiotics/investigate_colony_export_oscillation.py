'''Script for investigating why colony imports, then exports antibiotic

This import and export oscillation may not even occur anymore, which
this investigation will also reveal.
'''

import argparse
import os

from cell.analysis.analyze import Analyzer
from cell.plots.multibody_physics import plot_snapshots

from colony.constants import OUT_DIR
from colony.projects.antibiotics.investigate_utils import (
	filter_raw_data_by_time,
)


TIME_RANGE = (0.51, 1)
FIELDS = ['nitrocefin']


def main():
	parser = argparse.ArgumentParser()
	Analyzer.add_connection_args(parser)
	parser.add_argument(
		'experiment_id',
		type=str,
		help='ID of experiment to analyze.'
	)
	args = parser.parse_args()

	out_dir = os.path.join(OUT_DIR, args.experiment_id)
	if not os.path.exists(out_dir):
		os.makedirs(out_dir)

	data, environment_config = Analyzer.get_data(
		args, args.experiment_id)
	data = filter_raw_data_by_time(data, TIME_RANGE)
	snapshots_data = Analyzer.format_data_for_snapshots(
		data, environment_config)
	plot_config = {
		'out_dir': out_dir,
		'filename': 'snapshots_during_pulse',
		'include_fields': FIELDS,
	}
	plot_snapshots(snapshots_data, plot_config)


if __name__ == '__main__':
	main()
