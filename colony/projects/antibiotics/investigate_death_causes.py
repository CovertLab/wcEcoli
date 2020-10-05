'''Script for investigating what characteristics predict death.
'''


import argparse
import os

from vivarium_cell.analysis.analyze import Analyzer
from vivarium_cell.plots.expression_survival_dotplot import (
	plot_expression_survival)
from vivarium_cell.plots.multibody_physics import plot_tags
from vivarium.plots.agents_multigen import plot_agents_multigen

from colony.constants import OUT_DIR
from colony.compartments.antibiotics import (
	BETA_LACTAMASE_KEY,
	PUMP_KEY,
)
from colony.projects.antibiotics.investigate_utils import (
	filter_raw_data_by_time,
	split_raw_data_by_survival,
)


PUMP_PATH = (
	'boundary', 'bulk_molecule_concentrations', PUMP_KEY)
BETA_LACTAMASE_PATH = (
	'boundary', 'bulk_molecule_concentrations', BETA_LACTAMASE_KEY)
ANTIBIOTIC_TIME_RANGE = (0.5, 1)


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

	fig_pump = plot_expression_survival(
		data, PUMP_PATH,
		'Average AcrAB-TolC Concentration (mmol/L) over Cell Lifetime',
		ANTIBIOTIC_TIME_RANGE,
	)
	fig_pump.savefig(os.path.join(out_dir, 'pump'))
	fig_beta_lactamase = plot_expression_survival(
		data, BETA_LACTAMASE_PATH,
		'Average Beta-Lactamase Concentration (mmol/L) over Cell Lifetime',
		ANTIBIOTIC_TIME_RANGE,
	)
	fig_beta_lactamase.savefig(os.path.join(out_dir, 'beta_lactamase'))

	multigen_settings = {
		'include_paths': [
			PUMP_PATH,
			BETA_LACTAMASE_PATH,
			('boundary', 'bulk_molecules_report', PUMP_KEY),
			('boundary', 'bulk_molecules_report', BETA_LACTAMASE_KEY),
			('boundary', 'cytoplasm', 'nitrocefin_hydrolyzed'),
			('boundary', 'cytoplasm', 'nitrocefin'),
			('boundary', 'external', 'nitrocefin'),
			('boundary', 'dead'),
		],
	}
	filtered = filter_raw_data_by_time(data, ANTIBIOTIC_TIME_RANGE)
	survive_data, die_data = split_raw_data_by_survival(filtered)
	plot_agents_multigen(
		survive_data, multigen_settings, out_dir, 'survive')
	plot_agents_multigen(
		die_data, multigen_settings, out_dir, 'die')

	tag_path_name_map = {
		PUMP_PATH: 'AcrAB-TOlC',
		BETA_LACTAMASE_PATH: 'Beta-Lactamase',
	}
	tags_config = {
		'out_dir': out_dir,
		'tagged_molecules': tag_path_name_map.keys(),
		'filename': 'expression',
		'tag_path_name_map': tag_path_name_map,
	}
	tags_data = Analyzer.format_data_for_tags(data, environment_config)
	plot_tags(tags_data, tags_config)


if __name__ == '__main__':
	main()
