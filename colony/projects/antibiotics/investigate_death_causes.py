'''Script for investigating what characteristics predict death.
'''


import argparse
import copy
import os

from cell.analysis.analyze import Analyzer
from cell.plots.expression_survival_dotplot import (
	plot_expression_survival)
from vivarium.core.composition import plot_agents_multigen
from vivarium.core.experiment import get_in, assoc_path

from colony.constants import OUT_DIR
from colony.compartments.antibiotics import (
	BETA_LACTAMASE_KEY,
	PUMP_KEY,
)


PUMP_PATH = (
	'boundary', 'bulk_molecule_concentrations', PUMP_KEY)
BETA_LACTAMASE_PATH = (
	'boundary', 'bulk_molecule_concentrations', BETA_LACTAMASE_KEY)
ANTIBIOTIC_TIME_RANGE = (0.5, 1)
PATH_TO_AGENTS = ('agents',)
PATH_TO_DEAD = ('boundary', 'dead')


def filter_raw_data_by_time(raw_data, time_range):
	floor, ceil = time_range
	end = max(raw_data.keys())
	filtered = {
		time: time_data
		for time, time_data in raw_data.items()
		if floor * end <= time <= ceil * end
	}
	return filtered


def split_raw_data_by_survival(raw_data):
	# Establish which agents die
	agents_die = set()
	for time_data in raw_data.values():
		agents_data = get_in(time_data, PATH_TO_AGENTS)
		for agent, agent_data in agents_data.items():
			dead = get_in(agent_data, PATH_TO_DEAD, False)
			if dead:
				agents_die.add(agent)

	# Split the data
	survive_data = dict()
	for time in raw_data:
		path = (time,) + PATH_TO_AGENTS
		assoc_path(survive_data, path, dict())
	die_data = copy.deepcopy(survive_data)

	for time, time_data in raw_data.items():
		agents_data = get_in(time_data, PATH_TO_AGENTS)
		for agent, agent_data in agents_data.items():
			dest = die_data if agent in agents_die else survive_data
			path = (time,) + PATH_TO_AGENTS + (agent,)
			assoc_path(dest, path, agent_data)

	return survive_data, die_data


def main():
	if not os.path.exists(OUT_DIR):
		os.makedirs(OUT_DIR)

	parser = argparse.ArgumentParser()
	Analyzer.add_connection_args(parser)
	parser.add_argument(
		'experiment_id',
		type=str,
		help='ID of experiment to analyze.'
	)
	args = parser.parse_args()

	data, _ = Analyzer.get_data(
		args, args.experiment_id)

	fig_pump = plot_expression_survival(
		data, PUMP_PATH,
		'Average AcrAB-TolC Concentration (mmol/L) over Cell Lifetime',
		ANTIBIOTIC_TIME_RANGE,
	)
	fig_pump.savefig(os.path.join(OUT_DIR, 'pump'))
	fig_beta_lactamase = plot_expression_survival(
		data, PUMP_PATH,
		'Average Beta-Lactamase Concentration (mmol/L) over Cell Lifetime',
		ANTIBIOTIC_TIME_RANGE,
	)
	fig_beta_lactamase.savefig(os.path.join(OUT_DIR, 'beta_lactamase'))

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
		survive_data, multigen_settings, OUT_DIR, 'survive')
	plot_agents_multigen(
		die_data, multigen_settings, OUT_DIR, 'die')


if __name__ == '__main__':
	main()
