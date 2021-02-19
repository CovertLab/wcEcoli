from argparse import ArgumentParser
from datetime import datetime
import json
import sys

import numpy as np
from vivarium.core.experiment import get_in

from colony.experiments.antibiotics_scan import (
	simulate_one,
	SIMULATION_TIME,
	PULSE_CONCENTRATION,
	ANTIBIOTIC_THRESHOLD,
	AGENT_NAME,
)
import wholecell.utils.filepath as fp


LOGREG_BOUNDARY_M = -0.1824149289775941
LOGREG_BOUNDARY_B = 0.1343552138570287e-3  # mM
LOGREG_Y_ERROR = 0.05e-3  # mM
X_LOCATIONS = np.linspace(0, 0.18e-3, 10)  # mM
DESIRED_PRECISION = 1e-6  # mM
DEATH_PATH = ('agents', AGENT_NAME, 'boundary', 'dead')


def _logreg_estimate(x):
	y = LOGREG_BOUNDARY_M * x + LOGREG_BOUNDARY_B
	return y, LOGREG_Y_ERROR


def get_metadata():
	metadata = {
		'git_hash': fp.run_cmdline('git rev-parse HEAD'),
		'git_branch': fp.run_cmdline(
			'git symbolic-ref --short HEAD'),
		'git_status': fp.run_cmdline(
			'git status --porcelain').split('\n'),
		'time': datetime.utcnow().isoformat() + '+00:00',
		'python': sys.version.splitlines()[0],
	}
	return metadata


def evaluate_point(x, y):
	data = simulate_one(
		SIMULATION_TIME, PULSE_CONCENTRATION,
		ANTIBIOTIC_THRESHOLD, x, y)
	dead = get_in(data[max(data.keys())], DEATH_PATH)
	return dead


def search(estimate_func, x_values, precision):
	y_values = []
	points = []
	for x in x_values:
		y_estimate, y_error = estimate_func(x)
		upper = y_estimate + y_error
		if evaluate_point(x, upper):
			raise ValueError(
				f'Incorrectly assuming cell lives at ({x}, {upper}) mM')
		lower = y_estimate - y_error
		if not evaluate_point(x, lower):
			raise ValueError(
				f'Incorrectly assuming cell dies at ({x}, {lower}) mM')
		while upper - lower > precision:
			middle = (upper + lower) / 2
			dead = evaluate_point(x, middle)
			points.append((x, middle, dead))
			if dead:
				lower = middle
			else:
				upper = middle
		y_values.append((upper + lower) / 2)
	return y_values, points


def main(tokens=None):
	parser = ArgumentParser()
	parser.add_argument(
		'out_path', type=str, help='Path to output JSON file')
	args = parser.parse_args(tokens)
	y_values, points = search(_logreg_estimate, X_LOCATIONS, DESIRED_PRECISION)
	output = {
		'x_values': X_LOCATIONS.tolist(),
		'y_values': y_values,
		'points': points,
		'estimate_m': LOGREG_BOUNDARY_M,
		'estimate_b': LOGREG_BOUNDARY_B,
		'estiamte_error': LOGREG_Y_ERROR,
		'precision': DESIRED_PRECISION,
		'simulation_time': SIMULATION_TIME,
		'pulse_concentration': PULSE_CONCENTRATION,
		'antibiotic_threshold': ANTIBIOTIC_THRESHOLD,
		'metadata': get_metadata(),
	}
	with open(args.out_path, 'w') as f:
		json.dump(output, f, indent=4)


if __name__ == '__main__':
	main()
