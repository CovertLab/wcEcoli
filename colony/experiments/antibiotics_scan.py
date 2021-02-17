'''Simulate a colony of E. coli

Each E. coli cell is modeled with wcEcoli. They are placed into a shared
environment using Vivarium.
'''

import argparse
import datetime
import json
import sys

import numpy as np
from vivarium.core.composition import (
	agent_environment_experiment,
	simulate_experiment,
)
from vivarium.core.experiment import timestamp, get_in
from vivarium.library.units import units
from vivarium_cell.composites.lattice import Lattice
from vivarium_cell.processes.diffusion_field import make_gradient

from colony.compartments.antibiotics_simple import (
	ANTIBIOTIC_KEY,
	INITIAL_EXTERNAL_ANTIBIOTIC,
	SimpleAntibioticsCell,
)
from colony.experiments.antibiotics import (
	get_antibiotics_timeline,
	BOUNDS,
	N_BINS,
	NUM_EMISSIONS,
)
from wholecell.utils import filepath as fp


AGENT_NAME = 'simple_antibiotic_cell'
PUMP_CONCENTRATIONS = np.linspace(0, 0.000175, 10)
BETA_LACTAMASE_CONCENTRATIONS = np.linspace(0, 0.0003, 10)
ANTIBIOTIC_THRESHOLD = 0.02
SIMULATION_TIME = 48 * 60
PULSE_CONCENTRATION = 0.1239


def simulate_one(
	simulation_time, pulse_concentration, antibiotic_threshold,
	pump_concentration, beta_lactamase_concentration,
):
	recipe = {ANTIBIOTIC_KEY: INITIAL_EXTERNAL_ANTIBIOTIC}

	antibiotic_pulse = (
		simulation_time * 0.5,
		simulation_time * 0.5,
		pulse_concentration,
	)
	timeline_config = {
		'timeline': get_antibiotics_timeline(
			N_BINS, BOUNDS, [antibiotic_pulse], simulation_time),
	}
	process_config = {
		'_parallel': False,
		'death': {
			'detectors': {
				'antibiotic': {
					'antibiotic_threshold': antibiotic_threshold,
				},
			},
		},
		'antibiotic_transport': {
			'initial_pump': pump_concentration,
		},
		'hydrolysis': {
			'initial_catalyst': beta_lactamase_concentration,
		},
	}
	agent_ids = [AGENT_NAME]

	# These come from the start of a wcEcoli simulation with seed 1.
	# Specifically, they come from experiment 20210125.180927 at time 2.
	initial_agent_state = {
		'boundary': {
			'dry_mass': 351.4542634162847 * units.fg,
			'mass': 1171.8377856780676 * units.fg,
			'volume': 1.0653070778891522 * units.fL,
			'width': 1,
		},
	}
	settings = {
		'emitter': {'type': 'timeseries'},
		'emit_step': max(simulation_time // NUM_EMISSIONS, 1),
		'timeline': timeline_config,
		'experiment_id': timestamp(datetime.datetime.utcnow()),
		'display_info': False,
		'progress_bar': False,
	}
	agents_config = {
		'type': SimpleAntibioticsCell,
		'ids': agent_ids,
		'config': process_config,
	}
	environment_config = {
		'type': Lattice,
		'config': {
			'multibody': {
				'bounds': BOUNDS,
				'jitter_force': 1e-6,
			},
			'diffusion': {
				'bounds': BOUNDS,
				'n_bins': N_BINS,
				'molecules': list(recipe.keys()),
				'depth': 1000,  # Deep to avoid depleting local molecules
				'diffusion': 5e-1,
				'gradient': {
					'type': 'linear',
					'molecules': {
						molecule: {
							'center': (0, 0),
							'base': concentration,
							'slope': 0,
						}
						for molecule, concentration in recipe.items()
					},
				},
			},
		},
	}
	experiment = agent_environment_experiment(
		agents_config=agents_config,
		environment_config=environment_config,
		initial_state={},
		initial_agent_state=initial_agent_state,
		settings=settings,
	)
	settings = {
		'timestep': 1.0,
		'return_raw_data': True,
		'timeline': timeline_config,
	}
	raw_data = simulate_experiment(experiment, settings)
	experiment.end()
	return raw_data


def scan(pump_concentrations, beta_lactamase_concentrations):
	raw_datasets = []  # List of 3-tuples: (pump, beta lactamase, data)
	for pump in pump_concentrations:
		for beta_lactamase in beta_lactamase_concentrations:
			raw_data = simulate_one(
				simulation_time=SIMULATION_TIME,
				pulse_concentration=PULSE_CONCENTRATION,
				antibiotic_threshold=ANTIBIOTIC_THRESHOLD,
				pump_concentration=pump,
				beta_lactamase_concentration=beta_lactamase,
			)
			raw_datasets.append((pump, beta_lactamase, raw_data))
	return raw_datasets


def main():
	parser = argparse.ArgumentParser()
	parser.add_argument(
		'-o', '--out',
		type=str,
		default='antibiotics_scan_data.json',
		help='Output file for raw scan data.',
	)
	args = parser.parse_args()
	results = scan(PUMP_CONCENTRATIONS, BETA_LACTAMASE_CONCENTRATIONS)
	output = {
		'metadata': {
			'python': sys.version.splitlines()[0],
			'git_hash': fp.run_cmdline('git rev-parse HEAD'),
			'git_branch': fp.run_cmdline('git symbolic-ref --short HEAD'),
			'git_status': fp.run_cmdline(
				'git status --porcelain').split('\n'),
			'time': datetime.datetime.utcnow().isoformat() + '+00:00',
		},
		'data': results,
		'parameters': {
			'pump_concentrations': PUMP_CONCENTRATIONS.tolist(),
			'beta_lactamase_concentrations': (
				BETA_LACTAMASE_CONCENTRATIONS.tolist()),
			'agent_name': AGENT_NAME,
		},
	}
	with open(args.out, 'w') as f:
		json.dump(output, f, indent=4)


if __name__ == '__main__':
	main()
