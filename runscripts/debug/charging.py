#! /usr/bin/env python
"""
Tools and analysis to debug charging problems.
"""

import argparse
import os
import pickle
from typing import Tuple

import numpy as np

from models.ecoli.processes.polypeptide_elongation import (calculate_trna_charging,
	CONC_UNITS, get_charging_params)
from wholecell.io.tablereader import TableReader
from wholecell.utils import constants, scriptBase


class ChargingDebug(scriptBase.ScriptBase):
	def define_parameters(self, parser):
		super().define_parameters(parser)

		# Path args
		self.define_parameter_sim_dir(parser)
		self.define_parameter_variant_index(parser)
		parser.add_argument('-s', '--seed', type=int, default=0,
			help='The initial simulation number (int). The value will get'
				 ' formatted as a subdirectory name like "000000". Default = 0.')
		parser.add_argument('-g', '--generation', type=int, default=0,
			help='The generation number (int). The value will get formatted'
				 ' as a subdirectory name like "generation_000000". Default = 0.')
		parser.add_argument('-d', '--daughter', type=int, default=0,
			help='The daughter number (int). The value will get formatted as'
				 ' a subdirectory name like "000000". Default = 0.')

		# Sim options
		self.define_parameter_bool(parser, 'variable_elongation_translation',
			default_key='variable_elongation_translation',
			help='set if sims were run with variable_elongation_translation')

		# Debug options
		parser.add_argument('--validation', type=int, default=1,
			help='Number of time steps to run for validation. If < 0, will run all.')

	def update_args(self, args):
		super().update_args(args)

		# Extract data from args
		variant_dir_name, _, _ = args.variant_dir
		seed_str = '%06d' % (args.seed,)
		gen_str = 'generation_%06d' % (args.generation,)
		daughter_str = '%06d' % (args.daughter,)
		dirs = os.path.join(seed_str, gen_str, daughter_str)
		input_variant_directory = os.path.join(args.sim_path, variant_dir_name)

		# Set paths from args
		args.sim_data_file = os.path.join(input_variant_directory, 'kb',
			constants.SERIALIZED_SIM_DATA_MODIFIED)
		args.sim_out_dir = os.path.join(input_variant_directory, dirs, 'simOut')

	def load_data(self, sim_data_file: str, sim_out_dir: str) -> None:
		"""
		Loads sim_data and simulation output data from files and saves it as
		instance variables.

		Args:
			sim_data_file: path to the sim_data file for the simulation
			sim_out_dir: path to the simOut dir for the simulation
		"""

		with open(sim_data_file, 'rb') as f:
			self.sim_data = pickle.load(f)
		self.aa_from_trna = self.sim_data.process.transcription.aa_from_trna
		self.charging_params = get_charging_params(self.sim_data,
			variable_elongation=self.variable_elongation)

		# Listeners used
		growth_reader = TableReader(os.path.join(sim_out_dir, 'GrowthLimits'))
		main_reader = TableReader(os.path.join(sim_out_dir, 'Main'))

		# Load data
		self.synthetase_conc = CONC_UNITS * growth_reader.readColumn('synthetase_conc')[1:, :]
		self.uncharged_trna_conc = CONC_UNITS * growth_reader.readColumn('uncharged_trna_conc')[1:, :]
		self.charged_trna_conc = CONC_UNITS * growth_reader.readColumn('charged_trna_conc')[1:, :]
		self.aa_conc = CONC_UNITS * growth_reader.readColumn('aa_conc')[1:, :]
		self.ribosome_conc = CONC_UNITS * growth_reader.readColumn('ribosome_conc')[1:]
		self.fraction_aa_to_elongate = growth_reader.readColumn('fraction_aa_to_elongate')[1:, :]
		self.charging_fraction = growth_reader.readColumn('fraction_trna_charged')[1:, :]

		self.time_step_sizes = main_reader.readColumn('timeStepSec')[1:]
		self.n_time_steps = len(self.time_step_sizes)

	def solve_timestep(self, timestep: int) -> Tuple[np.ndarray, float]:
		"""
		Calculates charging and elongation rate for a given timestep.

		Args:
			timestep: simulation timestep to select data from
		"""

		return calculate_trna_charging(
			self.synthetase_conc[timestep, :],
			self.uncharged_trna_conc[timestep, :],
			self.charged_trna_conc[timestep, :],
			self.aa_conc[timestep, :],
			self.ribosome_conc[timestep],
			self.fraction_aa_to_elongate[timestep, :],
			self.charging_params,
			time_limit=self.time_step_sizes[timestep],
			)

	def validation(self, n_steps: int) -> None:
		"""
		Performs a validation check to makes sure solving the model from
		loaded data matches the objective from the original solution during
		the simulation.

		Args:
			n_steps: number of timesteps to check
				if 0: does not check
				if <0: runs all timepoints from the simulation
		"""

		if n_steps == 0:
			return
		elif n_steps < 0:
			n_steps = self.n_time_steps
		else:
			n_steps = min(self.n_time_steps, n_steps)

		print('Running validation to check output...')
		for timestep in range(n_steps):
			charging_fraction, _ = self.solve_timestep(timestep)
			if np.any(charging_fraction @ self.aa_from_trna != self.charging_fraction[timestep]):
				raise ValueError(f'Charging fraction does not match for time step {timestep}')
		print('All {} timesteps match the results from the whole-cell model.'.format(n_steps))

	def run(self, args: argparse.Namespace) -> None:
		self.variable_elongation = args.variable_elongation_translation
		self.load_data(args.sim_data_file, args.sim_out_dir)
		self.validation(args.validation)


if __name__ == '__main__':
	analysis = ChargingDebug()
	analysis.cli()
