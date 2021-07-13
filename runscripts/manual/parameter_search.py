#! /usr/bin/env python
"""
Run a parameter search using parca and simulation output.  See methods and
solvers in models/ecoli/sim/parameter_search and wholecell/optimization,
respectively.

TODO: Share more code with fw_queue.py and other runscripts.

Run with '-h' for command line help.
Set PYTHONPATH when running this.
"""

import os
import re
from typing import Tuple

from models.ecoli.sim.parameter_search import PARAMETER_METHODS
from wholecell.optimization import SOLVERS
from wholecell.utils import constants, scriptBase
import wholecell.utils.filepath as fp


SIM_DIR_PATTERN = r'({})__(.+)'.format(fp.TIMESTAMP_PATTERN)

# Command line arg defaults for solver options
DEFAULT_ITERATIONS = 5
DEFAULT_STARTING_ITERATION = 0
DEFAULT_LEARNING_RATE = 1
DEFAULT_PARAMETER_STEP = 0.1
DEFAULT_MAX_CHANGE = 0.1
DEFAULT_ALPHA = 0.1
DEFAULT_GAMMA = 0.1


def parse_timestamp_description(sim_path):
	# type: (str) -> Tuple[str, str]
	"""Parse `timestamp, description` from a sim_path that ends with a dir like
	'20190704.101500__Latest_sim_run' or failing that, return defaults.
	"""
	sim_dir = os.path.basename(sim_path)
	if not sim_dir:  # sim_path is empty or ends with '/'
		sim_dir = os.path.basename(os.path.dirname(sim_path))

	match = re.match(SIM_DIR_PATTERN, sim_dir)
	if match:
		timestamp = match.group(1)
		description = match.group(2).replace('_', ' ')
	else:
		timestamp = fp.timestamp()
		description = 'a manual run'

	return timestamp, description


class RunParameterSearch(scriptBase.ScriptBase):
	"""Runs various parameter search algorithms to optimize desired objectives."""

	def description(self):
		"""Describe the command line program."""
		return 'Whole Cell E. coli parameter search'

	def help(self):
		"""Return help text for the Command Line Interface."""
		return '''Run a {}.
				If the sim_path ends with a dir like
				"20190704.101500__Latest_sim_run", this will get the
				timestamp and description from the path to write into
				metadata.json.
				The command line option names are long but you can use any
				unambiguous prefix.'''.format(self.description())

	def define_parameters(self, parser):
		super().define_parameters(parser)

		# TODO: reuse params from other parsers?
		# self.define_parameter_sim_dir(parser)
		# self.define_sim_loop_options(parser, manual_script=True)
		# self.define_sim_options(parser)
		# self.define_elongation_options(parser)

		default_solver = list(SOLVERS.keys())[0]

		# NOTE: Don't name this arg sim_dir since that makes parse_args() look
		# for an existing sim_dir directory while here we aim to create one.
		parser.add_argument('sim_outdir', nargs='?', default='manual',
			help='The simulation "out/" subdirectory to write to.'
				 ' Default = "manual".')
		parser.add_argument('--timestamp', action='store_true',
			help='Timestamp the given `sim_outdir`, transforming e.g.'
				 ' "Fast run" to "20190514.135600__Fast_run".')

		parser.add_argument('--solver',
			default=default_solver,
			choices=SOLVERS.keys(),
			help=f'Solver for optimizing parameters (default: {default_solver}).')
		parser.add_argument('--method',
			required=True,
			choices=PARAMETER_METHODS.keys(),
			help='Class defining parameters and an objective for parameter search.')
		parser.add_argument('--iterations',
			default=DEFAULT_ITERATIONS,
			type=int,
			help=f'Number of iterations to update parameters before stopping (default: {DEFAULT_ITERATIONS}).')
		parser.add_argument('--starting-iteration',
			default=DEFAULT_STARTING_ITERATION,
			type=int,
			help=f'Number of iterations to update parameters before stopping (default: {DEFAULT_STARTING_ITERATION}).')
		parser.add_argument('--learning-rate',
			default=DEFAULT_LEARNING_RATE,
			type=float,
			help=f'Learning rate for updating parameters (default: {DEFAULT_LEARNING_RATE}).')
		parser.add_argument('--parameter-step',
			default=DEFAULT_PARAMETER_STEP,
			type=float,
			help=f'Fraction to update parameters by to determine the gradient (default: {DEFAULT_PARAMETER_STEP}).')
		parser.add_argument('--max-change',
			default=DEFAULT_MAX_CHANGE,
			type=float,
			help=f'Maximum fraction to update a parameter in a given time step (default: {DEFAULT_MAX_CHANGE}).')
		parser.add_argument('--cpus',
			default=1,
			type=int,
			help='Number of CPUs to use for running sims in parallel for a given iteration.')

		# SPSA options
		parser.add_argument('--alpha',
			type=float,
			default=DEFAULT_ALPHA,
			help=f'Set alpha parameter for SPSA solver to control learning rate decay (default: {DEFAULT_ALPHA}).')
		parser.add_argument('--gamma',
			type=float,
			default=DEFAULT_GAMMA,
			help=f'Set gamma parameter for SPSA solver to control parameter update decay (default: {DEFAULT_GAMMA}).')

	def parse_args(self):
		args = super().parse_args()

		args.time = fp.timestamp()
		args.description = args.sim_outdir.replace(' ', '_')

		if args.timestamp:
			args.sim_outdir = args.time + '__' + args.description

		args.sim_path = fp.makedirs(fp.ROOT_PATH, "out", args.sim_outdir)
		return args

	def run(self, args):
		fp.makedirs(args.sim_path, constants.KB_DIR)
		fp.makedirs(args.sim_path, 'metadata')  # TODO: write metadata?

		method = PARAMETER_METHODS[args.method]()
		solver = SOLVERS[args.solver](method, args)

		# TODO: address this issue by saving raw/sim data directly or saving the updates to apply
		if args.starting_iteration != 0:
			print('Warning: starting at a new iteration is not guaranteed to give identical results.'
				' The updated parameters should be saved after each iteration to restart at the same'
				' place in the algorithm.')

		for it in range(args.starting_iteration, args.starting_iteration + args.iterations):
			print(f'** Starting iteration {it} **')
			objective = solver.run_iteration()
			solver.print_update(objective)


if __name__ == '__main__':
	script = RunParameterSearch()
	script.cli()
