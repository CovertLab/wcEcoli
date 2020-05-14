#! /usr/bin/env python
"""
Run a parameter search using simulation output.

Prerequisite: Run the parameter calculator (runParca.py).

TODO: Share more code with fw_queue.py.

Run with '-h' for command line help.
Set PYTHONPATH when running this.
"""

from __future__ import absolute_import, division, print_function

import cPickle
import re
import os
from typing import Tuple

from wholecell.fireworks.firetasks import SimulationTask, VariantSimDataTask
from wholecell.utils import constants, data, scriptBase
import wholecell.utils.filepath as fp


SIM_DIR_PATTERN = r'({})__(.+)'.format(fp.TIMESTAMP_PATTERN)


# TODO: make base class to inherit from
class ppGpp():
	parameters = ('a', 'b')
	learning_rate = 0.001
	parameter_step = 0.9

	def get_parameter(self, iteration):
		idx = iteration % len(self.parameters)
		return self.parameters[idx]

	def update_parameter_value(self, sim_data, parameter):
		old_value = sim_data[parameter]
		new_value = old_value * self.parameter_step
		diff = new_value - old_value

		return new_value, diff

	def perturb_sim_data(self, sim_data_file, parameter):
		# with open(sim_data_file) as f:
		# 	sim_data = cPickle.load(f)
		sim_data = sim_data_file.copy()

		sim_data[parameter] = self.update_parameter_value(sim_data, parameter)[0]

		return sim_data

	def get_parameter_update(self, sim_data_file, parameter, old_objective, new_objective):
		# with open(sim_data_file) as f:
		# 	sim_data = cPickle.load(f)
		sim_data = sim_data_file.copy()

		update = self.learning_rate * (new_objective - old_objective) / self.update_parameter_value(sim_data, parameter)[1]

		return update

	def update_sim_data(self, sim_data_file, updates):
		# with open(sim_data_file) as f:
		# 	sim_data = cPickle.load(f)
		sim_data = sim_data_file.copy()

		for parameter, update in zip(self.parameters, updates):
			sim_data[parameter] -= update

		return sim_data

	def get_objective_value(self, sim_data_file, sim_out_dir):
		# with open(sim_data_file) as f:
		# 	sim_data = cPickle.load(f)

		# Load listeners
		y0 = 1
		y = 10*sim_data_file['a'] + 5*sim_data_file['a']*sim_data_file['b']
		objective = (y - y0)**2

		return objective

# TODO: move to another file
PARAMETER_METHODS = {
	'ppGpp': ppGpp,
	}


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


class RunSimulation(scriptBase.ScriptBase):
	"""Drives a simple simulation run."""

	def description(self):
		"""Describe the command line program."""
		return 'Whole Cell E. coli simulation'

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
		super(RunSimulation, self).define_parameters(parser)
		self.define_parameter_sim_dir(parser)
		self.define_sim_loop_options(parser, manual_script=True)
		self.define_sim_options(parser)
		self.define_elongation_options(parser)

		parser.add_argument('--method',
			required=True,
			choices=PARAMETER_METHODS.keys(),
			help='Method for updating parameters in parameter search.')
		parser.add_argument('--iterations',
			default=5,
			type=int,
			help='Number of iterations to update parameters before stopping.')

	def run(self, args):
		kb_directory = os.path.join(args.sim_path, 'kb')
		sim_data_file = os.path.join(kb_directory, constants.SERIALIZED_SIM_DATA_FILENAME)
		fp.verify_file_exists(sim_data_file, 'Run runParca?')

		timestamp, description = parse_timestamp_description(args.sim_path)

		variant_type = args.variant[0]
		variant_spec = (variant_type, int(args.variant[1]), int(args.variant[2]))

		cli_sim_args = data.select_keys(vars(args), scriptBase.SIM_KEYS)

		# Write the metadata file.
		metadata = data.select_keys(
			vars(args),
			scriptBase.METADATA_KEYS,
			git_hash=fp.run_cmdline("git rev-parse HEAD") or '--',
			git_branch=fp.run_cmdline("git symbolic-ref --short HEAD") or '--',
			description=description,
			time=timestamp,
			analysis_type=None,
			variant=variant_type,
			total_variants=str(variant_spec[2] + 1 - variant_spec[1]),
			total_gens=args.total_gens or args.generations)
		metadata_dir = fp.makedirs(args.sim_path, 'metadata')
		metadata_path = os.path.join(metadata_dir, constants.JSON_METADATA_FILE)
		fp.write_json_file(metadata_path, metadata)

		method = PARAMETER_METHODS[args.method]()
		sim_data = {'a': 1, 'b': 1}  # TODO: handle sim_data file
		for i in range(args.iterations):
			sim_out_dir = None  # TODO: run sim, handle out dirs
			objective = method.get_objective_value(sim_data, sim_out_dir)

			# Check all parameters before updating objective for faster convergence
			updates = []
			for parameter in method.parameters:
				perturbed_sim_data = method.perturb_sim_data(sim_data, parameter)

				sim_out_dir = None  # TODO: run sim, handle out dirs

				new_objective = method.get_objective_value(perturbed_sim_data, sim_out_dir)
				updates.append(method.get_parameter_update(sim_data, parameter, objective, new_objective))

			sim_data = method.update_sim_data(sim_data, updates)
			print('It {}: {}\n\t{}'.format(i, new_objective, sim_data))

		import ipdb; ipdb.set_trace()


		# args.sim_path is called INDIV_OUT_DIRECTORY in fw_queue.
		for i, subdir in fp.iter_variants(*variant_spec):
			variant_directory = os.path.join(args.sim_path, subdir)
			variant_sim_data_directory = os.path.join(variant_directory,
				VariantSimDataTask.OUTPUT_SUBDIR_KB)

			variant_sim_data_modified_file = os.path.join(
				variant_sim_data_directory, constants.SERIALIZED_SIM_DATA_MODIFIED)

			if args.require_variants:
				fp.verify_file_exists(
					variant_sim_data_modified_file, 'Run makeVariants?')
			else:
				variant_metadata_directory = os.path.join(variant_directory,
					VariantSimDataTask.OUTPUT_SUBDIR_METADATA)
				task = VariantSimDataTask(
					variant_function=variant_type,
					variant_index=i,
					input_sim_data=sim_data_file,
					output_sim_data=variant_sim_data_modified_file,
					variant_metadata_directory=variant_metadata_directory,
					)
				task.run_task({})

			for j in xrange(args.seed, args.seed + args.init_sims):  # init sim seeds
				seed_directory = fp.makedirs(variant_directory, "%06d" % j)

				for k in xrange(args.generations):  # generation number k
					gen_directory = fp.makedirs(seed_directory, "generation_%06d" % k)

					# l is the daughter number among all of this generation's cells,
					# which is 0 for single-daughters but would span range(2**k) if
					# each parent had 2 daughters.
					l = 0
					cell_directory = fp.makedirs(gen_directory, "%06d" % l)
					cell_sim_out_directory = fp.makedirs(cell_directory, "simOut")

					options = dict(cli_sim_args,
						input_sim_data=variant_sim_data_modified_file,
						output_directory=cell_sim_out_directory,
						)

					if k == 0:
						task = SimulationTask(seed=j, **options)
					else:
						parent_gen_directory = os.path.join(seed_directory, "generation_%06d" % (k - 1))
						parent_cell_directory = os.path.join(parent_gen_directory, "%06d" % (l // 2))
						parent_cell_sim_out_directory = os.path.join(parent_cell_directory, "simOut")
						daughter_state_path = os.path.join(
							parent_cell_sim_out_directory,
							constants.SERIALIZED_INHERITED_STATE % (l % 2 + 1))
						task = SimulationDaughterTask(
							seed=(j + 1) * ((2 ** k - 1) + l),
							inherited_state_path=daughter_state_path,
							**options
							)
					task.run_task({})


if __name__ == '__main__':
	script = RunSimulation()
	script.cli()
