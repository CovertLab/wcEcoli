#! /usr/bin/env python
"""
Run a parameter search using simulation output.

Prerequisite: Run the parameter calculator (runParca.py).

TODO: Share more code with fw_queue.py.

Run with '-h' for command line help.
Set PYTHONPATH when running this.
"""

import os
import pickle
import re
from typing import Tuple

import numpy as np

from wholecell.analysis.analysis_tools import read_bulk_molecule_counts
from wholecell.fireworks.firetasks import SimulationTask, SimulationDaughterTask, VariantSimDataTask
from wholecell.io.tablereader import TableReader
from wholecell.utils import constants, data, scriptBase, units
import wholecell.utils.filepath as fp


SIM_DIR_PATTERN = r'({})__(.+)'.format(fp.TIMESTAMP_PATTERN)

# Command line arg defaults for solver options
DEFAULT_ITERATIONS = 5
DEFAULT_LEARNING_RATE = 1e-3
DEFAULT_PARAMETER_STEP = 1e-2


class ParameterSearch(object):
	parameters = ()

	def __init__(self, sim_dir, lr=1e-3, step=0.99):
		if len(self.parameters) == 0:
			raise NotImplementedError('Must define parameters to search in a subclass.')

		self.learning_rate = lr
		self.parameter_step = step
		self.sim_dir = sim_dir

	@staticmethod
	def reader(out_dir, table):
		return TableReader(os.path.join(out_dir, table))

	def sim_data_path(self, variant):
		kb_dir = fp.makedirs(self.sim_dir, '{}_{:06n}'.format(
			self.__class__.__name__, variant), VariantSimDataTask.OUTPUT_SUBDIR_KB)
		return os.path.join(kb_dir, constants.SERIALIZED_SIM_DATA_MODIFIED)

	def parameter_summary(self, sim_data):
		summary = {}
		for p in self.parameters:
			summary['sim_data.{}'.format(p)] = self.get_attrs(sim_data, p)
		return summary

	@staticmethod
	def get_attrs(sim_data, attr):
		val = sim_data
		for a in attr.split('.'):
			val = getattr(val, a)
		return val

	@staticmethod
	def set_attrs(sim_data, attr, value):
		parent = sim_data
		attrs = attr.split('.')
		for a in attrs[:-1]:
			parent = getattr(parent, a)
		setattr(parent, attrs[-1], value)

	def update_parameter_value(self, sim_data, parameter):
		old_value = self.get_attrs(sim_data, parameter)
		new_value = old_value * (1 - self.parameter_step)
		diff = new_value - old_value

		return new_value, diff

	def perturb_sim_data(self, sim_data_file, parameter):
		with open(sim_data_file, 'rb') as f:
			sim_data = pickle.load(f)

		value = self.update_parameter_value(sim_data, parameter)[0]
		self.set_attrs(sim_data, parameter, value)

		return sim_data

	def perturb_sim_data_spsa(self, sim_data_file, iteration, alpha, gamma, direction):
		with open(sim_data_file, 'rb') as f:
			sim_data = pickle.load(f)

		_, deltas = self.get_spsa_parameters(iteration, alpha, gamma)
		for d, parameter in zip(deltas, self.parameters):
			value = self.get_attrs(sim_data, parameter) * (1 + direction * d)
			self.set_attrs(sim_data, parameter, value)

		return sim_data

	def get_parameter_update(self, sim_data_file, parameter, old_objective, new_objective):
		with open(sim_data_file, 'rb') as f:
			sim_data = pickle.load(f)

		update = (self.learning_rate * (new_objective - old_objective)
			/ self.update_parameter_value(sim_data, parameter)[1])

		return update

	def update_sim_data(self, sim_data_file, updates):
		with open(sim_data_file, 'rb') as f:
			sim_data = pickle.load(f)

		for parameter, update in zip(self.parameters, updates):
			old_val = self.get_attrs(sim_data, parameter)
			if units.hasUnit(old_val):
				unit = units.getUnit(old_val)
				old_val = old_val.asNumber(unit)
				update = update.asNumber(1/unit)
			else:
				unit = 1
			new_val = unit * (old_val - update)
			self.set_attrs(sim_data, parameter, new_val)

		return sim_data

	def update_sim_data_spsa(self, sim_data_file, objectives, iteration, alpha, gamma):
		with open(sim_data_file, 'rb') as f:
			sim_data = pickle.load(f)

		objective_difference = objectives[1] - objectives[0]

		at, deltas = self.get_spsa_parameters(iteration, alpha, gamma)
		updates = []
		# TODO: iterate at if different learning rate for each parameter
		for d, parameter in zip(deltas, self.parameters):
			updates.append(at * objective_difference / (2 * d * self.get_attrs(sim_data, parameter)))

		return self.update_sim_data(sim_data_file, updates)

	def get_spsa_parameters(self, iteration, alpha, gamma):
		it = iteration + 1
		np.random.seed(it)
		at = self.learning_rate / it**alpha
		ct = self.parameter_step / it**gamma
		deltas = ct * (np.random.rand(len(self.parameters)) * 2 - 1)

		return at, deltas

	def get_objective_value(self, sim_data_file, sim_out_dir):
		raise NotImplementedError('Need to implement this in a subclass.')

class BasicppGpp(ParameterSearch):
	parameters = ('constants.KD_RelA_ribosome', 'constants.k_RelA_ppGpp_synthesis')

	def get_objective_value(self, sim_data_file, sim_out_dir):
		with open(sim_data_file, 'rb') as f:
			sim_data = pickle.load(f)

		# Listeners used
		enzyme_kinetics_reader = self.reader(sim_out_dir, 'EnzymeKinetics')

		# Load data
		counts_to_molar = enzyme_kinetics_reader.readColumn('countsToMolar')
		ppgpp_counts, = read_bulk_molecule_counts(sim_out_dir, ([sim_data.moleculeIds.ppGpp],))
		ppgpp_conc = (counts_to_molar * ppgpp_counts)[1:].mean() * 1000

		# Calculate objective
		ppgpp_target = 45
		objective = (ppgpp_conc - ppgpp_target)**2

		return objective

class ppGpp(ParameterSearch):
	parameters = (
		'constants.synthetase_charging_rate',
		'constants.Km_synthetase_uncharged_trna',
		'constants.Km_synthetase_amino_acid',
		'constants.Kdissociation_charged_trna_ribosome',
		'constants.Kdissociation_uncharged_trna_ribosome',
		'constants.KD_RelA_ribosome',
		'constants.k_RelA_ppGpp_synthesis',
		'constants.k_SpoT_ppGpp_synthesis',
		'constants.k_SpoT_ppGpp_degradation',
		'constants.KI_SpoT_ppGpp_degradation',
		)

	def get_objective_value(self, sim_data_file, sim_out_dir):
		with open(sim_data_file, 'rb') as f:
			sim_data = pickle.load(f)

		# Listeners used
		enzyme_kinetics_reader = self.reader(sim_out_dir, 'EnzymeKinetics')

		# Load data
		counts_to_molar = enzyme_kinetics_reader.readColumn('countsToMolar')
		ppgpp_counts, = read_bulk_molecule_counts(sim_out_dir, ([sim_data.moleculeIds.ppGpp],))
		ppgpp_conc = (counts_to_molar * ppgpp_counts)[1:] * 1000

		# Calculate objective
		ppgpp_target = 45
		objective = np.linalg.norm(ppgpp_conc - ppgpp_target)

		return objective

# TODO: move to another file
PARAMETER_METHODS = {
	'ppGpp': ppGpp,
	'BasicppGpp': BasicppGpp,
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

def run_sim(cli_args, variant):
	cli_sim_args = data.select_keys(vars(cli_args), scriptBase.SIM_KEYS)

	variant_type = cli_args.method
	variant_directory = os.path.join(cli_args.sim_path, '{}_{:06n}'
		.format(variant_type, variant))
	variant_sim_data_directory = os.path.join(variant_directory,
		VariantSimDataTask.OUTPUT_SUBDIR_KB)

	variant_sim_data_modified_file = os.path.join(
		variant_sim_data_directory, constants.SERIALIZED_SIM_DATA_MODIFIED)

	for j in range(cli_args.seed,
			cli_args.seed + cli_args.init_sims):  # init sim seeds
		seed_directory = fp.makedirs(variant_directory, "%06d" % j)

		for k in range(cli_args.generations):  # generation number k
			gen_directory = fp.makedirs(seed_directory,
				"generation_%06d" % k)

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
				parent_gen_directory = os.path.join(seed_directory,
					"generation_%06d" % (k - 1))
				parent_cell_directory = os.path.join(parent_gen_directory,
					"%06d" % (l // 2))
				parent_cell_sim_out_directory = os.path.join(
					parent_cell_directory, "simOut")
				daughter_state_path = os.path.join(
					parent_cell_sim_out_directory,
					constants.SERIALIZED_INHERITED_STATE % (l % 2 + 1))
				task = SimulationDaughterTask(
					seed=(j + 1) * ((2 ** k - 1) + l),
					inherited_state_path=daughter_state_path,
					**options
					)
			task.run_task({})

			# TODO: currently only supports one sim cycle with out dir
			return cell_sim_out_directory

def gradient_descent(method, args, n_variants, sim_data_file, iteration):
	sim_out_dir = run_sim(args, n_variants)
	n_variants += 1

	objective = method.get_objective_value(sim_data_file, sim_out_dir)

	print(f'Initial objective: {objective:.3f}\n')

	# Check all parameters before updating objective for faster convergence
	updates = []
	for parameter in method.parameters:
		# Perturb parameter in sim_data
		perturbed_sim_data = method.perturb_sim_data(sim_data_file, parameter)

		# Save perturbed sim_data for variant sim
		perturbed_sim_data_file = method.sim_data_path(n_variants)
		with open(perturbed_sim_data_file, 'wb') as f:
			pickle.dump(perturbed_sim_data, f, protocol=pickle.HIGHEST_PROTOCOL)

		# Run sim with perturbed sim_data
		sim_out_dir = run_sim(args, n_variants)
		n_variants += 1

		# Calculate objective and resulting parameter update
		new_objective = method.get_objective_value(perturbed_sim_data_file, sim_out_dir)
		print(f'Updated {parameter}: objective = {new_objective:.3f}\n')
		updates.append(method.get_parameter_update(sim_data_file, parameter, objective, new_objective))

	# Apply all updates to sim_data
	sim_data = method.update_sim_data(sim_data_file, updates)

	return sim_data, n_variants

def spsa(method, args, n_variants, sim_data_file, iteration):
	"""Simultaneous perturbation stochastic approximation"""

	# TODO: pass these in as args and allow adjustment from the command line
	alpha = 0.1
	gamma = 0.1

	objectives = []
	for direction in [-1, 1]:
		# Perturb parameter in sim_data
		perturbed_sim_data = method.perturb_sim_data_spsa(sim_data_file, iteration, alpha, gamma, direction)

		# Save perturbed sim_data for variant sim
		perturbed_sim_data_file = method.sim_data_path(n_variants)
		with open(perturbed_sim_data_file, 'wb') as f:
			pickle.dump(perturbed_sim_data, f, protocol=pickle.HIGHEST_PROTOCOL)

		# Run sim with perturbed sim_data
		sim_out_dir = run_sim(args, n_variants)
		n_variants += 1

		# Calculate objective and resulting parameter update
		new_objective = method.get_objective_value(perturbed_sim_data_file, sim_out_dir)
		print(f'Updated parameter direction {direction}: objective = {new_objective:.3f}\n')
		objectives.append(new_objective)

	# Apply all updates to sim_data
	sim_data = method.update_sim_data_spsa(sim_data_file, objectives, iteration, alpha, gamma)

	return sim_data, n_variants


SOLVERS = {
	'spsa': spsa,
	'gradient-descent': gradient_descent,
	}


class RunParameterSearch(scriptBase.ScriptBase):
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
		super(RunParameterSearch, self).define_parameters(parser)
		self.define_parameter_sim_dir(parser)
		self.define_sim_loop_options(parser, manual_script=True)
		self.define_sim_options(parser)
		self.define_elongation_options(parser)

		default_solver = list(SOLVERS.keys())[0]

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
		parser.add_argument('--learning-rate',
			default=DEFAULT_LEARNING_RATE,
			type=float,
			help=f'Learning rate for updating parameters (default: {DEFAULT_LEARNING_RATE}).')
		parser.add_argument('--parameter-step',
			default=DEFAULT_PARAMETER_STEP,
			type=float,
			help=f'Fraction to update parameters by to determine the gradient (default: {DEFAULT_PARAMETER_STEP}).')


	def run(self, args):
		kb_directory = os.path.join(args.sim_path, 'kb')
		sim_data_file = os.path.join(kb_directory, constants.SERIALIZED_SIM_DATA_FILENAME)
		with open(sim_data_file, 'rb') as f:
			sim_data = pickle.load(f)
		fp.verify_file_exists(sim_data_file, 'Run runParca?')

		timestamp, description = parse_timestamp_description(args.sim_path)

		variant_type = args.variant[0]
		variant_spec = (variant_type, int(args.variant[1]), int(args.variant[2]))

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
			total_variants=str(variant_spec[2] + 1 - variant_spec[1]),  # TODO: arg for total iters
			total_gens=args.total_gens or args.generations)
		metadata_dir = fp.makedirs(args.sim_path, 'metadata')
		metadata_path = os.path.join(metadata_dir, constants.JSON_METADATA_FILE)
		fp.write_json_file(metadata_path, metadata)

		solver = SOLVERS[args.solver]
		method = PARAMETER_METHODS[args.method](args.sim_path,
			lr=args.learning_rate, step=args.parameter_step)
		n_variants = 0
		sim_data_file = method.sim_data_path(n_variants)
		with open(sim_data_file, 'wb') as f:
			pickle.dump(sim_data, f, protocol=pickle.HIGHEST_PROTOCOL)
		for i in range(args.iterations):
			print(f'** Starting iteration {i} **')

			sim_data, n_variants = solver(method, args, n_variants, sim_data_file, i)

			# Save updated sim_data
			sim_data_file = method.sim_data_path(n_variants)
			n_variants += 1
			with open(sim_data_file, 'wb') as f:
				pickle.dump(sim_data, f, protocol=pickle.HIGHEST_PROTOCOL)

			# Print status update
			print('New parameter values after this iteration:')
			for p, val in sorted(method.parameter_summary(sim_data).items()):
				print(f'   {p}: {val}')
			print('')


if __name__ == '__main__':
	script = RunParameterSearch()
	script.cli()
