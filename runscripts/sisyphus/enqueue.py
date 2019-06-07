#!/usr/bin/env python

"""
Create a workflow DAG of Sisyphus tasks and specifies their links and enqueue
it to run.
"""

from __future__ import absolute_import, division, print_function

import argparse
import pprint as pp
from typing import Any, Dict, List

from runscripts.sisyphus.workflow import Task, Workflow


DOCKER_IMAGE = 'gcr.io/allen-discovery-center-mcovert/wcm-code:latest'
STORAGE_PREFIX = 'sisyphus:data/'

VARIANT_OPTION = 'variant'  # this unique option takes 3 args

# CLI options for the workflow tasks.
PARCA_OPTIONS = {'ribosome_fitting', 'rnapoly_fitting'}
MAKE_VARIANT_OPTIONS = {VARIANT_OPTION}
SIM_OPTIONS = {VARIANT_OPTION, 'generations', 'total_gens', 'seed', 'init_sims',
	'timeline', 'length_sec', 'timestep_safety_frac', 'timestep_max',
	'timestep_update_freq', 'mass_distribution', 'growth_rate_noise',
	'd_period_division', 'translation_supply', 'trna_charging'}


def _add_options(tokens, args, *options):
	# type: (List[str], Dict[str, Any], *str) -> None
	"""Add the named options from `args` to command line `tokens`."""
	for option_name in options:
		value = args.get(option_name)
		if value is not None:
			if option_name == VARIANT_OPTION:
				# `value` is [VARIANT_TYPE, FIRST_INDEX, LAST_INDEX]
				assert int(value[1]) <= int(value[2])
				tokens += ['--' + option_name] + [str(v) for v in value]
			else:
				tokens += ['--' + option_name, value]

class WCM_Workflow(Workflow):
	"""A Workflow builder for the wcEcoli Whole Cell Model."""

	def __init__(self):
		# type: () -> None
		super(WCM_Workflow, self).__init__()

		# TODO(jerry): A developer-specific Docker image with their code.
		self.image = DOCKER_IMAGE

		# TODO(jerry): A developer- and experiment-specific storage_prefix.
		self.storage_prefix = STORAGE_PREFIX

		self.sim_outdir = 'workflow'  # subdir of './out' in the container

	def add_python_task(self, python_args, *upstream_tasks, **kwargs):
		# type: (List[str], *Task, **Any) -> Task
		"""Operator: Add a Python task to the workflow and return it."""
		kwargs['image'] = self.image

		# Sisyphus needs `-u` unbuffered output.
		kwargs['commands'] = [{'command': ['python', '-u'] + python_args}]

		return self.add_task(Task(*upstream_tasks, **kwargs))

	def add_parca_task(self, *upstream_tasks, **kwargs):
		# type: (*Task, **Any) -> Task
		"""Operator: Add a parameter-calc task to the workflow and return it."""
		kwargs.setdefault('name', 'parca')

		tokens = ['runscripts/manual/runParca.py', self.sim_outdir]
		_add_options(tokens, kwargs, *PARCA_OPTIONS)

		# Use `upstream_tasks` to set the inputs, then add the outputs.
		task = self.add_python_task(tokens, *upstream_tasks, **kwargs)
		task.set_output_mapping(self.storage_prefix, 'kb/')
		return task

	def add_variants_task(self, *upstream_tasks, **kwargs):
		# type: (*Task, **Any) -> Task
		"""Operator: Add a make-variants task to the workflow and return it."""
		kwargs.setdefault('name', 'makeVariants')

		tokens = ['runscripts/manual/makeVariants.py', self.sim_outdir]
		_add_options(tokens, kwargs, *MAKE_VARIANT_OPTIONS)

		# TODO(jerry): Compute the actual variant paths.
		task = self.add_python_task(tokens, *upstream_tasks, **kwargs)
		task.set_output_mapping(self.storage_prefix, 'wildtype_000000', 'kb/')
		task.set_output_mapping(self.storage_prefix, 'wildtype_000000', 'metadata/')
		return task

	def add_sim_task(self, *upstream_tasks, **kwargs):
		# type: (*Task, **Any) -> Task
		"""Operator: Add a run-sim task to the workflow and return it."""
		kwargs.setdefault('name', 'runSim')

		tokens = [
			'runscripts/manual/runSim.py', self.sim_outdir,
			'--require_variants', '1']
		_add_options(tokens, kwargs, *SIM_OPTIONS)

		# TODO(jerry): Compute the actual simOut paths.
		task = self.add_python_task(tokens, *upstream_tasks, **kwargs)
		task.set_output_mapping(self.storage_prefix, 'wildtype_000000',
			'000000', 'generation_000000', '000000', 'simOut/')
		return task


def wc_ecoli_workflow(args):
	# type: (Dict[str, Any]) -> WCM_Workflow
	"""Construct a workflow DAG for wcEcoli."""
	wf = WCM_Workflow()

	t_parca = wf.add_parca_task(**args)
	t_variants = wf.add_variants_task(t_parca, **args)

	# TODO(jerry): Loop over variants, seeds, & gens, adding sim & sim daughter
	#  tasks. Change args in kwargs so runSim will only run one iteration and
	#  get the seed, etc. Setup each sim's upstream tasks and simOut path.
	t_sim = wf.add_sim_task(t_variants, **args)

	# TODO(jerry): Add the analysis tasks...

	return wf

def enqueue_workflow(args):
	# type: (Dict[str, Any]) -> None
	"""Construct and enqueue a workflow."""
	wf = wc_ecoli_workflow(args)
	wf.enqueue()


class Default(object):
	def __repr__(self):
		return 'default'


def default(v):
	return Default() if v is None else v


def cli():
	"""Command line interpreter."""
	# This simple CLI parser passes the provided args to the runscripts,
	# letting their CLIs handle parameter types and default values.
	parser = argparse.ArgumentParser(description='''
		Construct a wcEcoli Whole Cell Model workflow and [TODO] queue it up.
		See runscripts/manual/*.py for help on the arguments and their defaults.
		''')
	parser.add_argument('--' + VARIANT_OPTION, nargs=3,
		metavar=('VARIANT_TYPE', 'FIRST_INDEX', 'LAST_INDEX'))

	options = (PARCA_OPTIONS | MAKE_VARIANT_OPTIONS | SIM_OPTIONS) - {
		VARIANT_OPTION}
	for option in options:
		parser.add_argument('--' + option)

	args = vars(parser.parse_args())
	pp.pprint({'Arguments': {k: default(v) for k, v in args.viewitems()}})
	print()

	enqueue_workflow(args)


if __name__ == '__main__':
	cli()
