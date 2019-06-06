#!/usr/bin/env python

"""
Create a workflow DAG of Sisyphus tasks and specifies their links and enqueue
it to run.
"""

from __future__ import absolute_import, division, print_function

from os import path as op
from typing import Any, Dict, List

from .workflow import Task, Workflow


DOCKER_IMAGE = 'gcr.io/allen-discovery-center-mcovert/wcm-code:latest'
STORAGE_PREFIX = 'sisyphus:data/'


def _add_options(tokens, args, *options):
	# type: (List[str], Dict[str, Any], *str) -> None
	"""Add the named options from `args` to command line `tokens`."""
	for option_name in options:
		value = args.get(option_name)
		if value:
			if option_name != 'variant':
				value = (value,)
				# else `value` is (VARIANT_TYPE, FIRST_INDEX, LAST_INDEX)
			tokens += ['--' + option_name] + [str(v) for v in value]

class WCM_Workflow(Workflow):
	"""A Workflow builder for the wcEcoli Whole Cell Model."""

	# TODO(jerry): Finish adding the I/O mappings.

	def __init__(self):
		# type: () -> None
		super(WCM_Workflow, self).__init__()

		# TODO(jerry): A developer-specific Docker image with their code.
		self.image = DOCKER_IMAGE

		# TODO(jerry): A developer- and experiment-specific storage_prefix.
		self.storage_prefix = STORAGE_PREFIX

		self.sim_outdir = 'workflow'  # subdir of './out' in the container

	def _set_path_mapping(self, args, key, *subdirs):
		# type: (Dict[str, Any], str, *str) -> None
		"""Set an I/O path mapping for the given `subdirs` of `sim_outdir`."""
		path = op.join(self.sim_outdir, *subdirs)
		args.setdefault(key, {})[op.join(self.storage_prefix, path)] = path

	def _set_input_mapping(self, args, *subdirs):
		# type: (Dict[str, Any], *str) -> None
		"""Set an input path mapping for the given `subdirs` of `sim_outdir`."""
		self._set_path_mapping(args, 'inputs', *subdirs)

	def _set_output_mapping(self, args, *subdirs):
		# type: (Dict[str, Any], *str) -> None
		"""Set an output path mapping for the given `subdirs` of `sim_outdir`."""
		self._set_path_mapping(args, 'outputs', *subdirs)

	def add_python_task(self, python_args, *upstream_tasks, **kwargs):
		# type: (List[str], *Task, **Any) -> Task
		"""Add a Python task to the workflow and return it."""
		kwargs['image'] = self.image

		# Sisyphus needs `-u` unbuffered output.
		kwargs['commands'] = [{'command': ['python', '-u'] + python_args}]

		return self.add_task(Task(*upstream_tasks, **kwargs))

	def add_parca_task(self, *upstream_tasks, **kwargs):
		# type: (*Task, **Any) -> Task
		kwargs.setdefault('name', 'parca')

		self._set_output_mapping(kwargs, 'kb')

		tokens = ['runscripts/manual/runParca.py', self.sim_outdir]
		_add_options(tokens, kwargs, 'ribosome_fitting', 'rnapoly_fitting')

		return self.add_python_task(tokens, *upstream_tasks, **kwargs)

	def add_variants_task(self, *upstream_tasks, **kwargs):
		# type: (*Task, **Any) -> Task
		kwargs.setdefault('name', 'makeVariants')

		# Set up the output directory mappings. `upstream_tasks` will add the
		# input mappings.
		# TODO(jerry): Compute all the requested variant paths.
		# TODO(jerry): Make the make-variant code output to just one subdir?
		self._set_output_mapping(kwargs, 'wildtype_000000', 'kb')
		self._set_output_mapping(kwargs, 'wildtype_000000', 'metadata')

		tokens = ['runscripts/manual/makeVariants.py', self.sim_outdir]
		_add_options(tokens, kwargs, 'variant')

		return self.add_python_task(tokens, *upstream_tasks, **kwargs)

	def add_sim_task(self, *upstream_tasks, **kwargs):
		# type: (*Task, **Any) -> Task
		kwargs.setdefault('name', 'runSim')

		# Set up the output directory mappings.
		# TODO(jerry): Compute all the requested sim-dir paths.
		#  f'{variant_type}_{i:06d}/{j:06d}/generation_{k:06d}/{l:06d}/simOut'
		#  for variant index i, seed (cohort) j, gen k, daughter l.
		#  Refactor that out from fw_queue, runSim, etc.
		self._set_output_mapping(kwargs, 'wildtype_000000', '000000',
			'generation_000000', '000000', 'simOut')

		tokens = [
			'runscripts/manual/runSim.py', self.sim_outdir,
			'--require_variants', '1']
		_add_options(tokens, kwargs,
			'variant', 'generations', 'total_gens', 'seed', 'init_sims',
			'timeline', 'length_sec', 'timestep_safety_frac', 'timestep_max',
			'timestep_update_freq', 'mass_distribution', 'growth_rate_noise',
			'd_period_division', 'translation_supply', 'trna_charging')

		return self.add_python_task(tokens, *upstream_tasks, **kwargs)


def wc_ecoli_workflow(**kwargs):
	"""Construct a workflow DAG for wcEcoli."""
	wf = WCM_Workflow()

	t_parca = wf.add_parca_task(**kwargs)
	t_variants = wf.add_variants_task(t_parca, **kwargs)

	# TODO(jerry): Loop over variants, seeds, & gens, adding sim & sim daughter
	#  tasks. Change args in kwargs so runSim will only run one iteration and
	#  get the seed, etc. Setup each sim's upstream tasks and simOut path.
	t_sim = wf.add_sim_task(t_variants, **kwargs)

	# TODO(jerry): Add the analysis tasks...

	return wf

def enqueue_workflow(**kwargs):
	"""Construct and enqueue a workflow."""
	wf = wc_ecoli_workflow(**kwargs)
	wf.enqueue()


# TODO(jerry): An argparse CLI.
