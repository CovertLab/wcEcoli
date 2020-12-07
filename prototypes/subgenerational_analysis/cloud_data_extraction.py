#!/usr/bin/env python

"""Run a custom command line in Google Cloud using Fireworks, with access to
the GCS files from a WCM workflow run.

Prerequisite: Run `cloud/build-wcm.sh` to build a Docker Image containing the
command line code to invoke in the workers. The code doesn't have to be checked
in to git, just in the directories that build-wcm.sh bundles up.

The owner_id, timestamp arg, and description CLI args must match those of a
previous WCM workflow run so this can fetch input files from that workflow's
outputs and store output files to the same directory.
"""

import os
import posixpath as pp

from borealis.util import gcp

from runscripts.cloud.util.workflow import Task
from runscripts.cloud.util.workflow_cli import WorkflowCLI


# SEEDS = sorted(set(range(100)) - {9, 10, 28, 30, 49, 57, 58, 74, 80, 94} | {128} | set(range(300, 310)))


class CustomWorkflow(WorkflowCLI):
	"""A workflow to run a custom command line in Google Cloud."""

	WORKFLOW_BASENAME = 'ComplexCounts'
	DEFAULT_TIMEOUT = 10 * 60  # add_task() default timeout, in seconds

	def __init__(self):
		super().__init__(internal_prefix=pp.join(pp.sep, 'wcEcoli', 'out', 'wf'))
		self.DOCKER_IMAGE = ''  # set it after parsing CLI parameters

	def add_analysis_task(self, seed, num_gens):
		# type: (int, int) -> Task
		"""Add an analysis task to the workflow."""
		def in_sim_dir(*path_elements):
			return pp.join(sim_dir, *path_elements)

		base = self.internal('wildtype_000000')
		seed_key = format(seed, '06')  # in 6-digit format
		inputs = []

		for generation in range(num_gens):
			generation_key = 'generation_{:06d}'.format(generation)
			sim_dir = pp.join(
				base, seed_key, generation_key, '000000', 'simOut')
			inputs.append(in_sim_dir('BulkMolecules', 'counts'))
			inputs.append(in_sim_dir('BulkMolecules', 'attributes.json'))
			inputs.append(in_sim_dir('Main', 'time'))
			inputs.append(in_sim_dir('Main', 'attributes.json'))

		return self.add_task(
			name='pull_complex_counts_cloud{}'.format(seed),  # unique task name
			command=['python',
					 '-u',
					 # 'prototypes/subgenerational_analysis/pull_complex_counts_cloud.py',
					 'runscripts/cloud/util/multigen_analysis_example.py',
					 pp.join(base, seed_key)],
			inputs=inputs,
			outputs=[pp.join(base, 'count_out', seed_key, 'complex', '')])

	def build(self, args):
		"""Build the workflow."""
		seeds = [int(seed) for seed in args.seed_list.split()]
		for seed in seeds:
			self.add_analysis_task(seed, args.generations)

	def define_parameters(self, parser):
		self.define_option(parser, 'generations', int, 1, flag='g',
			help='The number of generations to analyze.')
		self.define_option(parser, 'seed_list', str, '0', flag='s',
			help='The list of cell sim seed numbers to analyze, e.g. "0 1 2".')

		self.define_wf_name_parameters(parser)

		super().define_parameters(parser)

	def run(self, args):
		owner_id = args.id or os.environ.get('WF_ID', os.environ['USER'])
		setattr(args, 'owner_id', owner_id)

		# Fetch from and store to a WCM workflow's storage directory, assuming
		# the owner_id, timestamp arg, and description arg match.
		setattr(args, 'storage_basename', 'WCM')

		self.DOCKER_IMAGE = f'gcr.io/{gcp.project()}/{owner_id}-wcm-code'
		super().run(args)


if __name__ == '__main__':
	CustomWorkflow().cli()
