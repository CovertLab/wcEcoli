#!/usr/bin/env python

"""Run a custom command line in Google Cloud using Fireworks, with access to
the GCS files from a WCM workflow run.

Prerequisite: Run `cloud/build-wcm.sh` to build a Docker Image containing the
command line code to invoke in the workers. The code doesn't have to be checked
in to git, just in the directory that build-wcm.sh bundles up."""

import os
import posixpath as pp

from borealis.util import gcp

import wholecell.utils.filepath as fp
from runscripts.cloud.util.workflow import Task
from runscripts.cloud.util.workflow_cli import WorkflowCLI


# SEEDS = sorted(set(range(100)) - {9, 10, 28, 30, 49, 57, 58, 74, 80, 94} | {128} | set(range(300, 310)))


class CustomWorkflow(WorkflowCLI):
	"""A workflow to run a custom command line in Google Cloud."""

	DEFAULT_TIMEOUT = 60 * 60  # add_task() default timeout, in seconds
	DOCKER_IMAGE_TEMPLATE = 'gcr.io/{}/{}-wcm-code'

	def __init__(self):
		super().__init__()
		self.internal_prefix = pp.join(pp.sep, 'wcEcoli', 'out')
		self.DOCKER_IMAGE = ''

	def add_analysis_task(self, seed, num_gens):
		# type: (int, int) -> Task
		"""Add an analysis task to the workflow."""
		def in_sim_dir(*path_elements):
			return pp.join(sim_dir, *path_elements)

		base = self.internal('wildtype_000000')  # 'counts/wildtype_000000'???
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
					 pp.join('prototypes',
							 'subgenerational_analysis',
							 'pull_complex_counts_cloud.py'),
					 pp.join(base, seed_key)],
			inputs=inputs,
			outputs=[pp.join(base, 'count_out', seed_key, 'complex')])

	def build(self, args):
		"""Build the workflow."""
		for seed in range(args.seeds):
			self.add_analysis_task(seed, args.generations)

	def define_parameters(self, parser):
		self.define_option(parser, 'generations', int, 1, flag='g',
			help='The number of generations to analyze.')
		self.define_option(parser, 'seeds', int, 1, flag='s',
			help='The number of seeds to analyze.')  # TODO(jerry): Take a seed list?

		self.define_option(parser, 'description', str, '',
			help='A simulation description; part of the storage folder name.')
		self.define_option(parser, 'id', str, default=None,
			help='Workflow ID or owner ID such as a user name or a CI build'
				 ' name to combine with the timestamp to form the unique'
				 ' workflow name. Default = $WF_ID environment variable or'
				 ' else the $USER environment variable.')
		self.define_option(parser, 'timestamp', str, fp.timestamp(),
			help='Timestamp for this workflow. It gets combined with the'
				 ' Workflow ID to form the workflow name. Set this if you want'
				 ' to upload new steps for an existing workflow. Default ='
				 ' the current local date-time.')
		self.define_parameter_bool(parser, 'verbose', True,
			help='Verbose workflow builder logging.')

		super().define_parameters(parser)

	def run(self, args):
		owner_id = args.id or os.environ.get('WF_ID', os.environ['USER'])
		setattr(args, 'basename', 'WCM')
		setattr(args, 'owner_id', owner_id)
		self.DOCKER_IMAGE = self.DOCKER_IMAGE_TEMPLATE.format(
			gcp.project(), owner_id)
		super().run(args)


if __name__ == '__main__':
	CustomWorkflow().cli()
