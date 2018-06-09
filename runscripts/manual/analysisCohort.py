"""
Runs all cohort analysis plots for a given sim with data from one location to
another.

Run with '-h' for command line help.
"""

from __future__ import absolute_import
from __future__ import division

import os

from runscripts.manual.analysisBase import AnalysisBase
from wholecell.fireworks.firetasks.analysisCohort import AnalysisCohortTask


class AnalysisCohort(AnalysisBase):
	"""Runs all cohort analysis plots for a given sim."""

	def define_parameters(self, parser):
		super(AnalysisCohort, self).define_parameters(parser)
		parser.add_argument('--variant',
			help='simulation variant, e.g. "condition_000001"')

	def add_args(self, args):
		super(AnalysisCohort, self).add_args(args)

		if args.variant is None:  # defaulted
			args.variant = self.find_variant_dir(args.sim_path)

		metadata = args.metadata
		metadata['analysis_type'] = 'cohort'
		metadata['variant_function'] = args.variant
		metadata['variant_index'] = None

	def run(self, args):
		sim_path = args.sim_path
		variant = args.variant

		input_variant_directory = os.path.join(sim_path, variant)
		sim_data_modified = os.path.join(
			sim_path, variant, 'kb/simData_Modified.cPickle')
		# TODO(jerry): Load simData_Modified into metadata?
		output_dir = os.path.join(sim_path, variant, 'plotOut')

		task = AnalysisCohortTask(
			input_variant_directory=input_variant_directory,
			input_sim_data=sim_data_modified,
			input_validation_data=args.input_validation_data,
			output_plots_directory=output_dir,
			metadata=args.metadata)
		task.run_task({})


if __name__ == '__main__':
	analysis = AnalysisCohort()
	analysis.cli()
