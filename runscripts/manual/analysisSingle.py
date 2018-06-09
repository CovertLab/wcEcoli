"""
Runs all single analysis plots for a given sim w/optional variant.

Run with '-h' for command line help.
"""

from __future__ import absolute_import
from __future__ import division

import os

from runscripts.manual.analysisBase import AnalysisBase
from wholecell.fireworks.firetasks.analysisSingle import AnalysisSingleTask


SEED = '000000'
GEN = 'generation_000000'
DAUGHTER = '000000'
DIRS = os.path.join(SEED, GEN, DAUGHTER)


class AnalysisSingle(AnalysisBase):
	"""Runs all single analysis plots for a given sim w/optional variant."""

	def define_parameters(self, parser):
		super(AnalysisSingle, self).define_parameters(parser)
		parser.add_argument('--variant',
			help='simulation variant, e.g. "condition_000000"')

	def add_args(self, args):
		super(AnalysisSingle, self).add_args(args)

		if args.variant is None:  # defaulted
			args.variant = self.find_variant_dir(args.sim_path)

		metadata = args.metadata
		metadata['analysis_type'] = 'single'
		metadata['variant_function'] = args.variant
		metadata['variant_index'] = None
		metadata['seed'] = SEED
		metadata['gen'] = GEN

	def run(self, args):
		sim_path = args.sim_path
		variant = args.variant

		results_dir = os.path.join(sim_path, variant, DIRS, 'simOut')
		sim_data_modified = os.path.join(
			sim_path, variant, 'kb', 'simData_Modified.cPickle')
		# TODO(jerry): Load simData_Modified into metadata?
		output_dir = os.path.join(sim_path, variant, DIRS, 'plotOut')

		task = AnalysisSingleTask(
			input_results_directory=results_dir,
			input_sim_data=sim_data_modified,
			input_validation_data=args.input_validation_data,
			output_plots_directory=output_dir,
			metadata=args.metadata)
		task.run_task({})


if __name__ == '__main__':
	analysis = AnalysisSingle()
	analysis.cli()
