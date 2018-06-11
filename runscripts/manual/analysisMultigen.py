"""
Run all multigen analysis plots for a given sim.

Run with '-h' for command line help.

Set the environment variable WC_ANALYZE_FAST to run multiple analysis scripts
in parallel processes.
"""

from __future__ import absolute_import
from __future__ import division

import os

from wholecell.utils import filepath

from runscripts.manual.analysisBase import AnalysisBase
from wholecell.fireworks.firetasks.analysisMultiGen import AnalysisMultiGenTask


DIR = "000000"


class AnalysisMultigen(AnalysisBase):
	"""Runs all multigen analysis plots for a given sim."""

	def define_parameters(self, parser):
		super(AnalysisMultigen, self).define_parameters(parser)
		parser.add_argument('--variant',
			help='simulation variant, e.g. "geneKnockdown_000030"')

	def parse_args(self):
		args = super(AnalysisMultigen, self).parse_args()

		if args.variant is None:  # defaulted
			args.variant = self.find_variant_dir(args.sim_path)

		metadata = args.metadata
		metadata['analysis_type'] = 'cohort'
		metadata['variant_function'] = args.variant
		metadata['variant_index'] = None

	def run(self, args):
		sim_path = args.sim_path
		variant = args.variant

		input_path = os.path.join(sim_path, variant, DIR)
		sim_data_modified = os.path.join(
			sim_path, variant, 'kb', 'simData_Modified.cPickle')
		# TODO(jerry): Load simData_Modified into metadata?
		output_dir = filepath.makedirs(sim_path, variant, DIR, "plotOut")

		task = AnalysisMultiGenTask(
			input_seed_directory=input_path,
			input_sim_data=sim_data_modified,
			input_validation_data=args.input_validation_data,
			output_plots_directory=output_dir,
			metadata=args.metadata)
		task.run_task({})


if __name__ == '__main__':
	analysis = AnalysisMultigen()
	analysis.cli()
