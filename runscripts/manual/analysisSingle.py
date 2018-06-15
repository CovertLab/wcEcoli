"""
Runs all single analysis plots for a given variant of a given sim.

Run with '-h' for command line help.
"""

from __future__ import absolute_import
from __future__ import division

import os

from runscripts.manual.analysisBase import AnalysisBase
from wholecell.fireworks.firetasks.analysisSingle import AnalysisSingleTask
from wholecell.utils import constants
from wholecell.utils import filepath


SEED = '000000'
GEN = 'generation_000000'
DAUGHTER = '000000'
DIRS = os.path.join(SEED, GEN, DAUGHTER)


class AnalysisSingle(AnalysisBase):
	"""Runs all single analysis plots for a given sim w/optional variant."""

	def define_parameters(self, parser):
		super(AnalysisSingle, self).define_parameters(parser)
		self.define_parameter_variant_index(parser)

	def parse_args(self):
		args = super(AnalysisSingle, self).parse_args()

		metadata = args.metadata
		metadata['analysis_type'] = 'single'
		metadata['seed'] = SEED
		metadata['gen'] = GEN

		return args

	def run(self, args):
		sim_path = args.sim_path
		variant_dir_name = args.variant_dir_name

		input_variant_directory = os.path.join(sim_path, variant_dir_name)
		input_dir = os.path.join(input_variant_directory, DIRS, 'simOut')
		sim_data_modified = os.path.join(input_variant_directory, 'kb',
			constants.SERIALIZED_SIM_DATA_MODIFIED)
		output_dir = filepath.makedirs(input_variant_directory, DIRS, 'plotOut')

		task = AnalysisSingleTask(
			input_results_directory=input_dir,
			input_sim_data=sim_data_modified,
			input_validation_data=args.input_validation_data,
			output_plots_directory=output_dir,
			metadata=args.metadata)
		task.run_task({})


if __name__ == '__main__':
	analysis = AnalysisSingle()
	analysis.cli()
