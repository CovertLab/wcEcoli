"""
Runs all variant analysis plots for a given sim.

Run with '-h' for command line help.
"""

from __future__ import absolute_import
from __future__ import division

import os

from runscripts.manual.analysisBase import AnalysisBase
from wholecell.fireworks.firetasks.analysisVariant import AnalysisVariantTask


class AnalysisVariant(AnalysisBase):
	"""Runs all variant analysis plots for a given sim."""

	def add_args(self, args):
		super(AnalysisVariant, self).add_args(args)

		metadata = args.metadata
		metadata['analysis_type'] = 'variant'
		metadata['total_variants'] = None

	def run(self, args):
		output_dir = os.path.join(args.sim_path, 'plotOut')

		task = AnalysisVariantTask(
			input_directory=args.sim_path,
			input_validation_data=args.input_validation_data,
			output_plots_directory=output_dir,
			metadata=args.metadata)
		task.run_task({})


if __name__ == '__main__':
	analysis = AnalysisVariant()
	analysis.cli()
