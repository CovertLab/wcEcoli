"""
Runs all parca analysis plots for a given sim.

Run with '-h' for command line help.
"""

from __future__ import absolute_import, division, print_function

import os

from runscripts.manual.analysisBase import AnalysisBase
from wholecell.fireworks.firetasks.analysisParca import AnalysisParcaTask
from wholecell.fireworks.firetasks.parca import ParcaTask
from wholecell.utils import constants


class AnalysisParca(AnalysisBase):
	"""Runs some or all the ACTIVE parca analysis plots for a given sim."""

	def run(self, args):
		# TODO: where to put output?
		# Does kb/plotOut make sense? - Is this compatible with cloud runs?
		output_dir = os.path.join(args.sim_path, constants.KB_PLOT_OUTPUT_DIR)
		input_sim_data = os.path.join(args.sim_path,
			ParcaTask.OUTPUT_SUBDIR, constants.SERIALIZED_SIM_DATA_FILENAME)

		task = AnalysisParcaTask(
			input_directory=args.sim_path,
			input_sim_data=input_sim_data,
			input_validation_data=args.input_validation_data,
			output_plots_directory=output_dir,
			metadata=args.metadata,
			output_filename_prefix=args.output_prefix,
			**self.select_analysis_keys(args)
			)
		task.run_task({})


if __name__ == '__main__':
	analysis = AnalysisParca()
	analysis.cli()
