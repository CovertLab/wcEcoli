"""
Runs all parca analysis plots for a given sim.

Run with '-h' for command line help.
"""

from __future__ import absolute_import, division, print_function

import os

from runscripts.manual.analysisBase import AnalysisBase
from wholecell.fireworks.firetasks.analysisParca import AnalysisParcaTask
from wholecell.utils import constants


class AnalysisParca(AnalysisBase):
	"""Run some or all the ACTIVE parca analysis plots for a given sim_path on
	both its primary (kb/) dir and (if present) its secondary (kb-poly/) dir.
	"""

	def run(self, args):
		input_dir = os.path.join(args.sim_path, constants.KB_DIR)
		output_dir = os.path.join(args.sim_path, constants.KB_PLOT_OUTPUT_DIR)
		input_sim_data = os.path.join(
			input_dir, constants.SERIALIZED_SIM_DATA_FILENAME)
		input_validation_data = os.path.join(
			input_dir, constants.SERIALIZED_VALIDATION_DATA)

		params = dict(self.select_analysis_keys(args),
			input_directory=input_dir,
			input_sim_data=input_sim_data,
			input_validation_data=input_validation_data,
			output_plots_directory=output_dir,
			metadata=args.metadata,
			output_filename_prefix=args.output_prefix)

		input_dir2 = os.path.join(args.sim_path, constants.PKB_DIR)
		input_sim_data2 = os.path.join(
			input_dir2, constants.SERIALIZED_SIM_DATA_FILENAME)
		input_validation_data2 = os.path.join(
			input_dir2, constants.SERIALIZED_VALIDATION_DATA)

		if os.path.isfile(input_sim_data2):
			params.update(dict(
				input_directory2=input_dir2,
				input_sim_data2=input_sim_data2,
				input_validation_data2=input_validation_data2))

		task = AnalysisParcaTask(**params)
		task.run_task({})


if __name__ == '__main__':
	analysis = AnalysisParca()
	analysis.cli()
