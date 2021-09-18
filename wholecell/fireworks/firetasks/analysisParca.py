"""
Run the plots for a given list of `models.ecoli.analysis.parca` analyses, or
by default the ACTIVE plots listed in that package's `__init__.py`.

If the `DEBUG_GC` environment variable is true, enable memory leak detection.
"""

from __future__ import absolute_import, division, print_function

from fireworks import explicit_serialize

from wholecell.fireworks.firetasks.analysisBase import AnalysisBase
import models.ecoli.analysis.parca


@explicit_serialize
class AnalysisParcaTask(AnalysisBase):
	"""Run the Parca analysis plots on both the primary (kb/) dir and (if
	requested) the secondary (kb-poly/) dir.
	"""

	_fw_name = "AnalysisParcaTask"
	required_params = [
		"input_directory",
		"input_sim_data",
		"input_validation_data",
		"output_plots_directory",
		"output_filename_prefix",
		"metadata",

		"plot",
		"cpus",
		"compile",
		]
	optional_params = [
		"input_directory2",  # if specified, the other secondary params are required
		"input_sim_data2",
		"input_validation_data2",
		]
	MODULE_PATH = 'models.ecoli.analysis.parca'
	TAGS = models.ecoli.analysis.parca.TAGS

	def run_task(self, fw_spec):
		self["metadata"] = dict(self["metadata"], analysis_type="parca")

		# Run with the primary params.
		self["_params"] = ''
		super().run_task(fw_spec)

		# If requested, run with the secondary params.
		if self.get("input_directory2"):
			self["_params"] = '2'
			super().run_task(fw_spec)

	def plotter_args(self, module_filename):
		params = self.get("_params", '')  # select primary or secondary params
		prefix = 'polycistronic_' if params else ''

		return (
			self["input_directory" + params],
			self["output_plots_directory"],
			self["output_filename_prefix"] + prefix + module_filename[:-3],
			self["input_sim_data" + params],
			self["input_validation_data" + params],
			self["metadata"],
			)

	def description(self):
		params = self.get("_params", '')
		return f'{type(self).__name__} on {self["input_directory" + params]}'
