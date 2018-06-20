"""
Common code for analysis plots. The abstract base class AnalysisPlot defines a
plot() method for scripts to call.

TODO: Fill in for defaulted command line args.

TODO: Set the fallback for the matplotlib back end so we needn't hack its
installation. Set the CWD so matplotlib finds the matplotlibrc file even when
FireWorks doesn't launch it in the right directory.

TODO: In the abstract subclasses CohortAnalysisPlot et al, add methods to
instantiate and run a list of subclasses in a controlled order, with unneeded
ones commented out, to simplify the Firetasks. That also makes it easy to
adjust the order when you're working on particular ones.

TODO: Setup/reset matplotlib before each script and cleanup afterwards.

TODO: Enable future warnings, esp. for matplotlib.

TODO: Memory leak detection.

TODO: Move the run_plot() args to instance variables?

TODO: Other shared code to simplify the subclasses, e.g. make the output dir.
"""

from __future__ import absolute_import
from __future__ import division

import abc

from wholecell.utils import scriptBase


class AnalysisPlot(scriptBase.ScriptBase):
	"""Abstract Base Class for analysis plots.
	Call the cli() method to run the command line interface.
	"""

	def description(self):
		"""Describe the command line program."""
		return '{cls.__module__}.{cls.__name__}'.format(cls=type(self))

	def define_parameters(self, parser):
		super(AnalysisPlot, self).define_parameters(parser)
		self.define_parameter_sim_dir(parser)

		parser.add_argument("--plotOutDir",
			help="Directory for plot output (will get created if necessary).")
		parser.add_argument("--plotOutFileName", help="Plot output filename.")
		parser.add_argument("--simDataFile", help="KB file name.")

	def parse_args(self):
		"""Parse the command line args into an `argparse.Namespace`, including
		the `sim_dir` and `sim_path` args. Overrides should first call super().
		"""
		args = super(AnalysisPlot, self).parse_args()

		# TODO: Fill in for defaulted args. Default args.simDataFile to
		# constants.SERIALIZED_KB_DIR/constants.SERIALIZED_KB_MOST_FIT_FILENAME?

		return args

	@abc.abstractmethod
	def run_plot(self, inputDir, plotOutDir, plotOutFileName, simDataFile,
			validationDataFile, metadata):
		"""Inner method to write out a plot with the given arguments."""
		raise NotImplementedError("AnalysisPlot subclass must implement plot()")

	def plot(self, inputDir, plotOutDir, plotOutFileName, simDataFile,
			validationDataFile=None, metadata=None):
		"""Public method to set up, write out a plot, and cleanup."""
		# TODO: Setup.

		self.run_plot(inputDir, plotOutDir, plotOutFileName, simDataFile,
			validationDataFile, metadata)

		# TODO: Cleanup.

	def run(self, args):
		"""Public method to write out a plot with the given CLI arguments."""
		self.plot(
			args.sim_path,
			args.plotOutDir,
			args.plotOutFileName,
			args.simDataFile,
			getattr(args, 'validationDataFile', None),
			getattr(args, 'metadata', None),
		)

	@classmethod
	def main(cls, inputDir, plotOutDir, plotOutFileName, simDataFile,
			validationDataFile=None, metadata=None):
		"""Instantiate this (sub)class and run a plot."""
		instance = cls()
		instance.plot(inputDir, plotOutDir, plotOutFileName, simDataFile,
			validationDataFile, metadata)
