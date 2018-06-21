"""
Common code for analysis plots. The abstract base class AnalysisPlot defines a
plot() method for scripts to call.

TODO: Fill in for defaulted command line args. Use ScriptBase methods for
finding variant directories, etc.

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

TODO: Other shared code to simplify the subclasses, e.g. make plotOutDir,
check that `os.path.isdir(simOutDir)`, instantiate an AnalysisPaths (except for
SingleAnalysisPlot subclasses), etc.
"""

from __future__ import absolute_import
from __future__ import division

import abc

from wholecell.utils import scriptBase


class AnalysisPlot(scriptBase.ScriptBase):
	"""Abstract Base Class for analysis plots.

	Each analysis class must override do_plot().

	Call cli() to run the command line interface for one analysis class.

	Call main() to run an analysis plot for a Firetask.

	The abstract subclasses CohortAnalysisPlot et al will have methods to run
	all their current analyses in a controlled order.
	"""

	def description(self):
		"""Describe the command line program."""
		return '{cls.__module__}.{cls.__name__}'.format(cls=type(self))

	def define_parameters(self, parser):
		"""Define command line parameters to run a single analysis plot as a
		standalone script.
		"""
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
	def do_plot(self, inputDir, plotOutDir, plotOutFileName, simDataFile,
			validationDataFile, metadata):
		"""Inner method that each analysis class must override."""
		raise NotImplementedError("AnalysisPlot subclass must implement do_plot()")

	def plot(self, inputDir, plotOutDir, plotOutFileName, simDataFile,
			validationDataFile, metadata):
		"""Public method to set up, make a plot, and cleanup."""
		# TODO: Setup.

		self.do_plot(inputDir, plotOutDir, plotOutFileName, simDataFile,
			validationDataFile, metadata)

		# TODO: Cleanup.


	def run(self, args):
		"""Run an analysis plot with the given CLI arguments. Called by cli()."""
		self.plot(
			args.sim_path,
			args.plotOutDir,
			args.plotOutFileName,
			args.simDataFile,
			getattr(args, 'validationDataFile', None),
			getattr(args, 'metadata', None),
		)

	def main(self, inputDir, plotOutDir, plotOutFileName, simDataFile,
			validationDataFile=None, metadata=None):
		"""Run an analysis plot for a Firetask."""
		self.plot(inputDir, plotOutDir, plotOutFileName, simDataFile,
			validationDataFile, metadata)
