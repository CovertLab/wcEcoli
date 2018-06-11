"""
Common code for scripts that manually run analysis plots.

The command line interpreter has one positional parameter sim_dir to identify
the simulation's output directory, defaulting to the latest timestamped (or else
alphabetically first) subdirectory of wcEcoli/out/.

Run with '-h' for command line help.
"""

from __future__ import absolute_import
from __future__ import division

import cPickle
import errno
import os

from runscripts.manual import scriptBase


class AnalysisBase(scriptBase.ScriptBase):
	"""Abstract base class for scripts that manually run analysis plots in an
	existing sim_dir: Defines sim_dir as the first CLI positional parameter;
	uses it to find the sim_path and set the input_validation_data path and the
	metadata_path; loads the metadata; and provides the `find_variant_dir()`
	utility.

	run() is still an abstract method.
	"""

	def description(self):
		"""Describe the command line program."""
		return '{} plots'.format(type(self).__name__)

	def define_parameters(self, parser):
		"""Define command line parameters. Subclasses should call super to
		define the positional parameter "sim_dir", then call
		`parser.add_argument()` as needed. See the super docstring for examples.
		"""
		super(AnalysisBase, self).define_parameters(parser)

		parser.add_argument('sim_dir', nargs='?',
			help='The simulation "out/" subdirectory to read from (optionally'
				+ ' starting with "out/"), or an absolute directory name, or'
				+ ' default to the the most interesting subdirectory of "out/".'
				+ ' It must contain a "metadata/metadata.cPickle" file.')

	def parse_args(self):
		"""Parse the command line args into an `argparse.Namespace`, attach
		derived args for analysis scripts, and return the Namespace of args.

		This attaches a `sim_path` argument derived from the optional "sim_dir"
		argument, the `input_validation_data` path, the
		`metadata_path` path "<sim_path>/metadata/metadata.cPickle", and the
		`metadata` dict loaded from `metadata_path`.

		Overrides should first call super().
		"""
		args = super(AnalysisBase, self).parse_args()

		args.sim_path = scriptBase.find_sim_path(args.sim_dir)

		args.input_validation_data = os.path.join(
			args.sim_path, 'kb', 'validationData.cPickle')

		args.metadata_path = os.path.join(
			args.sim_path, 'metadata', 'metadata.cPickle')
		with open(args.metadata_path) as f:
			args.metadata = cPickle.load(f)

		return args

	def find_variant_dir(self, sim_path):
		"""Find a simulation variant dir in the given `sim_path` or raise an
		IOError.
		"""
		for subdir in os.listdir(sim_path):
			if '_' in subdir and os.path.isdir(os.path.join(sim_path, subdir)):
				return subdir

		raise IOError(errno.ENOENT, 'No simulation variant directory found')


class TestAnalysis(AnalysisBase):
	"""To test out the command line parser."""
	def run(self, analysis_args):
		print "[TEST] Analysis args:", analysis_args


if __name__ == '__main__':
	analysis = TestAnalysis()
	analysis.cli()
