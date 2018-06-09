"""
Common code for scripts that manually run analysis plots.

The command line interpreter takes one positional argument to identify the
simulation's output directory; and it defaults to the latest timestamped
(or else the alphabetically first) subdirectory of wcEcoli/out/.

Run with '-h' for command line help.
"""

from __future__ import absolute_import
from __future__ import division

import abc
import argparse
import cPickle
import datetime
import os
import time
import wholecell


# The wcEcoli project root path.
ROOT_PATH = os.path.dirname(os.path.dirname(os.path.abspath(wholecell.__file__)))


def default_wcecoli_out_subdir_path():
	"""Return an absolute path to the most interesting subdirectory of
	wcEcoli/out: the one that starts with the latest timestamp or else the
	alphabetically first subdirectory.
	"""
	out_dir = os.path.join(ROOT_PATH, 'out')
	fallback = None

	for directory in sorted(os.listdir(out_dir), reverse=True):
		path = os.path.join(out_dir, directory)
		if os.path.isdir(path):
			if directory[0].isdigit():
				return path
			fallback = path

	if fallback:
		return fallback

	raise IOError("{} has no subdirectories".format(out_dir))

def pick_sim_path(directory=None):
	"""Pick a simulation path, looking for the given directory name as an
	absolute path, or as a subdirectory of wcEcoli/out/, or as a subdirectory
	name that starts with out/, or (if None) default to
	default_wcecoli_out_subdir_path().
	"""
	if directory is None:
		input_dir = default_wcecoli_out_subdir_path()
	elif os.path.isabs(directory):
		input_dir = directory
	elif directory.startswith('out/'):
		input_dir = os.path.join(ROOT_PATH, directory)
	else:
		input_dir = os.path.join(ROOT_PATH, 'out', directory)

	if not os.path.isdir(input_dir):
		raise ValueError("{} is not a simulation path".format(input_dir))
	return input_dir


class AnalysisBase(object):
	"""Abstract base class for scripts that manually run analysis plots."""
	__metaclass__ = abc.ABCMeta

	def description(self):
		"""Describe the command line program."""
		return '{} plots'.format(type(self).__name__)

	def define_parameters(self, parser):
		"""Define any subclass-specific command line parameters to follow the
		standard positional parameter "sim_dir". Override this method to call
		parser.add_argument() as needed.

		Examples include positional arguments
			`parser.add_argument('variant', nargs='?',
			help='simulation variant')`
		options
			`parser.add_argument('--seed', default='000000',
			help='simulation seed')`.
		and flags
			`parser.add_argument('-v', '--verbose', action='store_true',
			help='set verbose logging')`.
		"""
		pass

	def parse_args(self):
		"""Parse the command line args: Construct an ArgumentParser, call
		`define_parameters()` to add any subclass-specific parameters, use it
		to parse the command line into an `argparse.Namespace`, attach a
		"sim_path" argument derived from the optional "sim_dir" argument, and
		return the Namespace of args.

		(A `Namespace` is an object with attributes and some methods like
		`__repr__()` and `__eq__()`. Call `vars(args)` to turn it into a dict.)
		"""
		parser = argparse.ArgumentParser(
			description='Run {}.'.format(self.description()))
		parser.add_argument('sim_dir', nargs='?',
			help='The simulation "out/" subdirectory to read from (optionally'
				+ ' starting with "out/"), or an absolute directory name, or'
				+ ' default to the the most interesting subdirectory of "out/".'
				+ ' It must contain a "metadata/metadata.cPickle" file.')

		self.define_parameters(parser)

		args = parser.parse_args()
		setattr(args, 'sim_path', pick_sim_path(args.sim_dir))
		return args

	def find_variant_dir(self, sim_path, default=''):
		"""Find a simulation variant dir in the given `sim_path`; return the
		given `default` if none are found."""
		for subdir in os.listdir(sim_path):
			if os.path.isdir(subdir) and '_' in subdir:
				return subdir

		return default

	def add_args(self, args):
		"""Add attributes to the Namespace of parsed args. This is an
		opportunity to translate defaulted args into useful values and compute
		derived args.

		This base method adds the `input_validation_data` path, the base
		`output_plots_directory` path, the `metadata_path` path
		"<sim_path>/metadata/metadata.cPickle" (not simData_Modified.cPickle),
		and `metadata` loaded from `metadata_path`.
		Overrides should first call super().

		SIDE EFFECTS: Modifies `args`.
		"""
		args.input_validation_data = os.path.join(
			args.sim_path, 'kb', 'validationData.cPickle')

		args.output_plots_directory = os.path.join(args.sim_path, 'plotOut')

		args.metadata_path = os.path.join(
			args.sim_path, 'metadata', 'metadata.cPickle')
		with open(args.metadata_path) as f:
			args.metadata = cPickle.load(f)

	@abc.abstractmethod
	def run(self, args):
		"""Run the analysis with the given arguments."""
		raise NotImplementedError("AnalysisBase subclass must implement run()")

	def cli(self):
		"""Command-Line Interpreter: parse_args(), add_args(), then run()."""
		args = self.parse_args()
		self.add_args(args)

		print('{}: {} from {}'.format(time.ctime(), self.description(), args.sim_path))
		start_sec = time.clock()
		self.run(args)
		end_sec = time.clock()
		elapsed = datetime.timedelta(seconds = (end_sec - start_sec))
		print "Analysis run in {}h {}m {}s total".format(*str(elapsed).split(':'))


class TestAnalysis(AnalysisBase):
	"""To test out the command line parser."""
	def run(self, analysis_args):
		print "[DEBUG] Analysis args:", analysis_args


if __name__ == '__main__':
	analysis = TestAnalysis()
	analysis.cli()
