"""
Common code for scripts that manually run simulation and analysis operations
outside of Fireworks workflows.

Run with '-h' for command line help.
"""

from __future__ import absolute_import
from __future__ import division

import abc
import argparse
import datetime
import errno
import os
import time

import wholecell


# The wcEcoli project root path.
ROOT_PATH = os.path.dirname(os.path.dirname(os.path.abspath(wholecell.__file__)))


def default_wcecoli_out_subdir_path():
	"""Return an absolute path to the most interesting subdirectory of
	wcEcoli/out/: the one that starts with the latest timestamp or else the
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

	raise IOError(errno.ENOENT, "{} has no subdirectories".format(out_dir))

def find_sim_path(directory=None):
	"""Find a simulation path, looking for the given directory name as an
	absolute path, or as a subdirectory of wcEcoli/out/, or as a subdirectory
	name that starts with out/, or (if None) call
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
		raise IOError(errno.ENOENT, "{} is not a simulation path".format(input_dir))
	return input_dir


class ScriptBase(object):
	"""Abstract base class for scripts. This defines a template where
	`description()` describes the script,
	`define_parameters()` defines its command line parameters,
	`parse_args()` parses the command line args,
	`run()` does the work,
	`cli()` is the driving Command-Line Interpreter.
	"""
	__metaclass__ = abc.ABCMeta

	def description(self):
		"""Describe the command line program."""
		return type(self).__name__

	def define_parameters(self, parser):
		"""Define command line parameters.

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
		`define_parameters()` to define parameters including subclass-specific
		parameters, use it to parse the command line into an
		`argparse.Namespace`, and return that.

		Overrides should first call super().

		(A `Namespace` is an object with attributes and some methods like
		`__repr__()` and `__eq__()`. Call `vars(args)` to turn it into a dict.)
		"""
		parser = argparse.ArgumentParser(
			description='Run {}.'.format(self.description()))

		self.define_parameters(parser)

		args = parser.parse_args()
		return args

	@abc.abstractmethod
	def run(self, args):
		"""Run the operation with the given arguments."""
		raise NotImplementedError("ScriptBase subclass must implement run()")

	def cli(self):
		"""Command Line Interpreter: parse_args() then run(). This also prints
		a starting message (including args.sim_path if defined) and an ending
		message (including the elapsed run time).
		"""
		args = self.parse_args()

		location = getattr(args, 'sim_path', '')
		if location:
			location = ' at ' + location

		print('{}: {}{}'.format(time.ctime(), self.description(), location))

		start_sec = time.clock()
		self.run(args)
		end_sec = time.clock()
		elapsed = datetime.timedelta(seconds = (end_sec - start_sec))

		print "Run in {}h {}m {}s total".format(*str(elapsed).split(':'))


class TestScript(ScriptBase):
	"""To test out the command line parser."""

	def define_parameters(self, parser):
		super(TestScript, self).define_parameters(parser)
		parser.add_argument('--seed', default='000001', help='simulation seed')

	def run(self, args):
		print "[TEST] Run args:", args


if __name__ == '__main__':
	script = TestScript()
	script.cli()
