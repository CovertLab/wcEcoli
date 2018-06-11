"""
Run the Fitter. The output goes into the named subdirectory of wcEcoli/out/,
defaulting to "manual".

TODO: Call the firetasks to do the work in the same way as in FireWorks.

Run with '-h' for command line help.
Set PYTHONPATH when running this.
"""

from __future__ import absolute_import
from __future__ import division

import cPickle
import time
import os

from reconstruction.ecoli.knowledge_base_raw import KnowledgeBaseEcoli
from reconstruction.ecoli.fit_sim_data_1 import fitSimData_1
from wholecell.utils import filepath
from runscripts.manual import scriptBase


class RunFitter(scriptBase.ScriptBase):
	"""Runs the Fitter, aka simulation parameter calculator."""

	def define_parameters(self, parser):
		super(RunFitter, self).define_parameters(parser)

		parser.add_argument('sim_dir', nargs='?', default='manual',
			help='The simulation "out/" subdirectory to write to.')
		parser.add_argument('-d', '--debug', action='store_true',
			help="Enable Fitter debugging. This fits only one"
				 " arbitrarily-chosen transcription factor for a faster debug"
				 " cycle. Don't use it for an actual simulation.")

	def parse_args(self):
		args = super(RunFitter, self).parse_args()
		args.sim_path = filepath.makedirs(
			scriptBase.ROOT_PATH, "out", args.sim_dir)
		return args

	def run(self, args):
		location = filepath.makedirs(args.sim_path, "kb")
		raw_data_file = os.path.join(location, "rawData.cPickle")
		sim_data_file = os.path.join(location, "simData_Fit_1.cPickle")

		print "{}: Loading raw".format(time.ctime())
		raw_data = KnowledgeBaseEcoli()

		print "{}: {}Fitting".format(time.ctime(),
			'DEBUG-' if args.debug else '')
		sim_data = fitSimData_1(raw_data, debug=args.debug)
		print "{}: Done fitting".format(time.ctime())

		with open(raw_data_file, "wb") as f:
			cPickle.dump(raw_data, f, protocol=cPickle.HIGHEST_PROTOCOL)
		with open(sim_data_file, "wb") as f:
			cPickle.dump(sim_data, f, protocol=cPickle.HIGHEST_PROTOCOL)
		print "{}: Wrote parameter files {}, {}".format(time.ctime(),
			raw_data_file, sim_data_file)


if __name__ == '__main__':
	script = RunFitter()
	script.cli()
