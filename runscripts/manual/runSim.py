"""
Run a simulation.
Only runs first gen at the moment.

TODO: Call the firetask to do the work in the same way as the FireWorks
workflow.

Run with '-h' for command line help.
Set PYTHONPATH when running this.
"""

from __future__ import absolute_import
from __future__ import division

import errno
import os

from models.ecoli.sim.simulation import EcoliSimulation
from runscripts.manual import scriptBase



class RunSimulation(scriptBase.ScriptBase):
	"""Drives a simple simulation run."""

	def define_parameters(self, parser):
		super(RunSimulation, self).define_parameters(parser)

		parser.add_argument('sim_dir', nargs='?',
			help='The simulation "out/" subdirectory to read from (optionally'
				+ ' starting with "out/"), or an absolute directory name, or'
				+ ' default to the the most interesting subdirectory of "out/".')

	def parse_args(self):
		args = super(RunSimulation, self).parse_args()
		args.sim_path = scriptBase.find_sim_path(args.sim_dir)
		return args

	def run(self, args):
		sim_data_file = os.path.join(args.sim_path, 'kb', 'simData_Fit_1.cPickle')
		if not os.path.isfile(sim_data_file):
			raise IOError(
				errno.ENOENT,
				'Missing "{}".  Run the Fitter.'.format(sim_data_file))

		options = {
			'simDataLocation': sim_data_file,
			'outputDir': os.path.join(args.sim_path, 'sim'),
			# 'lengthSec': 30,
		}

		print 'Simulation options: {}\n'.format(options)

		sim = EcoliSimulation(**options)
		sim.run()


if __name__ == '__main__':
	script = RunSimulation()
	script.cli()
