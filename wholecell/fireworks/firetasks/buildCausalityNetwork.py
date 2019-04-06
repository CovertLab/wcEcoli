"""
Run the analysis scripts that generate input files to the Causality Network
tool.
"""
from __future__ import absolute_import, division, print_function

import importlib
import os
import time
import traceback

from fireworks import FireTaskBase, explicit_serialize
from models.ecoli.analysis.causality_network.build_network import BuildNetwork
from models.ecoli.analysis.causality_network.network_components import NODELIST_JSON


@explicit_serialize
class BuildCausalityNetworkTask(FireTaskBase):

	_fw_name = "BuildCausalNetworkTask"
	required_params = [
		"input_results_directory",
		"input_sim_data",
		"output_network_directory",
		"output_dynamics_directory",
		"metadata",
		]
	optional_params = [
		"check_sanity",
		]

	READER_FILE_PATH = 'models.ecoli.analysis.causality_network.read_dynamics'

	def plotter_args(self):
		self["metadata"] = dict(self["metadata"], analysis_type = "causality_network")

		return (
			self["input_results_directory"],
			self["output_dynamics_directory"],
			None,
			self["input_sim_data"],
			self["node_list_file"],
			self["metadata"],
			)

	def run_task(self, fw_spec):
		startTime = time.time()
		print("\n{}: --- Starting {} ---".format(
			time.ctime(startTime), type(self).__name__))

		self['node_list_file'] = os.path.join(
			self["output_network_directory"], NODELIST_JSON)

		self["check_sanity"] = self.get("check_sanity", False)

		# Build network first if the node list file does not exist
		if not os.path.isfile(self['node_list_file']):
			print("{}: Building causality network".format(time.ctime()))

			causality_network = BuildNetwork(
				self["input_sim_data"], self["output_network_directory"],
				self["check_sanity"])
			causality_network.run()

		# Read dynamics from sim results if the dynamics directory is empty
		if len(os.listdir(self["output_dynamics_directory"])) == 0:
			mod = importlib.import_module(self.READER_FILE_PATH)
			args = self.plotter_args()

			print("{}: Reading simulation results for causality network"
				.format(time.ctime()))
			mod.Plot.main(*args)

		timeTotal = time.time() - startTime

		duration = time.strftime("%H:%M:%S", time.gmtime(timeTotal))
		print("{}: Completed building causality network in {}".format(
			time.ctime(), duration)
			)

		# Root directory of wcEcoli repo
		root_dir = os.path.dirname(os.path.dirname(
			os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

		# Launch app through causality repo
		print("{}: Launching causality network app".format(time.ctime()))
		os.system("%s %s %s" % (
			os.path.join(os.path.dirname(root_dir), "causality", "site", "server.py"),
			os.path.dirname(self["node_list_file"]),
			self["output_dynamics_directory"]))
