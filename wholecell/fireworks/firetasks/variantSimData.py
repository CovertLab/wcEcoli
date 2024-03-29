import os
import pickle
import sys

from fireworks import FiretaskBase, explicit_serialize
from models.ecoli.sim.variants import apply_variant
import wholecell.utils.filepath as fp


@explicit_serialize
class VariantSimDataTask(FiretaskBase):

	_fw_name = "VariantSimDataTask"
	required_params = [
		"variant_function",
		"variant_index",
		"input_sim_data",
		"output_sim_data",  # the output variant sim_data file (dir gets created)
		"variant_metadata_directory",  # the output variant metadata directory
		]

	def run_task(self, fw_spec):
		fp.makedirs(os.path.dirname(self["output_sim_data"]))
		fp.makedirs(self["variant_metadata_directory"])

		info, sim_data = apply_variant.apply_variant(
			self["input_sim_data"],
			self["variant_function"],
			self["variant_index"])

		sys.setrecursionlimit(4000)

		with open(self["output_sim_data"], "wb") as f:
			pickle.dump(sim_data, f, protocol = pickle.HIGHEST_PROTOCOL)

		with open(os.path.join(self["variant_metadata_directory"], "short_name"), "w") as h:
			h.write("%s\n" % info["shortName"])

		with open(os.path.join(self["variant_metadata_directory"], "description"), "w") as h:
			h.write("%s\n" % info["desc"])
