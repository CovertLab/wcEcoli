import cPickle
import time
import os

from fireworks import FireTaskBase, explicit_serialize
from models.ecoli.sim.variants import nameToFunctionMapping

@explicit_serialize
class VariantSimDataTask(FireTaskBase):

	_fw_name = "VariantSimDataTask"
	required_params = [
		"variant_function", "variant_index",
		"input_sim_data", "output_sim_data",
		"variant_metadata_directory",
		]

	def run_task(self, fw_spec):
		return