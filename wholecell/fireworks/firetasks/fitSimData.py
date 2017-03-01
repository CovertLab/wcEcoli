import cPickle
import time

from fireworks import FireTaskBase, explicit_serialize
from reconstruction.ecoli.fit_sim_data_1 import fitSimData_1
from reconstruction.ecoli.fit_sim_data_2 import fitSimData_2

@explicit_serialize
class FitSimDataTask(FireTaskBase):

	_fw_name = "FitSimDataTask"
	required_params = ["fit_level", "input_data", "output_data"]
	optional_params = ["sim_out_dir"]

	def run_task(self, fw_spec):
		return