import cPickle
import time

from fireworks import FireTaskBase, explicit_serialize
from validation.ecoli.validation_data import ValidationDataEcoli

@explicit_serialize
class InitValidationDataTask(FireTaskBase):

	_fw_name = "InitValidationDataTask"
	required_params = ["validation_data_input", "knowledge_base_raw", "output_data"]

	def run_task(self, fw_spec):
		return