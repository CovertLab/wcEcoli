import cPickle
import time

from fireworks import FireTaskBase, explicit_serialize
from validation.ecoli.validation_data_raw import ValidationDataRawEcoli

@explicit_serialize
class InitRawValidationDataTask(FireTaskBase):

	_fw_name = "InitRawValidationDataTask"
	required_params = ["output"]

	def run_task(self, fw_spec):
		return