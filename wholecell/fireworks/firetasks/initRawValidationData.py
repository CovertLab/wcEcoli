import pickle
import time

from fireworks import FiretaskBase, explicit_serialize
from validation.ecoli.validation_data_raw import ValidationDataRawEcoli


@explicit_serialize
class InitRawValidationDataTask(FiretaskBase):

	_fw_name = "InitRawValidationDataTask"
	required_params = ["output"]

	def run_task(self, fw_spec):
		print("%s: Instantiating validation_data_raw" % (time.ctime()))

		validation_data_raw = ValidationDataRawEcoli()

		print("%s: Saving validation_data_raw" % (time.ctime()))

		with open(self["output"], "wb") as fh:
			pickle.dump(validation_data_raw, fh, protocol=pickle.HIGHEST_PROTOCOL)
