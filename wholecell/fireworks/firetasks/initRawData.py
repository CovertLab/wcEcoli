from __future__ import absolute_import, division, print_function

import pickle
import time

from fireworks import FiretaskBase, explicit_serialize
from reconstruction.ecoli.knowledge_base_raw import KnowledgeBaseEcoli
from wholecell.utils.constants import DEFAULT_OPERON_OPTION


@explicit_serialize
class InitRawDataTask(FiretaskBase):

	_fw_name = "InitRawDataTask"
	required_params = ["output"]
	optional_params = [
		'operons',
		'remove_rrna_operons',
		'remove_rrff',
		]

	def run_task(self, fw_spec):
		operon_option = self.get('operons') or DEFAULT_OPERON_OPTION
		print(f"{time.ctime()}: Instantiating raw_data with operons={operon_option}")

		raw_data = KnowledgeBaseEcoli(
			operons_on=(operon_option == 'on'),
			remove_rrna_operons=self.get('remove_rrna_operons', False),
			remove_rrff=self.get('remove_rrff', False),
			)

		print(f"{time.ctime()}: Saving raw_data")

		with open(self["output"], "wb") as f:
			pickle.dump(raw_data, f, protocol = pickle.HIGHEST_PROTOCOL)
