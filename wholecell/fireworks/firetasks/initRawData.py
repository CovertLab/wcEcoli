from __future__ import absolute_import, division, print_function

from six.moves import cPickle
import time

from fireworks import FiretaskBase, explicit_serialize
from reconstruction.ecoli.knowledge_base_raw import KnowledgeBaseEcoli
from wholecell.utils.constants import DEFAULT_OPERON_OPTION
from wholecell.utils.constants import DEFAULT_NEW_GENES_OPTION


@explicit_serialize
class InitRawDataTask(FiretaskBase):

	_fw_name = "InitRawDataTask"
	required_params = ["output"]
	optional_params = ['operons','new_genes']

	def run_task(self, fw_spec):
		operon_option = self.get('operons') or DEFAULT_OPERON_OPTION
		print(f"{time.ctime()}: Instantiating raw_data with operons={operon_option}")

		new_gene_option = self.get('new_genes') or DEFAULT_NEW_GENES_OPTION
		print(f"{time.ctime()}: Instantiating raw_data with new_genes={new_gene_option}")

		raw_data = KnowledgeBaseEcoli(
			operons_on=(operon_option == 'on'),
			new_genes_on=(new_gene_option == 'on'))

		print(f"{time.ctime()}: Saving raw_data")

		with open(self["output"], "wb") as f:
			cPickle.dump(raw_data, f, protocol = cPickle.HIGHEST_PROTOCOL)
