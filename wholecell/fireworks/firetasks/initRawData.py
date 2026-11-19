import pickle
import time

from fireworks import FiretaskBase, explicit_serialize
from reconstruction.ecoli.knowledge_base_raw import KnowledgeBaseEcoli
from wholecell.utils.constants import DEFAULT_OPERON_OPTION
from wholecell.utils.constants import DEFAULT_NEW_GENES_OPTION
from wholecell.utils.constants import DEFAULT_PROTEIN_DEGRADATION_COMBO


@explicit_serialize
class InitRawDataTask(FiretaskBase):

	_fw_name = "InitRawDataTask"
	required_params = ["output"]
	optional_params = [
		'operons',
		'new_genes',
		'protein_degradation_combo',
		'remove_rrna_operons',
		'remove_rrff',
		'stable_rrna',
		]

	def run_task(self, fw_spec):
		operon_option = self.get('operons') or DEFAULT_OPERON_OPTION
		print(f"{time.ctime()}: Instantiating raw_data with operons={operon_option}")

		new_gene_option = self.get('new_genes') or DEFAULT_NEW_GENES_OPTION
		print(f"{time.ctime()}: Instantiating raw_data with new_genes={new_gene_option}")

		protein_degradation_combo = self.get('protein_degradation_combo') or DEFAULT_PROTEIN_DEGRADATION_COMBO
		print(f"{time.ctime()}: Instantiating raw_data with protein_degradation_combo={protein_degradation_combo}")

		raw_data = KnowledgeBaseEcoli(
			operons_on=(operon_option == 'on'),
			new_genes_option=new_gene_option,
			protein_degradation_combo_option=protein_degradation_combo,
			remove_rrna_operons=self.get('remove_rrna_operons', False),
			remove_rrff=self.get('remove_rrff', False),
			stable_rrna=self.get('stable_rrna', False),
			)

		print(f"{time.ctime()}: Saving raw_data")

		with open(self["output"], "wb") as f:
			pickle.dump(raw_data, f, protocol = pickle.HIGHEST_PROTOCOL)
