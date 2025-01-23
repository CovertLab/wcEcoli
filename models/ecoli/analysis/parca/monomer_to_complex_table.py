"""
This script generates a tsv file where every monomer in the model is reported along with which complexes it is found in,
the counts of the specified monomers per complex, monomer protease assignment, free monomer half life specified in
given parca and the source of this degradation rate.
"""

import pickle
import os

from matplotlib import pyplot as plt
# noinspection PyUnresolvedReferences
import numpy as np
import csv

from models.ecoli.analysis import parcaAnalysisPlot
from wholecell.utils import units

complexToMonomer = {
	"CPLX0-7620[c]": "PD00260[c]", # CPLX0-7620's monomer is EG10359-MONOMER, which is ID'ed as PD00260 (proteins.tsv)
	"CPLX0-8801[c]": "G6420-MONOMER[c]",
	"CPLX0-7677[c]": "EG11969-MONOMER[c]",
	"CPLX0-7702[c]": "G6263-MONOMER[c]",
	"CPLX0-7701[c]": "G6420-MONOMER[c]",
	}

monomerToTranslationMonomer = {
	"MONOMER0-1781[c]": "EG11519-MONOMER[c]", # MONOMER0-1781 is a complex, EG11519 is its monomer
	"EG10359-MONOMER[c]": "PD00260[c]", # EG10359-MONOMER is not the ID of fur monomer, it's PD00260 (proteins.tsv)
	}

class Plot(parcaAnalysisPlot.ParcaAnalysisPlot):

	def do_plot(self, input_dir, plot_out_dir, plot_out_filename, sim_data_file, validation_data_file, metadata):
		with open(sim_data_file, 'rb') as f:
			sim_data = pickle.load(f)
		monomer_ids = sim_data.process.translation.monomer_data['id']
		monomer_id_to_index = {
			monomer_id: i for (i, monomer_id)
			in enumerate(monomer_ids)
		}
		protease = sim_data.process.translation.monomer_data['protease_assignment']
		deg_rate_source = sim_data.process.translation.monomer_data['deg_rate_source']
		degradation_rates = sim_data.process.translation.monomer_data['deg_rate'].asNumber(1 / units.s)
		half_lives = np.log(2) / degradation_rates / 60  # in minutes

		# Get all protein ids reqiured
		ids_complexation_complexes = sim_data.process.complexation.ids_complexes  # Only complexes
		ids_two_comp_sys_complexes = list(
			sim_data.process.two_component_system.complex_to_monomer.keys())  # complexes, generally monomers that are phosphorylated
		ids_equilibrium_complexes = sim_data.process.equilibrium.ids_complexes  # Only complexes

		sim_complexation = sim_data.process.complexation
		sim_equilibrium = sim_data.process.equilibrium
		sim_two_comp = sim_data.process.two_component_system


		def find_monomers_mult_complex(complex_ids, sim_data_process):
			# Identify monomers that are in complexes
			monomersInvolvedInManyComplexes = []
			monomersInvolvedInComplexes = []
			monomer_complex_dict = {}

			for complexId in complex_ids:
				#account for number of the given monomer per complex
				subunitIds = sim_data_process.get_monomers(complexId)["subunitIds"]
				subunit_stoic = sim_data_process.get_monomers(complexId)["subunitStoich"]
				for i, subunitId in enumerate(subunitIds):
					monomers_per_complex = subunit_stoic[i]
					if subunitId not in monomer_id_to_index:
						if subunitId in monomerToTranslationMonomer:
							# couple monomers have different ID in ids_translation
							subunitId = monomerToTranslationMonomer[subunitId]
						elif "CPLX" in subunitId:
							# few transcription factors are complexed with ions
							if subunitId in complexToMonomer:
								subunitId = complexToMonomer[subunitId]
						elif "RNA" in subunitId:
							continue
					if subunitId in monomer_id_to_index:
						if subunitId in monomersInvolvedInComplexes:
							monomersInvolvedInManyComplexes.append(subunitId)
						if subunitId not in monomer_complex_dict:
							monomer_complex_dict[subunitId] = [(complexId, monomers_per_complex)]
						else:
							monomer_complex_dict[subunitId].append((complexId, monomers_per_complex))
						monomersInvolvedInComplexes.append(subunitId)

			return monomer_complex_dict

		cl_monomers_comp_dict = find_monomers_mult_complex(
			ids_complexation_complexes, sim_complexation)
		eq_monomers_comp_dict = find_monomers_mult_complex(
			ids_equilibrium_complexes, sim_equilibrium)
		two_comp_monomers_comp_dict = find_monomers_mult_complex(
			ids_two_comp_sys_complexes, sim_two_comp)


		# Some monomers are listed in more than one complex-associated dictionary
		def merge_dicts_repeat_vals(dict1, dict2):
			"""
			Merges two dictionaries with list values.
			If a key exists in both dictionaries, the values (which are lists) are appended.
			"""
			merged_dict = {}

			for key, value in dict1.items():
				if key in dict2:
					merged_dict[key] = value + dict2[key]
				else:
					merged_dict[key] = value

			for key, value in dict2.items():
				if key not in merged_dict:
					merged_dict[key] = value

			return merged_dict

		cl_eq_comp_dict = merge_dicts_repeat_vals(cl_monomers_comp_dict, eq_monomers_comp_dict)
		total_comp_dict = merge_dicts_repeat_vals(cl_eq_comp_dict, two_comp_monomers_comp_dict)

		# Remove any duplicate complexes in the values
		clean_comp_dict = {}
		for monomer, complex_list in total_comp_dict.items():
			clean_comp_dict[monomer] = list(set(complex_list))

		all_monomers_set = set(monomer_id_to_index.keys())

		monomers_not_in_complex_set = all_monomers_set.symmetric_difference(set(clean_comp_dict.keys()))

		# Add monomers that do not form complexes to the dictionary
		for monomer in monomers_not_in_complex_set:
			clean_comp_dict[monomer] = [(None, None)]


		def write_table_for_monomers_in_complexes(clean_comp_dict, monomer_id_to_index, protease, half_lives,
												  deg_rate_source):
			# Write data to table
			with open(os.path.join(plot_out_dir, plot_out_filename + '.tsv'),
					  'w') as f:
				writer = csv.writer(f, delimiter='\t')
				writer.writerow(
					['monomer_id', 'complex_id', 'monomers_per_complex', 'monomer_protease_assignment', 'total_monomer_half_life_(min)',
					 'deg_rate_source'])

				for monomer, complex_list in clean_comp_dict.items():
					monomer_index = monomer_id_to_index[monomer]
					protease_assign = protease[monomer_index]
					half_life = half_lives[monomer_index]
					source = deg_rate_source[monomer_index]
					for complex, monomer_number in complex_list:
						writer.writerow([
							monomer, complex, monomer_number, protease_assign, half_life, source
						])

		write_table_for_monomers_in_complexes(clean_comp_dict, monomer_id_to_index, protease, half_lives,
											  deg_rate_source)

if __name__ == "__main__":
	Plot().cli()
