from __future__ import absolute_import, division, print_function


from reconstruction import spreadsheets
from functools import partial
from collections import defaultdict

DIALECT = "excel-tab"
JsonReader = partial(spreadsheets.JsonReader, dialect = DIALECT)
JsonWriter = partial(spreadsheets.JsonWriter, dialect = DIALECT)

def get_operon_rna_list(pc_rnas_path):


    rna_list, fieldnames = parse_tsv(pc_rnas_path)

    pc_dict = {}
    for row in rna_list:
        operon_rna = '_'.join(row['transcription_units']) + '_RNA[c]'
        pc_dict[operon_rna] = [s + '_RNA[c]' for s in row['transcription_units']]

        monomer_rna_included = set(row['transcription_units']) - set(row['monomers_to_remove'])

        if len(monomer_rna_included) > 0:
            for mono_rna in monomer_rna_included:
                pc_dict[mono_rna + '_RNA[c]'] = mono_rna + '_RNA[c]'

    return pc_dict

def pc_mc_conversion_dicts(pc_rnas_path):
    rna_list, fieldnames = parse_tsv(pc_rnas_path)

    pc_to_mc_dict = {}
    for row in rna_list:
        operon_rna = '_'.join(row['transcription_units']) + '_RNA[c]'
        pc_to_mc_dict[operon_rna] = [s + '_RNA[c]' for s in row['transcription_units']]

    mc_to_pc_dict = defaultdict(list)

    for key in pc_to_mc_dict.keys():
        for v in pc_to_mc_dict[key]:
            mc_to_pc_dict[v].append(key)

    return pc_to_mc_dict, mc_to_pc_dict


def parse_tsv(tsv_file):
#Takes in a tsv file, and creates a list of lists of the rows
#contained within the TSV.
	tsv_list = []
	with open(tsv_file) as tsvfile:
		reader = JsonReader(tsvfile)
		fieldnames = reader.fieldnames
		for row in reader:
			tsv_list.append(row)
	return tsv_list, fieldnames

