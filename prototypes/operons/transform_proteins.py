import os
import numpy as np
from functools import partial
from reconstruction import spreadsheets
import csv


'''
Purpose:
Here we are simply taking in proteins.tsv and appending a column to the 
file called rna_sets which is intended to contain all the RNA's attached 
to a specific monomer in a list.

This file will take in data from operon_rnas.tsv and look for all 
RNAs that are assigned to a specific monomer through the monomer_sets column
of that file.

Currently saves to a new filename. This filename is manually changed after 
it has been verified to be generating the correct information, to prevent
overwriting the original data.

TODO:
In future need to make sure that operon_rnas.tsv is automatically rather than
semi manuallly generated as it is now. Some processing for that file
occurs in transform_rnas.
'''


DIALECT = "excel-tab"
JsonReader = partial(spreadsheets.JsonReader, dialect = DIALECT)
JsonWriter = partial(spreadsheets.JsonWriter, dialect = DIALECT)


FLAT_DIR = os.path.join("reconstruction", "ecoli", "flat")
RNA_FILE = os.path.join(FLAT_DIR, "operon_rnas.tsv")
PROTEIN_FILE = os.path.join(FLAT_DIR, "proteins_old.tsv")
#saving to a new file for now so that all manually input TUs in 
#operon_rnas.tsv are not overwritten.
output_file = os.path.join(FLAT_DIR, "proteins_1.tsv")
#


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


def make_collection():
	protein_info, protein_fieldnames = parse_tsv(PROTEIN_FILE)
	rna_info, rna_fieldnames = parse_tsv(RNA_FILE)

	#Go through monomerSet line by line. Find the matching monomers within
	#those lists then find the corresponding monomer in proteins.tsv.
	#Add the id from operon_rnas to the rnaSet list

	protein_index = {}
	for protein_row in protein_info:
		protein_row['rnaSet'] = []
		protein_index[protein_row['id']] = protein_row

	
	for rna_row in rna_info:
		for monomer in rna_row['monomerSet']:
			protein_row = protein_index[monomer]
			protein_row['rnaSet'].append(rna_row['id'])	

	
	#add fieldname for 'monomersets'
	protein_fieldnames.append('rnaSet')

	with open(output_file, "w") as f:
		writer = JsonWriter(f, protein_fieldnames)
		writer.writeheader()
		for protein_row in protein_info:
			writer.writerow(protein_row)


if __name__ == "__main__":
	make_collection()
