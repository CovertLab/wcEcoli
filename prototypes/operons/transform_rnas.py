import os
from functools import partial

from reconstruction import spreadsheets

'''
There is an error in this file when there is no monomerID the 
value is not saved.

Fix this

'''
import csv

DIALECT = "excel-tab"
JsonReader = partial(spreadsheets.JsonReader, dialect = DIALECT)
JsonWriter = partial(spreadsheets.JsonWriter, dialect = DIALECT)


FLAT_DIR = os.path.join("reconstruction", "ecoli", "flat")
RNA_FILE = os.path.join(FLAT_DIR, "rnas.tsv")
#saving to a new file for now so that all manually input TUs in 
#operon_rnas.tsv are not overwritten.
output_file = os.path.join(FLAT_DIR, "operon_rnas_3.tsv")


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
	#pull in rna_info as a list
	rna_info, fieldnames= parse_tsv(RNA_FILE)

	#look for instances where rnas do not have an assigned monomerId
	#replace with a empty list, else put monomerId(s) into a list.
	#Do the same for GeneIDs

	for rna_row in rna_info:
		if not rna_row['monomerId']:
			rna_row['monomerId'] = []
		else:
			rna_row['monomerId'] = [rna_row['monomerId']]

		#if not rna_row['geneId']:
			#rna_row['geneId'] = []
		#else:
			#rna_row['geneId'] = [rna_row['geneId']]

	with open(output_file, "w") as f:
		writer = JsonWriter(f, fieldnames)
		writer.writeheader()
		for rna_row in rna_info:
			writer.writerow(rna_row)


if __name__ == "__main__":
	make_collection()

