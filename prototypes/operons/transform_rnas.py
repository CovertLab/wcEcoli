import os
from functools import partial
from reconstruction import spreadsheets
import csv

DIALECT = "excel-tab"
JsonReader = partial(spreadsheets.JsonReader, dialect = DIALECT)
JsonWriter = partial(spreadsheets.JsonWriter, dialect = DIALECT)


FLAT_DIR = os.path.join("reconstruction", "ecoli", "flat")
RNA_FILE = os.path.join(FLAT_DIR, "rnas.tsv")
#saving to a new file for now so that all manually input TUs in 
#operon_rnas.tsv are not overwritten.
output_file = os.path.join(FLAT_DIR, "operon_rnas_3.tsv")
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
	#pull in rna_info as a list
	rna_info, fieldnames = parse_tsv(RNA_FILE)
	
	#Creating a new column for monomerSets. Which are the monomers for given
	#mRNA transcripts. Right now we are manually adding in multi-gene TUs, 
	#so Im not going to worry about anything except taking the current rnas.tsv
	#and creating an output that has a new column putting the monomerID into a list.
	#for a multigene TU, will have two IDs within this list.

	#look for instances where rnas do not have an assigned monomerId
	#replace with a empty list, else put monomerId(s) into a list.

	

	for rna_row in rna_info:
		if not rna_row['monomerId']:
			rna_row['monomerSet'] = []
		else:
			rna_row['monomerSet'] = [rna_row['monomerId']]

	#add fieldname for 'monomersets'
	fieldnames.append('monomerSet')

	with open(output_file, "w") as f:
		writer = JsonWriter(f, fieldnames)
		writer.writeheader()
		for rna_row in rna_info:
			writer.writerow(rna_row)


if __name__ == "__main__":
	make_collection()


