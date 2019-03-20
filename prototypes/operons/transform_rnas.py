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
output_file = os.path.join(FLAT_DIR, "operon_rnas_2.tsv")


def parse_tsv(tsv_file):
	
#Takes in a tsv file, and creates a list of lists of the rows 
#contained within the TSV.
	
	tsv_list = []
	with open(tsv_file) as tsvfile:
		reader = JsonReader(tsvfile)
		fieldnames = reader.fieldnames
		#reader.next()
		for row in reader:
			tsv_list.append(row)
	return tsv_list, fieldnames


def make_collection():
	#pull in rna_info as a list
	rna_info, fieldnames= parse_tsv(RNA_FILE)

	#look for instances where rnas do not have an assigned monomerId
	#replace with a empty list, else put monomerId(s) into a list.

	for rna_row in rna_info:
		if not rna_row['monomerId']:
			rna_row['monomerId'] = []
		else:
			rna_row['monomerId'] = [rna_row['monomerId']]


	#now need to write to a file:
	#Tried doing with Json writer but had an issue bringin in the header.
	with open (output_file, 'w') as f:
		writer = csv.writer(f, dialect = DIALECT, quoting=csv.QUOTE_MINIMAL)
		writer.writerow(fieldnames)	
		for rna_row in rna_info:
			rna_data = list(rna_row.values())
			writer.writerow(rna_data)

	#import ipdb; ipdb.set_trace()
'''


def make_collection():
	with open(RNA_FILE, "r") as f:
		reader = JsonReader(f)
		entries = {entry["monomerId"]: entry for entry in reader} #this is the line that is excluting information.
		fieldnames = reader.fieldnames
		for entry in entries.values():
			entry["monomerId"] = [entry["monomerId"]]
	with open(output_file, "w") as f:
		writer = JsonWriter(f, fieldnames)
		writer.writeheader()

		for key in sorted(entries.keys()):
			writer.writerow(entries[key])
'''
if __name__ == "__main__":
	make_collection()

