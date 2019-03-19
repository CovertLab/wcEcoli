import os
from functools import partial

from reconstruction import spreadsheets

'''
There is an error in this file when there is no monomerID the 
value is not saved.

Fix this error today.

'''

DIALECT = "excel-tab"
JsonReader = partial(spreadsheets.JsonReader, dialect = DIALECT)
JsonWriter = partial(spreadsheets.JsonWriter, dialect = DIALECT)


FLAT_DIR = os.path.join("reconstruction", "ecoli", "flat")
RNA_FILE = os.path.join(FLAT_DIR, "operon_rnas.tsv")
output_file = os.path.join(FLAT_DIR, "operon_rnas_2.tsv")

def make_collection():
	with open(RNA_FILE, "r") as f:
		reader = JsonReader(f)
		entries = {entry["monomerId"]: entry for entry in reader}
		fieldnames = reader.fieldnames
		for entry in entries.values():
			entry["monomerId"] = [entry["monomerId"]]

	with open(output_file, "w") as f:
		writer = JsonWriter(f, fieldnames)
		writer.writeheader()

		for key in sorted(entries.keys()):
			writer.writerow(entries[key])

if __name__ == "__main__":
	make_collection()