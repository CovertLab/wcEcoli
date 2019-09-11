
import numpy as np
import os
from functools import partial
from reconstruction import spreadsheets
from prototypes.operons import transform_rnas
import csv

'''
Purpose:
This file is analogous to relation.py, but will run from command line, without
having to go through the ParCa. 
Using as a scratch pad to rewrite the functions so that it will deal with
operon data.

TODO: 
-Rewrite file so it will run in the ParCa - ref: relation.py
'''

DIALECT = "excel-tab"
JsonReader = partial(spreadsheets.JsonReader, dialect = DIALECT)
JsonWriter = partial(spreadsheets.JsonWriter, dialect = DIALECT)


FLAT_DIR = os.path.join('reconstruction', 'ecoli', 'flat')
OPERON_RNA_DATA = os.path.join(FLAT_DIR, 'operon_rnas.tsv')
RNA_DATA = os.path.join(FLAT_DIR, 'rnas.tsv')
PROTEIN_DATA = os.path.join(FLAT_DIR, 'proteins.tsv')
output_file = os.path.join(FLAT_DIR, 'rnaid_to_monomerid.tsv')

rna_info = transform_rnas.parse_tsv(RNA_DATA)
operon_rna_info = transform_rnas.parse_tsv(OPERON_RNA_DATA)
protein_info = transform_rnas.parse_tsv(PROTEIN_DATA)

#find all operon_rnas that arent in the rnas.tsv file.




def buildRnaIndexToMonomerMapping():
	'''
	#Old Method: 

	rnaIndexToMonomerMapping = np.array([
		np.where(x == RNA_DATA["id"])[0][0] 
		for x in PROTEIN_DATA["rnaId"]])
	'''
	import pdb; pdb.set_trace()
	rnaIndexToMonomerMapping = np.array([
		np.where(x == rna_info["id"])[0][0] 
		for x in protein_info["rnaId"]])
	

def buildMonomerIndexToRnaMapping(self, raw_data, sim_data):
		monomerIndexToRnaMapping = np.array([
			np.where(x == PROTEIN_DATA["rnaId"])[0][0] 
			for x in RNA_DATA["id"] 
			if len(np.where(x == PROTEIN_DATA["rnaId"])[0])])

buildRnaIndexToMonomerMapping()



'''
from __future__ import division

import re
import numpy as np

# Unit imports
from wholecell.utils import units
from wholecell.utils.unit_struct_array import UnitStructArray

'''
