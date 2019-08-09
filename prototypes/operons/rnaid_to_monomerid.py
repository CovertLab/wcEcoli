
import os
from functools import partial
from reconstruction import spreadsheets
from prototypes.operons.transform_rnas import parse_tsv
import csv

DIALECT = "excel-tab"
JsonReader = partial(spreadsheets.JsonReader, dialect = DIALECT)
JsonWriter = partial(spreadsheets.JsonWriter, dialect = DIALECT)


FLAT_DIR = os.path.join('reconstruction', 'ecoli', 'flat')
OPERON_RNA_DATA = os.path.join(FLAT_DIR, 'operon_rnas.tsv')
RNA_DATA = os.path.join(FLAT_DIR, 'rnas.tsv')
PROTEIN_DATA = os.path.join(FLAT_DIR, 'proteins.tsv')
output_file = os.path.join(FLAT_DIR, 'rnaid_to_monomerid.tsv')

#find all operon_rnas that arent in the rnas.tsv file.


def buildRnaIndexToMonomerMapping(self, raw_data, sim_data):
	rnaIndexToMonomerMapping = np.array([
		np.where(x == RNA_DATA["id"])[0][0] 
		for x in PROTEIN_DATA["rnaId"]])

def buildMonomerIndexToRnaMapping(self, raw_data, sim_data):
		monomerIndexToRnaMapping = np.array([
			np.where(x == PROTEIN_DATA["rnaId"])[0][0] 
			for x in RNA_DATA["id"] 
			if len(np.where(x == PROTEIN_DATA["rnaId"])[0])])



#need to pull in the RNA file and

'''
from __future__ import division

import re
import numpy as np

# Unit imports
from wholecell.utils import units
from wholecell.utils.unit_struct_array import UnitStructArray

'''
