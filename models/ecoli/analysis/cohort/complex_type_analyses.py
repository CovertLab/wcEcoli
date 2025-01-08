"""
Template for cohort analysis plots
"""

import pickle
import os
import pandas as pd
import ast
from matplotlib import pyplot as plt
# noinspection PyUnresolvedReferences
import numpy as np

from models.ecoli.analysis import cohortAnalysisPlot
from wholecell.analysis.analysis_tools import (exportFigure,
	read_bulk_molecule_counts, read_stacked_bulk_molecules, read_stacked_columns)
from wholecell.io.tablereader import TableReader

# path to the complexes file:
complexes_path = '~/wcEcoli/reconstruction/ecoli/flat/complexation_reactions.tsv'

# read in the complexes file
complexes = pd.read_csv(complexes_path, skiprows=4, sep='\t')
import json

# Clean up the stoichiometry column
complexes['stoichiometry'] = complexes['stoichiometry'].str.replace(r'[^\x00-\x7F]+', '', regex=True)  # Remove non-ASCII characters
complexes['stoichiometry'] = complexes['stoichiometry'].str.strip('"')  # Remove surrounding quotes if present

# Try using json.loads to convert the string to a dictionary
complexes['stoichiometry'] = complexes['stoichiometry'].apply(lambda x: json.loads(x) if isinstance(x, str) else x)

# Check the result
print(complexes['stoichiometry'].head())
hi = 3


def generate_data(complexes_table):
	# extract the "stoichiometry" column from the complexes table:
	stoichiometry = complexes_table['stoichiometry']

	complex_names = []
	for row in stoichiometry:
		name = list(row.keys())[0]
		keys = list(row.keys())
		complex_names.append(name)

	monomers = []
	monomer_stoichiometry = []
	monomer_number = []
	for row in stoichiometry:
		current_monomers = []
		current_stoichiometry = []
		for key in list(row.keys())[1:]:
			current_monomers.append(key)
			current_stoichiometry.append(row[key])
		monomers.append(current_monomers)
		monomer_stoichiometry.append(current_stoichiometry)
		monomer_number.append(len(current_monomers))

	hi = 5



	return complexes_table


complexes_table = generate_data(complexes)

