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

# clean up the stoichiometry column to be parsed correctly
complexes['stoichiometry'] = complexes['stoichiometry'].str.replace(r'[^\x00-\x7F]+', '', regex=True)
complexes['stoichiometry'] = complexes['stoichiometry'].str.strip('"')

# use json.loads to convert the string to a dictionary
complexes['stoichiometry'] = complexes['stoichiometry'].apply(lambda x: json.loads(x) if isinstance(x, str) else x)

def generate_full_stoichiometry_table(complexes_table):
	# extract the "stoichiometry" column from the complexes table:
	stoichiometry = complexes_table['stoichiometry']

	complex_names = []
	for row in stoichiometry:
		name = list(row.keys())[0]
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

	# create a table with the complex names, monomers, and stoichiometry
	table = pd.DataFrame({'complex_name': complex_names, 'monomers': monomers, 'stoichiometry': monomer_stoichiometry, 'unique_monomers': monomer_number})

	# add an extra column to the table with the complex's classification:
	table['complex_type'] = ''

	# loop throught the complex name column and classify the complex type based on the number of monomers in it
	for i in range(len(table)):
		if table.loc[i, 'unique_monomers'] == 0:
			table.loc[i, 'complex_type'] = 'unknown'
		elif table.loc[i, 'unique_monomers'] == 1:
			table.loc[i, 'complex_type'] = 'homogeneous'
		else:
			table.loc[i, 'complex_type'] = 'heterogeneous'

	# check the result by printing the number of homogeneous and heterogeneous complexes
	print('Number of homogeneous complexes:', len(table[table['complex_type'] == 'homogeneous']))
	print('Number of heterogeneous complexes:', len(table[table['complex_type'] == 'heterogeneous']))
	print('Number of unknown complexes:', len(table[table['complex_type'] == 'unknown']))

	return table

# make a function that finds the monomers that exist only in a homogeneous complex
def generate_homogeneous_monomer_table(full_complexes_table):
	"""
	Find monomers that exist ONLY in a homogeneous complex(s)
	Args:
		full_complexes_table: A table consisting of each complex and its contents

	Returns:
		monomers_in_homogeneous_complexes: A table of monomers that exist only
		in a homogeneous complex, along with the identity of the complex and
		the number of monomers in the complex
	"""
	# create a copy of the full_complexes_table that contains only the homogeneous complexes
	homogeneous_complexes = full_complexes_table[full_complexes_table['complex_type'] == 'homogeneous'].copy()
	homogeneous_complexes.reset_index(drop=True, inplace=True)

	# create a new table to store monomers that are only in homogeneous complexes
	monomers_ids = []
	complex_ids = []
	stoichiometries = []
	for i in range(len(homogeneous_complexes)):
		# extract the monomers and stoichiometry for the current complex
		monomer = homogeneous_complexes.loc[i, 'monomers']
		complex_id = homogeneous_complexes.loc[i, 'complex_name']
		stoichiometry = homogeneous_complexes.loc[i, 'stoichiometry']
		# append the monomers, complex_id, and stoichiometry to the new table
		monomers_ids.append(monomer[0])
		complex_ids.append(complex_id)
		# get the stoichiometry value (some times it is "None")
		stoich = stoichiometry[0]
		if type(stoich) == int:
			stoich = -1 * stoich
			stoichiometries.append(stoich)
		else:
			stoichiometries.append(stoichiometry)

	# make the table
	monomers_in_homogeneous_complexes = pd.DataFrame(
			{'monomer_id': monomers_ids, 'complex_id': complex_ids,
			 'stoichiometry': stoichiometries})

	duplicates = monomers_in_homogeneous_complexes['monomer_id'].duplicated()
	print(monomers_in_homogeneous_complexes[duplicates])

	return monomers_in_homogeneous_complexes

# todo: generate a table of monomers in heterogeneous compelexes
# todo: generate a table of monomers in homogeneous compelexes
# todo: generate a table of monomers in both homogeneous and heterogeneous complexes
# todo: generate a table of monomers in heterogeneous compelexes only
# todo: generate a table of monomers in homogeneous compelexes only
# todo: add the added, change the modified, and remove the removed monomers from the starting table
# todo: ask nora where exactly the random complexed free monomers are stored?

full_complexes_table = generate_full_stoichiometry_table(complexes)
monomers_in_homogeneous_complexes = generate_homogeneous_monomer_table(full_complexes_table)


# save the tables as csv files

outpath = os.path.join('~/wcEcoli/out/', 'complex_classification_tables')
outpath = os.path.expanduser('~/wcEcoli/out/complex_classification_tables')
if not os.path.exists(outpath):
	os.makedirs(outpath)
	print('Directory created:', outpath)
full_complexes_table.to_csv(os.path.join(outpath, 'all_complexes.csv'))
monomers_in_homogeneous_complexes.to_csv(os.path.join(outpath, 'monomers_in_homogeneous_complexes.csv'))



