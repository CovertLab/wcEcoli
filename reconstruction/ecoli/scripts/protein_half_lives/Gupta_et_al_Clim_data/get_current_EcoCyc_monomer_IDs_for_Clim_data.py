# Use this file to covert the UniProt IDs in the Gupta et al. (2024) supplementary
# data files to the current corresponding EcoCyc monomer ID that will be recognized by
# the model using a UniProt ID conversion. This file specifically avoids assgining
# multiple monomer IDs to the same gene common name.

# Note: this file will take a while to run AND requires a number of user inputs.

import pandas as pd
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
import requests
import ast
import io
import time

import pickle
import os

from models.ecoli.analysis import parcaAnalysisPlot
from wholecell.analysis.analysis_tools import exportFigure
from wholecell.utils import constants
from wholecell.utils import units
from wholecell.io import tsv

# USER INPUTS:

# CHANGE THIS TO THE SPECIFIC FILE YOU WANT TO FIND THE MONOMER IDS FOR
# (in the 'Gupta_et_al_2024_data_files' folder)
file_to_convert = '41467_2024_49920_MOESM4_ESM_ST1.xlsx'
# CHANGE THIS TO THE CURRENT DATE (in the form of DDMMYYYY)
date = '11202024'
# CHANGE THIS TO YOUR ECOCYC USERNAME
username = 'miagrahn@stanford.edu' # typically a user's email address
# CHANGE THIS TO YOUR ECOCYC PASSWORD
password = 'RvKrA.%yUr'

# END OF USER INPUTS



# Directory for this file
FILE_LOCATION = os.path.realpath(os.path.dirname(__file__))

# Directory to supplemental data files and data file of interest
INPUT_FOLDER = os.path.join(FILE_LOCATION, 'Gupta_et_al_2024_data_files')
INPUT = os.path.join(INPUT_FOLDER, file_to_convert)

# Make the output file name and specify the output file location
OUTPUT_FILE_NAME = file_to_convert[:-5] + '_EcoCyc_monomer_ID_matches_'+ date + '.tsv'
OUTPUT_FLAT_FILE_PATH = os.path.join(
    'reconstruction/ecoli/scripts/protein_half_lives/Gupta_et_al_Clim_'
    'data/Clim_EcoCyc_monomer_ID_matches', OUTPUT_FILE_NAME)

# read in the data from the table:
FullTable = pd.read_excel(INPUT, skiprows=[0, 1, 2, 3])

# Get the protein IDs from the table
FullTable = FullTable[['Protein ID', 'Gene names ']]
InputTable = FullTable.copy(deep=True)

# Isolate the UniProt ID and gene name from the Protein ID column
InputTable['UniProt ID'] = (
    InputTable['Protein ID'].str.split('|', expand = True))[1]
InputTable['UniProt Gene Name'] = (
    InputTable['Protein ID'].str.split('|', expand = True))[2]

# Create a session for accessing the EcoCyc database
s = requests.Session()
s.post('https://websvc.biocyc.org/credentials/login/',
       data={'email':username, 'password':password})

# Function to grab the EcoCyc monomer ID for a given UniProt ID
def get_ecocyc_id(uniprot_id):
    monomer_id = None
    # Issue web service request:
    url = f"https://websvc.biocyc.org/ECOLI/foreignid?ids=UniProt:{uniprot_id}"
    response = s.get(url).text.split('\t')
    if response[1] == '1':
        monomer_id = response[2].replace("\n", "")
    return monomer_id

# Retrieve the EcoCyc monomer ID for each UniProt ID in the table
# NOTE: this will take a while to run!
InputTable['Monomer ID'] = InputTable['UniProt ID'].apply(get_ecocyc_id)

# Retrieve the common name from rnas.tsv for each monomer ID
# TODO: maybe consider using the gene common names from proteins.tsv instead
rnasTable = pd.read_csv(
    '~/wcEcoli/reconstruction/ecoli/flat/rnas.tsv', sep='\t',
                   skiprows=[0, 1, 2, 3])
rnas_common_names = rnasTable[['monomer_ids','common_name']]

# Function to evaluate strings as lists, as the format of the monomer_ids column in rnas.tsv is a string
def parse_monomer_ids(row):
    try:
        # Evaluate the string as a list
        return ast.literal_eval(row)
    except (ValueError, SyntaxError):
        return []

# Apply the parsing function to the monomer_ids column to make it a real list
rnas_data = rnas_common_names.copy(deep=True)
rnas_data['monomer_ids'] = rnas_data['monomer_ids'].apply(parse_monomer_ids)

# Determine which common names in rnas.tsv correspond to multiple monomer IDs:
common_names_with_multiple_monomer_ids = []
multiple_monomer_ids = []
for i in range(len(rnas_data)):
    if len(rnas_data.iloc[i]['monomer_ids']) > 1:
        multiple_monomer_ids.append(rnas_data.iloc[i]['monomer_ids'])
        common_names_with_multiple_monomer_ids.append(
            rnas_data.iloc[i]['common_name'])

# View which common names correspond to multiple monomer IDs:
multiple_monomer_ids_df = pd.DataFrame(
    {'common_name': common_names_with_multiple_monomer_ids, 'monomer_ids': multiple_monomer_ids})
print("Gene common names with multiple monomer IDs in rnas.tsv:", len(multiple_monomer_ids_df))
print(multiple_monomer_ids_df)

# Make a new table with the updated relevant info:
Temp_InputTable = InputTable[
    ['Protein ID','UniProt ID', 'UniProt Gene Name', 'Monomer ID']]
Interest_InputTable = Temp_InputTable.copy(deep=True)
Interest_InputTable.rename(columns = {'Protein ID':'Gupta et al. 2024 Protein ID'}, inplace = True)

# Add a common name column to be populated with the common name from rnas.tsv
Interest_InputTable['Common Name'] = None

# Save which monomer IDs in Interest_InputTable correspond to multiple
multiple_monomer_ids_with_gene_names = []
gene_names_for_multiple_monomer_ids = []

# Populate the common name column in Interest_InputTable with the common name from rnas.tsv based
# on the UniProt monomer ID from the Gupta et al. dataset:
for i in range(len(rnas_data['monomer_ids'])):
    monomer_id = rnas_data.iloc[i]['monomer_ids']
    if len(monomer_id) > 1:
        # check if the monomer_id is in the Interest_InputTable table
        for j in range(len(monomer_id)):
            if monomer_id[j] in Interest_InputTable['Monomer ID'].values:
                common_name = rnas_common_names.iloc[i]['common_name']
                Interest_InputTable.loc[Interest_InputTable['Monomer ID'] == monomer_id[j],
                'Common Name']  = common_name
                multiple_monomer_ids_with_gene_names.append(monomer_id[j])
                gene_names_for_multiple_monomer_ids.append(rnas_data.iloc[i]['common_name'])
    elif len(monomer_id) == 0:
        # some monomer_ids have no common name, so just skip them
        pass
    else:
        # cases where there is one monomer ID to one common name
        monomer_id = monomer_id[0]
        common_name = rnas_common_names.iloc[i]['common_name']
        if monomer_id in Interest_InputTable['Monomer ID'].values:
            Interest_InputTable.loc[Interest_InputTable['Monomer ID'] == monomer_id,
            'Common Name']  = common_name

# Save a copy of the common names in rnas.tsv that correspond to multiple monomer IDs that also show
# up in the Interest_InputTable monomer ID list:
multiple_monomer_ids_with_gene_names_df = (
    pd.DataFrame({'monomer_ids': multiple_monomer_ids_with_gene_names,
                  'common_name': gene_names_for_multiple_monomer_ids}))
print("Common names with multiple monomer IDs in rnas.tsv that match to monomer"
      " IDs in the Gupta et al. dataset:",
      len(multiple_monomer_ids_with_gene_names_df))
print(multiple_monomer_ids_with_gene_names_df)

# CHECK TO MAKE SURE ALL THE NUMBERS ADD UP AS EXPECTED!
print("Total number of monomer IDs in the dataset: ", len(Interest_InputTable['Monomer ID']))

# Check if any common_names in Interest_InputTable show up for more than one Monomer ID
non_unique_common_names = []
for common_name in range(len(Interest_InputTable['Common Name'])):
    common_name = Interest_InputTable.iloc[common_name]['Common Name']
    if len(Interest_InputTable[Interest_InputTable['Common Name'] == common_name]) > 1:
        non_unique_common_names.append(common_name)
# NOTE: this should come out to zero if the code is working correctly
print("# of common names that show up for more than one monomer ID in the dataset: ",
      len(non_unique_common_names)) # 0

# Check how monomer IDs do not share a common name with another monomer ID:
unique_common_names = []
for common_name in range(len(Interest_InputTable['Common Name'])):
    common_name = Interest_InputTable.iloc[common_name]['Common Name']
    if len(Interest_InputTable[Interest_InputTable['Common Name'] == common_name]) == 1:
        unique_common_names.append(common_name)
print("# of monomer IDs that do not share a common name with another monomer ID: ",
      len(unique_common_names)) # 2354

# Check how many monomer IDs do not have a common name:
none_common_names = 0
for common_name in range(len(Interest_InputTable['Common Name'])):
    common_name = Interest_InputTable.iloc[common_name]['Common Name']
    if common_name == None:
        none_common_names += 1
print(none_common_names) # 8

# Check that everything adds up to the total number of monomer IDs in the dataset
total = (
        len(non_unique_common_names) + len(unique_common_names) + none_common_names)
if total == len(Interest_InputTable['Monomer ID']):
    print("All the numbers add up to the total number of monomer IDs in the dataset: ",
          str(total))
else:
    print("The numbers do not add up, please check the files manually.")

# Save the Interest_InputTable to a .csv file:
print("Saving .csv file...")
with io.open(OUTPUT_FLAT_FILE_PATH, 'wb') as f:
    writer = tsv.writer(f, quotechar="'", lineterminator='\n')
    writer.writerow(['# Generated by {} on {}'.format(__file__, time.ctime())])
    writer.writerow(
        ['Gupta et al. 2024 Protein ID','UniProt ID', 'UniProt Gene Name',
         'Monomer ID', 'Common Name'])

    # save the Interest_InputTable as a .tsv file
    for i in range(len(Interest_InputTable)):
        row = Interest_InputTable.iloc[i]
        writer.writerow(
            [f'"{row[0]}"', f'{row[1]}', f'{row[2]}', f'{row[3]}', f'{row[4]}'])


#Interest_InputTable.to_csv(OUTPUT_FLAT_FILE_PATH, sep='\t', index=False)
print("Done. File saved to: ", OUTPUT_FLAT_FILE_PATH)








