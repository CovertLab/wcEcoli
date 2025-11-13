"""
Converts sorted files of protease assignments selected from supplementary table
 S2 of Gupta et al. 2024,
"Global Protein-Turnover Quantification in Escherichia coli Reveals Cytoplasmic
Recycling under Nitrogen Limitation"
(doi: https://doi.org/10.1038/s41467-024-49920-8) to flat files to be used in the
translation.py within the model.

NOTE:
After the new flat file is generated, it must also be added to the list in
knowledge_base_raw.py in order for it to be used in the model.

NOTE: This script uses the EcoCyc monomer IDs that correspond to the UniProt IDs
provided in the Gupta et. al. (2024) dataset. The EcoCyc monomer IDs may update
over time when the model is updated with a new version of EcoCyc, so it is
important to ensure that the EcoCyc monomer IDs are up-to-date from time to time
and generating new EcoCyc monomer ID match files as needed by using the script:
get_current_EcoCyc_monomer_IDs_for_Clim_data.py (in the
uniprot_id_to_ecocyc_monomer_id_matches folder).
"""

import io
import os
import time
import pandas as pd
from wholecell.io import tsv
from wholecell.utils.filepath import ROOT_PATH
os.chdir(os.path.expanduser('~/wcEcoli/'))

"""
USER INPUTS
"""

# PATH TO THE SORTED INPUT FILE (in the
# reconstruction/ecoli/flat/clim_half_life_data/raw_protease_assignment_files folder)
# Note: this file should have a column named "Protein ID" (formatted exactly as
# it is in the original file: 41467_2024_49920_MOESM4_ESM_ST1.xlsx)
# and a column named "half_life (units.min)"
protease_assignment_file = 'priority_proteases_10182024.xlsx'

# SPECIFY THE ECOCYC COMPARISON FILE (in the
# uniprot_id_to_ecocyc_monomer_id_matches folder)
# Note: if a new EcoCyc update has been released, use
# get_current_EcoCyc_monomer_IDs_for_protease_assignments.py.py
# to get an updated comparison file
EcoCyc_file = \
    '41467_2024_49920_MOESM5_ESM_ST2_EcoCyc_monomer_ID_matches_11122025.tsv'

# SPECIFY THE OUTPUT FILE NAME (to be placed in the flat file list)
# Note: this file name needs to be added to the list in "knowledge_base_raw.py"
# to be used in the model
OUTPUT_FILE_NAME = 'protease_assignments.tsv' # e.g. 'protease_assignments_Clim#.tsv'

"""
END OF USER INPUTS
"""

# file folder location:
raw_data_file_location = 'reconstruction/ecoli/flat/clim_half_life_data/'

# Path to the input file location:
INPUT_PATH = os.path.join(
    raw_data_file_location, 'raw_protease_assignment_files', protease_assignment_file)

# Path to the EcoCyc comparison file location:
ECOCYC_PATH = os.path.join(
    raw_data_file_location, 'uniprot_id_to_ecocyc_monomer_id_matches', EcoCyc_file)

# Path to the output file location:
OUTPUT_PATH = os.path.join(ROOT_PATH, 'reconstruction', 'ecoli', 'flat',
                           OUTPUT_FILE_NAME)

# duplicate the sort to post_flat_conversion_clim_sorts folder for comparisons later:
storage_name = OUTPUT_FILE_NAME[:-4] + EcoCyc_file[-13:]
post_flat_conversion_storage_path = os.path.join(
    raw_data_file_location, 'post_flat_conversion_protease_assignments', storage_name)

# Load the EcoCyc comparison file:
EcoCyc_file = pd.read_csv(ECOCYC_PATH, sep='\t', skiprows=[0])

# Load the input file:
protease_assignment_file = pd.read_excel(INPUT_PATH, )

# Function to find the monomer ID for each of the proteins in the input file:
def get_monomer_ID(EcoCyc_file, protease_assignment_file):
    # add a Monomer ID column to the input file
    protease_assignment_file['Monomer ID'] = None
    protease_assignment_file['Common Name'] = None
    for i in range(len(protease_assignment_file)):
        # find the monomer ID that matches the protein ID
        monomer_row = EcoCyc_file[EcoCyc_file['Gupta et al. 2024 Protein ID'] == protease_assignment_file['Protein ID'][i]]
        # add the monomer ID to the input file
        protease_assignment_file['Monomer ID'][i] = monomer_row['Monomer ID'].values[0]
        protease_assignment_file['Common Name'][i] = monomer_row['Common Name'].values[0]
    return protease_assignment_file

# Clean up the assignment column:
def clean_assignments(protease_assignment_file):
    # add a new protease assginment column to the input file
    protease_assignment_file['protease_assignment'] = None
    for i in range(len(protease_assignment_file)):
        # if the protease assignment is contains the word "only" in it,
        # then assign it to the current protease assignment column item
        if 'only' in protease_assignment_file['Protease assignment'][i]:
            protease_assignment_file['protease_assignment'][i] = protease_assignment_file['Protease assignment'][i]
        # if the assignment is not one singular protease, then list those included:
        else:
            proteases = protease_assignment_file['Protease assignment'][i] + ':'
            # check what proteases are included in the assignment:
            if protease_assignment_file['ClpP'][i] > 0:
                proteases += ' ClpP,'
                if protease_assignment_file['Lon'][i] > 0:
                    proteases += ' Lon'
                    if protease_assignment_file['HslV'][i] > 0:
                        proteases += ', HslV'
                else:
                    if protease_assignment_file['HslV'][i] > 0:
                        proteases += ' HslV'
            else:
                if protease_assignment_file['Lon'][i] > 0:
                    proteases += ' Lon'
                    if protease_assignment_file['HslV'][i] > 0:
                        proteases += ', HslV'

            # add the protease assignment to the protease assignment column item
            protease_assignment_file['protease_assignment'][i] = proteases

    return protease_assignment_file


# Convert the sort file to flat file format:
def make_flat_file(protease_assignment_file, output_path=OUTPUT_PATH):
    print("Generating flat file...")
    with io.open(output_path, 'wb') as f:
        print('Writing to {}'.format(f.name))
        writer = tsv.writer(f, quotechar="'", lineterminator='\n')
        writer.writerow(['# Generated by {} on {}'.format(__file__,
                                                          time.ctime())])
        writer.writerow(['"id"', '"common_name"', '"protease_assignment"', '"ClpP"', '"Lon"', '"HslV"', '"Unexplained"'])

        for i in range(len(protease_assignment_file)):
            writer.writerow([f'"{protease_assignment_file["Monomer ID"][i]}"',
                             f'"{protease_assignment_file["Common Name"][i]}"',
                             f'"{protease_assignment_file["protease_assignment"][i]}"',
                             f'{protease_assignment_file["ClpP"][i]}',
                             f'{protease_assignment_file["Lon"][i]}',
                             f'{protease_assignment_file["HslV"][i]}',
                             f'{protease_assignment_file["Unexplained"][i]}'])

# Get the monomer IDs for the proteins in the input file:
protease_assignment_file = get_monomer_ID(EcoCyc_file, protease_assignment_file)

# Clean up the assignment column:
protease_assignment_file = clean_assignments(protease_assignment_file)

# Make the flat file:
make_flat_file(protease_assignment_file)
make_flat_file(protease_assignment_file,
                output_path=post_flat_conversion_storage_path)

# Print the output file location:
print("Done. DON'T FORGET TO ADD THE FILE NAME TO knowledge_base_raw.py. "
      "File saved to: ", OUTPUT_PATH)
