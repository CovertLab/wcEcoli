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
get_current_EcoCyc_monomer_IDs_for_Clim_data.py (in the Gupta_et_al_Clim_data folder).
"""

import io
import os
import time
import pandas as pd
from wholecell.io import tsv
from wholecell.utils.filepath import ROOT_PATH

"""
USER INPUTS
"""

# PATH TO THE SORTED INPUT FILE (in the Clim_sort_files_pre_flat_conversion folder)
# Note: this file should have a column named "Protein ID" (formatted exactly as
# it is in the original file: 41467_2024_49920_MOESM4_ESM_ST1.xlsx)
# and a column named "half_life (units.min)"
protease_assignment_file = 'protease_assignments_test_file.xlsx'

# SPECIFY THE ECOCYC COMPARISON FILE (in the Clim_EcoCyc_monomer_ID_matches folder)
# Note: if a new EcoCyc update has been released, use
# get_current_EcoCyc_monomer_IDs_for_Clim_data.py to get an updated comparison file
EcoCyc_file = \
    '41467_2024_49920_MOESM4_ESM_ST1_EcoCyc_monomer_ID_matches_11202024.tsv'

# SPECIFY THE OUTPUT FILE NAME (to be placed in the flat file list)
# Note: this file name needs to be added to the list in "knowledge_base_raw.py"
# to be used in the model
OUTPUT_FILE_NAME = 'protease_assignments_Clim0_TEST.tsv' # e.g. 'protease_assignments_Clim#.tsv'

"""
END OF USER INPUTS
"""


# Current file location:
CURRENT_LOCATION = os.path.realpath(os.path.dirname(__file__))

# Path to the input file location:
INPUT_PATH = os.path.join(CURRENT_LOCATION, 'Gupta_et_al_Clim_data',
                     'protease_assignment_files_pre_flat_conversion',
                          protease_assignment_file)

# Path to the EcoCyc comparison file location:
ECOCYC_PATH = os.path.join(CURRENT_LOCATION, 'Gupta_et_al_Clim_data',
                           'Clim_EcoCyc_monomer_ID_matches', EcoCyc_file)

# Path to the output file location:
OUTPUT_PATH = os.path.join(ROOT_PATH, 'reconstruction', 'ecoli', 'flat',
                           OUTPUT_FILE_NAME)

# Load the EcoCyc comparison file:
EcoCyc_file = pd.read_csv(ECOCYC_PATH, sep='\t', skiprows=[0])

# Load the input file:
Clim_file = pd.read_excel(INPUT_PATH,)


# Function to find the monomer ID for each of the proteins in the input file:
def get_monomer_ID(EcoCyc_file, Clim_file):
    # add a Monomer ID column to the input file
    Clim_file['Monomer ID'] = None
    for i in range(len(Clim_file)):
        # find the monomer ID that matches the protein ID
        monomer_row = EcoCyc_file[EcoCyc_file['Gupta et al. 2024 Protein ID'] == Clim_file['Protein ID'][i]]
        # add the monomer ID to the input file
        Clim_file['Monomer ID'][i] = monomer_row['Monomer ID'].values[0]
    return Clim_file

# Clean up the assignment column:
def clean_assignments(Clim_file):
    # add a new protease assginment column to the input file
    Clim_file['protease_assignment'] = None
    for i in range(len(Clim_file)):
        # if the protease assignment is contains the word "only" in it,
        # then assign it to the current protease assignment column item
        if 'only' in Clim_file['Protease assignment'][i]:
            Clim_file['protease_assignment'][i] = Clim_file['Protease assignment'][i]
        # if the assignment is not one singular protease, then list those included:
        else:
            proteases = Clim_file['Protease assignment'][i] + ':'
            # check what proteases are included in the assignment:
            if Clim_file['ClpP'][i] > 0:
                proteases += ' ClpP,'
                if Clim_file['Lon'][i] > 0:
                    proteases += ' Lon'
                    if Clim_file['HslV'][i] > 0:
                        proteases += ', HslV'
                else:
                    if Clim_file['HslV'][i] > 0:
                        proteases += ' HslV'
            else:
                if Clim_file['Lon'][i] > 0:
                    proteases += ' Lon'
                    if Clim_file['HslV'][i] > 0:
                        proteases += ', HslV'

            # add the protease assignment to the protease assignment column item
            Clim_file['protease_assignment'][i] = proteases

    return Clim_file

# TODO: decide if the common name should also be included in the flat file

# Convert the sort file to flat file format:
def make_flat_file(Clim_file):
    print("Generating flat file...")
    with io.open(OUTPUT_PATH, 'wb') as f:
        print('Writing to {}'.format(f.name))
        writer = tsv.writer(f, quotechar="'", lineterminator='\n')
        writer.writerow(['# Generated by {} on {}'.format(__file__,
                                                          time.ctime())])
        writer.writerow(['id', 'protease_assignment'])

        for i in range(len(Clim_file)):
            writer.writerow([f'"{Clim_file["Monomer ID"][i]}"',
                             f'{Clim_file["protease_assignment"][i]}'])

# Get the monomer IDs for the proteins in the input file:
Clim_file = get_monomer_ID(EcoCyc_file, Clim_file)

# Clean up the assignment column:
Clim_file = clean_assignments(Clim_file)

# Make the flat file:
make_flat_file(Clim_file)

# Print the output file location:
print("Done. DON'T FORGET TO ADD THE FILE NAME TO knowledge_base_raw.py. "
      "File saved to: ", OUTPUT_PATH)
