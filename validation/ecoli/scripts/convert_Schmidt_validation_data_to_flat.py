"""
Converts sorted files of protein counts from Schmidt et al. 2016:
https://www.nature.com/articles/nbt.3418


NOTE:
After the new flat file is generated, it must also be added to the list in
validation_data_raw.py in order for it to be used in the model.

NOTE: This script uses the EcoCyc monomer IDs that correspond to the UniProt IDs
provided in the Schmidt et al. dataset. The EcoCyc monomer IDs may update
over time when the model is updated with a new version of EcoCyc, so it is
important to ensure that the EcoCyc monomer IDs are up-to-date from time to time
and generating new EcoCyc monomer ID match files as needed by using the script:
convert_Schmidt_UniProt_IDs_to_monomer_IDs.py (in the uniprot_data_conversion folder).

NOTE: you can plug in any old Ecocyc update comparison file into this to see
how the validation data changes over updates.

NOTE 2: "Common Name" here will refer to the common name of the that
corresponds to the Monomer ID listed in rnas.tsv. "Gene" will refer to the
gene name listed in the original Schmidt et al. dataset.
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

# PATH TO THE SORTED INPUT FILE
# (in validation/ecoli/flat/Schmidt_2016_uniprot_conversion_files/raw_files):
validation_file = 'Schmidt_2016_ST9_raw.csv'

# TABLE TYPE (6 or 9)
# specify either 6 or 9 depending on which table is being converted as the
# columns in the raw files differ. Other tables/columns can be added as well
# manually following the same format used in make_flat_file() below.
TABLE_TYPE = 9

# SPECIFY THE ECOCYC COMPARISON FILE
# (in validation/ecoli/flat/Schmidt_2016_uniprot_conversion_files/schmidt_uniprot_ids_to_monomer_ids)
# Note: if a new EcoCyc update has been released, use
# convert_Schmidt_UniProt_IDs_to_monomer_IDs.py to get an updated comparison file!
EcoCyc_file = \
    'Schmidt_2016_ST9_raw_EcoCyc_uniprot_ID_to_monomer_ID_matches_12052025.tsv'

# SPECIFY THE OUTPUT FILE NAME (to be placed in the flat file list)
# Note: this file name needs to be added to the list in "validation_data_raw.py"
# to be used in the model
OUTPUT_FILE_NAME = 'Schmidt_2016_ST9.tsv' # e.g. 'Schmidt_2016_ST6.tsv'


"""
END OF USER INPUTS
"""

# Current file location:
CURRENT_LOCATION = os.path.realpath(os.path.dirname(__file__))

# Path to the input file location:
INPUT_PATH = os.path.join(ROOT_PATH, 'validation', 'ecoli', 'flat',
                          'Schmidt_2016_uniprot_conversion_files', 'raw_files',
                          validation_file)

# Path to the EcoCyc comparison file location:
ECOCYC_PATH = os.path.join(ROOT_PATH, 'validation', 'ecoli', 'flat',
                          'Schmidt_2016_uniprot_conversion_files',
                           'schmidt_uniprot_ids_to_monomer_ids', EcoCyc_file)

# Path to the output file location:
OUTPUT_PATH = os.path.join(ROOT_PATH, 'validation', 'ecoli', 'flat',
                           OUTPUT_FILE_NAME)

# Also save a copy of the output file to the all_converted_schmidt_files folder:
relevant_date = EcoCyc_file[-13:-4]  # get the date from the EcoCyc file name
long_save_name = OUTPUT_FILE_NAME.replace('.tsv', f'{relevant_date}.tsv')
ALL_CONVERTED_SCHMIDT_PATH = os.path.join(ROOT_PATH, 'validation', 'ecoli',
                                          'flat', 'Schmidt_2016_uniprot_conversion_files',
                                          'all_converted_schmidt_files',
                                          long_save_name)

# Load the EcoCyc comparison file:
EcoCyc_file = pd.read_csv(ECOCYC_PATH, sep='\t', skiprows=[0])

# Load the input file:
validation_file = pd.read_csv(INPUT_PATH, header=0)

# Generate a dictionary mapping UniProt Accession to Monomer ID and another to Common Name:
UniProt_to_MonomerID = {}
UniProt_to_CommonName = {}
for i in range(len(EcoCyc_file)):
    UniProt_to_MonomerID[EcoCyc_file['Uniprot Accession'][i]] = EcoCyc_file['Monomer ID'][i]
    UniProt_to_CommonName[EcoCyc_file['Uniprot Accession'][i]] = EcoCyc_file['Common Name'][i]

# Function to find the monomer ID for each of the proteins in the input file:
def get_monomer_ID(EcoCyc_file, validation_file):
    # add a Monomer ID column to the input file
    validation_file['Monomer ID'] = None
    validation_file['Common Name'] = None
    for i in range(len(validation_file)):
        # find the monomer ID that matches the protein ID
        uniprot_id = EcoCyc_file['Uniprot Accession'][i]
        monomer_id = UniProt_to_MonomerID.get(uniprot_id)
        common_name = UniProt_to_CommonName.get(uniprot_id)
        # add the monomer ID to the input file
        validation_file['Monomer ID'][i] = monomer_id
        validation_file['Common Name'][i] = common_name
    return validation_file

# Convert the sort file to flat file format:
def make_flat_file(validation_file, table=6, output_path=OUTPUT_PATH):
    print("Generating flat file...")
    with io.open(output_path, 'wb') as f:
        print('Writing to {}'.format(f.name))
        writer = tsv.writer(f, quotechar="'", lineterminator='\n')
        writer.writerow(['# Generated by {} on {}'.format(__file__, time.ctime())])
        if table == 6:
            # relevant columns for table 6:
            writer.writerow(['Uniprot Accession','Description','Gene', 'Monomer ID',
                             'Common Name','Peptides.used.for.quantitation',
                             'Confidence.score','Molecular weight (Da)','Dataset',
                             'Glucose','LB','Glycerol + AA','Acetate','Fumarate',
                             'Glucosamine','Glycerol','Pyruvate','Chemostat µ=0.5',
                             'Chemostat µ=0.35','Chemostat µ=0.20','Chemostat µ=0.12',
                             'Stationary phase 1 day','Stationary phase 3 days',
                             'Osmotic-stress glucose','42°C glucose','pH6 glucose',
                             'Xylose','Mannose','Galactose','Succinate','Fructose',
                             'Glucose','LB','Glycerol + AA','Acetate','Fumarate'])

            for i in range(len(validation_file)):
                writer.writerow([f'"{validation_file["Uniprot Accession"][i]}"',
                                 f'"{validation_file["Description"][i]}"',
                                 f'"{validation_file["Gene"][i]}"',
                                 f'"{validation_file["Monomer ID"][i]}"',
                                 f'"{validation_file["Common Name"][i]}"',
                                 f'{validation_file["Peptides.used.for.quantitation"][i]}',
                                 f'{validation_file["Confidence.score"][i]}',
                                 f'{validation_file["Molecular weight (Da)"][i]}',
                                    f'"{validation_file["Dataset"][i]}"',
                                    f'{validation_file["Glucose"][i]}',
                                    f'{validation_file["LB"][i]}',
                                    f'{validation_file["Glycerol + AA"][i]}',
                                    f'{validation_file["Acetate"][i]}',
                                    f'{validation_file["Fumarate"][i]}',
                                    f'{validation_file["Glucosamine"][i]}',
                                    f'{validation_file["Glycerol"][i]}',
                                    f'{validation_file["Pyruvate"][i]}',
                                    f'{validation_file["Chemostat µ=0.5"][i]}',
                                    f'{validation_file["Chemostat µ=0.35"][i]}',
                                    f'{validation_file["Chemostat µ=0.20"][i]}',
                                    f'{validation_file["Chemostat µ=0.12"][i]}',
                                    f'{validation_file["Stationary phase 1 day"][i]}',
                                    f'{validation_file["Stationary phase 3 days"][i]}',
                                    f'{validation_file["Osmotic-stress glucose"][i]}',
                                    f'{validation_file["42°C glucose"][i]}',
                                    f'{validation_file["pH6 glucose"][i]}',
                                    f'{validation_file["Xylose"][i]}',
                                    f'{validation_file["Mannose"][i]}',
                                    f'{validation_file["Galactose "][i]}',
                                    f'{validation_file["Succinate"][i]}',
                                    f'{validation_file["Fructose"][i]}',
                                    f'{validation_file["Glucose.1"][i]}',
                                    f'{validation_file["LB.1"][i]}',
                                    f'{validation_file["Glycerol + AA.1"][i]}',
                                    f'{validation_file["Acetate.1"][i]}',
                                    f'{validation_file["Fumarate.1"][i]}'])
        elif table == 9:
            # relevant columns for table 9 (note: not all columns are included)
            writer.writerow(['Uniprot Accession', 'Description', 'Gene',
                             'Monomer ID', 'Common Name',
                             'Peptides.used.for.quantitation', 'Confidence.score',
                             'medianRatio_MG1655.LB_vs_BW25113.Glucose',
                             'medianRatio_NCM3722.LB_vs_BW25113.Glucose',
                             'medianRatio_MG1655.Glucose_vs_BW25113.Glucose',
                             'medianRatio_NCM3722.Glucose_vs_BW25113.Glucose',
                             'medianRatio_BW25113.LB_vs_BW25113.Glucose',
                             'Copies/Cell_MG1655.LB', 'Copies/Cell_MG1655.Glucose',
                             'Copies/Cell_NCM3722.LB', 'Copies/Cell_NCM3722.Glucose',
                             'Copies/Cell_BW25113.LB', 'Copies/Cell_BW25113.Glucose'])

            for i in range(len(validation_file)):
                writer.writerow([f'"{validation_file["Uniprot Accession"][i]}"',
                                 f'"{validation_file["Description"][i]}"',
                                 f'"{validation_file["Gene"][i]}"',
                                 f'"{validation_file["Monomer ID"][i]}"',
                                 f'"{validation_file["Common Name"][i]}"',
                                 f'{validation_file["Peptides.used.for.quantitation"][i]}',
                                 f'{validation_file["Confidence.score"][i]}',
                                 f'{validation_file["medianRatio_MG1655.LB_vs_BW25113.Glucose"][i]}',
                                 f'{validation_file["medianRatio_NCM3722.LB_vs_BW25113.Glucose"][i]}',
                                 f'{validation_file["medianRatio_MG1655.Glucose_vs_BW25113.Glucose"][i]}',
                                 f'{validation_file["medianRatio_NCM3722.Glucose_vs_BW25113.Glucose"][i]}',
                                 f'{validation_file["medianRatio_BW25113.LB_vs_BW25113.Glucose"][i]}',
                                 f'{validation_file["Copies/Cell_MG1655.LB"][i]}',
                                 f'{validation_file["Copies/Cell_MG1655.Glucose"][i]}',
                                 f'{validation_file["Copies/Cell_NCM3722.LB"][i]}',
                                 f'{validation_file["Copies/Cell_NCM3722.Glucose"][i]}',
                                 f'{validation_file["Copies/Cell_BW25113.LB"][i]}',
                                 f'{validation_file["Copies/Cell_BW25113.Glucose"][i]}'])

        else:
            raise ValueError("TABLE_TYPE must be either 6 or 9, or add new "
                             "table format manually.")


# Get the monomer IDs for the proteins in the input file:
validation_file = get_monomer_ID(EcoCyc_file, validation_file)

# Make the flat file:
make_flat_file(validation_file, table=TABLE_TYPE)

# Also save a copy of the output file to the all_converted_schmidt_files folder:
make_flat_file(validation_file, table=TABLE_TYPE, output_path=ALL_CONVERTED_SCHMIDT_PATH)

# Print the output file location:
print("Done. DON'T FORGET TO ADD THE FILE NAME TO validation_data_raw.py!!! "
      "File saved to: ", OUTPUT_PATH)
