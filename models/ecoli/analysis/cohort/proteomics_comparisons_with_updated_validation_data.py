# Plan for file:
# 1. load in the three raw validation proteomics datasets:
#    - Schmidt BW 2016 (ST6) [named SBWST6]
#    - Schmidt MG 2016 (ST9) [named SMGST9]
#    - Schmidt BW 2016 (ST9) [named SBWST9]
# 2. Load in the simulated proteomics data from a simulation AND validation data
#    - Sim Data from CLMLNE1 [named CLMLNE1]
#    - Validation BW Data from Schmidt 2015 [named VBW2015]
# 3. Using all the 4309 protein names in the simulated data, make a map of protein ids to gene ids
#    - gene_id_to_protein_id_dict
#   NOTE: consider looking at the function where every common name is searched? Not just the first one in RNAS.tsv
# 4. Then, make a map of gene ids to descriptive IDs, half life, and counts

# 5. Make a plot comparing the four validation data sets to the simulated data
# 6. Make plots comparing each validation data set to each other
# TODO: update description above to explain what this file does. explain where this pulls from and how to update the flat files based on the ECOCYC udpates

import pickle
import os
import ast
import matplotlib.pyplot as plt
import csv
import pandas as pd
import numpy as np
from models.ecoli.analysis import cohortAnalysisPlot
from wholecell.analysis.analysis_tools import (exportFigure,
    read_bulk_molecule_counts, read_stacked_bulk_molecules,read_stacked_columns)
from wholecell.io.tablereader import TableReader
from wholecell.utils.protein_counts import get_simulated_validation_counts
import plotly.graph_objects as go
from sklearn.metrics import r2_score
from scipy.stats import pearsonr

""" USER INPUTS """

# Indicate the number of generations to be ignored at the start of each seed:
IGNORE_FIRST_N_GENS = 2 # 2 for local, 14 for Sherlock (w/ 24 total gens)

# Indicate any proteins that you want to highlight on the plots (these can include the compartment tag or not, it does not matter):
PROTEINS_OF_INTEREST = ["EG10241-MONOMER", "1-ACYLGLYCEROL-3-P-ACYLTRANSFER-MONOMER", "EG10599-MONOMER", "EG10241-MONOMER", "EG10762-MONOMER", "G7320-MONOMER"]

# Indicate whether to highlight the proteins listed in PROTEINS_OF_INTEREST on the plots:
HIGHLIGHT_PROTEINS_OF_INTEREST = False

# Indicate whether to generate a plot with some manually mapped proteins included
# (not all Uniprot IDs successfully match to protein IDs via the BioCyc mapping
# completed in validation/ecoli/scripts/convert_Schmidt_UniProt_IDs_to_monomer_IDs.py,
# so the proteins that could not be mapped via their uniprot ID are manually
# mapped using their common name and finding the match within rnas.tsv.
# This option was added because it appeared that at least a few of these
# unmapped Uniprot IDs did match to proteins with common names that matched the
# one listed on the raw validation file when manually searched in EcoCyc.
# However, it should be noted that not all unmapped proteins have been checked
# for this, so keep that in mind when analyzing results.)
INCLUDE_MANUALLY_MAPPED_PROTEINS = False

# Indicate if comparisons of validation source to validation source should be plotted:
PLOT_VALIDATION_SOURCE_COMPARISONS = True

  # example Uniprot IDs
# TODO: check If you can read in uniprot from ecyocyc and not biocyc
# TODO: add option to output validation comparison plots as well
# TODO: add "Glucose" to the long name of the data
""" END USER INPUTS """


# PRINT NOTE ABOUT UPDATING FLAT FILES
print("NOTE: Check for recent EcoCyc updates for the most up-to-date plots. "
      "Please update the flat files in "
      "models/ecoli/data/ecoli/flat/ accordingly using the "
      "convert_Schmidt_UniProt_IDs_to_monomer_IDs.py script and "
      "convert_Schmidt_validation_data_to_flat.py scripts in "
      "validation/ecoli/scripts if needed.")


class Plot(cohortAnalysisPlot.CohortAnalysisPlot):
    def generate_data(self, simDataFile):
        """
        Extracts the average total monomer counts for each protein in the
        simulation
        Args:
            simDataFile: simulation data file

        Returns:
            total_protein_counts: average total monomer counts for all
            proteins in the simulation
            self.all_monomer_ids: list of all the monomer ids in the simulation
        """
        with open(simDataFile, 'rb') as f:
            sim_data = pickle.load(f)
        monomer_sim_data = (
            sim_data.process.translation.monomer_data.struct_array)

        # Extract monomer count data for each protein:
        self.all_monomer_ids = monomer_sim_data['id']
        all_cells = self.ap.get_cells(
            generation=np.arange(IGNORE_FIRST_N_GENS,self.n_total_gens),
            only_successful=True)

        self.total_cells = len(all_cells)

        # Get the average total protein counts for each monomer:
        total_counts = (
            read_stacked_columns(all_cells, 'MonomerCounts',
                                 'monomerCounts', ignore_exception=True))
        avg_total_counts = np.mean(total_counts, axis=0)

        # Get the average free protein counts for each monomer:
        (free_counts,) = read_stacked_bulk_molecules(
            all_cells, self.all_monomer_ids, ignore_exception=True)
        avg_free_counts = np.mean(free_counts, axis=0)

        # Get the average complex counts for each monomer:
        avg_counts_for_monomers_in_complexes = avg_total_counts - avg_free_counts

        return avg_total_counts, avg_free_counts, avg_counts_for_monomers_in_complexes

    def get_validation_data(self, simDataFile, validationDataFile):
        # adapted from multigen and single/proteinCountsValidation.py
        """
        Extract the protein counts for proteins that are present in either
        validation data set
        Args:
            simDataFile: Simulation data file
            validationDataFile: validation data file

        Returns: lists of proteins that are present in both the simulation and
        validation datasets, as well as the protein counts for each protein.
        """
        # Get simulation data:
        sim_data = self.read_pickle_file(simDataFile)
        sim_monomer_ids = sim_data.process.translation.monomer_data["id"]

        # Get validation data:
        validation_data = self.read_pickle_file(validationDataFile)
        wisniewski_ids = validation_data.protein.wisniewski2014Data["monomerId"]
        schmidt_ids = validation_data.protein.schmidt2015Data["monomerId"]
        wisniewski_counts = validation_data.protein.wisniewski2014Data["avgCounts"]
        schmidt_counts = validation_data.protein.schmidt2015Data["glucoseCounts"]

        # Get the simulation cell directories:
        all_cells = self.ap.get_cells(
            generation=np.arange(IGNORE_FIRST_N_GENS, self.n_total_gens),
            only_successful=True)

        # get all the simulation monomer counts:
        monomer_counts = (read_stacked_columns(all_cells,
                                               'MonomerCounts', 'monomerCounts'))

        # the simluation and validation datasets:
        sim_schmidt_counts, val_schmidt_counts, schmidt_overlap_ids = (
            get_simulated_validation_counts(
                schmidt_counts, monomer_counts, schmidt_ids, sim_monomer_ids))
        sim_wisniewski_counts, val_wisniewski_counts, wisniewski_overlap_ids = (
            get_simulated_validation_counts(wisniewski_counts,
                            monomer_counts, wisniewski_ids, sim_monomer_ids))

        return (sim_schmidt_counts, val_schmidt_counts,
                schmidt_overlap_ids, sim_wisniewski_counts,
                val_wisniewski_counts, wisniewski_overlap_ids)

    # Create a function preps the validation data saved with the sim data for the comparison:
    def prep_validation_data_for_comparison(self, simDataFile, validationDataFile):
        # Obtain the validation data saved with the sim data:
        (sim_schmidt_counts, val_schmidt_counts,
         schmidt_overlap_ids, sim_wisniewski_counts,
         val_wisniewski_counts,
         wisniewski_overlap_ids) = self.get_validation_data(
            simDataFile, validationDataFile)

        # Generate a validation data dictionary for easy access later:
        sim_validation_common_names_to_counts_dict = {}
        sim_validation_common_names_to_monomer_ids_dict = {}
        for i in range(len(schmidt_overlap_ids)):
            monomer_id = schmidt_overlap_ids[i][:-3]  # remove compartment tag
            common_name = self.sim_monomer_ids_to_common_names_dict[monomer_id]
            val_count = val_schmidt_counts[i]
            sim_validation_common_names_to_counts_dict[common_name] = val_count
            sim_validation_common_names_to_monomer_ids_dict[common_name] = monomer_id

        # TODO: consider also making it possible to compare Wisniewski if there appears to be a desire to do so.

        return (sim_validation_common_names_to_counts_dict,
                sim_validation_common_names_to_monomer_ids_dict)



    # Function that checks if an input is a vaild molecule in the simulation:
    def check_validity_and_get_compartment(self, simDataFile, molecule_list):
        with open(simDataFile, 'rb') as f:
            sim_data = pickle.load(f)
        revised_molecule_list = []
        for protein in molecule_list:
            if "[" in protein:
                protein = protein[:-3]  # remove compartment
            if sim_data.getter.is_valid_molecule(protein):
                revised_name = protein + sim_data.getter.get_compartment_tag(protein)
                revised_molecule_list.append(revised_name)

        return revised_molecule_list


    # Extract the data from the Schmidt et al. ST6 (BW25113) dataset:
    def get_schmidt_BW_ST6_data(self,):
        # NOTE: here, monomer IDs refer to the monomer ID mapped to by the
        # Uniprot IDs in the Schmidt et al. dataset and the corresponding
        # rnas.tsv common name for that monomer ID

        # Retrieve the common name from rnas.tsv for each monomer ID
        schmidt_ST6 = pd.read_csv(
            '~/wcEcoli/validation/ecoli/flat/Schmidt_2016_ST6.tsv',
            sep='\t', header=1)

        # Make a dictionary mapping the Uniprot Accession to Gene, Monomer ID,
        # and Glucose Counts:
        uniprot_IDs_to_schmidt_common_names = {}
        uniprot_IDs_to_monomer_IDs = {} # for now, do not map directly to monomer IDs as some are None
        uniprot_IDs_to_schmidt_glucose_counts = {}
        for i in range(len(schmidt_ST6)):
            uniprot_ID = schmidt_ST6.iloc[i]['Uniprot Accession']
            gene_symbol = schmidt_ST6.iloc[i]['Gene']
            monomer_ID = schmidt_ST6.iloc[i]['Monomer ID']
            glucose_counts = schmidt_ST6.iloc[i]['Glucose']

            #uniprot_IDs_to_schmidt_common_names[uniprot_ID] = gene_symbol
            #uniprot_IDs_to_monomer_IDs[uniprot_ID] = monomer_ID
            #uniprot_IDs_to_schmidt_glucose_counts[uniprot_ID] = glucose_counts

            if monomer_ID == None:
                print("WARNING: No simulation monomer ID mapped for Uniprot ID",
                      uniprot_ID, " (common name: ", gene_symbol,
                      ") within the Schmidt et al. 2016 dataset.")

            uniprot_IDs_to_schmidt_common_names[uniprot_ID] = gene_symbol
            uniprot_IDs_to_monomer_IDs[uniprot_ID] = monomer_ID
            uniprot_IDs_to_schmidt_glucose_counts[uniprot_ID] = glucose_counts

        return (uniprot_IDs_to_schmidt_common_names,
                uniprot_IDs_to_monomer_IDs,
                uniprot_IDs_to_schmidt_glucose_counts)


    # Obtain a dictionary of the common names associated with multiple monomer IDs:
    def check_common_name_to_gene_id_map(self, simDataFile):
        # Retrieve the common name from rnas.tsv for each monomer ID
        rnasTable = pd.read_csv(
            '~/wcEcoli/reconstruction/ecoli/flat/rnas.tsv',
            sep='\t', skiprows=[0, 1, 2, 3], header=1)
        rnas_common_names = rnasTable[['monomer_ids', 'common_name']]

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
        common_names_with_multiple_monomer_ids_dict = {}
        for i in range(len(rnas_data)):
            if len(rnas_data.iloc[i]['monomer_ids']) > 1:
                multiple_monomer_ids.append(rnas_data.iloc[i]['monomer_ids'])
                common_names_with_multiple_monomer_ids.append(
                    rnas_data.iloc[i]['common_name'])
                common_names_with_multiple_monomer_ids_dict[
                    rnas_data.iloc[i]['common_name']] = rnas_data.iloc[i]['monomer_ids']

        # View which common names correspond to multiple monomer IDs:
        multiple_monomer_ids_df = pd.DataFrame(
            {'common_name': common_names_with_multiple_monomer_ids,
             'monomer_ids': multiple_monomer_ids})
        print("Gene common names with multiple monomer IDs in rnas.tsv:",
              len(multiple_monomer_ids_df))
        print(multiple_monomer_ids_df)

        # Check which IDs end up in the simulation...
        print(
            "Checking which, if any, of the monomer IDs belonging to the same gene common name are in the simulation:")
        common_names_with_multiple_monomer_ids_narrowed_dict = {}
        for common_name in common_names_with_multiple_monomer_ids_dict.keys():
            monomer_ids = common_names_with_multiple_monomer_ids_dict[common_name]
            monomer_ids_in_sim = []
            print("Checking if any monomer IDs for", common_name, "are in the simulation...")
            for monomer_id in monomer_ids:
                monomer_ID = self.check_validity_and_get_compartment(simDataFile, [monomer_id])[0]
                if monomer_ID in self.all_monomer_ids:
                    monomer_ids_in_sim.append(
                        monomer_ID[:-3])  # remove compartment tag for clarity
                    print(monomer_ID, "is in the simulation.")
            # only append the monomer IDs that are actually in the simulation:
            common_names_with_multiple_monomer_ids_narrowed_dict[common_name] = monomer_ids_in_sim

        # make it official:
        self.common_names_with_multiple_monomer_ids_narrowed_dict = common_names_with_multiple_monomer_ids_narrowed_dict

        # return a dictionary of these values in case it is helpful later:
        return common_names_with_multiple_monomer_ids_dict

    # Create a dictionary mapping all monomer IDs in rnas.tsv to common names (and vice versa):
    def map_monomer_ids_to_common_names(self):
        # Retrieve the common name from rnas.tsv for each monomer ID
        rnasTable = pd.read_csv(
            '~/wcEcoli/reconstruction/ecoli/flat/rnas.tsv', sep='\t',
            skiprows=[0, 1, 2, 3], header=1)
        rnas_common_names = rnasTable[['monomer_ids', 'common_name']]

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

        # Create a dictionary mapping common names to monomer IDs and vice versa:
        print("Mapping simulation monomer IDs to common names...")
        common_name_to_monomer_id = {}
        monomer_id_to_common_name = {}
        for i in range(len(rnas_data)):
            common_name = rnas_data.iloc[i]['common_name']
            monomer_ids = rnas_data.iloc[i]['monomer_ids']

            # add the compartment tag to each monomer ID:
            # just kidding, do not do this, takes way too much time to do here.

            # if only one monomer ID maps to the common name, map it directly:
            if len(monomer_ids) == 1:
                common_name_to_monomer_id[common_name] = monomer_ids[0]
                monomer_id_to_common_name[monomer_ids[0]] = common_name

            # generate a warning if a common name maps to multiple monomer IDs:
            if len(monomer_ids) > 1:
                print("WARNING: ", common_name,
                      "maps to multiple monomer IDs:", monomer_ids)
                # find the monomer ID that is actually in the simulation:
                monomer_IDs_in_sim = self.common_names_with_multiple_monomer_ids_narrowed_dict[
                    common_name]
                if len(monomer_IDs_in_sim) == 1:
                    print("Mapping to the "
                          "monomer ID that appears in the simulation by default:",
                          monomer_IDs_in_sim)
                    # map to first monomer ID by default
                    common_name_to_monomer_id[common_name] = monomer_IDs_in_sim[0]
                    monomer_id_to_common_name[monomer_IDs_in_sim[0]] = common_name
                if len(monomer_IDs_in_sim) > 1:
                    print("Multiple monomer IDs for", common_name,
                          "are in the simulation:", monomer_IDs_in_sim,
                          "Mapping to the first one by default:", monomer_IDs_in_sim[0])
                    common_name_to_monomer_id[common_name] = [monomer_IDs_in_sim[0]]
                    monomer_id_to_common_name[monomer_IDs_in_sim[0]] = common_name

        return monomer_id_to_common_name, common_name_to_monomer_id

    # Create a function that can get the simuation descriptions from proteins.tsv:
    def get_monomer_descriptions_from_flat_file(self):
        # Retrieve the common name from rnas.tsv for each monomer ID
        proteinsTable = pd.read_csv(
            '~/wcEcoli/reconstruction/ecoli/flat/proteins.tsv', sep='\t',
            skiprows=[0, 1, 2, 3], header=1)
        proteins_data = proteinsTable[['id', 'common_name']]

        # Create a dictionary mapping monomer ids to descriptions:
        monomer_id_to_common_name_description = {}
        for i in range(len(proteins_data)):
            monomer_id = proteins_data.iloc[i]['id']
            common_name = proteins_data.iloc[i]['common_name']
            monomer_id_to_common_name_description[monomer_id] = common_name

        return monomer_id_to_common_name_description

    # Create a function that maps monomer IDs in the rnas.tsv file to their synonyms:
    def map_synonyms_to_monomer_ids(self):
        # Retrieve the synonyms from rnas.tsv for each monomer ID
        rnasTable = pd.read_csv(
            '~/wcEcoli/reconstruction/ecoli/flat/rnas.tsv', sep='\t',
            skiprows=[0, 1, 2, 3], header=1)
        rnas_synonyms = rnasTable[['monomer_ids', 'synonyms']]

        # Function to evaluate strings as lists, as the format of the synonyms and monomer_ids columns in rnas.tsv is a string
        def parse_rows(row):
            try:
                # Evaluate the string as a list
                return ast.literal_eval(row)
            except (ValueError, SyntaxError):
                return []

        # Apply the parsing function to the monomer_ids column to make it a real list
        rnas_data = rnas_synonyms.copy(deep=True)
        rnas_data['monomer_ids'] = rnas_data['monomer_ids'].apply(parse_rows)
        rnas_data['synonyms'] = rnas_data['synonyms'].apply(parse_rows)
        # TODO: might need to delete synonyms that match the gene_symbol itself if this occurs

        # Create a dictionary mapping monomer IDs to synonyms:
        synonyms_to_monomer_ids = {}
        monomer_ids_to_synonyms = {}
        for i in range(len(rnas_data)):
            monomer_ids = rnas_data.iloc[i]['monomer_ids']
            synonyms = rnas_data.iloc[i]['synonyms']
            for synonym in synonyms:
                synonyms_to_monomer_ids[synonym] = monomer_ids
            for monomer_id in monomer_ids:
                monomer_ids_to_synonyms[monomer_id] = synonyms

        return synonyms_to_monomer_ids, monomer_ids_to_synonyms

    # Create a function that checks if there is a synonym for a given common name:
    def check_for_synonym_match_in_sim_monomer_ids(self, simDataFile, common_name):
        synonyms_to_monomer_ids, monomer_ids_to_synonyms = self.map_synonyms_to_monomer_ids()
        monomer_IDs_in_sim = []
        if common_name in synonyms_to_monomer_ids.keys():
            monomer_ids = synonyms_to_monomer_ids[common_name]
            print("The common name", common_name,
                  "has synonyms that map to the following monomer IDs:",
                  monomer_ids)
            # check which (if any) of these monomer IDs are in the simulation:
            print("Checking if any of these monomer IDs are in the simulation...")
            for monomer_id in monomer_ids:
                monomer_ID = self.check_validity_and_get_compartment(
                    simDataFile, [monomer_id])[0]
                if monomer_ID in self.all_monomer_ids:
                    monomer_IDs_in_sim.append(monomer_ID[:-3])
                    print(monomer_ID, "is in the simulation.")

        return monomer_IDs_in_sim

    # Create a dictionary mapping the simulation monomer IDs to common names, descriptions, half lifes, average protein counts, and monomer names with compartment tags:
    def generate_sim_monomer_id_to_info_dicts(self, simDataFile, avg_total_counts, avg_free_counts, avg_counts_for_monomers_in_complexes):
        # to reduce the run time associated with using the self.check_validity_and_get_compartment() function,
        # instead just map all the monomer ids to their compartment tags here:
        monomer_id_to_monomer_id_with_compartment_tag_dict = {}
        for monomer_id in self.all_monomer_ids:
            base_monomer_id = monomer_id[:-3]
            monomer_id_to_monomer_id_with_compartment_tag_dict[base_monomer_id] = monomer_id

        # map simulation monomer IDs to common names (and vice versa):
        sim_monomer_ids_to_common_names_dict = {}
        sim_common_names_to_monomer_ids_dict = {}
        for monomer_id in self.all_monomer_ids:
            base_monomer_id = monomer_id[:-3]
            common_name = self.monomer_id_to_common_name[base_monomer_id]
            sim_monomer_ids_to_common_names_dict[base_monomer_id] = common_name
            sim_common_names_to_monomer_ids_dict[common_name] = base_monomer_id

        # map simulation monomer IDs to descriptions:
        sim_monomer_ids_to_descriptions_dict = {}
        # retrieve the monomer ID to common name description mapping from proteins.tsv:
        monomer_id_to_common_name_description = self.get_monomer_descriptions_from_flat_file()
        for monomer_id in self.all_monomer_ids:
            base_monomer_id = monomer_id[:-3]
            description = monomer_id_to_common_name_description[base_monomer_id]
            sim_monomer_ids_to_descriptions_dict[base_monomer_id] = description

        # map simulation monomer IDs to half lifes:
        with open(simDataFile, 'rb') as f:
            sim_data = pickle.load(f)
        deg_rates = sim_data.process.translation.monomer_data["deg_rate"]
        # convert deg_rate to half life in minutes:
        half_lives = np.log(2) / deg_rates.asNumber() / 60  # convert to minutes
        monomer_ids_wo_compartment = [id[:-3] for id in self.all_monomer_ids]
        sim_monomer_ids_to_half_lifes_dict = dict(zip(monomer_ids_wo_compartment, half_lives))

        # map simulation monomer IDs to average total protein counts:
        sim_monomer_ids_to_avg_total_protein_counts_dict = dict(zip(monomer_ids_wo_compartment, avg_total_counts))
        # map simulation monomer IDs to average free protein counts:
        sim_monomer_ids_to_avg_free_protein_counts_dict = dict(zip(monomer_ids_wo_compartment, avg_free_counts))
        # map simulation monomer IDs to average complex protein counts:
        sim_monomer_ids_to_avg_complex_protein_counts_dict = dict(zip(monomer_ids_wo_compartment, avg_counts_for_monomers_in_complexes))

        # make official:
        self.sim_monomer_ids_to_monomer_ids_with_compartment_tags_dict = monomer_id_to_monomer_id_with_compartment_tag_dict
        self.sim_monomer_ids_to_common_names_dict = sim_monomer_ids_to_common_names_dict
        self.sim_common_names_to_monomer_ids_dict = sim_common_names_to_monomer_ids_dict
        self.sim_monomer_ids_to_descriptions_dict = sim_monomer_ids_to_descriptions_dict
        self.sim_monomer_ids_to_half_lifes_dict = sim_monomer_ids_to_half_lifes_dict
        self.sim_monomer_ids_to_avg_total_protein_counts_dict = sim_monomer_ids_to_avg_total_protein_counts_dict
        self.sim_monomer_ids_to_avg_free_protein_counts_dict = sim_monomer_ids_to_avg_free_protein_counts_dict
        self.sim_monomer_ids_to_avg_complex_protein_counts_dict = sim_monomer_ids_to_avg_complex_protein_counts_dict


    # Function that matches protein IDs from a given validation dataset to simulation monomer IDs:
    def match_validation_dataset_monomer_IDs_to_simulation_monomer_IDs(
            self, uniprot_IDs_to_monomer_IDs, uniprot_IDs_to_schmidt_common_names):

        # Determine which Uniprot IDs correspond to an actual monomer ID in the simulation:
        validation_dataset_uniprot_ids_to_sim_monomer_ids_dict = {}
        for uniprot_ID in uniprot_IDs_to_monomer_IDs.keys():
            monomer_ID = uniprot_IDs_to_monomer_IDs[uniprot_ID]
            if monomer_ID is not None:
                # check if this monomer ID is valid:
                if monomer_ID in self.sim_monomer_ids_to_monomer_ids_with_compartment_tags_dict.keys():
                    validation_dataset_uniprot_ids_to_sim_monomer_ids_dict[uniprot_ID] = monomer_ID
            else:
                print("WARNING: No simulation monomer ID mapped to UniProt ID:",
                      uniprot_ID, " (listed as", uniprot_IDs_to_schmidt_common_names[uniprot_ID], "in Schmidt et al.)")


        return validation_dataset_uniprot_ids_to_sim_monomer_ids_dict

    # Function that finds monomer ids for uniprot ids that were not mapped directly:
    def find_monomer_ids_for_unmapped_uniprot_ids(self, uniprot_IDs_to_RVS_common_names,
                uniprot_IDs_to_monomer_IDs):
        # find the uniprot IDs that mapped to None monomer IDs:
        unmapped_uniprot_IDs = []
        for uniprot_ID in uniprot_IDs_to_monomer_IDs.keys():
            monomer_ID = uniprot_IDs_to_monomer_IDs[uniprot_ID]
            if pd.isna(monomer_ID):
                unmapped_uniprot_IDs.append(uniprot_ID)


        # Find the RVS common names for these unmapped uniprot IDs:
        unmapped_uniprot_IDs_to_RVS_common_names = {}
        for uniprot_ID in unmapped_uniprot_IDs:
            RVS_common_name = uniprot_IDs_to_RVS_common_names[uniprot_ID]
            unmapped_uniprot_IDs_to_RVS_common_names[uniprot_ID] = RVS_common_name

        # Use these RVS common names to try to map to simulation monomer IDs:
        monomer_id_to_common_name, common_name_to_monomer_id  = self.map_monomer_ids_to_common_names()

        # define a dictionary mapping the unmapped uniprot IDs to simulation monomer IDs:
        unmapped_uniprot_IDs_to_sim_monomer_IDs = {}
        unmapped_uniprot_IDs_to_sim_common_names = {}
        for uniprot_ID in unmapped_uniprot_IDs_to_RVS_common_names.keys():
            RVS_common_name = unmapped_uniprot_IDs_to_RVS_common_names[uniprot_ID]
            if RVS_common_name in common_name_to_monomer_id.keys():
                sim_monomer_ID = common_name_to_monomer_id[RVS_common_name]
                unmapped_uniprot_IDs_to_sim_monomer_IDs[uniprot_ID] = sim_monomer_ID
                unmapped_uniprot_IDs_to_sim_common_names[uniprot_ID] = RVS_common_name
                print("The unmapped UniProt ID", uniprot_ID,
                      "maps to simulation monomer ID", sim_monomer_ID,
                      "via the raw validation source's common name: ", RVS_common_name)
            else:
                print("The unmapped UniProt ID", uniprot_ID,
                      "with RVS common name", RVS_common_name,
                      "did not map to any simulation monomer ID via common names in rnas.tsv.")


        return unmapped_uniprot_IDs_to_sim_monomer_IDs, unmapped_uniprot_IDs_to_sim_common_names




    # Create a function that automatically maps a list of common names to simulation monomer IDs:
    def map_common_names_to_sim_monomer_ids(self, simDataFile, common_name_list, record_extra_mappings=False):
        # Need to be careful here becuase there is a small chance that a
        # common name in the dataset being compared here uses a common name for
        # a protein that is not listed as the default "common_name" in rnas.tsv
        # and instead is in the "synonyms" column of rnas.tsv

        # make this function run to get self.map_synonyms_to_monomer_ids defined:
        synonyms_to_monomer_ids, monomer_ids_to_synomnyms = self.map_synonyms_to_monomer_ids()

        # define a dictionary mapping the comparison dataset's common names to the simulation monomer IDs:
        comparison_dataset_common_name_to_sim_monomer_id_dict = {}

        # define a list of common names not directly found in the sim data common names:
        comparison_dataset_common_names_not_mapped_directly_to_sim_common_names_list = []

        # define a dictionary that shows what the mapping was for dataset common names that mapped
        # to a monomer ID via an rnas.tsv synonyms match:
        comparison_dataset_common_name_to_simulation_common_name_dict = {}

        # define a dictionary that records which common names could not be mapped to simulation monomer IDs:
        comparison_dataset_common_names_not_mapped_to_sim_ids_list = []

        print("Mapping common names from the validaiton comparison dataset to simulation monomer IDs...")
        for common_name in common_name_list:
            if common_name in self.sim_common_names_to_monomer_ids_dict.keys():
                monomer_id = self.sim_common_names_to_monomer_ids_dict[common_name]
                comparison_dataset_common_name_to_sim_monomer_id_dict[common_name] = monomer_id

            else:
                comparison_dataset_common_names_not_mapped_directly_to_sim_common_names_list.append(common_name)

        initial_mapping_count = len(comparison_dataset_common_name_to_sim_monomer_id_dict.keys())
        print(initial_mapping_count, "common names from the comparison dataset were mapped directly to simulation monomer IDs based on matching common names.")
        print(f"Now checking if any of the remaining {int(len(common_name_list)) - int(initial_mapping_count)} common names from the comparison dataset map to simulation monomer IDs via synonyms in rnas.tsv...")

        # For those common names that were not direct match with simulation
        # common name ids, check if they match with a synonym in rnas.tsv that maps
        # to a simulation monomer id (note: split up the for loops this way
        # because some "common names" in rnas.tsv are also listed in "synonyms"):
        for common_name in comparison_dataset_common_names_not_mapped_directly_to_sim_common_names_list:
            # if it is not in the sim_common_names_to_monomer_ids_dict, check if it is anywhere in the common_name column of the rnas.tsv file:
            if common_name in self.common_name_to_monomer_id:
                monomer_id = self.common_name_to_monomer_id[common_name]
                print("The monomer", monomer_id, "(common name:",
                    common_name, ") from the comparison dataset is in rnas.tsv "
                    "but is not a protein that matches to a simulation common "
                                 "name mapping, so it won't be included in the mapping.")
                comparison_dataset_common_names_not_mapped_to_sim_ids_list.append(common_name)

            # if it is not in the rnas.tsv "common name" column, check if it is in the synonyms column of rnas.tsv:
            elif common_name in synonyms_to_monomer_ids.keys():
                matching_proteins = self.check_for_synonym_match_in_sim_monomer_ids(simDataFile, common_name)
                if len(matching_proteins) == 1:
                    monomer_id = matching_proteins[0]
                    comparison_dataset_common_name_to_sim_monomer_id_dict[common_name] = monomer_id
                    corresponding_sim_common_name = self.sim_monomer_ids_to_common_names_dict[monomer_id]
                    comparison_dataset_common_name_to_simulation_common_name_dict[common_name] = corresponding_sim_common_name
                    print("The common name", common_name,
                          "from the validation comparison dataset maps to the "
                          "simulation monomer ID", monomer_id, "via synonyms in rnas.tsv, so it will be included in the mapping.")
                if len(matching_proteins) > 1:
                    print("The common name", common_name,
                          "from the comparison dataset maps to multiple "
                          "simulation monomer IDs (",matching_proteins, ") via synonyms in rnas.tsv. Matching the common name to the first monomer ID by default:",
                          matching_proteins[0])
                    monomer_id = matching_proteins[0]
                    comparison_dataset_common_name_to_sim_monomer_id_dict[common_name] = monomer_id
                    corresponding_sim_common_name = self.sim_monomer_ids_to_common_names_dict[
                        monomer_id]
                    comparison_dataset_common_name_to_simulation_common_name_dict[
                        common_name] = corresponding_sim_common_name

            else:
                comparison_dataset_common_names_not_mapped_to_sim_ids_list.append(common_name)
                print(common_name, "did not map to a simulation monomer name via the simulation common names or rnas.tsv synonyms, and thus, will not be included in the mapping.")

        # Share final stats:
        final_mapping_count = len(comparison_dataset_common_name_to_sim_monomer_id_dict)
        mapped_via_synonyms_count = len(comparison_dataset_common_name_to_simulation_common_name_dict.keys())
        print(mapped_via_synonyms_count, f"common names from the comparison dataset were mapped to simulation monomer IDs after checking synonyms in rnas.tsv ({final_mapping_count} mapped total).")
        print("A total of", len(comparison_dataset_common_names_not_mapped_to_sim_ids_list),
              "common names from the comparison dataset could not be mapped to simulation monomer IDs.")

        if record_extra_mappings:
            return (comparison_dataset_common_name_to_sim_monomer_id_dict,
                    comparison_dataset_common_name_to_simulation_common_name_dict,
                    comparison_dataset_common_names_not_mapped_to_sim_ids_list)
        else:
            return comparison_dataset_common_name_to_sim_monomer_id_dict

    # Create a function for mapping input proteins to simulation monomer IDs:
    def map_input_proteins_to_simulation_monomer_IDs(self, input_protein_list, full_dataframe):
        # map the input protein list common names to simulation monomer IDs:
        print("Mapping input protein list to simulation monomer IDs...")
        sim_monomer_id_with_compartment_to_input_protein_dict = {}
        for protein in input_protein_list:
            if "[" in protein:
                protein = protein[:-3]
            if protein in self.sim_monomer_ids_to_monomer_ids_with_compartment_tags_dict.keys():
                monomer_id = self.sim_monomer_ids_to_monomer_ids_with_compartment_tags_dict[protein]
                sim_monomer_id_with_compartment_to_input_protein_dict[monomer_id] = protein
            else:
                print("The input protein", protein,
                      "did not map to any simulation monomer IDs.")

        # find the relevant rows in the full dataframe based on which proteins were mapped (this will filter out proteins there is no validation data for):
        mapped_input_proteins_df = full_dataframe[
            full_dataframe['monomer_id'].isin(
                sim_monomer_id_with_compartment_to_input_protein_dict.keys())]

        print("Mapped", len(sim_monomer_id_with_compartment_to_input_protein_dict.keys()),
              f"input proteins (of {len(input_protein_list)} inputs total) to simulation monomer IDs."
              f" Of those proteins, only {len(mapped_input_proteins_df)} showed up in the validation data as well.")

        return mapped_input_proteins_df






    # Create a function that creates a table of relevant simulation protein info for a given list of monomer IDs:
    def create_simulation_protein_info_table(
            self, RVS_uniprot_ids_to_sim_monomer_ids_matches_dict,
            raw_validation_source_uniprot_ids_to_schmidt_common_names,
            raw_validation_source_uniprot_ids_to_counts_dict):

        # NOTE: specifically need the matches between the simulation monomer IDs
        # and the validation dataset monomer IDs to create this table!

        # Generate a table contiaining data for each monomer ID that overlaps
        # between the simulation and RVS:
        RVS_sim_data_table = []
        for uniprot_id in RVS_uniprot_ids_to_sim_monomer_ids_matches_dict.keys():
            monomer_id = RVS_uniprot_ids_to_sim_monomer_ids_matches_dict[uniprot_id]

            # extract the common name associated with the protein from rnas.tsv:
            sim_common_name = self.sim_monomer_ids_to_common_names_dict[monomer_id]
            monomer_id_with_compartment_tag = self.sim_monomer_ids_to_monomer_ids_with_compartment_tags_dict[monomer_id]
            description = self.sim_monomer_ids_to_descriptions_dict[monomer_id]
            half_life = self.sim_monomer_ids_to_half_lifes_dict[monomer_id]
            avg_total_count = self.sim_monomer_ids_to_avg_total_protein_counts_dict[monomer_id]
            avg_free_count = self.sim_monomer_ids_to_avg_free_protein_counts_dict[monomer_id]
            avg_complex_count = self.sim_monomer_ids_to_avg_complex_protein_counts_dict[monomer_id]

            # extract the common name associated with the protein directly from the validation dataset:
            RVS_common_name = raw_validation_source_uniprot_ids_to_schmidt_common_names[uniprot_id]
            RVS_count = raw_validation_source_uniprot_ids_to_counts_dict[uniprot_id]

            RVS_sim_data_table.append({
                'uniprot_id': uniprot_id,
                'common_name': sim_common_name,
                'monomer_id': monomer_id_with_compartment_tag,
                'description': description,
                'simulation_half_life': f'{half_life} mins, {half_life/60:.2f} hrs',
                'RVS_count': RVS_count,
                'avg_total_count': avg_total_count,
                'avg_free_count': avg_free_count,
                'avg_complex_count': avg_complex_count,
                'RVS_common_name': RVS_common_name
            })

        # convert to a dataframe:
        RVS_sim_data_df = pd.DataFrame(RVS_sim_data_table)

        return RVS_sim_data_df

    # Define a hovertext function:
    def generate_hovertext(self, dataframe, include_STDs=False, generation_averaged=False):
        hovertext = dataframe.apply(lambda row:
                                          f"Monomer ID: {row['monomer_id']}"
                                          f"<br>UniProt ID: {row['uniprot_id']}"
                                          f"<br>Simulation common name: {row['common_name']}"
                                          f"<br>Validation source common name: {row['RVS_common_name']}"
                                          f"<br>Description: {row['description']}"
                                          f"<br>Simulation half life: {row['simulation_half_life']}"
                                          f"<br>Validation source count: {row['RVS_count']}"
                                          f"<br>Avg. total simulation count: {row['avg_total_count']}"
                                          f"<br>Avg. complexed count: {row['avg_complex_count']}"
                                          f"<br>Avg. free count: {row['avg_free_count']}",
                                          axis=1)
        if include_STDs == True:
            hovertext = dataframe.apply(lambda row:
                                              f"Monomer ID: {row['monomer_id']}"
                                              f"<br>UniProt ID: {row['uniprot_id']}"
                                              f"<br>Simulation common name: {row['common_name']}"
                                              f"<br>Validation source common name: {row['RVS_common_name']}"
                                              f"<br>Description: {row['description']}"
                                              f"<br>Simulation half life: {row['simulation_half_life']}"
                                              f"<br>Validation source count: {row['RVS_count']}"
                                              f"<br>Avg. total simulation count: {row['avg_total_count']} ± {row['std_total_count']}"
                                              f"<br>Avg. complexed count: {row['avg_complex_count']} ± {row['std_complexed_count']}"
                                              f"<br>Avg. free count: {row['avg_free_count']} ± {row['std_free_count']}",
                                              axis=1)
        if generation_averaged == True:
            hovertext = dataframe.apply(lambda row:
                                        f"Monomer ID: {row['monomer_id']}"
                                        f"<br>UniProt ID: {row['uniprot_id']}"
                                        f"<br>Simulation common name: {row['common_name']}"
                                        f"<br>Validation source common name: {row['RVS_common_name']}"
                                        f"<br>Description: {row['description']}"
                                        f"<br>Simulation half life: {row['simulation_half_life']}"
                                        f"<br>Validation source count: {row['RVS_count']}"
                                        f"<br>Avg. total simulation count: {row['avg_total_count_gen']} ± {row['std_total_count_gen']}",
                                        axis=1)
        return hovertext

    # With the matching simulation and validation data obtained, make appropreite comparison plots:
    def compare_simulation_counts_to_raw_validation_source(
            self, plotOutDir, validation_source_name, validation_source_name_short,
            RVS_uniprot_ids_to_monomer_IDs,
            RVS_uniprot_ids_to_schmidt_common_names,
            RVS_uniprot_ids_to_counts_dict,
            plot_proteins_of_interest=False):

        print(f"Generating comparison plot for {self.sim_name} vs. {validation_source_name_short}.")

        # Obtain the mapping of uniprot IDs to simulation monomer IDs for the overlapping proteins:
        RVS_uniprot_ids_to_sim_monomer_ids_dict = self.match_validation_dataset_monomer_IDs_to_simulation_monomer_IDs(RVS_uniprot_ids_to_monomer_IDs, RVS_uniprot_ids_to_schmidt_common_names)


        # create a table of relevant simulation protein info for the overlapping proteins:
        RVS_sim_data_df = (
            self.create_simulation_protein_info_table(
                RVS_uniprot_ids_to_sim_monomer_ids_dict,
                RVS_uniprot_ids_to_schmidt_common_names,
                RVS_uniprot_ids_to_counts_dict))


        # Generate the plot:
        fig = go.Figure()

        # Compute log10 values for simulation and validation protein counts:
        x = np.log10(RVS_sim_data_df['RVS_count'].values + 1)
        y = np.log10(RVS_sim_data_df['avg_total_count'].values + 1)

        # Compute linear trendline
        z = np.polyfit(x, y, 1)
        p = np.poly1d(z)
        trendline_y = p(x)

        # Compute linear trendline for counts above log10(30) as done in Macklin et al. 2020:
        # NOTE: (+1 bc log(0) is undefined)
        above_30_idx = np.where((x > np.log10(30 + 1)) & (y > np.log10(30 + 1)))
        x_above_30 = x[above_30_idx]
        y_above_30 = y[above_30_idx]
        z_above_30 = np.polyfit(x_above_30, y_above_30, 1)
        p_above_30 = np.poly1d(z_above_30)
        trendline_y_above_30 = p_above_30(x)

        # Compute the pearson r, R2, and the coefficient of determination R2:
        r_value, p_val = pearsonr(x_above_30, y_above_30)
        pr2 = r_value ** 2
        COD_r2 = r2_score(x_above_30, y_above_30) # COD R2

        # Add total counts scatter data:
        hovertext = self.generate_hovertext(RVS_sim_data_df)
        fig.add_trace(
            go.Scatter(x=x, y=y, hovertext=hovertext, mode='markers',
                       name=f"Monomer Counts", marker=dict(color='lightseagreen'), opacity=0.7))

        protein_string = ''
        if plot_proteins_of_interest == True:
            if PROTEINS_OF_INTEREST != []:
                proteins_of_interest_df = self.map_input_proteins_to_simulation_monomer_IDs(
                    PROTEINS_OF_INTEREST, RVS_sim_data_df)
                # Add proteins of interest scatter data:
                hovertext_poi = self.generate_hovertext(proteins_of_interest_df)
                x_poi = np.log10(proteins_of_interest_df['RVS_count'].values + 1)
                y_poi = np.log10(proteins_of_interest_df['avg_total_count'].values + 1)
                fig.add_trace(
                    go.Scatter(x=x_poi, y=y_poi, hovertext=hovertext_poi, mode='markers',
                               name=f"Proteins of Interest (n={len(proteins_of_interest_df)})", marker=dict(color='red', size=6)))
                protein_string = '_'
                for p in range(len(proteins_of_interest_df)):
                    common_name = proteins_of_interest_df.iloc[p]["common_name"]
                    protein_string = protein_string + '_' + str(common_name)
            else:
                print("HIGHLIGHT_PROTEINS_OF_INTEREST is set to True but no "
                      "proteins were listed in PROTEINS_OF_INTEREST, "
                      "so individual proteins will not be highlighted.")


        # Add linear trendline:
        fig.add_trace(
            go.Scatter(x=x, y=trendline_y, mode='lines',
                       name=f'Linear fit (all data): {p}',
                       line=dict(color='green')))

        # Add linear trendline for counts above 30:
        fig.add_trace(
            go.Scatter(x=x, y=trendline_y_above_30, mode='lines',
                       name=f'Linear fit (counts > 30): {p_above_30}',
                       line=dict(color='pink')))

        # Add a y=x line:
        fig.add_trace(
            go.Scatter(x=[0, 6], y=[0, 6], mode="lines",
                       line=go.scatter.Line(color="black", dash="dash"),
                       opacity=0.2, name="y=x"))


        # Update layout
        fig.update_traces(marker_size=3)
        fig.update_layout(
            title=f"Simulation Protein Counts vs. Validation Protein Counts<br>"
                  f"Sim ID: {self.sim_name} (averaged over {self.total_cells} cells), "
                  f"Validation dataset: {validation_source_name}<br>"
                  f"Pearson R<sup>2</sup> for counts > 30: {round(pr2, 3)}, n={len(above_30_idx[0])} (of {len(x)} total plotted)",
            title_font=dict(size=8),
            xaxis_title="log₁₀(Validation Protein Counts + 1)",
            yaxis_title="log₁₀(Simulation Protein Counts + 1)",
            autosize=False,
            width=900,
            height=600,
            plot_bgcolor='white',  # Set the plot area background color to white
            paper_bgcolor='white'  # Set the entire graph background to white
        )

        # Define the text to display
        text = (
            f'Pearson R (counts > 30): {round(r_value, 3)}<br>'
            f'Pearson R<sup>2</sup> (counts > 30): {round(pr2, 3)}<br>'
            f'Coefficient of determination R<sup>2</sup> (counts > 30): {round(COD_r2, 3)}'
        )

        # Get the maximum x and minimum y to position the text in the bottom-right
        x_max = x.max()
        y_min = y.min()

        # Adjust x_max and y_min slightly outside the actual graph boundaries
        text_offset_x = 0.1
        text_offset_y = 0.05

        # Adding text annotation to the bottom right
        fig.add_annotation(
            x=x_max + text_offset_x,  # Move the x position slightly to the right
            y=y_min - text_offset_y,  # Move the y position slightly below
            text=text,
            showarrow=False,
            bgcolor='rgba(255, 255, 255, 0.8)',
            bordercolor='rgba(0, 0, 0, 0.5)',
            borderwidth=1,
            borderpad=4,
            align='right',
            font=dict(size=10, color='gray'),
            xref='x',
            yref='y',
        )

        # save the figure as an html:
        plot_name = (f"proteomics_comparison_to_raw_validation_source_"
                     f"{self.sim_name}_vs_{validation_source_name_short}{protein_string}.html")
        fig.write_html(os.path.join(plotOutDir, plot_name))


    # generate the same graph as above but allow for the manually mapped proteins to show up:
    def compare_simulation_counts_to_raw_validation_source_with_manual_mappings(
            self, plotOutDir, validation_source_name, validation_source_name_short,
            RVS_uniprot_ids_to_monomer_IDs,
            RVS_uniprot_ids_to_schmidt_common_names,
            RVS_uniprot_ids_to_counts_dict,
            unmapped_uniprot_ids_to_monomer_ids, plot_proteins_of_interest=False):

        print(f"Generating comparison plot for {self.sim_name} vs. {validation_source_name_short} (including manual mappings).")

        # Obtain the mapping of uniprot IDs to simulation monomer IDs for the overlapping proteins:
        RVS_uniprot_ids_to_sim_monomer_ids_dict = self.match_validation_dataset_monomer_IDs_to_simulation_monomer_IDs(RVS_uniprot_ids_to_monomer_IDs, RVS_uniprot_ids_to_schmidt_common_names)

        # create a table of relevant simulation protein info for the overlapping proteins:
        RVS_sim_data_df = (
            self.create_simulation_protein_info_table(
                RVS_uniprot_ids_to_sim_monomer_ids_dict,
                RVS_uniprot_ids_to_schmidt_common_names,
                RVS_uniprot_ids_to_counts_dict))

        # generate a table for the manually mapped proteins as well:
        RVS_unmapped_uniprot_ids_to_sim_monomer_ids_dict = self.match_validation_dataset_monomer_IDs_to_simulation_monomer_IDs(
            unmapped_uniprot_ids_to_monomer_ids, RVS_uniprot_ids_to_schmidt_common_names)

        RVS_extra_data_df = self.create_simulation_protein_info_table(RVS_unmapped_uniprot_ids_to_sim_monomer_ids_dict,
                                                  RVS_uniprot_ids_to_schmidt_common_names,
                                                  RVS_uniprot_ids_to_counts_dict)

        # append the RVS_extra_data_df to the RVS_sim_data_df
        RVS_sim_data_df = pd.concat([RVS_sim_data_df, RVS_extra_data_df], ignore_index=True)


        # Generate the plot:
        fig = go.Figure()

        # Compute log10 values for simulation and validation protein counts:
        x = np.log10(RVS_sim_data_df['RVS_count'].values + 1)
        y = np.log10(RVS_sim_data_df['avg_total_count'].values + 1)

        # Generate data with manual mappings highlighted
        x2 = np.log10(RVS_extra_data_df['RVS_count'].values + 1)
        y2 = np.log10(RVS_extra_data_df['avg_total_count'].values + 1)

        # Compute linear trendline
        z = np.polyfit(x, y, 1)
        p = np.poly1d(z)
        trendline_y = p(x)

        # Compute linear trendline for counts above log10(30) as done in Macklin et al. 2020:
        # NOTE: (+1 bc log(0) is undefined)
        above_30_idx = np.where((x > np.log10(30 + 1)) & (y > np.log10(30 + 1)))
        x_above_30 = x[above_30_idx]
        y_above_30 = y[above_30_idx]
        z_above_30 = np.polyfit(x_above_30, y_above_30, 1)
        p_above_30 = np.poly1d(z_above_30)
        trendline_y_above_30 = p_above_30(x)

        # Compute the pearson r, R2, and the coefficient of determination R2:
        r_value, p_val = pearsonr(x_above_30, y_above_30)
        pr2 = r_value ** 2
        COD_r2 = r2_score(x_above_30, y_above_30) # COD R2

        # Add total counts scatter data:
        hovertext = self.generate_hovertext(RVS_sim_data_df)
        fig.add_trace(
            go.Scatter(x=x, y=y, hovertext=hovertext, mode='markers',
                       name=f"Monomer Counts (n={len(x)})", marker=dict(color='lightseagreen'), opacity=0.7))

        protein_string = ''
        if plot_proteins_of_interest == True:
            if PROTEINS_OF_INTEREST != []:
                proteins_of_interest_df = self.map_input_proteins_to_simulation_monomer_IDs(
                    PROTEINS_OF_INTEREST, RVS_sim_data_df)
                # Add proteins of interest scatter data:
                hovertext_poi = self.generate_hovertext(proteins_of_interest_df)
                x_poi = np.log10(proteins_of_interest_df['RVS_count'].values + 1)
                y_poi = np.log10(proteins_of_interest_df['avg_total_count'].values + 1)
                fig.add_trace(
                    go.Scatter(x=x_poi, y=y_poi, hovertext=hovertext_poi, mode='markers',
                               name=f"Proteins of Interest (n={len(proteins_of_interest_df)})", marker=dict(color='red', size=6)))
                protein_string = '_'
                for p in range(len(proteins_of_interest_df)):
                    common_name = proteins_of_interest_df.iloc[p]["common_name"]
                    protein_string = protein_string + '_' + str(common_name)
            else:
                print("HIGHLIGHT_PROTEINS_OF_INTEREST is set to True but no "
                      "proteins were listed in PROTEINS_OF_INTEREST, "
                      "so individual proteins will not be highlighted.")


        # Add manually mapped proteins scatter data:
        manual_data_hovertext = self.generate_hovertext(RVS_extra_data_df)
        fig.add_trace(go.Scatter(x=x2, y=y2, hovertext=manual_data_hovertext, mode='markers',
                                 name=f'Manually Mapped Proteins (n={len(x2)})',
                                 marker=dict(color='orange', size=8, symbol='circle-open')))

        # Add linear trendline:
        fig.add_trace(
            go.Scatter(x=x, y=trendline_y, mode='lines',
                       name=f'Linear fit (all data): {p}',
                       line=dict(color='green')))

        # Add linear trendline for counts above 30:
        fig.add_trace(
            go.Scatter(x=x, y=trendline_y_above_30, mode='lines',
                       name=f'Linear fit (counts > 30): {p_above_30}',
                       line=dict(color='pink')))

        # Add a y=x line:
        fig.add_trace(
            go.Scatter(x=[0, 6], y=[0, 6], mode="lines",
                       line=go.scatter.Line(color="black", dash="dash"),
                       opacity=0.2, name="y=x"))


        # Update layout
        fig.update_traces(marker_size=3)
        fig.update_layout(
            title=f"Simulation Protein Counts vs. Validation Protein Counts<br>"
                  f"Sim ID: {self.sim_name} (averaged over {self.total_cells} cells), "
                  f"Validation dataset: {validation_source_name}<br>"
                  f"Pearson R<sup>2</sup> for counts > 30: {round(pr2, 3)}, n={len(above_30_idx[0])} (of {len(x)} total plotted)",
            title_font=dict(size=8),
            xaxis_title="log₁₀(Validation Protein Counts + 1)",
            yaxis_title="log₁₀(Simulation Protein Counts + 1)",
            autosize=False,
            width=900,
            height=600,
            plot_bgcolor='white',  # Set the plot area background color to white
            paper_bgcolor='white'  # Set the entire graph background to white
        )

        # Define the text to display
        text = (
            f'Pearson R (counts > 30): {round(r_value, 3)}<br>'
            f'Pearson R<sup>2</sup> (counts > 30): {round(pr2, 3)}<br>'
            f'Coefficient of determination R<sup>2</sup> (counts > 30): {round(COD_r2, 3)}'
        )

        # Get the maximum x and minimum y to position the text in the bottom-right
        x_max = x.max()
        y_min = y.min()

        # Adjust x_max and y_min slightly outside the actual graph boundaries
        text_offset_x = 0.1
        text_offset_y = 0.05

        # Adding text annotation to the bottom right
        fig.add_annotation(
            x=x_max + text_offset_x,  # Move the x position slightly to the right
            y=y_min - text_offset_y,  # Move the y position slightly below
            text=text,
            showarrow=False,
            bgcolor='rgba(255, 255, 255, 0.8)',
            bordercolor='rgba(0, 0, 0, 0.5)',
            borderwidth=1,
            borderpad=4,
            align='right',
            font=dict(size=10, color='gray'),
            xref='x',
            yref='y',
        )

        # save the figure as an html:
        plot_name = (f"proteomics_comparison_to_raw_validation_source_"
                     f"{self.sim_name}_vs_{validation_source_name_short}_with_manual_mappings{protein_string}.html")
        fig.write_html(os.path.join(plotOutDir, plot_name))


    # Create a function that makes a table of the validation data against another validation data source:
    def prep_validation_datasets_for_comparison_with_other_validation_data(
            self, common_name_to_counts_dict_1, common_name_to_sim_monomer_id_1,
            common_name_to_counts_dict_2, common_name_to_sim_monomer_id_2):
        # Temporarily build monomer_id_to_common_name dict for dataset 1:
        monomer_id_to_common_name_1 = {}
        for common_name in common_name_to_sim_monomer_id_1.keys():
            monomer_id = common_name_to_sim_monomer_id_1[common_name]
            monomer_id_to_common_name_1[monomer_id] = common_name

        # Temporarily build monomer_id_to_common_name dict for dataset 2:
        monomer_id_to_common_name_2 = {}
        for common_name in common_name_to_sim_monomer_id_2.keys():
            monomer_id = common_name_to_sim_monomer_id_2[common_name]
            monomer_id_to_common_name_2[monomer_id] = common_name

        # Obtain overlapping protein counts between the two validation datasets
        sim_monomer_ids_in_both_datasets = set(monomer_id_to_common_name_1.keys()).intersection(
            set(monomer_id_to_common_name_2.keys()))
        print("Number of shared monomer IDs between the two validation datasets:", len(sim_monomer_ids_in_both_datasets))

        # Create dictionaries mapping the counts and common names to counts for the overlapping proteins
        dataset1_common_names_to_counts_dict = {}
        dataset1_common_names_to_sim_monomer_ids_dict = {}
        dataset2_common_names_to_counts_dict = {}
        dataset2_common_names_to_sim_monomer_ids_dict = {}
        for monomer_id in sim_monomer_ids_in_both_datasets:
            common_name_1 = monomer_id_to_common_name_1[monomer_id]
            dataset1_common_names_to_counts_dict[common_name_1] = common_name_to_counts_dict_1[common_name_1]
            dataset1_common_names_to_sim_monomer_ids_dict[common_name_1] = monomer_id

            common_name_2 = monomer_id_to_common_name_2[monomer_id]
            dataset2_common_names_to_counts_dict[common_name_2] = common_name_to_counts_dict_2[common_name_2]
            dataset2_common_names_to_sim_monomer_ids_dict[common_name_2] = monomer_id

        # Make a table with all the data:
        validation_datasets_comparison_table = []
        for monomer_id in sim_monomer_ids_in_both_datasets:
            common_name_1 = monomer_id_to_common_name_1[monomer_id]
            common_name_2 = monomer_id_to_common_name_2[monomer_id]
            count_1 = common_name_to_counts_dict_1[common_name_1]
            count_2 = common_name_to_counts_dict_2[common_name_2]

            validation_datasets_comparison_table.append({
                'monomer_id': monomer_id,
                'dataset1_common_name': common_name_1,
                'dataset2_common_name': common_name_2,
                'dataset1_count': count_1,
                'dataset2_count': count_2
            })

        # convert to a dataframe:
        validation_datasets_comparison_df = pd.DataFrame(validation_datasets_comparison_table)

        return validation_datasets_comparison_df

    # Function for comparing validaiton dataset to validation dataset:
    def compare_validation_datasets(self, plotOutDir, validation_source_name_1, validation_source_name_short_1,
            uniprot_IDs_to_common_names_1, uniprot_IDs_to_monomer_IDs_1, uniprot_IDs_to_counts_1,
            validation_source_name_2, validation_source_name_short_2,
            uniprot_IDs_to_common_names_2, uniprot_IDs_to_monomer_IDs_2, uniprot_IDs_to_counts_2):

        print(f"Plotting validation source comparison between "
              f"{validation_source_name_1} and {validation_source_name_2}.")

        # Determine which proteins overlap between the datasets:
        overlapping_uniprot_ids = set(uniprot_IDs_to_common_names_1.keys()).intersection(uniprot_IDs_to_common_names_2.keys())

        # Determine if the monomer also shows up in the simulation data:
        uniprot_IDs_to_sim_monomer_IDs_existance = {}
        uniprot_IDs_to_sim_common_name = {}
        uniprot_IDs_to_sim_monomer_counts = {}
        uniprot_IDs_to_sim_monomer_description = {}
        for id in overlapping_uniprot_ids:
            monomer_id = uniprot_IDs_to_monomer_IDs_1[id]
            if monomer_id in self.all_monomer_ids:
                uniprot_IDs_to_sim_monomer_IDs_existance[id] = True
                uniprot_IDs_to_sim_common_name[id] = self.sim_monomer_ids_to_common_names_dict[monomer_id]
                uniprot_IDs_to_sim_monomer_counts[id] = self.sim_monomer_ids_to_avg_total_protein_counts_dict[monomer_id]
                uniprot_IDs_to_sim_monomer_description[id] = self.sim_monomer_ids_to_descriptions_dict[monomer_id]
            else:
                uniprot_IDs_to_sim_monomer_IDs_existance[id] = False
                uniprot_IDs_to_sim_common_name[id] = "NA"
                uniprot_IDs_to_sim_monomer_counts[id] = "NA"
                uniprot_IDs_to_sim_monomer_description[id] = "NA"

        # Generate a large table with all the relevant info:
        validation_dataset_comparison_table = []
        for id in overlapping_uniprot_ids:
            monomer_id_1 = uniprot_IDs_to_monomer_IDs_1[id]
            monomer_id_2 = uniprot_IDs_to_monomer_IDs_2[id]
            monomer_id_existance_in_sim_data = uniprot_IDs_to_sim_monomer_IDs_existance[id]
            common_name_1 = uniprot_IDs_to_common_names_1[id]
            common_name_2 = uniprot_IDs_to_common_names_2[id]
            sim_common_name = uniprot_IDs_to_sim_common_name[id]
            count_1 = uniprot_IDs_to_counts_1[id]
            count_2 = uniprot_IDs_to_counts_2[id]
            sim_count = uniprot_IDs_to_sim_monomer_counts[id]
            sim_description = uniprot_IDs_to_sim_monomer_description[id]

            validation_dataset_comparison_table.append({
                'uniprot_id': id,
                'monomer_id_1': monomer_id_1,
                'monomer_id_2': monomer_id_2,
                'monomer_id_in_sim': monomer_id_existance_in_sim_data,
                'dataset1_common_name': common_name_1,
                'dataset2_common_name': common_name_2,
                'sim_common_name': sim_common_name,
                'dataset1_count': count_1,
                'dataset2_count': count_2,
                'sim_count': sim_count,
                'sim_description': sim_description
            })

        # convert to a dataframe:
        validation_dataset_comparison_df = pd.DataFrame(validation_dataset_comparison_table)

        # Generate the plot:
        fig = go.Figure()

        # Compute log10 values for simulation validation and raw validation protein counts:
        x = np.log10(validation_dataset_comparison_df['dataset1_count'].values + 1)
        y = np.log10(validation_dataset_comparison_df['dataset2_count'].values + 1)

        # Compute linear trendline
        z = np.polyfit(x, y, 1)
        p = np.poly1d(z)
        trendline_y = p(x)

        # Compute linear trendline for counts above log10(30+1): (+1 bc log(0) is undefined)
        above_30_idx = np.where((x > np.log10(30 + 1)) & (y > np.log10(30 + 1)))
        x_above_30 = x[above_30_idx]
        y_above_30 = y[above_30_idx]
        z_above_30 = np.polyfit(x_above_30, y_above_30, 1)
        p_above_30 = np.poly1d(z_above_30)
        trendline_y_above_30 = p_above_30(x)

        # Compute the pearson r, R2, and the coefficient of determination R2:
        r_value, p_val = pearsonr(x_above_30, y_above_30)
        pr2 = r_value ** 2
        COD_r2 = r2_score(x_above_30, y_above_30)  # COD R2

        # Define the hover text for the plot output:
        hovertext = validation_dataset_comparison_df.apply(lambda
                                                                row:
                                                            f"Uniprot ID: {row['uniprot_id']}"
                                                            f"<br>{validation_source_name_short_1} Monomer ID: {row['monomer_id_1']}"
                                                            f"<br>{validation_source_name_short_2} Monomer ID: {row['monomer_id_2']}"
                                                            f"<br>Monomer found in {self.sim_name}: {row['monomer_id_in_sim']}"
                                                            f"<br>{validation_source_name_short_1} common name: {row['dataset1_common_name']}"
                                                            f"<br>{validation_source_name_short_2} common name: {row['dataset2_common_name']}"
                                                            f"<br>{self.sim_name} common name: {row['sim_common_name']}"
                                                            f"<br>{validation_source_name_short_1} count: {row['dataset1_count']}"
                                                            f"<br>{validation_source_name_short_2} count: {row['dataset2_count']}"
                                                            f"<br>{self.sim_name} average total count: {row['sim_count']}"
                                                            f"<br>{self.sim_name} monomer description: {row['sim_description']}",
                                                            axis=1)

        # Add total counts scatter data:
        fig.add_trace(
            go.Scatter(x=x, y=y, hovertext=hovertext, mode='markers',
                       name=f"Monomer Counts"))

        # Add linear trendline:
        fig.add_trace(
            go.Scatter(x=x, y=trendline_y, mode='lines',
                       name=f'Linear fit (all data): {p}',
                       line=dict(color='green')))

        # Add linear trendline for counts above 30:
        fig.add_trace(
            go.Scatter(x=x, y=trendline_y_above_30, mode='lines',
                       name=f'Linear fit (counts > 30): {p_above_30}',
                       line=dict(color='pink')))

        # Add a y=x line:
        fig.add_trace(
            go.Scatter(x=[0, 6], y=[0, 6], mode="lines",
                       line=go.scatter.Line(color="black", dash="dash"),
                       opacity=0.2, name="y=x"))

        # Update layout
        fig.update_traces(marker_size=3)
        fig.update_layout(
            title=f"{validation_source_name_1} ({validation_source_name_short_1}) vs. {validation_source_name_2} ({validation_source_name_short_2})<br>"
                  f"Sim ID: {self.sim_name} (averaged over {self.total_cells} cells), "
                  f"Pearson R<sup>2</sup> for counts > 30: {round(pr2, 3)}, n={len(above_30_idx[0])} (of {len(x)} total plotted)",
            title_font=dict(size=8),
            xaxis_title=f"log₁₀({validation_source_name_short_1} Protein Counts)",
            yaxis_title=f"log₁₀({validation_source_name_short_2} Protein Counts)",
            autosize=False,
            width=900,
            height=600,
            plot_bgcolor='white',  # Set the plot area background color to white
            paper_bgcolor='white'  # Set the entire graph background to white
        )

        # Define the text to display
        text = (
            f'Pearson R (counts > 30): {round(r_value, 3)}<br>'
            f'Pearson R<sup>2</sup> (counts > 30): {round(pr2, 3)}<br>'
            f'Coefficient of determination R<sup>2</sup> (counts > 30): {round(COD_r2, 3)}'
        )

        # Get the maximum x and minimum y to position the text in the bottom-right
        x_max = x.max()
        y_min = y.min()

        # Adjust x_max and y_min slightly outside the actual graph boundaries
        text_offset_x = 0.1
        text_offset_y = 0.05

        # Adding text annotation to the bottom right
        fig.add_annotation(
            x=x_max + text_offset_x,  # Move the x position slightly to the right
            y=y_min - text_offset_y,  # Move the y position slightly below
            text=text,
            showarrow=False,
            bgcolor='rgba(255, 255, 255, 0.8)',
            bordercolor='rgba(0, 0, 0, 0.5)',
            borderwidth=1,
            borderpad=4,
            align='right',
            font=dict(size=10, color='gray'),
            xref='x',
            yref='y',
        )

        # save the figure as an html:
        plot_name = (f"validation_to_validation_source_comparison_"
                     f"{validation_source_name_short_1}_vs_{validation_source_name_short_2}.html")
        fig.write_html(os.path.join(plotOutDir, plot_name))

    # Create the main function that generates the validation comparison plots:
    def compare_raw_validation_dataset_to_saved_simulation_validation_dataset(
            self, simDataFile, validationDataFile, plotOutDir,
            validation_source_name, validation_source_name_short, uniprot_IDs_to_common_names,
            uniprot_IDs_to_monomer_IDs, uniprot_IDs_to_counts):

        print(f"Plotting comparison plot of the validation data saved with the simulation data (Schmidt et al. 2015 BW25113 Glucose data) with raw {validation_source_name}")

        # Extract the validation data automatically saved with the simulation  (Schmidt et al. 2015 BW25113):
        (sim_validation_common_names_to_counts_dict,
         sim_validation_common_names_to_monomer_ids_dict) = self.prep_validation_data_for_comparison(
            simDataFile, validationDataFile)

        # Reverse map monomer IDs to respective common names:
        monomer_IDs_to_sim_validation_common_names = {}
        for common_name in sim_validation_common_names_to_monomer_ids_dict.keys():
            monomer_ID = sim_validation_common_names_to_monomer_ids_dict[common_name]
            monomer_IDs_to_sim_validation_common_names[monomer_ID] = common_name

        # Convert RVS dictionaries to be mapped on monomers:
        RVS_monomer_IDs_to_uniprot_IDs = {}
        for id in uniprot_IDs_to_monomer_IDs.keys():
            monomer_ID = uniprot_IDs_to_monomer_IDs[id]
            RVS_monomer_IDs_to_uniprot_IDs[monomer_ID] = id

        # check matching data:
        overlapping_monomer_IDs = set(RVS_monomer_IDs_to_uniprot_IDs.keys()).intersection(set(monomer_IDs_to_sim_validation_common_names.keys()))

        # Generate a table of the relevant information:
        validation_dataset_comparison_table = []
        for id in overlapping_monomer_IDs:
            sim_monomer_ID = self.sim_monomer_ids_to_monomer_ids_with_compartment_tags_dict[id]
            RVS_uniprot_ID = RVS_monomer_IDs_to_uniprot_IDs[id]
            sim_common_name = monomer_IDs_to_sim_validation_common_names[id]
            RVS_common_name = uniprot_IDs_to_common_names[RVS_uniprot_ID]
            sim_validation_counts = sim_validation_common_names_to_counts_dict[sim_common_name]
            RVS_counts = uniprot_IDs_to_counts[RVS_uniprot_ID]
            sim_counts = self.sim_monomer_ids_to_avg_total_protein_counts_dict[id]
            sim_monomer_description = self.sim_monomer_ids_to_descriptions_dict[id]

            # Since the validation data is saved with the simulation here,
            # we can assume any monomer ID matches have matches for common names
            # and counts in the simulation as well.
            validation_dataset_comparison_table.append({
                'monomer_id': sim_monomer_ID,
                'uniprot_id': RVS_uniprot_ID,
                'sim_common_name': sim_common_name,
                'RVS_common_name': RVS_common_name,
                'sim_validation_counts': sim_validation_counts,
                'RVS_counts': RVS_counts,
                'sim_counts': sim_counts,
                'sim_monomer_description': sim_monomer_description
            })

        # Convert to a dataframe:
        validation_dataset_comparison_df = pd.DataFrame(validation_dataset_comparison_table)

        # Generate the plot:
        fig = go.Figure()


        # Compute log10 values for simulation validation and raw validation protein counts:
        x = np.log10(validation_dataset_comparison_df['sim_validation_counts'].values + 1)
        y = np.log10(validation_dataset_comparison_df['RVS_counts'].values + 1)

        # Compute linear trendline
        z = np.polyfit(x, y, 1)
        p = np.poly1d(z)
        trendline_y = p(x)

        # Compute linear trendline for counts above log10(30+1): (+1 bc log(0) is undefined)
        above_30_idx = np.where((x > np.log10(30 + 1)) & (y > np.log10(30 + 1)))
        x_above_30 = x[above_30_idx]
        y_above_30 = y[above_30_idx]
        z_above_30 = np.polyfit(x_above_30, y_above_30, 1)
        p_above_30 = np.poly1d(z_above_30)
        trendline_y_above_30 = p_above_30(x)

        # Compute the pearson r, R2, and the coefficient of determination R2:
        r_value, p_val = pearsonr(x_above_30, y_above_30)
        pr2 = r_value ** 2
        COD_r2 = r2_score(x_above_30, y_above_30)  # COD R2

        # Define the hover text for the plot output:
        hovertext = validation_dataset_comparison_df.apply(lambda
                                                                   row:
                                                               f"{validation_source_name_short} Uniprot ID: {row['uniprot_id']}"
                                                               f"<br>Monomer ID: {row['monomer_id']}"
                                                               f"<br>{validation_source_name_short} common name: {row['RVS_common_name']}"
                                                               f"<br>{self.sim_name} common name: {row['sim_common_name']}"
                                                               f"<br>{validation_source_name_short} counts: {row['RVS_counts']}"
                                                               f"<br>Saved validation data (Schmidt2015_ST6_BW) counts: {row['sim_validation_counts']}"
                                                               f"<br>{self.sim_name} average total counts: {row['sim_counts']}"
                                                               f"<br>{self.sim_name} monomer description: {row['sim_monomer_description']}",
                                                               axis=1)

        # Add total counts scatter data:
        fig.add_trace(
            go.Scatter(x=x, y=y, hovertext=hovertext, mode='markers',
                           name=f"Monomer Counts"))

        # Add linear trendline:
        fig.add_trace(
            go.Scatter(x=x, y=trendline_y, mode='lines',
                           name=f'Linear fit (all data): {p}',
                           line=dict(color='green')))

        # Add linear trendline for counts above 30:
        fig.add_trace(
            go.Scatter(x=x, y=trendline_y_above_30, mode='lines',
                           name=f'Linear fit (counts > 30): {p_above_30}',
                           line=dict(color='pink')))

        # Add a y=x line:
        fig.add_trace(
            go.Scatter(x=[0, 6], y=[0, 6], mode="lines",
                           line=go.scatter.Line(color="black", dash="dash"),
                           opacity=0.2, name="y=x"))

        # Update layout
        fig.update_traces(marker_size=3)
        fig.update_layout(
                title=f"{validation_source_name} ({validation_source_name_short}) vs. Schmidt et al. 2015 BW25113 Glucose validation data saved with the simulation (Schmidt2015_ST6_BW)<br>"
                      f"Sim ID: {self.sim_name} (averaged over {self.total_cells} cells), "
                      f"Pearson R<sup>2</sup> for counts > 30: {round(pr2, 3)}, n={len(above_30_idx[0])} (of {len(x)} total plotted)",
                title_font=dict(size=8),
                xaxis_title=f"log₁₀(Schmidt2015_ST6_BW Protein Counts + 1)",
                yaxis_title=f"log₁₀({validation_source_name_short} Protein Counts + 1)",
                autosize=False,
                width=900,
                height=600,
                plot_bgcolor='white',  # Set the plot area background color to white
                paper_bgcolor='white'  # Set the entire graph background to white
            )

        # Define the text to display
        text = (
                f'Pearson R (counts > 30): {round(r_value, 3)}<br>'
                f'Pearson R<sup>2</sup> (counts > 30): {round(pr2, 3)}<br>'
                f'Coefficient of determination R<sup>2</sup> (counts > 30): {round(COD_r2, 3)}'
            )

        # Get the maximum x and minimum y to position the text in the bottom-right
        x_max = x.max()
        y_min = y.min()

        # Adjust x_max and y_min slightly outside the actual graph boundaries
        text_offset_x = 0.1
        text_offset_y = 0.05

        # Adding text annotation to the bottom right
        fig.add_annotation(
                x=x_max + text_offset_x,  # Move the x position slightly to the right
                y=y_min - text_offset_y,  # Move the y position slightly below
                text=text,
                showarrow=False,
                bgcolor='rgba(255, 255, 255, 0.8)',
                bordercolor='rgba(0, 0, 0, 0.5)',
                borderwidth=1,
                borderpad=4,
                align='right',
                font=dict(size=10, color='gray'),
                xref='x',
                yref='y',
            )

            # save the figure as an html:
        plot_name = (f"validation_to_validation_source_comparison_"
                         f"{validation_source_name_short}_vs_saved_Schmidt2015_ST6_BW.html")
        fig.write_html(os.path.join(plotOutDir, plot_name))


    # Get the standard deviations of the simulation data:
    def get_sim_data_standard_deviations(self, simDataFile, protein_df):
        # get the data:
        with open(simDataFile, 'rb') as f:
            sim_data = pickle.load(f)
        monomer_sim_data = (
            sim_data.process.translation.monomer_data.struct_array)

        # generate a dictionary mapping monomer IDs to their indices in the monomer_sim_data array:
        monomer_id_to_index = {}
        for index, monomer_id in enumerate(monomer_sim_data['id']):
            monomer_id_to_index[monomer_id] = index

        # find where each protein in the "monomer_id" column of protein_df is in the monomer_sim_data array:
        protein_indices_in_sim_data = []
        protein_indices_in_sim_data_dict = {}
        for monomer_id in protein_df['monomer_id']:
            if monomer_id in monomer_id_to_index:
                protein_indices_in_sim_data.append(monomer_id_to_index[monomer_id])
                protein_indices_in_sim_data_dict[monomer_id] = monomer_id_to_index[monomer_id]
            else:
                # only proteins that show up in the model should be present here so hopefully this isnt used
                protein_indices_in_sim_data.append(None)

        # Extract monomer count data for each protein:
        self.all_monomer_ids = monomer_sim_data['id']
        all_cells = self.ap.get_cells(
            generation=np.arange(IGNORE_FIRST_N_GENS,self.n_total_gens),
            only_successful=True)

        # AVERAGE OVER TOTAL SIMULATION DURATION, GIVING EACH TIME POINT EQUAL WEIGHT:
        # read in the data:
        total_counts = (
            read_stacked_columns(all_cells, 'MonomerCounts',
                                 'monomerCounts', ignore_exception=True))
        (free_counts,) = read_stacked_bulk_molecules(
            all_cells, self.all_monomer_ids, ignore_exception=True)
        complexed_counts = total_counts - free_counts

        # calculate the STDs and the averges:
        avg_total_counts = np.mean(total_counts, axis=0)
        std_total_counts = np.std(total_counts, axis=0)
        avg_free_counts = np.mean(free_counts, axis=0)
        std_free_counts = np.std(free_counts, axis=0)
        avg_complexed_counts = np.mean(complexed_counts, axis=0)
        std_complexed_counts = np.std(complexed_counts, axis=0)

        # add the above entries to the protein_df using the protein indices:
        avg_total_counts = avg_total_counts[[protein_indices_in_sim_data_dict[monomer_id] for monomer_id in protein_df['monomer_id']]]
        std_total_counts = std_total_counts[[protein_indices_in_sim_data_dict[monomer_id] for monomer_id in protein_df['monomer_id']]]
        avg_free_counts = avg_free_counts[[protein_indices_in_sim_data_dict[monomer_id] for monomer_id in protein_df['monomer_id']]]
        std_free_counts = std_free_counts[[protein_indices_in_sim_data_dict[monomer_id] for monomer_id in protein_df['monomer_id']]]
        avg_complexed_counts = avg_complexed_counts[[protein_indices_in_sim_data_dict[monomer_id] for monomer_id in protein_df['monomer_id']]]
        std_complexed_counts = std_complexed_counts[[protein_indices_in_sim_data_dict[monomer_id] for monomer_id in protein_df['monomer_id']]]

        # add to protein_df:
        protein_df['avg_total_count'] = avg_total_counts
        protein_df['std_total_count'] = std_total_counts
        protein_df['avg_free_count'] = avg_free_counts
        protein_df['std_free_count'] = std_free_counts
        protein_df['avg_complexed_count'] = avg_complexed_counts
        protein_df['std_complexed_count'] = std_complexed_counts


        # AVERAGE BY GENERATION, GIVING EACH CELL EQUAL WEIGHT:
        total_counts_per_gen = (
            read_stacked_columns(all_cells, 'MonomerCounts',
                                 'monomerCounts',
                                 fun=lambda x: np.mean(x[:], axis=0)))
        # TODO: when lisener-update is pushed, also track the free and complexed count reads

        # calculate the STDs and the averages:
        avg_total_counts_gen = np.mean(total_counts_per_gen, axis=0)
        std_total_counts_gen = np.std(total_counts_per_gen, axis=0)

        # use the dictionary of protein indices to get only the relevant proteins:
        avg_total_counts_gen = avg_total_counts_gen[[protein_indices_in_sim_data_dict[monomer_id] for monomer_id in protein_df['monomer_id']]]
        std_total_counts_gen = std_total_counts_gen[[protein_indices_in_sim_data_dict[monomer_id] for monomer_id in protein_df['monomer_id']]]

        # add to protein_df:
        protein_df['avg_total_count_gen'] = avg_total_counts_gen
        protein_df['std_total_count_gen'] = std_total_counts_gen

        return protein_df

    def plot_proteins_of_interest_averaged_over_all_time(self, simDataFile, plotOutDir,
                                           proteins_of_interest_df,
                                           validation_source_name,
                                           validation_source_name_short):

        protein_df = self.get_sim_data_standard_deviations(
            simDataFile, proteins_of_interest_df)

        # Extract the right transformation of the error bars:
        def log_error_bars(avg, std):
            avg_log = np.log10(np.clip(avg + 1, 1, None))  # Ensure input is at least 1

            # Verify there are no negative values before calculating lower bounds
            valid_mask = (avg - std) > 0

            # Calculate lower bounds while avoiding invalid log operations
            lower = np.where(valid_mask, np.log10(np.clip(avg + 1 - std, 1, None)),
                             avg_log)  # using clip to avoid taking log(0)
            upper = np.log10(np.clip(avg + 1 + std, 1, None))

            return avg_log, lower, upper


        # Get average logs and bounds for error bars
        total_avg_log, total_lower, total_upper = log_error_bars(protein_df['avg_total_count'].values,
                                                                 protein_df['std_total_count'].values)
        free_avg_log, free_lower, free_upper = log_error_bars(protein_df['avg_free_count'].values,
                                                              protein_df['std_free_count'].values)
        complexed_avg_log, complexed_lower, complexed_upper = log_error_bars(protein_df['avg_complexed_count'].values,
                                                                        protein_df['std_complexed_count'].values)

        # Create plot
        plt.figure(figsize=(8, 6), dpi=300)

        # Designate the x data as the validation counts and y data as the simulation counts:
        x_data = np.log10(protein_df['RVS_count'].values + 1)
        y_data = np.log10(protein_df['avg_total_count'].values + 1)

        # Now, calculate yerr for the error bars
        yerr_upper = total_upper - total_avg_log
        yerr_lower = total_avg_log - total_lower

        # Ensure yerr does not contain negative values
        yerr_lower = np.maximum(yerr_lower, 0)  # Replace negative errors with 0

        # Combine into a 2D array with shape (2, n)
        yerr = np.array([yerr_lower, yerr_upper])

        # Use plt.errorbar to create scatter plot with error bars
        plt.scatter(x_data, y_data, color='lightblue')
        plt.errorbar(x_data, y_data, yerr=yerr, fmt='o', color='lightblue', capsize=5)

        # Plot the common name next to each point
        for i in range(len(protein_df)):
            plt.text(x_data[i] + np.log10(1), y_data[i], protein_df['common_name'].values[i],
                     fontsize=6, ha='left', va='bottom')

        # Add labels and title
        plt.title(f"Simulation Protein Counts vs. Validation Protein Counts\n"
                  f"Sim ID: {self.sim_name} (averaged over all time points spanned by {self.total_cells} cells)\n "
                  f"Validation dataset: {validation_source_name}")
        plt.ylabel('Log10(Total Simulation Counts + 1)')
        plt.xlabel('Log10(Validation Counts + 1)')
        plt.grid(True)

        # Adjust y-axis to show log scale
        plt.ylim(0, max(total_avg_log) + 0.5)  # Add some space above the upper limit

        # add a y=x line:
        max_limit = max(max(x_data), max(y_data)) + 0.5
        plt.plot([0, max_limit], [0, max_limit], 'k--', lw=1)

        # obtain a string list of the common names to include in the title:
        protein_string = '_'
        for p in range(len(protein_df)):
            common_name = protein_df.iloc[p]["common_name"]
            protein_string = protein_string + '_' + str(common_name)

        # save the plot as a png:
        plot_name = (f"proteins_of_interest_time_averaged_scatter_"
                     f"{self.sim_name}_vs_{validation_source_name_short}{protein_string}.png")
        plt.savefig(os.path.join(plotOutDir, plot_name))

        # also generate a plotly version of the same plot:
        fig = go.Figure()

        hovertext = self.generate_hovertext(protein_df, include_STDs=True)

        fig.add_trace(
            go.Scatter(x=x_data, y=y_data, hovertext=hovertext, mode='markers',
                error_y=dict(
                    type='data',
                    symmetric=False,
                    array=yerr_upper,
                    arrayminus=yerr_lower,
                    thickness=1.5,
                    width=5,
                    color='lightseagreen'
                ),
                name='Total Counts',
                marker=dict(color='lightseagreen', size=10)
            )
        )


        # also add the free and complexed counts as separate points:
        y_data_free = np.log10(protein_df['avg_free_count'].values + 1)
        free_err_upper = free_upper - free_avg_log
        free_err_lower = free_avg_log - free_lower
        fig.add_trace(
            go.Scatter(x=x_data, y=y_data_free, hovertext=hovertext, mode='markers',
                error_y=dict(
                    type='data',
                    symmetric=False,
                    array=free_err_upper,
                    arrayminus=free_err_lower,
                    thickness=1,
                    width=2,
                    color='purple'
                ),
                name='Free Counts',
                marker=dict(color='purple', size=5)
            )
        )

        y_data_complexed = np.log10(protein_df['avg_complexed_count'].values + 1)
        complexed_err_upper = complexed_upper - complexed_avg_log
        complexed_err_lower = complexed_avg_log - complexed_lower
        fig.add_trace(
            go.Scatter(x=x_data, y=y_data_complexed, hovertext=hovertext, mode='markers',
                error_y=dict(
                    type='data',
                    symmetric=False,
                    array=complexed_err_upper,
                    arrayminus=complexed_err_lower,
                    thickness=1,
                    width=2,
                    color='orange'
                ),
                name='Complexed Counts',
                marker=dict(color='orange', size=5)
            )
        )

        # Update layout
        fig.update_layout(
            title=f"Simulation Protein Counts vs. Validation Protein Counts<br>"
                  f"Sim ID: {self.sim_name} (averaged over all total time points spanned by {self.total_cells} cells)<br> "
                  f"Validation dataset: {validation_source_name}",
            xaxis_title='Log10(Validation Counts + 1)',
            yaxis_title='Log10(Simulation Counts + 1)',
            autosize=False,
            width=800,
            height=600,
            plot_bgcolor='white',  # Set the plot area background color to white
            paper_bgcolor='white'  # Set the entire graph background to white
        )

        # add a y=x line:
        fig.add_trace(
            go.Scatter(x=[0, 6], y=[0, 6], mode="lines",
                       line=go.scatter.Line(color="black", dash="dash"),
                          opacity=0.2, name="y=x"))

        # save the figure as an html:
        plot_name = (f"proteins_of_interest_plotly_time_averaged_scatter_"
                     f"{self.sim_name}_vs_{validation_source_name_short}{protein_string}.html")
        fig.write_html(os.path.join(plotOutDir, plot_name))

    def plot_proteins_of_interest_averaged_over_generations(self, simDataFile, plotOutDir,
                                           proteins_of_interest_df,
                                           validation_source_name,
                                           validation_source_name_short):

        protein_df = self.get_sim_data_standard_deviations(
            simDataFile, proteins_of_interest_df)

        # Extract the right transformation of the error bars:
        def log_error_bars(avg, std):
            avg_log = np.log10(np.clip(avg + 1, 1, None))  # Ensure input is at least 1

            # Verify there are no negative values before calculating lower bounds
            valid_mask = (avg - std) > 0

            # Calculate lower bounds while avoiding invalid log operations
            lower = np.where(valid_mask, np.log10(np.clip(avg + 1 - std, 1, None)),
                             avg_log)  # using clip to avoid taking log(0)
            upper = np.log10(np.clip(avg + 1 + std, 1, None))

            return avg_log, lower, upper


        # Get average logs and bounds for error bars
        total_avg_log, total_lower, total_upper = log_error_bars(protein_df['avg_total_count_gen'].values,
                                                                 protein_df['std_total_count_gen'].values)
        # TODO: add free and complexed generation-averaged counts when lisener-update is pushed and edit hovertext function to have the right STDs and averages remported

        # Create plot
        plt.figure(figsize=(8, 6), dpi=300)

        # Designate the x data as the validation counts and y data as the simulation counts:
        x_data = np.log10(protein_df['RVS_count'].values + 1)
        y_data = np.log10(protein_df['avg_total_count'].values + 1)

        # Now, calculate yerr for the error bars
        yerr_upper = total_upper - total_avg_log
        yerr_lower = total_avg_log - total_lower

        # Ensure yerr does not contain negative values
        yerr_lower = np.maximum(yerr_lower, 0)  # Replace negative errors with 0

        # Combine into a 2D array with shape (2, n)
        yerr = np.array([yerr_lower, yerr_upper])

        # Use plt.errorbar to create scatter plot with error bars
        plt.scatter(x_data, y_data, color='lightblue')
        plt.errorbar(x_data, y_data, yerr=yerr, fmt='o', color='lightblue', capsize=5)

        # Plot the common name next to each point
        for i in range(len(protein_df)):
            plt.text(x_data[i] + np.log10(1), y_data[i], protein_df['common_name'].values[i],
                     fontsize=6, ha='left', va='bottom')

        # Add labels and title
        plt.title(f"Simulation Protein Counts vs. Validation Protein Counts\n"
                  f"Sim ID: {self.sim_name} (averaged over the average cell count of {self.total_cells} cells)\n "
                  f"Validation dataset: {validation_source_name}")
        plt.ylabel('Log10(Total Simulation Counts + 1)')
        plt.xlabel('Log10(Validation Counts + 1)')
        plt.grid(True)

        # Adjust y-axis to show log scale
        plt.ylim(0, max(total_avg_log) + 0.5)  # Add some space above the upper limit

        # add a y=x line:
        max_limit = max(max(x_data), max(y_data)) + 0.5
        plt.plot([0, max_limit], [0, max_limit], 'k--', lw=1)

        # obtain a string list of the common names to include in the title:
        protein_string = '_'
        for p in range(len(protein_df)):
            common_name = protein_df.iloc[p]["common_name"]
            protein_string = protein_string + '_' + str(common_name)

        # save the plot as a png:
        plot_name = (f"proteins_of_interest_generation_averaged_scatter_"
                     f"{self.sim_name}_vs_{validation_source_name_short}{protein_string}.png")
        plt.savefig(os.path.join(plotOutDir, plot_name))

        # also generate a plotly version of the same plot:
        fig = go.Figure()

        hovertext = self.generate_hovertext(protein_df, include_STDs=True, generation_averaged=True)

        fig.add_trace(
            go.Scatter(x=x_data, y=y_data, hovertext=hovertext, mode='markers',
                error_y=dict(
                    type='data',
                    symmetric=False,
                    array=yerr_upper,
                    arrayminus=yerr_lower,
                    thickness=1.5,
                    width=5,
                    color='lightseagreen'
                ),
                name='Total Counts',
                marker=dict(color='lightseagreen', size=10)
            )
        )


        # Update layout
        fig.update_layout(
            title=f"Simulation Protein Counts vs. Validation Protein Counts<br>"
                  f"Sim ID: {self.sim_name} (averaged over the average cell count of {self.total_cells} cells)<br> "
                  f"Validation dataset: {validation_source_name}",
            xaxis_title='Log10(Validation Counts + 1)',
            yaxis_title='Log10(Simulation Counts + 1)',
            autosize=False,
            width=800,
            height=600,
            plot_bgcolor='white',  # Set the plot area background color to white
            paper_bgcolor='white'  # Set the entire graph background to white
        )

        # add a y=x line:
        fig.add_trace(
            go.Scatter(x=[0, 6], y=[0, 6], mode="lines",
                       line=go.scatter.Line(color="black", dash="dash"),
                          opacity=0.2, name="y=x"))

        # save the figure as an html:
        plot_name = (f"proteins_of_interest_plotly_generation_averaged_scatter_"
                     f"{self.sim_name}_vs_{validation_source_name_short}{protein_string}.html")
        fig.write_html(os.path.join(plotOutDir, plot_name))


# Create a function that plots just the highlighted proteins of interest:
    def plot_proteins_of_interest_only(self, simDataFile, plotOutDir, validation_source_name, validation_source_name_short,
            RVS_uniprot_ids_to_monomer_IDs,
            RVS_uniprot_ids_to_schmidt_common_names,
            RVS_uniprot_ids_to_counts_dict):

        # Obtain the mapping of uniprot IDs to simulation monomer IDs for the overlapping proteins:
        RVS_uniprot_ids_to_sim_monomer_ids_dict = self.match_validation_dataset_monomer_IDs_to_simulation_monomer_IDs(
            RVS_uniprot_ids_to_monomer_IDs, RVS_uniprot_ids_to_schmidt_common_names)

        # Create a table of relevant simulation protein info for the overlapping proteins:
        RVS_sim_data_df = (
            self.create_simulation_protein_info_table(
                RVS_uniprot_ids_to_sim_monomer_ids_dict,
                RVS_uniprot_ids_to_schmidt_common_names,
                RVS_uniprot_ids_to_counts_dict))

        # If proteins of interest are defined in the simulation and validation data, plot them:
        if PROTEINS_OF_INTEREST != []:
            proteins_of_interest_df = self.map_input_proteins_to_simulation_monomer_IDs(
                PROTEINS_OF_INTEREST, RVS_sim_data_df)
            # generate the plot that averages over generations:
            generation_averaged_outpath = os.path.join(plotOutDir,
                                                              "generation_averaged_proteins_of_interest_plots")
            if not os.path.exists(generation_averaged_outpath):
                os.mkdir(generation_averaged_outpath)
            self.plot_proteins_of_interest_averaged_over_generations(simDataFile, generation_averaged_outpath,
                                                                     proteins_of_interest_df,
                                                                     validation_source_name,
                                                                     validation_source_name_short)

            # generate the plot with averages over all time points:
            time_averaged_outpath = os.path.join(plotOutDir,
                                                       "time_averaged_proteins_of_interest_plots")
            if not os.path.exists(time_averaged_outpath):
                os.mkdir(time_averaged_outpath)
            self.plot_proteins_of_interest_averaged_over_all_time(simDataFile, time_averaged_outpath,
                                               proteins_of_interest_df, validation_source_name,
                                               validation_source_name_short)

        else:
            print("HIGHLIGHT_PROTEINS_OF_INTEREST is set to TRUE but no proteins "
                  "were included in PROTEINS_OF_INTEREST, so plots highlighting "
                  "individual proteins will not be generated. ")


    def plot_validation_comparison(self, simDataFile, validationDataFile, plotOutDir, sim_name):
        # Generate sim data and get self.all_monomer_ids defined:
        avg_total_counts, avg_free_counts, avg_counts_for_monomers_in_complexes = (
            self.generate_data(simDataFile))

        # Obtain overlapping protein counts between the simulation and validation data
        self.common_names_with_multiple_monomer_ids_dict = self.check_common_name_to_gene_id_map(
            simDataFile)

        # Obtain the monomer ID to common name mapping dictionaries for later
        # use (Note: this takes in all monomer IDs and common names from rnas.tsv,
        # which may be more than what is actually in the simulation):
        self.monomer_id_to_common_name, self.common_name_to_monomer_id = self.map_monomer_ids_to_common_names()
        # TODO: make sure to check for monomer names that dont match when mapping to the sim by using the above dictionary!

        # Generate dictionaries mapping the simulation monomer IDs to common names,
        # descriptions, half lifes, average protein counts, and monomer names
        # with compartment tags:
        self.generate_sim_monomer_id_to_info_dicts(
            simDataFile, avg_total_counts, avg_free_counts, avg_counts_for_monomers_in_complexes)

        # Compare simulation data to the Schmidt et al. 2016 ST6 BW25113 proteomics data:
        (SBWST6_uniprot_IDs_to_schmidt_common_names,
         SBWST6_uniprot_IDs_to_monomer_IDs,
         SBWST6_uniprot_IDs_to_schmidt_glucose_counts) = self.get_schmidt_BW_ST6_data()

        # Compare simulation data to the Schmidt et al. 2016 ST9 MG1655 proteomics data:
        (SMGST9_uniprot_IDs_to_schmidt_common_names,
         SMGST9_uniprot_IDs_to_monomer_IDs,
         SMGST9_uniprot_IDs_to_schmidt_glucose_counts) = self.get_schmidt_BW_ST6_data()

        # Create a validation comparison output directory if it does not already exist:
        comparison_outpath = os.path.join(plotOutDir, "proteomics_validation_comparisons")
        if not os.path.exists(comparison_outpath):
            os.mkdir(comparison_outpath)

        # Create comparison plots between the simulation and the raw validation data:
        print("Plotting simulation vs raw validation data...")
        self.compare_simulation_counts_to_raw_validation_source(
            comparison_outpath, "Schmidt et al. 2016 ST6 BW25113 data",
            "Schmidt2016_ST6_BW",
            SBWST6_uniprot_IDs_to_monomer_IDs,
            SBWST6_uniprot_IDs_to_schmidt_common_names,
            SBWST6_uniprot_IDs_to_schmidt_glucose_counts)

        if INCLUDE_MANUALLY_MAPPED_PROTEINS == True:
            # Find unmapped uniprot IDs and attempt to map them manually:
            print("Plotting simulation vs raw validation data (including manually mapped proteins) ...")
            unmapped_uniprot_ids_to_simulation_monomer_ids, unmapped_uniprot_IDs_to_sim_common_names = self.find_monomer_ids_for_unmapped_uniprot_ids(SBWST6_uniprot_IDs_to_schmidt_common_names,
             SBWST6_uniprot_IDs_to_monomer_IDs)

            self.compare_simulation_counts_to_raw_validation_source_with_manual_mappings(
                comparison_outpath, "Schmidt et al. 2016 ST6 BW25113 data",
                "Schmidt2016_ST6_BW",
                SBWST6_uniprot_IDs_to_monomer_IDs,
                SBWST6_uniprot_IDs_to_schmidt_common_names,
                SBWST6_uniprot_IDs_to_schmidt_glucose_counts,
                unmapped_uniprot_ids_to_simulation_monomer_ids)

        if HIGHLIGHT_PROTEINS_OF_INTEREST == True:
            # Plot the the proteins of interest on the plot too:
            print("Plotting proteins of interest...")
            # Generate a file path for the plots:
            proteins_of_interest_plot_outpaths = os.path.join(comparison_outpath, "comparisons_with_proteins_of_interest_only")
            if not os.path.exists(proteins_of_interest_plot_outpaths):
                os.mkdir(proteins_of_interest_plot_outpaths)

            # Create comparison plots between the simulation and the raw validation data:
            self.compare_simulation_counts_to_raw_validation_source(
                comparison_outpath, "Schmidt et al. 2016 ST6 BW25113 data",
                "Schmidt2016_ST6_BW",
                SBWST6_uniprot_IDs_to_monomer_IDs,
                SBWST6_uniprot_IDs_to_schmidt_common_names,
                SBWST6_uniprot_IDs_to_schmidt_glucose_counts,
                plot_proteins_of_interest=True)

            # Create the comparison plots including manually mapped proteins:
            unmapped_uniprot_ids_to_simulation_monomer_ids, unmapped_uniprot_IDs_to_sim_common_names = self.find_monomer_ids_for_unmapped_uniprot_ids(
                SBWST6_uniprot_IDs_to_schmidt_common_names,
                SBWST6_uniprot_IDs_to_monomer_IDs)
            self.compare_simulation_counts_to_raw_validation_source_with_manual_mappings(
                comparison_outpath, "Schmidt et al. 2016 ST6 BW25113 data",
                "Schmidt2016_ST6_BW",
                SBWST6_uniprot_IDs_to_monomer_IDs,
                SBWST6_uniprot_IDs_to_schmidt_common_names,
                SBWST6_uniprot_IDs_to_schmidt_glucose_counts,
                unmapped_uniprot_ids_to_simulation_monomer_ids,
                plot_proteins_of_interest=True)

            # Generate plots of just the proteins of interest:
            self.plot_proteins_of_interest_only(simDataFile, proteins_of_interest_plot_outpaths, "Schmidt et al. 2016 ST6 BW25113 data",
                "Schmidt2016_ST6_BW",
                SBWST6_uniprot_IDs_to_monomer_IDs,
                SBWST6_uniprot_IDs_to_schmidt_common_names,
                SBWST6_uniprot_IDs_to_schmidt_glucose_counts)

    # TODO: add validation source to validation source comparison options
    # TODO: decide if it makes sense to have validation data compare against all validation data available or just the data that matches the simulation data?
        if PLOT_VALIDATION_SOURCE_COMPARISONS == True:
            print("Plotting validation source comparisons...")
            # generate a validation dataset outpath:
            validation_outpath = os.path.join(comparison_outpath, "validation_source_comparisons")
            if not os.path.exists(validation_outpath):
                os.mkdir(validation_outpath)
            # Compare the validation data generated from the sim (Schmidt et al. 2015 BW25113) to the Schmidt et al. 2016 ST9 MG1655 data:
            self.compare_validation_datasets(validation_outpath,
                                                                  "Schmidt et al. 2015 BW25113 Glucose data from saved simulation",
                                                                  "Schmidt2015_ST6_BW",
                                             SBWST6_uniprot_IDs_to_schmidt_common_names,
                                             SBWST6_uniprot_IDs_to_monomer_IDs,
                                             SBWST6_uniprot_IDs_to_schmidt_glucose_counts,
                                                             "Schmidt et al. 2016 ST9 MG1655 Glucose data",
                                                                  "Schmidt2016_ST9_MG",
                                             SMGST9_uniprot_IDs_to_schmidt_common_names,
                                             SMGST9_uniprot_IDs_to_monomer_IDs,
                                             SMGST9_uniprot_IDs_to_schmidt_glucose_counts)

            # plot the data saved with the sim
            self.compare_raw_validation_dataset_to_saved_simulation_validation_dataset(simDataFile, validationDataFile, validation_outpath, "Schmidt et al. 2016 ST9 MG1655 Glucose data",
                                                                  "Schmidt2015_ST6_BW", SBWST6_uniprot_IDs_to_schmidt_common_names,
                                             SBWST6_uniprot_IDs_to_monomer_IDs,
                                             SBWST6_uniprot_IDs_to_schmidt_glucose_counts)


    def do_plot(self, variantDir, plotOutDir, plotOutFileName, simDataFile,
                validationDataFile, metadata):

        # Generate the data for the simulation:
        self.sim_name = metadata["description"]
        self.n_total_gens = self.ap.n_generation
        self.plot_validation_comparison(simDataFile, validationDataFile, plotOutDir, self.sim_name)



if __name__ == '__main__':
    Plot().cli()

