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
import io
from wholecell.io import tsv
from wholecell.utils.filepath import ROOT_PATH


""" USER INPUTS """

# Indicate the number of generations to be ignored at the start of each seed:
IGNORE_FIRST_N_GENS = 2 # 2 for local, 14 for Sherlock (w/ 24 total gens)

""" END USER INPUTS """


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

        # # Calculate the average total monomer counts across all cells:
        # total_protein_counts = (read_stacked_columns(all_cells,
        # 	'MonomerCounts', 'monomerCounts')).mean(axis=0)
        #
        # # Make into a np.array:
        # total_protein_counts = np.array(total_protein_counts)
        # todo: delete the above if it matches below

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


    # extract the validation data from Schmidt et al. 2016 supplementary table 9:
    # two ways to do this: make a validation data set that is literally the same as the schmidt set by going
    def get_schmidt_MG_validation_data_from_ST9(self):
        """
        Extracts the protein counts from Schmidt et al. 2016 supplementary
        table 9
        Returns: a dictionary mapping gene symbols to their average protein
        counts
        """
        schmidt_st9_file = os.path.join(
            ROOT_PATH, 'reconstruction', 'ecoli',
                                     'flat', 'Schmidt_2016_ST9.csv')
        gene_symbol_to_avg_counts = {}
        with io.open(schmidt_st9_file, 'r') as f:
            reader = csv.reader(f, delimiter=',')
            headers = next(reader)
            gene_symbol_index = headers.index('Gene')
            avg_counts_index = headers.index('Copies/Cell_MG1655.Glucose')

            for line in reader:
                gene_symbol = line[gene_symbol_index]
                avg_counts = line[avg_counts_index]
                gene_symbol_to_avg_counts[gene_symbol] = avg_counts

        # Clean ST9_dict: convert counts to int and remove bogus entries
        for gene, count in list(gene_symbol_to_avg_counts.items()):
            if count == '' or count is None:
                del gene_symbol_to_avg_counts[gene]  # Remove bogus entries
            else:
                gene_symbol_to_avg_counts[gene] = int(count)  # Convert counts to integers

        return gene_symbol_to_avg_counts

    # get the ST9 BW counts dictonary:
    def get_schmidt_BW_validation_data_from_ST9(self):
        """
        Extracts the protein counts from Schmidt et al. 2016 supplementary
        table 9 data for BW25
        Returns: a dictionary mapping gene symbols to their average protein
        counts
        """
        schmidt_st9_file = os.path.join(
            ROOT_PATH, 'reconstruction', 'ecoli',
                                     'flat', 'Schmidt_2016_ST9.csv')
        gene_symbol_to_avg_counts = {}
        with io.open(schmidt_st9_file, 'r') as f:
            reader = csv.reader(f, delimiter=',')
            headers = next(reader)
            gene_symbol_index = headers.index('Gene')
            avg_counts_index = headers.index('Copies/Cell_BW25113.Glucose')

            for line in reader:
                gene_symbol = line[gene_symbol_index]
                avg_counts = line[avg_counts_index]
                gene_symbol_to_avg_counts[gene_symbol] = avg_counts

        # Clean ST9_dict: convert counts to int and remove bogus entries
        for gene, count in list(gene_symbol_to_avg_counts.items()):
            if count == '' or count is None:
                del gene_symbol_to_avg_counts[gene]  # Remove bogus entries
            else:
                gene_symbol_to_avg_counts[gene] = int(count)  # Convert counts to integers

        return gene_symbol_to_avg_counts


    # get the ST6 BW counts dictonary:
    def get_schmidt_BW_validation_data_from_ST6(self):
        """
        Extracts the protein counts from Schmidt et al. 2016 supplementary
        table 6 data for BW25
        Returns: a dictionary mapping gene symbols to their average protein
        counts
        """
        schmidt_st6_file = os.path.join(
            ROOT_PATH, 'validation', 'ecoli',
            'flat', 'schmidt2015_javier_table.tsv')
        gene_symbol_to_avg_counts = {}
        with io.open(schmidt_st6_file, 'r') as f:
            reader = csv.reader(f, delimiter='\t')
            headers = next(reader)
            gene_symbol_index = headers.index('GeneName')
            avg_counts_index = headers.index('Glucose')

            for line in reader:
                gene_symbol = line[gene_symbol_index]
                avg_counts = line[avg_counts_index]
                gene_symbol_to_avg_counts[gene_symbol] = avg_counts

        # Clean ST9_dict: convert counts to int and remove bogus entries
        for gene, count in list(gene_symbol_to_avg_counts.items()):
            if count == '' or count is None:
                del gene_symbol_to_avg_counts[gene]  # Remove bogus entries
            else:
                gene_symbol_to_avg_counts[gene] = int(count)  # Convert counts to integers

        return gene_symbol_to_avg_counts

    # Obtain a dictionary of the common names associated with multiple monomer IDs:
    def check_common_name_to_gene_id_map(self, simDataFile):
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
        common_names_with_multiple_monomer_ids_narrowed_dict = {}
        for common_name in common_names_with_multiple_monomer_ids_dict.keys():
            monomer_ids = common_names_with_multiple_monomer_ids_dict[common_name]
            monomer_ids_in_sim = []
            print("Checking which monomer IDs for", common_name, "are in the simulation...")
            for monomer_id in monomer_ids:
                monomer_ID = self.check_validity_and_get_compartment(simDataFile, [monomer_id])[0]
                if monomer_ID in self.all_monomer_ids:
                    monomer_ids_in_sim.append(monomer_ID)
                    print(monomer_ID, "is in the simulation.")
            # only append the monomer IDs that are actually in the simulation:
            common_names_with_multiple_monomer_ids_narrowed_dict[common_name] = monomer_ids_in_sim

        # make it official:
        self.common_names_with_multiple_monomer_ids_narrowed_dict = common_names_with_multiple_monomer_ids_narrowed_dict

        # return a dictionary of these values in case it is helpful later:
        return common_names_with_multiple_monomer_ids_dict

    # create a function that checks if something is a vaild molecule in the sim:
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

    # Create a dictionary mapping the monomer_IDs to common names (and vice versa):
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
        common_name_to_monomer_id = {}
        monomer_id_to_common_name = {}
        for i in range(len(rnas_data)):
            common_name = rnas_data.iloc[i]['common_name']
            monomer_ids = rnas_data.iloc[i]['monomer_ids']
            # if only one monomer ID maps to the common name, map it directly:
            if len(monomer_ids) == 1:
                common_name_to_monomer_id[common_name] = monomer_ids
                monomer_id_to_common_name[monomer_ids] = common_name

            # generate a warning if a common name maps to multiple monomer IDs:
            if len(monomer_ids) > 1:
                print("WARNING: ", common_name,
                      "maps to multiple monomer IDs:", monomer_ids)
                # find the monomer ID that is actually in the simulation:
                monomer_IDs_in_sim = self.common_names_with_multiple_monomer_ids_narrowed_dict[common_name]
                if len(monomer_IDs_in_sim) == 1:
                    print("This can mess up later mappings! Mapping to the "
                          "monomer ID that appears in the simulation by default:",
                          monomer_IDs_in_sim)
                    # map to first monomer ID by default
                    common_name_to_monomer_id[common_name] = monomer_IDs_in_sim
                if len(monomer_IDs_in_sim) > 1:
                    print("Multiple monomer IDs for", common_name,
                          "are in the simulation:", monomer_IDs_in_sim,
                          "Mapping to the first one by default:", monomer_IDs_in_sim[0])
                    common_name_to_monomer_id[common_name] = [monomer_IDs_in_sim[0]]

        return monomer_id_to_common_name, common_name_to_monomer_id




    def plot_validation_comparison(self, simDataFile, validationDataFile,plotOutDir, sim_name):
        # run this so that self.all_monomer_ids is defined:
        avg_total_counts, avg_free_counts, avg_counts_for_monomers_in_complexes = self.generate_data(
            simDataFile)

        # obtain overlapping protein counts between the simulation and validation data
        self.common_names_with_multiple_monomer_ids_dict = self.check_common_name_to_gene_id_map(simDataFile)


        self.monomer_id_to_common_name, self.common_name_to_monomer_id = self.map_monomer_ids_to_common_names()
        # TODO: make sure to check for monomer names that dont match when mapping to the sim by using the above dictionary!

        hi = 5







    def do_plot(self, variantDir, plotOutDir, plotOutFileName, simDataFile,
                validationDataFile, metadata):

        # Generate the data for the simulation:
        self.sim_name = metadata["description"]
        self.n_total_gens = self.ap.n_generation
        self.plot_validation_comparison(simDataFile, validationDataFile, plotOutDir, self.sim_name)



if __name__ == '__main__':
    Plot().cli()
