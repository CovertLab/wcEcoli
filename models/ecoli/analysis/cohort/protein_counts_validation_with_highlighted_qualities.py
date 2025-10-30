
import pickle
import os
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

# input proteins to highlight:
HIGHLIGHT_PROTEINS =[] #["TRYPSYN-APROTEIN[c]", 'TRYPSYN-BPROTEIN[c]', 'ANTHRANSYNCOMPI-MONOMER[c]', 'CYCLASE-MONOMER[c]', 'EG11274-MONOMER[m]', 'PRAI-IGPS[c]', 'ASPKINIHOMOSERDEHYDROGI-MONOMER[c]',
					  #'ATPPHOSRIBOSTRANS-MONOMER[c]', 'GLUTAMIDOTRANS-MONOMER[c]', 'HISTCYCLOPRATPPHOS[c]', 'HISTDEHYD-MONOMER[c]', 'HISTPHOSTRANS-MONOMER[c]','HOMOSERKIN-MONOMER[c]',
					  #'IMIDPHOSPHADEHYDHISTIDPHOSPHA-MONOMER[c]', 'PRIBFAICARPISOM-MONOMER[c]', 'THRESYN-MONOMER[c]']



# threshold for complex fraction:
COMPLEX_FRACTION_THRESHOLD = 0.9

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

	# make a function that makes a table with the protein IDs, gene symbols, descriptive names, schmidt validation data, and schmidt validation data from ST9:
	def make_schmidt_validation_table(self, simDataFile, validationDataFile):
		# Get the dictionary of gene symbols to avg counts from ST9:
		ST9_dict = self.get_schmidt_BW_validation_data_from_ST9()

		# Get the validation data:
		(sim_schmidt_counts, val_schmidt_counts,
		schmidt_overlap_ids, sim_wisniewski_counts,
		val_wisniewski_counts, wisniewski_overlap_ids) = self.get_validation_data(simDataFile, validationDataFile)

		# Make a dictionary of Schmidt ST6 gene_symbols to protein IDs:
		ST6_dict = {}

		# Create a dictonary of each ST6 gene Id to its descriptive name:
		ST6_descriptive_dict = {}

		# Create a dictonary of each ST6 gene ID to its simulation counts:
		ST6_sim_counts_dict = {}

		# Create a dictonary of each ST6 gene ID to the validation counts:
		ST6_val_counts_dict = {}

		# Create all dictionaries:
		#NOTE: this is not directly taking from ST6, rather it is taking from the validationDataFile proteins that "overlap" with schmidt data.
		for idx, protein_id in enumerate(schmidt_overlap_ids):
			gene_symbol = self.get_common_name(protein_id)
			ST6_dict[gene_symbol] = protein_id

			descriptive_name = self.get_descriptive_name(protein_id)
			ST6_descriptive_dict[gene_symbol] = descriptive_name

			sim_counts = sim_schmidt_counts[idx]
			ST6_sim_counts_dict[gene_symbol] = sim_counts

			val_counts = val_schmidt_counts[idx]
			ST6_val_counts_dict[gene_symbol] = val_counts

		# Create a DataFrame to hold the combined data
		SDF = []

		for gene in ST6_dict.keys():
			if gene in ST9_dict:
				row = {
					'Gene_Symbol': gene,
					'Protein_ID': ST6_dict[gene],
					'Descriptive_Name': ST6_descriptive_dict[gene],
					'Sim_Counts': ST6_sim_counts_dict[gene],
					'Val_BW_Counts': ST6_val_counts_dict[gene],
					'Val_MG_Counts': ST9_dict[gene]
				}
				SDF.append(row)

		# Convert the combined data into a DataFrame
		SDF_df = pd.DataFrame(SDF)

		return SDF_df

	# make a Schmit comparison data for validation file schmidt data vs ST6 BW data:
	def make_schmidt_validation_table_comparing_val_BW_to_direct_ST6(self, simDataFile, validationDataFile):
		# Get the dictionary of gene symbols to avg counts from ST9:
		ST6_dict = self.get_schmidt_BW_validation_data_from_ST6()

		# Get the validation data:
		(sim_schmidt_counts, val_schmidt_counts,
		schmidt_overlap_ids, sim_wisniewski_counts,
		val_wisniewski_counts, wisniewski_overlap_ids) = self.get_validation_data(simDataFile, validationDataFile)

		# Make a dictionary of Schmidt ST6 gene_symbols to protein IDs:
		val_dict = {}

		# Create a dictonary of each ST6 gene Id to its descriptive name:
		val_descriptive_dict = {}

		# Create a dictonary of each ST6 gene ID to its simulation counts:
		val_sim_counts_dict = {}

		# Create a dictonary of each ST6 gene ID to the validation counts:
		val_val_counts_dict = {}

		# Create all dictionaries:
		#NOTE: this is not directly taking from ST6, rather it is taking from the validationDataFile proteins that "overlap" with schmidt data.
		for idx, protein_id in enumerate(schmidt_overlap_ids):
			gene_symbol = self.get_common_name(protein_id)
			val_dict[gene_symbol] = protein_id

			descriptive_name = self.get_descriptive_name(protein_id)
			val_descriptive_dict[gene_symbol] = descriptive_name

			sim_counts = sim_schmidt_counts[idx]
			val_sim_counts_dict[gene_symbol] = sim_counts

			val_counts = val_schmidt_counts[idx]
			val_val_counts_dict[gene_symbol] = val_counts

		# Create a DataFrame to hold the combined data
		SDF = []

		for gene in val_dict.keys():
			if gene in ST6_dict:
				row = {
					'Gene_Symbol': gene,
					'Protein_ID': val_dict[gene],
					'Descriptive_Name': val_descriptive_dict[gene],
					'Sim_Counts': val_sim_counts_dict[gene],
					'Val_BW_Counts': val_val_counts_dict[gene],
					'ST6_BW_Counts': ST6_dict[gene]
				}
				SDF.append(row)

		# Convert the combined data into a DataFrame
		SDF_df = pd.DataFrame(SDF)

		return SDF_df

	# make comparison between ST9 and validatoin file data:
	def make_schmidt_validation_table_comparing_val_BW_to_direct_ST9(self, simDataFile, validationDataFile):
		# Get the dictionary of gene symbols to avg counts from ST9:
		ST9_dict = self.get_schmidt_BW_validation_data_from_ST9()

		# Get the validation data:
		(sim_schmidt_counts, val_schmidt_counts,
		schmidt_overlap_ids, sim_wisniewski_counts,
		val_wisniewski_counts, wisniewski_overlap_ids) = self.get_validation_data(simDataFile, validationDataFile)

		# Make a dictionary of Schmidt ST6 gene_symbols to protein IDs:
		val_dict = {}

		# Create a dictonary of each ST6 gene Id to its descriptive name:
		val_descriptive_dict = {}

		# Create a dictonary of each ST6 gene ID to its simulation counts:
		val_sim_counts_dict = {}

		# Create a dictonary of each ST6 gene ID to the validation counts:
		val_val_counts_dict = {}

		# Create all dictionaries:
		#NOTE: this is not directly taking from ST6, rather it is taking from the validationDataFile proteins that "overlap" with schmidt data.
		for idx, protein_id in enumerate(schmidt_overlap_ids):
			gene_symbol = self.get_common_name(protein_id)
			val_dict[gene_symbol] = protein_id

			descriptive_name = self.get_descriptive_name(protein_id)
			val_descriptive_dict[gene_symbol] = descriptive_name

			sim_counts = sim_schmidt_counts[idx]
			val_sim_counts_dict[gene_symbol] = sim_counts

			val_counts = val_schmidt_counts[idx]
			val_val_counts_dict[gene_symbol] = val_counts

		# Create a DataFrame to hold the combined data
		SDF = []

		for gene in val_dict.keys():
			if gene in ST9_dict:
				row = {
					'Gene_Symbol': gene,
					'Protein_ID': val_dict[gene],
					'Descriptive_Name': val_descriptive_dict[gene],
					'Sim_Counts': val_sim_counts_dict[gene],
					'Val_BW_Counts': val_val_counts_dict[gene],
					'ST9_BW_Counts': ST9_dict[gene]
				}
				SDF.append(row)

		# Convert the combined data into a DataFrame
		SDF_df = pd.DataFrame(SDF)

		return SDF_df





	# hi = 5
	# function to match gene symbols to monomer ids
	def get_gene_symbols_for_monomer_ids(self):
		"""
		Extracts the gene symbols for each monomer id in the model.

		Returns: a dictionary mapping monomer ids to gene symbols.
		Code adapted from convert_to_flat.py.
		"""
		RNAS_FILE = os.path.join(ROOT_PATH, 'reconstruction', 'ecoli',
									 'flat', 'rnas.tsv')
		with (io.open(RNAS_FILE, 'rb') as f):
			reader = tsv.reader(f, delimiter='\t')
			headers = next(reader)
			while headers[0].startswith('#'):
				headers = next(reader)

			# extract relevant information
			gene_symbol_index = headers.index('common_name')
			protein_id_index = headers.index('monomer_ids')
			monomer_ids_to_gene_symbols = {}
			for line in reader:
				gene_symbol = line[gene_symbol_index]
				protein_id = list(
						line[protein_id_index][2:-2].split('", "'))[0]  # not sure what this does
				monomer_ids_to_gene_symbols[protein_id] = gene_symbol

		return monomer_ids_to_gene_symbols

	def get_descriptive_name_for_monomer_ids(self):
		PROTEINS_FILE = os.path.join(ROOT_PATH, 'reconstruction', 'ecoli',
										 'flat', 'proteins.tsv')
		with (io.open(PROTEINS_FILE, 'rb') as f):
			reader = tsv.reader(f, delimiter='\t')
			headers = next(reader)
			while headers[0].startswith('#'):
				headers = next(reader)

			# extract relevant information
			gene_description_index = headers.index('common_name')
			protein_id_index = headers.index('id')
			monomer_ids_to_descriptive_names = {}
			for line in reader:
				gene_description = line[gene_description_index]
				protein_id = line[protein_id_index]
				monomer_ids_to_descriptive_names[protein_id] = gene_description

		return monomer_ids_to_descriptive_names

	def get_common_name(self, protein_id):
		"""
		Obtains the common names for each protein of interest
		Args:
		    protein_id: the name of the protein(s) of interest

		Returns: the common name for the protein(s) of interest
		"""
		if protein_id == 'NG-GFP-MONOMER[c]':
			return 'GFP'

		else:
			if "[" in protein_id:
				protein = protein_id[:-3]  # subtract the compartment
			else:
				protein = protein_id
			common_name = self.get_gene_symbols_for_monomer_ids()[protein]

		return common_name

	def get_descriptive_name(self, protein_id):
		if protein_id == 'NG-GFP-MONOMER[c]':
			return 'GFP'

		else:
			if "[" in protein_id:
				protein = protein_id[:-3]  # subtract the compartment
			else:
				protein = protein_id
			descriptive_name = self.get_descriptive_name_for_monomer_ids()[protein]

		return descriptive_name

	def get_ids(self, monomer_idx_dict, protein_idxs):
		"""
        Obtain the protein ids for each protein based on its respective index
        Args:
            monomer_idx_dict: a dictionary that maps protein names to their idx
            protein_idxs: an array of indices for each protein, respectively

        Returns: the corresponding id for each respective protein index
        """
		inv_monomer_idx_dict = {idx: i for i, idx in monomer_idx_dict.items()}
		protein_ids = [inv_monomer_idx_dict.get(monomer_id) for
					   monomer_id in protein_idxs]

		return protein_ids

	def get_half_lives(self, simDataFile):
		"""
		Extract the half lives for each protein in the simulation
		Args:
			simDataFile: simulation data file

		Returns: a dictionary of the half lives for each protein
		"""
		with open(simDataFile, 'rb') as f:
			sim_data = pickle.load(f)

		# Extract half lives for each protein:
		half_lives = sim_data.process.translation.monomer_data["deg_rate"]

		# Create a dictionary of the half lives for each protein:
		total_protein_counts, _, _ = self.generate_data(simDataFile)
		monomer_to_half_life = dict(zip(self.all_monomer_ids, half_lives))

		return monomer_to_half_life

	def get_LogData(self, protein_idxs, interest_protein_counts, index_vals=[]):
		"""
        Covert monomer count data to a log10 scale
        Args:
            protein_idxs: an array of the indices for proteins to be plotted
             (this should be smaller than interest_protein_counts if the data
              has been filtered)
            interest_protein_counts: the full data structure of all proteins
            and their respective counts (no filter applied)
            (usually size variants by # of proteins), either filtered or not
            index_vals: if the protein idxs are not in sequential order (usually
            happens after filtering the data), include a list of the original
            indices for each protein in this variable.

        Returns: the log10 values for the total average counts for each protein
        """
		avg_log_interest_proteins = np.zeros(len(protein_idxs))
		for idx in range(len(protein_idxs)):
			if len(index_vals) == 0:
				index = idx
			else:
				index = index_vals[idx]
			avg_log_interest_proteins[idx] = (
				np.log10(interest_protein_counts[index] + 1))
		return avg_log_interest_proteins


	def transpose_and_reshape(self, data):
		"""
		Transpose and reshape the data to obtain the desired format for saving
		Args:
			data: data to be transposed and reshaped

		Returns: the transposed and reshaped data
		"""
		data = np.transpose(np.array(data)); data = data.reshape(-1, 1)
		return data


	def generate_validation_plotly(self, simDataFile, plotOutDir, simulationCounts, validationCounts,
								 overlapIDs, sim_name, val_name):
		fig = go.Figure()

		# Compute log10 values
		x = self.get_LogData(overlapIDs, validationCounts)
		y = self.get_LogData(overlapIDs, simulationCounts)

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

		# compute the rsquared value:
		r_squared_30_above = r2_score(x_above_30, y_above_30)

		# have the monomer IDs be the overlap text
		monomer_counts_table, monomer_id_to_complex_fraction, monomer_id_to_complex_counts = self.determine_fraction_table(
			simDataFile)
		monomer_to_half_life = self.get_half_lives(simDataFile)

		# retreive the common names and the descriptions:
		common_names = [self.get_common_name(protein_id) for protein_id in overlapIDs]
		descriptive_names = [self.get_descriptive_name(protein_id) for protein_id in overlapIDs]

		# create a dataframe of the protein ids and their half life sources:
		protein_df = pd.DataFrame({"protein_id": overlapIDs,
								   "common_name": common_names,
								   "descriptive_name": descriptive_names,
								   'simulation_protein_counts': y,
								   'validation_protein_counts': x})
		protein_df['fraction_in_complex'] = protein_df['protein_id'].map(
			monomer_id_to_complex_fraction)
		protein_df['complex_counts'] = protein_df['protein_id'].map(monomer_id_to_complex_counts)
		protein_df['half_life'] = protein_df['protein_id'].map(monomer_to_half_life)
		hovertext = self.hover_text_info(protein_df)

		# Add scatter trace
		fig.add_trace(
			go.Scatter(x=x, y=y, hovertext=hovertext, mode='markers',
					   name=f"Counts (R^2 for counts > 30:{ round(r_squared_30_above,3)})"))


		# Add trendline trace
		fig.add_trace(
			go.Scatter(x=x, y=trendline_y, mode='lines',
					   name=f'Linear fit: {p}',
					   line=dict(color='green')))
		fig.add_trace(
			go.Scatter(x=x, y=trendline_y_above_30, mode='lines',
					   name=f'Linear fit (counts > 30): {p_above_30}',
					   line=dict(color='pink')))

		# Update layout
		fig.update_traces(marker_size=3)
		fig.update_layout(
			title=f"Simulation Protein Counts ({sim_name}) "
				  f"vs. Validation Protein Counts ({val_name} et al.)",
			xaxis_title="log10(Validation Protein Counts)",
			yaxis_title=f"log10(Simulation Protein Counts)",
			autosize=False, width=900, height=600)

		# add a y=x line
		fig.add_trace(
			go.Scatter(x=[0, 6], y=[0, 6], mode="lines",
					line=go.scatter.Line(color="black", dash="dash"),
					opacity=0.2, name="y=x"));

		# save the figure as an html:
		plot_name = f"proteinCountsValidation_cohortPlot_{sim_name}_vs_{val_name}.html"
		fig.write_html(os.path.join(plotOutDir, plot_name))


	def generate_validation_plot(self, plotOutDir, simulationCounts, validationCounts,
								 overlapIDs, sim_name, val_name):
		plt.figure(figsize=(9, 6), dpi=200)

		# Compute log10 values
		x = self.get_LogData(overlapIDs, validationCounts)
		y = self.get_LogData(overlapIDs, simulationCounts)

		# Compute linear trendline for counts above log10(30):
		above_30_idx = np.where((x > np.log10(30 + 1)) & (y > np.log10(30 + 1)))
		x_above_30 = x[above_30_idx]
		y_above_30 = y[above_30_idx]
		z_above_30 = np.polyfit(x_above_30, y_above_30, 1)
		p_above_30 = np.poly1d(z_above_30)
		trendline_y_above_30 = p_above_30(x)

		# compute the pearsonr and R2:
		r_value, p_val = pearsonr(x_above_30, y_above_30)
		pr2 = r_value ** 2
		# compute the coefficent of determination rsquared value:
		r_squared_30_above = r2_score(x_above_30, y_above_30)

		# Add scatter trace
		plt.scatter(x=x, y=y, s=5, label=f"Counts (n={len(x)})",
					alpha=0.5, color='lightseagreen')

		# Compute linear trendline
		z = np.polyfit(x, y, 1)
		p = np.poly1d(z)
		trendline_y = p(x)

		# Add trendline trace
		plt.plot(x, trendline_y,
					   label=f'Linear fit: {p}',
					   color='orange')
		plt.plot(x,trendline_y_above_30,
					   label=f'Linear fit (counts > 30): {p_above_30}',
					   color='pink')
		# add a y=x line
		plt.plot([0, 6], [0, 6],
				color="black", linestyle="dashed",
				 alpha=0.2, label="y=x");


		# Update layout

		plt.title(f"Simulation Protein Counts ({sim_name}) "
				  f"vs.\n Validation Protein Counts ({val_name} et al.)")
		plt.xlabel("log10(Validation Protein Counts)")
		plt.ylabel(f"log10(Simulation Protein Counts)")
		plt.xlim(0, 6)
		plt.ylim(0, 6)
		plt.gca().set_aspect('equal', adjustable='box')

		plt.legend(bbox_to_anchor=(1.45,1.), loc='upper right', fontsize=8, )
		plt.text(0.96, 0.03, f'Pearson R (counts > 30): {round(r_value,3)}\n Pearson $R^2$ (counts > 30): {round(pr2,3)}\nCoefficent of determination $R^2$ (counts > 30): {round(r_squared_30_above,3)}',
				 ha='right', va='bottom',
				 transform=plt.gca().transAxes,
				 fontsize=6, color='gray')



		# save the figure as an html:
		plot_name = f"proteinCountsValidation_cohortPlot_{sim_name}_vs_{val_name}_matplotlib.pdf"
		plt.savefig(os.path.join(plotOutDir, plot_name))

	def determine_fraction_table(self, simDataFile):
		"""
		Determines the fraction of proteins in complex form for each monomer
		Args:
			simDataFile: simulation data file
		Returns:
			complex_fraction: fraction of proteins in complex form for each monomer
		"""
		# Obtain counts data:
		avg_total_counts, avg_free_counts, avg_complex_counts = self.generate_data(
			simDataFile)
		# Remove the last three characters from each value:
		#monomer_ids = [id[:-3] for id in self.all_monomer_ids]
		monomer_ids = self.all_monomer_ids

		# Calculate the faction in complex and in free from:
		complex_fraction = np.zeros(len(monomer_ids))
		free_fraction = np.zeros(len(monomer_ids))
		for i in range(len(monomer_ids)):
			if avg_total_counts[i] == 0:
				complex_fraction[i] = 0 # set to zero so it is still accounted for
				free_fraction[i] = 0
				print('Monomer ID:', monomer_ids[i], 'has no total counts')
			else:
				complex_fraction[i] = avg_complex_counts[i] / avg_total_counts[i]
				free_fraction[i] = avg_free_counts[i] / avg_total_counts[i]


		# make a table of the monomers and their counts and their fractions:
		monomer_counts_table = pd.DataFrame({'Monomer ID': monomer_ids,
									   'Total Counts': avg_total_counts,
									   'Free Counts': avg_free_counts,
									   'Free Fraction': free_fraction,
									   'Complex Counts': avg_complex_counts,
									   'Complex Fraction': complex_fraction})

		monomer_id_to_complex_fraction = dict(zip(monomer_ids, complex_fraction))
		monomer_id_to_complex_counts = dict(zip(monomer_ids, avg_complex_counts))

		return monomer_counts_table, monomer_id_to_complex_fraction, monomer_id_to_complex_counts

	def hover_text_info(self, dataframe):
		hovertext = dataframe.apply(lambda
					row: f"Monomer ID: {row['protein_id']}<br>Common Name: {row['common_name']}<br>Description: {row['descriptive_name']}<br>HL Value: {row['half_life']}<br>validation AMC: {10 ** (row['validation_protein_counts'])}<br>Simulation AMC: {10 ** (row['simulation_protein_counts'])}<br>Avg. Complexed Monomer Counts: {row['complex_counts']} Complexed Fraction: {row['fraction_in_complex']}<br>",
									axis=1)
		return hovertext
	def plot_by_complex_fraction_plotly(self, simDataFile, plotOutDir, simulationCounts, validationCounts,
								 overlapIDs, sim_name, val_name):
		"""
		Extract the complex count fractions for each protein in the simulation
		Args:
			simDataFile: simulation data file

		Returns: a dictionary of the complex count fractions for each protein
		"""


		monomer_counts_table, monomer_id_to_complex_fraction, monomer_id_to_complex_counts = self.determine_fraction_table(simDataFile)

		# Compute log10 values
		x = self.get_LogData(overlapIDs, validationCounts)
		y = self.get_LogData(overlapIDs, simulationCounts)

		monomer_to_half_life = self.get_half_lives(simDataFile)

		# retreive the common names and the descriptions:
		common_names = [self.get_common_name(protein_id) for protein_id in overlapIDs]
		descriptive_names = [self.get_descriptive_name(protein_id) for protein_id in overlapIDs]

		# create a dataframe of the protein ids and their half life sources:
		protein_df = pd.DataFrame({"protein_id": overlapIDs,
								   "common_name": common_names,
								   "descriptive_name": descriptive_names,
								   'simulation_protein_counts': y,
								   'validation_protein_counts': x})
		protein_df['fraction_in_complex'] = protein_df['protein_id'].map(monomer_id_to_complex_fraction)
		protein_df['complex_counts'] = protein_df['protein_id'].map(monomer_id_to_complex_counts)
		protein_df['half_life'] = protein_df['protein_id'].map(monomer_to_half_life)

		# make a yes or no column for whether the protein is in complex form:
		protein_df['in_complex'] = protein_df['fraction_in_complex'].apply(lambda x: 'Yes' if x > 0 else 'No')


		# split up the data frames by complex fraction:
		complex_fraction_types = protein_df['in_complex'].unique()
		complex_fraction_dfs = {}
		for type in complex_fraction_types:
			complex_fraction_dfs[type] = protein_df[protein_df['in_complex'] == type]

		# create a color map for each complex fraction:
		color_map = {
			"Yes": "lightseagreen",
			"No": "yellowgreen"}

		name_map = {
			"Yes": "Has Complexed Counts",
			"No": "Free Monomer Only"}

		# create a scatter plot for each half life source:
		plt.figure(figsize=(9, 6), dpi=200)
		fig = go.Figure()
		for type, df in complex_fraction_dfs.items():
			print(type)
			c = color_map[type]; name = name_map[type]
			hovertext = self.hover_text_info(df)
			fig.add_trace(go.Scatter(x=df['validation_protein_counts'], y=df['simulation_protein_counts'], hovertext=hovertext, mode='markers',
									 name=f"{name} (n={len(df)})",
									 marker=dict(color=c, size=.7, opacity=.5)))

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

		# compute the rsquared value:
		r_squared_30_above = r2_score(x_above_30, y_above_30)

		# Add trendline trace
		fig.add_trace(
			go.Scatter(x=x, y=trendline_y, mode='lines',
					   name=f'Linear fit: {p}',
					   line=dict(color='green')))
		fig.add_trace(
			go.Scatter(x=x, y=trendline_y_above_30, mode='lines',
					   name=f'Linear fit (counts > 30): {p_above_30}',
					   line=dict(color='pink')))

		# Update layout
		fig.update_traces(marker_size=3)
		fig.update_layout(
			title=f"Simulation Protein Counts ({sim_name}) "
				  f"vs. Validation Protein Counts ({val_name} et al.),<br> $R^2$ > 30: {round(r_squared_30_above,3)}, n={len(above_30_idx[0])} (of {len(overlapIDs)} total)",
			xaxis_title="log10(Validation Protein Counts+1)",
			yaxis_title=f"log10(Simulation Protein Counts+1)",
			autosize=False, width=900, height=600)

		# add a y=x line
		fig.add_trace(
			go.Scatter(x=[0, 6], y=[0, 6], mode="lines",
					   line=go.scatter.Line(color="black", dash="dash"),
					   opacity=0.2, name="y=x"));

		# save the figure as an html:
		plot_name = f"proteinCountsValidation_cohortPlot_{sim_name}_vs_{val_name}_complex_fraction_highlighted.html"
		fig.write_html(os.path.join(plotOutDir, plot_name))

	hi = 5
	def plot_by_complex_fraction_with_proteins_highlighted_plotly(self, simDataFile, plotOutDir, simulationCounts, validationCounts,
								 overlapIDs, sim_name, val_name):
		"""
		Extract the complex count fractions for each protein in the simulation
		Args:
			simDataFile: simulation data file

		Returns: a dictionary of the complex count fractions for each protein
		"""

		monomer_counts_table, monomer_id_to_complex_fraction, monomer_id_to_complex_counts = self.determine_fraction_table(simDataFile)

		# Compute log10 values
		x = self.get_LogData(overlapIDs, validationCounts)
		y = self.get_LogData(overlapIDs, simulationCounts)

		monomer_to_half_life = self.get_half_lives(simDataFile)

		# retreive the common names and the descriptions:
		common_names = [self.get_common_name(protein_id) for protein_id in overlapIDs]
		descriptive_names = [self.get_descriptive_name(protein_id) for protein_id in overlapIDs]


		# create a dataframe of the protein ids and their half life sources:
		protein_df = pd.DataFrame({"protein_id": overlapIDs,
								   "common_name": common_names,
								   "descriptive_name": descriptive_names,
								   'simulation_protein_counts': y,
								   'validation_protein_counts': x})
		protein_df['fraction_in_complex'] = protein_df['protein_id'].map(monomer_id_to_complex_fraction)
		protein_df['complex_counts'] = protein_df['protein_id'].map(monomer_id_to_complex_counts)
		protein_df['half_life'] = protein_df['protein_id'].map(monomer_to_half_life)

		# make a yes or no column for whether the protein is in complex form:
		protein_df['in_complex'] = protein_df['fraction_in_complex'].apply(lambda x: 'Yes' if x > COMPLEX_FRACTION_THRESHOLD else 'No')

		# Edit the in complex fraction to say "HIGHLIGHTED_PROTEINS" if the protein is in the highlight list:
		protein_df.loc[protein_df['protein_id'].isin(
			HIGHLIGHT_PROTEINS), 'in_complex'] = 'HIGHLIGHTED_PROTEINS'

		# split up the data frames by complex fraction:
		complex_fraction_types = protein_df['in_complex'].unique()
		complex_fraction_dfs = {}
		for type in complex_fraction_types:
			complex_fraction_dfs[type] = protein_df[protein_df['in_complex'] == type]

		# create a color map for each complex fraction:
		color_map = {
			"Yes": "lightseagreen",
			"No": "yellowgreen"}
		color_map['HIGHLIGHTED_PROTEINS'] = 'purple'

		name_map = {
			"Yes": f"Complex fraction above {COMPLEX_FRACTION_THRESHOLD}",
			"No": f"Complex fraction below {COMPLEX_FRACTION_THRESHOLD}"}
		name_map['HIGHLIGHTED_PROTEINS'] = 'Highlighted Proteins'

		# create a scatter plot for each half life source:
		plt.figure(figsize=(9, 6), dpi=200)
		fig = go.Figure()
		for type, df in complex_fraction_dfs.items():
			print(type)
			c = color_map[type]; name = name_map[type]
			hovertext = self.hover_text_info(df)
			fig.add_trace(go.Scatter(x=df['validation_protein_counts'], y=df['simulation_protein_counts'], hovertext=hovertext, mode='markers',
									 name=f"{name} (n={len(df)})",
									 marker=dict(color=c, size=.7, opacity=.5)))

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

		# compute the rsquared value:
		r_squared_30_above = r2_score(x_above_30, y_above_30)

		# Add trendline trace
		fig.add_trace(
			go.Scatter(x=x, y=trendline_y, mode='lines',
					   name=f'Linear fit: {p}',
					   line=dict(color='green')))
		fig.add_trace(
			go.Scatter(x=x, y=trendline_y_above_30, mode='lines',
					   name=f'Linear fit (counts > 30): {p_above_30}',
					   line=dict(color='pink')))

		# Update layout
		fig.update_traces(marker_size=3)
		fig.update_layout(
			title=f"Simulation Protein Counts ({sim_name}) "
				  f"vs. Validation Protein Counts ({val_name} et al.),<br> R^2 > 30: {round(r_squared_30_above,3)}, n={len(above_30_idx[0])} (of {len(overlapIDs)} total)",
			xaxis_title="log10(Validation Protein Counts+1)",
			yaxis_title=f"log10(Simulation Protein Counts+1)",
			autosize=False, width=900, height=600)

		# add a y=x line
		fig.add_trace(
			go.Scatter(x=[0, 6], y=[0, 6], mode="lines",
					   line=go.scatter.Line(color="black", dash="dash"),
					   opacity=0.2, name="y=x"));

		# save the figure as an html:
		plot_name = f"proteinCountsValidation_cohortPlot_{sim_name}_vs_{val_name}_complex_fraction_above_{COMPLEX_FRACTION_THRESHOLD}_highlighted.html"
		fig.write_html(os.path.join(plotOutDir, plot_name))

	def validation_data_source_comparison_plot(self, sim_schmidt_counts, val_schmidt_counts,
		 schmidt_overlap_ids, schmidt_name, sim_wisniewski_counts,
		 val_wisniewski_counts, wisniewski_overlap_ids, wisniewski_name, plotOutDir):

		# # find the validation data values that overlap with simulation proteins from both validation sets:
		# do not use set, as it does not guarentee the same order every iteration:
		validation_overlap = [n for n in schmidt_overlap_ids if n in wisniewski_overlap_ids]

		# create a map for each dataset:
		schmidt_map = {n: i for i, n in enumerate(schmidt_overlap_ids)}
		wisniewski_map = {n: i for i, n in enumerate(wisniewski_overlap_ids)}

		# map to the indicies of the original overlapping proteins (with the sim data) datasets for each validation source:
		schmidt_idxs = [schmidt_map[n] for n in validation_overlap]
		wisniewski_idxs = [wisniewski_map[n] for n in validation_overlap]


		# extract the overlapping counts:
		x = self.get_LogData(list(range(len(validation_overlap))), val_schmidt_counts, schmidt_idxs)
		y = self.get_LogData(list(range(len(validation_overlap))), val_wisniewski_counts, wisniewski_idxs)

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

		# compute the rsquared value:
		r_squared_30_above = r2_score(x_above_30, y_above_30)

		hovertext = np.array(validation_overlap)

		# Add scatter trace
		fig = go.Figure()
		fig.add_trace(
			go.Scatter(x=x, y=y, hovertext=hovertext, mode='markers',
					   name=f"Counts (R^2 for counts > 30:{round(r_squared_30_above, 3)})"))

		# Add trendline trace
		fig.add_trace(
			go.Scatter(x=x, y=trendline_y, mode='lines',
					   name=f'Linear fit: {p}',
					   line=dict(color='green')))
		fig.add_trace(
			go.Scatter(x=x, y=trendline_y_above_30, mode='lines',
					   name=f'Linear fit (counts > 30): {p_above_30}',
					   line=dict(color='pink')))

		# Update layout
		fig.update_traces(marker_size=3)
		fig.update_layout(
			title=f"Validation Protein Counts Comparison: {wisniewski_name} et al. "
				  f"vs. {schmidt_name} et al.<br> (n={len(validation_overlap)} plotted)",
			xaxis_title=f"log10({schmidt_name} et al. Counts + 1))",
			yaxis_title=f"log10({wisniewski_name} et al. Counts + 1))",
			autosize=False, width=900, height=600)

		# add a y=x line
		fig.add_trace(
			go.Scatter(x=[0, 6], y=[0, 6], mode="lines",
					   line=go.scatter.Line(color="black", dash="dash"),
					   opacity=0.2, name="y=x"));

		# save the figure as an html:
		plot_name = f"proteinCountsValidation_cohortPlot_validation_source_comparison_{schmidt_name}_vs_{wisniewski_name}.html"
		fig.write_html(os.path.join(plotOutDir, plot_name))

	# compare the Schmidt validation sources to each other:
	def compare_Schmidt_BW_to_MG(self, simDataFile, validationDataFile, plotOutDir):
		# Obtain Schmit validation table:
		SDF = self.make_schmidt_validation_table(simDataFile, validationDataFile)

		# obtain the protein IDs:
		schmidt_IDs = SDF['Protein_ID'].values

		# obtain the BW and MG counts:
		val_schmidt_counts_bw = SDF['Val_BW_Counts'].values
		val_schmidt_counts_mg = SDF['Val_MG_Counts'].values
		sim_counts = SDF['Sim_Counts'].values


		# plot the comparison of the validation counts:
		self.validation_data_source_comparison_plot(sim_counts, val_schmidt_counts_bw,
													schmidt_IDs, "BW25113 Schmidt from validation file",
													sim_counts,
													val_schmidt_counts_mg, schmidt_IDs,
													"MG1655 Schmidt", plotOutDir)


	# Make a plotly of the ST9 counts vs the sim counts:
	def compare_ST9_MG_to_simulation(self, simDataFile, plotOutDir):
		# Obtain ST9 validation table:
		ST9_dict = self.get_schmidt_MG_validation_data_from_ST9()

		# Generate the simulation data:
		avg_total_counts, avg_free_counts, avg_counts_for_monomers_in_complexes = self.generate_data(
			simDataFile)

		# Make dictonaries for mapping:
		sim_genes_to_counts = {}
		sim_genes_to_protein_ids = {}
		sim_genes_to_descriptions = {}
		for i, monomer_id in enumerate(self.all_monomer_ids):
			gene_name = self.get_common_name(monomer_id)
			sim_genes_to_counts[gene_name] = avg_total_counts[i]
			sim_genes_to_protein_ids[gene_name] = monomer_id
			sim_genes_to_descriptions[gene_name] = self.get_descriptive_name(monomer_id)

		# Create a dataframe for ST9 data with simulation counts:
		SDF = []
		for gene in sim_genes_to_protein_ids.keys():
			if gene in ST9_dict:
				row = {
					'Gene_Symbol': gene,
					'Protein_ID': sim_genes_to_protein_ids[gene],
					'Descriptive_Name': sim_genes_to_descriptions[gene],
					'Sim_Counts': sim_genes_to_counts[gene],
					'Val_MG_Counts': ST9_dict[gene]
				}
				SDF.append(row)

		# Convert the combined data into a DataFrame
		SDF_df = pd.DataFrame(SDF)

		# create a plotly comparing the ST9 counts (MG1655 Schmidt et al.) to the simulation counts:
		simulationCounts = SDF_df['Sim_Counts'].values
		validationCounts = SDF_df['Val_MG_Counts'].values
		overlapIDs = SDF_df['Protein_ID'].values
		sim_name = "Simulation"
		val_name = "ST9 MG1655 Glucose Schmidt et al."

		# Generate the plotly:
		fig = go.Figure()

		# Compute log10 values
		x = self.get_LogData(overlapIDs, validationCounts)
		y = self.get_LogData(overlapIDs, simulationCounts)

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

		# compute the pearsonr and R2:
		r_value, p_val = pearsonr(x_above_30, y_above_30)
		pr2 = r_value ** 2
		# compute the coefficent of determination rsquared value:
		r_squared_30_above = r2_score(x_above_30, y_above_30)

		# have the monomer IDs be the overlap text
		monomer_counts_table, monomer_id_to_complex_fraction, monomer_id_to_complex_counts = self.determine_fraction_table(
			simDataFile)
		monomer_to_half_life = self.get_half_lives(simDataFile)

		# retreive the common names and the descriptions:
		common_names = [self.get_common_name(protein_id) for protein_id in overlapIDs]
		descriptive_names = [self.get_descriptive_name(protein_id) for protein_id in overlapIDs]

		# create a dataframe of the protein ids and their half life sources:
		protein_df = pd.DataFrame({"protein_id": overlapIDs,
								   "common_name": common_names,
								   "descriptive_name": descriptive_names,
								   'simulation_protein_counts': y,
								   'validation_protein_counts': x})
		protein_df['fraction_in_complex'] = protein_df['protein_id'].map(
			monomer_id_to_complex_fraction)
		protein_df['complex_counts'] = protein_df['protein_id'].map(monomer_id_to_complex_counts)
		protein_df['half_life'] = protein_df['protein_id'].map(monomer_to_half_life)
		hovertext = protein_df.apply(lambda
										 row: f"Monomer ID: {row['protein_id']}<br>Common Name: {row['common_name']}<br>Description: {row['descriptive_name']}<br>HL Value: {row['half_life']}<br>Validation AMC: {10 ** (row['validation_protein_counts'])}<br>Simulation AMC: {10 ** (row['simulation_protein_counts'])}<br>Avg. Complexed Monomer Counts: {row['complex_counts']} Complexed Fraction: {row['fraction_in_complex']}<br>",
									 axis=1)

		# Add scatter trace
		fig.add_trace(
			go.Scatter(x=x, y=y, hovertext=hovertext, mode='markers',
					   name=f"Average Monomer Counts"))

		# Add trendline trace
		fig.add_trace(
			go.Scatter(x=x, y=trendline_y, mode='lines',
					   name=f'Linear fit: {p}',
					   line=dict(color='green')))
		fig.add_trace(
			go.Scatter(x=x, y=trendline_y_above_30, mode='lines',
					   name=f'Linear fit (counts > 30): {p_above_30}',
					   line=dict(color='pink')))

		# Update layout
		fig.update_traces(marker_size=3)
		fig.update_traces(marker_size=3)
		fig.update_layout(
			title=f"Simulation Protein Counts ({sim_name}) "
				  f"vs. Validation Protein Counts ({val_name}) <br> Pearson R<sup>2</sup> counts > 30: {round(pr2, 3)}, n={len(above_30_idx[0])} (of {len(overlapIDs)} total)",
			title_font=dict(size=12),  # Make the title font smaller
			xaxis_title="log10(Validation Protein Counts)",
			yaxis_title="log10(Simulation Protein Counts)",
			autosize=False,
			width=900,  # Set equal width and height for a square graph
			height=600,
			plot_bgcolor='white',  # Set the plot area background color to white
			paper_bgcolor='white'  # Set the entire graph background to white
		)

		# add a y=x line
		fig.add_trace(
			go.Scatter(x=[0, 6], y=[0, 6], mode="lines",
					   line=go.scatter.Line(color="black", dash="dash"),
					   opacity=0.2, name="y=x"));

		# Define the text to display
		text = (
			f'Pearson R (counts > 30): {round(r_value, 3)}<br>'
			f'Pearson R<sup>2</sup> (counts > 30): {round(pr2, 3)}<br>'
			f'Coefficient of determination R<sup>2</sup> (counts > 30): {round(r_squared_30_above, 3)}'
		)

		# Get the maximum x and minimum y to position the text in the bottom-right
		x_max = x.max()
		y_min = y.min()

		# Adding text annotation just outside the graph
		# Adjust x_max and y_min slightly outside the actual graph boundaries
		text_offset_x = 0.05  # Offset to place the text box outside the graph
		text_offset_y = 0.05  # Adjust according to your needs

		# Adding text annotation to the bottom right
		fig.add_annotation(
			x=x_max + text_offset_x,  # Move the x position slightly to the right
			y=y_min - text_offset_y,  # Move the y position slightly below
			text=text,  # Your custom multi-line text
			showarrow=False,  # No arrow
			bgcolor='rgba(255, 255, 255, 0.8)',
			# Slightly less transparent white background for better visibility
			bordercolor='rgba(0, 0, 0, 0.5)',  # Optional border color
			borderwidth=1,  # Optional border width
			borderpad=4,  # Padding around the text
			align='right',  # Align the text to the right
			font=dict(size=10, color='gray'),  # Font properties
			xref='x',  # Reference for x-coordinate
			yref='y',  # Reference for y-coordinate
		)

		# save the figure as an html:
		plot_name = f"proteinCountsValidation_cohortPlot_{sim_name}_vs_{val_name}.html"
		fig.write_html(os.path.join(plotOutDir, plot_name))


	# compare ST9 BW25113 to sim data:
	def compare_ST9_BW_to_simulation(self, simDataFile, plotOutDir):
		# Obtain ST9 validation table:
		ST9_BW_dict = self.get_schmidt_BW_validation_data_from_ST9()

		# Generate the simulation data:
		avg_total_counts, avg_free_counts, avg_counts_for_monomers_in_complexes = self.generate_data(
			simDataFile)

		# Make dictonaries for mapping:
		sim_genes_to_counts = {}
		sim_genes_to_protein_ids = {}
		sim_genes_to_descriptions = {}
		for i, monomer_id in enumerate(self.all_monomer_ids):
			gene_name = self.get_common_name(monomer_id)
			sim_genes_to_counts[gene_name] = avg_total_counts[i]
			sim_genes_to_protein_ids[gene_name] = monomer_id
			sim_genes_to_descriptions[gene_name] = self.get_descriptive_name(monomer_id)

		# Create a dataframe for ST9 data with simulation counts:
		SDF = []
		for gene in sim_genes_to_protein_ids.keys():
			if gene in ST9_BW_dict:
				row = {
					'Gene_Symbol': gene,
					'Protein_ID': sim_genes_to_protein_ids[gene],
					'Descriptive_Name': sim_genes_to_descriptions[gene],
					'Sim_Counts': sim_genes_to_counts[gene],
					'Val_BW_Counts': ST9_BW_dict[gene]
				}
				SDF.append(row)

		# Convert the combined data into a DataFrame
		SDF_df = pd.DataFrame(SDF)

		# create a plotly comparing the ST9 counts (MG1655 Schmidt et al.) to the simulation counts:
		simulationCounts = SDF_df['Sim_Counts'].values
		validationCounts = SDF_df['Val_BW_Counts'].values
		overlapIDs = SDF_df['Protein_ID'].values
		sim_name = self.sim_name
		val_name = "ST9 BW25113 Glucose Schmidt et al."

		# Generate the plotly:
		fig = go.Figure()

		# Compute log10 values
		x = self.get_LogData(overlapIDs, validationCounts)
		y = self.get_LogData(overlapIDs, simulationCounts)

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

		# compute the pearsonr and R2:
		r_value, p_val = pearsonr(x_above_30, y_above_30)
		pr2 = r_value ** 2
		# compute the coefficent of determination rsquared value:
		r_squared_30_above = r2_score(x_above_30, y_above_30)

		# have the monomer IDs be the overlap text
		monomer_counts_table, monomer_id_to_complex_fraction, monomer_id_to_complex_counts = self.determine_fraction_table(
			simDataFile)
		monomer_to_half_life = self.get_half_lives(simDataFile)

		# retreive the common names and the descriptions:
		common_names = [self.get_common_name(protein_id) for protein_id in overlapIDs]
		descriptive_names = [self.get_descriptive_name(protein_id) for protein_id in overlapIDs]

		# create a dataframe of the protein ids and their half life sources:
		protein_df = pd.DataFrame({"protein_id": overlapIDs,
								   "common_name": common_names,
								   "descriptive_name": descriptive_names,
								   'simulation_protein_counts': y,
								   'validation_protein_counts': x})
		protein_df['fraction_in_complex'] = protein_df['protein_id'].map(
			monomer_id_to_complex_fraction)
		protein_df['complex_counts'] = protein_df['protein_id'].map(monomer_id_to_complex_counts)
		protein_df['half_life'] = protein_df['protein_id'].map(monomer_to_half_life)
		hovertext = protein_df.apply(lambda
					row: f"Monomer ID: {row['protein_id']}<br>Common Name: {row['common_name']}<br>Description: {row['descriptive_name']}<br>HL Value: {row['half_life']}<br>Validation AMC: {10 ** (row['validation_protein_counts'])}<br>Simulation AMC: {10 ** (row['simulation_protein_counts'])}<br>Avg. Complexed Monomer Counts: {row['complex_counts']} Complexed Fraction: {row['fraction_in_complex']}<br>",
									axis=1)

		# Add scatter trace
		fig.add_trace(
			go.Scatter(x=x, y=y, hovertext=hovertext, mode='markers',
					   name=f"Average Monomer Counts"))

		# Add trendline trace
		fig.add_trace(
			go.Scatter(x=x, y=trendline_y, mode='lines',
					   name=f'Linear fit: {p}',
					   line=dict(color='green')))
		fig.add_trace(
			go.Scatter(x=x, y=trendline_y_above_30, mode='lines',
					   name=f'Linear fit (counts > 30): {p_above_30}',
					   line=dict(color='pink')))

		# Update layout
		fig.update_traces(marker_size=3)
		fig.update_traces(marker_size=3)
		fig.update_layout(
			title=f"Simulation Protein Counts ({sim_name}) "
				  f"vs. Validation Protein Counts ({val_name}) <br> Pearson R<sup>2</sup> counts > 30: {round(pr2, 3)}, n={len(above_30_idx[0])} (of {len(overlapIDs)} total)",
			title_font=dict(size=12),  # Make the title font smaller
			xaxis_title="log10(Validation Protein Counts)",
			yaxis_title="log10(Simulation Protein Counts)",
			autosize=False,
			width=900,  # Set equal width and height for a square graph
			height=600,
			plot_bgcolor='white',  # Set the plot area background color to white
			paper_bgcolor='white'  # Set the entire graph background to white
		)

		# add a y=x line
		fig.add_trace(
			go.Scatter(x=[0, 6], y=[0, 6], mode="lines",
					   line=go.scatter.Line(color="black", dash="dash"),
					   opacity=0.2, name="y=x"));

		# Define the text to display
		text = (
			f'Pearson R (counts > 30): {round(r_value, 3)}<br>'
			f'Pearson R<sup>2</sup> (counts > 30): {round(pr2, 3)}<br>'
			f'Coefficient of determination R<sup>2</sup> (counts > 30): {round(r_squared_30_above, 3)}'
		)

		# Get the maximum x and minimum y to position the text in the bottom-right
		x_max = x.max()
		y_min = y.min()

		# Adding text annotation just outside the graph
		# Adjust x_max and y_min slightly outside the actual graph boundaries
		text_offset_x = 0.05  # Offset to place the text box outside the graph
		text_offset_y = 0.05  # Adjust according to your needs

		# Adding text annotation to the bottom right
		fig.add_annotation(
			x=x_max + text_offset_x,  # Move the x position slightly to the right
			y=y_min - text_offset_y,  # Move the y position slightly below
			text=text,  # Your custom multi-line text
			showarrow=False,  # No arrow
			bgcolor='rgba(255, 255, 255, 0.8)',
			# Slightly less transparent white background for better visibility
			bordercolor='rgba(0, 0, 0, 0.5)',  # Optional border color
			borderwidth=1,  # Optional border width
			borderpad=4,  # Padding around the text
			align='right',  # Align the text to the right
			font=dict(size=10, color='gray'),  # Font properties
			xref='x',  # Reference for x-coordinate
			yref='y',  # Reference for y-coordinate
		)

		# save the figure as an html:
		plot_name = f"proteinCountsValidation_cohortPlot_{sim_name}_vs_{val_name}.html"
		fig.write_html(os.path.join(plotOutDir, plot_name))

	# compare ST6 BW Schmidt data to sim data:
	def compare_ST6_BW_to_simulation(self, simDataFile, plotOutDir):
		# Obtain ST9 validation table:
		ST6_BW_dict = self.get_schmidt_BW_validation_data_from_ST6()

		# Generate the simulation data:
		avg_total_counts, avg_free_counts, avg_counts_for_monomers_in_complexes = self.generate_data(
			simDataFile)

		# Make dictonaries for mapping:
		sim_genes_to_counts = {}
		sim_genes_to_protein_ids = {}
		sim_genes_to_descriptions = {}
		for i, monomer_id in enumerate(self.all_monomer_ids):
			gene_name = self.get_common_name(monomer_id)
			sim_genes_to_counts[gene_name] = avg_total_counts[i]
			sim_genes_to_protein_ids[gene_name] = monomer_id
			sim_genes_to_descriptions[gene_name] = self.get_descriptive_name(monomer_id)

		# Create a dataframe for ST9 data with simulation counts:
		SDF = []
		for gene in sim_genes_to_protein_ids.keys():
			if gene in ST6_BW_dict:
				row = {
					'Gene_Symbol': gene,
					'Protein_ID': sim_genes_to_protein_ids[gene],
					'Descriptive_Name': sim_genes_to_descriptions[gene],
					'Sim_Counts': sim_genes_to_counts[gene],
					'Val_BW_Counts': ST6_BW_dict[gene]
				}
				SDF.append(row)

		# Convert the combined data into a DataFrame
		SDF_df = pd.DataFrame(SDF)

		# create a plotly comparing the ST9 counts (MG1655 Schmidt et al.) to the simulation counts:
		simulationCounts = SDF_df['Sim_Counts'].values
		validationCounts = SDF_df['Val_BW_Counts'].values
		overlapIDs = SDF_df['Protein_ID'].values
		sim_name = self.sim_name
		val_name = "ST6 BW25113 Glucose Schmidt et al."

		# Generate the plotly:
		fig = go.Figure()

		# Compute log10 values
		x = self.get_LogData(overlapIDs, validationCounts)
		y = self.get_LogData(overlapIDs, simulationCounts)

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

		# compute the pearsonr and R2:
		r_value, p_val = pearsonr(x_above_30, y_above_30)
		pr2 = r_value ** 2
		# compute the coefficent of determination rsquared value:
		r_squared_30_above = r2_score(x_above_30, y_above_30)

		# have the monomer IDs be the overlap text
		monomer_counts_table, monomer_id_to_complex_fraction, monomer_id_to_complex_counts = self.determine_fraction_table(
			simDataFile)
		monomer_to_half_life = self.get_half_lives(simDataFile)

		# retreive the common names and the descriptions:
		common_names = [self.get_common_name(protein_id) for protein_id in overlapIDs]
		descriptive_names = [self.get_descriptive_name(protein_id) for protein_id in overlapIDs]

		# create a dataframe of the protein ids and their half life sources:
		protein_df = pd.DataFrame({"protein_id": overlapIDs,
								   "common_name": common_names,
								   "descriptive_name": descriptive_names,
								   'simulation_protein_counts': y,
								   'validation_protein_counts': x})
		protein_df['fraction_in_complex'] = protein_df['protein_id'].map(
			monomer_id_to_complex_fraction)
		protein_df['complex_counts'] = protein_df['protein_id'].map(monomer_id_to_complex_counts)
		protein_df['half_life'] = protein_df['protein_id'].map(monomer_to_half_life)
		hovertext = protein_df.apply(lambda
					row: f"Monomer ID: {row['protein_id']}<br>Common Name: {row['common_name']}<br>Description: {row['descriptive_name']}<br>HL Value: {row['half_life']}<br>Validation AMC: {10 ** (row['validation_protein_counts'])}<br>Simulation AMC: {10 ** (row['simulation_protein_counts'])}<br>Avg. Complexed Monomer Counts: {row['complex_counts']} Complexed Fraction: {row['fraction_in_complex']}<br>",
									axis=1)

		# Add scatter trace
		fig.add_trace(
			go.Scatter(x=x, y=y, hovertext=hovertext, mode='markers',
					   name=f"Average Monomer Counts"))

		# Add trendline trace
		fig.add_trace(
			go.Scatter(x=x, y=trendline_y, mode='lines',
					   name=f'Linear fit: {p}',
					   line=dict(color='green')))
		fig.add_trace(
			go.Scatter(x=x, y=trendline_y_above_30, mode='lines',
					   name=f'Linear fit (counts > 30): {p_above_30}',
					   line=dict(color='pink')))

		# Update layout
		fig.update_traces(marker_size=3)
		fig.update_layout(
			title=f"Simulation Protein Counts ({sim_name}) "
				  f"vs. Validation Protein Counts ({val_name}) <br> Pearson R<sup>2</sup> counts > 30: {round(pr2, 3)}, n={len(above_30_idx[0])} (of {len(overlapIDs)} total)",
			title_font=dict(size=12),  # Make the title font smaller
			xaxis_title="log10(Validation Protein Counts)",
			yaxis_title="log10(Simulation Protein Counts)",
			autosize=False,
			width=900,  # Set equal width and height for a square graph
			height=600,
			plot_bgcolor='white',  # Set the plot area background color to white
			paper_bgcolor='white'  # Set the entire graph background to white
		)

		# add a y=x line
		fig.add_trace(
			go.Scatter(x=[0, 6], y=[0, 6], mode="lines",
					   line=go.scatter.Line(color="black", dash="dash"),
					   opacity=0.2, name="y=x"));

		# Define the text to display
		text = (
			f'Pearson R (counts > 30): {round(r_value, 3)}<br>'
			f'Pearson R<sup>2</sup> (counts > 30): {round(pr2, 3)}<br>'
			f'Coefficient of determination R<sup>2</sup> (counts > 30): {round(r_squared_30_above, 3)}'
		)

		# Get the maximum x and minimum y to position the text in the bottom-right
		x_max = x.max()
		y_min = y.min()

		# Adding text annotation just outside the graph
		# Adjust x_max and y_min slightly outside the actual graph boundaries
		text_offset_x = 0.05  # Offset to place the text box outside the graph
		text_offset_y = 0.05  # Adjust according to your needs

		# Adding text annotation to the bottom right
		fig.add_annotation(
			x=x_max + text_offset_x,  # Move the x position slightly to the right
			y=y_min - text_offset_y,  # Move the y position slightly below
			text=text,  # Your custom multi-line text
			showarrow=False,  # No arrow
			bgcolor='rgba(255, 255, 255, 0.8)',
			# Slightly less transparent white background for better visibility
			bordercolor='rgba(0, 0, 0, 0.5)',  # Optional border color
			borderwidth=1,  # Optional border width
			borderpad=4,  # Padding around the text
			align='right',  # Align the text to the right
			font=dict(size=10, color='gray'),  # Font properties
			xref='x',  # Reference for x-coordinate
			yref='y',  # Reference for y-coordinate
		)

		# save the figure as an html:
		plot_name = f"proteinCountsValidation_cohortPlot_{sim_name}_vs_{val_name}.html"
		fig.write_html(os.path.join(plotOutDir, plot_name))


	# finally, compare ST6 BW to ST9 BW:
	def compare_ST6_BW_to_ST9_BW(self, validationDataFile, plotOutDir):
		# Obtain ST9 and ST6 BW validation tables:
		ST9_BW_dict = self.get_schmidt_BW_validation_data_from_ST9()
		ST6_BW_dict = self.get_schmidt_BW_validation_data_from_ST6()

		# find overlapping genes:
		overlap_genes = set(ST9_BW_dict.keys()).intersection(set(ST6_BW_dict.keys()))

		# obtain the counts for the overlapping genes:
		ST9_counts = [ST9_BW_dict[gene] for gene in overlap_genes]
		ST6_counts = [ST6_BW_dict[gene] for gene in overlap_genes]

		# plot the comparison of the validation counts:
		self.validation_data_source_comparison_plot(ST9_counts, ST6_counts,
													list(overlap_genes), "ST9 BW25113 Schmidt",
													ST9_counts,
													ST6_counts, list(overlap_genes),
													"ST6 BW25113 Schmidt", plotOutDir)


	# compare the validation BW output to the direct read in of ST6 BW data:
	def compare_validation_file_data_to_direct_ST6_BW(self, simDataFile, validationDataFile, plotOutDir):
		# Obtain Schmit validation table:
		SDF = self.make_schmidt_validation_table_comparing_val_BW_to_direct_ST6(simDataFile, validationDataFile)

		# obtain the protein IDs:
		schmidt_IDs = SDF['Protein_ID'].values

		# obtain the BW and MG counts:
		val_schmidt_counts_bw = SDF['Val_BW_Counts'].values
		val_schmidt_counts_bw_from_ST6 = SDF['ST6_BW_Counts'].values
		sim_counts = SDF['Sim_Counts'].values


		# plot the comparison of the validation counts:
		self.validation_data_source_comparison_plot(sim_counts, val_schmidt_counts_bw,
													schmidt_IDs, "BW25113 Schmidt from validation file",
													sim_counts,
													val_schmidt_counts_bw_from_ST6, schmidt_IDs,
													"BW25113 Schmidt from ST6", plotOutDir)

	# plot the ST9 BW data against the validation file BW data:
	def compare_validation_file_data_to_direct_ST9_BW(self, simDataFile, validationDataFile, plotOutDir):
		# Obtain Schmit validation table:
		SDF = self.make_schmidt_validation_table_comparing_val_BW_to_direct_ST9(simDataFile, validationDataFile)

		# obtain the protein IDs:
		schmidt_IDs = SDF['Protein_ID'].values

		# obtain the BW and MG counts:
		val_schmidt_counts_bw = SDF['Val_BW_Counts'].values
		val_schmidt_counts_bw_from_ST9 = SDF['ST9_BW_Counts'].values
		sim_counts = SDF['Sim_Counts'].values


		# plot the comparison of the validation counts:
		self.validation_data_source_comparison_plot(sim_counts, val_schmidt_counts_bw,
													schmidt_IDs, "BW25113 Schmidt from validation file",
													sim_counts,
													val_schmidt_counts_bw_from_ST9, schmidt_IDs,
													"BW25113 Schmidt from ST9", plotOutDir)




	def plot_validation_comparison(self, simDataFile, validationDataFile,plotOutDir, sim_name):

		# obtain overlapping protein counts between the simulation and validation data
		(sim_schmidt_counts, val_schmidt_counts,
		 schmidt_overlap_ids, sim_wisniewski_counts,
		 val_wisniewski_counts, wisniewski_overlap_ids) = (
			self.get_validation_data(simDataFile, validationDataFile))

		# plot the comparison of the validation counts:
		self.validation_data_source_comparison_plot(sim_schmidt_counts, val_schmidt_counts,
													schmidt_overlap_ids, "Schmidt",
													sim_wisniewski_counts,
													val_wisniewski_counts, wisniewski_overlap_ids,
													"Wisniewski", plotOutDir)

		#self.compare_Schmidt_BW_to_MG(simDataFile, validationDataFile, plotOutDir)
		#self.compare_ST9_BW_to_simulation(simDataFile, plotOutDir)
		#self.compare_ST9_MG_to_simulation(simDataFile, plotOutDir)
		#self.compare_ST6_BW_to_simulation(simDataFile, plotOutDir)
		#self.compare_validation_file_data_to_direct_ST6_BW(simDataFile, validationDataFile, plotOutDir)
		#self.compare_validation_file_data_to_direct_ST9_BW(simDataFile, validationDataFile, plotOutDir)
		#self.compare_ST6_BW_to_ST9_BW(validationDataFile, plotOutDir)


		# generate interactive validation plotlys:
		self.generate_validation_plotly(simDataFile, plotOutDir, sim_schmidt_counts,
			val_schmidt_counts, schmidt_overlap_ids, sim_name, "Schmidt")
		self.generate_validation_plotly(simDataFile, plotOutDir, sim_wisniewski_counts,
			val_wisniewski_counts, wisniewski_overlap_ids, sim_name, "Wisniewski")


		# generate matplotlib validation plots:
		self.generate_validation_plot(plotOutDir, sim_schmidt_counts,
			val_schmidt_counts, schmidt_overlap_ids, sim_name, "Schmidt")
		self.generate_validation_plot(plotOutDir, sim_wisniewski_counts,
			val_wisniewski_counts, wisniewski_overlap_ids, sim_name, "Wisniewski")

		# generate plotly validation plots with complex fraction highlighted:
		self.plot_by_complex_fraction_plotly(simDataFile, plotOutDir, sim_schmidt_counts,
			val_schmidt_counts, schmidt_overlap_ids, sim_name, "Schmidt")
		self.plot_by_complex_fraction_plotly(simDataFile, plotOutDir, sim_wisniewski_counts,
			val_wisniewski_counts, wisniewski_overlap_ids, sim_name, "Wisniewski")

		# generate plotly validation plots with complex fraction highlighted:
		self.plot_by_complex_fraction_with_proteins_highlighted_plotly(simDataFile, plotOutDir, sim_schmidt_counts,
			val_schmidt_counts, schmidt_overlap_ids, sim_name, "Schmidt")
		self.plot_by_complex_fraction_with_proteins_highlighted_plotly(simDataFile, plotOutDir, sim_wisniewski_counts,
			val_wisniewski_counts, wisniewski_overlap_ids, sim_name, "Wisniewski")
		hi = 6







	def do_plot(self, variantDir, plotOutDir, plotOutFileName, simDataFile,
				validationDataFile, metadata):

		# Generate the data for the simulation:
		self.sim_name = metadata["description"]
		self.n_total_gens = self.ap.n_generation
		self.plot_validation_comparison(simDataFile, validationDataFile, plotOutDir, self.sim_name)



if __name__ == '__main__':
	Plot().cli()
