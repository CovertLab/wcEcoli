
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

""" USER INPUTS """

# Indicate the number of generations to be ignored at the start of each seed:
IGNORE_FIRST_N_GENS = 2 # 2 for local, 14 for Sherlock (w/ 24 total gens)

# input proteins to highlight:
HIGHLIGHT_PROTEINS = ['G6890-MONOMER[c]',
 					   'PD03938[c]',
 					   'G6737-MONOMER[c]',
 					   'RPOD-MONOMER[c]',
 					   'PD02936[c]',
 					   'RED-THIOREDOXIN2-MONOMER[c]',
  						"EG10542-MONOMER[c]"]

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
		avg_counts_for_monomers_in_complexs = avg_total_counts - avg_free_counts

		return avg_total_counts, avg_free_counts, avg_counts_for_monomers_in_complexs

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

		# Initialize lists to store data for overlapping protein counts:
		sim_schmidt_counts_multigen = []
		sim_wisniewski_counts_multigen = []

		# Extract data for proteins that are present in the validation data:
		for simDir in all_cells:
			simOutDir = os.path.join(simDir, "simOut")
			monomer_counts_reader = (
				TableReader(os.path.join(simOutDir, "MonomerCounts")))
			monomer_counts = monomer_counts_reader.readColumn("monomerCounts")

			# Obtain the protein counts for protiens that are present in both
			# the simluation and validation datasets:
			sim_schmidt_counts, val_schmidt_counts, schmidt_overlap_ids = (
				get_simulated_validation_counts(
				schmidt_counts, monomer_counts, schmidt_ids, sim_monomer_ids))
			sim_wisniewski_counts, val_wisniewski_counts, wisniewski_overlap_ids = (
				get_simulated_validation_counts(wisniewski_counts,
					monomer_counts, wisniewski_ids, sim_monomer_ids))

			# Append the protein counts for the current cell:
			sim_schmidt_counts_multigen.append(sim_schmidt_counts)
			sim_wisniewski_counts_multigen.append(sim_wisniewski_counts)

		# Average over all the cells:
		sim_schmidt_counts = (
			(np.array(sim_schmidt_counts_multigen)).mean(axis=0))
		sim_wisniewski_counts = (
			(np.array(sim_wisniewski_counts_multigen)).mean(axis=0))

		return (sim_schmidt_counts, val_schmidt_counts,
				schmidt_overlap_ids, sim_wisniewski_counts,
				val_wisniewski_counts, wisniewski_overlap_ids)

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
		half_life_source = sim_data.process.translation.monomer_data["deg_rate_source"]

		# Create a dictionary of the half lives for each protein:
		total_protein_counts, _, _ = self.generate_data(simDataFile)
		monomer_to_half_life = dict(zip(self.all_monomer_ids, half_lives))
		monomer_to_half_life_source = dict(zip(self.all_monomer_ids, half_life_source))

		return monomer_to_half_life, monomer_to_half_life_source



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


	def generate_validation_plotly(self, plotOutDir, simulationCounts, validationCounts,
								 overlapIDs, sim_name, val_name):
		fig = go.Figure()
		# have the monomer IDs be the overlap text
		hovertext = overlapIDs

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

		# compute the rsquared value:
		r_squared_30_above = r2_score(x_above_30, y_above_30)

		# Add scatter trace
		plt.scatter(x=x, y=y, s=5, label=f"Counts ($R^2$ counts > 30: {round(r_squared_30_above,3)})",
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



		# save the figure as an html:
		plot_name = f"proteinCountsValidation_cohortPlot_{sim_name}_vs_{val_name}_matplotlib.pdf"
		plt.savefig(os.path.join(plotOutDir, plot_name))


	# function that colors the counts by their half live value
	def plot_with_color_by_half_life(self, simDataFile, plotOutDir, simulationCounts, validationCounts,
								 overlapIDs, sim_name, val_name):
		"""
		Color the points in the plot by their half life value
		Args:
			protein_ids: the ids of the proteins
			protein_counts: the counts of the proteins
			half_lives: the half lives of the proteins

		Returns:
		"""
		# Compute log10 values
		x = self.get_LogData(overlapIDs, validationCounts)
		y = self.get_LogData(overlapIDs, simulationCounts)

		monomer_to_half_life, monomer_to_half_life_source = self.get_half_lives(simDataFile)

		# create a dataframe of the protein ids and their half life sources:
		protein_df = pd.DataFrame({"protein_id": overlapIDs,
								   'simulation_protein_counts': y,
								   'validation_protein_counts': x})
		protein_df['half_life_source'] = protein_df['protein_id'].map(monomer_to_half_life_source)

		# split up the data frames by half life source:
		half_life_sources = protein_df['half_life_source'].unique()
		half_life_dfs = {}
		for source in half_life_sources:
			half_life_dfs[source] = protein_df[protein_df['half_life_source'] == source]

		# create a color map for each half life source:
		color_map = {
			"N_end_rule": "lightseagreen",
			"Gupta_et_al_MS_2024": "yellowgreen",
			"CL_measured_deg_rates_2020": "orange",
		}

		size_map = {
			"N_end_rule": 5,
			"Gupta_et_al_MS_2024": 2,
			"CL_measured_deg_rates_2020": 5,
		}

		alpha_map = {
			"N_end_rule": 0.5,
			"Gupta_et_al_MS_2024": 0.3,
			"CL_measured_deg_rates_2020": 0.9,
		}

		# create a scatter plot for each half life source:
		plt.figure(figsize=(9, 6), dpi=200)
		for source, df in half_life_dfs.items():
			c = color_map[source]; si = size_map[source]; a = alpha_map[source]
			print(source)
			plt.scatter(df['validation_protein_counts'], df['simulation_protein_counts'],
						label=source, color=c, alpha=a, s=si)

		# Compute linear trendline for counts above log10(30):
		above_30_idx = np.where((x > np.log10(30 + 1)) & (y > np.log10(30 + 1)))
		x_above_30 = x[above_30_idx]
		y_above_30 = y[above_30_idx]
		z_above_30 = np.polyfit(x_above_30, y_above_30, 1)
		p_above_30 = np.poly1d(z_above_30)
		trendline_y_above_30 = p_above_30(x)

		# compute the rsquared value:
		r_squared_30_above = r2_score(x_above_30, y_above_30)

		# Compute linear trendline
		z = np.polyfit(x, y, 1)
		p = np.poly1d(z)
		trendline_y = p(x)

		# Add trendline trace
		plt.plot(x, trendline_y,
				 label=f'Linear fit: {p}',
				 color='plum')
		plt.plot(x, trendline_y_above_30,
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
		plt.text(0.95, 0.02, f'$R^2$ for counts > 30: {round(r_squared_30_above,2)}',
				 ha='right', va='bottom',
				 transform=plt.gca().transAxes,
				 fontsize=8, color='gray')

		plt.legend(bbox_to_anchor=(1.45, 1.), loc='upper right', fontsize=8, )

		# save the figure as an html:
		plot_name = f"proteinCountsValidation_cohortPlot_{sim_name}_vs_{val_name}_half_life_source_highlighted_matplotlib.pdf"
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
										row: f"Monomer ID: {row['protein_id']}<br>HL Value: {row['half_life']}<br>HL Source: {row['half_life_source']}<br>validation AMC: {10 ** (row['validation_protein_counts'])}<br>Simulation AMC: {10 ** (row['simulation_protein_counts'])}<br>Avg. Complexed Monomer Counts: {row['complex_counts']} Complexed Fraction: {row['fraction_in_complex']}<br>",
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

		monomer_to_half_life, monomer_to_half_life_source = self.get_half_lives(simDataFile)

		# create a dataframe of the protein ids and their half life sources:
		protein_df = pd.DataFrame({"protein_id": overlapIDs,
								   'simulation_protein_counts': y,
								   'validation_protein_counts': x})
		protein_df['fraction_in_complex'] = protein_df['protein_id'].map(monomer_id_to_complex_fraction)
		protein_df['complex_counts'] = protein_df['protein_id'].map(monomer_id_to_complex_counts)
		protein_df['half_life_source'] = protein_df['protein_id'].map(monomer_to_half_life_source)
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
				  f"vs. Validation Protein Counts ({val_name} et al.),<br> R^2 > 30: {round(r_squared_30_above,3)}",
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

		# get the

		monomer_to_half_life, monomer_to_half_life_source = self.get_half_lives(simDataFile)

		# create a dataframe of the protein ids and their half life sources:
		protein_df = pd.DataFrame({"protein_id": overlapIDs,
								   'simulation_protein_counts': y,
								   'validation_protein_counts': x})
		protein_df['fraction_in_complex'] = protein_df['protein_id'].map(monomer_id_to_complex_fraction)
		protein_df['complex_counts'] = protein_df['protein_id'].map(monomer_id_to_complex_counts)
		protein_df['half_life_source'] = protein_df['protein_id'].map(monomer_to_half_life_source)
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
				  f"vs. Validation Protein Counts ({val_name} et al.),<br> R^2 > 30: {round(r_squared_30_above,3)}",
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




	def plot_validation_comparison(self, simDataFile, validationDataFile,plotOutDir, sim_name):

		# obtain overlapping protein counts between the simulation and validation data
		(sim_schmidt_counts, val_schmidt_counts,
		 schmidt_overlap_ids, sim_wisniewski_counts,
		 val_wisniewski_counts, wisniewski_overlap_ids) = (
			self.get_validation_data(simDataFile, validationDataFile))

		# generate interactive validation plotlys:
		self.generate_validation_plotly(plotOutDir, sim_schmidt_counts,
			val_schmidt_counts, schmidt_overlap_ids, sim_name, "Schmidt")
		self.generate_validation_plotly(plotOutDir, sim_wisniewski_counts,
			val_wisniewski_counts, wisniewski_overlap_ids, sim_name, "Wisniewski")


		# generate matplotlib validation plots:
		self.generate_validation_plot(plotOutDir, sim_schmidt_counts,
			val_schmidt_counts, schmidt_overlap_ids, sim_name, "Schmidt")
		self.generate_validation_plot(plotOutDir, sim_wisniewski_counts,
			val_wisniewski_counts, wisniewski_overlap_ids, sim_name, "Wisniewski")

		# generate matplotlib validation plots with half life source highlighted:
		self.plot_with_color_by_half_life(simDataFile, plotOutDir, sim_schmidt_counts,
			val_schmidt_counts, schmidt_overlap_ids, sim_name, "Schmidt")
		self.plot_with_color_by_half_life(simDataFile, plotOutDir, sim_wisniewski_counts,
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




	def do_plot(self, variantDir, plotOutDir, plotOutFileName, simDataFile,
				validationDataFile, metadata):

		# Generate the data for the simulation:
		sim_name = metadata["description"]
		self.n_total_gens = self.ap.n_generation
		self.plot_validation_comparison(simDataFile, validationDataFile, plotOutDir, sim_name)



if __name__ == '__main__':
	Plot().cli()
