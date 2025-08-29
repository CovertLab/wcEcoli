
import pickle
import os
import matplotlib.pyplot as plt
import csv
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

		# Calculate the average total monomer counts across all cells:
		total_protein_counts = (read_stacked_columns(all_cells,
			'MonomerCounts', 'monomerCounts')).mean(axis=0)

		# Make into a np.array:
		total_protein_counts = np.array(total_protein_counts)

		return total_protein_counts

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

		# Add scatter trace
		fig.add_trace(
			go.Scatter(x=x, y=y, hovertext=hovertext, mode='markers', name="Counts"))

		# Compute linear trendline
		z = np.polyfit(x, y, 1)
		p = np.poly1d(z)
		trendline_y = p(x)

		# Compute linear trendline for counts above log10(30+1): (+1 bc log(0) is undefined)
		above_30_idx = np.where((x > np.log10(30+1)) & (y > np.log10(30+1)))
		x_above_30 = x[above_30_idx]
		y_above_30 = y[above_30_idx]
		z_above_30 = np.polyfit(x_above_30, y_above_30, 1)
		p_above_30 = np.poly1d(z_above_30)
		trendline_y_above_30 = p_above_30(x)


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
		plt.figure(figsize=(8, 6), dpi=200)

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
		plt.scatter(x=x, y=y, s=5, label=f"Counts",
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
				  f"vs.\n Validation Protein Counts ({val_name} et al., 2016)")
		plt.xlabel("log10(Validation Counts)")
		plt.ylabel(f"log10(Simulation Counts)")
		plt.xlim(0, 6)
		plt.ylim(0, 6)
		plt.gca().set_aspect('equal', adjustable='box')

		plt.text(0.95, 0.02, f'$R^2$ for counts > 30: {round(r_squared_30_above, 2)}',
				 ha='right', va='bottom',
				 transform=plt.gca().transAxes,
				 fontsize=8, color='gray')
		plt.legend(bbox_to_anchor=(1.3,1.), loc='upper right', fontsize=6, )





		# save the figure as an html:
		plot_name = f"proteinCountsValidation_cohortPlot_{sim_name}_vs_{val_name}_matplotlib.pdf"
		plt.savefig(os.path.join(plotOutDir, plot_name))



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




	def do_plot(self, variantDir, plotOutDir, plotOutFileName, simDataFile,
				validationDataFile, metadata):

		# Generate the data for the simulation:
		sim_name = metadata["description"]
		self.n_total_gens = self.ap.n_generation
		self.plot_validation_comparison(simDataFile, validationDataFile, plotOutDir, sim_name)


if __name__ == '__main__':
	Plot().cli()
