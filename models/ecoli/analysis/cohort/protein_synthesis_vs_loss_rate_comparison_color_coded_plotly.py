"""
Template for cohort analysis plots
"""

import pickle
import os

from matplotlib import pyplot as plt
# noinspection PyUnresolvedReferences
import numpy as np
import plotly.graph_objects as go
from models.ecoli.analysis import cohortAnalysisPlot
from wholecell.analysis.analysis_tools import (exportFigure,
	read_bulk_molecule_counts, read_stacked_bulk_molecules, read_stacked_columns)
from wholecell.io.tablereader import TableReader

# TODO: make a separate plot with just the highlighted proteins and their standard deviations
# TODO: add more information to the hoverdata
# TODO: change plot symbols based on complex fraction


HIGHLIGHT_IN_RED = []#['EG10863-MONOMER[c]','DETHIOBIOTIN-SYN-MONOMER[c]','DCUR-MONOMER[c]']
HIGHLIGHT_IN_BLUE = []#['CARBPSYN-SMALL[c]', 'CDPDIGLYSYN-MONOMER[i]','EG10743-MONOMER[c]','GLUTCYSLIG-MONOMER[c]']
HIGHLIGHT_IN_PURPLE = ['G6890-MONOMER[c]','PD03938[c]','G6737-MONOMER[c]','RPOD-MONOMER[c]','PD02936[c]','RED-THIOREDOXIN2-MONOMER[c]']
class Plot(cohortAnalysisPlot.CohortAnalysisPlot):
	def do_plot(self, seedOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		with open(simDataFile, 'rb') as f:
			sim_data = pickle.load(f)
		with open(validationDataFile, 'rb') as f:
			validation_data = pickle.load(f)
		sim_id = metadata['description']
		cell_paths = self.ap.get_cells()
		sim_dir = cell_paths[0]
		simOutDir = os.path.join(sim_dir, 'simOut')

		# Extract protein indexes for each new gene
		monomer_counts_reader = TableReader(
			os.path.join(simOutDir, "MonomerCounts"))
		monomer_idx_dict = {monomer: i for i, monomer in
							enumerate(monomer_counts_reader.readAttribute(
								'monomerIds'))}

		# extract the numbers of interest
		monomerIds = monomer_counts_reader.readAttribute(
			"monomerIds")  # this is the one that matches the indexing  I used earlier to construct the listeners!


		# extract the data over all generations:
		# Load data
		time = read_stacked_columns(cell_paths, 'Main', 'time')
		(free_monomer_counts,) = read_stacked_bulk_molecules(
			cell_paths, monomerIds)

		# doubling time function from nora:
		def extract_doubling_times(cell_paths):
			# Load data
			time = read_stacked_columns(cell_paths, 'Main', 'time').squeeze()
			# Determine doubling time
			doubling_times = read_stacked_columns(cell_paths, 'Main', 'time',
												  fun=lambda x: (x[-1] - x[0])).squeeze().astype(
				int)
			end_generation_times = np.cumsum(doubling_times) + time[0]  #
			start_generation_indices = np.searchsorted(time, end_generation_times[:-1],
													   side='left').astype(int)
			start_generation_indices = np.insert(start_generation_indices, 0, 0) + np.arange(
				len(doubling_times))
			end_generation_indices = start_generation_indices + doubling_times
			return time, doubling_times, end_generation_times, start_generation_indices, end_generation_indices

		# extract the doubling times:
		time, doubling_times, end_generation_times, start_generation_indices, end_generation_indices = extract_doubling_times(
			cell_paths)

		# find how many proteins were removed via dilution for each doubling time:
		diluted_counts = np.zeros(((len(doubling_times) - 1), len(monomerIds)))
		diluted_counts_over_time = np.zeros(((len(time)), len(monomerIds)))
		for i in range(
				len(doubling_times) - 1):  # -1 is to account for not doing the last generation
			end_gen = end_generation_indices[i]
			start_gen = start_generation_indices[i + 1]  # the first start is zero, so skip that
			print(end_gen, start_gen)

			# find the protein counts at the end of the generation:
			monomer_counts_at_gen_end = free_monomer_counts[end_gen,:]  # get this for each protein

			# find the protein counts at the start of the next generation:
			monomer_counts_at_gen_start = free_monomer_counts[start_gen, :]

			# find the difference between the two:
			protein_counts_removed = monomer_counts_at_gen_end - monomer_counts_at_gen_start
			diluted_counts[i, :] = protein_counts_removed
			diluted_counts_over_time[start_gen,:] = protein_counts_removed  # put it at the start of the next gen in terms of time

		print(diluted_counts)

		# Get the free monomer counts:
		fmcs = read_stacked_columns(cell_paths, 'MonomerCounts', "freeMonomerCounts")

		# Compute how many proteins were removed via degradation over the entire sim length:
		degraded_counts = read_stacked_columns(cell_paths, 'MonomerCounts', "monomersDegraded")
		# take the average:
		avg_degraded_counts = np.mean(degraded_counts, axis = 0)

		# compute how many proteins were removed via dilution over the entire sim length:
		# todo: make sure this is the right way to compute the average (when you only have one dilution timepoint (becuase there was that one graph I did where it was different depending on the size of the array)
		total_diluted_counts = np.sum(diluted_counts, axis = 0)
		avg_diluted_counts = total_diluted_counts / len(time) # divide by the number of timesteps to get the average per timestep

		# compute the average loss rate for each protein:
		avg_loss_rate = avg_degraded_counts + avg_diluted_counts
		log_avg_loss_rate = np.log10(avg_loss_rate) # todo: consider adding 1 to this to avoid log(0) and all


		# compute how many counts were added via elongation over the entire sim length:
		elongated_counts = read_stacked_columns(cell_paths, 'MonomerCounts', "monomersElongated")
		avg_elongated_counts = np.mean(elongated_counts, axis = 0)
		log_avg_production_rate = np.log10(avg_elongated_counts) # todo: consider adding 1 to this to avoid log(0) and all and being able to see all proteins


		# plot the loss rate vs the production rate:
		# Create figure
		fig = go.Figure()

		# Scatter plot for all proteins (grey)
		fig.add_trace(go.Scatter(
			x=log_avg_production_rate,
			y=log_avg_loss_rate,
			mode='markers',
			hovertext=monomerIds,
			marker=dict(size=5, color='lightseagreen', opacity=0.3),
			name='All Proteins'
		))

		if len(HIGHLIGHT_IN_RED) > 0:
			# Scatter plot for red proteins
			red_protein_indices = [monomer_idx_dict[protein] for protein in HIGHLIGHT_IN_RED]
			fig.add_trace(go.Scatter(
				x=log_avg_production_rate[red_protein_indices],
				y=log_avg_loss_rate[red_protein_indices],
				mode='markers',
				hovertext=HIGHLIGHT_IN_RED,
				marker=dict(size=5, color='red', opacity=1),
				name='Red Proteins'
		))
		if len(HIGHLIGHT_IN_BLUE) > 0:
			# Scatter plot for blue proteins
			blue_protein_indices = [monomer_idx_dict[protein] for protein in HIGHLIGHT_IN_BLUE]
			fig.add_trace(go.Scatter(
				x=log_avg_production_rate[blue_protein_indices],
				y=log_avg_loss_rate[blue_protein_indices],
				mode='markers',
				hovertext=HIGHLIGHT_IN_BLUE,
				marker=dict(size=5, color='blue', opacity=1),
				name='Blue Proteins'
			))

		if len(HIGHLIGHT_IN_PURPLE) > 0:
			# Scatter plot for purple proteins
			purple_protein_indices = [monomer_idx_dict[protein] for protein in HIGHLIGHT_IN_PURPLE]
			fig.add_trace(go.Scatter(
				x=log_avg_production_rate[purple_protein_indices],
				y=log_avg_loss_rate[purple_protein_indices],
				mode='markers',
				hovertext=HIGHLIGHT_IN_PURPLE,
				marker=dict(size=5, color='hotpink', opacity=1),
				name='Purple Proteins'
			))
		# Scater plot for purple proteins:

		# Generate line data
		x = np.linspace(-4, 3, 100)

		# y = x line (black)
		fig.add_trace(go.Scatter(
			x=x, y=x, mode='lines', line=dict(color='black', width=2,),
			name='y = x'
		))

		# y = x + 1 line (green, dashed)
		fig.add_trace(go.Scatter(
			x=x, y=x + 1, mode='lines',
			line=dict(color='green', width=2, dash='dash'), name='y = x + 1'
		))

		# y = x - 1 line (green, dashed)
		fig.add_trace(go.Scatter(
			x=x, y=x - 1, mode='lines',
			line=dict(color='green', width=2, dash='dash'), name='y = x - 1'
		))

		# Layout settings
		fig.update_layout(
			title=f"Log10 Average Loss Rate vs Log10 Average Production Rate<br>Simulation ID: {sim_id}, Seed: {metadata['seed']}",
			xaxis_title="Log10 Average Production Rate",
			yaxis_title="Log10 Average Loss Rate",
			width=700, height=700,
			showlegend=True,
		)

		#save the plot:
		seed = metadata['seed']
		plot_name = plotOutFileName +"_"+ sim_id + "_seed_"+ seed +"_red_" + str(HIGHLIGHT_IN_RED) + "_blue_" + str(HIGHLIGHT_IN_BLUE) + "_purple_" + str(HIGHLIGHT_IN_PURPLE) + ".html"
		fig.write_html(os.path.join(plotOutDir, plot_name))
		#exportFigure(plt, plotOutDir, plotOutFileName, metadata)


		# compute the average synthesis rate for each protein:
		avg_synthesis_rate = avg_elongated_counts


if __name__ == '__main__':
	Plot().cli()
