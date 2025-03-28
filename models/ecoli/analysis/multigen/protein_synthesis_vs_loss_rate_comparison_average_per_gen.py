"""
Template for multigen analysis plots
"""

import pickle
import os

from matplotlib import pyplot as plt
# noinspection PyUnresolvedReferences
import numpy as np

from models.ecoli.analysis import multigenAnalysisPlot
from wholecell.analysis.analysis_tools import (exportFigure,
	read_bulk_molecule_counts, read_stacked_bulk_molecules, read_stacked_columns)
from wholecell.io.tablereader import TableReader


class Plot(multigenAnalysisPlot.MultigenAnalysisPlot):
	def do_plot(self, seedOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		with open(simDataFile, 'rb') as f:
			sim_data = pickle.load(f)
		with open(validationDataFile, 'rb') as f:
			validation_data = pickle.load(f)

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

		# Load data
		time = read_stacked_columns(cell_paths, 'Main', 'time')
		(free_monomer_counts,) = read_stacked_bulk_molecules(
			cell_paths, monomerIds)

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
			monomer_counts_at_gen_end = free_monomer_counts[end_gen,
										:]  # get this for each protein

			# find the protein counts at the start of the next generation:
			monomer_counts_at_gen_start = free_monomer_counts[start_gen, :]

			# find the difference between the two:
			protein_counts_removed = monomer_counts_at_gen_end - monomer_counts_at_gen_start
			diluted_counts[i, :] = protein_counts_removed
			diluted_counts_over_time[start_gen,
			:] = protein_counts_removed  # put it at the start of the next gen in terms of time


		# compute how many proteins were removed via degradation over the entire sim length:
		degraded_counts = read_stacked_columns(cell_paths, 'MonomerCounts', "protein_deg_ES1_counts") # take from evolve state counter!

		# find the average degradation for each protein over all generations:
		generation_averages = np.zeros((len(cell_paths), len(monomerIds)))
		for i in range(len(end_generation_indices)):
			# retrieve the generation start and end time:
			end_gen = end_generation_indices[i]
			start_gen = start_generation_indices[i]

			# find the protein counts at all time points in the generation:
			monomer_counts_at_gen = degraded_counts[start_gen:end_gen, :]

			# take the average:
			generation_averages[i, :] = np.mean(monomer_counts_at_gen, axis = 0)

		# find the average degradation for each protein over all generations:
		avg_degraded_counts = np.mean(generation_averages, axis = 0)


		# compute how many proteins were removed via dilution over the entire sim length:
		# todo: make sure this is the right way to compute the average (when you only have one dilution timepoint (becuase there was that one graph I did
		total_diluted_counts = np.sum(diluted_counts, axis = 0)
		avg_diluted_counts = total_diluted_counts / len(time) # divide by the number of timesteps to get the average per timestep

		# compute the average loss rate for each protein:
		# todo: decide if I should include dilution rate or not
		avg_loss_rate = avg_degraded_counts

		# compute how many counts were added via elongation over the entire sim length:
		elongated_counts = read_stacked_columns(cell_paths, 'MonomerCounts', "peptide_elongate_ES1_counts")
		# find the average elongation for each protein over all generations:
		generation_averages = np.zeros((len(cell_paths), len(monomerIds)))
		for i in range(len(end_generation_indices)):
			# retrieve the generation start and end time:
			end_gen = end_generation_indices[i]
			start_gen = start_generation_indices[i]

			# find the protein counts at all time points in the generation:
			monomer_counts_at_gen = elongated_counts[start_gen:end_gen, :]

			# take the average:
			generation_averages[i, :] = np.mean(monomer_counts_at_gen, axis=0)

		# find the average elongation for each protein over all generations:
		avg_elongated_counts = np.mean(generation_averages, axis=0)


		# determine the minimum value:
		# find the smallest nonzero value in the array:
		min_avg_loss_rate_indices = np.nonzero(avg_loss_rate)
		min_avg_loss_rate = avg_loss_rate[min_avg_loss_rate_indices]
		min_avg_elongated_counts_indices = np.nonzero(avg_elongated_counts)
		min_avg_elongated_counts = avg_elongated_counts[min_avg_elongated_counts_indices]

		min_values = [np.min(min_avg_loss_rate), np.min(min_avg_elongated_counts)]
		min_value = np.min(min_values)
		log_factor = min_value * 0.1 # add this to avoid negative and zero log values

		# compute the log of the loss and production rates
		log_avg_loss_rate = np.log10(avg_loss_rate + log_factor)
		log_avg_production_rate = np.log10(avg_elongated_counts + log_factor)

		# plot the loss rate vs the production rate:
		plt.figure()

		from sklearn.metrics import r2_score
		#rsquared = r2_score( true, predicted)
		r_squared = r2_score(log_avg_loss_rate, log_avg_production_rate
							 )

		plt.scatter(log_avg_production_rate, log_avg_loss_rate, s=5, alpha=0.3, color = 'grey', label=f"$R^2$ = " + str(r_squared))

		# add an y=x line:
		x = np.linspace(-4, 3, 100)
		plt.plot(x, x, color = 'black', alpha = 0.5)

		# plot a line at 1 magnitude above the y=x line:
		plt.plot(x, x + 1, color = 'green', alpha = 0.5, linestyle='--')

		# plot a line at 1 magnitude below the y=x line:
		plt.plot(x, x - 1, color = 'green', alpha = 0.5, linestyle='--')




		plt.xlabel("Log10 Average Production Rate")
		plt.ylabel("Log10 Average Loss Rate")
		plt.legend()
		plt.title("Log10 Average Loss Rate vs Log10 Average Production Rate")
		plt.axis('square')
		plt.tight_layout()



		#save the plot:
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)


		# compute the average synthesis rate for each protein:
		avg_synthesis_rate = avg_elongated_counts


if __name__ == '__main__':
	Plot().cli()
