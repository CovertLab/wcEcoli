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

		# Load data
		## Simple stacking functions for data from all cells
		names = ['ATP[c]']  # Replace with desired list of names
		time = read_stacked_columns(cell_paths, 'Main', 'time')
		(counts,) = read_stacked_bulk_molecules(cell_paths, (names,))

		# Extract protein indexes for each new gene
		monomer_counts_reader = TableReader(
			os.path.join(simOutDir, "MonomerCounts"))
		monomer_idx_dict = {monomer: i for i, monomer in
							enumerate(monomer_counts_reader.readAttribute(
								'monomerIds'))}

		bulkMolecules = TableReader(os.path.join(simOutDir, "BulkMolecules"))
		bulkMoleculeCounts = bulkMolecules.readColumn("counts")
		moleculeIds = bulkMolecules.readAttribute("objectNames")

		# other items from bulk molecules:
		bulkMolecules_time = monomer_counts_reader.readColumn(
			"time")  # equal to the right number of seconds? based on the length of the other attributes below
		bulkMolecules_simulationStep = monomer_counts_reader.readColumn("simulationStep")
		bulkMolecules_monomerCounts = monomer_counts_reader.readColumn("monomerCounts")

		# extract the numbers of interest
		monomerIds = monomer_counts_reader.readAttribute(
			"monomerIds")  # this is the one that matches the indexing  I used earlier to construct the listeners!

		# these likely are the actual counts to be degraded at the timestep
		bulkMolecules_pd_CR2_counts = monomer_counts_reader.readColumn(
			"protein_deg_CR2_counts")  # this is likely still zero
		bulkMolecules_pd_ES1_counts = monomer_counts_reader.readColumn(
			"protein_deg_ES1_counts")  # of the two, this is likely the one where counts is updated to the correct number

		# total counts from protein deg:
		bulkMolecules_pd_CR2__TC = monomer_counts_reader.readColumn(
			"protein_deg_CR2__totalCount")  # this is likely still zero
		bulkMolecules_pd_CR1__TC = monomer_counts_reader.readColumn(
			"protein_deg_CR1__totalCount")  # this is likely still zero

		# peptide elongation:
		bulkMolecules_pe_ES1_counts = monomer_counts_reader.readColumn(
			"peptide_elongate_ES1_counts")

		# peptide elgonation total counts:
		bulkMolecules_pe_ES1__TC = monomer_counts_reader.readColumn(
			"peptide_elongate_ES1__totalCount")

		# extract the data over all generations:
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

		print(diluted_counts)

		# compute how many proteins were removed via degradation over the entire sim length:
		degraded_counts = read_stacked_columns(cell_paths, 'MonomerCounts', "protein_deg_ES1_counts") # take from evolve state counter!
		# take the average:
		avg_degraded_counts = np.mean(degraded_counts, axis = 0)

		# compute how many proteins were removed via dilution over the entire sim length:
		# todo: make sure this is the right way to compute the average (when you only have one dilution timepoint (becuase there was that one graph I did
		total_diluted_counts = np.sum(diluted_counts, axis = 0)
		avg_diluted_counts = total_diluted_counts / len(time) # divide by the number of timesteps to get the average per timestep

		# compute the average loss rate for each protein:
		avg_loss_rate = avg_degraded_counts + avg_diluted_counts
		log_avg_loss_rate = np.log10(avg_loss_rate)


		# compute how many counts were added via elongation over the entire sim length:
		elongated_counts = read_stacked_columns(cell_paths, 'MonomerCounts', "peptide_elongate_ES1_counts")
		avg_elongated_counts = np.mean(elongated_counts, axis = 0)
		log_avg_production_rate = np.log10(avg_elongated_counts)

		# plot the loss rate vs the production rate:
		plt.figure()

		from sklearn.metrics import r2_score

		# rsquared = r2_score( true, predicted)
		#r_squared = r2_score(log_avg_production_rate,
							 #log_avg_loss_rate)

		plt.scatter(log_avg_production_rate, log_avg_loss_rate, s=5, alpha=0.3, color = 'grey', )

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
