"""
Template for multigen analysis plots
"""

import pickle
import os
from scipy.optimize import curve_fit
from matplotlib import pyplot as plt
# noinspection PyUnresolvedReferences
import numpy as np

from models.ecoli.analysis import multigenAnalysisPlot
from wholecell.analysis.analysis_tools import (exportFigure,
	read_bulk_molecule_counts, read_stacked_bulk_molecules, read_stacked_columns)
from wholecell.io.tablereader import TableReader


PLOT_PROTEINS = ['G6890-MONOMER[c]',
					   'PD03938[c]',
						'G6737-MONOMER[c]',
					   	'RPOD-MONOMER[c]',
					   	'PD02936[c]',
					   	'RED-THIOREDOXIN2-MONOMER[c]',
 						"PD03867[c]"]

						#"EG11734-MONOMER[c]",]
#"PD03867[c]"



class Plot(multigenAnalysisPlot.MultigenAnalysisPlot):
	def do_plot(self, seedOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		with open(simDataFile, 'rb') as f:
			sim_data = pickle.load(f)
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
		diluted_counts = np.zeros(((len(doubling_times)-1), len(monomerIds)))
		diluted_counts_over_time = np.zeros(((len(time)), len(monomerIds)))
		for i in range(len(doubling_times) -1): # -1 is to account for not doing the last generation
			end_gen = end_generation_indices[i]
			start_gen = start_generation_indices[i+1] # the first start is zero, so skip that
			print(end_gen, start_gen)

			# find the protein counts at the end of the generation:
			monomer_counts_at_gen_end = free_monomer_counts[end_gen,:] # get this for each protein

			# find the protein counts at the start of the next generation:
			monomer_counts_at_gen_start = free_monomer_counts[start_gen,:]

			# find the difference between the two:
			protein_counts_removed = monomer_counts_at_gen_end - monomer_counts_at_gen_start
			diluted_counts[i,:] = protein_counts_removed
			diluted_counts_over_time[start_gen,:] = protein_counts_removed # put it at the start of the next gen in terms of time

		diluted_counts



		# compute how many proteins were removed via degradation over the entire sim length:
		degraded_counts = read_stacked_columns(cell_paths, 'MonomerCounts', "protein_deg_ES1_counts") # take from evolve state counter!
		# take the average:
		avg_degraded_counts = np.mean(degraded_counts, axis = 0)

		# compute how many proteins were removed via dilution over the entire sim length:
		# todo: make sure this is the right way to compute the average (when you only have one dilution timepoint (becuase there was that one graph I did
		total_diluted_counts = np.sum(diluted_counts, axis = 0)
		avg_diluted_counts = total_diluted_counts / len(time) # divide by the number of timesteps to get the average per timestep

		# compute the average loss rate for each protein:
		#avg_loss_rate = avg_degraded_counts + avg_diluted_counts
		#log_avg_loss_rate = np.log10(avg_loss_rate)


		# compute how many counts were added via elongation over the entire sim length:
		elongated_counts = read_stacked_columns(cell_paths, 'MonomerCounts', "peptide_elongate_ES1_counts")
		avg_elongated_counts = np.mean(elongated_counts, axis = 0)
		log_avg_production_rate = np.log10(avg_elongated_counts)


		# reorganize the data:
		degraded_counts = degraded_counts * -1
		diluted_counts_over_time = diluted_counts_over_time * -1

		# calculate the net change per time:
		net_rate = degraded_counts + diluted_counts_over_time + elongated_counts



		# plot the loss rate and the production rate:
		for protein in PLOT_PROTEINS:
			protein_idx = monomer_idx_dict[protein]
			protein_FMC = free_monomer_counts[:, protein_idx]
			proteins_degraded = (degraded_counts[:, protein_idx])*-1
			monomer_data = sim_data.process.translation.monomer_data[protein_idx]
			deg_rate = monomer_data["deg_rate"]
			measured_HL = (np.log(2) / deg_rate) / 60



			fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, figsize=(10, 6))

			# Monomer Counts plot
			ax1.plot(time, protein_FMC, color='lightseagreen', label='Free Monomer Counts')
			for i in range(len(end_generation_times)):
				dt = end_generation_times[i]
				ax1.axvline(x=dt, linestyle='--', color="yellowgreen")


				# get the effective HL:
				start_time= start_generation_indices[i]
				end_time = end_generation_indices[i]
				gen_time = time[start_time:end_time]
				hi = 5
				gen_protein_counts = protein_FMC[start_time:end_time]
				gen_degraded_counts = proteins_degraded[start_time:end_time]
				# todo: double check that the protein counts are being sliced correctly

				# calc effective deg rate
				non_zero_counts = np.nonzero(gen_protein_counts)[0]
				k_eff = gen_degraded_counts[non_zero_counts]/gen_protein_counts[non_zero_counts]
				avg_k_eff = np.mean(k_eff) # todo: is it fine to have the mean happen here?
				avg_half_life = (np.log(2) / avg_k_eff)/60
				print(avg_half_life) # Pd03938 has a divide by issue error in the line above. how?
				hi = 5 # PD029 had a single non nan


				# make a plot of the effective half lives
				def effective_HL(t, C0, k):
					return C0 * k * t
					#return C0 * np.exp(-k * t)

				k_avg = np.log(2) / (avg_half_life*60) # todo is this an issue units wise agh. does it need to be multipled by 60
				C0_fit = gen_protein_counts
				# todo: fix this! I think gen_time is what should not be in here. it should be like 1 to end_gen_time!
				time_for_graph = gen_time - np.ones(len(gen_time))*gen_time[0]
				y_data = effective_HL(time_for_graph, C0_fit, k_avg)
				hi = 5
				ax1.plot(gen_time, y_data,
						 color='gray', linestyle='--',)

				y_pos = y_data[0] * .8
				ax1.text(dt-time_for_graph[-1]/2, y_pos, f"HL ≈ {avg_half_life:.1f} min", color="black",
						 ha="center", va="bottom", fontsize=8, rotation=0)

				# also plot the measured half life as a comparison:
				k_measured = np.log(2)/(measured_HL*60)
				y_data_measured = effective_HL(time_for_graph, C0_fit, k_measured)
				ax1.plot(gen_time, y_data_measured,
						 color='orange', linestyle=':', alpha=0.5 )

			sim_name = metadata["description"]


			ax1.text(0,protein_FMC.max()*.95, f"measured HL= {measured_HL:.1f} min", color="orange")
			ax1.set_ylabel("Free Monomer Counts")
			ax1.set_title(f"Monomer counts and influx/efflux behavior over time for {protein}\n {sim_name}")

			# rates
			ax2.plot(time, degraded_counts[:, protein_idx], color='red', alpha=0.5,
					 label='Degradation')
			ax2.plot(time, diluted_counts_over_time[:, protein_idx], color='blue', alpha=0.5,
					 label='Dilution')
			ax2.plot(time, elongated_counts[:, protein_idx], color='green', alpha=0.5,
					 label='Elongation')
			ax2.plot(time, net_rate[:, protein_idx], color='black', linestyle='--', alpha=1,
					 linewidth=0.2, label='Net change')

			ax2.set_xlabel("Time (s)")
			ax2.set_ylabel("Rate")
			ax2.set_ylim([-20, 20])
			ax2.legend()

			plt.tight_layout()








			#save the plot:
			file_name = plotOutFileName + "_synthesis_vs_loss_rate_over_time_with_counts_instantanious_"+protein
			exportFigure(plt, plotOutDir, file_name, metadata)


if __name__ == '__main__':
	Plot().cli()
