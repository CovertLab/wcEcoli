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
					   	'RED-THIOREDOXIN2-MONOMER[c]', 'EG10542-MONOMER[c]'] # ["PD00196", "EG11111-MONOMER", "EG11545-MONOMER", "EG10320-MONOMER", "ADHP-MONOMER", "EG10580-MONOMER", "YJCQ-MONOMER", "G6988-MONOMER", "MOTB-FLAGELLAR-MOTOR-STATOR-PROTEIN"]


PLOT_COMPLEXES = ["CPLX0-2881"]

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





		def check_validity_and_get_compartment(protein_list):
			revised_protein_list = []
			for protein in protein_list:
				if "[" in protein:
					protein = protein[:-3] # remove compartment
				if sim_data.getter.is_valid_molecule(protein):
					revised_name = protein + sim_data.getter.get_compartment_tag(protein)
					revised_protein_list.append(revised_name)

			return revised_protein_list

		PLOT_PROTEINS_revised = check_validity_and_get_compartment(PLOT_PROTEINS)
		PLOT_COMPLEXES_revised = check_validity_and_get_compartment(PLOT_COMPLEXES)

		# GENERATE PLOT
		plt.figure(figsize=(14, 10))


		# EXRACT THE FREE MONOMER COUNTS:
		# Extract monomer indexes for each protein of interest
		monomer_counts_reader = TableReader(os.path.join(simOutDir,
														 'MonomerCounts'))
		monomer_idx_dict = {monomer: i for i, monomer in enumerate(
			monomer_counts_reader.readAttribute('monomerIds'))}


		# Load the time data
		time = read_stacked_columns(cell_paths, 'Main', 'time',
									ignore_exception=True)

		# Get the free monomer counts for each protein
		(free_monomer_counts,) = read_stacked_bulk_molecules(
			cell_paths, monomerIds, ignore_exception=True)

		# reshape free_monomer_counts if only one protein is being plotted:
		if len(PLOT_PROTEINS_revised) == 1:
			free_monomer_counts = free_monomer_counts.reshape(-1, 1)





		# get the complex counts:
		for complex in PLOT_COMPLEXES_revised:
			complex_counts = read_stacked_bulk_molecules(cell_paths, [complex])
			# reshape:
			complex_counts = np.array(complex_counts)
			complex_counts = complex_counts.reshape(-1, 1)
			complex_counts_average = complex_counts.mean(axis=0)
			# plot the complex counts
			name = complex + f" (average ~ {complex_counts_average})"
			plt.plot(time, complex_counts, label=name, linestyle=":")





		# plot the loss rate and the production rate:
		for protein in PLOT_PROTEINS_revised:
			protein_idx = monomer_idx_dict[protein]
			protein_FMC = free_monomer_counts[:, protein_idx]
			# get the average counts over all the gens:
			protein_FMC_average = protein_FMC.mean(axis=0)
			# Monomer Counts plot
			name = protein + f" (average ~ {protein_FMC_average})"
			plt.plot(time, protein_FMC, label=name)

		# plot the generation lines
		for i in range(len(end_generation_times)):
			dt = end_generation_times[i]
			plt.axvline(x=dt, linestyle='--', color="yellowgreen")



		sim_name = metadata["description"]


		plt.ylabel("Counts")
		plt.xlabel("Time (s)")
		plt.title(f"Selected Monomer Counts and Complex Counts Over Time \n {sim_name}")


		plt.legend()

		plt.tight_layout()








		#save the plot:
		file_name = plotOutFileName + "_proteins_" + str(PLOT_PROTEINS_revised) +"_complexes_"+str(PLOT_COMPLEXES_revised)+ ".png"
		exportFigure(plt, plotOutDir, file_name, metadata)


if __name__ == '__main__':
	Plot().cli()