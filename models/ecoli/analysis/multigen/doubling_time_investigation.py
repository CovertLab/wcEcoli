"""
Plot doubling time, growth rate,
RNAP and ribosome counts, and ppGpp concentration.
"""

import pickle
import os

from matplotlib import pyplot as plt
import matplotlib as mpl
# noinspection PyUnresolvedReferences
import numpy as np
from numpy import inf

from models.ecoli.analysis import multigenAnalysisPlot
from models.ecoli.sim.variants.new_gene_internal_shift import determine_new_gene_ids_and_indices
from wholecell.analysis.analysis_tools import (exportFigure,
	read_stacked_bulk_molecules, read_stacked_columns, read_bulk_molecule_counts)
from wholecell.io.tablereader import TableReader
from wholecell.utils import units

class Plot(multigenAnalysisPlot.MultigenAnalysisPlot):
	def get_mRNA_ids_from_monomer_ids(self, sim_data, target_monomer_ids):
		"""
		Map monomer ids back to the mRNA ids that they were translated from.

		Args:
			target_monomer_ids: ids of the monomers to map to mRNA ids

		Returns: set of mRNA ids
		"""
		# Map protein ids to cistron ids
		monomer_ids = sim_data.process.translation.monomer_data['id']
		cistron_ids = sim_data.process.translation.monomer_data[
			'cistron_id']
		monomer_to_cistron_id_dict = {
			monomer_id: cistron_ids[i] for i, monomer_id in
			enumerate(monomer_ids)}
		target_cistron_ids = [
			monomer_to_cistron_id_dict.get(RNAP_monomer_id) for
			RNAP_monomer_id in target_monomer_ids]
		# Map cistron ids to RNA indexes
		target_RNA_indexes = [
			sim_data.process.transcription.cistron_id_to_rna_indexes(
				RNAP_cistron_id) for RNAP_cistron_id in
			target_cistron_ids]
		# Map RNA indexes to RNA ids
		RNA_ids = sim_data.process.transcription.rna_data['id']
		target_RNA_ids = set()
		for i in range(len(target_RNA_indexes)):
			for index in target_RNA_indexes[i]:
				target_RNA_ids.add(RNA_ids[index])
		return target_RNA_ids

	def get_mRNA_indexes_from_monomer_ids(self, sim_data, all_cells, target_monomer_ids, index_type):
		"""
		Retrieve new gene indexes of a given type.

		Args:
			all_cells: Paths to all cells to read data from (directories should
				contain a simOut/ subdirectory), typically the return from
				AnalysisPaths.get_cells()
			target_monomer_ids: ids of the monomers to map to mRNA indexes
			index_type: Type of indexes to extract, currently supported options
				are 'mRNA' and 'monomer'

		Returns:
			List of requested indexes
		"""
		sim_dir = all_cells[0]
		simOutDir = os.path.join(sim_dir, 'simOut')

		if index_type == 'mRNA':
			# Map protein ids to RNA ids
			target_RNA_ids = self.get_mRNA_ids_from_monomer_ids(
				sim_data, target_monomer_ids)
			# Get index of those RNA ids in the output
			mRNA_counts_reader = TableReader(os.path.join(simOutDir, 'RNACounts'))
			mRNA_idx_dict = {
				rna: i for i, rna in
				enumerate(mRNA_counts_reader.readAttribute('mRNA_ids'))}
			output_indexes = [
				mRNA_idx_dict.get(mRNA_id) for mRNA_id in target_RNA_ids]

		elif index_type == 'monomer':
			# Get index of those monomer ids in the output
			monomer_counts_reader = TableReader(
				os.path.join(simOutDir, "MonomerCounts"))
			monomer_idx_dict = {
				monomer: i for i, monomer in enumerate(
					monomer_counts_reader.readAttribute('monomerIds'))}
			output_indexes = [
				monomer_idx_dict.get(monomer_id) for monomer_id in
				target_monomer_ids]

		else:
			raise Exception(
				"Index type " + index_type +
				" has no instructions for data extraction.")

		return output_indexes

	def do_plot(self, seedOutDir, plotOutDir, plotOutFileName, simDataFile,
				validationDataFile, metadata):
		with open(simDataFile, 'rb') as f:
			sim_data = pickle.load(f)
		with open(validationDataFile, 'rb') as f:
			validation_data = pickle.load(f)

		cell_paths = self.ap.get_cells()
		sim_dir = cell_paths[0]
		simOutDir = os.path.join(sim_dir, 'simOut')

		# Determine new gene ids
		with open(simDataFile, 'rb') as f:
			sim_data = pickle.load(f)

		# Load data
		time = read_stacked_columns(
			cell_paths, 'Main', 'time', ignore_exception=True)
		time_no_first = read_stacked_columns(
			cell_paths, 'Main', 'time', remove_first=True, ignore_exception=True)
		all_mRNA_stacked_counts = read_stacked_columns(
			cell_paths, 'RNACounts', 'mRNA_counts', ignore_exception=True)

		plot_suffixes = ["", "_standard_axes_y"]
		total_plots = 30 # TODO Modularize and get rid of this magic number

		for i in range(len(plot_suffixes)):

			plot_suffix = plot_suffixes[i]

			# Plotting
			mpl.rcParams['axes.spines.right'] = False
			mpl.rcParams['axes.spines.top'] = False
			plt.figure(figsize = (8.5, 60))
			plot_num = 1
			ax1 = plt.subplot(total_plots, 1, plot_num)

			# Growth Rate
			growth_rate = np.ravel(read_stacked_columns(
				cell_paths, "Mass", "instantaneous_growth_rate",
				ignore_exception=True))
			moving_window = min(500, len(growth_rate))
			convolution_array = (np.ones(moving_window) / moving_window)
			growth_rate_convolved = np.convolve(
				convolution_array, growth_rate, mode='same')
			plt.plot(time.flatten() / 60., growth_rate_convolved)
			if plot_suffix == "_standard_axes_both" or plot_suffix == "_standard_axes_y":
				plt.ylim(0,.0004)
			else:
				plt.ylim(bottom=0)
			plt.xlabel("Time (min)")
			plt.ylabel("Growth Rate", fontsize="small")
			plt.title("Growth Rate")
			plot_num += 1

			# Mass
			plt.subplot(total_plots, 1, plot_num, sharex=ax1)
			mass = read_stacked_columns(
				cell_paths, "Mass", "cellMass", ignore_exception=True)
			plt.plot(time / 60., mass)
			if plot_suffix == "_standard_axes_both" or plot_suffix == "_standard_axes_y":
				plt.ylim(0,4000)
			else:
				plt.ylim(bottom=0)
			plt.xlabel("Time (min)")
			plt.ylabel("Mass (fg)", fontsize="small")
			plt.title("Cell Mass")
			plot_num += 1

			# Doubling time
			plt.subplot(total_plots, 1, plot_num, sharex=ax1)
			dt = read_stacked_columns(
				cell_paths, 'Main', 'time',
				fun=lambda x: (x[-1] - x[0]) / 60.).squeeze()
			num_time_steps = read_stacked_columns(
				cell_paths, 'Main', 'time',
				fun=lambda x: len(x)).squeeze()
			# Create a new numpy array where each dt[i] is repeated num_time_steps[i] times
			dt_to_plot = np.repeat(dt, num_time_steps)
			plt.plot(time / 60., dt_to_plot)
			if plot_suffix == "_standard_axes_both" or plot_suffix == "_standard_axes_y":
				plt.ylim((0, 80))
			plt.xlabel("Time (min)")
			plt.ylabel("Doubling Time (min)", fontsize="small")
			plt.title("Doubling Time")
			plot_num += 1

			# Total Ribosome Counts
			plt.subplot(total_plots, 1, plot_num, sharex=ax1)
			# Inactive
			complex_id_30s = [sim_data.molecule_ids.s30_full_complex]
			complex_id_50s = [sim_data.molecule_ids.s50_full_complex]
			(complex_counts_30s, complex_counts_50s) = read_stacked_bulk_molecules(
				cell_paths, (complex_id_30s, complex_id_50s), ignore_exception=True)
			inactive_ribosome_counts = np.minimum(
				complex_counts_30s, complex_counts_50s)
			# Active
			unique_molecule_counts_table = TableReader(
				os.path.join(simOutDir, "UniqueMoleculeCounts"))
			ribosome_index = unique_molecule_counts_table.readAttribute(
				"uniqueMoleculeIds").index('active_ribosome')
			active_ribosome_counts = read_stacked_columns(
				cell_paths, 'UniqueMoleculeCounts',
				'uniqueMoleculeCounts', ignore_exception=True)[:, ribosome_index]
			# Total
			total_ribosome_counts = active_ribosome_counts + inactive_ribosome_counts
			plt.plot(time / 60., total_ribosome_counts)
			if plot_suffix == "_standard_axes_both" or plot_suffix == "_standard_axes_y":
				plt.ylim((0,40000))
			else:
				plt.ylim(bottom=0)
			plt.xlabel("Time (min)")
			plt.ylabel("Total Ribosome Counts", fontsize="small")
			plt.title("Total Ribosome Counts")
			plot_num += 1

			# Total RNAP Counts
			plt.subplot(total_plots, 1, plot_num, sharex=ax1)
			# Inactive
			rnap_id = [sim_data.molecule_ids.full_RNAP]
			(inactive_rnap_counts,) = read_stacked_bulk_molecules(
				cell_paths, (rnap_id,), ignore_exception=True)
			# Active
			uniqueMoleculeCounts = TableReader(
				os.path.join(simOutDir, "UniqueMoleculeCounts"))
			active_rnap_index = uniqueMoleculeCounts.readAttribute(
				"uniqueMoleculeIds").index('active_RNAP')
			active_rnap_counts = read_stacked_columns(
				cell_paths, 'UniqueMoleculeCounts',
				'uniqueMoleculeCounts',
				ignore_exception=True)[:, active_rnap_index]
			# Total
			total_rnap_counts = inactive_rnap_counts + active_rnap_counts
			plt.plot(time / 60., total_rnap_counts)
			if plot_suffix == "_standard_axes_both" or plot_suffix == "_standard_axes_y":
				plt.ylim((0,10000))
			else:
				plt.ylim(bottom=0)
			plt.xlabel("Time (min)")
			plt.ylabel("Total RNAP Counts", fontsize="small")
			plt.title("Total RNAP Counts")
			plot_num += 1

			# ppGpp
			plt.subplot(total_plots, 1, plot_num, sharex=ax1)
			ppGpp_concentration = read_stacked_columns(
				cell_paths, "GrowthLimits", "ppgpp_conc", remove_first=True,
				ignore_exception=True)
			plt.plot(time_no_first / 60., ppGpp_concentration)
			if plot_suffix == "_standard_axes_both" or plot_suffix == "_standard_axes_y":
				plt.ylim((0,150))
			else:
				plt.ylim(bottom=0)
			plt.xlabel("Time (min)")
			plt.ylabel("ppGpp Concentration (uM)", fontsize="small")
			plt.title("ppGpp Concentration")
			plot_num += 1

			# Glucose Comsumption Rate
			# TODO: extend to other carbon sources
			GLUCOSE_ID = "GLC[p]"
			FLUX_UNITS = units.mmol / units.g / units.h
			MASS_UNITS = units.fg
			GROWTH_UNITS = units.fg / units.s
			fba_results = TableReader(os.path.join(simOutDir, "FBAResults"))
			external_molecule_ids = np.array(
				fba_results.readAttribute("externalMoleculeIDs"))
			fba_results.close()
			if GLUCOSE_ID not in external_molecule_ids:
				print("This plot only runs when glucose is the carbon source.")
				return
			glucose_idx = np.where(external_molecule_ids == GLUCOSE_ID)[0][0]
			glucose_flux = FLUX_UNITS * read_stacked_columns(
				cell_paths, "FBAResults", "externalExchangeFluxes",
				ignore_exception=True, remove_first=True)[:, glucose_idx]
			glucose_mw = sim_data.getter.get_mass(GLUCOSE_ID)
			cell_dry_mass = MASS_UNITS * read_stacked_columns(
				cell_paths, "Mass", "dryMass", ignore_exception=True,
				remove_first=True).flatten()
			glucose_mass_flux = glucose_flux * glucose_mw * cell_dry_mass
			plt.subplot(total_plots, 1, plot_num, sharex=ax1)
			moving_window = min(300, len(glucose_mass_flux.asNumber()))
			convolution_array = (np.ones(moving_window) / moving_window)
			glucose_mass_flux_convolved = np.convolve(
				convolution_array, glucose_mass_flux.asNumber(), mode='same')
			plt.plot(time_no_first / 60., -glucose_mass_flux_convolved)
			if plot_suffix == "_standard_axes_both" or plot_suffix == "_standard_axes_y":
				plt.ylim((0,1800))
			plt.xlabel("Time (min)")
			plt.ylabel("Glucose Comsumption Rate (fg/h)", fontsize='x-small')
			plt.title("Glucose Consumption Rate")
			plot_num += 1

			# Active RNAP Portion Allocation
			plt.subplot(total_plots, 1, plot_num, sharex=ax1)
			# Active RNAP Counts
			uniqueMoleculeCounts = TableReader(
				os.path.join(simOutDir, "UniqueMoleculeCounts"))
			active_rnap_index = uniqueMoleculeCounts.readAttribute(
				"uniqueMoleculeIds").index('active_RNAP')
			active_rnap_counts = read_stacked_columns(
				cell_paths, 'UniqueMoleculeCounts',
				'uniqueMoleculeCounts',
				ignore_exception=True)[:, active_rnap_index]
			# rRNA RNAP Portion
			rrna_rnap_counts = np.sum(read_stacked_columns(
				cell_paths, "RNACounts", "partial_rRNA_counts",
				ignore_exception=True), axis = 1).flatten()
			rrna_rnap_portion = rrna_rnap_counts / active_rnap_counts
			# RNAP Subunit RNAP Portion
			RNAP_subunit_monomer_ids = sim_data.molecule_groups.RNAP_subunits
			rnap_subunit_mRNA_indexes = self.get_mRNA_indexes_from_monomer_ids(
				sim_data, cell_paths, RNAP_subunit_monomer_ids, "mRNA")
			rnap_subunit_rnap_counts = np.sum(read_stacked_columns(
				cell_paths, "RNACounts", "partial_mRNA_counts",
				ignore_exception=True)[:, rnap_subunit_mRNA_indexes], axis=1).flatten()
			rnap_subunit_rnap_portion = rnap_subunit_rnap_counts / active_rnap_counts
			# Ribosomal Proteins RNAP Portion
			ribosomal_monomer_ids = sim_data.molecule_groups.ribosomal_proteins
			ribosomal_mRNA_indexes = self.get_mRNA_indexes_from_monomer_ids(
				sim_data, cell_paths, ribosomal_monomer_ids, "mRNA")
			ribosomal_rnap_counts = np.sum(read_stacked_columns(
				cell_paths, "RNACounts", "partial_mRNA_counts",
				ignore_exception=True)[:, ribosomal_mRNA_indexes], axis=1).flatten()
			ribosomal_rnap_portion = ribosomal_rnap_counts / active_rnap_counts
			# Plot
			plt.plot(time / 60., rnap_subunit_rnap_portion, label="RNAP Subunit")
			plt.plot(time / 60., ribosomal_rnap_portion, label="Ribo. Prot.")
			plt.plot(time / 60., rrna_rnap_portion, label="rRNA")
			if plot_suffix == "_standard_axes_both" or plot_suffix == "_standard_axes_y":
				plt.ylim((-0.1, 1.1))
			plt.xlabel("Time (min)")
			plt.ylabel("Portion of Active RNAPs", fontsize='x-small')
			plt.title("Allocation of Active RNAPs")
			plt.legend(fontsize="x-small")
			plot_num += 1

			# Active Ribosome Portion Allocation
			plt.subplot(total_plots, 1, plot_num, sharex=ax1)
			# Active Ribosome Counts
			unique_molecule_counts_table = TableReader(
				os.path.join(simOutDir, "UniqueMoleculeCounts"))
			ribosome_index = unique_molecule_counts_table.readAttribute(
				"uniqueMoleculeIds").index('active_ribosome')
			active_ribosome_counts = read_stacked_columns(
				cell_paths, 'UniqueMoleculeCounts',
				'uniqueMoleculeCounts', ignore_exception=True)[:, ribosome_index]
			# RNAP Subunits Ribosome Portion
			RNAP_subunit_monomer_ids = sim_data.molecule_groups.RNAP_subunits
			rnap_subunit_monomer_indexes = self.get_mRNA_indexes_from_monomer_ids(
				sim_data, cell_paths, RNAP_subunit_monomer_ids, "monomer")
			rnap_subunit_ribosome_counts = np.sum(read_stacked_columns(
				cell_paths, "RibosomeData", "n_ribosomes_per_transcript",
				ignore_exception=True)[:, rnap_subunit_monomer_indexes], axis=1).flatten()
			rnap_subunit_ribosome_portion = rnap_subunit_ribosome_counts / active_ribosome_counts
			# Ribosomal Proteins Ribosome Portion
			ribosomal_monomer_ids = sim_data.molecule_groups.ribosomal_proteins
			ribosomal_monomer_indexes = self.get_mRNA_indexes_from_monomer_ids(
				sim_data, cell_paths, ribosomal_monomer_ids, "monomer")
			ribosomal_ribosome_counts = np.sum(read_stacked_columns(
				cell_paths, "RibosomeData", "n_ribosomes_per_transcript",
				ignore_exception=True)[:, ribosomal_monomer_indexes], axis=1).flatten()
			ribosomal_ribosome_portion = ribosomal_ribosome_counts / active_ribosome_counts
			# Plot
			plt.plot(time / 60., rnap_subunit_ribosome_portion, label="RNAP Subunit")
			plt.plot(time / 60., ribosomal_ribosome_portion, label="Ribo. Prot.")
			if plot_suffix == "_standard_axes_both" or plot_suffix == "_standard_axes_y":
				plt.ylim((-0.1, 1.1))
			plt.xlabel("Time (min)")
			plt.ylabel("Portion of Active Ribosomes", fontsize='x-small')
			plt.title("Allocation of Active Ribosomes")
			plt.legend(fontsize="x-small")
			plot_num += 1

			# Amino acid concentrations
			aa_ids_of_interest = sim_data.molecule_groups.amino_acids
			targets = np.array(
				[sim_data.process.metabolism.conc_dict[key].asNumber(units.mmol / units.L) for key
					in aa_ids_of_interest])
			aa_counts = read_stacked_bulk_molecules(
				cell_paths, aa_ids_of_interest, remove_first=True, ignore_exception=True)[0]
			counts_to_molar = read_stacked_columns(
				cell_paths, 'EnzymeKinetics', 'countsToMolar',
				remove_first=True, ignore_exception=True)
			aa_conc = aa_counts * counts_to_molar

			for ii, aa in enumerate(aa_ids_of_interest):
				plt.subplot(total_plots, 1, plot_num, sharex=ax1)
				plt.plot(time_no_first / 60., aa_conc[:, aa_ids_of_interest.index(aa)])
				if plot_suffix == "_standard_axes_both" or plot_suffix == "_standard_axes_y":
					plt.ylim((0, 0.1))
				plt.axhline(targets[ii], color='r', linestyle='--')
				plt.xlabel("Time (min)")
				plt.ylabel("Concentration (mmol/L)", fontsize='x-small')
				plt.title("Amino Acid Concentration: " + aa)
				plot_num += 1

			plt.subplots_adjust(hspace = 0.7, top = 0.95, bottom = 0.05)
			exportFigure(plt, plotOutDir, plotOutFileName + plot_suffix, metadata)
			plt.close("all")

if __name__ == '__main__':
	Plot().cli()
