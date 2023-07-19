"""
Plot one value per index via heatmap for
new_gene_expression_and_translation_efficiency variant.

Possible Plots:
- Percent of sims that successfully reached a given generation number
- Average doubling time
- Average cell volume, mass, dry cell mass, mRNA mass, protein mass
- Average new gene mRNA count
- Average new gene mRNA mass fraction
- Average new gene NTP mass fraction
- Average new gene protein count
- Average new gene protein mass fraction
- Average new gene initialization rate for RNAP and Ribosomes
- Average fraction of time new gene is overcrowded by RNAP and Ribosomes
- Average number of overcrowded genes for RNAP and Ribosomes
- Average number of ribosomes
- Average number of RNA polymerases
- Average ppGpp concentration
"""

import numpy as np
from matplotlib import pyplot as plt
import math

from wholecell.io.tablereader import TableReader
from models.ecoli.analysis import variantAnalysisPlot
from wholecell.analysis.analysis_tools import exportFigure, \
	read_stacked_columns, read_stacked_bulk_molecules, \
	stacked_cell_identification
from wholecell.analysis.plotting_tools import heatmap
from unum.units import fg

import os.path
import pickle

FONT_SIZE=9

"""
Dashboard Flag
0: Separate Only (Each plot is its own file)
1: Dashboard Only (One file with all plots)
2: Both Dashboard and Separate
"""
DASHBOARD_FLAG = 2

"""
Count number of sims that reach this generation (remember index 7 
corresponds to generation 8)
"""
COUNT_INDEX = 15

"""
Plot data from generations [MIN_CELL_INDEX, MAX_CELL_INDEX)
Note that early generations may not be representative of dynamics 
due to how they are initialized
"""
MIN_CELL_INDEX = 4
MAX_CELL_INDEX = 16

"""
Specify which subset of heatmaps should be made
Completed_gens heatmap is always made, because it is used to
create the other heatmaps, and should not be included here.
The order listed here will be the order of the heatmaps in the
dashboard plot.
"""
HEATMAPS_TO_MAKE_LIST = ["doubling_times_heatmap",
						 "cell_mass_heatmap",
						 "cell_dry_mass_heatmap",
						 "cell_volume_heatmap",
						 "ppgpp_concentration_heatmap",
						 "rnap_crowding_heatmap",
						 "ribosome_crowding_heatmap",
						 "cell_mRNA_mass_heatmap",
						 "cell_protein_mass_heatmap",
						 "rnap_counts_heatmap",
						 "ribosome_counts_heatmap",
						 "new_gene_mRNA_counts_heatmap",
						 "new_gene_monomer_counts_heatmap",
						 "new_gene_rnap_init_rate_heatmap",
						 "new_gene_ribosome_init_rate_heatmap",
						 "new_gene_mRNA_mass_fraction_heatmap",
						 "new_gene_monomer_mass_fraction_heatmap",
						 "new_gene_rnap_time_overcrowded_heatmap",
						 "new_gene_ribosome_time_overcrowded_heatmap",
						 "new_gene_mRNA_NTP_fraction_heatmap",
						 ]

class Plot(variantAnalysisPlot.VariantAnalysisPlot):
	def save_heatmap_data(
			self, h, initial_index, trl_eff_index, exp_index, curr_heatmap_data,
			cell_mask):
		"""
		Applies cell_mask and saves average data value across seeds and generations
		to the appropriate location in data structure for heatmap h.

		Args:
			h: heatmap identifier
			initial_index: 0 for non new gene heatmaps, otherwise the relative
				index of the new gene
			trl_eff_index: New gene translation efficiency value index for this
				variant
			exp_index: New gene expression value index for this variant
			curr_heatmap_data: Data to save
			cell_mask: Should be same size as curr_heatmap_data, typically used
			 	to filter based on generations
		"""
		self.heatmap_data[h][initial_index, trl_eff_index, exp_index] = round(
			np.mean(curr_heatmap_data[cell_mask]),
			self.heatmap_details[h]['num_digits_rounding'])


	# Functions for extracting heatmap data
	def extract_heatmap_data(
			self, all_cells, h, trl_eff_index, exp_index,
			cell_mask):
		"""
		Extracts and saves data associated with heatmap h to heatmap_data for a
		single variant/parameter combination. Extraction is done based on the
		information in heatmap_details (if a standard, non new gene heatmap) or
		using functions designed to handle new genes and other special cases.

		Args:
			all_cells: paths to all cells to read data from (directories should
				contain a simOut/ subdirectory), typically the return from
				AnalysisPaths.get_cells()
			h: heatmap identifier
			trl_eff_index: New gene translation efficiency value index for this
				variant
			exp_index: New gene expression value index for this variant
			cell_mask: Should be same size as curr_heatmap_data, typically used
			 	to filter based on generations
		"""
		if not self.heatmap_details[h]['is_new_gene_heatmap']:
			if not self.heatmap_details[h]['is_nonstandard_data_retrieval']:
				curr_heatmap_data = read_stacked_columns(
					all_cells, self.heatmap_details[h]['data_table'],
					self.heatmap_details[h]['data_column'],
					remove_first = self.heatmap_details[h]['remove_first'],
					fun=self.heatmap_details[h]['function_to_apply'])
				self.save_heatmap_data(
					h, 0, trl_eff_index, exp_index, curr_heatmap_data, cell_mask)
			elif h == "rnap_counts_heatmap":
				self.extract_rnap_counts_data(
					all_cells, h, trl_eff_index, exp_index, cell_mask)
			elif h == "ribosome_counts_heatmap":
				self.extract_ribosome_counts_data(
					all_cells, h, trl_eff_index, exp_index, cell_mask)
			elif h == "rnap_crowding_heatmap":
				self.extract_crowding_heatmap_data(
					all_cells, h, trl_eff_index, exp_index, cell_mask,
					'actual_rna_synth_prob', 'target_rna_synth_prob')
			elif h == "ribosome_crowding_heatmap":
				self.extract_crowding_heatmap_data(
					all_cells, h, trl_eff_index, exp_index, cell_mask,
					'actual_prob_translation_per_transcript',
					'target_prob_translation_per_transcript')
			else:
				raise Exception(
					"Heatmap " + h + " is neither a standard heatmap nor a"
					" nonstandard heatmap that has specific instructions for"
					" data extraction.")

		else: # New gene heatmaps
			if h == "new_gene_mRNA_counts_heatmap":
				self.extract_new_gene_counts_heatmap_data(
					all_cells, h, trl_eff_index,
					exp_index, cell_mask, 'RNACounts', 'mRNA_counts',
					'mRNA')
			elif h == "new_gene_monomer_counts_heatmap":
				self.extract_new_gene_counts_heatmap_data(
					all_cells, h, trl_eff_index,
					exp_index, cell_mask, 'MonomerCounts', 'monomerCounts',
					'monomer')
			elif h == "new_gene_mRNA_mass_fraction_heatmap":
				self.extract_new_gene_mass_fraction_heatmap_data(
					all_cells, h, trl_eff_index, exp_index, cell_mask,
					'RNACounts', 'mRNA_counts', 'Mass', 'mRnaMass',
					'mRNA', self.new_gene_mRNA_ids)
			elif h == "new_gene_monomer_mass_fraction_heatmap":
				self.extract_new_gene_mass_fraction_heatmap_data(
					all_cells, h, trl_eff_index, exp_index, cell_mask,
					'MonomerCounts', 'monomerCounts', 'Mass', 'proteinMass',
					'monomer', self.new_gene_monomer_ids)
			elif h == "new_gene_mRNA_NTP_fraction_heatmap":
				self.extract_new_gene_mRNA_NTP_fraction_heatmap_data(
					all_cells, h, trl_eff_index, exp_index, cell_mask)
			elif h == "new_gene_rnap_init_rate_heatmap":
				self.extract_new_gene_rnap_init_rate_heatmap_data(
					all_cells, h, trl_eff_index, exp_index, cell_mask)
			elif h == "new_gene_ribosome_init_rate_heatmap":
				self.extract_new_gene_ribosome_init_rate_heatmap_data(
					all_cells, h, trl_eff_index, exp_index, cell_mask)
			elif h == "new_gene_rnap_time_overcrowded_heatmap":
				self.extract_new_gene_time_overcrowded_data(
					all_cells, h, trl_eff_index, exp_index, cell_mask,
					'RnaSynthProb', 'tu_is_overcrowded',
					'RNA')
			elif h == "new_gene_ribosome_time_overcrowded_heatmap":
				self.extract_new_gene_time_overcrowded_data(
					all_cells, h, trl_eff_index, exp_index, cell_mask,
					'RibosomeData', 'mRNA_is_overcrowded',
					'monomer')
			else:
				raise Exception(
					"Heatmap " + h + " has no instructions for"
					" data extraction.")

	def extract_rnap_counts_data(
			self, all_cells, h, trl_eff_index, exp_index, cell_mask):
		"""
		Special function to handle extraction and saving of RNAP counts heatmap
			data.

		Args:
			all_cells: paths to all cells to read data from (directories should
				contain a simOut/ subdirectory), typically the return from
				AnalysisPaths.get_cells()
			h: heatmap identifier
			trl_eff_index: New gene translation efficiency value index for this
				variant
			exp_index: New gene expression value index for this variant
			cell_mask: Should be same size as curr_heatmap_data, typically used
				to filter based on generations
		"""
		rnap_id = [self.sim_data.molecule_ids.full_RNAP]
		(rnapCountsBulk,) = read_stacked_bulk_molecules(
			all_cells, (rnap_id,))
		cell_id_vector = stacked_cell_identification(all_cells, 'Main', 'time')
		cell_ids, idx, cell_total_timesteps = np.unique(
			cell_id_vector, return_inverse=True, return_counts=True)
		sum_rnap_counts = np.bincount(idx, weights=rnapCountsBulk)
		avg_rnap_counts = (sum_rnap_counts / cell_total_timesteps)

		self.save_heatmap_data(
			h, 0, trl_eff_index, exp_index, avg_rnap_counts, cell_mask)

	def extract_ribosome_counts_data(
			self, all_cells, h, trl_eff_index, exp_index, cell_mask):
		"""
		Special function to handle extraction and saving of ribosome counts
		heatmap data.

		Args:
			all_cells: paths to all cells to read data from (directories should
				contain a simOut/ subdirectory), typically the return from
				AnalysisPaths.get_cells()
			h: heatmap identifier
			trl_eff_index: New gene translation efficiency value index for this
			 	variant
			exp_index: New gene expression value index for this variant
			cell_mask: Should be same size as curr_heatmap_data, typically used
				to filter based on generations
		"""
		# Determine ribosome index
		sim_dir = all_cells[0]
		simOutDir = os.path.join(sim_dir, 'simOut')
		uniqueMoleculeCounts = TableReader(os.path.join(simOutDir, "UniqueMoleculeCounts"))
		ribosome_index = uniqueMoleculeCounts.readAttribute(
			"uniqueMoleculeIds").index('active_ribosome')

		curr_heatmap_data = read_stacked_columns(
			all_cells, self.heatmap_details[h]['data_table'],
			self.heatmap_details[h]['data_column'],
			remove_first=self.heatmap_details[h]['remove_first'],
			fun=lambda x: np.mean(x[:, ribosome_index], axis=0))
		self.save_heatmap_data(
			h, 0, trl_eff_index, exp_index, curr_heatmap_data, cell_mask)

	def extract_crowding_heatmap_data(
			self, all_cells, h, trl_eff_index, exp_index, cell_mask,
			actual_probs_column, target_probs_column):
		"""
		Special function to handle extraction and saving of RNAP and ribosome
		crowding counts heatmap data.

		Args:
			all_cells: paths to all cells to read data from (directories should
				contain a simOut/ subdirectory), typically the return from
				AnalysisPaths.get_cells()
			h: heatmap identifier
			trl_eff_index: New gene translation efficiency value index for this
				variant
			exp_index: New gene expression value index for this variant
			cell_mask: Should be same size as curr_heatmap_data, typically used
				to filter based on generations
			actual_probs_column: Data column which contains actual
				probabilities of RNA synthesis or translation
			target_probs_column: Data column which contains target
				probabilities of RNA synthesis or translation
		"""
		avg_actual_prob = read_stacked_columns(
			all_cells, self.heatmap_details[h]['data_table'],
			actual_probs_column, fun=lambda x: np.mean(x, axis=0))
		avg_target_prob = read_stacked_columns(
			all_cells, self.heatmap_details[h]['data_table'],
			target_probs_column, fun=lambda x: np.mean(x, axis=0))
		# Get indexes that on average were overcrowded in any generation for
		# any seed
		num_overcrowded_indexes = np.array([len(np.where(sum(
			(avg_actual_prob < avg_target_prob)[cell_mask, :]) > 0)[0])])
		self.save_heatmap_data(
			h, 0, trl_eff_index, exp_index, num_overcrowded_indexes, np.array([True]))

	def get_new_gene_indexes(self, all_cells, index_type):
		"""
		Retrieve new gene indexes of a given type.

		Args:
			all_cells: Paths to all cells to read data from (directories should
				contain a simOut/ subdirectory), typically the return from
				AnalysisPaths.get_cells()
			index_type: Type of indexes to extract, currently supported options
				are 'cistron', 'RNA', 'mRNA', and 'monomer'

		Returns:
			List of requested indexes
		"""
		sim_dir = all_cells[0]
		simOutDir = os.path.join(sim_dir, 'simOut')
		if index_type == 'cistron':
			# Extract cistron indexes for each new gene
			rnap_reader = TableReader(os.path.join(simOutDir, 'RnapData'))
			cistron_idx_dict = {
				cis: i for i, cis in enumerate(rnap_reader.readAttribute('cistron_ids'))}
			new_gene_indexes = [
				cistron_idx_dict.get(mRNA_id) for mRNA_id in self.new_gene_mRNA_ids]
		elif index_type == 'RNA':
			# Extract RNA indexes for each new gene
			rnap_reader = TableReader(os.path.join(simOutDir, 'RnaSynthProb'))
			RNA_idx_dict = {
				rna[:-3]: i for i, rna in enumerate(rnap_reader.readAttribute('rnaIds'))}
			new_gene_indexes = [
				RNA_idx_dict.get(mRNA_id) for mRNA_id in self.new_gene_mRNA_ids]
		elif index_type == 'mRNA':
			# Extract mRNA indexes for each new gene
			mRNA_counts_reader = TableReader(os.path.join(simOutDir, 'RNACounts'))
			mRNA_idx_dict = {
				rna[:-3]: i for i, rna in enumerate(mRNA_counts_reader.readAttribute('mRNA_ids'))}
			new_gene_indexes = [
				mRNA_idx_dict.get(mRNA_id) for mRNA_id in self.new_gene_mRNA_ids]
		elif index_type == 'monomer':
			# Extract protein indexes for each new gene
			monomer_counts_reader = TableReader(os.path.join(simOutDir, "MonomerCounts"))
			monomer_idx_dict = {
				monomer: i for i, monomer in enumerate(
				monomer_counts_reader.readAttribute('monomerIds'))}
			new_gene_indexes = [
				monomer_idx_dict.get(monomer_id) for monomer_id in self.new_gene_monomer_ids]
		else:
			raise Exception(
				"Index type " + index_type + " has no instructions for data extraction.")

		return new_gene_indexes

	def get_avg_new_gene_counts(
			self, all_cells, data_table, data_column, new_gene_indexes):
		"""
		Retreives average counts of new gene mRNAs or proteins, which are needed
		for multiple heatmaps.

		Args:
			all_cells: paths to all cells to read data from (directories should
				contain a simOut/ subdirectory), typically the return from
				AnalysisPaths.get_cells()
			data_table: Table to find data that needs to be retrieved
			data_column: Column to find data that needs to be retreived
			new_gene_indexes: Global indexes of the new genes within data_table

		Returns:
			Average counts of new gene mRNAs or proteins.
		"""
		return (read_stacked_columns(
				all_cells, data_table, data_column, fun=lambda
				x: np.mean( x[:, new_gene_indexes], axis=0)))

	def extract_new_gene_counts_heatmap_data(
			self, all_cells, h, trl_eff_index, exp_index, cell_mask,
			data_table, data_column, new_gene_index_type):
		"""
		Special function to handle extraction and saving of new gene mRNA and
		protein counts heatmap data.

		Args:
			all_cells: paths to all cells to read data from (directories should
				contain a simOut/ subdirectory), typically the return from
				AnalysisPaths.get_cells()
			h: heatmap identifier
			trl_eff_index: New gene translation efficiency value index for this
				variant
			exp_index: New gene expression value index for this variant
			cell_mask: Should be same size as curr_heatmap_data, typically used
				to filter based on generations
			data_table: Table to find data that needs to be retrieved
			data_column: Column to find data that needs to be retreived
			new_gene_index_type: Index type to use for the data table
		"""
		new_gene_indexes = self.get_new_gene_indexes(all_cells, new_gene_index_type)
		avg_new_gene_counts = self.get_avg_new_gene_counts(all_cells, data_table,
			data_column, new_gene_indexes)
		for i in range(len(new_gene_indexes)):
			self.save_heatmap_data(h, i, trl_eff_index, exp_index,
				np.log10(avg_new_gene_counts[:, i] + 1), cell_mask)

	def extract_new_gene_mass_fraction_heatmap_data(
			self, all_cells, h, trl_eff_index, exp_index, cell_mask,
			counts_data_table, counts_data_column, mass_data_table,
			mass_data_column, new_gene_index_type, new_gene_ids):
		"""
		Special function to handle extraction and saving of new gene mRNA and
		protein mass fraction heatmap data.

		Args:
			all_cells: paths to all cells to read data from (directories should
				contain a simOut/ subdirectory), typically the return from
				AnalysisPaths.get_cells()
			h: heatmap identifier
			trl_eff_index: New gene translation efficiency value index for this
				variant
			exp_index: New gene expression value index for this variant
			cell_mask: Should be same size as curr_heatmap_data, typically used
				to filter based on generations
			counts_data_table: Table to find the counts data that needs to be
				retrieved
			counts_data_column: Column to find the counts data that needs to be
				retreived
			mass_data_table: Table to find the mass data that needs to be
				retrieved
			mass_data_column: Column to find the mass data that needs to be
				retreived
			new_gene_index_type: Index type to use for the data table
			new_gene_ids: Ids of new genes in sim_data
		"""
		new_gene_indexes = self.get_new_gene_indexes(all_cells, new_gene_index_type)
		# Get mass for each new gene
		new_gene_masses = [1 for id in new_gene_ids]
		for i in range(len(new_gene_ids)):
			new_gene_masses[i] = (
				self.sim_data.getter.get_mass(
				new_gene_ids[i])/self.sim_data.constants.n_avogadro).asNumber(fg)
		# Get counts for each new gene
		avg_new_gene_counts = self.get_avg_new_gene_counts(
			all_cells, counts_data_table, counts_data_column, new_gene_indexes)
		# Get average mass for all genes
		avg_mass = read_stacked_columns(
					all_cells, mass_data_table, mass_data_column, fun=lambda x: np.mean(x))
		# Determine mass fractions for each new gene
		for i in range(len(self.new_gene_mRNA_ids)):
			new_gene_mass_fraction = (avg_new_gene_counts[:, i] *
				new_gene_masses[i]) / avg_mass[:, i]
			self.save_heatmap_data(
				h, i, trl_eff_index, exp_index, new_gene_mass_fraction,
				cell_mask)

	def extract_new_gene_mRNA_NTP_fraction_heatmap_data(
			self, all_cells, h, trl_eff_index, exp_index, cell_mask):
		"""
		Special function to handle extraction and saving of new gene mRNA and
		NTP fraction heatmap data.

		Args:
			all_cells: paths to all cells to read data from (directories should
				contain a simOut/ subdirectory), typically the return from
				AnalysisPaths.get_cells()
			h: heatmap identifier
			trl_eff_index: New gene translation efficiency value index for this
				variant
			exp_index: New gene expression value index for this variant
			cell_mask: Should be same size as curr_heatmap_data, typically used
				to filter based on generations
		"""
		# Determine number of NTPs per new gene mRNA
		new_gene_mRNA_ntp_counts = [{} for id in self.new_gene_mRNA_ids]
		all_rna_counts_ACGU = self.sim_data.process.transcription.rna_data[
			'counts_ACGU'].asNumber()
		rna_ids = self.sim_data.process.transcription.rna_data['id']
		rna_id_to_index_mapping = {rna[:-3]: i for i, rna in enumerate(rna_ids)}
		for i in range(len(self.new_gene_mRNA_ids)):
			new_gene_mRNA_index = rna_id_to_index_mapping[self.new_gene_mRNA_ids[i]]
			for ntp_index in range(len(self.ntp_ids)):
				new_gene_mRNA_ntp_counts[i][self.ntp_ids[ntp_index]] = (
					all_rna_counts_ACGU[new_gene_mRNA_index, ntp_index])
		# All mRNA
		all_mRNA_counts_ACGU = (
			all_rna_counts_ACGU[self.sim_data.process.transcription.rna_data["is_mRNA"]])
		avg_mRNA_counts = read_stacked_columns(
			all_cells, 'RNACounts', 'mRNA_counts', fun=lambda x: np.mean(x, axis=0))
		# New gene mRNA
		new_gene_mRNA_indexes = self.get_new_gene_indexes(all_cells, 'mRNA')
		avg_new_gene_mRNA_counts = self.get_avg_new_gene_counts(
			all_cells, 'RNACounts', 'mRNA_counts', new_gene_mRNA_indexes)
		# Compute new gene NTP fractions
		all_mRNA_ntp_totals = {}
		for ntp_index in range(len(self.ntp_ids)):
			ntp_id = self.ntp_ids[ntp_index]
			all_mRNA_ntp_totals[ntp_id] = \
				(avg_mRNA_counts @ all_mRNA_counts_ACGU[:, ntp_index])
			for i in range(len(self.new_gene_mRNA_ids)):
				self.heatmap_data[h][ntp_id][i, trl_eff_index, exp_index] = round(
					np.mean((avg_new_gene_mRNA_counts[:, i][cell_mask] *
					new_gene_mRNA_ntp_counts[i][ntp_id]) / all_mRNA_ntp_totals[
					ntp_id][cell_mask]),
					self.heatmap_details[h]['num_digits_rounding'])

	def extract_new_gene_rnap_init_rate_heatmap_data(
			self, all_cells, h, trl_eff_index, exp_index, cell_mask):
		"""
		Special function to handle extraction and saving of RNAP new gene
		initialization rate heatmap data.

		Args:
			all_cells: paths to all cells to read data from (directories should
				contain a simOut/ subdirectory), typically the return from
				AnalysisPaths.get_cells()
			h: heatmap identifier
			trl_eff_index: New gene translation efficiency value index for this
				variant
			exp_index: New gene expression value index for this variant
			cell_mask: Should be same size as curr_heatmap_data, typically used
				to filter based on generations
		"""
		new_gene_cistron_indexes = self.get_new_gene_indexes(all_cells, 'cistron')
		avg_new_gene_copy_number = (read_stacked_columns(
			all_cells, 'RnaSynthProb', 'gene_copy_number',
			fun=lambda x: np.mean(x[:, new_gene_cistron_indexes], axis=0)))
		avg_new_gene_rnap_init_rates = (read_stacked_columns(
			all_cells, 'RnapData', 'rna_init_event_per_cistron',
			fun=lambda x: np.mean(x[:, new_gene_cistron_indexes],
			axis=0))) / avg_new_gene_copy_number
		for i in range(len(self.new_gene_mRNA_ids)):
			self.save_heatmap_data(
				h, i, trl_eff_index, exp_index,
				avg_new_gene_rnap_init_rates[:, i], cell_mask)

	def extract_new_gene_ribosome_init_rate_heatmap_data(
			self, all_cells, h, trl_eff_index, exp_index, cell_mask):
		"""
		Special function to handle extraction and saving of ribosome new gene
		initialization rate heatmap data.

		Args:
			all_cells: paths to all cells to read data from (directories should
				contain a simOut/ subdirectory), typically the return from
				AnalysisPaths.get_cells()
			h: heatmap identifier
			trl_eff_index: New gene translation efficiency value index for this
				variant
			exp_index: New gene expression value index for this variant
			cell_mask: Should be same size as curr_heatmap_data, typically used
				to filter based on generations
		"""
		new_gene_mRNA_indexes = self.get_new_gene_indexes(all_cells, 'mRNA')
		new_gene_monomer_indexes = self.get_new_gene_indexes(all_cells, 'monomer')
		avg_new_gene_mRNA_counts = self.get_avg_new_gene_counts(
			all_cells, 'RNACounts', 'mRNA_counts', new_gene_mRNA_indexes)
		avg_new_gene_ribosome_init_rates = (read_stacked_columns(
			all_cells, 'RibosomeData', 'ribosome_init_event_per_monomer',
			fun=lambda x: np.mean(x[:, new_gene_monomer_indexes],
			axis=0))) / avg_new_gene_mRNA_counts
		for i in range(len(self.new_gene_mRNA_ids)):
			self.save_heatmap_data(
				h, i, trl_eff_index, exp_index,
				avg_new_gene_ribosome_init_rates[:, i], cell_mask)

	def extract_new_gene_time_overcrowded_data(
			self, all_cells, h, trl_eff_index, exp_index, cell_mask,
			data_table, data_column, new_gene_index_type):
		"""
		Special function to handle extraction and saving of RNAP and ribosome
		new gene time overcrowded heatmap data.

		Args:
			all_cells: paths to all cells to read data from (directories should
				contain a simOut/ subdirectory), typically the return from
				AnalysisPaths.get_cells()
			h: heatmap identifier
			trl_eff_index: New gene translation efficiency value index for this
				variant
			exp_index: New gene expression value index for this variant
			cell_mask: Should be same size as curr_heatmap_data, typically used
				to filter based on generations
			data_table: Table to find data that needs to be retrieved
			data_column: Column to find data that needs to be retreived
			new_gene_index_type: Index type to use for the data table
		"""
		new_gene_indexes = self.get_new_gene_indexes(all_cells, new_gene_index_type)
		# Average fraction of time steps that overcrowding occurs for new genes
		# per generation
		new_gene_num_time_steps_overcrowded = (read_stacked_columns(
			all_cells, data_table, data_column,
			fun=lambda x: np.sum(x[:, new_gene_indexes], axis=0)/ (
			x[:, new_gene_indexes].shape[0])))
		for i in range(len(new_gene_indexes)):
			self.save_heatmap_data(
				h, i, trl_eff_index, exp_index,
				new_gene_num_time_steps_overcrowded[:,i], cell_mask)


	# Functions for plotting heatmaps
	def plot_heatmaps(
			self, is_dashboard, variant_mask, heatmap_x_label, heatmap_y_label,
			new_gene_expression_factors, new_gene_translation_efficiency_values,
			figsize_x, figsize_y, plotOutDir, plot_suffix):
		"""
		Plots all heatmaps in order given by HEATMAPS_TO_MAKE_LIST.

		Args:
			is_dashboard: Boolean flag for whether we are creating a dashboard
				of heatmaps or a number of individual heatmaps
			variant_mask: np.array of dimension
				(len(new_gene_translation_efficiency_values),
				len(new_gene_expression_factors)) with entries set to True if
				variant was run, False otherwise.
			heatmap_x_label: Label for x axis of heatmap
			heatmap_y_label: Label for y axis of heatmap
			new_gene_expression_factors: New gene expression factors used in
				these variants
			new_gene_translation_efficiency_values: New gene translation
				efficiency values used in these variants
			figsize_x: Horizontal size of each heatmap
			figsize_y: Vertical size of each heatmap
			plotOutDir: Output directory for plots
			plot_suffix: Suffix to add to plot file names, usually specifying
				which generations were plotted
		"""
		if is_dashboard:
			# Determine dashboard layout
			if len(HEATMAPS_TO_MAKE_LIST) > 3:
				dashboard_ncols = 4
				dashboard_nrows = math.ceil((len(HEATMAPS_TO_MAKE_LIST) + 1) / dashboard_ncols)
			else:
				dashboard_ncols = len(HEATMAPS_TO_MAKE_LIST) + 1
				dashboard_nrows = 1
			fig, axs = plt.subplots(nrows=dashboard_nrows,
				ncols=dashboard_ncols,
				figsize=(figsize_y * dashboard_ncols,figsize_x * dashboard_nrows),
				layout='constrained'
				)
			if dashboard_nrows == 1:
				axs = np.reshape(axs, (1, dashboard_ncols))

			# Percent Completion Heatmap
			heatmap(
				self, axs[0,0], variant_mask,
				self.heatmap_data["completed_gens_heatmap"][0,:,:],
				self.heatmap_data["completed_gens_heatmap"][0,:,:],
				new_gene_expression_factors,
				new_gene_translation_efficiency_values,
				heatmap_x_label,
				heatmap_y_label,
				"Percentage of Sims That Reached Generation " + str(COUNT_INDEX + 1))
			row_ax = 0
			col_ax = 1

			for h in HEATMAPS_TO_MAKE_LIST:
				if not self.heatmap_details[h]["is_nonstandard_plot"]:
					stop_index = 1
					title_addition = ""
					if self.heatmap_details[h]["is_new_gene_heatmap"]:
						stop_index = len(self.new_gene_mRNA_ids)
					for i in range(stop_index):
						if self.heatmap_details[h]["is_new_gene_heatmap"]:
							title_addition = ": " + self.new_gene_mRNA_ids[i][:-4]
						self.make_single_heatmap(
							h, axs[row_ax, col_ax], variant_mask,
							heatmap_x_label, heatmap_y_label, i,
							new_gene_expression_factors,
							new_gene_translation_efficiency_values,
							title_addition)
						col_ax += 1
						if (col_ax == dashboard_ncols):
							col_ax = 0
							row_ax += 1
				elif h == "new_gene_mRNA_NTP_fraction_heatmap":
					for i in range(len(self.new_gene_mRNA_ids)):
						for ntp_id in self.ntp_ids:
							self.make_new_gene_mRNA_NTP_fraction_heatmap(
								h, axs[row_ax, col_ax], variant_mask,
								heatmap_x_label, heatmap_y_label, i,
								new_gene_expression_factors,
								new_gene_translation_efficiency_values, ntp_id)
							fig.tight_layout()
							plt.show()
							exportFigure(
								plt, plotOutDir,
								'new_gene_mRNA_' + ntp_id[:-3] + '_fraction_heatmap'
								+ "_" + self.new_gene_mRNA_ids[i][:-4] +
								plot_suffix)
							col_ax += 1
							if (col_ax == dashboard_ncols):
								col_ax = 0
								row_ax += 1
				else:
					raise Exception(
						"Heatmap " + h + " is neither a standard plot nor a"
						" nonstandard plot that has specific instructions for"
						" plotting.")
			fig.tight_layout()
			exportFigure(plt, plotOutDir,
				"new_gene_exp_trl_eff_dashboard" + plot_suffix)
			plt.close("all")

		else: # individual plots
			# Plot percent completion heatmap
			fig, ax = plt.subplots(1, 1, figsize=(figsize_x, figsize_y))
			heatmap(
				self, ax, variant_mask,
				self.heatmap_data["completed_gens_heatmap"][0, :, :],
				self.heatmap_data["completed_gens_heatmap"][0, :, :],
				new_gene_expression_factors,
				new_gene_translation_efficiency_values,
				heatmap_x_label,
				heatmap_y_label,
				"Percentage of Sims that Reached Generation " \
					+ str(COUNT_INDEX + 1))
			fig.tight_layout()
			plt.show()
			exportFigure(plt, plotOutDir, 'completed_gens_heatmap')

			for h in HEATMAPS_TO_MAKE_LIST:
				if not self.heatmap_details[h]["is_nonstandard_plot"]:
					stop_index = 1
					title_addition = ""
					filename_addition = ""
					if self.heatmap_details[h]["is_new_gene_heatmap"]:
						stop_index = len(self.new_gene_mRNA_ids)
					for i in range(stop_index):
						if self.heatmap_details[h]["is_new_gene_heatmap"]:
							title_addition = ": " + self.new_gene_mRNA_ids[i][:-4]
							filename_addition = "_" + self.new_gene_mRNA_ids[i][:-4]
						fig, ax = plt.subplots(1, 1, figsize=(figsize_x, figsize_y))
						self.make_single_heatmap(
							h, ax, variant_mask, heatmap_x_label, heatmap_y_label,
							i, new_gene_expression_factors,
							new_gene_translation_efficiency_values, title_addition)
						fig.tight_layout()
						plt.show()
						exportFigure(plt, plotOutDir, h + filename_addition +
							plot_suffix)
						plt.close()
				elif h == "new_gene_mRNA_NTP_fraction_heatmap":
					for i in range(len(self.new_gene_mRNA_ids)):
						for ntp_id in self.ntp_ids:
							fig, ax = plt.subplots(1, 1, figsize=(figsize_x,
								figsize_y))
							self.make_new_gene_mRNA_NTP_fraction_heatmap(
								h, ax, variant_mask, heatmap_x_label,
								heatmap_y_label, i, new_gene_expression_factors,
								new_gene_translation_efficiency_values, ntp_id)
							fig.tight_layout()
							plt.show()
							exportFigure(
								plt, plotOutDir,
								'new_gene_mRNA_' + ntp_id[:-3] + '_fraction_heatmap'
								+ "_" + self.new_gene_mRNA_ids[i][:-4] +
								plot_suffix)
				else:
					raise Exception(
						"Heatmap " + h + " is neither a standard plot nor a"
						" nonstandard plot that has specific instructions for"
						" plotting.")

	def make_single_heatmap(
			self, h, ax, variant_mask, heatmap_x_label, heatmap_y_label,
			initial_index, new_gene_expression_factors,
			new_gene_translation_efficiency_values, title_addition):
		"""
		Creates a heatmap for h.

		Args:
			h: Heatmap identifier
			ax: Axes to plot on
			variant_mask: np.array of dimension
				(len(new_gene_translation_efficiency_values),
				len(new_gene_expression_factors)) with entries set to True if
				variant was run, False otherwise.
			heatmap_x_label: Label for x axis of heatmap
			heatmap_y_label: Label for y axis of heatmap
			initial_index: 0 for non new gene heatmaps, otherwise the relative
				index of the new gene
			new_gene_expression_factors: New gene expression factors used in
				these variants
			new_gene_translation_efficiency_values: New gene translation
				efficiency values used in these variants
			title_addition: Any string that needs to be added to the title of
				the heatmap, e.g. a new gene id
		"""
		heatmap(
			self, ax, variant_mask,
			self.heatmap_data[h][initial_index, :, :],
			self.heatmap_data["completed_gens_heatmap"][0, :, :],
			new_gene_expression_factors, new_gene_translation_efficiency_values,
			heatmap_x_label, heatmap_y_label,
			self.heatmap_details[h]['plot_title'] + title_addition,
			self.heatmap_details[h]['box_text_size'])

	def make_new_gene_mRNA_NTP_fraction_heatmap(
			self, h, ax, variant_mask, heatmap_x_label, heatmap_y_label,
			initial_index, new_gene_expression_factors,
			new_gene_translation_efficiency_values, ntp_id):
		"""
		Special function that creates a new gene mRNA NTP fraction heatmap for
		one new gene and one NTP.

		Args:
			h: Heatmap identifier
			ax: Axes to plot on
			variant_mask: np.array of dimension
				(len(new_gene_translation_efficiency_values),
				len(new_gene_expression_factors)) with entries set to True if
				variant was run, False otherwise.
			heatmap_x_label: Label for x axis of heatmap
			heatmap_y_label: Label for y axis of heatmap
			initial_index: 0 for non new gene heatmaps, otherwise the relative
				index of the new gene
			new_gene_expression_factors: New gene expression factors used in
				these variants
			new_gene_translation_efficiency_values: New gene translation
				efficiency values used in these variants
			ntp_id: Id of NTP to plot
		"""
		heatmap(
			self, ax, variant_mask,
			self.heatmap_data[h][ntp_id][initial_index, :, :],
			self.heatmap_data["completed_gens_heatmap"][0, :, :],
			new_gene_expression_factors, new_gene_translation_efficiency_values,
			heatmap_x_label, heatmap_y_label,
			self.heatmap_details[h]['plot_title'] + " " + ntp_id[:-3] +
				" Fraction: " + self.new_gene_mRNA_ids[initial_index][:-4],
			self.heatmap_details[h]['box_text_size'])


	def do_plot(self, inputDir, plotOutDir, plotOutFileName, simDataFile,
				validationDataFile, metadata):
		heatmaps_to_make = set(HEATMAPS_TO_MAKE_LIST)
		with open(simDataFile, 'rb') as f:
			self.sim_data = pickle.load(f)

		"""
		Details needed to create all possible heatmaps
			key (string): Heatmap identifier, will also be used in file name if
				plots are saved separately
			is_new_gene_heatmap (bool): If True, one heatmap will be made 
				for each new gene
			is_nonstandard_data_retrieval (bool): False if only one column needs
			 	to be read from one table to extract data for a non new gene
			 	heatmap. True in all other cases.
			is_nonstandard_plotting (bool): False if only one plot (or one plot
				per new gene) needs to be made. True in all other cases.
			data_table (string): Table to get data from.
			data_column (string): Column in table to get data from.
			default_value (int): Value to use in heatmap if no data is 
				extracted for this parameter combination.
			remove_first (bool): If True, removes the first column of data 
				from each cell (which might be set to a default value in 
				some cases)
			function_to_apply (lambda): Function to apply to data in each generation
				(eg. np.mean will return and array with the mean value for 
				each generation instead of each time point)
			num_digits_rounding (int): Specifies how to round the number 
				displayed in each sequare of the heatmap
			box_text_size (string): Specifies font size of number displayed 
				in each square of the heatmap
			plot_title (string): Title of heatmap to display
		"""
		# Defaults - unless otherwise specified, these values will be
		# used for plotting
		default_is_nonstandard_data_retrieval = False
		default_is_nonstandard_plot = False
		default_value = -1
		default_remove_first = False
		default_function_to_apply = lambda x: np.mean(x)
		default_num_digits_rounding = 2
		default_box_text_size = 'medium'
		# Specify unique fields and non-default values here
		self.heatmap_details = {
			"doubling_times_heatmap" :
				{'data_table': 'Main',
				 'data_column': 'time',
				 'function_to_apply': lambda x: (x[-1] - x[0]) / 60.,
				 'num_digits_rounding': 0,
				 'plot_title': 'Doubling Time (minutes)',
				 },
			"cell_volume_heatmap":
				{'data_table': 'Mass',
				 'data_column': 'cellVolume',
				 'plot_title': 'Cell Volume (fL)',
				 },
			"cell_mass_heatmap":
				{'data_table': 'Mass',
				 'data_column': 'cellMass',
				 'num_digits_rounding': 0,
				 'box_text_size': 'x-small',
				 'plot_title': 'Cell Mass (fg)',
				 },
			"cell_dry_mass_heatmap":
				{'data_table': 'Mass',
				 'data_column': 'dryMass',
				 'num_digits_rounding': 0,
				 'box_text_size': 'x-small',
				 'plot_title': 'Dry Cell Mass (fg)',
				 },
			"cell_mRNA_mass_heatmap":
				{'data_table': 'Mass',
				 'data_column': 'mRnaMass',
				 'plot_title': 'Total mRNA Mass (fg)',
				 },
			"cell_protein_mass_heatmap":
				{'data_table': 'Mass',
				 'data_column': 'proteinMass',
				 'num_digits_rounding': 0,
				 'box_text_size': 'x-small',
				 'plot_title': 'Total Protein Mass (fg)',
				 },
			"ppgpp_concentration_heatmap":
				{'data_table': 'GrowthLimits',
				 'data_column': 'ppgpp_conc',
				 'remove_first': True,
				 'num_digits_rounding': 1,
				 'plot_title': 'ppGpp Concentration (uM)',
				 },
			"rnap_counts_heatmap":
				{'is_nonstandard_data_retrieval': True,
				 'num_digits_rounding': 0,
				 'box_text_size': 'x-small',
				 'plot_title': 'RNA Polymerase (RNAP) Counts',
				 },
			"ribosome_counts_heatmap":
				{'is_nonstandard_data_retrieval': True,
				 'data_table': 'UniqueMoleculeCounts',
				 'data_column': 'uniqueMoleculeCounts',
				 'num_digits_rounding': 0,
				 'box_text_size': 'x-small',
				 'plot_title': 'Ribosome Counts',
				 },
			"rnap_crowding_heatmap":
				{'is_nonstandard_data_retrieval': True,
				 'data_table': 'RnaSynthProb',
				 'function_to_apply': lambda x: np.mean(x, axis = 0),
				 'num_digits_rounding': 0,
				 'plot_title': 'RNAP Crowding: # of TUs',
				 },
			"ribosome_crowding_heatmap":
				{'is_nonstandard_data_retrieval': True,
				 'data_table': 'RibosomeData',
				 'function_to_apply': lambda x: np.mean(x, axis = 0),
				 'num_digits_rounding': 0,
				 'plot_title': 'Ribosome Crowding: # of Monomers',
				 },
			"new_gene_mRNA_counts_heatmap":
				{'plot_title': 'Log(New Gene mRNA Counts+1)'},
			"new_gene_monomer_counts_heatmap":
				{'plot_title': 'Log(New Gene Protein Counts+1)'},
			"new_gene_mRNA_mass_fraction_heatmap":
				{'plot_title': 'New Gene mRNA Mass Fraction'},
			"new_gene_mRNA_NTP_fraction_heatmap":
				{'is_nonstandard_plot': True,
				 'num_digits_rounding': 4,
				 'box_text_size': 'x-small',
				 'plot_title': 'New Gene',
				 },
			"new_gene_monomer_mass_fraction_heatmap":
				{'plot_title': 'New Gene Protein Mass Fraction'},
			"new_gene_rnap_init_rate_heatmap":
				{'plot_title': 'New Gene RNAP Initialization Rate'},
			"new_gene_ribosome_init_rate_heatmap":
				{'plot_title': 'New Gene Ribosome Initalization Rate'},
			"new_gene_rnap_time_overcrowded_heatmap":
				{'plot_title': 'Fraction of Time RNAP Overcrowded New Gene'},
			"new_gene_ribosome_time_overcrowded_heatmap":
				{'plot_title': 'Fraction of Time Ribosome Overcrowded New Gene'},
		}
		assert "completed_gens_heatmap" not in heatmaps_to_make, \
			"the completed_gens_heatmap is run by default, do not include in heatmaps_to_make"
		# Check validity of requested heatmaps and fill in default values where needed
		for h in heatmaps_to_make:
			assert h in self.heatmap_details, "Heatmap " + h + " is not an option"
			self.heatmap_details[h]['is_new_gene_heatmap'] = h.startswith("new_gene_")
			self.heatmap_details[h].setdefault(
				'is_nonstandard_data_retrieval',
				default_is_nonstandard_data_retrieval)
			self.heatmap_details[h].setdefault(
				'is_nonstandard_plot',default_is_nonstandard_plot)
			self.heatmap_details[h].setdefault(
				'default_value', default_value)
			self.heatmap_details[h].setdefault(
				'remove_first', default_remove_first)
			self.heatmap_details[h].setdefault(
				'function_to_apply', default_function_to_apply)
			self.heatmap_details[h].setdefault(
				'num_digits_rounding', default_num_digits_rounding)
			self.heatmap_details[h].setdefault(
				'box_text_size', default_box_text_size)

		# Map variant indices to expression factors and translation efficiency
		# values
		if 'new_gene_expression_factors' not in metadata or \
				'new_gene_translation_efficiency_values' not in metadata:
			print("This plot is intended to be run on simulations where the"
				  " new gene expression-translation efficiency variant was "
				  "enabled, but no parameters for this variant were found.")
			return
		new_gene_expression_factors = metadata['new_gene_expression_factors']
		new_gene_translation_efficiency_values = metadata[
			'new_gene_translation_efficiency_values']
		separator = len(new_gene_translation_efficiency_values)
		variants = self.ap.get_variants()
		variant_index_to_values = {}
		variant_index_to_list_indices = {}
		variant_mask = np.zeros(( # Track whether we ran this sim
			len(new_gene_translation_efficiency_values),
			len(new_gene_expression_factors)), dtype=bool)
		for index in variants:
			if index == 0:
				expression_list_index = 0
				trl_eff_list_index = len(
					new_gene_translation_efficiency_values) - 1
				expression_variant_index = 0
				# Note: this value should not matter since gene is knocked out
				trl_eff_value = 0
			else:
				trl_eff_list_index = index % separator
				if trl_eff_list_index == 0:
					expression_list_index = index // separator
				else:
					expression_list_index = index // separator + 1

				expression_variant_index = new_gene_expression_factors[expression_list_index]
				trl_eff_value = new_gene_translation_efficiency_values[trl_eff_list_index]
			variant_index_to_values[index] = np.array([
				expression_variant_index, trl_eff_value])
			variant_index_to_list_indices[index] = np.array([
				expression_list_index, trl_eff_list_index])
			variant_mask[trl_eff_list_index, expression_list_index] = True

		# Determine new gene ids
		mRNA_sim_data = self.sim_data.process.transcription.cistron_data.struct_array
		monomer_sim_data = self.sim_data.process.translation.monomer_data.struct_array
		self.new_gene_mRNA_ids = mRNA_sim_data[mRNA_sim_data['is_new_gene']]['id'].tolist()
		mRNA_monomer_id_dict = dict(
			zip(monomer_sim_data['cistron_id'], monomer_sim_data['id']))
		self.new_gene_monomer_ids = [
			mRNA_monomer_id_dict.get(mRNA_id) for mRNA_id in self.new_gene_mRNA_ids]
		if len(self.new_gene_mRNA_ids) == 0:
			print("This plot is intended to be run on simulations where the"
				  " new gene option was enabled, but no new gene mRNAs were "
				  "found.")
			return
		if len(self.new_gene_monomer_ids) == 0:
			print("This plot is intended to be run on simulations where the "
				  "new gene option was enabled, but no new gene proteins "
				  "were found.")
			return
		assert len(self.new_gene_monomer_ids) == len(self.new_gene_mRNA_ids),\
			'number of new gene monomers and mRNAs should be equal'

		# Create data structures that to use for the heatmaps
		self.heatmap_data = {}
		self.heatmap_data["completed_gens_heatmap"] = np.zeros((
			1, len(new_gene_translation_efficiency_values),
			len(new_gene_expression_factors)))
		for h in heatmaps_to_make:
			if not self.heatmap_details[h]['is_new_gene_heatmap']:
				self.heatmap_data[h] = np.zeros((
					1, len(new_gene_translation_efficiency_values),
					len(new_gene_expression_factors))
					) + self.heatmap_details[h]['default_value']
			else:
				if h == "new_gene_mRNA_NTP_fraction_heatmap":
					self.ntp_ids = list(
						self.sim_data.ntp_code_to_id_ordered.values())
					self.heatmap_data[
						"new_gene_mRNA_NTP_fraction_heatmap"] = {}
					for ntp_id in self.ntp_ids:
						self.heatmap_data["new_gene_mRNA_NTP_fraction_heatmap"][
							ntp_id] = np.zeros((len(self.new_gene_mRNA_ids),
							len(new_gene_translation_efficiency_values),
							len(new_gene_expression_factors))
							) + self.heatmap_details[h]['default_value']
				else:
					self.heatmap_data[h] = np.zeros((
						len(self.new_gene_mRNA_ids),
						len(new_gene_translation_efficiency_values),
						len(new_gene_expression_factors))
						) + self.heatmap_details[h]['default_value']

		# Data extraction
		print("---Data Extraction---")
		reached_count_gen = {}
		generations = {}

		variants = self.ap.get_variants()
		min_variant = min(variants)
		for variant in variants:

			print("Variant: ",variant)
			all_cells = self.ap.get_cells(variant=[variant], only_successful=True)
			exp_index, trl_eff_index = variant_index_to_list_indices[variant]
			if len(all_cells) == 0:
				continue

			# Determine mask to apply based on generations
			all_cells_gens = np.array([
				int(os.path.basename(os.path.dirname(cell_path))[-6:])
				for cell_path in all_cells])
			generations[variant] = all_cells_gens
			cell_mask = np.logical_and(
				(generations[variant] >=MIN_CELL_INDEX),
				(generations[variant] < MAX_CELL_INDEX))
			if len(cell_mask) == 1:
				cell_mask = cell_mask.reshape(1)
			if sum(cell_mask) < 1:
				continue

			# Completed Gens Heatmap: Count the number of simulations that
			# reach gen COUNT_INDEX + 1
			num_count_gen = len(self.ap.get_cells(
				variant=[variant], generation=[COUNT_INDEX], only_successful=True))
			num_zero_gen = len(self.ap.get_cells(
				variant=[variant], generation=[0], only_successful=True))
			reached_count_gen[variant] = num_count_gen / num_zero_gen
			self.heatmap_data["completed_gens_heatmap"][0, trl_eff_index,
				exp_index]	= round(reached_count_gen[variant], 2)

			# Extract data for each heatmap
			for h in heatmaps_to_make:
				self.extract_heatmap_data(
					all_cells, h, trl_eff_index, exp_index, cell_mask)

		# Plotting
		print("---Plotting---")
		plot_suffix = "_gens_" + str(MIN_CELL_INDEX) + "_through_" + str(MAX_CELL_INDEX)
		heatmap_x_label = "Expression Variant"
		heatmap_y_label = "Translation Efficiency Value"
		figsize_x =  2 + len(new_gene_expression_factors)/2
		figsize_y = 0.5 + len(new_gene_translation_efficiency_values)/2

		# Create dashboard plot
		if DASHBOARD_FLAG == 1 or DASHBOARD_FLAG == 2:
			self.plot_heatmaps(
				True, variant_mask, heatmap_x_label, heatmap_y_label,
				new_gene_expression_factors,
				new_gene_translation_efficiency_values,
				figsize_x, figsize_y, plotOutDir, plot_suffix)

		# Create separate plots
		if DASHBOARD_FLAG == 0 or DASHBOARD_FLAG == 2:
			self.plot_heatmaps(
				False, variant_mask, heatmap_x_label,heatmap_y_label,
				new_gene_expression_factors, new_gene_translation_efficiency_values,
				figsize_x, figsize_y, plotOutDir, plot_suffix)

if __name__ == "__main__":
	Plot().cli()
