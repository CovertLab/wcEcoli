"""
Plot one value per index via heatmap for
new_gene_expression_and_translation_efficiency variant.

Possible Plots:
- Percent of sims that successfully reached a given generation number
- Average doubling time
- Average cell volume, mass, dry cell mass, mRNA mass, protein mass
- Average translation efficiency, weighted by cistron count
- Average mRNA count, monomer count, mRNA mass fraction, protein mass fraction,
	RNAP portion, and ribosome portion for a capacity gene to measure burden on
	overall host expression
- Average new gene mRNA count
- Average new gene mRNA mass fraction
- Average new gene mRNA counts fraction
- Average new gene NTP mass fraction
- Average new gene protein count
- Average new gene protein mass fraction
- Average new gene protein counts fraction
- Average new gene initialization rate for RNAP and Ribosomes
- Average new gene initialization probabilities for RNAP and Ribosomes
- Average number and proportion of RNAP on new genes at a given time step
- Average number and proportion of ribosomes on new gene mRNAs at a given time
	step
- Average number and proportion of RNAP making rRNAs at a given time step
- Average proportion of RNAP and ribosomes making RNAP subunits at a given time
	step
- Average proportion of RNAP and ribosomes making ribosomal proteins at a given
 	time step
- Average fraction of time new gene is overcrowded by RNAP and Ribosomes
- Average number of overcrowded genes for RNAP and Ribosomes
- Average number of total and free ribosomes
- Average number of total and free RNA polymerases
- Average ppGpp concentration
- Average rate of glucose consumption
- Average new gene monomer yields - per hour and per fg of glucose
"""

import numpy as np
from matplotlib import pyplot as plt
import math
from unum.units import fg

from models.ecoli.analysis import variantAnalysisPlot
from models.ecoli.sim.variants.new_gene_internal_shift import (
	get_new_gene_expression_factor_and_translation_efficiency,
	determine_new_gene_ids_and_indices)
from wholecell.analysis.analysis_tools import exportFigure, \
	read_stacked_columns, read_stacked_bulk_molecules, \
	stacked_cell_identification
from wholecell.analysis.plotting_tools import heatmap
from wholecell.io.tablereader import TableReader
from wholecell.utils import units

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
Standard Deviations Flag
True: Plot an additional copy of all plots with standard deviation displayed
	insted of the average
False: Plot no additional plots
"""
STD_DEV_FLAG = True

"""
Count number of sims that reach this generation (remember index 7 
corresponds to generation 8)
"""
# COUNT_INDEX = 23
COUNT_INDEX = 2 ### TODO: revert back after developing plot locally

"""
Plot data from generations [MIN_CELL_INDEX, MAX_CELL_INDEX)
Note that early generations may not be representative of dynamics 
due to how they are initialized
"""
# MIN_CELL_INDEX = 16
MIN_CELL_INDEX = 1 ### TODO: revert back after developing plot locally
MAX_CELL_INDEX = 24

"""
Specify which subset of heatmaps should be made
Completed_gens heatmap is always made, because it is used to
create the other heatmaps, and should not be included here.
The order listed here will be the order of the heatmaps in the
dashboard plot.
"""
HEATMAPS_TO_MAKE_LIST = [
		"doubling_times_heatmap",
		"cell_mass_heatmap",
		"cell_dry_mass_heatmap",
		"cell_volume_heatmap",
		"ppgpp_concentration_heatmap",
		# # "rnap_crowding_heatmap",
		# # "ribosome_crowding_heatmap",
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
		"new_gene_mRNA_counts_fraction_heatmap",
		"new_gene_monomer_counts_fraction_heatmap",
		"new_gene_rnap_counts_heatmap",
		"new_gene_rnap_portion_heatmap",
		"rrna_rnap_counts_heatmap",
		"rrna_rnap_portion_heatmap",
		"rnap_subunit_rnap_portion_heatmap",
		"rnap_subunit_ribosome_portion_heatmap",
		"ribosomal_protein_rnap_portion_heatmap",
		"ribosomal_protein_ribosome_portion_heatmap",
		"new_gene_ribosome_counts_heatmap",
		"new_gene_ribosome_portion_heatmap",
		# # "weighted_avg_translation_efficiency_heatmap",
		"new_gene_target_protein_init_prob_heatmap",
		"new_gene_actual_protein_init_prob_heatmap",
		"new_gene_target_rna_synth_prob_heatmap",
		"new_gene_actual_rna_synth_prob_heatmap",
		"capacity_gene_mRNA_counts_heatmap",
		"capacity_gene_monomer_counts_heatmap",
		"capacity_gene_rnap_portion_heatmap",
		"capacity_gene_ribosome_portion_heatmap",
		"capacity_gene_mRNA_mass_fraction_heatmap",
		"capacity_gene_monomer_mass_fraction_heatmap",
		"capacity_gene_mRNA_counts_fraction_heatmap",
		"capacity_gene_monomer_counts_fraction_heatmap",
		"free_rnap_counts_heatmap",
		"free_ribosome_counts_heatmap",
		"rnap_ribosome_counts_ratio_heatmap",
		"new_gene_yield_per_glucose",
		"new_gene_yield_per_hour",
		"glucose_consumption_rate",
		"new_gene_mRNA_NTP_fraction_heatmap",
	]

### TODO map id to common name, don't hardcode, add error checking?
capacity_gene_monomer_id = "EG10544-MONOMER[m]"
capacity_gene_common_name = "lpp"
# capacity_gene_monomer_id = "EG11036-MONOMER[c]"
# capacity_gene_common_name = "tufA"

class Plot(variantAnalysisPlot.VariantAnalysisPlot):
	"""
	General Functions
	"""
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
		self.heatmap_data[h]["mean"][
			initial_index, trl_eff_index, exp_index] = round(
			np.mean(curr_heatmap_data[cell_mask]),
			self.heatmap_details[h]['num_digits_rounding'])
		self.heatmap_data[h]["std_dev"][
			initial_index, trl_eff_index, exp_index] = round(
			np.std(curr_heatmap_data[cell_mask]),
			self.heatmap_details[h]['num_digits_rounding'])

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
			elif h == "rnap_ribosome_counts_ratio_heatmap":
				self.extract_rnap_ribosome_counts_ratio_data(
					all_cells, h, trl_eff_index, exp_index, cell_mask)
			elif h == "free_rnap_counts_heatmap":
				self.extract_rnap_counts_data(
					all_cells, h, trl_eff_index, exp_index, cell_mask,
					'inactive')
			elif h == "free_ribosome_counts_heatmap":
				self.extract_ribosome_counts_data(
					all_cells, h, trl_eff_index, exp_index, cell_mask,
					'inactive')
			elif h == "rnap_crowding_heatmap":
				self.extract_crowding_heatmap_data(
					all_cells, h, trl_eff_index, exp_index, cell_mask,
					'actual_rna_synth_prob', 'target_rna_synth_prob')
			elif h == "ribosome_crowding_heatmap":
				self.extract_crowding_heatmap_data(
					all_cells, h, trl_eff_index, exp_index, cell_mask,
					'actual_prob_translation_per_transcript',
					'target_prob_translation_per_transcript')
			elif h == "weighted_avg_translation_efficiency_heatmap":
				self.extract_trl_eff_weighted_avg_heatmap_data(
					all_cells, h, trl_eff_index, exp_index, cell_mask)
			elif h == "rrna_rnap_counts_heatmap":
				self.extract_rrna_rnap_counts_heatmap_data(
					all_cells, h, trl_eff_index, exp_index, cell_mask)
			elif h == "rrna_rnap_portion_heatmap":
				self.extract_rrna_rnap_portion_heatmap_data(
					all_cells, h, trl_eff_index, exp_index, cell_mask)
			elif h == "rnap_subunit_rnap_portion_heatmap":
				self.extract_rnap_subunits_rnap_portion_heatmap_data(
					all_cells, h, trl_eff_index, exp_index, cell_mask)
			elif h == "rnap_subunit_ribosome_portion_heatmap":
				self.extract_rnap_subunits_ribosome_portion_heatmap_data(
					all_cells, h, trl_eff_index, exp_index, cell_mask)
			elif h == "ribosomal_protein_rnap_portion_heatmap":
				self.extract_ribosomal_protein_rnap_portion_heatmap_data(
					all_cells, h, trl_eff_index, exp_index, cell_mask)
			elif h == "ribosomal_protein_ribosome_portion_heatmap":
				self.extract_ribosomal_protein_ribosome_portion_heatmap_data(
					all_cells, h, trl_eff_index, exp_index, cell_mask)
			elif h == "capacity_gene_mRNA_counts_heatmap":
				self.extract_capacity_gene_counts_heatmap_data(
					all_cells, h, trl_eff_index,
					exp_index, cell_mask, 'RNACounts', 'mRNA_counts',
					'mRNA')
			elif h == "capacity_gene_monomer_counts_heatmap":
				self.extract_capacity_gene_counts_heatmap_data(
					all_cells, h, trl_eff_index,
					exp_index, cell_mask, 'MonomerCounts', 'monomerCounts',
					'monomer')
			elif h == "capacity_gene_rnap_portion_heatmap":
				self.extract_capacity_gene_rnap_portion_heatmap_data(
					all_cells, h, trl_eff_index, exp_index, cell_mask)
			elif h == "capacity_gene_ribosome_portion_heatmap":
				self.extract_capacity_gene_ribosome_portion_heatmap_data(
					all_cells, h, trl_eff_index, exp_index, cell_mask)
			elif h == "capacity_gene_mRNA_mass_fraction_heatmap":
				self.extract_capacity_gene_mass_fraction_heatmap_data(
					all_cells, h, trl_eff_index, exp_index, cell_mask,
					'RNACounts', 'mRNA_counts', 'Mass', 'mRnaMass',
					'mRNA')
			elif h == "capacity_gene_monomer_mass_fraction_heatmap":
				self.extract_capacity_gene_mass_fraction_heatmap_data(
					all_cells, h, trl_eff_index, exp_index, cell_mask,
					'MonomerCounts', 'monomerCounts', 'Mass', 'proteinMass',
					'monomer')
			elif h == "capacity_gene_mRNA_counts_fraction_heatmap":
				self.extract_capacity_gene_counts_fraction_heatmap_data(
					all_cells, h, trl_eff_index, exp_index, cell_mask,
					'RNACounts', 'mRNA_counts', 'mRNA')
			elif h == "capacity_gene_monomer_counts_fraction_heatmap":
				self.extract_capacity_gene_counts_fraction_heatmap_data(
					all_cells, h, trl_eff_index, exp_index, cell_mask,
					'MonomerCounts', 'monomerCounts', 'monomer')
			elif h == "glucose_consumption_rate":
				self.extract_glucose_consumption_heatmap_data(
					all_cells, h, trl_eff_index, exp_index, cell_mask)
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
			elif h == "new_gene_mRNA_counts_fraction_heatmap":
				self.extract_new_gene_counts_fraction_heatmap_data(
					all_cells, h, trl_eff_index, exp_index, cell_mask,
					'RNACounts', 'mRNA_counts', 'mRNA')
			elif h == "new_gene_monomer_counts_fraction_heatmap":
				self.extract_new_gene_counts_fraction_heatmap_data(
					all_cells, h, trl_eff_index, exp_index, cell_mask,
					'MonomerCounts', 'monomerCounts', 'monomer')
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
			elif h == "new_gene_target_protein_init_prob_heatmap":
				self.extract_new_gene_init_prob_data(
					all_cells, h, trl_eff_index, exp_index, cell_mask,
					'RibosomeData', 'target_prob_translation_per_transcript',
					'monomer')
			elif h == "new_gene_actual_protein_init_prob_heatmap":
				self.extract_new_gene_init_prob_data(
					all_cells, h, trl_eff_index, exp_index, cell_mask,
					'RibosomeData', 'actual_prob_translation_per_transcript',
					'monomer')
			elif h == "new_gene_target_rna_synth_prob_heatmap":
				self.extract_new_gene_rna_synth_prob_data(
					all_cells, h, trl_eff_index, exp_index, cell_mask,
					'RnaSynthProb', 'target_rna_synth_prob',
					'RNA')
			elif h == "new_gene_actual_rna_synth_prob_heatmap":
				self.extract_new_gene_rna_synth_prob_data(
					all_cells, h, trl_eff_index, exp_index, cell_mask,
					'RnaSynthProb', 'actual_rna_synth_prob',
					'RNA')
			elif h == "new_gene_rnap_counts_heatmap":
				self.extract_new_gene_rnap_counts_heatmap_data(
					all_cells, h, trl_eff_index, exp_index, cell_mask)
			elif h == "new_gene_rnap_portion_heatmap":
				self.extract_new_gene_rnap_portion_heatmap_data(
					all_cells, h, trl_eff_index, exp_index, cell_mask)
			elif h == "new_gene_ribosome_counts_heatmap":
				self.extract_new_gene_ribosome_counts_heatmap_data(
					all_cells, h, trl_eff_index, exp_index, cell_mask)
			elif h == "new_gene_ribosome_portion_heatmap":
				self.extract_new_gene_ribosome_portion_heatmap_data(
					all_cells, h, trl_eff_index, exp_index, cell_mask)
			elif h == "new_gene_yield_per_glucose":
				self.extract_new_gene_yield_per_glucose_heatmap_data(
					all_cells, h, trl_eff_index, exp_index, cell_mask,
					self.new_gene_monomer_ids)
			elif h == "new_gene_yield_per_hour":
				self.extract_new_gene_yield_per_hour_heatmap_data(
					all_cells, h, trl_eff_index, exp_index, cell_mask,
					self.new_gene_monomer_ids)
			else:
				raise Exception(
					f"Heatmap {h} has no instructions for"
					f" data extraction.")

	def get_mRNA_ids_from_monomer_ids(self, target_monomer_ids):
		"""
		Map monomer ids back to the mRNA ids that they were translated from.

		Args:
			target_monomer_ids: ids of the monomers to map to mRNA ids

		Returns: set of mRNA ids
		"""
		# Map protein ids to cistron ids
		monomer_ids = self.sim_data.process.translation.monomer_data['id']
		cistron_ids = self.sim_data.process.translation.monomer_data[
			'cistron_id']
		monomer_to_cistron_id_dict = {
			monomer_id: cistron_ids[i] for i, monomer_id in
			enumerate(monomer_ids)}
		target_cistron_ids = [
			monomer_to_cistron_id_dict.get(RNAP_monomer_id) for
			RNAP_monomer_id in target_monomer_ids]
		# Map cistron ids to RNA indexes
		target_RNA_indexes = [
			self.sim_data.process.transcription.cistron_id_to_rna_indexes(
				RNAP_cistron_id) for RNAP_cistron_id in
			target_cistron_ids]
		# Map RNA indexes to RNA ids
		RNA_ids = self.sim_data.process.transcription.rna_data['id']
		target_RNA_ids = set()
		for i in range(len(target_RNA_indexes)):
			for index in target_RNA_indexes[i]:
				target_RNA_ids.add(RNA_ids[index])
		return target_RNA_ids

	def get_mRNA_indexes_from_monomer_ids(
			self, all_cells, target_monomer_ids, index_type):
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
			target_RNA_ids = self.get_mRNA_ids_from_monomer_ids(target_monomer_ids)
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
				f"Index type {index_type} has no instructions for data"
				f" extraction.")

		return output_indexes

	def extract_trl_eff_weighted_avg_heatmap_data(
			self, all_cells, h, trl_eff_index, exp_index, cell_mask):
		"""
		Special function to handle extraction and saving of average translation
		efficiency, weighted by full cistron counts.

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
		# Get normalized translation efficiency for all mRNAs
		trl_effs = self.sim_data.process.translation.translation_efficiencies_by_monomer
		# Get avg counts for all mRNAs
		mRNA_cistron_counts = read_stacked_columns(
			all_cells, 'RNACounts', 'full_mRNA_cistron_counts',
			fun=lambda x: np.mean(x, axis=0))[cell_mask, :]
		total_mRNA_cistron_count = np.expand_dims(
			np.sum(mRNA_cistron_counts,axis = 1), axis = 1)

		sim_dir = all_cells[0]
		simOutDir = os.path.join(sim_dir, 'simOut')
		mRNA_counts_reader = TableReader(os.path.join(simOutDir, 'RNACounts'))
		mRNA_cistron_idx_dict = {
			rna: i for i, rna in
			enumerate(mRNA_counts_reader.readAttribute('mRNA_cistron_ids'))}
		trl_eff_ids = self.sim_data.process.translation.monomer_data['cistron_id']
		trl_eff_id_mapping = np.array([
			mRNA_cistron_idx_dict[id] for id in trl_eff_ids])

		# Compute average translation efficiency, weighted by mRNA counts
		weighted_avg_trl_eff = np.array([
			np.sum(mRNA_cistron_counts / total_mRNA_cistron_count
			* trl_effs[np.argsort(trl_eff_id_mapping)], axis = 1)])

		all_true_mask = np.ones_like(weighted_avg_trl_eff, dtype=bool)
		self.save_heatmap_data(
			h, 0, trl_eff_index, exp_index, weighted_avg_trl_eff, all_true_mask)


	"""
	Shared RNA Polymerase and Ribosome Functions
	"""
	# Crowding
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
		# Get counts of overcrowded genes for each simulation
		num_overcrowded_indexes = np.sum((
				avg_actual_prob < avg_target_prob), axis = 1)
		self.save_heatmap_data(
			h, 0, trl_eff_index, exp_index, num_overcrowded_indexes, cell_mask)

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
		new_gene_indexes = self.get_new_gene_indexes(
			all_cells, new_gene_index_type)
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

	def extract_rnap_ribosome_counts_ratio_data(
			self, all_cells, h, trl_eff_index, exp_index, cell_mask,
			category='total'):
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
			category: 'active', 'inactive', or 'total'
		"""
		if category == 'inactive':
			avg_rnap_counts = self.get_avg_inactive_rnap_counts(all_cells)
			avg_ribosome_counts = self.get_avg_inactive_ribosome_counts(all_cells)
		elif category == 'active':
			avg_rnap_counts = self.get_avg_active_rnap_counts(all_cells)
			avg_ribosome_counts = self.get_avg_active_ribosome_counts(all_cells)
		elif category == 'total':
			avg_rnap_counts = self.get_avg_total_rnap_counts(all_cells)
			avg_ribosome_counts = self.get_avg_total_ribosome_counts(all_cells)
		else:
			raise Exception("The only supported categories are 'active',"
							" 'inactive', and 'total.'")
		avg_counts_ratio = avg_rnap_counts / avg_ribosome_counts

		self.save_heatmap_data(
			h, 0, trl_eff_index, exp_index, avg_counts_ratio, cell_mask)


	"""
	RNA Polymerase Functions
	"""
	# RNA Polymerase Counts
	def get_avg_inactive_rnap_counts(self, all_cells):
		"""
		Retrieve inactive RNAP counts.

		Args:
			all_cells: paths to all cells to read data from (directories should
				contain a simOut/ subdirectory), typically the return from
				AnalysisPaths.get_cells()

		Returns:
			Average counts of RNAPs.
		"""
		rnap_id = [self.sim_data.molecule_ids.full_RNAP]
		(rnapCountsBulk,) = read_stacked_bulk_molecules(all_cells, (rnap_id,))
		cell_id_vector = stacked_cell_identification(all_cells, 'Main', 'time')
		cell_ids, idx, cell_total_timesteps = np.unique(
			cell_id_vector, return_inverse=True, return_counts=True)
		sum_rnap_counts = np.bincount(idx, weights=rnapCountsBulk)
		avg_inactive_rnap_counts = (sum_rnap_counts / cell_total_timesteps)
		return avg_inactive_rnap_counts

	def get_active_rnap_counts(self, all_cells):
		"""
		Retrieve active RNAP counts for all simulations.

		Args:
			all_cells: paths to all cells to read data from (directories should
				contain a simOut/ subdirectory), typically the return from
				AnalysisPaths.get_cells()

		Returns:
			Counts of RNAPs for all simulations.
		"""
		# Determine active RNAP index
		sim_dir = all_cells[0]
		simOutDir = os.path.join(sim_dir, 'simOut')
		uniqueMoleculeCounts = TableReader(
			os.path.join(simOutDir, "UniqueMoleculeCounts"))
		active_rnap_index = uniqueMoleculeCounts.readAttribute(
			"uniqueMoleculeIds").index('active_RNAP')
		active_rnap_counts = read_stacked_columns(
			all_cells, 'UniqueMoleculeCounts',
			'uniqueMoleculeCounts')[:, active_rnap_index]
		return active_rnap_counts

	def get_avg_active_rnap_counts(self, all_cells):
		"""
		Retrieve average active RNAP counts for each simulation.

		Args:
			all_cells: paths to all cells to read data from (directories should
				contain a simOut/ subdirectory), typically the return from
				AnalysisPaths.get_cells()

		Returns:
			Average counts of RNAPs for each simulation.
		"""
		# Determine active RNAP index
		sim_dir = all_cells[0]
		simOutDir = os.path.join(sim_dir, 'simOut')
		uniqueMoleculeCounts = TableReader(
			os.path.join(simOutDir, "UniqueMoleculeCounts"))
		active_rnap_index = uniqueMoleculeCounts.readAttribute(
			"uniqueMoleculeIds").index('active_RNAP')
		avg_active_rnap_counts = read_stacked_columns(
			all_cells, 'UniqueMoleculeCounts',
			'uniqueMoleculeCounts',
			fun=lambda x: np.mean(x[:, active_rnap_index], axis=0))
		return avg_active_rnap_counts

	def get_avg_total_rnap_counts(self, all_cells):
		"""
		Retrieve total (active + inactive) RNAP counts.

		Args:
			all_cells: paths to all cells to read data from (directories should
				contain a simOut/ subdirectory), typically the return from
				AnalysisPaths.get_cells()

		Returns:
			Average counts of RNAPs.
		"""
		return (self.get_avg_inactive_rnap_counts(all_cells) +
			self.get_avg_active_rnap_counts(all_cells))

	def extract_rnap_counts_data(
			self, all_cells, h, trl_eff_index, exp_index, cell_mask,
			rnap_category='total'):
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
			rnap_category: 'active', 'inactive', or 'total'
		"""
		if rnap_category == 'inactive':
			avg_rnap_counts = self.get_avg_inactive_rnap_counts(all_cells)
		elif rnap_category == 'active':
			avg_rnap_counts = self.get_avg_active_rnap_counts(all_cells)
		elif rnap_category == 'total':
			avg_rnap_counts = self.get_avg_total_rnap_counts(all_cells)
		else:
			raise Exception("The only supported RNAP categories are 'active',"
							" 'inactive', and 'total.'")
		self.save_heatmap_data(
			h, 0, trl_eff_index, exp_index, avg_rnap_counts, cell_mask)

	def extract_new_gene_rnap_counts_heatmap_data(
			self, all_cells, h, trl_eff_index, exp_index, cell_mask):
		"""
		Special function to handle extraction and saving average counts of RNAP
		that are on new genes at a time heatmap data.

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
		avg_new_gene_rnap_counts = read_stacked_columns(
			all_cells, "RNACounts", "partial_mRNA_counts",
			fun=lambda x: np.mean(x[:, new_gene_mRNA_indexes], axis=0))

		for i in range(len(self.new_gene_mRNA_ids)):
			self.save_heatmap_data(
				h, i, trl_eff_index, exp_index,
				avg_new_gene_rnap_counts[:, i], cell_mask)

	def extract_rrna_rnap_counts_heatmap_data(
			self, all_cells, h, trl_eff_index, exp_index, cell_mask):
		"""
		Special function to handle extraction and saving average counts of RNAP
		that are making rRNAs at a time heatmap data.

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
		avg_rrna_rnap_counts = np.sum(read_stacked_columns(
			all_cells, "RNACounts", "partial_rRNA_counts",
			fun=lambda x: np.mean(x, axis=0)), axis = 1)

		self.save_heatmap_data(
			h, 0, trl_eff_index, exp_index, avg_rrna_rnap_counts, cell_mask)

	# RNA Polymerase Initialization Rate
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
		new_gene_cistron_indexes = self.get_new_gene_indexes(
			all_cells, 'cistron')
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

	# RNA Polymerase Portion
	def extract_new_gene_rnap_portion_heatmap_data(
			self, all_cells, h, trl_eff_index, exp_index, cell_mask):
		"""
		Special function to handle extraction and saving average portion of RNAP
		that are on new genes at a time heatmap data.

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
		new_gene_rnap_counts = read_stacked_columns(
			all_cells, "RNACounts", "partial_mRNA_counts",
			fun=lambda x: x[:, new_gene_mRNA_indexes])
		active_rnap_counts = self.get_active_rnap_counts(all_cells)
		active_rnap_counts = np.expand_dims(active_rnap_counts, axis = 1)
		avg_new_gene_rnap_portion = np.mean(
			new_gene_rnap_counts / active_rnap_counts, axis = 0, keepdims = True)

		# TODO: fix for standard deviation?
		for i in range(len(self.new_gene_mRNA_ids)):
			self.save_heatmap_data(
				h, i, trl_eff_index, exp_index,
				avg_new_gene_rnap_portion[:, i], [True])

	def extract_rrna_rnap_portion_heatmap_data(
			self, all_cells, h, trl_eff_index, exp_index, cell_mask):
		"""
		Special function to handle extraction and saving average portion of RNAP
		that are making rRNAs at a time heatmap data.

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
		rrna_rnap_counts = np.sum(read_stacked_columns(
			all_cells, "RNACounts", "partial_rRNA_counts"), axis=1)
		active_rnap_counts = self.get_active_rnap_counts(all_cells)
		avg_rrna_rnap_portion = np.mean(
			rrna_rnap_counts / active_rnap_counts, axis=0, keepdims=True)

		# TODO: fix standard deviation?
		self.save_heatmap_data(
			h, 0, trl_eff_index, exp_index, avg_rrna_rnap_portion, [True])

	def extract_rnap_subunits_rnap_portion_heatmap_data(
			self, all_cells, h, trl_eff_index, exp_index, cell_mask):
		"""
		Special function to handle extraction and saving average portion of RNAP
		that are making RNAP subunits at a time heatmap data.

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
		RNAP_subunit_monomer_ids = self.sim_data.molecule_groups.RNAP_subunits
		rnap_subunit_mRNA_indexes = self.get_mRNA_indexes_from_monomer_ids(
			all_cells, RNAP_subunit_monomer_ids, "mRNA")
		rnap_subunit_rnap_counts = np.sum(read_stacked_columns(
			all_cells, "RNACounts", "partial_mRNA_counts",
			fun=lambda x: x[:, rnap_subunit_mRNA_indexes]), axis=1)
		active_rnap_counts = self.get_active_rnap_counts(all_cells)
		avg_rnap_subunit_rnap_portion = np.mean(
			rnap_subunit_rnap_counts / active_rnap_counts, axis=0,
			keepdims=True)
		# TODO: fix standard deviation
		self.save_heatmap_data(
			h, 0, trl_eff_index, exp_index, avg_rnap_subunit_rnap_portion,
			[True])

	def extract_ribosomal_protein_rnap_portion_heatmap_data(
			self, all_cells, h, trl_eff_index, exp_index, cell_mask):
		"""
		Special function to handle extraction and saving average portion of RNAP
		that are making ribosomal proteins at a time heatmap data.

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
		ribosomal_monomer_ids = self.sim_data.molecule_groups.ribosomal_proteins
		ribosomal_mRNA_indexes = self.get_mRNA_indexes_from_monomer_ids(
			all_cells, ribosomal_monomer_ids, "mRNA")
		ribosomal_rnap_counts = np.sum(read_stacked_columns(
			all_cells, "RNACounts", "partial_mRNA_counts",
			fun=lambda x: x[:, ribosomal_mRNA_indexes]), axis=1)
		active_rnap_counts = self.get_active_rnap_counts(all_cells)
		avg_ribosomal_rnap_portion = np.mean(
			ribosomal_rnap_counts / active_rnap_counts, axis=0, keepdims=True)

		# TODO: fix standard deviation
		self.save_heatmap_data(
			h, 0, trl_eff_index, exp_index, avg_ribosomal_rnap_portion,
			[True])

	def extract_capacity_gene_rnap_portion_heatmap_data(
			self, all_cells, h, trl_eff_index, exp_index, cell_mask):
		"""
		Special function to handle extraction and saving average portion of RNAP
		that are making capacity gene mRNAs at a time heatmap data.

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
		gene_mRNA_indexes = self.get_mRNA_indexes_from_monomer_ids(
			all_cells, [capacity_gene_monomer_id], "mRNA")
		gene_rnap_counts = np.sum(read_stacked_columns(
			all_cells, "RNACounts", "partial_mRNA_counts",
			fun=lambda x: x[:, gene_mRNA_indexes]), axis=1)
		active_rnap_counts = self.get_active_rnap_counts(all_cells)
		avg_gene_rnap_portion = np.mean(
			gene_rnap_counts / active_rnap_counts, axis=0, keepdims=True)

		# TODO: fix standard deviation?
		self.save_heatmap_data(
			h, 0, trl_eff_index, exp_index, avg_gene_rnap_portion,
			[True])


	"""
	Ribosome Functions
	"""
	# Ribosome Counts
	def get_avg_inactive_ribosome_counts(self, all_cells):
		"""
		Retrieve inactive ribosome counts.

		Args:
			all_cells: paths to all cells to read data from (directories should
				contain a simOut/ subdirectory), typically the return from
				AnalysisPaths.get_cells()

		Returns:
			Average counts of ribosomes.
		"""
		complex_id_30s = [self.sim_data.molecule_ids.s30_full_complex]
		complex_id_50s = [self.sim_data.molecule_ids.s50_full_complex]

		(complex_counts_30s, complex_counts_50s) = read_stacked_bulk_molecules(
			all_cells, (complex_id_30s, complex_id_50s))
		cell_id_vector = stacked_cell_identification(all_cells, 'Main', 'time')
		cell_ids, idx, cell_total_timesteps = np.unique(
			cell_id_vector, return_inverse=True, return_counts=True)

		inactive_ribosome_counts = np.minimum(
			complex_counts_30s, complex_counts_50s)

		sum_inactive_ribosome_counts = np.bincount(
			idx, weights=inactive_ribosome_counts)
		avg_inactive_ribosome_counts = (sum_inactive_ribosome_counts / cell_total_timesteps)

		return avg_inactive_ribosome_counts

	def get_avg_active_ribosome_counts(self, all_cells):
		"""
		Retrieve active ribosome counts.

		Args:
			all_cells: paths to all cells to read data from (directories should
				contain a simOut/ subdirectory), typically the return from
				AnalysisPaths.get_cells()

		Returns:
			Average counts of ribosomes.
		"""
		# Determine ribosome index
		sim_dir = all_cells[0]
		simOutDir = os.path.join(sim_dir, 'simOut')
		uniqueMoleculeCounts = TableReader(
			os.path.join(simOutDir, "UniqueMoleculeCounts"))
		ribosome_index = uniqueMoleculeCounts.readAttribute(
			"uniqueMoleculeIds").index('active_ribosome')

		avg_ribosome_counts = read_stacked_columns(
			all_cells, 'UniqueMoleculeCounts',
			'uniqueMoleculeCounts',
			fun=lambda x: np.mean(x[:, ribosome_index], axis=0))

		return avg_ribosome_counts

	def get_avg_total_ribosome_counts(self, all_cells):
		"""
		Retrieve ribosome counts.

		Args:
			all_cells: paths to all cells to read data from (directories should
				contain a simOut/ subdirectory), typically the return from
				AnalysisPaths.get_cells()

		Returns:
			Average counts of ribosomes.
		"""

		return (self.get_avg_inactive_ribosome_counts(all_cells) +
			self.get_avg_active_ribosome_counts(all_cells))

	def extract_ribosome_counts_data(
			self, all_cells, h, trl_eff_index, exp_index, cell_mask,
			ribosome_category='total'):
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
			ribosome_category: 'active', 'inactive', or 'total'
		"""
		if ribosome_category == 'inactive':
			avg_ribosome_counts = self.get_avg_inactive_ribosome_counts(all_cells)
		elif ribosome_category == 'active':
			avg_ribosome_counts = self.get_avg_active_ribosome_counts(all_cells)
		elif ribosome_category == 'total':
			avg_ribosome_counts = self.get_avg_total_ribosome_counts(all_cells)
		else:
			raise Exception("The only supported ribosome categories are 'active',"
							" 'inactive', and 'total.'")
		self.save_heatmap_data(
			h, 0, trl_eff_index, exp_index, avg_ribosome_counts, cell_mask)

	def extract_new_gene_ribosome_counts_heatmap_data(
			self, all_cells, h, trl_eff_index, exp_index, cell_mask):
		"""
		Special function to handle extraction and saving average counts of
		ribosomes that are on new genes at a time heatmap data.

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
		new_gene_monomer_indexes = self.get_new_gene_indexes(all_cells, 'monomer')
		avg_new_gene_ribosome_counts = read_stacked_columns(
			all_cells, "RibosomeData", "n_ribosomes_per_transcript",
			fun=lambda x: np.mean(x[:, new_gene_monomer_indexes], axis=0))

		for i in range(len(self.new_gene_monomer_ids)):
			self.save_heatmap_data(
				h, i, trl_eff_index, exp_index,
				avg_new_gene_ribosome_counts[:, i], cell_mask)

	# Ribosome Initialization Rate
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
		avg_new_gene_mRNA_counts = self.get_avg_gene_counts(
			all_cells, 'RNACounts', 'mRNA_counts', new_gene_mRNA_indexes)
		avg_new_gene_ribosome_init_rates = (read_stacked_columns(
			all_cells, 'RibosomeData', 'ribosome_init_event_per_monomer',
			fun=lambda x: np.mean(x[:, new_gene_monomer_indexes],
			axis=0))) / avg_new_gene_mRNA_counts
		for i in range(len(self.new_gene_mRNA_ids)):
			self.save_heatmap_data(
				h, i, trl_eff_index, exp_index,
				avg_new_gene_ribosome_init_rates[:, i], cell_mask)

	# Ribosome Portion
	def extract_new_gene_ribosome_portion_heatmap_data(
			self, all_cells, h, trl_eff_index, exp_index, cell_mask):
		"""
		Special function to handle extraction and saving average portion of
		ribosomes that are on new genes at a time heatmap data.

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
		new_gene_monomer_indexes = self.get_new_gene_indexes(all_cells, 'monomer')
		avg_new_gene_ribosome_counts = read_stacked_columns(
			all_cells, "RibosomeData", "n_ribosomes_per_transcript",
			fun=lambda x: np.mean(x[:, new_gene_monomer_indexes], axis=0))
		avg_ribosome_counts = self.get_avg_active_ribosome_counts(all_cells)
		avg_new_gene_ribosome_portion = avg_new_gene_ribosome_counts/avg_ribosome_counts

		for i in range(len(self.new_gene_monomer_ids)):
			self.save_heatmap_data(
				h, i, trl_eff_index, exp_index,
				avg_new_gene_ribosome_portion[:, i], cell_mask)

	def extract_rnap_subunits_ribosome_portion_heatmap_data(
			self, all_cells, h, trl_eff_index, exp_index, cell_mask):
		"""
		Special function to handle extraction and saving average portion of ribosomes
		that are making RNAP subunits at a time heatmap data.

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
		RNAP_subunit_monomer_ids = self.sim_data.molecule_groups.RNAP_subunits
		rnap_subunit_monomer_indexes = self.get_mRNA_indexes_from_monomer_ids(
			all_cells, RNAP_subunit_monomer_ids, "monomer")
		avg_rnap_subunit_ribosome_counts = np.sum(read_stacked_columns(
			all_cells, "RibosomeData", "n_ribosomes_per_transcript",
			fun=lambda x: np.mean(x[:, rnap_subunit_monomer_indexes], axis=0)),
			axis = 1)
		avg_ribosome_counts = self.get_avg_active_ribosome_counts(all_cells)
		avg_rnap_subunit_ribosome_portion = avg_rnap_subunit_ribosome_counts/avg_ribosome_counts
		self.save_heatmap_data(
			h, 0, trl_eff_index, exp_index, avg_rnap_subunit_ribosome_portion,
			cell_mask)

	def extract_ribosomal_protein_ribosome_portion_heatmap_data(
			self, all_cells, h, trl_eff_index, exp_index, cell_mask):
		"""
		Special function to handle extraction and saving average portion of ribosomes
		that are making ribosomal proteins at a time heatmap data.

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
		ribosomal_monomer_ids = self.sim_data.molecule_groups.ribosomal_proteins
		ribosomal_monomer_indexes = self.get_mRNA_indexes_from_monomer_ids(
			all_cells, ribosomal_monomer_ids, "monomer")
		avg_ribosomal_ribosome_counts = np.sum(read_stacked_columns(
			all_cells, "RibosomeData", "n_ribosomes_per_transcript",
			fun=lambda x: np.mean(x[:, ribosomal_monomer_indexes], axis=0)),
			axis=1)
		avg_ribosome_counts = self.get_avg_active_ribosome_counts(all_cells)
		avg_ribosomal_ribosome_portion = avg_ribosomal_ribosome_counts / avg_ribosome_counts
		self.save_heatmap_data(
			h, 0, trl_eff_index, exp_index, avg_ribosomal_ribosome_portion,
			cell_mask)

	def extract_capacity_gene_ribosome_portion_heatmap_data(
			self, all_cells, h, trl_eff_index, exp_index, cell_mask):
		"""
		Special function to handle extraction and saving average portion of ribosomes
		that are making capacity gene proteins at a time heatmap data.

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
		gene_monomer_indexes = self.get_mRNA_indexes_from_monomer_ids(
			all_cells, [capacity_gene_monomer_id], "monomer")
		avg_gene_ribosome_counts = np.sum(read_stacked_columns(
			all_cells, "RibosomeData", "n_ribosomes_per_transcript",
			fun=lambda x: np.mean(x[:, gene_monomer_indexes], axis=0)),
			axis = 1)
		avg_ribosome_counts = self.get_avg_active_ribosome_counts(all_cells)
		avg_gene_ribosome_portion = avg_gene_ribosome_counts/avg_ribosome_counts
		self.save_heatmap_data(
			h, 0, trl_eff_index, exp_index, avg_gene_ribosome_portion,
			cell_mask)


	"""
	New Gene Functions
	"""
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
				cis: i for i, cis in
				enumerate(rnap_reader.readAttribute('cistron_ids'))}
			new_gene_indexes = [
				cistron_idx_dict.get(mRNA_id) for mRNA_id in
				self.new_gene_mRNA_ids]
		elif index_type == 'RNA':
			# Extract RNA indexes for each new gene
			rnap_reader = TableReader(os.path.join(simOutDir, 'RnaSynthProb'))
			RNA_idx_dict = {
				rna[:-3]: i for i, rna in
				enumerate(rnap_reader.readAttribute('rnaIds'))}
			new_gene_indexes = [
				RNA_idx_dict.get(mRNA_id) for mRNA_id in self.new_gene_mRNA_ids]
		elif index_type == 'mRNA':
			# Extract mRNA indexes for each new gene
			mRNA_counts_reader = TableReader(os.path.join(simOutDir, 'RNACounts'))
			mRNA_idx_dict = {
				rna[:-3]: i for i, rna in
				enumerate(mRNA_counts_reader.readAttribute('mRNA_ids'))}
			new_gene_indexes = [
				mRNA_idx_dict.get(mRNA_id) for mRNA_id in self.new_gene_mRNA_ids]
		elif index_type == 'monomer':
			# Extract protein indexes for each new gene
			monomer_counts_reader = TableReader(
				os.path.join(simOutDir, "MonomerCounts"))
			monomer_idx_dict = {
				monomer: i for i, monomer in enumerate(
				monomer_counts_reader.readAttribute('monomerIds'))}
			new_gene_indexes = [
				monomer_idx_dict.get(monomer_id) for monomer_id in
				self.new_gene_monomer_ids]
		else:
			raise Exception(
				"Index type " + index_type +
				" has no instructions for data extraction.")

		return new_gene_indexes

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
		avg_new_gene_counts = self.get_avg_gene_counts(all_cells, data_table,
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
		new_gene_indexes = self.get_new_gene_indexes(
			all_cells, new_gene_index_type)
		# Get mass for each new gene
		new_gene_masses = [1 for id in new_gene_ids]
		for i in range(len(new_gene_ids)):
			new_gene_masses[i] = (
				self.sim_data.getter.get_mass(
				new_gene_ids[i])/self.sim_data.constants.n_avogadro).asNumber(fg)
		# Get counts for each new gene
		avg_new_gene_counts = self.get_avg_gene_counts(
			all_cells, counts_data_table, counts_data_column, new_gene_indexes)
		# Get average mass for all genes
		avg_mass = read_stacked_columns(
			all_cells, mass_data_table, mass_data_column,
			fun=lambda x: np.mean(x))
		# Determine mass fractions for each new gene
		for i in range(len(self.new_gene_mRNA_ids)):
			new_gene_mass_fraction = (avg_new_gene_counts[:, i] *
				new_gene_masses[i]) / avg_mass
			self.save_heatmap_data(
				h, i, trl_eff_index, exp_index, new_gene_mass_fraction,
				cell_mask)

	def extract_new_gene_counts_fraction_heatmap_data(
			self, all_cells, h, trl_eff_index, exp_index, cell_mask,
			counts_data_table, counts_data_column, new_gene_index_type):
		"""
		Special function to handle extraction and saving of new gene mRNA and
		protein counts fraction heatmap data.

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
			new_gene_index_type: Index type to use for the data table
		"""
		new_gene_indexes = self.get_new_gene_indexes(
			all_cells, new_gene_index_type)
		# Get counts for each new gene
		avg_new_gene_counts = self.get_avg_gene_counts(
			all_cells, counts_data_table, counts_data_column, new_gene_indexes)
		# Get total avg counts for all genes
		total_counts = np.sum(read_stacked_columns(
			all_cells, counts_data_table, counts_data_column,
			fun=lambda x: np.mean(x, axis=0)), axis = 1)

		# Determine count fractions for each new gene
		for i in range(len(self.new_gene_mRNA_ids)):
			new_gene_counts_fraction = avg_new_gene_counts[:, i] / total_counts
			self.save_heatmap_data(
				h, i, trl_eff_index, exp_index, new_gene_counts_fraction,
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
		avg_new_gene_mRNA_counts = self.get_avg_gene_counts(
			all_cells, 'RNACounts', 'mRNA_counts', new_gene_mRNA_indexes)
		# Compute new gene NTP fractions
		all_mRNA_ntp_totals = {}
		for ntp_index in range(len(self.ntp_ids)):
			ntp_id = self.ntp_ids[ntp_index]
			all_mRNA_ntp_totals[ntp_id] = \
				(avg_mRNA_counts @ all_mRNA_counts_ACGU[:, ntp_index])
			for i in range(len(self.new_gene_mRNA_ids)):
				self.heatmap_data[h]["mean"][ntp_id][
					i, trl_eff_index, exp_index] = round(
					np.mean((avg_new_gene_mRNA_counts[:, i][cell_mask] *
					new_gene_mRNA_ntp_counts[i][ntp_id]) / all_mRNA_ntp_totals[
					ntp_id][cell_mask]),
					self.heatmap_details[h]['num_digits_rounding'])
				self.heatmap_data[h]["std_dev"][ntp_id][
					i, trl_eff_index, exp_index] = round(
					np.std((avg_new_gene_mRNA_counts[:, i][cell_mask] *
							 new_gene_mRNA_ntp_counts[i][ntp_id]) /
							all_mRNA_ntp_totals[
								ntp_id][cell_mask]),
					self.heatmap_details[h]['num_digits_rounding'])

	def extract_new_gene_init_prob_data(
			self, all_cells, h, trl_eff_index, exp_index, cell_mask,
			data_table, data_column, new_gene_index_type):
		"""
		Special function to handle extraction and saving of target and actual
		initiation probabilities for new genes.

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
		new_gene_indexes = self.get_new_gene_indexes(
			all_cells, new_gene_index_type)
		# Average init probability for each new gene
		new_gene_init_probs = read_stacked_columns(
			all_cells, data_table, data_column,
			fun=lambda x: np.mean(x[:, new_gene_indexes], axis=0))
		for i in range(len(new_gene_indexes)):
			self.save_heatmap_data(
				h, i, trl_eff_index, exp_index,
				new_gene_init_probs[:,i], cell_mask)

	def extract_new_gene_rna_synth_prob_data(
			self, all_cells, h, trl_eff_index, exp_index, cell_mask,
			data_table, data_column, new_gene_index_type):
		"""
		Special function to handle extraction and saving of target and actual
		RNA synthesis probabilities for new genes.

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
		new_gene_indexes = self.get_new_gene_indexes(
			all_cells, new_gene_index_type)
		# Average init probability for each new gene
		new_gene_rna_synth_probs = read_stacked_columns(
			all_cells, data_table, data_column,
			fun=lambda x: np.mean(x[:, new_gene_indexes], axis=0))
		for i in range(len(new_gene_indexes)):
			self.save_heatmap_data(
				h, i, trl_eff_index, exp_index,
				new_gene_rna_synth_probs[:,i], cell_mask)

	def get_gen_new_gene_monomer_counts_diff(self, all_cells, new_gene_indexes):
		"""
		Retrieves the difference between the number of proteins at the start
		and end of each generation.

		Args:
			all_cells: paths to all cells to read data from (directories should
				contain a simOut/ subdirectory), typically the return from
				AnalysisPaths.get_cells()
			new_gene_indexes: Global indexes of the new gene proteins

		Returns:
			Number of new gene proteins at final time step - number of new gene
				proteins at initial time step for each cell in all_cells
		"""
		return (read_stacked_columns(
				all_cells, "MonomerCounts", "monomerCounts", fun=lambda
				x: (x[-1, new_gene_indexes] - x[0, new_gene_indexes])))

	def get_avg_glucose_consumption_rate(self, all_cells):
		"""
		Computes average glucose consumption rate (fg/h) for each cell.

		Args:
			all_cells: paths to all cells to read data from (directories should
				contain a simOut/ subdirectory), typically the return from
				AnalysisPaths.get_cells()

		Returns:
			Average glucose consumption rate (fg/h) for each cell
		"""
		GLUCOSE_ID = "GLC[p]"
		FLUX_UNITS = units.mmol / units.g / units.h
		MASS_UNITS = units.fg

		# Determine glucose index in exchange fluxes
		sim_dir = all_cells[0]
		simOutDir = os.path.join(sim_dir, 'simOut')
		fba_results = TableReader(os.path.join(simOutDir, "FBAResults"))
		external_molecule_ids = np.array(
			fba_results.readAttribute("externalMoleculeIDs"))
		fba_results.close()
		if GLUCOSE_ID not in external_molecule_ids:
			print("This plot only runs when glucose is the carbon source.")
			return
		glucose_idx = np.where(external_molecule_ids == GLUCOSE_ID)[0][0]

		glucose_flux = FLUX_UNITS * read_stacked_columns(
			all_cells, "FBAResults", "externalExchangeFluxes",
			ignore_exception=True, remove_first=True,
			fun=lambda x: np.mean(x[:, glucose_idx]))
		glucose_mw = self.sim_data.getter.get_mass(GLUCOSE_ID)
		cell_dry_mass = MASS_UNITS * read_stacked_columns(
			all_cells, "Mass", "dryMass", ignore_exception=True,
			remove_first=True, fun=lambda x: np.mean(x))

		glucose_mass_flux = glucose_flux * glucose_mw * cell_dry_mass

		return -(glucose_mass_flux.asNumber())

	# TODO: make this more flexible for other carbon sources
	def extract_new_gene_yield_per_glucose_heatmap_data(
			self, all_cells, h, trl_eff_index, exp_index, cell_mask,
			new_gene_ids):
		"""
		Special function to handle extraction and saving of fg new gene
		monomer produced per fg of glucose uptake.

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
			new_gene_ids: Ids of new gene monomers in sim_data
		"""
		new_gene_indexes = self.get_new_gene_indexes(all_cells, "monomer")

		# Get mass for each new gene
		new_gene_masses = [1 for id in new_gene_ids]
		for i in range(len(new_gene_ids)):
			new_gene_masses[i] = (
				self.sim_data.getter.get_mass(
				new_gene_ids[i])/self.sim_data.constants.n_avogadro).asNumber(fg)

		# Determine how much new gene mass was acquired each gen
		new_gene_monomer_diffs = self.get_gen_new_gene_monomer_counts_diff(
			all_cells, new_gene_indexes)
		new_gene_monomer_mass_diffs = new_gene_monomer_diffs * new_gene_masses

		# Doubling time of each gen
		dt_h = "doubling_times_heatmap"
		dt = read_stacked_columns(
			all_cells, self.heatmap_details[dt_h]["data_table"],
			self.heatmap_details[dt_h]["data_column"],
			fun=self.heatmap_details[dt_h]["function_to_apply"])

		# Average glucose consumption rate fg/h
		avg_glucose_consumption_rate = self.get_avg_glucose_consumption_rate(
			all_cells)

		# Avg new gene monomer mass yield
		# Numerator: fg new gene monomer mass gained this gen
		# Denominator: (fg/hr glucose * dt in hours) approximates amount of
		# glucose consumed this gen
		avg_new_gene_yield = (
			new_gene_monomer_mass_diffs / (avg_glucose_consumption_rate * dt/60.))

		for i in range(len(new_gene_indexes)):
			self.save_heatmap_data(
				h, i, trl_eff_index, exp_index, avg_new_gene_yield[:, i],
				cell_mask)

	def extract_glucose_consumption_heatmap_data(
			self, all_cells, h, trl_eff_index, exp_index, cell_mask):
		"""
		Special function to handle extraction and saving of fg/hr of glucose uptake.

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
		# Average glucose consumption rate fg/h
		avg_glucose_consumption_rate = self.get_avg_glucose_consumption_rate(all_cells)

		self.save_heatmap_data(
			h, 0, trl_eff_index, exp_index, avg_glucose_consumption_rate,
			cell_mask)

	def extract_new_gene_yield_per_hour_heatmap_data(
			self, all_cells, h, trl_eff_index, exp_index, cell_mask,
			new_gene_ids):
		"""
		Special function to handle extraction and saving of fg new gene
		monomer produced per hour.

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
			new_gene_ids: Ids of new gene monomers in sim_data
		"""
		new_gene_indexes = self.get_new_gene_indexes(all_cells, "monomer")

		# Get mass for each new gene
		new_gene_masses = [1 for id in new_gene_ids]
		for i in range(len(new_gene_ids)):
			new_gene_masses[i] = (
				self.sim_data.getter.get_mass(
				new_gene_ids[i])/self.sim_data.constants.n_avogadro).asNumber(fg)

		# Determine how much new gene mass was acquired each gen
		new_gene_monomer_diffs = self.get_gen_new_gene_monomer_counts_diff(
			all_cells, new_gene_indexes)
		new_gene_monomer_mass_diffs = new_gene_monomer_diffs * new_gene_masses

		# Doubling time of each gen
		dt_h = "doubling_times_heatmap"
		dt = read_stacked_columns(
			all_cells, self.heatmap_details[dt_h]["data_table"],
			self.heatmap_details[dt_h]["data_column"],
			fun=self.heatmap_details[dt_h]["function_to_apply"])

		# Avg new gene monomer yield rate (fg new gene monomer mass / hr)
		# Numerator: (fg) new gene monomer mass gained this gen
		# Denominator: (dt in hours) length of gen
		avg_new_gene_yield_rate = (
				new_gene_monomer_mass_diffs / (dt/60.))

		for i in range(len(new_gene_indexes)):
			self.save_heatmap_data(
				h, i, trl_eff_index, exp_index, avg_new_gene_yield_rate[:, i],
				cell_mask)


	"""
	Capacity and General Gene Functions
	"""
	def get_capacity_gene_indexes(
			self, all_cells, index_type, capacity_gene_monomer_ids):
		"""
		Retrieve capacity gene indexes of a given type.

		Args:
			all_cells: Paths to all cells to read data from (directories should
				contain a simOut/ subdirectory), typically the return from
				AnalysisPaths.get_cells()
			index_type: Type of indexes to extract, currently supported options
				are 'cistron', 'RNA', 'mRNA', and 'monomer'
			capacity_gene_monomer_ids: monomer ids of capacity gene we need
				indexes for

		Returns:
			List of requested indexes
		"""
		sim_dir = all_cells[0]
		simOutDir = os.path.join(sim_dir, 'simOut')

		if index_type != "monomer":
			capacity_gene_mRNA_ids = self.get_mRNA_ids_from_monomer_ids(
				capacity_gene_monomer_ids)

		if index_type == 'cistron':
			# Extract cistron indexes for each new gene
			rnap_reader = TableReader(os.path.join(simOutDir, 'RnapData'))
			cistron_idx_dict = {
				cis: i for i, cis in
				enumerate(rnap_reader.readAttribute('cistron_ids'))}
			capacity_gene_indexes = [
				cistron_idx_dict.get(mRNA_id) for mRNA_id in
				capacity_gene_mRNA_ids]
		elif index_type == 'RNA':
			# Extract RNA indexes for each new gene
			rnap_reader = TableReader(os.path.join(simOutDir, 'RnaSynthProb'))
			RNA_idx_dict = {
				rna[:-3]: i for i, rna in
				enumerate(rnap_reader.readAttribute('rnaIds'))}
			capacity_gene_indexes = [
				RNA_idx_dict.get(mRNA_id) for mRNA_id in capacity_gene_mRNA_ids]
		elif index_type == 'mRNA':
			# Extract mRNA indexes for each new gene
			mRNA_counts_reader = TableReader(os.path.join(simOutDir, 'RNACounts'))
			mRNA_idx_dict = {
				rna[:-3]: i for i, rna in
				enumerate(mRNA_counts_reader.readAttribute('mRNA_ids'))}
			capacity_gene_indexes = [
				mRNA_idx_dict.get(mRNA_id[:-3]) for mRNA_id in capacity_gene_mRNA_ids]
		elif index_type == 'monomer':
			# Extract protein indexes for each new gene
			monomer_counts_reader = TableReader(
				os.path.join(simOutDir, "MonomerCounts"))
			monomer_idx_dict = {
				monomer: i for i, monomer in enumerate(
				monomer_counts_reader.readAttribute('monomerIds'))}
			capacity_gene_indexes = [
				monomer_idx_dict.get(monomer_id) for monomer_id in
				capacity_gene_monomer_ids]
		else:
			raise Exception(
				"Index type " + index_type +
				" has no instructions for data extraction.")

		return capacity_gene_indexes

	def get_avg_gene_counts(
			self, all_cells, data_table, data_column, capacity_gene_indexes):
		"""
		Retreives average counts of gene mRNAs or proteins, which are needed
		for multiple heatmaps.

		Args:
			all_cells: paths to all cells to read data from (directories should
				contain a simOut/ subdirectory), typically the return from
				AnalysisPaths.get_cells()
			data_table: Table to find data that needs to be retrieved
			data_column: Column to find data that needs to be retreived
			capacity_gene_indexes: Global indexes of the genes within data_table

		Returns:
			Average counts of gene mRNAs or proteins.
		"""
		return (read_stacked_columns(
				all_cells, data_table, data_column, fun=lambda
				x: np.mean( x[:, capacity_gene_indexes], axis=0)))

	def extract_capacity_gene_counts_heatmap_data(
			self, all_cells, h, trl_eff_index, exp_index, cell_mask,
			data_table, data_column, capacity_gene_index_type):
		"""
		Special function to handle extraction and saving of capacity gene mRNA
		and protein counts heatmap data.

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
			capacity_gene_index_type: Index type to use for the data table
		"""
		capacity_gene_indexes = self.get_capacity_gene_indexes(
			all_cells, capacity_gene_index_type, [capacity_gene_monomer_id])

		avg_capacity_gene_counts = self.get_avg_gene_counts(
			all_cells, data_table, data_column, capacity_gene_indexes)

		self.save_heatmap_data(h, 0, trl_eff_index, exp_index,
			np.log10(np.sum(avg_capacity_gene_counts, axis = 1) + 1), cell_mask)

	def extract_capacity_gene_mass_fraction_heatmap_data(
			self, all_cells, h, trl_eff_index, exp_index, cell_mask,
			counts_data_table, counts_data_column, mass_data_table,
			mass_data_column, capacity_gene_index_type):
		"""
		Special function to handle extraction and saving of capacity gene mRNA and
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
			capacity_gene_index_type: Index type to use for the data table
		"""
		capacity_gene_indexes = self.get_capacity_gene_indexes(
			all_cells, capacity_gene_index_type, [capacity_gene_monomer_id])

		if capacity_gene_index_type == 'monomer':
			capacity_gene_ids = [capacity_gene_monomer_id]
		else:
			capacity_gene_ids = list(
				self.get_mRNA_ids_from_monomer_ids([capacity_gene_monomer_id]))

		# Get mass for capacity gene (mRNAs or monomer)
		capacity_gene_masses = [1 for id in capacity_gene_ids]
		for i in range(len(capacity_gene_ids)):
			capacity_gene_masses[i] = (
				self.sim_data.getter.get_mass(
				capacity_gene_ids[i])/self.sim_data.constants.n_avogadro).asNumber(fg)

		# Get counts for each capacity gene (mRNAs or monomer)
		avg_capacity_gene_counts = self.get_avg_gene_counts(
			all_cells, counts_data_table, counts_data_column, capacity_gene_indexes)

		# Get average mass for all genes
		avg_mass = read_stacked_columns(
			all_cells, mass_data_table, mass_data_column,
			fun=lambda x: np.mean(x)).squeeze()

		# Determine mass fractions for each capacity gene (mRNAs or monomer)
		capacity_gene_mass_fractions = np.zeros_like((avg_capacity_gene_counts))
		for i in range(len(capacity_gene_ids)):
			capacity_gene_mass_fractions[:, i] = (avg_capacity_gene_counts[:, i] *
				capacity_gene_masses[i]) / avg_mass
		if len(capacity_gene_ids) > 1:
			capacity_gene_mass_fractions = np.sum(
				capacity_gene_mass_fractions, axis = 1)
		self.save_heatmap_data(
			h, 0, trl_eff_index, exp_index, capacity_gene_mass_fractions,
			cell_mask)

	def extract_capacity_gene_counts_fraction_heatmap_data(
			self, all_cells, h, trl_eff_index, exp_index, cell_mask,
			counts_data_table, counts_data_column, capacity_gene_index_type):
		"""
		Special function to handle extraction and saving of capacity gene mRNA and
		protein counts fraction heatmap data.

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
			capacity_gene_index_type: Index type to use for the data table
		"""
		capacity_gene_indexes = self.get_capacity_gene_indexes(
			all_cells, capacity_gene_index_type, [capacity_gene_monomer_id])

		if capacity_gene_index_type == 'monomer':
			capacity_gene_ids = [capacity_gene_monomer_id]
		else:
			capacity_gene_ids = list(
				self.get_mRNA_ids_from_monomer_ids([capacity_gene_monomer_id]))

		# Get counts for each new gene
		avg_capacity_gene_counts = self.get_avg_gene_counts(
			all_cells, counts_data_table, counts_data_column, capacity_gene_indexes)
		if len(capacity_gene_ids) > 1:
			avg_capacity_gene_counts = np.sum(
				avg_capacity_gene_counts, axis = 1)

		# Get total avg counts for all genes
		total_counts = np.sum(read_stacked_columns(
			all_cells, counts_data_table, counts_data_column,
			fun=lambda x: np.mean(x, axis=0)), axis=1)

		capacity_gene_counts_fraction = avg_capacity_gene_counts / total_counts
		self.save_heatmap_data(
			h, 0, trl_eff_index, exp_index, capacity_gene_counts_fraction,
			cell_mask)

	"""
	Plotting Functions
	"""
	def plot_heatmaps(
			self, is_dashboard, variant_mask, heatmap_x_label, heatmap_y_label,
			new_gene_expression_factors, new_gene_translation_efficiency_values,
			summary_statistic, figsize_x, figsize_y, plotOutDir, plot_suffix):
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
			summary_statistic: Specifies whether average ('mean') or
				standard deviation ('std_dev') should be displayed on the
				heatmaps
			figsize_x: Horizontal size of each heatmap
			figsize_y: Vertical size of each heatmap
			plotOutDir: Output directory for plots
			plot_suffix: Suffix to add to plot file names, usually specifying
				which generations were plotted
		"""
		if summary_statistic == 'std_dev':
			plot_suffix = plot_suffix + "_std_dev"
		elif summary_statistic != 'mean':
			raise Exception(
				"'mean' and 'std_dev' are the only currently supported"
				" summary statistics")

		if is_dashboard:
			# Determine dashboard layout
			if self.total_heatmaps_to_make > 3:
				dashboard_ncols = 4
				dashboard_nrows = math.ceil((self.total_heatmaps_to_make + 1) / dashboard_ncols)
			else:
				dashboard_ncols = self.total_heatmaps_to_make + 1
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
				f"Percentage of Sims That Reached Generation {COUNT_INDEX + 1}")
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
							title_addition = f": {self.new_gene_mRNA_ids[i][:-4]}"
						self.make_single_heatmap(
							h, axs[row_ax, col_ax], variant_mask,
							heatmap_x_label, heatmap_y_label, i,
							new_gene_expression_factors,
							new_gene_translation_efficiency_values,
							summary_statistic, title_addition)
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
								new_gene_translation_efficiency_values,
								summary_statistic,
								ntp_id)
							fig.tight_layout()
							plt.show()
							exportFigure(
								plt, plotOutDir,
								f'new_gene_mRNA_{ntp_id[:-3]}_fraction_heatmap'
								f'_{self.new_gene_mRNA_ids[i][:-4]}'
								f'{plot_suffix}')
							col_ax += 1
							if (col_ax == dashboard_ncols):
								col_ax = 0
								row_ax += 1
				else:
					raise Exception(
						f"Heatmap {h} is neither a standard plot nor a"
						f" nonstandard plot that has specific instructions for"
						f" plotting.")
			fig.tight_layout()
			exportFigure(plt, plotOutDir,
				f"00_new_gene_exp_trl_eff_dashboard{plot_suffix}")
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
				f"Percentage of Sims that Reached Generation {COUNT_INDEX + 1}")
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
							title_addition = f": {self.new_gene_mRNA_ids[i][:-4]}"
							filename_addition = f"_{self.new_gene_mRNA_ids[i][:-4]}"
						fig, ax = plt.subplots(1, 1, figsize=(figsize_x, figsize_y))
						self.make_single_heatmap(
							h, ax, variant_mask, heatmap_x_label, heatmap_y_label,
							i, new_gene_expression_factors,
							new_gene_translation_efficiency_values,
							summary_statistic, title_addition)
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
								new_gene_translation_efficiency_values,
								summary_statistic, ntp_id)
							fig.tight_layout()
							plt.show()
							exportFigure(
								plt, plotOutDir,
								f'new_gene_mRNA_{ntp_id[:-3]}_fraction_heatmap'
								f'_{self.new_gene_mRNA_ids[i][:-4]}{plot_suffix}')
				else:
					raise Exception(
						f"Heatmap {h} is neither a standard plot nor a"
						f" nonstandard plot that has specific instructions for"
						f" plotting.")

	def make_single_heatmap(
			self, h, ax, variant_mask, heatmap_x_label, heatmap_y_label,
			initial_index, new_gene_expression_factors,
			new_gene_translation_efficiency_values,
			summary_statistic, title_addition):
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
			summary_statistic: Specifies whether average ('mean') or
				standard deviation ('std_dev') should be displayed on the
				heatmaps
			title_addition: Any string that needs to be added to the title of
				the heatmap, e.g. a new gene id
		"""
		title = self.heatmap_details[h]['plot_title'] + title_addition
		if summary_statistic == "std_dev":
			title = f"Std Dev: {title}"

		heatmap(
			self, ax, variant_mask,
			self.heatmap_data[h][summary_statistic][initial_index, :, :],
			self.heatmap_data["completed_gens_heatmap"][0, :, :],
			new_gene_expression_factors, new_gene_translation_efficiency_values,
			heatmap_x_label, heatmap_y_label,
			title,
			self.heatmap_details[h]['box_text_size'])

	def make_new_gene_mRNA_NTP_fraction_heatmap(
			self, h, ax, variant_mask, heatmap_x_label, heatmap_y_label,
			initial_index, new_gene_expression_factors,
			new_gene_translation_efficiency_values, summary_statistic, ntp_id):
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
			summary_statistic: Specifies whether average ('mean') or
				standard deviation ('std_dev') should be displayed on the
				heatmaps
			ntp_id: Id of NTP to plot
		"""
		title = (f"{self.heatmap_details[h]['plot_title']} {ntp_id[:-3]}"
				 f" Fraction: {self.new_gene_mRNA_ids[initial_index][:-4]}")
		if summary_statistic == 'std_dev':
			title = f"Std Dev: {title}"

		heatmap(
			self, ax, variant_mask,
			self.heatmap_data[h][summary_statistic][ntp_id][initial_index, :, :],
			self.heatmap_data["completed_gens_heatmap"][0, :, :],
			new_gene_expression_factors, new_gene_translation_efficiency_values,
			heatmap_x_label, heatmap_y_label,
			title,
			self.heatmap_details[h]['box_text_size'])

	"""
	Orchestrating Function
	"""
	def do_plot(self, inputDir, plotOutDir, plotOutFileName, simDataFile,
				validationDataFile, metadata):

		with open(simDataFile, 'rb') as f:
			self.sim_data = pickle.load(f)

		# Determine new gene mRNA and monomer ids
		(self.new_gene_mRNA_ids, self.new_gene_indices, self.new_gene_monomer_ids,
			self.new_gene_monomer_indices) = determine_new_gene_ids_and_indices(self.sim_data)

		# Map variant indexes to parameters used for new genes
		assert ('new_gene_expression_factors' in metadata and
				'new_gene_translation_efficiency_values' in metadata), (
			"This plot is intended to be run on simulations where the"
			" new gene expression-translation efficiency variant was "
			"enabled, but no parameters for this variant were found.")
		new_gene_expression_factors = metadata['new_gene_expression_factors']
		new_gene_translation_efficiency_values = metadata[
			'new_gene_translation_efficiency_values']

		variants = self.ap.get_variants()
		variant_index_to_values = {}
		variant_index_to_list_indices = {}
		variant_mask = np.zeros((  # Track whether we ran this sim
			len(new_gene_translation_efficiency_values),
			len(new_gene_expression_factors)), dtype=bool)
		for index in variants:
			condition_index = index // 1000
			index_remainder = index - condition_index * 1000

			(expression_list_index, trl_eff_list_index, expression_variant_index,
			trl_eff_value) = get_new_gene_expression_factor_and_translation_efficiency(
				self.sim_data, index_remainder)

			expression_variant_index = new_gene_expression_factors[expression_list_index]
			trl_eff_value = new_gene_translation_efficiency_values[trl_eff_list_index]

			variant_index_to_values[index] = np.array([
				expression_variant_index, trl_eff_value])
			variant_index_to_list_indices[index] = np.array([
				expression_list_index, trl_eff_list_index])
			variant_mask[trl_eff_list_index, expression_list_index] = True

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
			"free_rnap_counts_heatmap":
				{'is_nonstandard_data_retrieval': True,
				 'num_digits_rounding': 0,
				 'box_text_size': 'x-small',
				 'plot_title': 'Free RNA Polymerase (RNAP) Counts',
				 },
			"free_ribosome_counts_heatmap":
				{'is_nonstandard_data_retrieval': True,
				 'data_table': 'UniqueMoleculeCounts',
				 'data_column': 'uniqueMoleculeCounts',
				 'num_digits_rounding': 0,
				 'box_text_size': 'x-small',
				 'plot_title': 'Free Ribosome Counts',
				 },
			"rnap_ribosome_counts_ratio_heatmap":
				{'is_nonstandard_data_retrieval': True,
				 'num_digits_rounding': 4,
				 'box_text_size': 'x-small',
				 'plot_title': 'RNAP Counts / Ribosome Counts',
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
			"weighted_avg_translation_efficiency_heatmap":
				{'is_nonstandard_data_retrieval': True,
				 'num_digits_rounding': 3,
				 'plot_title': 'Translation Efficiency (Weighted Average)',
				},
			"capacity_gene_mRNA_counts_heatmap":
				{'is_nonstandard_data_retrieval': True,
				 'plot_title': 'Log(Capacity Gene mRNA Counts+1): '
					+ capacity_gene_common_name},
			"capacity_gene_monomer_counts_heatmap":
				{'is_nonstandard_data_retrieval': True,
				 'plot_title': 'Log(Capacity Gene Protein Counts+1): '
					+ capacity_gene_common_name},
			"capacity_gene_mRNA_mass_fraction_heatmap":
				{'is_nonstandard_data_retrieval': True,
				 'num_digits_rounding': 3,
				 'plot_title': 'Capacity Gene mRNA Mass Fraction: '
					+ capacity_gene_common_name},
			"capacity_gene_monomer_mass_fraction_heatmap":
				{'is_nonstandard_data_retrieval': True,
				 'num_digits_rounding': 3,
				 'plot_title': 'Capacity Gene Protein Mass Fraction: '
					+ capacity_gene_common_name},
			"capacity_gene_mRNA_counts_fraction_heatmap":
				{'is_nonstandard_data_retrieval': True,
				 'num_digits_rounding': 3,
				 'plot_title': 'Capacity Gene mRNA Counts Fraction: '
							   + capacity_gene_common_name},
			"capacity_gene_monomer_counts_fraction_heatmap":
				{'is_nonstandard_data_retrieval': True,
				 'num_digits_rounding': 3,
				 'plot_title': 'Capacity Gene Protein Counts Fraction: '
							   + capacity_gene_common_name},
			"new_gene_mRNA_counts_heatmap":
				{'plot_title': 'Log(New Gene mRNA Counts+1)'},
			"new_gene_monomer_counts_heatmap":
				{'plot_title': 'Log(New Gene Protein Counts+1)'},
			"new_gene_mRNA_mass_fraction_heatmap":
				{'plot_title': 'New Gene mRNA Mass Fraction'},
			"new_gene_mRNA_counts_fraction_heatmap":
				{'plot_title': 'New Gene mRNA Counts Fraction'},
			"new_gene_mRNA_NTP_fraction_heatmap":
				{'is_nonstandard_plot': True,
				 'num_digits_rounding': 4,
				 'box_text_size': 'x-small',
				 'plot_title': 'New Gene',
				},
			"new_gene_monomer_mass_fraction_heatmap":
				{'plot_title': 'New Gene Protein Mass Fraction'},
			"new_gene_monomer_counts_fraction_heatmap":
				{'plot_title': 'New Gene Protein Counts Fraction'},
			"new_gene_rnap_init_rate_heatmap":
				{'plot_title': 'New Gene RNAP Initialization Rate'},
			"new_gene_ribosome_init_rate_heatmap":
				{'plot_title': 'New Gene Ribosome Initalization Rate'},
			"new_gene_rnap_time_overcrowded_heatmap":
				{'plot_title': 'Fraction of Time RNAP Overcrowded New Gene'},
			"new_gene_ribosome_time_overcrowded_heatmap":
				{'plot_title': 'Fraction of Time Ribosome Overcrowded New Gene'},
			"new_gene_actual_protein_init_prob_heatmap":
				{'plot_title': 'New Gene Actual Protein Init Prob',
				 'num_digits_rounding': 4},
			"new_gene_target_protein_init_prob_heatmap":
				{'plot_title': 'New Gene Target Protein Init Prob',
				 'num_digits_rounding': 4},
			"new_gene_actual_rna_synth_prob_heatmap":
				{'plot_title': 'New Gene Actual RNA Synth Prob',
				 'num_digits_rounding': 4},
			"new_gene_target_rna_synth_prob_heatmap":
				{'plot_title': 'New Gene Target RNA Synth Prob',
				 'num_digits_rounding': 4},
			"new_gene_rnap_counts_heatmap":
				{'box_text_size': 'x-small',
				 'num_digits_rounding': 0,
				 'plot_title': 'New Gene RNAP Counts',
				},
			"new_gene_rnap_portion_heatmap":
				{'plot_title': 'New Gene RNAP Portion',
				 'num_digits_rounding': 3,
				},
			"rrna_rnap_counts_heatmap":
				{'is_nonstandard_data_retrieval': True,
				 'plot_title': 'rRNA RNAP Counts',
				 'num_digits_rounding': 0,
				},
			"rrna_rnap_portion_heatmap":
				{'is_nonstandard_data_retrieval': True,
				 'plot_title': 'rRNA RNAP Portion',
				 'num_digits_rounding': 3,
				 },
			"rnap_subunit_rnap_portion_heatmap":
				{'is_nonstandard_data_retrieval': True,
				 'plot_title': 'RNAP Subunit RNAP Portion',
				 'num_digits_rounding': 3,
				},
			"rnap_subunit_ribosome_portion_heatmap":
				{'is_nonstandard_data_retrieval': True,
				 'plot_title': 'RNAP Subunit Ribosome Portion',
				 'num_digits_rounding': 3,
				 },
			"ribosomal_protein_rnap_portion_heatmap":
				{'is_nonstandard_data_retrieval': True,
				 'plot_title': 'Ribosomal Protein RNAP Portion',
				 'num_digits_rounding': 3,
				 },
			"ribosomal_protein_ribosome_portion_heatmap":
				{'is_nonstandard_data_retrieval': True,
				 'plot_title': 'Ribosomal Protein Ribosome Portion',
				 'num_digits_rounding': 3,
				 },
			"new_gene_ribosome_counts_heatmap":
				{'box_text_size': 'x-small',
				 'num_digits_rounding': 0,
				 'plot_title': 'New Gene Ribosome Counts',
				},
			"new_gene_ribosome_portion_heatmap":
				{'plot_title': 'New Gene Ribosome Portion',
				 'num_digits_rounding': 3,
				},
			"capacity_gene_rnap_portion_heatmap":
				{'is_nonstandard_data_retrieval': True,
				 'plot_title':
				 'Capacity Gene RNAP Portion: ' + capacity_gene_common_name,
				 'num_digits_rounding': 4,
				 },
			"capacity_gene_ribosome_portion_heatmap":
				{'is_nonstandard_data_retrieval': True,
				 'plot_title':
				 'Capacity Gene Ribosome Portion: ' + capacity_gene_common_name,
				 'num_digits_rounding': 4,
				 },
			"new_gene_yield_per_glucose":
				{'is_nonstandard_data_retrieval': True,
				 'plot_title': 'New Gene fg Protein Yield per fg Glucose',
				 'num_digits_rounding': 3,
				 },
			"new_gene_yield_per_hour":
				{'is_nonstandard_data_retrieval': True,
				 'plot_title': 'New Gene fg Protein Yield per Hour',
				 'num_digits_rounding': 2,
				 },
			"glucose_consumption_rate":
				{'is_nonstandard_data_retrieval': True,
				 'plot_title': 'Average Glucose Consumption Rate (fg/hr)',
				 'num_digits_rounding': 1,
				 },
		}

		# Check validity of requested heatmaps and fill in default values where needed
		heatmaps_to_make = set(HEATMAPS_TO_MAKE_LIST)
		assert "completed_gens_heatmap" not in heatmaps_to_make, \
			"the completed_gens_heatmap is run by default, do not include in heatmaps_to_make"
		self.total_heatmaps_to_make = 0
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
			if not h.startswith("new_gene_"):
				self.total_heatmaps_to_make += 1
			elif h == "new_gene_mRNA_NTP_fraction_heatmap":
				self.ntp_ids = list(
					self.sim_data.ntp_code_to_id_ordered.values())
				self.total_heatmaps_to_make += len(self.ntp_ids)
			else:
				self.total_heatmaps_to_make += len(self.new_gene_mRNA_ids)

		# Create data structures to use for the heatmaps
		self.heatmap_data = {}
		self.heatmap_data["completed_gens_heatmap"] = np.zeros((
			1, len(new_gene_translation_efficiency_values),
			len(new_gene_expression_factors)))
		for h in heatmaps_to_make:
			if not self.heatmap_details[h]['is_new_gene_heatmap']:
				self.heatmap_data[h] = {}
				self.heatmap_data[h]["mean"] = np.zeros((
					1, len(new_gene_translation_efficiency_values),
					len(new_gene_expression_factors))
					) + self.heatmap_details[h]['default_value']
				self.heatmap_data[h]["std_dev"] = np.zeros((
					1, len(new_gene_translation_efficiency_values),
					len(new_gene_expression_factors))
				) + self.heatmap_details[h]['default_value']
			else:
				if h == "new_gene_mRNA_NTP_fraction_heatmap":
					self.heatmap_data[h] = {}
					self.heatmap_data[
						"new_gene_mRNA_NTP_fraction_heatmap"]["mean"] = {}
					self.heatmap_data[
						"new_gene_mRNA_NTP_fraction_heatmap"]["std_dev"] = {}
					for ntp_id in self.ntp_ids:
						self.heatmap_data[
							"new_gene_mRNA_NTP_fraction_heatmap"][
							"mean"][ntp_id] = np.zeros(
							(len(self.new_gene_mRNA_ids),
							 len(new_gene_translation_efficiency_values),
							 len(new_gene_expression_factors))
							) + self.heatmap_details[h]['default_value']
						self.heatmap_data[
							"new_gene_mRNA_NTP_fraction_heatmap"][
							"std_dev"][ntp_id] = np.zeros(
							(len(self.new_gene_mRNA_ids),
							 len(new_gene_translation_efficiency_values),
							 len(new_gene_expression_factors))
							) + self.heatmap_details[h]['default_value']
				else:
					self.heatmap_data[h] = {}
					self.heatmap_data[h]["mean"] = np.zeros((
						len(self.new_gene_mRNA_ids),
						len(new_gene_translation_efficiency_values),
						len(new_gene_expression_factors))
						) + self.heatmap_details[h]['default_value']
					self.heatmap_data[h]["std_dev"] = np.zeros((
						len(self.new_gene_mRNA_ids),
						len(new_gene_translation_efficiency_values),
						len(new_gene_expression_factors))
					) + self.heatmap_details[h]['default_value']

		# Data extraction
		print("---Data Extraction---")
		reached_count_gen = {}
		generations = {}
		variants = self.ap.get_variants()
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
		figsize_x =  2 + 2*len(new_gene_expression_factors)/3
		figsize_y = 2*len(new_gene_translation_efficiency_values)/2

		# Create dashboard plot
		if DASHBOARD_FLAG == 1 or DASHBOARD_FLAG == 2:
			self.plot_heatmaps(
				True, variant_mask, heatmap_x_label, heatmap_y_label,
				new_gene_expression_factors,
				new_gene_translation_efficiency_values, 'mean', figsize_x,
				figsize_y, plotOutDir, plot_suffix)

			if STD_DEV_FLAG:
				self.plot_heatmaps(
					True, variant_mask, heatmap_x_label, heatmap_y_label,
					new_gene_expression_factors,
					new_gene_translation_efficiency_values, 'std_dev',
					figsize_x, figsize_y, plotOutDir, plot_suffix)

		# Create separate plots
		if DASHBOARD_FLAG == 0 or DASHBOARD_FLAG == 2:
			self.plot_heatmaps(
				False, variant_mask, heatmap_x_label,heatmap_y_label,
				new_gene_expression_factors,
				new_gene_translation_efficiency_values, 'mean', figsize_x,
				figsize_y, plotOutDir, plot_suffix)

			if STD_DEV_FLAG:
				self.plot_heatmaps(
					False, variant_mask, heatmap_x_label, heatmap_y_label,
					new_gene_expression_factors,
					new_gene_translation_efficiency_values, 'std_dev',
					figsize_x, figsize_y, plotOutDir, plot_suffix)

if __name__ == "__main__":
	Plot().cli()
