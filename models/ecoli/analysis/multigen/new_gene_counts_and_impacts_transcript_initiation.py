"""
Plot mRNA and protein counts for new genes across multiple generations, as well
as plots to analyze the impact of new gene expression, including growth rate,
RNAP and ribosome counts, and ppGpp concentration.
"""
# TODO: update file header comment

import pickle
import os

from matplotlib import pyplot as plt
import matplotlib as mpl
# noinspection PyUnresolvedReferences
import mplcursors
import numpy as np
from numpy import inf

from models.ecoli.analysis import multigenAnalysisPlot
from models.ecoli.sim.variants.new_gene_internal_shift import determine_new_gene_ids_and_indices
from wholecell.analysis.analysis_tools import (exportFigure,
	read_stacked_bulk_molecules, read_stacked_columns, read_bulk_molecule_counts)
from wholecell.io.tablereader import TableReader
from wholecell.utils import units

LINE_COLOR = (66/255, 170/255, 154/255)

# OTHER_GENES_OF_INTEREST = [
# 	"EG11024 ", # trpA, TRYPSYN-APROTEIN[c]
# 	"EG11025", # trpB, TRYPSYN-BPROTEIN[c]
# 	"EG10447", # hisD, HISTDEHYD-MONOMER[c]
# 	"EG11000 ", # thrC, THRESYN-MONOMER[c]
# ]

OTHER_MONOMERS_OF_INTEREST = [
	"TRYPSYN-APROTEIN[c]",
	"TRYPSYN-BPROTEIN[c]",
	"HISTDEHYD-MONOMER[c]",
	"THRESYN-MONOMER[c]",
]

OTHER_GENE_COMMON_NAMES = [ # TODO: extract these using sim_data
	"trpA",
	"trpB",
	"hisD",
	"thrC",
]

OTHER_COMPLEXES_OF_INTEREST = ["TRYPSYN[c]"]

# TODO: also plot tryptophan synthase complex? ID: TRYPSYN

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
		mRNA_sim_data = sim_data.process.transcription.cistron_data.struct_array
		monomer_sim_data = sim_data.process.translation.monomer_data.struct_array
		new_gene_mRNA_ids = mRNA_sim_data[mRNA_sim_data['is_new_gene']]['id'].tolist()
		mRNA_monomer_id_dict = dict(zip(
			monomer_sim_data['cistron_id'], monomer_sim_data['id']))
		monomer_to_mRNA_id_dict = dict(zip(
			monomer_sim_data['id'],monomer_sim_data['cistron_id']))
		new_gene_monomer_ids = [mRNA_monomer_id_dict.get(mRNA_id)
								for mRNA_id in new_gene_mRNA_ids]
		if len(new_gene_mRNA_ids) == 0:
			print("This plot is intended to be run on simulations where the"
				  " new gene option was enabled, but no new gene mRNAs were "
				  "found.")
			return
		if len(new_gene_monomer_ids) == 0:
			print("This plot is intended to be run on simulations where the "
				  "new gene option was enabled, but no new gene proteins "
				  "were "
				  "found.")
			return
		assert len(new_gene_monomer_ids) == len(new_gene_mRNA_ids), \
			'number of new gene monomers and mRNAs should be equal'

		# TODO: make a option for this plot that you could run without new genes

		other_gene_mRNA_ids = {}
		other_gene_mRNA_indexes = {}
		for monomer in OTHER_MONOMERS_OF_INTEREST:
			other_gene_mRNA_ids[monomer] = list(self.get_mRNA_ids_from_monomer_ids(
				sim_data, [monomer]))
			other_gene_mRNA_indexes[monomer] = self.get_mRNA_indexes_from_monomer_ids(
				sim_data, cell_paths, [monomer], 'mRNA')

		# determine cistron ids from otehr_monomers_of_interest
		other_gene_cistron_ids = [
			monomer_to_mRNA_id_dict.get(monomer_id) for monomer_id
			in OTHER_MONOMERS_OF_INTEREST]
		new_gene_cistron_ids = [
			monomer_to_mRNA_id_dict.get(monomer_id) for monomer_id in new_gene_monomer_ids]
		# Extract cistron indexes for each gene of interest
		cistron_counts_reader = TableReader(os.path.join(simOutDir, 'RNACounts'))
		cistron_idx_dict = {cistron: i for i, cistron in enumerate(
			cistron_counts_reader.readAttribute('mRNA_cistron_ids'))}
		other_gene_cistron_indexes = [
			cistron_idx_dict.get(cistron_id) for cistron_id in other_gene_cistron_ids]
		cistron_counts_reader.close()

		# Extract mRNA indexes for each new gene
		mRNA_counts_reader = TableReader(os.path.join(simOutDir,
													  'RNACounts'))
		mRNA_idx_dict = {rna[:-3]: i for i, rna in enumerate(
			mRNA_counts_reader.readAttribute('mRNA_ids'))}
		new_gene_mRNA_indexes = [mRNA_idx_dict.get(mRNA_id) for mRNA_id in
								 new_gene_mRNA_ids]
		mRNA_counts_reader.close()

		# Extract protein indexes for each new gene
		monomer_counts_reader = TableReader(
			os.path.join(simOutDir, "MonomerCounts"))
		monomer_idx_dict = {monomer: i for i, monomer in
							enumerate(monomer_counts_reader.readAttribute(
								'monomerIds'))}
		new_gene_monomer_indexes = [monomer_idx_dict.get(monomer_id) for
									monomer_id in new_gene_monomer_ids]
		other_gene_monomer_indexes = [monomer_idx_dict.get(monomer_id) for
									  monomer_id in OTHER_MONOMERS_OF_INTEREST]
		monomer_counts_reader.close()

		# Extract RNA indexes for each new gene
		rnap_reader = TableReader(os.path.join(simOutDir, 'RnaSynthProb'))
		RNA_idx_dict = {
			rna[:-3]: i for i, rna in
			enumerate(rnap_reader.readAttribute('rnaIds'))}
		new_gene_RNA_indexes = [
			RNA_idx_dict.get(mRNA_id) for mRNA_id in new_gene_mRNA_ids]
		rna_ids = rnap_reader.readAttribute('rnaIds')
		rnap_reader.close()

		# Load data
		time = read_stacked_columns(
			cell_paths, 'Main', 'time', ignore_exception=True)
		time_no_first = read_stacked_columns(
			cell_paths, 'Main', 'time', remove_first=True, ignore_exception=True)
		(new_gene_monomer_counts,) = read_stacked_bulk_molecules(
			cell_paths, new_gene_monomer_ids, ignore_exception=True)
		all_mRNA_stacked_counts = read_stacked_columns(
			cell_paths, 'RNACounts', 'mRNA_counts', ignore_exception=True)
		new_gene_mRNA_counts = all_mRNA_stacked_counts[:,new_gene_mRNA_indexes]
		new_gene_promoter_copy_numbers = read_stacked_columns(
			cell_paths, 'RnaSynthProb', 'promoter_copy_number',
			ignore_exception=True)[:,new_gene_RNA_indexes]

		# Other genes of interest
		other_gene_monomer_counts = read_stacked_columns(
			cell_paths, 'MonomerCounts', 'monomerCounts', ignore_exception=True)[:,other_gene_monomer_indexes]
		all_mRNA_cistron_stacked_counts = read_stacked_columns(
			cell_paths, 'RNACounts', 'mRNA_cistron_counts', ignore_exception=True)
		other_gene_cistron_counts = all_mRNA_cistron_stacked_counts[:,other_gene_cistron_indexes]

		rna_synth_prob_max_p = read_stacked_columns(
			cell_paths, 'RnaSynthProb', 'max_p',
			ignore_exception=True, remove_first=True)

		# Transcript Initiation
		basal_probs_orig = read_stacked_columns(
			cell_paths, 'RnaSynthProb', 'basal_prob_orig',
			ignore_exception=True, remove_first=True)
		basal_probs_updated = read_stacked_columns(
			cell_paths, 'RnaSynthProb', 'basal_prob_updated',
			ignore_exception=True, remove_first=True)

		basal_probs_orig_no_gfp = read_stacked_columns(
			[cell_paths[0]], 'RnaSynthProb', 'basal_prob_orig',
			ignore_exception=True, remove_first=True)
		basal_probs_orig_yes_gfp = read_stacked_columns(
			[cell_paths[1]], 'RnaSynthProb', 'basal_prob_orig',
			ignore_exception=True, remove_first=True)
		basal_probs_updated_no_gfp = read_stacked_columns(
			[cell_paths[0]], 'RnaSynthProb', 'basal_prob_updated',
			ignore_exception=True, remove_first=True)
		basal_probs_updated_yes_gfp = read_stacked_columns(
			[cell_paths[1]], 'RnaSynthProb', 'basal_prob_updated',
			ignore_exception=True, remove_first=True)

		# if basal_probs_updated_no_gfp and basal_probs_updated_yes_gfp are not the same shape, trim rows off of the longer one
		if basal_probs_updated_no_gfp.shape[0] > basal_probs_updated_yes_gfp.shape[0]:
			basal_probs_updated_no_gfp = basal_probs_updated_no_gfp[:basal_probs_updated_yes_gfp.shape[0], :]
		elif basal_probs_updated_no_gfp.shape[0] < basal_probs_updated_yes_gfp.shape[0]:
			basal_probs_updated_yes_gfp = basal_probs_updated_yes_gfp[:basal_probs_updated_no_gfp.shape[0], :]

		# find all columns where the basal_prob_updated_yes_gfp is greater than the basal_prob_updated_no_gfp for any row
		threshold = 0.70
		rna_idx_where_updated_yes_gfp_greater = np.where(np.any(basal_probs_updated_yes_gfp >= basal_probs_updated_no_gfp * threshold, axis=0))[0]
		counts_rna_idx_where_updated_yes_gfp_greater = np.sum(basal_probs_updated_yes_gfp >= basal_probs_updated_no_gfp * threshold, axis=0)

		# get rna ids for all rna_idx_where_updated_yes_gfp_greater where counts_rna_idx_where_updated_yes_gfp_greater is greater than 2000
		rna_ids_where_updated_yes_gfp_greater = [sim_data.process.transcription.rna_data['id'][rna_idx] for rna_idx in rna_idx_where_updated_yes_gfp_greater if counts_rna_idx_where_updated_yes_gfp_greater[rna_idx] > 2000]

		# print(rna_idx_where_updated_yes_gfp_greater)
		# print(rna_ids_where_updated_yes_gfp_greater)

		# # make a dictionary
		# yes_gfp_greater_dict = {}
		# for j in range(len(rna_idx_where_updated_yes_gfp_greater)):
		# 	yes_gfp_greater_dict[rna_ids_where_updated_yes_gfp_greater[j]] = rna_idx_where_updated_yes_gfp_greater[j]

		# print(yes_gfp_greater_dict)
		# import ipdb; ipdb.set_trace()

		# map from rna_ids to cistrons
		cistron_ids_where_updated_yes_gfp_greater = np.array([])
		for rna_id in rna_ids_where_updated_yes_gfp_greater:
			cistron_indices = sim_data.process.transcription.rna_id_to_cistron_indexes(rna_id)
			cistron_ids_where_updated_yes_gfp_greater = np.append(cistron_ids_where_updated_yes_gfp_greater, sim_data.process.transcription.cistron_data[cistron_indices]['id'])

		# map these rna_ids to protein ids
		cistron_ids_to_monomer_ids = {monomer_data['cistron_id']: monomer_data['id'] for monomer_data in sim_data.process.translation.monomer_data}

		# remove RNA0-141 from cistron_ids_where_updated_yes_gfp_greater
		cistron_ids_where_updated_yes_gfp_greater = cistron_ids_where_updated_yes_gfp_greater[cistron_ids_where_updated_yes_gfp_greater != 'RNA0-141']

		monomer_ids_where_updated_yes_gfp_greater = []
		for cistron_id in cistron_ids_where_updated_yes_gfp_greater:
			if cistron_id in cistron_ids_to_monomer_ids.keys():
				monomer_ids_where_updated_yes_gfp_greater.append(cistron_ids_to_monomer_ids[cistron_id])

		epsilon = 1E-12
		basal_probs_updated_no_gfp_to_plot = basal_probs_updated_no_gfp.mean(axis=0) + epsilon
		basal_probs_updated_yes_gfp_to_plot = basal_probs_updated_yes_gfp.mean(axis=0) + epsilon
		# Remove GFP
		basal_probs_updated_no_gfp_to_plot = basal_probs_updated_no_gfp_to_plot[:-1]
		basal_probs_updated_yes_gfp_to_plot = basal_probs_updated_yes_gfp_to_plot[:-1]

		fold_changes = np.log10(
			basal_probs_updated_yes_gfp_to_plot / basal_probs_updated_no_gfp_to_plot)

		plt.figure(figsize=(16, 16))
		plt.scatter(np.log10(basal_probs_updated_no_gfp_to_plot), fold_changes)
		plt.xlabel('log10(basal_prob_updated_no_gfp)')
		plt.ylabel('log10 Fold Change')
		plt.title('Basal Transcription Probability: No GFP vs Yes GFP')

		# add the rna id labels
		for i, txt in enumerate(rna_ids[:-1]):
			if txt in ["TU3[c]", "TU00285[c]", "TU00178[c]"]:
				plt.annotate(txt, (np.log10(basal_probs_updated_no_gfp_to_plot[i]), fold_changes[i]))

		exportFigure(plt, plotOutDir, plotOutFileName, metadata)

		import plotly.express as px
		import plotly.io as pio

		# Create interactive plot with Plotly
		df = {
			'log10(basal_prob_updated_no_gfp)': np.log10(basal_probs_updated_no_gfp_to_plot),
			'log10 Fold Change': fold_changes,
			'basal_prob_updated_no_gfp': basal_probs_updated_no_gfp.mean(axis=0)[:-1],
			'basal_prob_updated_yes_gfp': basal_probs_updated_yes_gfp.mean(axis=0)[:-1],
			'rna_ids': rna_ids[:-1],
			}
		fig = px.scatter(df, x='log10(basal_prob_updated_no_gfp)', y='log10 Fold Change',
						 title='Basal Probability: No GFP vs Yes GFP',
						 labels={
							 'log10(basal_prob_updated_no_gfp)': 'log10(basal_prob_updated_no_gfp)',
							 'log10 Fold Change': 'log10 Fold Change'},
						 hover_data={
							 'basal_prob_updated_no_gfp': True,
							 'basal_prob_updated_yes_gfp': True,
							 'rna_ids': True,
						})

		# Save as HTML
		pio.write_html(fig, os.path.join(plotOutDir, plotOutFileName + '_interactive_scatterplot.html'))

if __name__ == '__main__':
	Plot().cli()
