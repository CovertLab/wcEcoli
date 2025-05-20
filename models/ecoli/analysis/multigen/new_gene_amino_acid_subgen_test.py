import pickle
import os
from cmath import phase

from matplotlib import pyplot as plt
import matplotlib as mpl
from matplotlib.collections import PatchCollection
import numpy as np

from models.ecoli.analysis import multigenAnalysisPlot
from wholecell.analysis.analysis_tools import (exportFigure,
	read_stacked_bulk_molecules, read_stacked_columns)
from wholecell.io.tablereader import TableReader

START_GEN = 0
END_GEN = 1

LINE_COLOR = (66/255, 170/255, 154/255)



genes_of_interest = {}

genes_of_interest["argH"] = {
	"gene_name": "argH",
	"cistron_of_interest": "EG11190_RNA",
	"reaction_ids_of_interest": ['ARGSUCCINLYA-RXN'],
	"reaction_product_of_interest": "ARG[c]",}

genes_of_interest["leuD"] = {
	"gene_name": "leuD",
	"cistron_of_interest": "EG11575_RNA",
	"reaction_ids_of_interest": ['RXN-13163'],
	"reaction_product_of_interest": "3-CARBOXY-3-HYDROXY-ISOCAPROATE[c]",}
genes_of_interest["argA"] = {
	"gene_name": "argA",
	"cistron_of_interest": "EG10063_RNA",
	"reaction_ids_of_interest": ['N-ACETYLTRANSFER-RXN'],
	"reaction_product_of_interest": "ACETYL-GLU[c]",
	"supplementary_gene_name": "argB",
	"supplementary_cistron_of_interest": "EG10064_RNA",
	"supplementary_reaction_ids_of_interest": ["ACETYLGLUTKIN-RXN"],
	"supplementary_reaction_product_of_interest": "N-ACETYL-GLUTAMYL-P[c]",
}

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

		cell_paths = self.ap.get_cells(generation=np.arange(START_GEN, END_GEN))
		sim_dir = cell_paths[0]
		simOutDir = os.path.join(sim_dir, 'simOut')

		reaction_dict = sim_data.process.metabolism.reaction_stoich

		for gene, gene_data in genes_of_interest.items():
			gene_name = gene_data["gene_name"]
			cistron_of_interest = gene_data["cistron_of_interest"]
			reaction_ids_of_interest = gene_data["reaction_ids_of_interest"]

			reaction_product_of_interest = gene_data["reaction_product_of_interest"]
			if "supplementary_gene_name" in gene_data.keys():
				supplementary_gene_name = gene_data["supplementary_gene_name"]
				supplementary_cistron_of_interest = gene_data["supplementary_cistron_of_interest"]
				supplementary_reaction_ids_of_interest = gene_data["supplementary_reaction_ids_of_interest"]
				supplementary_reaction_product_of_interest = gene_data["supplementary_reaction_product_of_interest"]

			# Determine gene ids
			with open(simDataFile, 'rb') as f:
				sim_data = pickle.load(f)
			monomer_sim_data = sim_data.process.translation.monomer_data.struct_array

			gene_of_interest_cistron_ids = [cistron_of_interest]
			mRNA_monomer_id_dict = dict(zip(monomer_sim_data['cistron_id'],
											monomer_sim_data['id']))
			gene_of_interest_monomer_ids = [
				mRNA_monomer_id_dict.get(mRNA_id) for mRNA_id in gene_of_interest_cistron_ids]

			# Extract mRNA indexes for each gene of interest
			mRNA_counts_reader = TableReader(os.path.join(simOutDir, 'RNACounts'))
			mRNA_idx_dict = {rna: i for i, rna in enumerate(
				mRNA_counts_reader.readAttribute('mRNA_cistron_ids'))}
			gene_of_interest_cistron_indexes = [
				mRNA_idx_dict.get(mRNA_id) for mRNA_id in gene_of_interest_cistron_ids]
			if "supplementary_cistron_of_interest" in gene_data.keys():
				supplementary_cistron_of_interest_ids = [supplementary_cistron_of_interest]
				supplementary_gene_of_interest_cistron_indexes = [
					mRNA_idx_dict.get(mRNA_id) for mRNA_id in supplementary_cistron_of_interest_ids]
			mRNA_counts_reader.close()

			# Extract protein indexes for each gene of interest
			monomer_counts_reader = TableReader(
				os.path.join(simOutDir, "MonomerCounts"))
			monomer_idx_dict = {
				monomer: i for i, monomer in enumerate(
				monomer_counts_reader.readAttribute('monomerIds'))}
			gene_of_interest_monomer_indexes = [
				monomer_idx_dict.get(monomer_id) for monomer_id in gene_of_interest_monomer_ids]
			if "supplementary_cistron_of_interest" in gene_data.keys():
				supplementary_gene_of_interest_monomer_ids = [
					mRNA_monomer_id_dict.get(mRNA_id) for mRNA_id in supplementary_cistron_of_interest_ids]
				supplementary_gene_of_interest_monomer_indexes = [
					monomer_idx_dict.get(monomer_id) for monomer_id in supplementary_gene_of_interest_monomer_ids]
			monomer_counts_reader.close()

			if 'supplementary_reaction_product_of_interest' in gene_data.keys():
				# List of all reactions that produce or consume the supplementary
				# product of interest
				reactions_for_supplementary_product = []
				# Specifies whether the reaction consumes, produces, or both the
				# supplementary product of interest
				handling_of_supplementary_product = []
				# Involves "ACETYL-COA[c]"
				involves_acetyl_CoA = []
				# Set of all metabolites that contain CoA
				metabolites_containing_CoA = set()
				for reaction in reaction_dict.keys():
					if supplementary_reaction_product_of_interest in reaction_dict[reaction].keys():
						reactions_for_supplementary_product.append(reaction)
						if reaction_dict[reaction][supplementary_reaction_product_of_interest] < 0:
							handling_of_supplementary_product.append("consumes")
						elif reaction_dict[reaction][supplementary_reaction_product_of_interest] > 0:
							handling_of_supplementary_product.append("produces")
						else:
							handling_of_supplementary_product.append("both")

						for key in reaction_dict[reaction].keys():
							if "COA" in key:
								metabolites_containing_CoA.add(key)

						if "ACETYL-COA[c]" in reaction_dict[reaction].keys():
							involves_acetyl_CoA.append(True)
						else:
							involves_acetyl_CoA.append(False)

			# Extract reaction index for each reaction of interest
			fba_results_reader = TableReader(os.path.join(simOutDir, 'FBAResults'))
			base_reaction_ids = fba_results_reader.readAttribute('base_reaction_ids')
			fba_reaction_ids = fba_results_reader.readAttribute('reactionIDs')
			base_reaction_idx_dict = {reaction: i for i, reaction in enumerate(base_reaction_ids)}
			all_reaction_idx_dict = {reaction: i for i, reaction in enumerate(fba_reaction_ids)}
			base_reaction_of_interest_indexes = [
				base_reaction_idx_dict.get(reaction_id) for reaction_id in reaction_ids_of_interest]
			if "supplementary_reaction_ids_of_interest" in gene_data.keys():
				base_supplementary_reaction_of_interest_indexes = [
					base_reaction_idx_dict.get(reaction_id) for reaction_id in supplementary_reaction_ids_of_interest]
			if 'supplementary_reaction_product_of_interest' in gene_data.keys():
				# First, map to base reaction names
				fba_reaction_id_to_base_reaction_id = sim_data.process.metabolism.reaction_id_to_base_reaction_id
				base_reactions_for_supplementary_product = [
					fba_reaction_id_to_base_reaction_id.get(reaction) for reaction
					in reactions_for_supplementary_product]
				base_supplementary_product_reaction_of_interest_indexes = [
					base_reaction_idx_dict.get(reaction_id) for reaction_id in base_reactions_for_supplementary_product]


			# Load data
			time = read_stacked_columns(
				cell_paths, 'Main', 'time', ignore_exception=True)
			time_no_first = read_stacked_columns(
				cell_paths, 'Main', 'time', remove_first=True, ignore_exception=True)

			time = time - time[0]
			time_no_first = time_no_first - time_no_first[0]

			all_monomer_stacked_counts = read_stacked_columns(
				cell_paths, 'MonomerCounts', 'monomerCounts', ignore_exception=True)
			gene_of_interest_monomer_counts = all_monomer_stacked_counts[:,gene_of_interest_monomer_indexes]
			if "supplementary_gene_name" in gene_data.keys():
				supplementary_gene_of_interest_monomer_counts = (
					all_monomer_stacked_counts[:,supplementary_gene_of_interest_monomer_indexes])
			when_is_monomer_nonexistent = time[gene_of_interest_monomer_counts == 0] / 60.

			# Find patches of time where monomer is nonexistent to highlight on plot
			monomer_nonexistent_indexes = np.where(np.array(gene_of_interest_monomer_counts) == 0)[0]
			patch_start_index = []
			patch_end_index = []
			if len(monomer_nonexistent_indexes):
				prev = monomer_nonexistent_indexes[0]
				patch_start_index.append(prev)
				for i in monomer_nonexistent_indexes:
					if np.abs(i - prev) > 1:
						patch_start_index.append(i)
						patch_end_index.append(prev)
					prev = i
				patch_end_index.append(prev)

			all_mRNA_stacked_counts = read_stacked_columns(
				cell_paths, 'RNACounts', 'mRNA_cistron_counts', ignore_exception=True)
			gene_of_interest_mRNA_counts = all_mRNA_stacked_counts[:,gene_of_interest_cistron_indexes]
			if "supplementary_gene_name" in gene_data.keys():
				supplementary_gene_of_interest_mRNA_counts = (
					all_mRNA_stacked_counts[:,supplementary_gene_of_interest_cistron_indexes])

			all_reaction_fluxes = read_stacked_columns(
				cell_paths, 'FBAResults', 'reactionFluxes', ignore_exception=True)
			all_base_reaction_fluxes = read_stacked_columns(
				cell_paths, 'FBAResults', 'base_reaction_fluxes', ignore_exception=True)
			reaction_of_interest_flux = all_base_reaction_fluxes[:,base_reaction_of_interest_indexes]
			if "supplementary_gene_name" in gene_data.keys():
				supplementary_reaction_of_interest_flux = (
					all_base_reaction_fluxes[:,base_supplementary_reaction_of_interest_indexes])
			if 'supplementary_reaction_product_of_interest' in gene_data.keys():
				supplementary_product_reaction_of_interest_flux = (
					all_base_reaction_fluxes[:, base_supplementary_product_reaction_of_interest_indexes])
			'''
			if phase2_reaction_ids_of_interest:
				phase2_reaction_indexes = [
					all_reaction_idx_dict.get(reaction_id) for reaction_id in phase2_reaction_ids_of_interest]
				phase2_reaction_flux = all_reaction_fluxes[:, phase2_reaction_indexes]'''


			(reaction_product_of_interest_counts,) = read_stacked_bulk_molecules(
					cell_paths, [reaction_product_of_interest], ignore_exception=True)
			if "supplementary_gene_name" in gene_data.keys():
				(supplementary_reaction_product_of_interest_counts,) = read_stacked_bulk_molecules(
					cell_paths, [supplementary_reaction_product_of_interest], ignore_exception=True)



			plot_suffixes = [""]
			standard_xlim = (0, 1500)
			total_plots = 18  # TODO Modularize and get rid of this magic number

			mpl.rcParams['axes.spines.right'] = False
			mpl.rcParams['axes.spines.top'] = False

			for i in range(len(plot_suffixes)):

				plot_suffix = plot_suffixes[i]

				# Plotting
				plt.figure(figsize=(10, 22 + 1.5 * total_plots))
				plot_num = 1

				# Mass
				ax1 = plt.subplot(total_plots, 1, plot_num)
				mass = read_stacked_columns(
					cell_paths, "Mass", "cellMass", ignore_exception=True)
				plt.plot(time / 60., mass, color=LINE_COLOR)
				if plot_suffix == "_standard_axes_both" or plot_suffix == "_standard_axes_y":
					plt.ylim(0, 4000)
				if plot_suffix == "_standard_axes_both" or plot_suffix == "_standard_axes_x":
					plt.xlim(standard_xlim)
				plt.xlabel("Time (min)")
				plt.ylabel("Cell Mass (fg)", fontsize="small")
				ax1_patches = []
				for j in range(len(patch_start_index)):
					width = (time[patch_end_index[j]] / 60. - time[patch_start_index[j]] / 60.)[0]
					if width <= 0.1:
						continue
					height = ax1.get_ylim()[1] - ax1.get_ylim()[0]
					ax1_patches.append(
						mpl.patches.Rectangle(
							((time[patch_start_index[j]] / 60.)[0], ax1.get_ylim()[0]),
							width, height, color='gray', alpha=0.15,
							linewidth=0.))
				ax1.add_collection(PatchCollection(ax1_patches, match_original=True))
				plot_num += 1
				# mRNA Cistron Counts
				ax2 = plt.subplot(total_plots, 1, plot_num, sharex=ax1)
				plt.plot(time / 60., gene_of_interest_mRNA_counts, color=LINE_COLOR)
				plt.xlabel("Time (min)")
				plt.ylabel("mRNA Cistron Counts: " + gene_name, fontsize="small")
				ax2_patches = []
				for j in range(len(patch_start_index)):
					width = (time[patch_end_index[j]] / 60. - time[patch_start_index[j]] / 60.)[0]
					if width <= 0.1:
						continue
					height = ax2.get_ylim()[1] - ax2.get_ylim()[0]
					ax2_patches.append(
						mpl.patches.Rectangle(
							((time[patch_start_index[j]] / 60.)[0], ax2.get_ylim()[0]),
							width, height, color='gray', alpha=0.15,
							linewidth=0.))
				ax2.add_collection(PatchCollection(ax2_patches, match_original=True))
				plot_num += 1

				# Protein Counts
				ax3 = plt.subplot(total_plots, 1, plot_num, sharex=ax1)
				plt.plot(time / 60., gene_of_interest_monomer_counts, color=LINE_COLOR)
				plt.xlabel("Time (min)")
				plt.ylabel("Monomer Counts: " + gene_name, fontsize="small")
				ax3_patches = []
				for j in range(len(patch_start_index)):
					width = (time[patch_end_index[j]] / 60. - time[patch_start_index[j]] / 60.)[0]
					if width <= 0.1:
						continue
					height = ax3.get_ylim()[1] - ax3.get_ylim()[0]
					ax3_patches.append(
						mpl.patches.Rectangle(
							((time[patch_start_index[j]] / 60.)[0], ax3.get_ylim()[0]),
							width, height, color='gray', alpha=0.15,
							linewidth=0.))
				ax3.add_collection(PatchCollection(ax3_patches, match_original=True))
				plot_num += 1

				# Reaction flux
				ax4 = plt.subplot(total_plots, 1, plot_num, sharex=ax1)
				for i in range(len(reaction_ids_of_interest)):
					plt.plot(
						time / 60., reaction_of_interest_flux[:, i],
						label=reaction_ids_of_interest[i][0:20] + "...")
				plt.legend(fontsize="x-small")
				plt.xlabel("Time (min)")
				plt.ylabel("Reaction flux (mmol/g DCW/hr)", fontsize="small")
				plt.ylim(-0.0001, 0.003)
				ax4_patches = []
				for j in range(len(patch_start_index)):
					width = (time[patch_end_index[j]] / 60. - time[patch_start_index[j]] / 60.)[0]
					if width <= 0.1:
						continue
					height = ax4.get_ylim()[1] - ax4.get_ylim()[0]
					ax4_patches.append(
						mpl.patches.Rectangle(
							((time[patch_start_index[j]] / 60.)[0], ax4.get_ylim()[0]),
							width, height, color='gray', alpha=0.15,
							linewidth=0.))
				ax4.add_collection(PatchCollection(ax4_patches, match_original=True))
				plot_num += 1

				# Reaction product counts
				ax5 = plt.subplot(total_plots, 1, plot_num, sharex=ax1)
				plt.plot(
					time / 60., reaction_product_of_interest_counts,
					color=LINE_COLOR)
				plt.xlabel("Time (min)")
				plt.ylabel(
					"" + reaction_product_of_interest + " counts",
					fontsize="small")
				ax5_patches = []
				for j in range(len(patch_start_index)):
					width = (time[patch_end_index[j]] / 60. - time[patch_start_index[j]] / 60.)[0]
					if width <= 0.1:
						continue
					height = ax5.get_ylim()[1] - ax5.get_ylim()[0]
					ax5_patches.append(
						mpl.patches.Rectangle(
							((time[patch_start_index[j]] / 60.)[0], ax5.get_ylim()[0]),
							width, height, color='gray', alpha=0.15,
							linewidth=0.))
				ax5.add_collection(PatchCollection(ax5_patches, match_original=True))
				plot_num += 1

				if "supplementary_gene_name" in gene_data.keys():
					# Plot counts of enzyme for next step in pathway

					# mRNA Cistron Counts
					ax2s = plt.subplot(total_plots, 1, plot_num, sharex=ax1)
					plt.plot(time / 60., supplementary_gene_of_interest_mRNA_counts, color=LINE_COLOR)
					plt.xlabel("Time (min)")
					plt.ylabel("mRNA Cistron Counts: " + supplementary_gene_name, fontsize="small")
					ax2s_patches = []
					for j in range(len(patch_start_index)):
						width = (time[patch_end_index[j]] / 60. - time[patch_start_index[j]] / 60.)[
							0]
						if width <= 0.1:
							continue
						height = ax2s.get_ylim()[1] - ax2s.get_ylim()[0]
						ax2s_patches.append(
							mpl.patches.Rectangle(
								((time[patch_start_index[j]] / 60.)[0], ax2s.get_ylim()[0]),
								width, height, color='gray', alpha=0.15,
								linewidth=0.))
					ax2s.add_collection(PatchCollection(ax2s_patches, match_original=True))
					plot_num += 1

					# Protein Counts
					ax3s = plt.subplot(total_plots, 1, plot_num, sharex=ax1)
					plt.plot(time / 60., supplementary_gene_of_interest_monomer_counts, color=LINE_COLOR)
					plt.xlabel("Time (min)")
					plt.ylabel("Monomer Counts: " + supplementary_gene_name, fontsize="small")
					ax3s_patches = []
					for j in range(len(patch_start_index)):
						width = (time[patch_end_index[j]] / 60. - time[patch_start_index[j]] / 60.)[
							0]
						if width <= 0.1:
							continue
						height = ax3s.get_ylim()[1] - ax3s.get_ylim()[0]
						ax3s_patches.append(
							mpl.patches.Rectangle(
								((time[patch_start_index[j]] / 60.)[0], ax3s.get_ylim()[0]),
								width, height, color='gray', alpha=0.15,
								linewidth=0.))
					ax3s.add_collection(PatchCollection(ax3s_patches, match_original=True))
					plot_num += 1

					# Reaction flux
					ax4s = plt.subplot(total_plots, 1, plot_num, sharex=ax1)
					for i in range(len(supplementary_reaction_ids_of_interest)):
						plt.plot(
							time / 60., supplementary_reaction_of_interest_flux[:, i],
							label=supplementary_reaction_ids_of_interest[i][0:20] + "...")
					plt.legend(fontsize="x-small")
					plt.xlabel("Time (min)")
					plt.ylabel("Reaction flux (mmol/g DCW/hr)", fontsize="small")
					plt.ylim(-0.0001, 0.003)
					ax4s_patches = []
					for j in range(len(patch_start_index)):
						width = (time[patch_end_index[j]] / 60. - time[patch_start_index[j]] / 60.)[
							0]
						if width <= 0.1:
							continue
						height = ax4s.get_ylim()[1] - ax4s.get_ylim()[0]
						ax4s_patches.append(
							mpl.patches.Rectangle(
								((time[patch_start_index[j]] / 60.)[0], ax4s.get_ylim()[0]),
								width, height, color='gray', alpha=0.15,
								linewidth=0.))
					ax4s.add_collection(PatchCollection(ax2s_patches, match_original=True))
					plot_num += 1

					# Reaction product counts
					ax5s = plt.subplot(total_plots, 1, plot_num, sharex=ax1)
					plt.plot(
						time / 60., supplementary_reaction_product_of_interest_counts,
						color=LINE_COLOR)
					plt.xlabel("Time (min)")
					plt.ylabel(
						"" + supplementary_reaction_product_of_interest + " counts",
						fontsize="small")
					ax5s_patches = []
					for j in range(len(patch_start_index)):
						width = (time[patch_end_index[j]] / 60. - time[patch_start_index[j]] / 60.)[
							0]
						if width <= 0.1:
							continue
						height = ax5s.get_ylim()[1] - ax5s.get_ylim()[0]
						ax5s_patches.append(
							mpl.patches.Rectangle(
								((time[patch_start_index[j]] / 60.)[0], ax5s.get_ylim()[0]),
								width, height, color='gray', alpha=0.15,
								linewidth=0.))
					ax5s.add_collection(PatchCollection(ax5s_patches, match_original=True))
					plot_num += 1


			plt.subplots_adjust(hspace = 0.7, top = 0.95, bottom = 0.05)
			exportFigure(
					plt, plotOutDir, plotOutFileName + plot_suffix + "_"
					+ str(START_GEN) + "_" + str(END_GEN) + "_" + gene_name,
					metadata)
			plt.close("all")



if __name__ == '__main__':
	Plot().cli()