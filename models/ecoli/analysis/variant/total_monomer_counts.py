"""
Plot mRNA and total monomer counts for proteins across multiple generations
and variants.

Use the --variant-range and  --seed-range commands to specify the range of
variants to plot over. By default, all variants will be plotted. In sherlock,
you can specify the specific variants you would like to plot by setting the
experimental_vars variable to a list of the variant numbers you would like to
plot specifically (if not all).

Note: this code works best with simulations containting four or fewer seeds.
"""

import pickle
import os

from matplotlib import pyplot as plt
# noinspection PyUnresolvedReferences
import numpy as np

from models.ecoli.analysis import variantAnalysisPlot
from wholecell.analysis.analysis_tools import (exportFigure,
	read_stacked_bulk_molecules, read_stacked_columns)
from wholecell.io.tablereader import TableReader

""" USER INPUTS"""

# Indicate the proteins to be visualized here (can be multiple):
interest_proteins = np.array([
	#'EG11854-MONOMER[c]',
	'TRYPSYN-APROTEIN[c]',
	'TRYPSYN-BPROTEIN[c]',
	#"NG-GFP-MONOMER[c]",


])

# If you know how many variants you have and want to specify which you would
# like to plot specifically, you can do so here. Otherwise, leave this as 0.
EXPERIMENTAL_VARS = [12, 16, 17, 19]
# example input for selecting specific variants (typically used in sherlock):
# experimental_vars = [16, 17, 18, 19, 20] # don't include the control (0) here

""" END USER INPUTS """

class Plot(variantAnalysisPlot.VariantAnalysisPlot):
	def seed_plot(self, control_var, experimental_vars, monomer_idx):
		"""
		Plot the protein and mRNA counts for the specified proteins across all
		variants and seeds.
		Args:
			control_var: variant with no modifications (control)
			experimental_vars: variants with modifications (such as new gene
			expression turned on) to compare with the control

		Returns: A plot of the total protein counts and mRNA counts for each of
		the specified proteins.
		"""
		# Retrieve the seed data
		seeds = self.ap.get_seeds()

		# Specify color options and line styles for the plot:
		colors = [["deepskyblue", "yellowgreen", "mediumpurple", "deeppink",
				   "darkturquoise", "turquoise","lightcoral", "gold",
				   "darkorange", "indigo", "darkred", "darkgreen", "darkblue",
				   "darkviolet", "cornflowerblue", "darkcyan", "darkmagenta",
				   "darkorange", "green", "blue", "pink", "orange"],]
		LS = ['-', ':', '-.', '--']

		# Generate Protein Counts Plot
		fig, ax = plt.subplots(nrows=2, ncols=1, figsize=(12, 18))

		# iterate through the seeds and plot the data for each:
		for seed in range(len(seeds)):
			# obtain the control variant data and plot it:
			cseed_dir = self.ap.get_cells(
				variant=[control_var], seed=[seed], only_successful=True)
			ctime = read_stacked_columns(
				cseed_dir, 'Main', 'time', ignore_exception=True)
			c_monomer_counts = (read_stacked_columns(
				cseed_dir, 'MonomerCounts', 'monomerCounts',
				ignore_exception=True))[:, monomer_idx]
			print("Plotting control variant seed " + str(seed))

			# Plot the experimental variants:
			if len(self.monomer_id) == 1:
				name = 'seed ' + str(seed) + ' cntrl var '
				ax[0].plot(ctime / 60., c_monomer_counts,
						 linestyle=LS[seed],
						 label=name, color="#FF796C", linewidth=.5)
				var_num = 0
				for variant in experimental_vars:
					print("Plotting variant "+str(variant)+", seed "+str(seed))
					seed_dir = self.ap.get_cells(
						variant=[variant], seed=[seed], only_successful=True)
					# Load data for the seed
					time = read_stacked_columns(seed_dir, 'Main',
						'time', ignore_exception=True)
					monomer_counts = (read_stacked_columns(
						seed_dir, 'MonomerCounts', 'monomerCounts',
						ignore_exception=True))[:, monomer_idx]
					name = 'seed ' + str(seed) + ' exp  var ' + str(variant)
					ax[0].plot(time / 60., monomer_counts, label=name,
						color=colors[0][var_num],linestyle=LS[seed],linewidth=.5)
					var_num = var_num + 1

				# plot specs
				ax[0].set_title(
					f"Total Monomer Counts for {self.monomer_id[0]} in\nthe "
					f"control variant and {len(experimental_vars)} experimental"
					f" variant(s)")

		# format the protein counts plot
		ax[0].set_xlabel("Time (min)")
		ax[0].set_ylabel("Total Monomer Counts")
		ax[0].legend(loc='lower center', bbox_to_anchor=(0.5, -.15),
					 ncols=len(seeds), fontsize='small')

		# Generate the mRNA Counts Plot
		for seed in range(len(seeds)):
			cseed_dir = self.ap.get_cells(
				variant=[control_var], seed=[seed], only_successful=True)
			ctime = read_stacked_columns(
				cseed_dir, 'Main', 'time',  ignore_exception=True)
			c_mRNA_counts = read_stacked_columns(
				cseed_dir, 'RNACounts', 'mRNA_cistron_counts',
				ignore_exception=True)[:, self.mRNA_index]
			if len(self.cistron_id) == 1:
				name = 'seed ' + str(seed) + ' cntrl var '
				ax[1].plot(ctime / 60., c_mRNA_counts, linestyle=LS[seed],
						 label=name, color="#FF796C", linewidth=.5)
				var_num = 0
				for variant in experimental_vars:
					seed_dir = self.ap.get_cells(
						variant=[variant], seed=[seed], only_successful=True)
					# Load data for the seed
					time = read_stacked_columns(seed_dir, 'Main',
						'time', ignore_exception=True)
					mRNA_counts = read_stacked_columns(seed_dir,
						'RNACounts', 'mRNA_cistron_counts',
						ignore_exception=True)[:, self.mRNA_index]
					name = 'seed ' + str(seed) + ' exp  var ' + str(variant)
					ax[1].plot(time / 60., mRNA_counts, label=name,
						color=colors[0][var_num],linestyle=LS[seed],linewidth=.5)
					var_num = var_num + 1

				# plot specs
				ax[1].set_title(f"mRNA Counts for {self.cistron_id[0]} in the "
					f"\ncontrol variant and {len(experimental_vars)} "
					f"experimental variant(s)")

		# format the mRNA counts plot:
		ax[1].set_xlabel("Time (min)")
		ax[1].set_ylabel("Cistron Counts")
		ax[1].legend(loc='lower center', bbox_to_anchor=(0.5, -.15),
					 ncols=len(seeds), fontsize='small')

		# format the final plot:
		plt.tight_layout()
		plt.subplots_adjust(hspace=0.3, top=0.95, bottom=0.15)

	def do_plot(self, simOutDir, plotOutDir, plotOutFileName, simDataFile,
				validationDataFile, metadata):
		with open(simDataFile, 'rb') as f:
			sim_data = pickle.load(f)
		monomer_sim_data = (
			sim_data.process.translation.monomer_data.struct_array)

		# extract info about the protein(s) from the monomer data:
		for protein in interest_proteins:
			monomer_idx = np.where(monomer_sim_data['id'] == protein)
			monomer_idx = monomer_idx[0][0]
			monomer_data = monomer_sim_data[monomer_idx]
			self.cistron_id = [monomer_data['cistron_id']]
			self.monomer_id = [protein]

			# get the simulation output directory:
			cell_paths = self.ap.get_cells()
			sim_dir = cell_paths[0] # arbitrary cell path
			simOutDir = os.path.join(sim_dir, 'simOut')

			# Extract mRNA indexes for the protein (which is not always the
			# same as the monomer index):
			mRNA_counts_reader = TableReader(os.path.join(simOutDir,'RNACounts'))
			mRNA_idx_dict = {rna: i for i, rna in enumerate(
				mRNA_counts_reader.readAttribute('mRNA_cistron_ids'))}
			self.mRNA_index = [
				mRNA_idx_dict.get(mRNA_id) for mRNA_id in self.cistron_id]

			# get the data for each variant:
			all_variants = self.ap.get_variants()
			# specifiy the control variant:
			control_var = all_variants[0]
			# if experimental_vars list is zero, plot all other variants
			if EXPERIMENTAL_VARS == 0:
				experimental_variants = all_variants[1:] # remove the control
			else:
				experimental_variants = EXPERIMENTAL_VARS

			# plot the data:
			print("Plotting for " + str(protein))
			self.seed_plot(control_var, experimental_variants, monomer_idx)

			# export the plot:
			if EXPERIMENTAL_VARS == 0:
				exportFigure(plt, plotOutDir, plotOutFileName +
							 f'_variantPlot_geneID_{str(protein)}', metadata)
			else:
				exportFigure(plt, plotOutDir, plotOutFileName + f'_variantPlot'
					f'_geneID_{str(protein)}_variants_{EXPERIMENTAL_VARS}', metadata)
			plt.close("all")

if __name__ == '__main__':
	Plot().cli()
