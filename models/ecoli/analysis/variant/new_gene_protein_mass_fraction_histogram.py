"""
Plot histogram of new gene protein mass fraction, colored by variant, for all generations, early generations, and late (i.e. not early) generations
"""

import numpy as np
from matplotlib import pyplot as plt

import os
from wholecell.io.tablereader import TableReader
from models.ecoli.analysis import variantAnalysisPlot
from wholecell.analysis.analysis_tools import exportFigure, read_stacked_columns, stacked_cell_identification
from wholecell.analysis.plotting_tools import DEFAULT_MATPLOTLIB_COLORS as COLORS


FONT_SIZE=9
MAX_CELL_LENGTH = 180
#MAX_CELL_LENGTH += 1 # comment out this line to filter sims that reach the max time of 180 min
MIN_LATE_CELL_INDEX = 4 # generations before this may not be representative of dynamics due to how they are initialized


class Plot(variantAnalysisPlot.VariantAnalysisPlot):
	def hist(self, ax, data, xlabel, bin_width=1., xlim=None, sf=1):
		if xlim:
			bins = np.histogram(range(xlim[0],xlim[1]+1), bins=int(np.ceil((xlim[1]-xlim[0])/bin_width)))[1]

		for variant, variant_data in data.items():
			color = COLORS[variant % len(COLORS)]
			if not xlim:
				bins = max(1, int(np.ceil((variant_data.max() - variant_data.min()) / bin_width)))
			mean = variant_data.mean()
			std = variant_data.std()
			ax.hist(variant_data, bins, color=color, alpha=0.5,
				label=f'Var {variant}: {mean:.{sf}f} +/- {std:.{sf+1}f}')
			ax.axvline(mean, color=color, linestyle='--', linewidth=1)

		if xlim:
			ax.set_xlim(xlim)
		self.remove_border(ax)
		ax.set_xlabel(xlabel, fontsize=FONT_SIZE)
		ax.tick_params(labelsize=FONT_SIZE)
		ax.legend()

	def do_plot(self, inputDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		total_protein_mass = {}

		# Data extraction
		print("---Data Extraction---")
		variants = self.ap.get_variants()
		min_variant = min(variants)
		for variant in variants:

			if variant >= 10:
				continue

			print("Variant: ",variant)
			all_cells = self.ap.get_cells(variant=[variant], only_successful=True)
			if len(all_cells) == 0:
				continue

			if len(all_cells) >= MIN_LATE_CELL_INDEX:
				all_cells_gens = [int(c.split("/")[-2][-6:]) for c in all_cells]
				early_cell_index = [i for i, v in enumerate(all_cells_gens) if v < MIN_LATE_CELL_INDEX]
				late_cell_index = [i for i, v in enumerate(all_cells_gens) if v >= MIN_LATE_CELL_INDEX]

			if variant == min_variant: ### TODO flag new gene mRNAs and proteins more efficiently
				sim_dir = all_cells[0]
				simOutDir = os.path.join(sim_dir, 'simOut')

				# Extract protein indexes for each new gene
				monomer_counts_reader = TableReader(os.path.join(simOutDir, "MonomerCounts"))
				monomer_idx = {monomer: i for i, monomer in
							   enumerate(monomer_counts_reader.readAttribute('monomerIds'))}
				new_gene_monomer_ids = [k for k, v in monomer_idx.items() if k.startswith('NG')]
				new_gene_monomer_indexes = [v for k, v in monomer_idx.items() if k.startswith('NG')]
				assert len(new_gene_monomer_ids) != 0, 'no new gene proteins found'

				new_gene_monomer_counts = [{} for id in new_gene_monomer_ids]
				new_gene_monomer_counts_early_gens = [{} for id in new_gene_monomer_ids]
				new_gene_monomer_counts_late_gens = [{} for id in new_gene_monomer_ids]
				new_gene_monomer_mass_fraction = [{} for id in new_gene_monomer_ids]
				new_gene_monomer_mass_fraction_early_gens = [{} for id in new_gene_monomer_ids]
				new_gene_monomer_mass_fraction_late_gens = [{} for id in new_gene_monomer_ids]

			all_monomer_counts = read_stacked_columns(all_cells, 'MonomerCounts', 'monomerCounts')

			cell_id_vector = stacked_cell_identification(all_cells, 'Main', 'time')
			cell_ids, idx, cell_total_timesteps = np.unique(cell_id_vector, return_inverse=True, return_counts=True)

			# Get total protein mass
			protein_mass_stacked = read_stacked_columns(all_cells, 'Mass','proteinMass')
			protein_mass_stacked = protein_mass_stacked[:,0]
			sum_protein_mass = np.bincount(idx, weights=protein_mass_stacked)
			avg_protein_mass = sum_protein_mass / cell_total_timesteps
			total_protein_mass[variant] = avg_protein_mass

			# TODO: Get molecular weight of new gene proteins
			molecular_weight = 27 * 0.0000016605402  # hardcode mw of GFP for now, 27 kD -> convert to fg


			# Get new gene protein counts
			for i in range(len(new_gene_monomer_ids)):
				new_gene_monomer_counts_var = all_monomer_counts[:, new_gene_monomer_indexes[i]]

				sum_new_gene_monomer_counts = np.bincount(idx, weights=new_gene_monomer_counts_var)
				avg_new_gene_monomer_counts = sum_new_gene_monomer_counts / cell_total_timesteps

				new_gene_monomer_counts[i][variant] = np.log10(avg_new_gene_monomer_counts + 1)
				new_gene_monomer_mass_fraction[i][variant] = (avg_new_gene_monomer_counts * molecular_weight)/avg_protein_mass

				if len(all_cells) >= MIN_LATE_CELL_INDEX:
					new_gene_monomer_counts_early_gens[i][variant] = new_gene_monomer_counts[i][variant][early_cell_index]
					new_gene_monomer_counts_late_gens[i][variant] = new_gene_monomer_counts[i][variant][late_cell_index]

					new_gene_monomer_mass_fraction_early_gens[i][variant] = new_gene_monomer_mass_fraction[i][variant][
						early_cell_index]
					new_gene_monomer_mass_fraction_late_gens[i][variant] = new_gene_monomer_mass_fraction[i][variant][late_cell_index]

		# Plotting
		print("---Plotting---")
		# ALL GENS
		for i in range(len(new_gene_monomer_ids)):
			_, axes = plt.subplots(2, 1, figsize=(10, 10))
			self.hist(axes[0], new_gene_monomer_mass_fraction[i], '' + new_gene_monomer_ids[i][:-3] +' Mass Fraction', bin_width=0.05, sf=2,xlim=[0,1])
			plt.tight_layout()
			exportFigure(plt, plotOutDir, plotOutFileName+'_all_gens_'+new_gene_monomer_ids[i][:-3], metadata)

		if len(all_cells) >= MIN_LATE_CELL_INDEX:
			# EARLY GENS
			for i in range(len(new_gene_monomer_ids)):
				_, axes = plt.subplots(2, 1, figsize=(10, 10))
				self.hist(axes[0], new_gene_monomer_mass_fraction_early_gens[i], '' + new_gene_monomer_ids[i][:-3] +' Mass Fraction', bin_width=0.05, sf=2,xlim=[0,1])
				plt.tight_layout()
				exportFigure(plt, plotOutDir, plotOutFileName + '_early_gens_' + new_gene_monomer_ids[i][:-3], metadata)

			# LATE GENS
			for i in range(len(new_gene_monomer_ids)):
				_, axes = plt.subplots(2, 1, figsize=(10, 10))
				self.hist(axes[0], new_gene_monomer_mass_fraction_late_gens[i],'' + new_gene_monomer_ids[i][:-3] +' Mass Fraction', bin_width=0.05, sf=2,xlim=[0,1])
				plt.tight_layout()
				exportFigure(plt, plotOutDir, plotOutFileName + '_late_gens_' + new_gene_monomer_ids[i][:-3], metadata)

		plt.close("all")


if __name__ == "__main__":
	Plot().cli()
