"""
Plot histogram of new gene protein mass fraction, colored by variant, for all generations, early generations, and/or late (i.e. not early) generations
"""

import numpy as np
from matplotlib import pyplot as plt

import pickle
import os
from wholecell.io.tablereader import TableReader
from models.ecoli.analysis import variantAnalysisPlot
from wholecell.analysis.analysis_tools import exportFigure, read_stacked_columns
from wholecell.analysis.plotting_tools import DEFAULT_MATPLOTLIB_COLORS as COLORS
from unum.units import g, mol

exclude_timeout_cells = 1 # 1 to exclude cells that took full MAX_CELL_LENGTH, 0 otherwise
exclude_early_gens = 1 # 1 to plot early (before MIN_LATE_CELL_INDEX), and late generationss in addition to all generations

FONT_SIZE=9
MAX_VARIANT = 10 # do not include any variant >= this index
MIN_LATE_CELL_INDEX = 4 # generations before this may not be representative of dynamics due to how they are initialized
MAX_CELL_LENGTH = 180
if exclude_timeout_cells:
	MAX_CELL_LENGTH += 1000

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
		print("Running analysis script with exclude_timeout_cells=", exclude_timeout_cells, " and exclude_early_gens=", exclude_early_gens)

		# Determine new gene indexes
		with open(simDataFile, 'rb') as f:
			sim_data = pickle.load(f)
		mRNA_sim_data = sim_data.process.transcription.cistron_data.struct_array
		monomer_sim_data = sim_data.process.translation.monomer_data.struct_array
		new_gene_mRNA_ids = mRNA_sim_data[mRNA_sim_data['is_new_gene']]['id'].tolist()
		mRNA_monomer_id_dict = dict(zip(monomer_sim_data['cistron_id'], monomer_sim_data['id']))
		new_gene_monomer_ids = [mRNA_monomer_id_dict.get(mRNA_id) for mRNA_id in new_gene_mRNA_ids]
		assert len(new_gene_mRNA_ids) != 0, 'no new gene mRNAs found'
		assert len(new_gene_monomer_ids) != 0, 'no new gene proteins found'
		assert len(new_gene_monomer_ids) == len(new_gene_mRNA_ids), 'number of new gene monomers and mRNAs should be equal'


		# Data extraction
		print("---Data Extraction---")
		doubling_times = {}
		doubling_times_early_gens = {}
		doubling_times_late_gens = {}
		new_gene_monomer_counts = [{} for id in new_gene_monomer_ids]
		new_gene_monomer_mass_fraction = [{} for id in new_gene_monomer_ids]
		if exclude_early_gens:
			new_gene_monomer_counts_early_gens = [{} for id in new_gene_monomer_ids]
			new_gene_monomer_counts_late_gens = [{} for id in new_gene_monomer_ids]

			new_gene_monomer_mass_fraction_early_gens = [{} for id in new_gene_monomer_ids]
			new_gene_monomer_mass_fraction_late_gens = [{} for id in new_gene_monomer_ids]

			new_gene_monomer_masses = [1 for id in new_gene_monomer_ids]
			for i in range(len(new_gene_monomer_ids)):
				new_gene_monomer_masses[i] = float(
					(sim_data.getter.get_mass(new_gene_monomer_ids[i]) / 1000 * 0.0000016605402) / (1 * g / mol))  # convert from g/mol to fg

		variants = self.ap.get_variants()
		min_variant = min(variants)
		for variant in variants:
			if variant >= MAX_VARIANT:
				continue

			print("Variant: ",variant)
			all_cells = self.ap.get_cells(variant=[variant], only_successful=True)
			if len(all_cells) == 0:
				continue

			# Doubling times
			dt = read_stacked_columns(all_cells, 'Main', 'time',
				fun=lambda x: (x[-1] - x[0]) / 60.).squeeze()
			doubling_times[variant] = dt[dt < MAX_CELL_LENGTH]

			if exclude_early_gens:
				all_cells_gens = np.array([int(c.split("/")[-2][-6:]) for c in all_cells])
				included_gens = all_cells_gens[dt < MAX_CELL_LENGTH]
				early_mask = included_gens < MIN_LATE_CELL_INDEX
				late_mask = ~early_mask

				doubling_times_early_gens[variant] = doubling_times[variant][early_mask]
				doubling_times_late_gens[variant] = doubling_times[variant][late_mask]

			# New gene mRNA and monomer counts
			if variant == min_variant:
				sim_dir = all_cells[0]
				simOutDir = os.path.join(sim_dir, 'simOut')

				# Extract protein indexes for each new gene
				monomer_counts_reader = TableReader(os.path.join(simOutDir, "MonomerCounts"))
				monomer_idx_dict = {monomer: i for i, monomer in
							   enumerate(monomer_counts_reader.readAttribute('monomerIds'))}
				new_gene_monomer_indexes = [monomer_idx_dict.get(monomer_id) for monomer_id in new_gene_monomer_ids]

			# Get new gene protein counts
			avg_new_gene_monomer_counts = read_stacked_columns(all_cells, 'MonomerCounts', 'monomerCounts',fun=lambda x: np.mean(x[:,new_gene_monomer_indexes],axis=0))
			avg_new_gene_monomer_counts = avg_new_gene_monomer_counts[dt < MAX_CELL_LENGTH]

			# Get total protein mass
			avg_protein_mass = read_stacked_columns(all_cells, 'Mass', 'proteinMass',fun=lambda x: np.mean(x))

			for i in range(len(new_gene_mRNA_ids)):
				new_gene_monomer_counts[i][variant] = np.log10(avg_new_gene_monomer_counts[:,i] + 1)

				new_gene_monomer_mass_fraction[i][variant] = (avg_new_gene_monomer_counts.flatten() * new_gene_monomer_masses[i] )/avg_protein_mass.flatten()

				if exclude_early_gens:
					new_gene_monomer_counts_early_gens[i][variant] = new_gene_monomer_counts[i][variant][early_mask]
					new_gene_monomer_counts_late_gens[i][variant] = new_gene_monomer_counts[i][variant][late_mask]

					new_gene_monomer_mass_fraction_early_gens[i][variant] = new_gene_monomer_mass_fraction[i][variant][early_mask]
					new_gene_monomer_mass_fraction_late_gens[i][variant] = new_gene_monomer_mass_fraction[i][variant][late_mask]


		# Plotting
		print("---Plotting---")
		std_xlim = [0,1]
		std_sf = 2
		std_bin_width = 0.05

		# ALL GENS
		for i in range(len(new_gene_monomer_ids)):
			std_xlab = '' + new_gene_monomer_ids[i][:-3] +' Mass Fraction'

			_, axes = plt.subplots(1, 1, figsize=(10, 5))
			self.hist(axes, new_gene_monomer_mass_fraction[i], std_xlab, bin_width=std_bin_width, sf=std_sf,xlim=std_xlim)
			plt.tight_layout()
			exportFigure(plt, plotOutDir, plotOutFileName+'_all_gens_'+new_gene_monomer_ids[i][:-3], metadata)

			if exclude_early_gens:
				# EARLY GENS
				_, axes = plt.subplots(1, 1, figsize=(10, 5))
				self.hist(axes, new_gene_monomer_mass_fraction_early_gens[i], std_xlab, bin_width=std_bin_width, sf=std_sf,xlim=std_xlim)
				plt.tight_layout()
				exportFigure(plt, plotOutDir, plotOutFileName + '_early_gens_' + new_gene_monomer_ids[i][:-3], metadata)

				# LATE GENS
				_, axes = plt.subplots(1, 1, figsize=(10, 5))
				self.hist(axes, new_gene_monomer_mass_fraction_late_gens[i], std_xlab, bin_width=std_bin_width, sf=std_sf,xlim=std_xlim)
				plt.tight_layout()
				exportFigure(plt, plotOutDir, plotOutFileName + '_late_gens_' + new_gene_monomer_ids[i][:-3], metadata)

		plt.close("all")


if __name__ == "__main__":
	Plot().cli()
