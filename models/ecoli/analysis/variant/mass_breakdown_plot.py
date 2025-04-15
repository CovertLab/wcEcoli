"""
Analyze the average mass breakdowns for each variant index as
bar charts.
"""

import pickle

import numpy as np
import matplotlib.pyplot as plt

from models.ecoli.analysis import variantAnalysisPlot
from wholecell.analysis.analysis_tools import exportFigure, read_stacked_columns

# Remove first N gens from plot
IGNORE_FIRST_N_GENS = 16

poster_colors = {
    "light_gray": (0.75, 0.75, 0.75),
    "poster_green": (66/255, 170/255, 154/255),
    "poster_blue": (27/255, 132/255, 198/255),
    "poster_purple": (188/255, 140/255, 191/255),
    "poster_gold": (221/255, 203/255, 119/255),
    "poster_light_blue": (136/255, 205/255, 240/255),
    "poster_red": (202/255, 0/255, 32/255),
}

class Plot(variantAnalysisPlot.VariantAnalysisPlot):
	def do_plot(self, inputDir, plotOutDir, plotOutFileName, simDataFile,
				validationDataFile, metadata):

		with open(simDataFile, 'rb') as f:
			sim_data = pickle.load(f)
		with open(validationDataFile, 'rb') as f:
			validation_data = pickle.load(f)

		n_total_gens = self.ap.n_generation

		selected_variant_indexes = self.ap.get_variants()

		n_variants = len(selected_variant_indexes)

		# Initialize everything
		all_cells = self.ap.get_cells(
			variant=[selected_variant_indexes[0]],
			generation=np.arange(IGNORE_FIRST_N_GENS, n_total_gens),
			only_successful=True)

		avg_mass = np.zeros((n_variants, 1))
		avg_DNA_mass = np.zeros((n_variants, 1))
		avg_mRNA_mass = np.zeros((n_variants, 1))
		avg_protein_mass = np.zeros((n_variants, 1))
		avg_rRNA_mass = np.zeros((n_variants, 1))
		avg_tRNA_mass = np.zeros((n_variants, 1))
		avg_membrane_mass = np.zeros((n_variants, 1))
		avg_water_mass = np.zeros((n_variants, 1))
		avg_dry_mass = np.zeros((n_variants, 1))
		avg_small_molecule_mass = np.zeros((n_variants, 1))

		# Loop through variant indexes
		for i, variant_index in enumerate(selected_variant_indexes):
			# Get all cells (within the generation range) of this variant index
			all_cells = self.ap.get_cells(
				variant=[variant_index],
				generation=np.arange(IGNORE_FIRST_N_GENS, n_total_gens),
				only_successful=True)

			if len(all_cells) == 0:
				continue

			# Load averages for all mass categories
			avg_mass[i] = np.mean(read_stacked_columns(
				all_cells, 'Mass', 'cellMass'
			))
			avg_DNA_mass[i] = np.mean(read_stacked_columns(
				all_cells, 'Mass', 'dnaMass'
			))
			avg_mRNA_mass[i] = np.mean(read_stacked_columns(
				all_cells, 'Mass', 'mRnaMass'
			))
			avg_protein_mass[i] = np.mean(read_stacked_columns(
				all_cells, 'Mass', 'proteinMass'
			))
			avg_rRNA_mass[i] = np.mean(read_stacked_columns(
				all_cells, 'Mass', 'rRnaMass'
			))
			avg_tRNA_mass[i] = np.mean(read_stacked_columns(
				all_cells, 'Mass', 'tRnaMass'
			))
			avg_membrane_mass[i] = np.mean(read_stacked_columns(
				all_cells, 'Mass', 'membrane_mass'
			))
			avg_water_mass[i] = np.mean(read_stacked_columns(
				all_cells, 'Mass', 'waterMass'
			))
			avg_dry_mass[i] = np.mean(read_stacked_columns(
				all_cells, 'Mass', 'dryMass'
			))
			avg_small_molecule_mass[i] = np.mean(read_stacked_columns(
				all_cells, 'Mass', 'smallMoleculeMass'
			))

		# Create a list of the mass categories
		mass_categories = ['Dry Mass', 'DNA Mass', 'mRNA Mass', 'Protein Mass',
						  'rRNA Mass', 'tRNA Mass', 'Membrane Mass', 'Small Molecule Mass']
		# Create a list of the average masses for each category
		avg_masses = [avg_dry_mass, avg_DNA_mass, avg_mRNA_mass, avg_protein_mass,
					  avg_rRNA_mass, avg_tRNA_mass, avg_membrane_mass, avg_small_molecule_mass]
		# Create a list of the colors for each category
		mass_colors = [poster_colors['poster_blue'], poster_colors['poster_green'],
					   poster_colors['poster_purple'], poster_colors['poster_gold'],
					   poster_colors['poster_light_blue'], poster_colors['poster_red'],
					   poster_colors['light_gray'], 'magenta']
		mass_colors_dict = {}
		for i, mass_category in enumerate(mass_categories):
			mass_colors_dict[mass_category] = mass_colors[i]

		# Sort the variants by descending cell mass
		sorted_indexes = np.argsort(avg_dry_mass[:, 0])[::-1]
		sorted_avg_dry_mass = avg_dry_mass[sorted_indexes]
		sorted_avg_DNA_mass = avg_DNA_mass[sorted_indexes]
		sorted_avg_mRNA_mass = avg_mRNA_mass[sorted_indexes]
		sorted_avg_protein_mass = avg_protein_mass[sorted_indexes]
		sorted_avg_rRNA_mass = avg_rRNA_mass[sorted_indexes]
		sorted_avg_tRNA_mass = avg_tRNA_mass[sorted_indexes]
		sorted_avg_membrane_mass = avg_membrane_mass[sorted_indexes]
		sorted_avg_small_molecule_mass = avg_small_molecule_mass[sorted_indexes]
		sorted_selected_variant_indexes = np.array(selected_variant_indexes)[sorted_indexes]
		# Create a list of the average masses for each category
		sorted_avg_masses = [sorted_avg_dry_mass, sorted_avg_DNA_mass, sorted_avg_mRNA_mass,
							 sorted_avg_protein_mass, sorted_avg_rRNA_mass, sorted_avg_tRNA_mass,
							 sorted_avg_membrane_mass, sorted_avg_small_molecule_mass]
		# Create a list of the colors for each category
		sorted_avg_masses_by_category = {}
		for i, mass_category in enumerate(mass_categories):
			sorted_avg_masses_by_category[mass_category] = sorted_avg_masses[i]

		x = np.arange(len(sorted_selected_variant_indexes)) * 20
		width = 1.5
		multiplier = 0
		fig, ax = plt.subplots(figsize=(14, 6))
		for i, (mass_category, avg_mass_for_category) in enumerate(
				sorted_avg_masses_by_category.items()):
			offset = width * multiplier
			bars = ax.bar(
				x + offset, avg_mass_for_category.squeeze(), width,
				label=mass_category, color=mass_colors_dict[mass_category])
			multiplier += 1
		ax.set_ylabel('Mass (fg)')
		ax.set_title('Mass Breakdown by Variant', fontsize=14, weight='bold')
		ax.set_xticks(x)
		ax.set_xticklabels(sorted_selected_variant_indexes, rotation=45, ha='right')
		ax.set_xlabel('Variant Index')
		ax.legend(loc='upper right', frameon=False, ncol=3)

		plt.tight_layout()
		exportFigure(plt, plotOutDir, f"{plotOutFileName}_dry_sorted", metadata)
		plt.close('all')

		# Make a bar plot where instead of plotting mass per category we plot
		# the mass per category as a percentage of the total mass
		sorted_avg_masses_by_category_as_percentage = {}
		for i, mass_category in enumerate(mass_categories):
			sorted_avg_masses_by_category_as_percentage[mass_category] = (
				sorted_avg_masses[i] / sorted_avg_dry_mass) * 100
		x = np.arange(len(sorted_selected_variant_indexes)) * 20
		width = 1.5
		multiplier = 0
		fig, ax = plt.subplots(figsize=(14, 6))
		for mass_category, avg_mass_for_category in sorted_avg_masses_by_category_as_percentage.items():
			if mass_category == 'Dry Mass' or mass_category == 'Water Mass':
				continue
			offset = width * multiplier
			bars = ax.bar(
				x + offset, avg_mass_for_category.squeeze(), width, label=mass_category,
				color=mass_colors_dict[mass_category])
			multiplier += 1
		ax.set_ylabel('Mass (%)')
		ax.set_title('Mass Breakdown by Variant as a Percentage of Total Dry Mass', fontsize=14, weight='bold')
		ax.set_xticks(x)
		ax.set_xticklabels(sorted_selected_variant_indexes, rotation=45, ha='right')
		ax.set_xlabel('Variant Index')
		ax.legend(loc='upper left', frameon=False, ncol=3)
		exportFigure(plt, plotOutDir, f"{plotOutFileName}_percentage_dry_sorted", metadata)
		plt.close('all')

		# Make a separate plot for each mass category
		for mass_category, avg_mass_for_category in sorted_avg_masses_by_category.items():
			fig, ax = plt.subplots()
			bars = ax.bar(
				x, avg_mass_for_category.squeeze(), width, label=mass_category,
				color=mass_colors_dict[mass_category])
			ax.set_ylabel('Mass (fg)')
			ax.set_title(f'{mass_category} by Variant')
			ax.set_xticks(x, sorted_selected_variant_indexes)
			ax.set_xticklabels(sorted_selected_variant_indexes)
			ax.legend()
			ax.set_xlabel('Variant Index')
			exportFigure(plt, plotOutDir, f"{plotOutFileName}_{mass_category}_sorted", metadata)
			plt.close('all')

		# Make a separate mass % plot for each mass category
		for mass_category, avg_mass_for_category in sorted_avg_masses_by_category_as_percentage.items():
			fig, ax = plt.subplots()
			bars = ax.bar(
				x, avg_mass_for_category.squeeze(), width, label=mass_category,
				color=mass_colors_dict[mass_category])
			ax.set_ylabel('Mass (%)')
			ax.set_title(f'{mass_category} by Variant as Percentage of Dry Mass')
			ax.set_xticks(x, sorted_selected_variant_indexes)
			ax.set_xticklabels(sorted_selected_variant_indexes)
			ax.legend()
			ax.set_xlabel('Variant Index')
			exportFigure(plt, plotOutDir, f"{plotOutFileName}_{mass_category}_percentage_dry_sorted", metadata)
			plt.close('all')


if __name__ == "__main__":
	Plot().cli()
