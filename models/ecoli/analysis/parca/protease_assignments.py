"""
Plot protein half lives and their distributions across sources of protein degradation rates.
"""

import pickle
import os

from matplotlib import pyplot as plt
from matplotlib import cm
import matplotlib.patches as mpatches

# noinspection PyUnresolvedReferences
import numpy as np

from models.ecoli.analysis import parcaAnalysisPlot
from wholecell.analysis.analysis_tools import exportFigure
from wholecell.utils import constants
from wholecell.utils import units


class Plot(parcaAnalysisPlot.ParcaAnalysisPlot):
	def do_plot(self, input_dir, plot_out_dir, plot_out_filename, sim_data_file, validation_data_file, metadata):

		with open(sim_data_file, 'rb') as f:
			sim_data = pickle.load(f)

			protein_ids = sim_data.process.translation.monomer_data['id']
			protease_assignment = sim_data.process.translation.monomer_data['protease_assignment']
			ClpP_contribution = sim_data.process.translation.monomer_data['ClpP_fraction']
			Lon_contribution = sim_data.process.translation.monomer_data['Lon_fraction']
			HslV_contribution = sim_data.process.translation.monomer_data['HslV_fraction']
			Unexplained_contribution = sim_data.process.translation.monomer_data['Unexplained_fraction']

			deg_rate_source = sim_data.process.translation.monomer_data['deg_rate_source']
			degradation_rates = sim_data.process.translation.monomer_data['deg_rate'].asNumber(1 / units.s)
			half_lives = np.log(2) / degradation_rates / 60

			# Compare deg rates across protease assignment
			indices_w_proteases = [i for i, x in enumerate(protease_assignment) if x != 'None']

			assigned_proteases = protease_assignment[indices_w_proteases]

			# get gist of protein counts that have proteases assigned
			protease_categories, protein_counts = np.unique(assigned_proteases, return_counts=True)

			# Map categories to unique integers
			protease_map = {protease: i for i, protease in enumerate(set(assigned_proteases))}

			# Create a colormap
			cmap = cm.get_cmap('tab20')

			protease_data = {}
			for protease in np.unique(assigned_proteases):
				indices_w_protease = [i for i, x in enumerate(protease_assignment) if x == protease]
				assigned_proteins_p = protein_ids[indices_w_protease]
				assigned_half_lives_p = half_lives[indices_w_protease]
				assigned_protease_idx = protease_assignment[indices_w_protease]
				protease_color_idx = [protease_map[protease] for protease in assigned_protease_idx]
				i = protease_map[protease]
				colors = [cmap(i) for i in protease_color_idx]
				protease_data[protease] = [assigned_proteins_p, assigned_half_lives_p, colors]


			def extract_dir(path):
				dir_parts = path.split('/')
				try:
					index_of_out = dir_parts.index('out')
					return dir_parts[index_of_out + 1]
				except ValueError:
					return None

			dir_name = extract_dir(input_dir)

			fig = plt.figure(figsize=(20, 10))

			def bar_plot(ax, protease_categories, protein_counts):
				# Create a bar plot
				ax.bar(protease_categories, protein_counts)

				# Add labels and title
				ax.set_xlabel('Protease assignments')
				ax.tick_params(axis='x', labelrotation=90)
				ax.set_ylabel('# of proteins')
				ax.set_title('Proteins with protease assignments')

			def plot_distributions(ax, protease_data):

				# Calculate positions for each category
				def calculate_positions(x, y, jitter=0.2):
					n = len(y)
					positions = np.linspace(x - jitter, x + jitter, n)
					np.random.shuffle(positions)
					return positions

				for i, category in enumerate(protease_data):
					half_lives = protease_data[category][1]
					y_positions = calculate_positions(i, half_lives)
					ax.scatter(half_lives, y_positions, c=protease_data[category][2], alpha=0.6, s=5)

				ax.set_xlabel('Half life (min)')
				ax.set_ylabel('Protease assignment')
				ax.set_yticks(range(len(protease_data)), protease_data.keys())
				ax.set_title('Distribution of protein half lives')


			ax1 = plt.subplot(1, 3, 1)
			ax2 = plt.subplot(1, 3, 2)

			bar_plot(ax1, protease_categories, protein_counts)
			plot_distributions(ax2, protease_data)

			# Create a list of patches to represent each category
			handles = [mpatches.Patch(color=cmap(i), label=protease) for protease, i in protease_map.items()]
			fig.legend(handles=handles, bbox_to_anchor=(1, 0.7))

			plt.tight_layout()
			exportFigure(plt, plot_out_dir, plot_out_filename, metadata)
			plt.close('all')


			import ipdb;
			ipdb.set_trace()



if __name__ == "__main__":
	Plot().cli()
