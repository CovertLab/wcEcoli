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

			# get gist of protein counts that have proteases assigned
			protease_categories, protein_counts = np.unique(protease_assignment, return_counts=True)

			# Compare deg rates across protease assignment
			indices_w_proteases = [i for i, x in enumerate(protease_assignment) if x != 'None']
			assigned_proteins = protein_ids[indices_w_proteases]
			assigned_half_lives = half_lives[indices_w_proteases]
			assigned_proteases = protease_assignment[indices_w_proteases]

			# Map categories to unique integers
			protease_map = {protease: i for i, protease in enumerate(set(assigned_proteases))}
			protease_indices = [protease_map[protease] for protease in assigned_proteases]

			# Create a colormap
			cmap = cm.get_cmap('Set1')
			colors = [cmap(i) for i in protease_indices]


			def extract_dir(path):
				dir_parts = path.split('/')
				try:
					index_of_out = dir_parts.index('out')
					return dir_parts[index_of_out + 1]
				except ValueError:
					return None

			dir_name = extract_dir(input_dir)

			fig = plt.figure(figsize=(12, 4))

			def bar_plot(ax, protease_categories, protein_counts)


			def plot_distributions(ax, half_lives, colors, dir_name):

				# Add jitter to the y-positions
				x_positions = np.random.randn(len(half_lives)) * 0.2

				ax.scatter(half_lives, x_positions, c = colors, alpha= 0.3, s = 5)

				ax.set_title('Distribution of protein half lives')
				ax.set_xlabel('Half life (min)')
				ax.set_ylabel('ParCa: '+ dir_name)
				ax.set_yticks([0])
				ax.set_yticklabels([])
				ax.set_ylim(-2, 2)

			ax1 = plt.subplot(1, 3, 1)
			ax2 = plt.subplot(1, 3, 2)
			plot_distributions(ax2, assigned_half_lives, colors, dir_name)

			# Create a list of patches to represent each category
			handles = [mpatches.Patch(color=cmap(i), label=source) for source, i in source_map.items()]
			fig.legend(handles=handles, bbox_to_anchor=(1, 0.7))

			plt.tight_layout()
			exportFigure(plt, plot_out_dir, plot_out_filename, metadata)
			plt.close('all')


			import ipdb;
			ipdb.set_trace()
			


if __name__ == "__main__":
	Plot().cli()
