"""
Plot protein half lives and their distributions across sources of protein degradation rates.
"""

import pickle
import os

from matplotlib import pyplot as plt
from matplotlib import cm
import matplotlib.patches as mpatches

# noinspection PyUnresolvedReferences
import csv
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
			deg_rate_source = sim_data.process.translation.monomer_data['deg_rate_source']
			degradation_rates = sim_data.process.translation.monomer_data['deg_rate'].asNumber(1 / units.s)
			half_lives = np.log(2) / degradation_rates / 60

			# yibQ's protein deg rate is adjusted, half life is set to almost 50000 min
			indices = [i for i, x in enumerate(protein_ids) if x == "EG12298-MONOMER[c]"]

			protein_ids_wo_yibQ = np.delete(protein_ids, indices[0])
			half_lives_wo_yibQ = np.delete(half_lives, indices[0])
			deg_rate_source_wo_yibQ = np.delete(deg_rate_source, indices[0])

			def extract_dir(path):
				dir_parts = path.split('/')
				try:
					index_of_out = dir_parts.index('out')
					return dir_parts[index_of_out + 1]
				except ValueError:
					return None

			dir_name = extract_dir(input_dir)

			fig = plt.figure(figsize=(12, 4))

			# Map categories to unique integers
			source_map = {source: i for i, source in enumerate(set(deg_rate_source_wo_yibQ))}
			source_indices = [source_map[source] for source in deg_rate_source_wo_yibQ]

			# Create a colormap
			cmap = cm.get_cmap('Set1')
			colors = [cmap(i) for i in source_indices]

			def plot_half_lives(ax, protein_ids, half_lives, colors, dir_name):

				sorted_idx = sorted(range(len(half_lives)), key=lambda i: half_lives[i])
				sorted_protein_ids = [protein_ids[i] for i in sorted_idx]
				sorted_half_lives = [half_lives[i] for i in sorted_idx]
				sorted_colors = [colors[i] for i in sorted_idx]

				ax.scatter(sorted_protein_ids, sorted_half_lives, c = sorted_colors, alpha= 0.3, s = 5)
				ax.set_xlabel('Proteins')
				ax.set_xticklabels([])
				ax.set_xticks([])
				ax.set_ylabel('Half lives (min)')
				ax.set_title('Protein half lives in '+ dir_name)

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
			plot_half_lives(ax1, protein_ids_wo_yibQ, half_lives_wo_yibQ, colors, dir_name)
			plot_distributions(ax2, half_lives_wo_yibQ, colors, dir_name)

			# Create a list of patches to represent each category
			handles = [mpatches.Patch(color=cmap(i), label=source) for source, i in source_map.items()]
			fig.legend(handles=handles, bbox_to_anchor=(1, 0.7))

			plt.tight_layout()
			exportFigure(plt, plot_out_dir, plot_out_filename, metadata)
			plt.close('all')


			# Write data to table
			with open(os.path.join(plot_out_dir, plot_out_filename + '.tsv'), 'w') as f:
				writer = csv.writer(f, delimiter='\t')
				writer.writerow([
					'monomer_id', 'degradation_rate(1/s)', 'half_life_(min)', 'degradation_rate_source'
				])

				for i, monomer in enumerate(protein_ids):
					writer.writerow([
						protein_ids[i][:-3], degradation_rates[i], half_lives[i],
						deg_rate_source[i]
					])



if __name__ == "__main__":
	Plot().cli()
