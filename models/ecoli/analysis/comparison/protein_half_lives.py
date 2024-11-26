"""
Template for comparison analysis plots
"""

from typing import Tuple

from matplotlib import pyplot as plt
from matplotlib import cm
import matplotlib.patches as mpatches
# noinspection PyUnresolvedReferences
import numpy as np

from models.ecoli.analysis import comparisonAnalysisPlot
from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from reconstruction.ecoli.simulation_data import SimulationDataEcoli
from validation.ecoli.validation_data import ValidationDataEcoli
from wholecell.analysis.analysis_tools import exportFigure
from wholecell.utils import units
# noinspection PyUnresolvedReferences
from wholecell.io.tablereader import TableReader


class Plot(comparisonAnalysisPlot.ComparisonAnalysisPlot):
	def do_plot(self, reference_sim_dir, plotOutDir, plotOutFileName, input_sim_dir, unused, metadata):
		# manual/analysisComparison.py can compare any two sim dirs.

		# noinspection PyUnusedLocal
		ap1, sim_data1, validation_data1 = self.setup(reference_sim_dir)
		# noinspection PyUnusedLocal
		ap2, sim_data2, validation_data2 = self.setup(input_sim_dir)

		if ap1.n_generation <= 2 or ap2.n_generation <= 2:
			print('Skipping analysis -- not enough sims run.')
			return

		def get_protein_data(sim_data, remove_yibQ):
			"""remove_yibQ is a boolean"""
			protein_ids= sim_data.process.translation.monomer_data['id']
			deg_rate_source = sim_data.process.translation.monomer_data['deg_rate_source']
			degradation_rates = sim_data.process.translation.monomer_data['deg_rate'].asNumber(1 / units.s)
			half_lives = np.log(2) / degradation_rates / 60  # in minutes
			if remove_yibQ:
				indices = [i for i, x in enumerate(protein_ids) if x == "EG12298-MONOMER[c]"]
				protein_ids_wo_yibQ = np.delete(protein_ids, indices[0])
				half_lives_wo_yibQ = np.delete(half_lives, indices[0])
				deg_rate_source_wo_yibQ = np.delete(deg_rate_source, indices[0])
				return protein_ids_wo_yibQ, half_lives_wo_yibQ, deg_rate_source_wo_yibQ
			else:
				return protein_ids, half_lives, deg_rate_source

		def extract_dir(path):
			dir_parts = path.split('/')
			try:
				index_of_out = dir_parts.index('out')
				return dir_parts[index_of_out + 1]
			except ValueError:
				return None

		ap1_dir_name = extract_dir(ap1._path_data['path'][0])
		ap2_dir_name = extract_dir(ap2._path_data['path'][0])
		protein_ids_1, half_lives_1, deg_rate_source_1 = get_protein_data(sim_data1, True)
		protein_ids_2, half_lives_2, deg_rate_source_2 = get_protein_data(sim_data2, True)

		# Create a colormap
		cmap = cm.get_cmap('Set1')

		unique_source_1 = set(deg_rate_source_1)
		unique_source_2 = set(deg_rate_source_2)
		deg_sources_combined = unique_source_1.union(unique_source_2)
		source_map = {source: i for i, source in enumerate(deg_sources_combined)}

		def color_map(deg_rate_source, source_map):
			source_indices = [source_map[source] for source in deg_rate_source]
			colors = [cmap(i) for i in source_indices]
			return colors

		protein_data = {ap1_dir_name:[protein_ids_1, half_lives_1, color_map(deg_rate_source_1, source_map)],
						ap2_dir_name: [protein_ids_2, half_lives_2, color_map(deg_rate_source_2, source_map)]}

		fig = plt.figure(figsize=(12, 4))

		def plot_distributions(ax, protein_data):

			# Calculate positions for each category
			def calculate_positions(x, y, jitter=0.2):
				n = len(y)
				positions = np.linspace(x - jitter, x + jitter, n)
				np.random.shuffle(positions)
				return positions

			for i,category in enumerate(protein_data):
				half_lives = protein_data[category][1]
				y_positions = calculate_positions(i, half_lives)
				ax.scatter(half_lives, y_positions, c = protein_data[category][2], alpha= 0.3, s = 5)

			ax.set_xlabel('Half life (min)')
			ax.set_ylabel('ParCa')
			ax.set_yticks(range(len(protein_data)), protein_data.keys())
			ax.set_title('Distribution of protein half lives')

		def plot_half_lives(ax, protein_data):
			category_names = list(protein_data.keys())
			protein_ids_1 = protein_data[category_names[0]][0]
			protein_ids_2 = protein_data[category_names[1]][0]
			category_1_set = set(protein_ids_1)
			category_2_set = set(protein_ids_2)

			common_proteins = {protein for protein in category_1_set} & {protein for protein in category_2_set}
			category_1_x_values = []
			category_2_y_values = []
			for protein in common_proteins:
				index_1 = [i for i, x in enumerate(protein_ids_1) if x == protein][0]
				index_2 = [i for i, x in enumerate(protein_ids_2) if x == protein][0]
				x_value = protein_data[category_names[0]][1][index_1]
				y_value = protein_data[category_names[1]][1][index_2]
				category_1_x_values.append(x_value)
				category_2_y_values.append(y_value)

			ax.scatter(category_1_x_values , category_2_y_values, color = '#555555', alpha= 0.3, s = 5)
			ax.set_xlabel(category_names[0] + ' protein half lives (min)')
			ax.set_ylabel(category_names[1] + ' protein half lives (min)')
			ax.set_title('Protein half lives')

		ax1 = plt.subplot(1, 3, 1)
		ax2 = plt.subplot(1, 3, 2)

		plot_distributions(ax1, protein_data)
		plot_half_lives(ax2, protein_data)

		# Create a list of patches to represent each category
		handles = [mpatches.Patch(color=cmap(i), label=source) for source, i in source_map.items()]
		fig.legend(handles=handles, bbox_to_anchor=(1, 0.7))

		plt.tight_layout()
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close('all')


	def setup(self, inputDir: str) -> Tuple[
				AnalysisPaths, SimulationDataEcoli, ValidationDataEcoli]:
			"""Return objects used for analyzing multiple sims."""
			ap = AnalysisPaths(inputDir, variant_plot=True)
			sim_data = self.read_sim_data_file(inputDir)
			validation_data = self.read_validation_data_file(inputDir)
			return ap, sim_data, validation_data

if __name__ == "__main__":
	Plot().cli()
