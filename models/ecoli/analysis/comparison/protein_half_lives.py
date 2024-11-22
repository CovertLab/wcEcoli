"""
Template for comparison analysis plots
"""

from typing import Tuple

from matplotlib import pyplot as plt
# noinspection PyUnresolvedReferences
import seaborn as sns
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

		protein_ids_1 = sim_data1.process.translation.monomer_data['id']
		protein_ids_2 = sim_data2.process.translation.monomer_data['id']
		degradation_rates_1 = sim_data1.process.translation.monomer_data['deg_rate'].asNumber(1 / units.s)
		degradation_rates_2 = sim_data2.process.translation.monomer_data['deg_rate'].asNumber(1 / units.s)
		half_lives_1 = np.log(2)/degradation_rates_1/60 # in minutes
		half_lives_2 = np.log(2)/degradation_rates_2/60 # in minutes


		def extract_dir(path):
			dir_parts = path.split('/')
			try:
				index_of_out = dir_parts.index('out')
				return dir_parts[index_of_out + 1]
			except ValueError:
				return None

		ap1_dir_name = extract_dir(ap1._path_data['path'][0])
		ap2_dir_name = extract_dir(ap2._path_data['path'][0])

		categories = {ap1_dir_name: dict(zip(protein_ids_1, half_lives_1)),
					  ap2_dir_name: dict(zip(protein_ids_2, half_lives_2))}


		# Calculate positions for each category
		def calculate_positions(x, y, jitter=0.2):
			n = len(y)
			positions = np.linspace(x - jitter, x + jitter, n)
			np.random.shuffle(positions)
			return positions

		plt.figure(figsize=(8, 6))


		for i,category in enumerate(categories):
			values = categories[category].values()
			x_positions = calculate_positions(i, values)
			plt.scatter(values, x_positions, label=category)


		plt.ylabel('ParCa')
		plt.xlabel('Half life (min)')
		plt.title('Distribution of protein half lives')
		plt.yticks(range(len(categories)), categories.keys())
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
