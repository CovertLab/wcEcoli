"""
Plots proteomics/fluxomics validation plots for two sets of simulations.
"""

from typing import Tuple

from matplotlib import pyplot as plt
# noinspection PyUnresolvedReferences
import numpy as np
import os
from scipy.stats import pearsonr
import pickle
from models.ecoli.analysis import comparisonAnalysisPlot
from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from models.ecoli.processes.metabolism import (
	COUNTS_UNITS, VOLUME_UNITS, TIME_UNITS, MASS_UNITS)
from reconstruction.ecoli.simulation_data import SimulationDataEcoli
from validation.ecoli.validation_data import ValidationDataEcoli
from wholecell.analysis.analysis_tools import exportFigure, read_stacked_columns
# noinspection PyUnresolvedReferences
from wholecell.analysis.analysis_tools import (exportFigure,
	read_bulk_molecule_counts, read_stacked_bulk_molecules, read_stacked_columns)
from wholecell.io.tablereader import TableReader
from wholecell.utils.protein_counts import get_simulated_validation_counts
from wholecell.utils import units

''' USER INPUTS '''
IGNORE_FIRST_N_GENS = 2
# note: the reference_sim is listed second when the function is called, the input_sim is listed first


class Plot(comparisonAnalysisPlot.ComparisonAnalysisPlot):
	def setup(self, inputDir: str) -> Tuple[
			AnalysisPaths, SimulationDataEcoli, ValidationDataEcoli]:
		"""Return objects used for analyzing multiple sims."""
		ap = AnalysisPaths(inputDir, variant_plot=True)
		sim_data = self.read_sim_data_file(inputDir)
		validation_data = self.read_validation_data_file(inputDir)
		return ap, sim_data, validation_data
	def generate_data(self, sim_dir):
		"""
		# function from analyze_complex_counts.py
        Generates csv files of protein count data from simulations
        Args:
            simDataFile: simulation data file
        Returns:
            #TODO: update this description
            protein_counts: protein count (PC) data for all proteins (originally
             present on the E.coli chromosome) in the simulation for each variant
              (the PC for each protein is averaged over all the generations)
            self.total_protein_counts: the original PCs and new gene (NG) PCs
            in one variable
            self.new_gene_monomer_ids: protein ids for new genes inserted into
            the E.coli genome
            self.original_gene_ids: protein ids for the original proteins on
            the E.coli genome
            self.all_monomer_ids: list of all the monomer ids (NG protein ids
            and orginal proteins' gene ids)
        """
		# Read data from sim_data
		ap, sim_data, validation_data = self.setup(sim_dir)

		monomer_sim_data = (
			sim_data.process.translation.monomer_data.struct_array)
		self.all_monomer_ids = monomer_sim_data['id']

		# Get the paths for all cells:
		n_total_gens = ap.n_generation
		# todo: change this depending on which sim you are using (dont always have to have a variant)
		all_cells = ap.get_cells(
			generation=np.arange(IGNORE_FIRST_N_GENS, n_total_gens),
			only_successful=True)

		# Get the average total protein counts for each monomer:
		total_counts = (
			read_stacked_columns(all_cells, 'MonomerCounts',
								 'monomerCounts', ignore_exception=True))
		avg_total_counts = np.mean(total_counts, axis=0)

		# Get the average free protein counts for each monomer:
		(free_counts,) = read_stacked_bulk_molecules(
			all_cells, self.all_monomer_ids, ignore_exception=True)
		avg_free_counts = np.mean(free_counts, axis=0)

		# Get the average complex counts for each monomer:
		avg_counts_for_monomers_in_complexs = avg_total_counts - avg_free_counts

		# get the counts data for specific complexes:
		complex_sim_data = (
			sim_data.process.complexation.ids_complexes)

		self.all_complex_ids = np.array(complex_sim_data)

		# get the counts of the complexes:
		(complex_counts,) = read_stacked_bulk_molecules(
			all_cells, self.all_complex_ids, ignore_exception=True)
		avg_complex_counts = np.mean(complex_counts, axis=0)



		return avg_total_counts, avg_free_counts, avg_counts_for_monomers_in_complexs, avg_complex_counts


	def do_plot(self, reference_sim_dir, plotOutDir, plotOutFileName, input_sim_dir, unused, metadata):
		# From fw_queue, reference_sim_dir has operons="off"; input_sim_dir has
		# operons="on".
		# manual/analysisComparison.py can compare any two sim dirs.
		# sim_data1.operons_on and sim_data2.operons_on indicate operons on/off.

		# todo: this actually lets you get the name of the input data !!!!!
		reference_sim_name = os.path.basename(os.path.normpath(reference_sim_dir))
		input_sim_name = os.path.basename(os.path.normpath(input_sim_dir))




		# noinspection PyUnusedLocal
		ap1, sim_data1, validation_data1 = self.setup(reference_sim_dir)
		# noinspection PyUnusedLocal
		ap2, sim_data2, _ = self.setup(input_sim_dir)

		if ap1.n_generation <= 2 or ap2.n_generation <= 2:
			print('Skipping analysis -- not enough sims run.')
			return


		avg_total_counts1, avg_free_counts1, avg_count_for_monomers_in_complexes1, avg_complex_counts1 = self.generate_data(reference_sim_dir)
		avg_total_counts2, avg_free_counts2, avg_count_for_monomers_in_complexes2, avg_complex_counts2 = self.generate_data(input_sim_dir)

		# correct for the free monomer:
		if avg_free_counts1.shape > avg_free_counts2.shape:
			diff = avg_free_counts1.shape[0] - avg_free_counts2.shape[0]
			counts_to_add = np.zeros(diff)
			avg_total_counts2 = np.append(avg_total_counts2, counts_to_add)
			avg_free_counts2 = np.append(avg_free_counts2, counts_to_add)
			avg_count_for_monomers_in_complexes2 = np.append(avg_count_for_monomers_in_complexes2, counts_to_add)

		if avg_free_counts2.shape > avg_free_counts1.shape:
			diff = avg_free_counts2.shape[0] - avg_free_counts1.shape[0]
			counts_to_add = np.zeros(diff)
			avg_total_counts1 = np.append(avg_total_counts1, counts_to_add)
			avg_free_counts1 = np.append(avg_free_counts1, counts_to_add)
			avg_count_for_monomers_in_complexes1 = np.append(avg_count_for_monomers_in_complexes1,
															 counts_to_add)



		protein_pearson1 = pearsonr(
			np.log10(avg_total_counts1 + 1),
			np.log10(avg_total_counts2 + 1))
		protein_pearson2 = pearsonr(
			np.log10(avg_free_counts1 + 1),
			np.log10(avg_free_counts2 + 1))

		protein_pearson3 = pearsonr(
			np.log10(avg_count_for_monomers_in_complexes1 + 1),
			np.log10(avg_count_for_monomers_in_complexes2 + 1))

		protein_pearson4 = pearsonr(
			np.log10(avg_complex_counts1 + 1),
			np.log10(avg_complex_counts2 + 1))



		plt.figure(figsize=(8, 8.4))


		def plot_total_protein_counts(ax, reference_count, input_count, pearson):
			ax.plot([0, 6], [0, 6], ls='--', lw=1, c='k', alpha=0.05)
			ax.scatter(
				np.log10(reference_count + 1),
				np.log10(input_count + 1), clip_on=False,
				c='#555555', edgecolor='none', s=20, alpha=0.25,
			)
			ax.set_title(f"$R^2$ = {pearson[0] ** 2:.3f}, n = {len(reference_count)}")
			ax.set_xlabel(
				f"log10(Mean simulated protein counts + 1),\n (total monomer counts) {reference_sim_name}")
			ax.set_ylabel(
				f"log10(Mean simulated protein counts + 1),\n (total monomer counts) {input_sim_name}")
			ax.spines["top"].set_visible(False)
			ax.spines["right"].set_visible(False)
			ax.spines["bottom"].set_position(("outward", 15))
			ax.spines["left"].set_position(("outward", 15))
			ax.set_xlim([0, 6])
			ax.set_ylim([0, 6])

		def plot_free_monomer_counts(ax, reference_count, input_count, pearson):
			ax.plot([0, 6], [0, 6], ls='--', lw=1, c='k', alpha=0.05)
			ax.scatter(
				np.log10(reference_count + 1),
				np.log10(input_count + 1), clip_on=False,
				c='#555555', edgecolor='none', s=20, alpha=0.25,
			)
			ax.set_title(f"$R^2$ = {pearson[0] ** 2:.3f}, n = {len(reference_count)}")
			ax.set_xlabel(
				f"log10(Mean simulated protein counts + 1),\n (free monomer counts) {reference_sim_name}")
			ax.set_ylabel(
				f"log10(Mean simulated protein counts + 1),\n (free monomer counts) {input_sim_name}")
			ax.spines["top"].set_visible(False)
			ax.spines["right"].set_visible(False)
			ax.spines["bottom"].set_position(("outward", 15))
			ax.spines["left"].set_position(("outward", 15))
			ax.set_xlim([0, 6])
			ax.set_ylim([0, 6])

		def plot_complexed_monomer_counts(ax, reference_count, input_count, pearson):
			ax.plot([0, 6], [0, 6], ls='--', lw=1, c='k', alpha=0.05)
			ax.scatter(
				np.log10(reference_count + 1),
				np.log10(input_count + 1), clip_on=False,
				c='#555555', edgecolor='none', s=20, alpha=0.25,
			)
			ax.set_title(f"$R^2$ = {pearson[0] ** 2:.3f}, n = {len(reference_count)}")
			ax.set_xlabel(
				f"log10(Mean simulated protein counts + 1),\n (complexed monomer counts) {reference_sim_name}")
			ax.set_ylabel(
				f"log10(Mean simulated protein counts + 1),\n (complexed monomer counts) {input_sim_name}")
			ax.spines["top"].set_visible(False)
			ax.spines["right"].set_visible(False)
			ax.spines["bottom"].set_position(("outward", 15))
			ax.spines["left"].set_position(("outward", 15))
			ax.set_xlim([0, 6])
			ax.set_ylim([0, 6])

		def plot_complex_counts(ax, reference_count, input_count, pearson):
			ax.plot([0, 6], [0, 6], ls='--', lw=1, c='k', alpha=0.05)
			ax.scatter(
				np.log10(reference_count + 1),
				np.log10(input_count + 1), clip_on=False,
				c='#555555', edgecolor='none', s=20, alpha=0.25,
			)
			ax.set_title(f"$R^2$ = {pearson[0] ** 2:.3f}, n = {len(reference_count)}")
			ax.set_xlabel(
				f"log10(Mean simulated complex counts + 1),\n {reference_sim_name}")
			ax.set_ylabel(
				f"log10(Mean simulated complex counts + 1),\n {input_sim_name}")
			ax.spines["top"].set_visible(False)
			ax.spines["right"].set_visible(False)
			ax.spines["bottom"].set_position(("outward", 15))
			ax.spines["left"].set_position(("outward", 15))
			ax.set_xlim([0, 6])
			ax.set_ylim([0, 6])




		ax1 = plt.subplot(2, 2, 1)
		ax2 = plt.subplot(2, 2, 2)
		plot_total_protein_counts(ax1, avg_total_counts1, avg_total_counts2, protein_pearson1)
		plot_free_monomer_counts(ax2, avg_free_counts1, avg_free_counts2, protein_pearson2)
		ax3 = plt.subplot(2, 2, 3)
		ax4 = plt.subplot(2, 2, 4)
		plot_complexed_monomer_counts(ax3, avg_count_for_monomers_in_complexes1, avg_count_for_monomers_in_complexes2, protein_pearson3)
		plot_complex_counts(ax4, avg_complex_counts1, avg_complex_counts2, protein_pearson4)




		out_file_name = \
			f"{plotOutFileName}_{reference_sim_name}_vs_{input_sim_name}_protein_count_comparisons"
		plt.tight_layout()
		exportFigure(plt, plotOutDir, out_file_name, metadata)
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
