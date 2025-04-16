"""
Plots proteomics/fluxomics validation plots for two sets of simulations.
"""

from typing import Tuple

from matplotlib import pyplot as plt
# noinspection PyUnresolvedReferences
import numpy as np
import os
from scipy.stats import pearsonr

from models.ecoli.analysis import comparisonAnalysisPlot
from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from models.ecoli.processes.metabolism import (
	COUNTS_UNITS, VOLUME_UNITS, TIME_UNITS, MASS_UNITS)
from reconstruction.ecoli.simulation_data import SimulationDataEcoli
from validation.ecoli.validation_data import ValidationDataEcoli
from wholecell.analysis.analysis_tools import exportFigure, read_stacked_columns
# noinspection PyUnresolvedReferences
from wholecell.io.tablereader import TableReader
from wholecell.utils.protein_counts import get_simulated_validation_counts
from wholecell.utils import units
SKIP_GENS = 2

class Plot(comparisonAnalysisPlot.ComparisonAnalysisPlot):
	def do_plot(self, reference_sim_dir, plotOutDir, plotOutFileName, input_sim_dir, unused, metadata):
		# From fw_queue, reference_sim_dir has operons="off"; input_sim_dir has
		# operons="on".
		# manual/analysisComparison.py can compare any two sim dirs.
		# sim_data1.operons_on and sim_data2.operons_on indicate operons on/off.


		reference_sim_name = os.path.basename(os.path.normpath(reference_sim_dir))
		input_sim_name = os.path.basename(os.path.normpath(input_sim_dir))


		# noinspection PyUnusedLocal
		ap1, sim_data1, validation_data1 = self.setup(reference_sim_dir)
		# noinspection PyUnusedLocal
		ap2, sim_data2, _ = self.setup(input_sim_dir)

		if ap1.n_generation <= 2 or ap2.n_generation <= 2:
			print('Skipping analysis -- not enough sims run.')
			return

		# Read data from sim_data
		monomer_ids = sim_data1.process.translation.monomer_data["id"]
		cell_density = sim_data1.constants.cell_density

		# Read validation protein counts
		val_monomer_ids = validation_data1.protein.schmidt2015Data["monomerId"]
		val_monomer_counts = validation_data1.protein.schmidt2015Data["glucoseCounts"]

		# Read validation fluxes
		val_rxn_ids = list(
			validation_data1.reactionFlux.toya2010fluxes["reactionID"])
		rxn_id_to_val_flux_mean = dict(zip(
			val_rxn_ids,
			[flux.asNumber(units.mmol / units.g / units.h) for flux
			in validation_data1.reactionFlux.toya2010fluxes["reactionFlux"]]
			))
		rxn_id_to_val_flux_std = dict(zip(
			val_rxn_ids,
			[flux.asNumber(units.mmol / units.g / units.h) for flux
			in validation_data1.reactionFlux.toya2010fluxes["reactionFluxStdev"]]
			))

		# Get list of relevant reaction IDs and indexes
		fba_reader = TableReader(
			os.path.join(ap1.get_cells()[0], 'simOut', 'FBAResults'))
		sim_rxn_ids = fba_reader.readAttribute('base_reaction_ids')
		sim_rxn_id_to_index = {
			rxn_id: i for (i, rxn_id) in enumerate(sim_rxn_ids)
		}

		def read_sim_protein_counts(ap):
			# Ignore data from first two gens
			cell_paths = ap.get_cells(generation=np.arange(SKIP_GENS, ap.n_generation))

			monomer_counts = read_stacked_columns(
				cell_paths, 'MonomerCounts', 'monomerCounts', ignore_exception=True)

			sim_monomer_counts, val_monomer_counts_filtered, overlapping_ids = get_simulated_validation_counts(
				val_monomer_counts, monomer_counts, val_monomer_ids, monomer_ids)

			# average the monomer counts over the generations
			all_sim_monomer_counts = monomer_counts.mean(axis=0)

			return all_sim_monomer_counts, sim_monomer_counts, val_monomer_counts_filtered


		sim1_all_monomer_counts, sim1_monomer_counts, val1_monomer_counts = read_sim_protein_counts(ap1)
		sim2_all_monomer_counts, sim2_monomer_counts, val2_monomer_counts = read_sim_protein_counts(ap2)



		protein_pearson1 = pearsonr(
			np.log10(val1_monomer_counts + 1),
			np.log10(sim1_monomer_counts + 1))
		protein_pearson2 = pearsonr(
			np.log10(val2_monomer_counts + 1),
			np.log10(sim2_monomer_counts + 1))

		protein_pearson3 = pearsonr(
			np.log10(sim1_monomer_counts + 1),
			np.log10(sim2_monomer_counts + 1))
		protein_pearson4 = pearsonr(
			np.log10(sim1_all_monomer_counts + 1),
			np.log10(sim2_all_monomer_counts + 1))





		plt.figure(figsize=(8, 8.4))

		def plot_sim_and_validation_counts(ax, val, sim, pearson, name):
			ax.plot([0, 6], [0, 6], ls='--', lw=1, c='k', alpha=0.05)
			ax.scatter(
				np.log10(val + 1),
				np.log10(sim + 1), clip_on=False,
				c='#555555', edgecolor='none', s=20, alpha=0.25,
				)
			ax.set_title(f"$R^2$ = {pearson[0] ** 2:.2f}, n = {len(val)}")
			ax.set_xlabel("log10(Schmidt 2015 protein counts + 1)")
			ax.set_ylabel(f"log10(Mean {name} simulation protein counts + 1)")
			ax.spines["top"].set_visible(False)
			ax.spines["right"].set_visible(False)
			ax.spines["bottom"].set_position(("outward", 15))
			ax.spines["left"].set_position(("outward", 15))
			ax.set_xlim([0, 6])
			ax.set_ylim([0, 6])

		ax1 = plt.subplot(2, 2, 1)
		ax2 = plt.subplot(2, 2, 2)
		plot_sim_and_validation_counts(
			ax1, val1_monomer_counts, sim1_monomer_counts, protein_pearson1,
			reference_sim_name)
		plot_sim_and_validation_counts(
			ax2, val2_monomer_counts, sim2_monomer_counts, protein_pearson2,
			input_sim_name)

		hi = 5
		def plot_validation_data_overlap_protein_counts(ax, reference_count, input_count, pearson):
			ax.plot([0, 6], [0, 6], ls='--', lw=1, c='k', alpha=0.05)
			ax.scatter(
				np.log10(reference_count + 1),
				np.log10(input_count + 1), clip_on=False,
				c='#555555', edgecolor='none', s=20, alpha=0.25,
			)
			ax.set_title(f"$R^2$ = {pearson[0] ** 2:.3f}, n = {len(reference_count)}")
			ax.set_xlabel(
				f"log10(Mean simulated protein counts + 1),\n (for proteins that overlap with validation data) {reference_sim_name}")
			ax.set_ylabel(
				f"log10(Mean simulated protein counts + 1),\n (for proteins that overlap with validation data) {input_sim_name}")
			ax.spines["top"].set_visible(False)
			ax.spines["right"].set_visible(False)
			ax.spines["bottom"].set_position(("outward", 15))
			ax.spines["left"].set_position(("outward", 15))
			ax.set_xlim([0, 6])
			ax.set_ylim([0, 6])

		def plot_all_protein_counts(ax, reference_count, input_count, pearson):
			ax.plot([0, 6], [0, 6], ls='--', lw=1, c='k', alpha=0.05)
			ax.scatter(
				np.log10(reference_count + 1),
				np.log10(input_count + 1), clip_on=False,
				c='#555555', edgecolor='none', s=20, alpha=0.25,
			)
			ax.set_title(f"$R^2$ = {pearson[0] ** 2:.3f}, n = {len(reference_count)}")
			ax.set_xlabel(
				f"log10(Mean simulated protein counts + 1),\n (for all proteins in the WCM) {reference_sim_name}")
			ax.set_ylabel(
				f"log10(Mean simulated protein counts + 1),\n (for all proteins in the WCM) {input_sim_name}")
			ax.spines["top"].set_visible(False)
			ax.spines["right"].set_visible(False)
			ax.spines["bottom"].set_position(("outward", 15))
			ax.spines["left"].set_position(("outward", 15))
			ax.set_xlim([0, 6])
			ax.set_ylim([0, 6])



		ax3 = plt.subplot(2, 2, 3)
		plot_validation_data_overlap_protein_counts(ax3, sim1_monomer_counts, sim2_monomer_counts, protein_pearson3)
		ax4 = plt.subplot(2, 2, 4)
		plot_all_protein_counts(ax4, sim1_all_monomer_counts, sim2_all_monomer_counts, protein_pearson4)



		out_file_name = \
			f"{plotOutFileName}_{reference_sim_name}_vs_{input_sim_name}_proteomics_validation_comparisons"
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
