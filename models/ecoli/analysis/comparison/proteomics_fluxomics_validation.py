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


class Plot(comparisonAnalysisPlot.ComparisonAnalysisPlot):
	def do_plot(self, reference_sim_dir, plotOutDir, plotOutFileName, input_sim_dir, unused, metadata):
		# From fw_queue, reference_sim_dir has operons="off"; input_sim_dir has
		# operons="on".
		# manual/analysisComparison.py can compare any two sim dirs.
		# sim_data1.operons_on and sim_data2.operons_on indicate operons on/off.

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
			cell_paths = ap.get_cells(generation=np.arange(2, ap.n_generation))

			monomer_counts = read_stacked_columns(
				cell_paths, 'MonomerCounts', 'monomerCounts', ignore_exception=True)

			sim_monomer_counts, val_monomer_counts_filtered = get_simulated_validation_counts(
				val_monomer_counts, monomer_counts, val_monomer_ids, monomer_ids)

			return sim_monomer_counts, val_monomer_counts_filtered


		sim1_monomer_counts, val1_monomer_counts = read_sim_protein_counts(ap1)
		sim2_monomer_counts, val2_monomer_counts = read_sim_protein_counts(ap2)


		def read_sim_fluxes(ap):
			# Ignore data from first two gens
			cell_paths = ap.get_cells(generation=np.arange(2, ap.n_generation))

			# Calculate coefficients to be used to convert flux units from
			# mM/s to mmol/gCDW/h
			cell_mass = read_stacked_columns(
				cell_paths, 'Mass', 'cellMass', ignore_exception=True)
			dry_mass = read_stacked_columns(
				cell_paths, 'Mass', 'dryMass', ignore_exception=True)
			conversion_coeffs = (
				dry_mass / cell_mass
				* cell_density.asNumber(MASS_UNITS / VOLUME_UNITS)
			)

			# Get flux for each reaction at each timestep
			sim_fluxes = (
				(COUNTS_UNITS / MASS_UNITS / TIME_UNITS)
				* (read_stacked_columns(cell_paths, 'FBAResults', 'base_reaction_fluxes', ignore_exception=True) / conversion_coeffs)
				).asNumber(units.mmol / units.g / units.h)

			# Add up all fluxes that contribute to each value in validation data
			rxn_id_to_sim_flux_mean = {}
			rxn_id_to_sim_flux_std = {}

			for val_rxn_id in val_rxn_ids:
				if val_rxn_id in sim_rxn_id_to_index:
					fluxes = sim_fluxes[:, sim_rxn_id_to_index[val_rxn_id]]
					rxn_id_to_sim_flux_mean[val_rxn_id] = fluxes.mean()
					rxn_id_to_sim_flux_std[val_rxn_id] = fluxes.std()

			return rxn_id_to_sim_flux_mean, rxn_id_to_sim_flux_std


		rxn_id_to_sim1_flux_mean, rxn_id_to_sim1_flux_std = read_sim_fluxes(ap1)
		rxn_id_to_sim2_flux_mean, rxn_id_to_sim2_flux_std = read_sim_fluxes(ap2)

		# Plot reactions that exist in both sets of sims and validation data
		plotted_rxn_ids = sorted(list(
			rxn_id_to_val_flux_mean.keys()
			& rxn_id_to_sim1_flux_mean.keys() & rxn_id_to_sim2_flux_mean.keys()
			))
		val_flux_mean = np.array([
			rxn_id_to_val_flux_mean[rxn_id] for rxn_id in plotted_rxn_ids])
		sim1_flux_mean = np.array([
			rxn_id_to_sim1_flux_mean[rxn_id] for rxn_id in plotted_rxn_ids])
		sim2_flux_mean = np.array([
			rxn_id_to_sim2_flux_mean[rxn_id] for rxn_id in plotted_rxn_ids])
		val_flux_std = np.array([
			rxn_id_to_val_flux_std[rxn_id] for rxn_id in plotted_rxn_ids])
		sim1_flux_std = np.array([
			rxn_id_to_sim1_flux_std[rxn_id] for rxn_id in plotted_rxn_ids])
		sim2_flux_std = np.array([
			rxn_id_to_sim2_flux_std[rxn_id] for rxn_id in plotted_rxn_ids])

		protein_pearson1 = pearsonr(
			np.log10(val1_monomer_counts + 1),
			np.log10(sim1_monomer_counts + 1))
		protein_pearson2 = pearsonr(
			np.log10(val2_monomer_counts + 1),
			np.log10(sim2_monomer_counts + 1))
		flux_pearson1 = pearsonr(val_flux_mean, sim1_flux_mean)
		flux_pearson2 = pearsonr(val_flux_mean, sim2_flux_mean)

		plt.figure(figsize=(8, 8.4))

		def plot_protein_counts(ax, val, sim, pearson):
			ax.plot([0, 6], [0, 6], ls='--', lw=1, c='k', alpha=0.05)
			ax.scatter(
				np.log10(val + 1),
				np.log10(sim + 1), clip_on=False,
				c='#555555', edgecolor='none', s=20, alpha=0.25,
				)
			ax.set_title(f"$R^2$ = {pearson[0] ** 2:.2f}")
			ax.set_xlabel("log10(Schmidt 2015 protein counts + 1)")
			ax.set_ylabel("log10(Mean simulated protein counts + 1)")
			ax.spines["top"].set_visible(False)
			ax.spines["right"].set_visible(False)
			ax.spines["bottom"].set_position(("outward", 15))
			ax.spines["left"].set_position(("outward", 15))
			ax.set_xlim([0, 6])
			ax.set_ylim([0, 6])

		ax1 = plt.subplot(2, 2, 1)
		ax2 = plt.subplot(2, 2, 2)
		plot_protein_counts(ax1, val1_monomer_counts, sim1_monomer_counts, protein_pearson1)
		plot_protein_counts(ax2, val2_monomer_counts, sim2_monomer_counts, protein_pearson2)

		def plot_flux(ax, val, sim, val_std, sim_std, pearson):
			ax.plot([-15, 25], [-15, 25], ls='--', lw=1, c='k', alpha=0.05)
			ax.errorbar(
				val, sim,
				xerr=val_std, yerr=sim_std,
				c='#555555', ms=4, fmt="o", ecolor="#cccccc")
			ax.set_title(f"$R^2$ = {pearson[0] ** 2:.2f}, p = {pearson[1]:.2e}")
			ax.set_xlabel("Toya 2010 flux [mmol/g/hr]")
			ax.set_ylabel("Mean simulated flux [mmol/g/hr]")
			ax.spines["top"].set_visible(False)
			ax.spines["right"].set_visible(False)
			ax.spines["bottom"].set_position(("outward", 15))
			ax.spines["left"].set_position(("outward", 15))
			ax.set_xticks([-15, -10, -5, 0, 5, 10, 15, 20, 25])
			ax.set_xlim([-15, 25])
			ax.set_ylim([-15, 25])

		ax3 = plt.subplot(2, 2, 3)
		ax4 = plt.subplot(2, 2, 4)
		plot_flux(ax3, val_flux_mean, sim1_flux_mean, val_flux_std, sim1_flux_std, flux_pearson1)
		plot_flux(ax4, val_flux_mean, sim2_flux_mean, val_flux_std, sim2_flux_std, flux_pearson2)

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
