"""
Analyze which fluxes differ the most across variants.
Originally created to debug high kinetic weight sims failing
in initial generation.
"""

import pickle
import os

import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import pandas as pd

from models.ecoli.analysis import variantAnalysisPlot
from wholecell.analysis.analysis_tools import (
	exportFigure, read_stacked_columns, read_stacked_bulk_molecules)
from wholecell.io.tablereader import TableReader


class Plot(variantAnalysisPlot.VariantAnalysisPlot):
	def do_plot(self, inputDir, plotOutDir, plotOutFileName, simDataFile,
				validationDataFile, metadata):

		with open(simDataFile, 'rb') as f:
			sim_data = pickle.load(f)
		plotOutDir = os.path.join(plotOutDir, plotOutFileName)

		if not os.path.exists(plotOutDir):
			os.mkdir(plotOutDir)

		var0_cells = self.ap.get_cells(variant=[0], generation=[0])
		fba_results_reader = TableReader(
			os.path.join(var0_cells[0], "simOut", "FBAResults"))
		reaction_ids = fba_results_reader.readAttribute("reactionIDs")
		fba_results_reader.close()
		avg_var0_fluxes = read_stacked_columns(
			var0_cells, "fbaResults",
			"reactionFluxes").mean(axis=0)

		avg_fluxes_all_var = avg_var0_fluxes.copy()

		for variant in self.ap.get_variants():
			if variant == 0:
				continue
			var_cells = self.ap.get_cells(variant=[variant], generation=[0])
			variant_fluxes = read_stacked_columns(
				var_cells, "fbaResults",
				"reactionFluxes").mean(axis=0)
			avg_fluxes_all_var = np.vstack(
				[avg_fluxes_all_var, variant_fluxes])

		variant_ids = self.ap.get_variants()
		variant_ids_str = [str(vid) for vid in variant_ids]
		avg_fluxes_all_var_df = pd.DataFrame(
			data=avg_fluxes_all_var,
			index=variant_ids_str,
			columns=reaction_ids).T

		df = avg_fluxes_all_var_df.copy()
		success_conditions = ["0", "5"]
		failure_conditions = ["7", "8"]

		df['mean_success'] = df[success_conditions].mean(axis=1)
		df['mean_failure'] = df[failure_conditions].mean(axis=1)
		df['delta'] = df['mean_failure'] - df['mean_success']
		df['abs_delta'] = df['delta'].abs()
		top_changed = df.sort_values('abs_delta', ascending=False)

		success_nonzero = (df[success_conditions].abs().sum(axis=1) > 0)
		failure_nonzero = (df[failure_conditions].abs().sum(axis=1) > 0)
		turned_off = df[success_nonzero & ~failure_nonzero].copy()
		turned_off = turned_off.sort_values('abs_delta', ascending=False)
		turned_on = df[~success_nonzero & failure_nonzero].copy()
		turned_on = turned_on.sort_values('abs_delta', ascending=False)

		sign_flip = df[
			np.sign(df['mean_success']) != np.sign(df['mean_failure'])
		].sort_values('abs_delta', ascending=False)

		top_changed_file = os.path.join(
			plotOutDir, "top_changed_fluxes.csv")
		top_changed.to_csv(top_changed_file)
		turned_off_file = os.path.join(
			plotOutDir, "turned_off_fluxes.csv")
		turned_off.to_csv(turned_off_file)
		turned_on_file = os.path.join(
			plotOutDir, "turned_on_fluxes.csv")
		turned_on.to_csv(turned_on_file)
		sign_flip_file = os.path.join(
			plotOutDir, "sign_flip_fluxes.csv")
		sign_flip.to_csv(sign_flip_file)


if __name__ == "__main__":
	Plot().cli()
