"""
Analyzes parameters sensitivity from running variant param_sensitivity

@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 5/17/19
"""

from __future__ import absolute_import
from __future__ import division

import cPickle
from matplotlib import pyplot as plt
import numpy as np
import os

from models.ecoli.analysis import variantAnalysisPlot
from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from wholecell.analysis.analysis_tools import exportFigure
from wholecell.io.tablereader import TableReader
from wholecell.utils import filepath


class Plot(variantAnalysisPlot.VariantAnalysisPlot):
	def do_plot(self, inputDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		if metadata.get('variant', '') != 'param_sensitivity':
			print 'This plot only runs for the param_sensitivity variant.'
			return

		if not os.path.isdir(inputDir):
			raise Exception, 'inputDir does not currently exist as a directory'

		filepath.makedirs(plotOutDir)

		ap = AnalysisPaths(inputDir, variant_plot=True)
		variants = ap.get_variants()

		initialize_from_sim_data = True
		for variant in variants:
			with open(ap.get_variant_kb(variant), 'rb') as f:
				sim_data = cPickle.load(f)

			# Populate values from the first sim_data object read
			if initialize_from_sim_data:
				initialize_from_sim_data = False
				rna_deg_increase_growth_rate = np.zeros(len(sim_data.process.transcription.rnaData['degRate']))
				rna_deg_decrease_growth_rate = np.zeros(len(sim_data.process.transcription.rnaData['degRate']))
				rna_deg_increase_counts = np.zeros(len(sim_data.process.transcription.rnaData['degRate']))
				rna_deg_decrease_counts = np.zeros(len(sim_data.process.transcription.rnaData['degRate']))
				trans_eff_increase_growth_rate = np.zeros(len(sim_data.process.translation.translationEfficienciesByMonomer))
				trans_eff_decrease_growth_rate = np.zeros(len(sim_data.process.translation.translationEfficienciesByMonomer))
				trans_eff_increase_counts = np.zeros(len(sim_data.process.translation.translationEfficienciesByMonomer))
				trans_eff_decrease_counts = np.zeros(len(sim_data.process.translation.translationEfficienciesByMonomer))
				protein_deg_growth_rate = None

			# TODO: save indices separately from sim_data for faster load
			rna_deg_increase_indices = sim_data.increase_rna_deg_indices
			# sim_data.increase_protein_deg_indices
			trans_eff_increase_indices = sim_data.increase_trans_eff_indices
			# sim_data.increase_synth_prob_indices
			rna_deg_decrease_indices = sim_data.decrease_rna_deg_indices
			# sim_data.decrease_protein_deg_indices
			trans_eff_decrease_indices = sim_data.decrease_trans_eff_indices
			# sim_data.decrease_synth_prob_indices
			if len(rna_deg_increase_indices) == 0:
				continue

			for sim_dir in ap.get_cells(variant=[variant]):
				simOutDir = os.path.join(sim_dir, "simOut")

				# Listeners used
				try:
					main_reader = TableReader(os.path.join(simOutDir, 'Main'))
					mass_reader = TableReader(os.path.join(simOutDir, 'Mass'))

					# Load data
					time = main_reader.readColumn('time')
					cell_mass = mass_reader.readColumn('cellMass')
					growth_rate = np.nanmean(mass_reader.readColumn('instantaniousGrowthRate'))
				except:
					continue

				mass_diff = (cell_mass[-1] - cell_mass[0]) / (time[-1] - time[0])

				rna_deg_increase_growth_rate[rna_deg_increase_indices] += growth_rate
				rna_deg_decrease_growth_rate[rna_deg_decrease_indices] += growth_rate
				rna_deg_increase_counts[rna_deg_increase_indices] += 1
				rna_deg_decrease_counts[rna_deg_decrease_indices] += 1
				trans_eff_increase_growth_rate[trans_eff_increase_indices] += growth_rate
				trans_eff_decrease_growth_rate[trans_eff_decrease_indices] += growth_rate
				trans_eff_increase_counts[trans_eff_increase_indices] += 1
				trans_eff_decrease_counts[trans_eff_decrease_indices] += 1

		plt.figure()

		### Create Plot ###

		rna_deg_growth_rate = rna_deg_increase_growth_rate / rna_deg_increase_counts - rna_deg_decrease_growth_rate / rna_deg_decrease_counts
		trans_eff_growth_rate = trans_eff_increase_growth_rate / trans_eff_increase_counts - trans_eff_decrease_growth_rate / trans_eff_decrease_counts

		data = trans_eff_growth_rate
		mean = data.mean()
		std = data.std()

		plt.bar(range(len(data)), np.sort(data))
		plt.axhline(mean + 3*std, color='r')
		plt.axhline(mean - 3*std, color='r')

		exportFigure(plt, plotOutDir, plotOutFileName, metadata)

		sorted_idx = np.argsort(data)
		ids = sim_data.process.translation.monomerData['id']
		print(np.sum(data > mean + 3*std))
		print(np.sum(data < mean - 3*std))
		print(ids[sorted_idx[:10]])
		print(ids[sorted_idx[-10:]])
		import ipdb; ipdb.set_trace()

		plt.close('all')


if __name__ == "__main__":
	Plot().cli()
