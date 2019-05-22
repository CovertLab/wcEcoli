"""
Analyzes parameters sensitivity from running variant param_sensitivity

@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 5/17/19
"""

from __future__ import absolute_import
from __future__ import division

import cPickle
from matplotlib import pyplot as plt
from multiprocessing import Pool
import numpy as np
import os

from models.ecoli.analysis import variantAnalysisPlot
from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from models.ecoli.sim.variants.param_sensitivity import number_params, param_indices, split_indices
from wholecell.analysis.analysis_tools import exportFigure
from wholecell.io.tablereader import TableReader
from wholecell.utils import filepath, parallelization


def analyze_variant((variant, ap)):
	'''

	'''

	increase_indices, decrease_indices = split_indices(sim_data, variant)
	(rna_deg_increase_indices,
		protein_deg_increase_indices,
		trans_eff_increase_indices,
		synth_prob_increase_indices) = param_indices(sim_data, increase_indices)
	(rna_deg_decrease_indices,
		protein_deg_decrease_indices,
		trans_eff_decrease_indices,
		synth_prob_decrease_indices) = param_indices(sim_data, decrease_indices)

	(n_rna_deg_rates,
		n_protein_deg_rates,
		n_translation_efficiencies,
		n_synth_prob) = number_params(sim_data)

	rna_deg_increase_growth_rate = np.zeros(n_rna_deg_rates)
	rna_deg_decrease_growth_rate = np.zeros(n_rna_deg_rates)
	rna_deg_increase_counts = np.zeros(n_rna_deg_rates)
	rna_deg_decrease_counts = np.zeros(n_rna_deg_rates)
	protein_deg_increase_growth_rate = np.zeros(n_protein_deg_rates)
	protein_deg_decrease_growth_rate = np.zeros(n_protein_deg_rates)
	protein_deg_increase_counts = np.zeros(n_protein_deg_rates)
	protein_deg_decrease_counts = np.zeros(n_protein_deg_rates)
	trans_eff_increase_growth_rate = np.zeros(n_translation_efficiencies)
	trans_eff_decrease_growth_rate = np.zeros(n_translation_efficiencies)
	trans_eff_increase_counts = np.zeros(n_translation_efficiencies)
	trans_eff_decrease_counts = np.zeros(n_translation_efficiencies)
	synth_prob_increase_growth_rate = np.zeros(n_synth_prob)
	synth_prob_decrease_growth_rate = np.zeros(n_synth_prob)
	synth_prob_increase_counts = np.zeros(n_synth_prob)
	synth_prob_decrease_counts = np.zeros(n_synth_prob)

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

			mass_diff = (cell_mass[-1] - cell_mass[0]) / (time[-1] - time[0])

		except:
			continue

		rna_deg_increase_growth_rate[rna_deg_increase_indices] += growth_rate
		rna_deg_decrease_growth_rate[rna_deg_decrease_indices] += growth_rate
		rna_deg_increase_counts[rna_deg_increase_indices] += 1
		rna_deg_decrease_counts[rna_deg_decrease_indices] += 1
		protein_deg_increase_growth_rate[protein_deg_increase_indices] += growth_rate
		protein_deg_decrease_growth_rate[protein_deg_decrease_indices] += growth_rate
		protein_deg_increase_counts[protein_deg_increase_indices] += 1
		protein_deg_decrease_counts[protein_deg_decrease_indices] += 1
		trans_eff_increase_growth_rate[trans_eff_increase_indices] += growth_rate
		trans_eff_decrease_growth_rate[trans_eff_decrease_indices] += growth_rate
		trans_eff_increase_counts[trans_eff_increase_indices] += 1
		trans_eff_decrease_counts[trans_eff_decrease_indices] += 1
		synth_prob_increase_growth_rate[synth_prob_increase_indices] += growth_rate
		synth_prob_decrease_growth_rate[synth_prob_decrease_indices] += growth_rate
		synth_prob_increase_counts[synth_prob_increase_indices] += 1
		synth_prob_decrease_counts[synth_prob_decrease_indices] += 1

	return (rna_deg_increase_growth_rate, rna_deg_decrease_growth_rate, rna_deg_increase_counts, rna_deg_decrease_counts,
		protein_deg_increase_growth_rate, protein_deg_decrease_growth_rate, protein_deg_increase_counts, protein_deg_decrease_counts,
		trans_eff_increase_growth_rate, trans_eff_decrease_growth_rate, trans_eff_increase_counts, trans_eff_decrease_counts,
		synth_prob_increase_growth_rate, synth_prob_decrease_growth_rate, synth_prob_increase_counts, synth_prob_decrease_counts)

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
		n_variants = len(variants)

		# Load one instance of sim_data to get number of parameters and ids
		global sim_data
		with open(ap.get_variant_kb(variants[0])) as f:
			sim_data = cPickle.load(f)

		pool = Pool(processes=min(8, parallelization.plotter_cpus()))
		args = zip(
			variants,
			[ap] * n_variants,
			)
		results = pool.map(analyze_variant, args)
		pool.close()
		pool.join()
		initialize = True
		for i, result in enumerate(results):
			(_rna_deg_increase_growth_rate,
				_rna_deg_decrease_growth_rate,
				_rna_deg_increase_counts,
				_rna_deg_decrease_counts,
				_protein_deg_increase_growth_rate,
				_protein_deg_decrease_growth_rate,
				_protein_deg_increase_counts,
				_protein_deg_decrease_counts,
				_trans_eff_increase_growth_rate,
				_trans_eff_decrease_growth_rate,
				_trans_eff_increase_counts,
				_trans_eff_decrease_counts,
				_synth_prob_increase_growth_rate,
				_synth_prob_decrease_growth_rate,
				_synth_prob_increase_counts,
				_synth_prob_decrease_counts) = result

			if initialize:
				initialize = False
				rna_deg_increase_growth_rate = _rna_deg_increase_growth_rate
				rna_deg_decrease_growth_rate = _rna_deg_decrease_growth_rate
				rna_deg_increase_counts = _rna_deg_increase_counts
				rna_deg_decrease_counts = _rna_deg_decrease_counts
				protein_deg_increase_growth_rate = _protein_deg_increase_growth_rate
				protein_deg_decrease_growth_rate = _protein_deg_decrease_growth_rate
				protein_deg_increase_counts = _protein_deg_increase_counts
				protein_deg_decrease_counts = _protein_deg_decrease_counts
				trans_eff_increase_growth_rate = _trans_eff_increase_growth_rate
				trans_eff_decrease_growth_rate = _trans_eff_decrease_growth_rate
				trans_eff_increase_counts = _trans_eff_increase_counts
				trans_eff_decrease_counts = _trans_eff_decrease_counts
				synth_prob_increase_growth_rate = _synth_prob_increase_growth_rate
				synth_prob_decrease_growth_rate = _synth_prob_decrease_growth_rate
				synth_prob_increase_counts = _synth_prob_increase_counts
				synth_prob_decrease_counts = _synth_prob_decrease_counts
			else:
				rna_deg_increase_growth_rate += _rna_deg_increase_growth_rate
				rna_deg_decrease_growth_rate += _rna_deg_decrease_growth_rate
				rna_deg_increase_counts += _rna_deg_increase_counts
				rna_deg_decrease_counts += _rna_deg_decrease_counts
				protein_deg_increase_growth_rate += _protein_deg_increase_growth_rate
				protein_deg_decrease_growth_rate += _protein_deg_decrease_growth_rate
				protein_deg_increase_counts += _protein_deg_increase_counts
				protein_deg_decrease_counts += _protein_deg_decrease_counts
				trans_eff_increase_growth_rate += _trans_eff_increase_growth_rate
				trans_eff_decrease_growth_rate += _trans_eff_decrease_growth_rate
				trans_eff_increase_counts += _trans_eff_increase_counts
				trans_eff_decrease_counts += _trans_eff_decrease_counts
				synth_prob_increase_growth_rate += _synth_prob_increase_growth_rate
				synth_prob_decrease_growth_rate += _synth_prob_decrease_growth_rate
				synth_prob_increase_counts += _synth_prob_increase_counts
				synth_prob_decrease_counts += _synth_prob_decrease_counts


		plt.figure()

		### Create Plot ###

		rna_deg_growth_rate = rna_deg_increase_growth_rate / rna_deg_increase_counts - rna_deg_decrease_growth_rate / rna_deg_decrease_counts
		protein_deg_growth_rate = protein_deg_increase_growth_rate / protein_deg_increase_counts - protein_deg_decrease_growth_rate / protein_deg_decrease_counts
		trans_eff_growth_rate = trans_eff_increase_growth_rate / trans_eff_increase_counts - trans_eff_decrease_growth_rate / trans_eff_decrease_counts
		synth_prob_growth_rate = synth_prob_increase_growth_rate / synth_prob_increase_counts - synth_prob_decrease_growth_rate / synth_prob_decrease_counts


		data = np.hstack((rna_deg_growth_rate, protein_deg_growth_rate, trans_eff_growth_rate, synth_prob_growth_rate))
		mean = data[np.isfinite(data)].mean()
		std = data[np.isfinite(data)].std()

		plt.bar(range(len(data)), np.sort(data))
		plt.axhline(mean + 3*std, color='r')
		plt.axhline(mean - 3*std, color='r')

		exportFigure(plt, plotOutDir, plotOutFileName, metadata)

		sorted_idx = np.argsort(data)
		# ids = sim_data.process.translation.monomerData['id']
		ids = sim_data.process.transcription.rnaData['id']
		print(np.sum(data > mean + 3*std))
		print(np.sum(data < mean - 3*std))
		# print(ids[sorted_idx[:10]])
		# print(ids[sorted_idx[-10:]])
		import ipdb; ipdb.set_trace()

		plt.close('all')


if __name__ == "__main__":
	Plot().cli()
