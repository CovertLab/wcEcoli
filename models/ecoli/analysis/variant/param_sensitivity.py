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
from wholecell.utils import constants, filepath, parallelization


def analyze_variant((variant, ap)):
	'''
	Method to map each variant to for parallel analysis.

	Args:
		variant (int): variant index
		ap (AnalysisPaths object): variant plot analysis paths object to get
			simulation output directories for the variant

	Returns:
		16 ndarray[float]: growth rate sums associated to a parameter being increased
			or decreased and count of number of times each parameter was increased
			or decreased
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

		try:
			# Listeners used
			mass_reader = TableReader(os.path.join(simOutDir, 'Mass'))

			# Load data
			growth_rate = np.nanmean(mass_reader.readColumn('instantaniousGrowthRate')[-5:])
		except:
			# Exclude failed sims
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

def rna_mapping(sim_data, monomer_to_gene):
	'''
	Maps RNA ids to gene symbols.

	Args:
		sim_data (SimulationData object)
		monomer_to_gene ({monomer id (str): gene id (str)}: mapping of monomer
			id to gene id, with no location tags for keys and values

	Returns:
		{rna with location tag (str): gene symbol without location tag (str)}:
			mapping of RNA to gene symbol
	'''

	gene_data = sim_data.process.replication.geneData
	rna_data = sim_data.process.transcription.rnaData
	monomer_data = sim_data.process.translation.monomerData

	rna_to_monomer = {rna: monomer for rna, monomer in zip(monomer_data['rnaId'], monomer_data['id'])}
	rna_to_gene_name = {rna: gene for rna, gene in zip(gene_data['rnaId'], gene_data['name'])}

	return {rna: monomer_to_gene.get(rna_to_monomer.get(rna, '')[:-3],
		rna_to_gene_name.get(rna[:-3], rna))
		for rna in rna_data['id']}


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
		with open(os.path.join(inputDir, 'kb', constants.SERIALIZED_FIT1_FILENAME)) as f:
			sim_data = cPickle.load(f)

		pool = Pool(processes=parallelization.plotter_cpus())
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


		rna_deg_growth_rate = rna_deg_increase_growth_rate / rna_deg_increase_counts - rna_deg_decrease_growth_rate / rna_deg_decrease_counts
		protein_deg_growth_rate = protein_deg_increase_growth_rate / protein_deg_increase_counts - protein_deg_decrease_growth_rate / protein_deg_decrease_counts
		trans_eff_growth_rate = trans_eff_increase_growth_rate / trans_eff_increase_counts - trans_eff_decrease_growth_rate / trans_eff_decrease_counts
		synth_prob_growth_rate = synth_prob_increase_growth_rate / synth_prob_increase_counts - synth_prob_decrease_growth_rate / synth_prob_decrease_counts

		data = np.hstack((rna_deg_growth_rate, protein_deg_growth_rate, trans_eff_growth_rate, synth_prob_growth_rate))
		mean = data[np.isfinite(data)].mean()
		std = data[np.isfinite(data)].std()

		n_std = 4
		upper_threshold = mean + n_std*std
		lower_threshold = mean - n_std*std

		plt.figure()

		plt.bar(range(len(data)), np.sort(data))
		plt.axhline(upper_threshold , color='r')
		plt.axhline(lower_threshold, color='r')

		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close('all')

		print('Number of params above threshold: {}'.format(np.sum(data > upper_threshold)))
		print('Number of params below threshold: {}'.format(np.sum(data < lower_threshold)))

		rna_ids = sim_data.process.transcription.rnaData['id']
		monomer_ids = sim_data.process.translation.monomerData['id']

		monomer_to_gene = sim_data.moleculeGroups.frameIDGeneSymbol_Dict
		rna_to_gene = rna_mapping(sim_data, monomer_to_gene)  # TODO: create dict in sim_data from raw_data

		print('Positive correlation between RNA deg and growth:')
		mask = rna_deg_growth_rate > upper_threshold
		for rna_id, effect in zip(rna_ids[mask], rna_deg_growth_rate[mask]):
			print('\t{}: {:.2f}'.format(rna_to_gene.get(rna_id, rna_id), (effect - mean) / std))
		print('Negative correlation between RNA deg and growth:')
		mask = rna_deg_growth_rate < lower_threshold
		for rna_id, effect in zip(rna_ids[mask], rna_deg_growth_rate[mask]):
			print('\t{}: {:.2f}'.format(rna_to_gene.get(rna_id, rna_id), (effect - mean) / std))

		print('Positive correlation between protein deg and growth:')
		mask = protein_deg_growth_rate > upper_threshold
		for monomer_id, effect in zip(monomer_ids[mask], protein_deg_growth_rate[mask]):
			print('\t{}: {:.2f}'.format(monomer_to_gene.get(monomer_id[:-3], monomer_id), (effect - mean) / std))
		print('Negative correlation between protein deg and growth:')
		mask = protein_deg_growth_rate < lower_threshold
		for monomer_id, effect in zip(monomer_ids[mask], protein_deg_growth_rate[mask]):
			print('\t{}: {:.2f}'.format(monomer_to_gene.get(monomer_id[:-3], monomer_id), (effect - mean) / std))

		print('Positive correlation between translation efficiency and growth:')
		mask = trans_eff_growth_rate > upper_threshold
		for monomer_id, effect in zip(monomer_ids[mask], trans_eff_growth_rate[mask]):
			print('\t{}: {:.2f}'.format(monomer_to_gene.get(monomer_id[:-3], monomer_id), (effect - mean) / std))
		print('Negative correlation between translation efficiency and growth:')
		mask = trans_eff_growth_rate < lower_threshold
		for monomer_id, effect in zip(monomer_ids[mask], trans_eff_growth_rate[mask]):
			print('\t{}: {:.2f}'.format(monomer_to_gene.get(monomer_id[:-3], monomer_id), (effect - mean) / std))

		print('Positive correlation between synthesis probability and growth:')
		mask = synth_prob_growth_rate > upper_threshold
		for rna_id, effect in zip(rna_ids[mask], synth_prob_growth_rate[mask]):
			print('\t{}: {:.2f}'.format(rna_to_gene.get(rna_id, rna_id), (effect - mean) / std))
		print('Negative correlation between synthesis probability and growth:')
		mask = synth_prob_growth_rate < lower_threshold
		for rna_id, effect in zip(rna_ids[mask], synth_prob_growth_rate[mask]):
			print('\t{}: {:.2f}'.format(rna_to_gene.get(rna_id, rna_id), (effect - mean) / std))


if __name__ == "__main__":
	Plot().cli()
