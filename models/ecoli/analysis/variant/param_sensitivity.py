"""
Analyzes parameters sensitivity from running variant param_sensitivity

@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 5/17/19
"""

from __future__ import absolute_import
from __future__ import division

import cPickle
from multiprocessing import Pool
import os
import re

from matplotlib import pyplot as plt
import numpy as np
from scipy import stats

from models.ecoli.analysis import variantAnalysisPlot
from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from models.ecoli.processes.metabolism import COUNTS_UNITS, MASS_UNITS, TIME_UNITS, VOLUME_UNITS
from models.ecoli.sim.variants.param_sensitivity import number_params, split_indices
from wholecell.analysis.analysis_tools import exportFigure
from wholecell.io.tablereader import TableReader
from wholecell.utils import constants, filepath, parallelization, sparkline, units


def analyze_variant((variant, total_params)):
	'''
	Method to map each variant to for parallel analysis.

	Requires global variables:
		sim_data (SimulationData object)
		validation_data (ValidationData object)
		ap (AnalysisPaths object)

	Args:
		variant (int): variant index
		total_params (int): total number of parameters that are changed

	Returns:
		ndarray[float]: number of times each parameter was increased
		ndarray[float]: number of times each parameter was decreased
		ndarray[float]: average growth rate for each parameter when increaed
		ndarray[float]: average growth rate for each parameter when decreased
		ndarray[float]: average flux correlation for each parameter when increaed
		ndarray[float]: average flux correlation for each parameter when decreased
	'''

	increase_indices, decrease_indices = split_indices(sim_data, variant)

	increase_params_counts = np.zeros(total_params)
	decrease_params_counts = np.zeros(total_params)
	increase_params_growth_rate = np.zeros(total_params)
	decrease_params_growth_rate = np.zeros(total_params)
	increase_params_flux_correlation = np.zeros(total_params)
	decrease_params_flux_correlation = np.zeros(total_params)

	cell_density = sim_data.constants.cellDensity
	flux_units = units.mmol / units.g / units.h

	# Validation data
	toya_reactions = validation_data.reactionFlux.toya2010fluxes["reactionID"]
	toya_fluxes = [f.asNumber(flux_units) for f in validation_data.reactionFlux.toya2010fluxes["reactionFlux"]]

	if variant > 0:  # exclude control with no parameters changed
		for sim_dir in ap.get_cells(variant=[variant]):
			simOutDir = os.path.join(sim_dir, "simOut")

			try:
				# Listeners used
				mass_reader = TableReader(os.path.join(simOutDir, 'Mass'))
				fba_results_reader = TableReader(os.path.join(simOutDir, "FBAResults"))

				# Load data
				## Growth rate
				growth_rate = np.nanmean(mass_reader.readColumn('instantaniousGrowthRate')[-5:])

				## Central carbon flux
				dry_mass = mass_reader.readColumn('dryMass')[-5:]
				cell_mass = mass_reader.readColumn('cellMass')[-5:]
				coefficient = dry_mass / cell_mass * cell_density.asNumber(MASS_UNITS / VOLUME_UNITS)

				reaction_ids = np.array(fba_results_reader.readAttribute("reactionIDs"))
				reaction_fluxes = (COUNTS_UNITS / MASS_UNITS / TIME_UNITS) * (fba_results_reader.readColumn("reactionFluxes")[-5:].T / coefficient).T
			# Exclude failed sims
			except Exception as e:
				print(e)
				continue

			# Extract fluxes in Toya data set from simulation output
			model_fluxes = np.zeros_like(toya_fluxes)
			for i, toya_reaction in enumerate(toya_reactions):
				flux_time_course = []

				for rxn in reaction_ids:
					if re.findall(toya_reaction, rxn):
						reverse = 1
						if re.findall("(reverse)", rxn):
							reverse = -1

						if len(flux_time_course):
							flux_time_course += reverse * reaction_fluxes[:, np.where(reaction_ids == rxn)]
						else:
							flux_time_course = reverse * reaction_fluxes[:, np.where(reaction_ids == rxn)]

				if len(flux_time_course):
					model_fluxes[i] = np.mean(flux_time_course).asNumber(flux_units)

			flux_r, _ = stats.pearsonr(toya_fluxes, model_fluxes)

			increase_params_counts[increase_indices] += 1
			decrease_params_counts[decrease_indices] += 1
			increase_params_growth_rate[increase_indices] += growth_rate
			decrease_params_growth_rate[decrease_indices] += growth_rate
			increase_params_flux_correlation[increase_indices] += flux_r
			decrease_params_flux_correlation[decrease_indices] += flux_r

	return (increase_params_counts, decrease_params_counts,
		increase_params_growth_rate, decrease_params_growth_rate,
		increase_params_flux_correlation, decrease_params_flux_correlation)

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

		global ap
		ap = AnalysisPaths(inputDir, variant_plot=True)
		variants = ap.get_variants()
		n_variants = len(variants)

		# Load one instance of sim_data to get number of parameters and ids
		global sim_data
		global validation_data
		with open(os.path.join(inputDir, 'kb', constants.SERIALIZED_FIT1_FILENAME), 'rb') as f:
			sim_data = cPickle.load(f)
		with open(validationDataFile, 'rb') as f:
			validation_data = cPickle.load(f)

		# sim_data information
		total_params = np.sum(number_params(sim_data))
		monomer_to_gene = sim_data.moleculeGroups.frameIDGeneSymbol_Dict
		rna_to_gene = rna_mapping(sim_data, monomer_to_gene)  # TODO: create dict in sim_data from raw_data
		rna_ids = sim_data.process.transcription.rnaData['id']
		monomer_ids = sim_data.process.translation.monomerData['id']

		# IDs must match order from param_indices() from param_sensitivity.py variant
		param_ids = np.array(
			['{} RNA deg rate'.format(rna_to_gene[rna]) for rna in rna_ids]
			+ ['{} protein deg rate'.format(monomer_to_gene[monomer[:-3]]) for monomer in monomer_ids]
			+ ['{} translation eff'.format(monomer_to_gene[monomer[:-3]]) for monomer in monomer_ids]
			+ ['{} synth prob'.format(rna_to_gene[rna]) for rna in rna_ids])
		if len(param_ids) != total_params:
			raise ValueError('Number of adjusted parameters and list of ids do not match.')

		pool = Pool(processes=parallelization.plotter_cpus())
		args = zip(
			variants,
			[total_params] * n_variants,
			)
		results = pool.map(analyze_variant, args)
		pool.close()
		pool.join()
		initialize = True
		for i, result in enumerate(results):
			(_increase_params_counts,
				_decrease_params_counts,
				_increase_params_growth_rate,
				_decrease_params_growth_rate,
				_increase_params_flux_correlation,
				_decrease_params_flux_correlation) = result

			if initialize:
				initialize = False
				increase_params_counts = _increase_params_counts
				decrease_params_counts = _decrease_params_counts
				increase_params_growth_rate = _increase_params_growth_rate
				decrease_params_growth_rate = _decrease_params_growth_rate
				increase_params_flux_correlation = _increase_params_flux_correlation
				decrease_params_flux_correlation = _decrease_params_flux_correlation
			else:
				increase_params_counts += _increase_params_counts
				decrease_params_counts += _decrease_params_counts
				increase_params_growth_rate += _increase_params_growth_rate
				decrease_params_growth_rate += _decrease_params_growth_rate
				increase_params_flux_correlation += _increase_params_flux_correlation
				decrease_params_flux_correlation += _decrease_params_flux_correlation

		data = increase_params_growth_rate / increase_params_counts - decrease_params_growth_rate / decrease_params_counts
		mean = data[np.isfinite(data)].mean()
		std = data[np.isfinite(data)].std()
		z_score = (data - mean) / std

		# Plot figure
		plt.figure()

		## Plot data
		n_std = 4
		plt.bar(range(len(z_score)), np.sort(z_score))
		plt.axhline(n_std , color='r')
		plt.axhline(-n_std, color='r')

		## Format axes
		sparkline.whitePadSparklineAxis(plt.gca(), xAxis=False)
		plt.xticks([])
		plt.yticks([-n_std, 0, n_std])
		plt.xlabel('Sorted Parameters')
		plt.ylabel('Z score\nparameter effect on growth rate')

		## Save figure
		plt.tight_layout()
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close('all')

		# Print analysis summary
		upper_threshold = mean + n_std*std
		lower_threshold = mean - n_std*std
		print('Number of params above threshold: {}'.format(np.sum(data > upper_threshold)))
		print('Number of params below threshold: {}'.format(np.sum(data < lower_threshold)))

		mask = (data > upper_threshold) | (data < lower_threshold)
		print('Significant correlation between parameter and growth rate:')
		for param_id, z in sorted(zip(param_ids[mask], z_score[mask]), key=lambda v: v[1], reverse=True):
			print('\t{}: {:.2f}'.format(param_id, z))


if __name__ == "__main__":
	Plot().cli()
