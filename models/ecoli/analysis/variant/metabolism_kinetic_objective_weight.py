'''
Analyze results from metabolism_kinetic_objective_weight variant

@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 6/7/18'
'''

from __future__ import division
from __future__ import absolute_import

import argparse
import os
import re

import numpy as np
from matplotlib import pyplot as plt
import cPickle

from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from wholecell.io.tablereader import TableReader
from wholecell.analysis.analysis_tools import exportFigure
from models.ecoli.analysis import variantAnalysisPlot
from wholecell.utils import filepath

from wholecell.utils import units
from models.ecoli.processes.metabolism import COUNTS_UNITS, VOLUME_UNITS, TIME_UNITS

MODEL_FLUX_UNITS = COUNTS_UNITS / VOLUME_UNITS / TIME_UNITS
DCW_FLUX_UNITS = units.mmol / units.g / units.h

FRAC_CONC_OFF_AXIS = 0.05
FRAC_FLUX_OFF_AXIS = 0.05

OUTLIER_REACTIONS = [
	'ISOCITDEH-RXN',
	'SUCCINATE-DEHYDROGENASE-UBIQUINONE-RXN-SUC/UBIQUINONE-8//FUM/CPD-9956.31.',
	]

def get_average_values(array):
	'''
	Input:
		array (numpy array of lists) - array of data to be averaged

	Returns numpy array of floats containing the average of each list in array
	'''

	return np.array([np.mean(x) for x in array])

class Plot(variantAnalysisPlot.VariantAnalysisPlot):
	def do_plot(self, inputDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		if not os.path.isdir(inputDir):
			raise Exception, 'inputDir does not currently exist as a directory'

		ap = AnalysisPaths(inputDir, variant_plot=True)
		variants = ap.get_variants()
		n_variants = len(variants)

		if n_variants <= 1:
			print('This plot only runs for multiple variants'.format(__name__))
			return

		filepath.makedirs(plotOutDir)
		validation_data = cPickle.load(open(validationDataFile, 'rb'))
		toya_reactions = validation_data.reactionFlux.toya2010fluxes['reactionID']
		toya_fluxes = np.array([x.asNumber(DCW_FLUX_UNITS) for x in validation_data.reactionFlux.toya2010fluxes['reactionFlux']])
		outlier_filter = [False if rxn in OUTLIER_REACTIONS else True for rxn in toya_reactions]

		# Arrays to populate for plots
		lambdas = []
		growth_rates = np.empty(n_variants, dtype=object)
		conc_correlation = np.empty(n_variants, dtype=object)
		n_conc_off_axis = np.empty(n_variants, dtype=object)
		flux_correlation = np.empty(n_variants, dtype=object)
		nonzero_flux_correlation = np.empty(n_variants, dtype=object)
		n_flux_above_0 = np.empty(n_variants, dtype=object)
		n_flux_off_axis = np.empty(n_variants, dtype=object)
		correlation_coefficient = np.empty(n_variants, dtype=object)
		filtered_correlation_coefficient = np.empty(n_variants, dtype=object)
		homeostatic_objective_value = np.empty(n_variants, dtype=object)
		kinetic_objective_value = np.empty(n_variants, dtype=object)
		homeostatic_objective_std = np.empty(n_variants, dtype=object)
		kinetic_objective_std = np.empty(n_variants, dtype=object)

		# Pull information from sim data and listeners
		for i, variant in enumerate(variants):
			# Load sim_data attributes for the given variant
			sim_data = cPickle.load(open(ap.get_variant_kb(variant), 'rb'))
			cell_density = sim_data.constants.cellDensity
			n_avogadro = sim_data.constants.nAvogadro
			constrained_reactions = sim_data.process.metabolism.constrainedReactionList
			disabled_constraints = sim_data.process.metabolism.constraintsToDisable
			lambdas.append(sim_data.process.metabolism.kinetic_objective_weight)

			toya_model_fluxes = {}
			for rxn in toya_reactions:
				toya_model_fluxes[rxn] = []

			# Setup lists to store values for each cell in the variant
			growth_rates[i] = []
			conc_correlation[i] = []
			n_conc_off_axis[i] = []
			flux_correlation[i] = []
			nonzero_flux_correlation[i] = []
			n_flux_above_0[i] = []
			n_flux_off_axis[i] = []
			correlation_coefficient[i] = []
			filtered_correlation_coefficient[i] = []
			homeostatic_objective_value[i] = []
			kinetic_objective_value[i] = []
			homeostatic_objective_std[i] = []
			kinetic_objective_std[i] = []

			for sim_dir in ap.get_cells(variant=[variant]):
				sim_out_dir = os.path.join(sim_dir, 'simOut')

				# Create readers for data
				try:
					mass_reader = TableReader(os.path.join(sim_out_dir, 'Mass'))
					fba_results_reader = TableReader(os.path.join(sim_out_dir, 'FBAResults'))
					enzyme_kinetics_reader = TableReader(os.path.join(sim_out_dir, 'EnzymeKinetics'))
					bulk_reader = TableReader(os.path.join(sim_out_dir, 'BulkMolecules'))
				except Exception as e:
					print(e)
					continue

				# Mass related values
				try:
					cell_mass = units.fg * mass_reader.readColumn('cellMass')
					dry_mass = units.fg * mass_reader.readColumn('dryMass')
				except Exception as e:
					print(e)
					continue
				dcw_to_volume = cell_density * (dry_mass / cell_mass).asNumber()
				volume = cell_mass / cell_density

				# Growth rates
				# Growth rate stored in units of per second and first value will be nan
				growth_rate = mass_reader.readColumn('instantaniousGrowthRate')
				if growth_rate.size <= 1:
					continue
				growth_rates[i].append(np.mean(3600 * growth_rate[1:]))

				# Metabolite comparison
				metabolite_ids = fba_results_reader.readAttribute('homeostaticTargetMolecules')
				bulk_ids = bulk_reader.readAttribute('objectNames')

				bulk_idxs = [bulk_ids.index(id) for id in metabolite_ids]
				actual_counts = bulk_reader.readColumn('counts')[:, bulk_idxs]
				actual_conc = np.mean((1. / n_avogadro / volume * actual_counts.T).asNumber(COUNTS_UNITS / VOLUME_UNITS), axis=1)
				target_conc = np.nanmean(fba_results_reader.readColumn('targetConcentrations')[1:,:], axis=0)
				# actual_conc = np.nanmean(enzyme_kinetics_reader.readColumn('metaboliteConcentrations')[1:,:], axis=0)

				conc_correlation[i].append(np.corrcoef(actual_conc, target_conc)[0, 1])
				n_conc_off_axis[i].append(np.sum(np.abs((target_conc - actual_conc)/target_conc) > FRAC_CONC_OFF_AXIS))

				# Flux target comparison
				# Nonzero includes fluxes at 0 if target is also 0
				# Flux for first recorded step is 0
				target_fluxes = MODEL_FLUX_UNITS * enzyme_kinetics_reader.readColumn('targetFluxes').T
				actual_fluxes = MODEL_FLUX_UNITS * enzyme_kinetics_reader.readColumn('actualFluxes').T

				target_fluxes = (target_fluxes / dcw_to_volume).asNumber(DCW_FLUX_UNITS)
				actual_fluxes = (actual_fluxes / dcw_to_volume).asNumber(DCW_FLUX_UNITS)

				target_ave = np.nanmean(target_fluxes[:,1:], axis=1)
				actual_ave = np.nanmean(actual_fluxes[:,1:], axis=1)

				flux_correlation[i].append(np.corrcoef(actual_ave, target_ave)[0, 1])
				n_flux_off_axis[i].append(np.sum(np.abs((target_ave - actual_ave)/target_ave) > FRAC_FLUX_OFF_AXIS))
				mask = (actual_ave != 0)
				nonzero_flux_correlation[i].append(np.corrcoef(actual_ave[mask], target_ave[mask])[0, 1])
				n_flux_above_0[i].append(np.sum(actual_ave > 0) + np.sum((actual_ave == 0) & (target_ave == 0)))

				# Toya comparison
				# Toya units read in as mmol/g/hr
				reaction_ids = np.array(fba_results_reader.readAttribute('reactionIDs'))
				reaction_fluxes = MODEL_FLUX_UNITS * fba_results_reader.readColumn('reactionFluxes').T
				reaction_fluxes = (reaction_fluxes / dcw_to_volume).asNumber(DCW_FLUX_UNITS).T

				for toya_reaction_id in toya_reactions:
					flux_time_course = []

					for rxn in reaction_ids:
						if re.findall(toya_reaction_id, rxn):
							reverse = 1
							if re.findall('(reverse)', rxn):
								reverse = -1

							if len(flux_time_course):
								flux_time_course += reverse * reaction_fluxes[1:, np.where(reaction_ids == rxn)]
							else:
								flux_time_course = reverse * reaction_fluxes[1:, np.where(reaction_ids == rxn)]

					if len(flux_time_course):
						flux_ave = np.mean(flux_time_course)
						toya_model_fluxes[toya_reaction_id].append(flux_ave)

				# Objective values
				# Need to filter nan and inf for kinetic
				enabled_idx = [idx for idx, rxn in enumerate(constrained_reactions) if rxn not in disabled_constraints]
				kinetic_objective_values = np.abs(1 - actual_fluxes[enabled_idx, :] / target_fluxes[enabled_idx, :])
				filter_idx = ~np.isfinite(kinetic_objective_values)
				kinetic_objective_values[filter_idx] = 0
				kinetic_objective = np.sum(kinetic_objective_values, axis=0)
				kinetic_objective_value[i].append(np.mean(kinetic_objective))
				kinetic_objective_std[i].append(np.std(kinetic_objective))
				homeostatic_objective_values = np.sum(fba_results_reader.readColumn('homeostaticObjectiveValues'), axis=1)
				homeostatic_objective_value[i].append(np.mean(homeostatic_objective_values))
				homeostatic_objective_std[i].append(np.std(homeostatic_objective_values))

			ave_toya_model = np.array([np.mean(toya_model_fluxes[rxn]) for rxn in toya_reactions])
			# dist_to_diag = (ave_toya_model - toya_fluxes) / np.sqrt(2)
			# dist_std = np.std(dist_to_diag)
			# outlier_filter = np.abs(dist_to_diag) < 3 * dist_std
			correlation_coefficient[i].append(np.corrcoef(ave_toya_model, toya_fluxes)[0, 1])
			filtered_correlation_coefficient[i].append(np.corrcoef(ave_toya_model[outlier_filter], toya_fluxes[outlier_filter])[0, 1])

		n_metabolites = len(actual_conc)
		n_fluxes = len(actual_ave)
		n_sims = np.array([len(x) for x in growth_rates])
		lambdas = [np.log10(x) if x != 0 else np.nanmin(np.log10(lambdas[lambdas != 0]))-1 for x in lambdas]
		growth_rates = get_average_values(growth_rates)
		conc_correlation = get_average_values(conc_correlation)
		n_conc_off_axis = get_average_values(n_conc_off_axis)
		flux_correlation = get_average_values(flux_correlation)
		nonzero_flux_correlation = get_average_values(nonzero_flux_correlation)
		n_flux_above_0 = get_average_values(n_flux_above_0)
		n_flux_off_axis = get_average_values(n_flux_off_axis)
		correlation_coefficient = get_average_values(correlation_coefficient)
		filtered_correlation_coefficient = get_average_values(filtered_correlation_coefficient)
		homeostatic_objective_value = get_average_values(homeostatic_objective_value)
		kinetic_objective_value = get_average_values(kinetic_objective_value)
		homeostatic_objective_std = get_average_values(homeostatic_objective_std)
		kinetic_objective_std = get_average_values(kinetic_objective_std)

		plt.figure(figsize = (8.5, 22))
		subplots = 8

		# Growth rates
		plt.subplot(subplots, 1, 1)
		plt.style.use('seaborn-deep')
		plt.bar(lambdas, growth_rates - growth_rates[0], align='center')
		plt.ylim([-1, 1])
		plt.title('Growth rate deviation from no kinetics')
		plt.ylabel('Deviation (1/hr)')

		# Metabolite comparisons
		plt.subplot(subplots, 1, 2)
		plt.bar(lambdas, conc_correlation, align='center')
		plt.ylim([0, 1])
		plt.ylabel('PCC')
		plt.title('Concentration correlation')

		plt.subplot(subplots, 1, 3)
		plt.bar(lambdas, n_conc_off_axis / n_metabolites, align='center')
		plt.ylim([0, 1])
		plt.ylabel('Fraction of concentrations')
		plt.title('Concentrations off axis (>{:.0f}%)'.format(FRAC_CONC_OFF_AXIS*100))

		# Flux target comparisons
		plt.subplot(subplots, 1, 4)
		plt.bar(lambdas, nonzero_flux_correlation, align='center', color='r')
		plt.bar(lambdas, flux_correlation, align='center')
		plt.ylim([0, 1])
		plt.ylabel('PCC')
		plt.title('Flux correlation')

		plt.subplot(subplots, 1, 5)
		plt.bar(lambdas, n_flux_above_0 / n_fluxes, align='center')
		plt.ylim([0, 1])
		plt.ylabel('Fraction of fluxes')
		plt.title('Flux above 0')

		plt.subplot(subplots, 1, 6)
		plt.bar(lambdas, n_flux_off_axis / n_fluxes, align='center')
		plt.ylim([0, 1])
		plt.ylabel('Fraction of fluxes')
		plt.title('Fluxes off axis (>{:.0f}%)'.format(FRAC_FLUX_OFF_AXIS*100))

		# Toya comparison
		plt.subplot(subplots, 1, 7)
		plt.bar(lambdas, filtered_correlation_coefficient, align='center', color='r')
		plt.bar(lambdas, correlation_coefficient, align='center')
		plt.ylim([0, 1])
		plt.ylabel('PCC')
		plt.title('Central carbon flux correlation')

		# Viable sims
		plt.subplot(subplots, 1, 8)
		plt.bar(lambdas, n_sims, align='center')
		plt.ylabel('Number of sims')
		plt.title('Complete sims')

		plt.xlabel('lambda')

		exportFigure(plt, plotOutDir, plotOutFileName, metadata)

		plt.figure()
		ax = plt.gca()
		ax.set_xscale("log", nonposx='clip')
		ax.set_yscale("log", nonposy='clip')
		plt.errorbar(homeostatic_objective_value, kinetic_objective_value, xerr=homeostatic_objective_std, yerr=kinetic_objective_std, fmt='o')
		for i in range(len(lambdas)):
			plt.text(homeostatic_objective_value[i], 0.8*kinetic_objective_value[i], i, horizontalalignment='center', verticalalignment='center')
		plt.xlabel('Homeostatic Objective Value')
		plt.ylabel('Kinetics Objective Value')
		exportFigure(plt, plotOutDir, '{}_obj'.format(plotOutFileName), metadata)

		plt.close('all')


if __name__ == "__main__":
	Plot().cli()
