#!/usr/bin/env python
'''
Analyze results from metabolism_lambda_weight variant

@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 6/7/18'
'''

from __future__ import division

import argparse
import os
import re

import numpy as np
from matplotlib import pyplot as plt
import cPickle

from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from wholecell.io.tablereader import TableReader
import wholecell.utils.constants
from wholecell.utils import filepath
from wholecell.analysis.analysis_tools import exportFigure

from wholecell.utils import units
from models.ecoli.processes.metabolism import COUNTS_UNITS, VOLUME_UNITS, TIME_UNITS

MODEL_FLUX_UNITS = COUNTS_UNITS / VOLUME_UNITS / TIME_UNITS
DCW_FLUX_UNITS = units.mmol / units.g / units.h

FRAC_CONC_OFF_AXIS = 0.05
FRAC_FLUX_OFF_AXIS = 0.05

def main(inputDir, plotOutDir, plotOutFileName, validationDataFile, metadata = None):
	if not os.path.isdir(inputDir):
		raise Exception, 'inputDir does not currently exist as a directory'

	ap = AnalysisPaths(inputDir, variant_plot=True)
	variants = sorted(ap._path_data['variant'].tolist())

	if len(variants) <= 1:
		return

	all_cells = sorted(ap.get_cells(seed=[0], generation=[0]))
	variants = ap.get_variants()

	filepath.makedirs(plotOutDir)
	validation_data = cPickle.load(open(validationDataFile, 'rb'))

	# Arrays to populate for plots
	lambdas = []
	growth_rates = []
	conc_correlation = []
	n_conc_off_axis = []
	flux_correlation = []
	nonzero_flux_correlation = []
	n_flux_above_0 = []
	n_flux_off_axis = []
	correlation_coefficient = []
	homeostatic_objective_value = []
	kinetic_objective_value = []
	homeostatic_objective_std = []
	kinetic_objective_std = []

	# Pull information from sim data and listeners
	for variant, sim_dir in zip(variants, all_cells):
		sim_out_dir = os.path.join(sim_dir, 'simOut')
		sim_data = cPickle.load(open(ap.get_variant_kb(variant), 'rb'))
		cell_density = sim_data.constants.cellDensity
		n_avogadro = sim_data.constants.nAvogadro
		constrained_reactions = sim_data.process.metabolism.constrainedReactionList
		disabled_constraints = sim_data.process.metabolism.constraintsToDisable

		mass_reader = TableReader(os.path.join(sim_out_dir, 'Mass'))
		fba_results_reader = TableReader(os.path.join(sim_out_dir, 'FBAResults'))
		enzyme_kinetics_reader = TableReader(os.path.join(sim_out_dir, 'EnzymeKinetics'))
		bulk_reader = TableReader(os.path.join(sim_out_dir, 'BulkMolecules'))

		# Mass related values
		cell_mass = units.fg * mass_reader.readColumn('cellMass')
		dry_mass = units.fg * mass_reader.readColumn('dryMass')
		dcw_to_volume = cell_density * (dry_mass / cell_mass).asNumber()
		avg_dcw_to_volume = np.mean(dcw_to_volume)
		volume = cell_mass / cell_density

		# lambda values
		lambdas.append(sim_data.constants.metabolismKineticObjectiveWeight)

		# Growth rates
		# Growth rate stored in units of per second and first value will be nan
		growth_rate = 3600 * mass_reader.readColumn('instantaniousGrowthRate')[1:]
		growth_rates.append(np.mean(growth_rate))

		# Metabolite comparison
		metabolite_ids = fba_results_reader.readAttribute('homeostaticTargetMolecules')
		bulk_ids = bulk_reader.readAttribute('objectNames')

		bulk_idxs = [bulk_ids.index(id) for id in metabolite_ids]
		actual_counts = bulk_reader.readColumn('counts')[:, bulk_idxs]
		actual_conc = np.mean((1. / n_avogadro / volume * actual_counts.T).asNumber(COUNTS_UNITS / VOLUME_UNITS), axis=1)
		target_conc = np.nanmean(fba_results_reader.readColumn('targetConcentrations')[1:,:], axis=0)
		# actual_conc = np.nanmean(enzyme_kinetics_reader.readColumn('metaboliteConcentrations')[1:,:], axis=0)

		conc_correlation.append(np.corrcoef(actual_conc, target_conc)[0, 1])
		n_conc_off_axis.append(np.sum(np.abs((target_conc - actual_conc)/target_conc) > FRAC_CONC_OFF_AXIS))

		# Flux target comparison
		# Nonzero includes fluxes at 0 if target is also 0
		target_fluxes = MODEL_FLUX_UNITS * enzyme_kinetics_reader.readColumn('targetFluxes').T
		actual_fluxes = MODEL_FLUX_UNITS * enzyme_kinetics_reader.readColumn('actualFluxes').T

		target_fluxes = (target_fluxes / dcw_to_volume).asNumber(DCW_FLUX_UNITS)
		actual_fluxes = (actual_fluxes / dcw_to_volume).asNumber(DCW_FLUX_UNITS)

		target_ave = np.nanmean(target_fluxes[:,1:], axis=1)
		actual_ave = np.nanmean(actual_fluxes[:,1:], axis=1)

		flux_correlation.append(np.corrcoef(actual_ave, target_ave)[0, 1])
		n_flux_off_axis.append(np.sum(np.abs((target_ave - actual_ave)/target_ave) > FRAC_FLUX_OFF_AXIS))
		mask = (actual_ave != 0)
		nonzero_flux_correlation.append(np.corrcoef(actual_ave[mask], target_ave[mask])[0, 1])
		n_flux_above_0.append(np.sum(actual_ave > 0) + np.sum((actual_ave == 0) & (target_ave == 0)))

		# Toya comparison
		# Toya units read in as mmol/g/hr
		reaction_ids = np.array(fba_results_reader.readAttribute('reactionIDs'))
		reaction_fluxes = fba_results_reader.readColumn('reactionFluxes')

		toya_reactions = validation_data.reactionFlux.toya2010fluxes['reactionID']
		toya_fluxes = np.array([(avg_dcw_to_volume * x).asNumber(MODEL_FLUX_UNITS) for x in validation_data.reactionFlux.toya2010fluxes['reactionFlux']])
		toya_stdev = np.array([(avg_dcw_to_volume * x).asNumber(MODEL_FLUX_UNITS) for x in validation_data.reactionFlux.toya2010fluxes['reactionFluxStdev']])
		toya_fluxes_dict = dict(zip(toya_reactions, toya_fluxes))
		toya_stdev_dict = dict(zip(toya_reactions, toya_stdev))

		toya_vs_reaction_ave = []
		toya_order = []
		for toya_reaction_id, toya_flux in toya_fluxes_dict.iteritems():
			flux_time_course = []

			for rxn in reaction_ids:
				if re.findall(toya_reaction_id, rxn):
					reverse = 1
					if re.findall('(reverse)', rxn):
						reverse = -1

					if len(flux_time_course):
						flux_time_course += reverse * reaction_fluxes[:, np.where(reaction_ids == rxn)]
					else:
						flux_time_course = reverse * reaction_fluxes[:, np.where(reaction_ids == rxn)]

			if len(flux_time_course):
				# flip sign if negative
				adjustment = 1
				if toya_flux < 0:
					adjustment = -1

				flux_ave = np.mean(flux_time_course)
				flux_stdev = np.std(flux_time_course)
				toya_vs_reaction_ave.append((adjustment*flux_ave, adjustment*toya_flux, flux_stdev, toya_stdev_dict[toya_reaction_id]))
				toya_order.append(toya_reaction_id)

		toya_vs_reaction_ave = np.array(toya_vs_reaction_ave)
		correlation_coefficient.append(np.corrcoef(toya_vs_reaction_ave[:, 0], toya_vs_reaction_ave[:, 1])[0, 1])

		# Objective values
		# Need to filter nan and inf for kinetic
		enabled_idx = [i for i, rxn in enumerate(constrained_reactions) if rxn not in disabled_constraints]
		kinetic_objective_values = np.abs(1 - actual_fluxes[enabled_idx, :] / target_fluxes[enabled_idx, :])
		filter_idx = ~np.isfinite(kinetic_objective_values)
		kinetic_objective_values[filter_idx] = 0
		kinetic_objective = np.sum(kinetic_objective_values, axis=0)
		kinetic_objective_value.append(np.mean(kinetic_objective))
		kinetic_objective_std.append(np.std(kinetic_objective))
		homeostatic_objective_values = np.sum(fba_results_reader.readColumn('homeostaticObjectiveValues'), axis=1)
		homeostatic_objective_value.append(np.mean(homeostatic_objective_values))
		homeostatic_objective_std.append(np.std(homeostatic_objective_values))

	n_metabolites = len(actual_conc)
	n_fluxes = len(actual_ave)
	lambdas = [np.log10(x) if x != 0 else np.nanmin(np.log10(lambdas[lambdas != 0]))-1 for x in lambdas]

	plt.figure(figsize = (8.5, 22))
	subplots = 7

	# Growth rates
	plt.subplot(subplots, 1, 1)
	plt.style.use('seaborn-deep')
	plt.bar(lambdas, growth_rates, align='center')
	plt.ylim([0, 2])
	plt.title('Growth rate')
	plt.ylabel('Growth rate (1/hr)')

	# Metabolite comparisons
	plt.subplot(subplots, 1, 2)
	plt.bar(lambdas, conc_correlation, align='center')
	plt.ylim([0, 1])
	plt.ylabel('PCC')
	plt.title('Concentration correlation')

	plt.subplot(subplots, 1, 3)
	plt.bar(lambdas, np.array(n_conc_off_axis) / n_metabolites, align='center')
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
	plt.bar(lambdas, np.array(n_flux_above_0) / n_fluxes, align='center')
	plt.ylim([0, 1])
	plt.ylabel('Fraction of fluxes')
	plt.title('Flux above 0')

	plt.subplot(subplots, 1, 6)
	plt.bar(lambdas, np.array(n_flux_off_axis) / n_fluxes, align='center')
	plt.ylim([0, 1])
	plt.ylabel('Fraction of fluxes')
	plt.title('Fluxes off axis (>{:.0f}%)'.format(FRAC_FLUX_OFF_AXIS*100))

	# Toya comparison
	plt.subplot(subplots, 1, 7)
	plt.bar(lambdas, correlation_coefficient, align='center')
	plt.ylim([0, 1])
	plt.ylabel('PCC')
	plt.title('Central carbon flux correlation')

	plt.xlabel('lambda')

	exportFigure(plt, plotOutDir, plotOutFileName, metadata)

	plt.figure()
	ax = plt.gca()
	ax.set_xscale("log", nonposx='clip')
	ax.set_yscale("log", nonposy='clip')
	plt.errorbar(homeostatic_objective_value, kinetic_objective_value, xerr=homeostatic_objective_std, yerr=kinetic_objective_std, fmt='o')
	plt.xlabel('Homeostatic Objective Value')
	plt.ylabel('Kinetics Objective Value')
	exportFigure(plt, plotOutDir, '{}_obj'.format(plotOutFileName), metadata)

	plt.close('all')


if __name__ == '__main__':
	defaultSimDataFile = os.path.join(
			wholecell.utils.constants.SERIALIZED_KB_DIR,
			wholecell.utils.constants.SERIALIZED_KB_MOST_FIT_FILENAME
			)

	parser = argparse.ArgumentParser()
	parser.add_argument('simOutDir', help = 'Directory containing simulation output', type = str)
	parser.add_argument('plotOutDir', help = 'Directory containing plot output (will get created if necessary)', type = str)
	parser.add_argument('plotOutFileName', help = 'File name to produce', type = str)
	parser.add_argument('--simDataFile', help = 'KB file name', type = str, default = defaultSimDataFile)

	args = parser.parse_args().__dict__

	main(args['simOutDir'], args['plotOutDir'], args['plotOutFileName'], args['simDataFile'])
