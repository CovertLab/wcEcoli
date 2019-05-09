from __future__ import absolute_import, division, print_function

import os
import csv
import json
import matplotlib.pyplot as plt
import numpy as np

import scipy.constants as constants
from reconstruction.spreadsheets import JsonReader
from itertools import ifilter

import environment.kinetic_rate_laws.kinetic_rate_laws as rate_laws


LOOKUP_DIR = os.path.join('environment', 'condition', 'look_up_tables')
TRANSPORT_REACTIONS_FILE = os.path.join(LOOKUP_DIR, 'transport_reactions.tsv')
CONC_LOOKUP_MINIMAL = os.path.join(LOOKUP_DIR, 'avg_concentrations', 'minimal.tsv')

# external molecules file to make a map between
KINETIC_PARAMETERS_FILE = os.path.join('environment', 'kinetic_rate_laws', 'parameters', 'glt.json')

# set output directory and files
OUTPUT_DIR = os.path.join('environment', 'kinetic_rate_laws', 'out')
OUTPUT_PARAM_TEMPLATE = os.path.join(OUTPUT_DIR, 'parameter_template.json')


ANALYZE_RATE_LAWS = True
SAVE_RATE_LAWS_CONFIG = False


def test_rate_laws():
	'''
	Tests the rate law configuration defined in KINETIC_PARAMETERS_FILE, save analysis plot

	'''

	# Make dict of transport reactions
	all_transport_reactions = {}
	with open(TRANSPORT_REACTIONS_FILE, 'rU') as csvfile:
		reader = JsonReader(
			ifilter(lambda x: x.lstrip()[0] != "#", csvfile),  # Strip comments
			dialect=csv.excel_tab)
		for row in reader:
			reaction_id = row['reaction id']
			stoichiometry = row['stoichiometry']
			reversible = row['is reversible']
			catalyzed = row['catalyzed by']
			all_transport_reactions[reaction_id] = {
				'stoichiometry': stoichiometry,
				'is reversible': reversible,
				'catalyzed by': catalyzed,
			}

	wcm_concs_minimal = {}
	with open(CONC_LOOKUP_MINIMAL, 'rU') as csvfile:
		reader = JsonReader(
			ifilter(lambda x: x.lstrip()[0] != "#", csvfile),  # Strip comments
			dialect=csv.excel_tab)
		for row in reader:
			molecule_id = row['molecule id']
			avg_conc = row['average concentration mmol/L']
			wcm_concs_minimal[molecule_id] = avg_conc

	# load dict of saved parameters
	with open(KINETIC_PARAMETERS_FILE, 'r') as fp:
		kinetic_parameters = json.load(fp)

	# make a dict of reactions that will be configured with the parameters
	make_reaction_ids = kinetic_parameters.keys()
	make_reactions = {
		reaction_id: specs
		for reaction_id, specs in all_transport_reactions.iteritems()
		if reaction_id in make_reaction_ids}

	# Make the kinetic model
	kinetic_rate_laws = rate_laws.KineticFluxModel(make_reactions, kinetic_parameters)

	# Get list of molecule_ids used by kinetic rate laws
	molecule_ids = kinetic_rate_laws.molecule_ids

	# get concentrations from wcm (mmol/L)
	concentrations = {mol_id: wcm_concs_minimal[mol_id] for mol_id in molecule_ids}

	# get reaction fluxes
	reaction_fluxes = kinetic_rate_laws.get_fluxes(concentrations)

	# run analyses and save output
	analyze_rate_laws(kinetic_rate_laws, concentrations)


def analyze_rate_laws(kinetic_rate_laws, baseline_concentrations):
	'''
	Args:
		kinetic_rate_laws (object): a configured kinetic_rate_law object
		baseline_concentrations (dict): concentrations for all molecules required for the rate laws

	Function:
		Runs an analysis of all rate laws in kinetic_rate_laws and saves the output

	'''

	test_transporter = True
	test_cofactor = True
	test_competitor = True

	reactions = kinetic_rate_laws.reactions
	kinetic_parameters = kinetic_rate_laws.kinetic_parameters
	rate_law_configuration = kinetic_rate_laws.rate_law_configuration

	## Plot analysis

	columns = 1 + sum([test_transporter, test_cofactor, test_competitor])
	n_samples = 100
	n_samples_shown = 10
	n_rxns = len(reactions)
	rows = 2*n_rxns + 2  # extra row for each reaction header

	cmap = plt.cm.get_cmap('Spectral')
	colors = [cmap(float(idx) / n_samples_shown) for idx in range(n_samples_shown)]

	plt.figure(figsize=(6*columns, 3*rows))
	plot_number = 1
	row_number = 0

	for reaction_id, specs in reactions.iteritems():
		transporters = specs.get('catalyzed by')
		stoich = specs.get('stoichiometry')
		parameters = kinetic_parameters.get(reaction_id)

		reactants = [mol for mol, coeff in stoich.iteritems() if coeff < 0]
		products = [mol for mol, coeff in stoich.iteritems() if coeff > 0]

		plt.subplot(rows, columns, plot_number)
		plt.text(0.02, 0.6, 'reaction: ' + reaction_id, weight='bold')
		plt.text(0.02, 0.45, 'reactants: %s' % reactants)
		plt.text(0.02, 0.3, 'products: %s' % products)
		plt.text(0.02, 0.15, 'transporters: %s' % transporters)
		plt.text(0.02, 0.0, 'parameters: %s' % parameters)
		plt.axis('off')
		plot_number += columns
		row_number += 1

		# test michaelis menten by sampling substrate concentrations
		for transporter in transporters:

			# TODO select the molecule more smartly
			for reactant in reactants:
				if parameters[transporter][reactant] is not None:
					a1 = reactant

			# get cofactor
			b1 = None
			if len(reactants) > 1:
				for reactant in reactants:
					if parameters[transporter][reactant] is not None and reactant is not a1:
						b1 = reactant

			concentrations = baseline_concentrations.copy()
			conc_samples = np.logspace(-9, 0, num=n_samples, endpoint=True, base=10)

			flux_values = np.empty_like(conc_samples)
			for idx, conc in enumerate(conc_samples):
				concentrations[a1] = conc
				reaction_fluxes = kinetic_rate_laws.get_fluxes(concentrations)
				flux_values[idx] = reaction_fluxes[reaction_id]

			# plot M-M curve for this reaction
			plt.subplot(rows, columns, plot_number)
			plt.plot(conc_samples, flux_values)

			# plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
			plt.xscale('log')
			plt.xlabel(a1 + ' concentration (M)')
			plt.ylabel('flux (M/s)')
			plt.title('transporter: %s' % transporter)

			plot_number += 1

			if test_transporter:
				concentrations = baseline_concentrations.copy()
				conc_samples = np.logspace(-8, 1, num=n_samples, endpoint=True, base=10)
				transporter_concs = np.logspace(-4, 1, num=n_samples_shown, endpoint=True, base=10)

				plt.subplot(rows, columns, plot_number)
				for index, transporter_conc in enumerate(transporter_concs):
					concentrations[transporter] = transporter_conc

					flux_values = np.empty_like(conc_samples)
					for idx, conc in enumerate(conc_samples):

						concentrations[a1] = conc
						reaction_fluxes = kinetic_rate_laws.get_fluxes(concentrations)
						flux_values[idx] = reaction_fluxes[reaction_id]

					# plot M-M curve for this reaction
					plt.plot(conc_samples, flux_values,
						color = colors[index],
						label = ('conc = %.2e' % (transporter_conc)),
						)

				plt.legend(loc='center left', title=transporter, bbox_to_anchor=(1.15, 0.5))
				plt.xscale('log')
				plt.xlabel(a1 + ' concentration (M)')
				plt.ylabel('flux (M/s)')
				plt.title('test transporter')

				plot_number += 1

			if test_cofactor:

				concentrations = baseline_concentrations.copy()
				conc_samples = np.logspace(-8, 1, num=n_samples, endpoint=True, base=10)
				cofactor_concs = np.logspace(-8, 1, num=n_samples_shown, endpoint=True, base=10)

				if b1 is not None:
					plt.subplot(rows, columns, plot_number)
					for index, cofactor_conc in enumerate(cofactor_concs):
						concentrations[b1] = cofactor_conc

						flux_values = np.empty_like(conc_samples)
						for idx, conc in enumerate(conc_samples):

							concentrations[a1] = conc
							reaction_fluxes = kinetic_rate_laws.get_fluxes(concentrations)
							flux_values[idx] = reaction_fluxes[reaction_id]

						# plot M-M curve for this reaction
						plt.plot(conc_samples, flux_values,
										color = colors[index],
										label = ('conc = %.2e' % (cofactor_conc)),
										)

					plt.legend(loc='center left', title=b1, bbox_to_anchor=(1.15, 0.5))
					plt.xscale('log')
					plt.xlabel(a1 + ' concentration (M)')
					plt.ylabel('flux (M/s)')
					plt.title('test cofactor')

				plot_number += 1

			if test_competitor:
				# get competitor
				rxns_transporter = rate_law_configuration[transporter]['reaction_cofactors'].keys()
				competing_rxns = [trpr for trpr in rxns_transporter if trpr not in reaction_id]
				competitor = None
				for rx in competing_rxns:
					competitor_candidates = rate_law_configuration[transporter]['reaction_cofactors'][rx]
					competitors = [mol for mol in competitor_candidates[0] if mol != a1 and mol != b1]
					if competitors:
						competitor = competitors[0]



				if competitor is not None:
					concentrations = baseline_concentrations.copy()
					conc_samples = np.logspace(-8, 1, num=n_samples, endpoint=True, base=10)
					competitor_concs = np.logspace(-8, 1, num=n_samples_shown, endpoint=True, base=10)

					plt.subplot(rows, columns, plot_number)
					for index, competitor_conc in enumerate(competitor_concs):
						concentrations[competitor] = competitor_conc

						flux_values = np.empty_like(conc_samples)
						for idx, conc in enumerate(conc_samples):

							concentrations[a1] = conc
							reaction_fluxes = kinetic_rate_laws.get_fluxes(concentrations)
							flux_values[idx] = reaction_fluxes[reaction_id]

						# plot M-M curve for this reaction
						plt.plot(conc_samples, flux_values,
												color = colors[index],
												label = ('conc = %.2e' % (competitor_conc)),
												)

					plt.legend(loc='center left', title=competitor, bbox_to_anchor=(1.15, 0.5))
					plt.xscale('log')
					plt.xlabel(a1 + ' concentration (M)')
					plt.ylabel('flux (M/s)')
					plt.title('test competitor')

				plot_number += 1

			row_number += 1

		plot_number = row_number * columns + 1

	plt.subplots_adjust(hspace=0.8, wspace=1.1)

	if not os.path.exists(OUTPUT_DIR):
		os.mkdir(OUTPUT_DIR)
	fig_name = ('rate_law_analysis')
	plt.savefig(os.path.join(OUTPUT_DIR, fig_name), bbox_inches='tight')

	print('rate law analysis plot saved')


def save_rate_law_configuration_template(reactions):
	'''
	Create an empty parameter file for convenience kinetics rate laws structured by reactions

	Args:
		reactions (dict): a dict of all reactions that are to be configured
	'''

	# get the rate law configuration for the set of reactions
	rate_law_configuration = rate_laws.make_configuration(reactions)

	# make a parameter template
	parameter_template = rate_laws.get_parameter_template(reactions, rate_law_configuration)

	with open(OUTPUT_PARAM_TEMPLATE, 'w') as fp:
		json.dump(parameter_template, fp, sort_keys=True, indent=2)

	print('rate law parameter template saved')


if SAVE_RATE_LAWS_CONFIG:
	'''
	make reaction dictionary for a set of exchange molecules and pass to save_rate_law_configuration_template
	'''

	# set up reactions
	exchange_molecules = ['GLT[p]']

	# Make dict of transport reactions
	all_reactions = {}
	with open(TRANSPORT_REACTIONS_FILE, 'rU') as csvfile:
		reader = JsonReader(
			ifilter(lambda x: x.lstrip()[0] != "#", csvfile),  # Strip comments
			dialect=csv.excel_tab)
		for row in reader:
			reaction_id = row['reaction id']
			stoichiometry = row['stoichiometry']
			reversible = row['is reversible']
			catalyzed = row['catalyzed by']
			all_reactions[reaction_id] = {
				'stoichiometry': stoichiometry,
				'is reversible': reversible,
				'catalyzed by': catalyzed,
			}

	# get a list of all reactions with exchange_molecules
	reactions_list = rate_laws.get_reactions_from_exchange(all_reactions, exchange_molecules)

	# make a dict of the given reactions using specs from all_reactions
	reactions = {reaction_id: all_reactions[reaction_id] for reaction_id in reactions_list}

	save_rate_law_configuration_template(reactions)


if ANALYZE_RATE_LAWS:
	'''test the rate law specified in KINETIC_PARAMETERS_FILE'''

	# Run test and save output
	test_rate_laws()
