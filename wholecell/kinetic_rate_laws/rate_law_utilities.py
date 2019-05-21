'''
Rate law utilities

This collection of functions can assist in constructing and analyzing kinetic rate laws
that are generated in wholecell.kinetic_rate_laws.kinetic_rate_laws.

Functions include:
	- load_reactions(): returns a dict of all_reactions, with location tags added to the enzymes
	- get_reactions_from_exchange: provided a dict of reactions and exchange molecules, this returns a list of all reactions for those exchange molecules
	- get_molecules_from_reactions: given a dict of reactions, returns all the relevant molecules -- substrates and enzymes

The RateLawUtilities module can be called with:
> python -m wholecell.kinetic_rate_laws.rate_law_utilities

'''

from __future__ import absolute_import, division, print_function

import os
import csv
import json
import matplotlib.pyplot as plt
import numpy as np
import argparse

from reconstruction.spreadsheets import JsonReader
from itertools import ifilter

from environment.condition.look_up_tables.look_up import LookUp
import wholecell.kinetic_rate_laws.kinetic_rate_laws as rate_laws

TSV_DIALECT = csv.excel_tab

REACTIONS_FILE = os.path.join("reconstruction", "ecoli", "flat", "reactions.tsv")
PROTEINS_FILE = os.path.join("reconstruction", "ecoli", "flat", "proteins.tsv")
COMPLEXATION_FILE = os.path.join("reconstruction", "ecoli", "flat", "complexationReactions.tsv")
KINETIC_PARAMETERS_PATH = os.path.join('wholecell', 'kinetic_rate_laws', 'parameters')
OUTPUT_DIR = os.path.join('wholecell', 'kinetic_rate_laws', 'out')


def analyze_rate_laws(kinetic_rate_laws, baseline_concentrations):
	'''
	Args:
		kinetic_rate_laws (object): a configured kinetic_rate_law object
		baseline_concentrations (dict): concentrations for all molecules required for the rate laws

	Function:
		Runs an analysis of all rate laws in kinetic_rate_laws and saves the output in a plot

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

			# TODO select the example molecule more smartly
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
	plt.savefig(os.path.join(OUTPUT_DIR, 'rate_law_analysis'), bbox_inches='tight')

	print('rate law analysis plot saved')

def load_reactions():
	'''
	Load all reactions, including the locations of enzymes into each reaction's 'catalyzed by' key

	Returns:
		all_reactions(dict)
	'''

	# get protein locations
	proteins_locations = {}
	with open(PROTEINS_FILE, 'rU') as tsvfile:
		reader = JsonReader(
			ifilter(lambda x: x.lstrip()[0] != "#", tsvfile), # Strip comments
			dialect = TSV_DIALECT)
		for row in reader:
			molecule_id = row["id"]
			location = row["location"]
			molecule_loc = ['{}[{}]'.format(molecule_id, loc) for loc in location]
			proteins_locations[molecule_id] = molecule_loc

	# get complex locations
	with open(COMPLEXATION_FILE, 'rU') as tsvfile:
		reader = JsonReader(
			ifilter(lambda x: x.lstrip()[0] != "#", tsvfile), # Strip comments
			dialect = TSV_DIALECT)
		for row in reader:
			stoichiometry = row["stoichiometry"]
			for stoich in stoichiometry:
				molecule_id =  stoich['molecule']
				location = stoich["location"]
				molecule_loc = ['{}[{}]'.format(molecule_id, loc) for loc in location]
				proteins_locations[molecule_id] = molecule_loc

	# make dict of all reactions
	all_reactions = {}
	with open(REACTIONS_FILE, 'rU') as tsvfile:
		reader = JsonReader(
			ifilter(lambda x: x.lstrip()[0] != "#", tsvfile), # Strip comments
			dialect = TSV_DIALECT)
		for row in reader:
			reaction_id = row["reaction id"]
			stoichiometry = row["stoichiometry"]
			reversible = row["is reversible"]
			enzymes = row["catalyzed by"]

			# get location
			enzymes_loc = []
			for enzyme in enzymes:
				if enzyme in proteins_locations.keys():
					enzymes_loc.extend(proteins_locations[enzyme])

			all_reactions[reaction_id] = {
				"stoichiometry": stoichiometry,
				"is reversible": reversible,
				"catalyzed by": enzymes_loc,
			}

	return all_reactions

def get_reactions_from_exchange(all_reactions, include_exchanges):
	'''
	Args:
		all_reactions (dict): all reactions with stoichiometry, reversibility, enzymes
		include_exchanges (list): molecules whose reactions are of interest

	Returns:
		include_reactions (list): all the reactions for molecules listed in include_exchanges

	'''
	include_reactions = []
	for reaction_id, specs in all_reactions.iteritems():
		reaction_molecules = specs['stoichiometry'].keys()
		for exchange in include_exchanges:
			if exchange in reaction_molecules:
				include_reactions.append(reaction_id)
	return include_reactions

def get_molecules_from_reactions(reactions):
	'''
	Inputs:
		   reaction_ids - a list of all reaction ids that will be used by transport
	Returns:
		   self.molecule_ids - a list of all molecules used by these reactions
	'''
	molecule_ids = []
	for reaction_id, specs in reactions.iteritems():
		stoichiometry = specs['stoichiometry']
		substrates = stoichiometry.keys()
		enzymes = specs['catalyzed by']
		# Add all relevant molecules_ids
		molecule_ids.extend(substrates)
		molecule_ids.extend(enzymes)
	return list(set(molecule_ids))


class RateLawUtilities(object):

	def __init__(self):

		parser = argparse.ArgumentParser(description='analyze rate laws')
		parser = self.add_arguments(parser)
		args = parser.parse_args()

		# load all reactions
		self.all_reactions = load_reactions()

		# load dict of saved parameters
		parameter_file = args.path
		with open(parameter_file, 'r') as fp:
			kinetic_parameters = json.load(fp)

		# make a dict of reactions that will be configured with the parameters
		make_reaction_ids = kinetic_parameters.keys()
		make_reactions = {
			reaction_id: specs
			for reaction_id, specs in self.all_reactions.iteritems()
			if reaction_id in make_reaction_ids}

		# Make the kinetic model
		self.kinetic_rate_laws = rate_laws.KineticFluxModel(make_reactions, kinetic_parameters)

		# Get list of molecule_ids used by kinetic rate laws
		self.molecule_ids = self.kinetic_rate_laws.molecule_ids

		# make look up object and get saved concentrations from wcm (mmol/L)
		self.look_up = LookUp()
		self.concentrations = self.look_up.look_up('average', args.media, self.molecule_ids)

		if args.analyze:
			self.run_analysis()

		if args.template:
			reactions_list = args.template.split(',')
			self.template_from_reactions(reactions_list)

		if args.exchange_template:
			molecule_list = args.exchange_template.split(',')
			self.template_from_exchange(molecule_list)

	def add_arguments(self, parser):

		parser.add_argument(
			'-a', '--analyze',
			action='store_true',
			default=False,
			help='run analysis on parameter files specified by path')

		parser.add_argument(
			'-t', '--template',
			type=str,
			default='',
			help='A list of reactions for making an empty parameter template, formatted as "reaction_id_1, reaction_id_2"')

		parser.add_argument(
			'-ex', '--exchange_template',
			type=str,
			default='',
			help='A list of exchange molecules for making an empty parameter template, formatted as "molecule_id_1, molecule_id_2"')

		parser.add_argument(
			'-m', '--media',
			type=str,
			default='minimal',
			help='The environment media')

		parser.add_argument(
			'--path',
			type=str,
			default=os.path.join(KINETIC_PARAMETERS_PATH, 'example_parameters.json') ,
			help='the path to the parameter file to be analyzed. parameters available in {}'.format(KINETIC_PARAMETERS_PATH))

		return parser

	def run_analysis(self):

		# run analyses and save output
		analyze_rate_laws(self.kinetic_rate_laws, self.concentrations)

	def template_from_exchange(self, exchange_molecules):
		'''
		saves a rate law parameter template for the list of exchange molecules

		'''

		# get a list of all reactions with exchange_molecules
		reactions_list = get_reactions_from_exchange(self.all_reactions, exchange_molecules)

		self.template_from_reactions(reactions_list)

	def template_from_reactions(self, reactions_list):
		'''
		saves a rate law parameter template for the list of reactions

		'''
		# make a dict of the given reactions using specs from all_reactions
		reactions = {reaction_id: self.all_reactions[reaction_id] for reaction_id in reactions_list}

		# make a parameter template
		parameter_template = self.get_parameter_template(reactions)

		output_name = os.path.join(OUTPUT_DIR, 'parameter_template.json')
		with open(output_name, 'w') as fp:
			json.dump(parameter_template, fp, sort_keys=True, indent=2)

		print('rate law parameter template saved')

	def get_parameter_template(self, reactions):
		'''
		Given a list of reactions, return a template for required parameters

		Args:
			reactions (dict): a reaction network, with
			 {reaction_id: {'catalyzed by': (list), 'is reversible': (bool), 'stoichiometry': (dict)}}

		Returns:
			parameter_template (dict): a template for all parameters required by this rate_law_configuration,
				filled with values of 0.0.

		'''

		rate_law_configuration = rate_laws.make_configuration(reactions)

		parameter_template = {}
		for enzyme_id, configuration in rate_law_configuration.iteritems():
			reaction_cofactors = configuration['reaction_cofactors']
			partition = configuration['partition']

			for reaction_id, cofactors in reaction_cofactors.iteritems():

				# check if reaction is already in the template
				if reaction_id not in parameter_template:
					parameter_template[reaction_id] = {}

				parameter_template[reaction_id][enzyme_id] = {}
				parameter_template[reaction_id][enzyme_id]['kcat_f'] = None

				reversible = reactions[reaction_id]['is reversible']
				if reversible:
					parameter_template[reaction_id][enzyme_id]['kcat_r'] = None

				all_bound_molecules = [mol_id for set in partition for mol_id in set]

				for molecule_id in all_bound_molecules:
					parameter_template[reaction_id][enzyme_id][molecule_id] = None

		return parameter_template


if __name__ == '__main__':
	command = RateLawUtilities()
