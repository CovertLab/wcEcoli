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


def analyze_rate_laws(kinetic_rate_laws, baseline_concentrations, output_filename):
	'''
	Args:
		kinetic_rate_laws (object): a configured kinetic_rate_law object
		baseline_concentrations (dict): concentrations for all molecules required for the rate laws

	Function:
		Runs an analysis of all rate laws in kinetic_rate_laws and saves the output in a plot

	'''

	plot_text = True
	test_transporter = True
	test_cofactor = True
	test_competitor = True

	reactions = kinetic_rate_laws.reactions
	kinetic_parameters = kinetic_rate_laws.kinetic_parameters
	rate_law_configuration = kinetic_rate_laws.rate_law_configuration

	## Plot analysis
	n_samples = 100
	n_samples_shown = 10
	n_rxns = len(reactions)

	# get concentrations array for testing and plotting
	conc_samples = np.logspace(-8, 1, num=n_samples, endpoint=True, base=10)
	conc_shown = np.flip(np.logspace(-8, 1, num=n_samples_shown, endpoint=True, base=10), 0) # use reverse ordering for legend order
	transporter_concs_shown = np.flip(np.logspace(-4, 1, num=n_samples_shown, endpoint=True, base=10), 0) # use reverse ordering for legend order

	# TODO -- should count number of transporters in each, because each gets an extra row
	rows = sum([plot_text, test_cofactor, test_transporter, test_competitor]) * n_rxns + n_rxns
	columns = 3


	cmap = plt.cm.get_cmap('Spectral')
	colors = [cmap(float(idx) / n_samples_shown) for idx in reversed(range(n_samples_shown))]  # use reverse ordering for legend order

	plt.figure(figsize=(6*columns, 3*rows))
	row_number = 0
	col_number = 0

	for reaction_id, specs in reactions.iteritems():
		transporters = specs.get('catalyzed by')
		stoich = specs.get('stoichiometry')
		parameters = kinetic_parameters.get(reaction_id)

		reactants = [mol for mol, coeff in stoich.iteritems() if coeff < 0]
		products = [mol for mol, coeff in stoich.iteritems() if coeff > 0]

		if plot_text:
			plt.subplot(rows, columns, row_number * columns + col_number + 1)
			plt.text(0.02, 0.6, 'reaction: %s' % reaction_id, weight='bold', fontsize=14)
			plt.text(0.02, 0.45, 'reactants: %s' % reactants, fontsize=14)
			plt.text(0.02, 0.3, 'products: %s' % products, fontsize=14)
			plt.text(0.02, 0.15, 'transporters: %s' % transporters, fontsize=14)
			plt.text(0.02, 0.0, 'parameters: %s' % parameters, fontsize=14, wrap=True)
			plt.axis('off')

			row_number += 1

		# test rate law by sampling substrate concentrations
		for transporter in transporters:
			if transporter not in parameters:
				continue

			# get cofactors list
			analyze_cofactors = []
			for reactant in reactants:
				if parameters[transporter][reactant] is not None:
					analyze_cofactors.append(reactant)

			# use first listed cofactor as basal TODO -- pass in cofactor of interest
			a1 = analyze_cofactors[0]

			# get competitors list
			rxns_transporter = rate_law_configuration[transporter]['reaction_cofactors'].keys()
			competing_rxns = [trpr for trpr in rxns_transporter if trpr not in reaction_id]

			analyze_competitors = []
			for rx in competing_rxns:
				competitor_candidates = rate_law_configuration[transporter]['reaction_cofactors'][rx]
				competing_mols = [mol for mol in competitor_candidates[0] if mol not in analyze_cofactors]
				analyze_competitors.extend(competing_mols)

			if test_cofactor:
				for cofactor in analyze_cofactors:

					concentrations = baseline_concentrations.copy()

					if cofactor is a1:
						# plot a1 alone

						flux_values = scan_conc(kinetic_rate_laws, concentrations, reaction_id, a1, conc_samples)

						# plot M-M curve for this reaction
						plt.subplot(rows, columns, row_number * columns + col_number + 1)
						plt.plot(conc_samples, flux_values)

						plt.xscale('log')
						plt.xlabel('%s concentration (M)' % cofactor)
						plt.ylabel('flux (M/s)')
						plt.title('cofactor: %s' % cofactor)

						col_number += 1

					else:
						# plot affect on a1 of varying cofactor
						plt.subplot(rows, columns, row_number * columns + col_number + 1)
						for index, cofactor_conc in enumerate(conc_shown):
							concentrations[cofactor] = cofactor_conc
							flux_values = scan_conc(kinetic_rate_laws, concentrations, reaction_id, a1, conc_samples)

							# plot M-M curve for this reaction
							plt.plot(conc_samples, flux_values,
								 color=colors[index],
								 label=('%.2e' % (cofactor_conc)))

						plt.legend(loc='center left', title=cofactor, bbox_to_anchor=(1.15, 0.5), prop={'size': 7})
						plt.xscale('log')
						plt.xlabel('%s concentration (M)' % a1)
						plt.ylabel('flux (M/s)')
						plt.title('cofactor: %s' % cofactor)

						col_number += 1

				col_number = 0
				row_number += 1


			if test_transporter:
				concentrations = baseline_concentrations.copy()

				plt.subplot(rows, columns, row_number * columns + col_number + 1)
				for index, transporter_conc in enumerate(transporter_concs_shown):
					concentrations[transporter] = transporter_conc
					flux_values = scan_conc(kinetic_rate_laws, concentrations, reaction_id, a1, conc_samples)

					# plot M-M curve for this reaction
					plt.plot(conc_samples, flux_values,
						color = colors[index],
						label = ('%.2e' % (transporter_conc)))

				plt.legend(loc='center left', title=transporter, bbox_to_anchor=(1.15, 0.5), prop={'size': 7})
				plt.xscale('log')
				plt.xlabel('%s concentration (M)' % a1)
				plt.ylabel('flux (M/s)')
				plt.title('transporter: %s' % transporter)

				col_number += 1


				# Test kcat
				kcat_f = parameters[transporter].get('kcat_f')

				# set concentrations
				concentrations = baseline_concentrations.copy()
				concentrations[transporter] = 1

				cofactor_rel_conc_kms = [1, 2, 100]

				plt.subplot(rows, columns, row_number * columns + col_number + 1)
				for rel_conc in cofactor_rel_conc_kms:

					for mol_index, mol_id in enumerate(analyze_cofactors):
						km = parameters[transporter].get(mol_id)
						concentrations[mol_id] = km * rel_conc  # cofactors set at km * rel_conc

					for mol_id in analyze_competitors:
						concentrations[mol_id] = 0  # remove competitors

					reaction_fluxes = kinetic_rate_laws.get_fluxes(concentrations)
					flux_value = reaction_fluxes[reaction_id]
					plt.axhline(y=flux_value, color=np.random.rand(3,), label='[cofactors] = %i*k_m' % (rel_conc))

				plt.axhline(y=kcat_f, color='r', linestyle='dashed', label='k_cat = %.6e' % kcat_f)
				plt.legend(loc='center left', title=transporter, bbox_to_anchor=(1.15, 0.5), prop={'size': 7})
				plt.ylim((0, 2*kcat_f))
				plt.title('[%s] = 1' % transporter)
				plt.tick_params(
					axis='x',  # changes apply to the x-axis
					which='both',  # both major and minor ticks are affected
					bottom=False,  # ticks along the bottom edge are off
					top=False,  # ticks along the top edge are off
					labelbottom=False)  # labels along the bottom edge are off

				col_number = 0
				row_number += 1

			if test_competitor:

				for competitor in analyze_competitors:
					concentrations = baseline_concentrations.copy()

					plt.subplot(rows, columns, row_number * columns + col_number + 1)
					for index, competitor_conc in enumerate(conc_shown):
						concentrations[competitor] = competitor_conc
						flux_values = scan_conc(kinetic_rate_laws, concentrations, reaction_id,	a1,	conc_samples)

						# plot M-M curve for this reaction
						plt.plot(conc_samples, flux_values,
							color = colors[index],
							label = ('%.2e' % (competitor_conc)))

					plt.legend(loc='center left', title=competitor, bbox_to_anchor=(1.15, 0.5), prop={'size': 7})
					plt.xscale('log')
					plt.xlabel(a1 + ' concentration (M)')
					plt.ylabel('flux (M/s)')
					plt.title('competitor: %s' % competitor)

					col_number += 1

				if analyze_competitors:
					row_number += 1

			col_number = 0

	plt.subplots_adjust(hspace=0.8, wspace=1.1)

	if not os.path.exists(OUTPUT_DIR):
		os.mkdir(OUTPUT_DIR)
	plt.savefig(os.path.join(OUTPUT_DIR, output_filename), bbox_inches='tight')

	print('rate law analysis plot saved')


def scan_conc(rate_law, concentrations, reaction_id, molecule_id, conc_samples):
	# scan the concentrations of molecule_id for a given rate law, and return flux values for the given reaction_id
	flux_values = np.empty_like(conc_samples)
	for idx, conc in enumerate(conc_samples):
		concentrations[molecule_id] = conc
		reaction_fluxes = rate_law.get_fluxes(concentrations)
		flux_values[idx] = reaction_fluxes[reaction_id]

	return flux_values

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

		# load dict of saved parameters and save parameter file name as param_id
		parameter_file = args.path
		with open(parameter_file, 'r') as fp:
			kinetic_parameters = json.load(fp)
		self.param_id = parameter_file.split(os.path.sep)[-1].split(".")[0]

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
		analyze_rate_laws(self.kinetic_rate_laws, self.concentrations, self.param_id)

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
