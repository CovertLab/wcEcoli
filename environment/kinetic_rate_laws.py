from __future__ import absolute_import, division, print_function

import scipy.constants as constants
import numpy as np
import matplotlib.pyplot as plt

mM_to_M = 1E-3 # convert mmol/L to mol/L

class KineticFluxModel(object):
	'''
	Args:
		make_reactions (list): these are the reactions that will be parameterized

		kinetic_parameters (dict): a dictionary of parameters a nested format:
			{reaction_id: {
				transporter_id : {
					param_id: param_value}}}

		all_reactions (dict): all metabolic reactions, with
			{reaction_id: {
				'catalyzed by': list,
				'is reversible': bool,
				'stoichiometry': dict,
				}}

	Attributes:
		rate_laws: a dict, with a key for each reaction id, and then subdictionaries with each reaction's transporters
			and their rate law function. These rate laws are used directly from within this dictionary
	'''

	def __init__(self, reactions, kinetic_parameters):

		self.avogadro = constants.Avogadro
		self.kinetic_parameters = kinetic_parameters

		self.reactions = reactions
		self.reaction_ids = self.reactions.keys()

		# TODO -- check that all reaction_ids are in kinetic_parameters

		# get all relevant molecule ids
		self.molecule_ids = get_molecules(self.reactions)

		convenience_kinetics_rate_laws = True
		minimal_rate_laws = False

		# make the rate laws
		if convenience_kinetics_rate_laws:
			self.rate_law_configuration = make_configuration(self.reactions)

			self.rate_laws = make_rate_laws(
				self.reactions,
				self.rate_law_configuration,
				self.kinetic_parameters)

		elif minimal_rate_laws:
			self.rate_laws = make_rate_laws_minimal(self.reaction_ids, self.kinetic_parameters)


	def get_fluxes(self, concentrations_dict):
		'''
		Use rate law functions to calculate flux

		Args:
			concentrations_dict (dict) {molecule_id: concentration} - a dictionary of all relevant
		molecules and their concentrations, in mmol/L.

		Returns:
			reaction_fluxes (dict) - with fluxes for all reactionss
		'''

		# Initialize reaction_fluxes and exchange_fluxes dictionaries
		reaction_fluxes = {reaction_id: 0.0 for reaction_id in self.reaction_ids}

		for reaction_id, transporters in self.rate_laws.iteritems():
			for transporter, rate_law in transporters.iteritems():
				flux = self.rate_laws[reaction_id][transporter](concentrations_dict)
				reaction_fluxes[reaction_id] += flux

		return reaction_fluxes


	def analyze_rate_laws(self, all_parameters):

		test_transporter = True
		test_cofactor = True
		test_competitor = True

		columns = 1 + sum([test_transporter, test_cofactor, test_competitor])

		n_vary = 10
		n_samples = 100
		n_rxns = len(self.parameter_indices)
		rows = 2*n_rxns + 2  # extra row for each reaction header

		cmap = plt.cm.get_cmap('Spectral')
		colors = [cmap(float(idx) / n_vary) for idx in range(n_vary)]

		baseline_concentrations = self.kinetic_model.baseline_concentrations


		plt.figure(figsize=(6*columns, 3*rows))
		plot_number = 1
		row_number = 0
		for reaction_id, specs in self.reactions.iteritems():
			transporters = specs['transporters']
			stoich = specs['stoichiometry']
			parameters = self.parameter_indices[reaction_id]

			# TODO -- set a1 to amino acid... or show all?
			reactants = [mol for mol, coeff in stoich.iteritems() if coeff < 0]
			products = [mol for mol, coeff in stoich.iteritems() if coeff > 0]

			a1_set = False
			for mol in self.exchange_molecules:
				if mol in reactants:
					a1 = mol
					a1_set = True

			if not a1_set:
				a1 = reactants[0]

			# get cofactor
			b1 = None
			if len(reactants) > 1:
				cofactors = [x for x in reactants if x != a1]
				b1 = cofactors[0]

			# plot info in whole row
			param_values = {}
			for trans, params in self.parameter_indices[reaction_id].iteritems():
				param_values[trans] = {}
				param_values[trans]['kms'] = {}
				for param_type, params in params.iteritems():
					if 'km' in param_type:
						for param, idx in params.iteritems():
							param_values[trans]['kms'][param] = all_parameters[idx]
					else:
						param_values[trans][param_type] = all_parameters[params]

			plt.subplot(rows, columns, plot_number)
			plt.text(0.02, 0.6, 'reaction: ' + reaction_id, weight='bold')
			plt.text(0.02, 0.45, 'reactants: %s' % reactants)
			plt.text(0.02, 0.3, 'products: %s' % products)
			plt.text(0.02, 0.15, 'transporters: %s' % transporters)
			plt.text(0.02, 0.0, 'parameters: %s' % param_values[transporters[0]])
			plt.axis('off')
			plot_number += columns
			row_number += 1

			# michaelis menten by sampling substrate concentrations
			for transporter in transporters:

				concentrations = baseline_concentrations.copy()
				conc_values = np.logspace(-9, 0, num=n_samples, endpoint=True, base=10)

				flux_values = np.empty_like(conc_values)
				for idx, conc in enumerate(conc_values):
					concentrations[a1] = conc
					reaction_fluxes, exchange_fluxes = self.kinetic_model.get_fluxes(all_parameters, concentrations)
					flux_values[idx] = reaction_fluxes[reaction_id]

				# plot M-M curve for this reaction
				plt.subplot(rows, columns, plot_number)
				plt.plot(conc_values, flux_values)

				# plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
				plt.xscale('log')
				plt.xlabel(a1 + ' concentration (M)')
				plt.ylabel('flux (M/s)')
				plt.title('transporter: %s' % transporter)

				plot_number += 1

				if test_transporter:
					concentrations = baseline_concentrations.copy()
					conc_values = np.logspace(-8, 1, num=n_samples, endpoint=True, base=10)
					transporter_concs = np.logspace(-4, 1, num=n_vary, endpoint=True, base=10)

					plt.subplot(rows, columns, plot_number)
					for index, transporter_conc in enumerate(transporter_concs):
						concentrations[transporter] = transporter_conc

						flux_values = np.empty_like(conc_values)
						for idx, conc in enumerate(conc_values):

							concentrations[a1] = conc
							reaction_fluxes, exchange_fluxes = self.kinetic_model.get_fluxes(all_parameters, concentrations)
							flux_values[idx] = reaction_fluxes[reaction_id]

						# plot M-M curve for this reaction
						plt.plot(conc_values, flux_values,
												color = colors[index],
												label = ('conc = %.2e' % (transporter_conc)),
												)

					plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
					plt.xscale('log')
					plt.xlabel(a1 + ' concentration (M)')
					plt.ylabel('flux (M/s)')
					plt.title('transporter: %s' % transporter)

					plot_number += 1

				if test_cofactor:

					concentrations = baseline_concentrations.copy()
					conc_values = np.logspace(-8, 1, num=n_samples, endpoint=True, base=10)
					cofactor_concs = np.logspace(-8, 1, num=n_vary, endpoint=True, base=10)

					if b1 is not None:
						plt.subplot(rows, columns, plot_number)
						for index, cofactor_conc in enumerate(cofactor_concs):
							concentrations[b1] = cofactor_conc

							flux_values = np.empty_like(conc_values)
							for idx, conc in enumerate(conc_values):

								concentrations[a1] = conc
								reaction_fluxes, exchange_fluxes = self.kinetic_model.get_fluxes(all_parameters, concentrations)
								flux_values[idx] = reaction_fluxes[reaction_id]

							# plot M-M curve for this reaction
							plt.plot(conc_values, flux_values,
											color = colors[index],
											label = ('conc = %.2e' % (cofactor_conc)),
											)

						plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
						plt.xscale('log')
						plt.xlabel(a1 + ' concentration (M)')
						plt.ylabel('flux (M/s)')
						plt.title('cofactor: %s' % b1)

					plot_number += 1

				if test_competitor:
					# get competitor
					rxns_transporter = self.rate_law_configuration[transporter]['reaction_cofactors'].keys()
					competing_rxns = [trpr for trpr in rxns_transporter if trpr not in reaction_id]
					competitor = None
					for rx in competing_rxns:
						competitors = self.rate_law_configuration[transporter]['reaction_cofactors'][rx]
						competitor = competitors[0][0]

					if competitor is not None:
						concentrations = baseline_concentrations.copy()
						conc_values = np.logspace(-8, 1, num=n_samples, endpoint=True, base=10)
						competitor_concs = np.logspace(-8, 1, num=n_vary, endpoint=True, base=10)

						plt.subplot(rows, columns, plot_number)
						for index, competitor_conc in enumerate(competitor_concs):
							concentrations[competitor] = competitor_conc

							flux_values = np.empty_like(conc_values)
							for idx, conc in enumerate(conc_values):

								concentrations[a1] = conc
								reaction_fluxes, exchange_fluxes = self.kinetic_model.get_fluxes(all_parameters, concentrations)
								flux_values[idx] = reaction_fluxes[reaction_id]

							# plot M-M curve for this reaction
							plt.plot(conc_values, flux_values,
													color = colors[index],
													label = ('conc = %.2e' % (competitor_conc)),
													)

						plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
						plt.xscale('log')
						plt.xlabel(a1 + ' concentration (M)')
						plt.ylabel('flux (M/s)')
						plt.title('competitor: %s' % (competitor))

					plot_number += 1

				row_number += 1

			plot_number = row_number * columns + 1

		plt.subplots_adjust(hspace=0.5, wspace=1.5)

		if not os.path.exists(self.out_dir):
			os.mkdir(self.out_dir)
		fig_name = ('MM_' + self.replicate_id)
		plt.savefig(os.path.join(self.out_dir, fig_name), bbox_inches='tight')

		print('rate law analysis plot saved')







def make_configuration(reactions):
	'''
	Make the rate law configuration, which tells the parameters where to be placed.

	Args:
		reactions: A dictionary with all reactions that will be made into rate laws.

	Returns:
		rate_law_configuration: a dictionary with partition and reaction_cofactor entries for each reaction
	'''

	rate_law_configuration = {}
	# gets all potential interactions between the reactions
	for reaction_id, specs in reactions.iteritems():
		transporters = specs['catalyzed by']
		# initialize all transporters' entries
		for transporter in transporters:
			if transporter not in rate_law_configuration:
				rate_law_configuration[transporter] = {
					'partition': [],
					'reaction_cofactors': {},
				}

	# identify parameters for reactions
	for reaction_id, specs in reactions.iteritems():
		stoich = specs.get('stoichiometry')
		transporters = specs.get('catalyzed by', None)
		reversibility = specs.get('is reversible', False)

		# get sets of cofactors driving this reaction
		forward_cofactors = [mol for mol, coeff in stoich.iteritems() if coeff < 0]
		cofactors = [forward_cofactors]

		if reversibility:
			reverse_cofactors = [mol for mol, coeff in stoich.iteritems() if coeff > 0]
			cofactors.append(reverse_cofactors)

		# get partition, reactions, and parameter indices for each transporter, and save to rate_law_configuration dictionary
		for transporter in transporters:

			# get competition for this transporter from all other reactions
			competing_reactions = [rxn for rxn, specs2 in reactions.iteritems() if
					(rxn is not reaction_id) and (transporter in specs2['catalyzed by'])]

			competitors = []
			for reaction2 in competing_reactions:
				stoich2 = reactions[reaction2]['stoichiometry']
				reactants2 = [mol for mol, coeff in stoich2.iteritems() if coeff < 0]
				competitors.append(reactants2)

			# partition includes both competitors and cofactors.
			partition = competitors + cofactors
			rate_law_configuration[transporter]['partition'] = partition
			rate_law_configuration[transporter]['reaction_cofactors'][reaction_id] = cofactors

	return rate_law_configuration



def get_molecules(reactions):
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




# Make rate laws

# def make_rate_laws_minimal(reactions, kinetic_parameters):
def make_rate_laws(reactions, rate_law_configuration, kinetic_parameters):

	# make rate law for each reaction
	rate_laws = {reaction_id: {} for reaction_id in reactions.iterkeys()}
	for reaction_id, specs in reactions.iteritems():
		stoichiometry = specs.get('stoichiometry')
		# reversible = specs.get('is reversible')
		transporters = specs.get('catalyzed by')

		# rate law for each transporter
		# TODO -- make sure that transporter is in the rate law configuration
		for transporter in transporters:
			cofactors_sets = rate_law_configuration[transporter]["reaction_cofactors"][reaction_id]
			partition = rate_law_configuration[transporter]["partition"]

		rate_law = construct_rate_law(
			stoichiometry,
			transporter,
			cofactors_sets,
			partition,
			kinetic_parameters[reaction_id][transporter]
		)

		# save the rate law for each transporter in this reaction
		rate_laws[reaction_id][transporter] = rate_law

	return rate_laws

def cofactor_numerator(concentration, km):
	def term():
		return concentration / km if km else 0

	return term

def cofactor_denominator(concentration, km):
	def term():
		return 1 + concentration / km if km else 1

	return term

def construct_rate_law(stoichiometry, transporter, cofactors_sets, partition, parameters):
	'''
	Args:
		cofactors: a list with the required cofactors , each pair needs a kcat.
		partition: a list of lists. each sublist is the set of cofactors for a given partition.
			[[C1, C2],[C3, C4], [C5]

	Returns:
		a kinetic rate law for the reaction, with arguments for concentrations and parameters,
		and returns flux.
	'''

	kcat_f = parameters.get('kcat_f')
	kcat_r = parameters.get('kcat_r')




	import ipdb;
	ipdb.set_trace()




	def rate_law(concentrations):

		# construct numerator
		transporter_concentration = concentrations[transporter]

		numerator = 0
		for cofactors in cofactors_sets:
			# if reversible, determine direction by looking at stoichiometry
			if kcat_r:
				coeffs = [stoichiometry[mol] for mol in cofactors]  # coeffs should be of the same sign
				if coeffs[0] > 0:
					kcat = -kcat_r  # use reverse rate
				else:
					kcat = kcat_f
			else:
				kcat = kcat_f

			# multiply the affinities of all cofactors
			term = np.prod([
				cofactor_numerator(
					concentrations[molecule],
					parameters[molecule]  # km of molecule
				)() for molecule in cofactors
			])
			numerator += kcat * term
		numerator *= transporter_concentration

		# construct denominator, with all competing terms in the partition
		# denominator starts at +1 for the unbound state
		denominator = 1
		for cofactors_set in partition:
			# multiply the affinities of all cofactors in this partition
			term = np.prod([
				cofactor_denominator(
					concentrations[molecule],
					parameters[molecule]
				)() for molecule in cofactors_set
			])
			denominator += term - 1
		flux = numerator / denominator

		return flux

	return rate_law






## Piecewise linear rate laws
def make_rate_laws_minimal(reactions, kinetic_parameters):
	'''
	Make rate laws for all reactions

	Returns:
		rate_laws (dict) {'reaction id': {transporter: rate_law}} - A dictionary with
			rate laws saved as functions, in nested.
	'''

	# Make rate laws for each reaction
	rate_laws = {reaction_id: {} for reaction_id in reactions}
	for reaction_id in reactions:
		transporters = kinetic_parameters[reaction_id]

		# Rate law for each transporter
		for transporter, params in transporters.iteritems():

			if transporter == 'None':
				del params['k_avg']
				# If there are remaining parameters (max_conc), build rate law
				if params:
					rate_law = construct_rate_law_no_transporter(params)

					# Save the rate law in this reaction under transporter 'None'
					rate_laws[reaction_id][transporter] = rate_law

					print('{}, has no transporter'.format(reaction_id))
				else:
					rate_law = construct_empty_rate_law()
					# Save the rate law in this reaction under transporter 'None'
					rate_laws[reaction_id][transporter] = rate_law

					print('{}, has no transporter, no parameters'.format(reaction_id))

			else:
				rate_law = construct_rate_law_piecewise_linear(
					transporter,
					params
				)

				# Save the rate law for each transporter in this reaction
				rate_laws[reaction_id][transporter] = rate_law

	return rate_laws

def construct_empty_rate_law():
	def rate_law(concentrations):
		return 0.0
	return rate_law

def construct_rate_law_piecewise_linear(transporter, parameters):
	'''
	Makes a single piecewise linear rate law.

	Inputs:
		transporter (str)
		parameters (dict) {parameter: value} - requires a 'k_avg' entry, and max_conc entries for all
			relevant substrates.

	Returns:
		rate_law - a function for the kinetic rate law of the given transporter, which takes a
			dictionary of concentrations as input, and which returns flux as output
	'''

	k_avg = parameters.get('k_avg')
	max_conc = dict(parameters)
	del max_conc['k_avg']

	# add units. max_conc is in mmol/L --> convert to mol/L
	max_conc = {mol_id : value * mM_to_M for mol_id, value in max_conc.iteritems()}

	# Construct rate law for this transporter
	def rate_law(concentrations):
		transporter_concentration = concentrations[transporter]

		substrates_normalized = 1.0
		for substrate, substrate_max_conc in max_conc.iteritems():
			# If concentration is above the max, assume saturated
			if concentrations[substrate] > substrate_max_conc:
				continue
			else:
				substrates_normalized *= concentrations[substrate] / substrate_max_conc

		flux = k_avg * transporter_concentration * substrates_normalized
		return flux

	return rate_law

def construct_rate_law_no_transporter(parameters):
	'''
	Makes a single piecewise linear rate law, based entirely on substrate concentrations.

	Inputs:
		parameters (dict) {parameter: value} - requires a 'k_avg' entry, and max_conc
			entries for all relevant substrates.

	Returns:
		rate_law - a function for the kinetic rate law, which takes a
			dictionary of concentrations as input, and which returns flux as output
	'''

	max_conc = dict(parameters)

	# add units. max_conc is in mmol/L --> convert to mol/L
	max_conc = {mol_id : value * mM_to_M for mol_id, value in max_conc.iteritems()}

	# Construct rate law for this transporter
	def rate_law(concentrations):
		substrates_normalized = 1.0
		for substrate, substrate_max_conc in max_conc.iteritems():
			# If concentration is above the max, assume saturated
			if concentrations[substrate] > substrate_max_conc:
				continue
			else:
				substrates_normalized *= concentrations[substrate] / substrate_max_conc

		flux = substrates_normalized

		# TODO -- we don't want a flux of 1... need a different rate law
		# return flux
		return 0.0

	return rate_law



















import os
import csv
from reconstruction.spreadsheets import JsonReader
from itertools import ifilter
import json

CSV_DIALECT = csv.excel_tab
TRANSPORT_REACTIONS_FILE = os.path.join('environment', 'condition', 'look_up_tables', 'transport_reactions.tsv')
KINETIC_PARAMETERS_FILE = os.path.join('environment', 'condition', 'parameters', 'glt_family.tsv')
EXTERNAL_MOLECULES_FILE = os.path.join('environment', 'condition', 'environment_molecules.tsv')
WCM_SIMDATA_FILE = os.path.join('environment', 'condition', 'look_up_tables', 'wcm_sim_data.json')


# self.baseline_concentrations = initialize_state(self.set_baseline)

def test_rate_laws():
	# Make dict of transport reactions
	all_transport_reactions = {}
	with open(TRANSPORT_REACTIONS_FILE, 'rU') as csvfile:
		reader = JsonReader(
			ifilter(lambda x: x.lstrip()[0] != "#", csvfile),  # Strip comments
			dialect=CSV_DIALECT)
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

	# Make kinetic_parameters in a nested format: {reaction_id: {transporter_id : {param_id: param_value}}}
	kinetic_parameters = {}
	with open(KINETIC_PARAMETERS_FILE, 'rU') as csvfile:
		reader = JsonReader(
			ifilter(lambda x: x.lstrip()[0] != "#", csvfile),  # Strip comments
			dialect=CSV_DIALECT)
		for row in reader:
			reaction_id = row['reaction id']
			transporter = row['transporter']
			k_avg = float(row['k_avg'])
			json_acceptable_max_conc = row['max_conc'].replace("'", "\"")
			max_conc = json.loads(json_acceptable_max_conc)

			# Combine kinetics into dictionary
			k_param = {'k_avg': k_avg}
			k_param.update(max_conc)
			transporter_kinetics = {transporter: k_param}

			# Add to kinetic_parameters dict
			if reaction_id in kinetic_parameters:
				kinetic_parameters[reaction_id].update(transporter_kinetics)
			else:
				kinetic_parameters[reaction_id] = transporter_kinetics

	# make a dict of reactions that will be configured with the parameters
	make_reaction_ids = kinetic_parameters.keys()
	make_reactions = {
		reaction_id: specs
		for reaction_id, specs in all_transport_reactions.iteritems()
		if reaction_id in make_reaction_ids}


	# Make the kinetic model
	kinetic_rate_laws = KineticFluxModel(make_reactions, kinetic_parameters)

	# Get list of molecule_ids used by kinetic rate laws
	molecule_ids = kinetic_rate_laws.molecule_ids

	# initialize concentrations and get fluxes
	concentrations = {molecule_id: 1e-2 for molecule_id in molecule_ids}
	reaction_fluxes = kinetic_rate_laws.get_fluxes(concentrations)

	import ipdb; ipdb.set_trace()


test_rate_laws()
