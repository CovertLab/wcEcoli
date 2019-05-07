from __future__ import absolute_import, division, print_function

import numpy as np

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

		self.kinetic_parameters = kinetic_parameters

		self.reactions = reactions
		self.reaction_ids = self.reactions.keys()

		# TODO -- check that all reaction_ids are in kinetic_parameters

		# get all relevant molecule ids
		self.molecule_ids = get_molecules(self.reactions)

		# make the rate laws
		self.rate_law_configuration = {}
		self.rate_law_configuration = make_configuration(self.reactions)

		self.rate_laws = make_rate_laws(
			self.reactions,
			self.rate_law_configuration,
			self.kinetic_parameters)


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

def get_parameter_template(reactions, rate_law_configuration):
	'''
	Given a rate law configuration, return a template for required parameters

	Args:
		reactions:
		rate_law_configuration:

	Returns:
		parameter_template (dict): a template for all parameters required by this rate_law_configuration,
			filled with values of 0.0.

	'''
	parameter_template = {}
	for enzyme_id, configuration in rate_law_configuration.iteritems():
		reaction_cofactors = configuration['reaction_cofactors']
		partition = configuration['partition']

		for reaction_id, cofactors in reaction_cofactors.iteritems():

			# check if reaction is already in the template
			if reaction_id not in parameter_template:
				parameter_template[reaction_id] = {}

			parameter_template[reaction_id][enzyme_id] = {}
			parameter_template[reaction_id][enzyme_id]['kcat_f'] = 0.0

			reversible = reactions[reaction_id]['is reversible']
			if reversible:
				parameter_template[reaction_id][enzyme_id]['kcat_r'] = 0.0

			all_bound_molecules = [mol_id for set in partition for mol_id in set]

			for molecule_id in all_bound_molecules:
				parameter_template[reaction_id][enzyme_id][molecule_id] = 0.0

	return parameter_template


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



## Make rate laws

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

		rate_law = construct_convenience_rate_law(
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

def construct_convenience_rate_law(stoichiometry, transporter, cofactors_sets, partition, parameters):
	'''
	Make a convenience kinetics rate law for one transporter

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

	# remove km parameters with None as their value
	for parameter, value in parameters.iteritems():
		if 'kcat' not in parameter:
			if value is None:
				for part in partition:
					if parameter in part:
						part.remove(parameter)
				for set in cofactors_sets:
					if parameter in set:
						set.remove(parameter)

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
