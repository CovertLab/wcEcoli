from __future__ import absolute_import, division, print_function

import scipy.constants as constants
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

	def __init__(self, make_reactions, kinetic_parameters, all_reactions):

		self.avogadro = constants.Avogadro
		self.kinetic_parameters = kinetic_parameters
		self.reaction_ids = make_reactions


		# TODO -- check that all reaction_ids are in kinetic_parameters

		# get all relevant molecule ids
		self.molecule_ids = get_molecules(all_reactions, self.reaction_ids)

		# make the rate laws
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

def get_molecules(all_reactions, reaction_ids):
	'''
	Inputs:
		reaction_ids - a list of all reaction ids that will be used by transport
	Returns:
		self.molecule_ids - a list of all molecules used by these reactions
	'''
	molecule_ids = []
	for reaction_id in reaction_ids:
		stoichiometry = all_reactions[reaction_id]['stoichiometry']
		substrates = stoichiometry.keys()
		enzymes = all_reactions[reaction_id]['catalyzed by']

		# Add all relevant molecules_ids
		molecule_ids.extend(substrates)
		molecule_ids.extend(enzymes)

	return list(set(molecule_ids))
