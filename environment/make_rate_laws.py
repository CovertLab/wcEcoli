from __future__ import absolute_import, division, print_function

import scipy.constants as constants
import numpy as np

mM_to_M = 1E-3 # convert mmol/L to mol/L

class KineticFluxModel(object):
	'''
	Attributes:
		parameter_indices: a dictionary, organized by transporter -- each with 'reaction_cofactors' for
			each reaction that it is a part of, 'partition' for all of its bounded forms, and 'parameter_indices'
			for the indices at which each of its parameters can be assigned in an array.

		rate_laws: a dict, with a key for each reaction id, and then subdictionaries with each reaction's transporters
			and their rate law function. These rate laws are used directly from within this dictionary
	'''

	def __init__(self, kinetic_parameters):

		self.avogadro = constants.Avogadro
		self.kinetic_parameters = kinetic_parameters
		self.reactions = self.kinetic_parameters.keys()
		self.rate_laws = self.make_rate_laws(self.reactions, self.kinetic_parameters) #, self.transport_configuration, self.parameter_indices)


	def get_fluxes(self, concentrations_dict):
		'''
		Use rate law functions to calculate flux

		Args:
			concentrations_dict (dict) {molecule_id: concentration} - a dictionary of all relevant
		molecules and their concentrations, in mmol/L.

		Returns:
			reaction_fluxes_array - an array with fluxes for all reactions in reaction_id
			exchange_fluxes_array - an array with exchange fluxes for all molecules in molecule_id

		'''

		# Initialize reaction_fluxes and exchange_fluxes dictionaries
		reaction_fluxes = {reaction["reaction id"]: 0.0 for reaction in self.transport_reactions}
		exchange_fluxes = {mol_id: 0.0 for mol_id in self.all_molecule_ids}

		# Save reaction flux and add to exchange flux.
		for reaction in self.transport_reactions:
			reaction_id = reaction["reaction id"]
			stoichiometry = reaction["stoichiometry"]

			# Get transporter ids with compartments, and find the flux through each
			transporters = self.kinetic_parameters[reaction_id].keys()
			for transporter in transporters:

				flux = self.rate_laws[reaction_id][transporter](concentrations_dict)
				reaction_fluxes[reaction_id] += flux

				for molecule, coefficient in stoichiometry.iteritems():
					# Only use molecules in molecule_ids
					if molecule in self.all_molecule_ids:
						if coefficient < 0:
							exchange_fluxes[molecule] -= flux
						elif coefficient > 0:
							exchange_fluxes[molecule] += flux

		reaction_fluxes_array = np.asarray([reaction_fluxes[reaction_id] for reaction_id in self.reaction_ids])
		exchange_fluxes_array = np.asarray([exchange_fluxes[molecule_id] for molecule_id in self.all_molecule_ids])

		return reaction_fluxes_array, exchange_fluxes_array


	def make_rate_laws(self, reactions, kinetic_parameters):
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
						rate_law = self.construct_rate_law_no_transporter(params)

						# Save the rate law in this reaction under transporter 'None'
						rate_laws[reaction_id][transporter] = rate_law

						print('{}, has no transporter'.format(reaction_id))
					else:
						rate_law = self.construct_empty_rate_law()
						# Save the rate law in this reaction under transporter 'None'
						rate_laws[reaction_id][transporter] = rate_law

						print('{}, has no transporter, no parameters'.format(reaction_id))

				else:
					rate_law = self.construct_rate_law_piecewise_linear(
						transporter,
						params
					)

					# Save the rate law for each transporter in this reaction
					rate_laws[reaction_id][transporter] = rate_law

		return rate_laws

	def construct_empty_rate_law(self):
		def rate_law(concentrations):
			return 0.0
		return rate_law

	def construct_rate_law_piecewise_linear(self, transporter, parameters):
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

	def construct_rate_law_no_transporter(self, parameters):
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






	## Configure reactions
	def get_molecules(self, reaction_ids):
		'''
		Inputs:
			reaction_ids - a list of all reaction ids that will be used by transport
		Returns:
			self.molecule_ids - a list of all molecules used by these reactions
		'''

		molecule_ids = []
		for reaction_id, specs in self.all_transport_reactions.iteritems():
			if reaction_id in reaction_ids:
				stoichiometry = specs["stoichiometry"]
				substrates = stoichiometry.keys()
				transporters = self.kinetic_parameters[reaction_id].keys()

				# Add all relevant molecules_ids
				molecule_ids.extend(substrates)
				molecule_ids.extend(transporters)

		# Save unique molecule ids
		return list(set(molecule_ids))
