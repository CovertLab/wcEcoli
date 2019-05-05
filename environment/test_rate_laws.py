# Testing

from __future__ import absolute_import, division, print_function

import os
import csv
import json
import matplotlib.pyplot as plt

from reconstruction.spreadsheets import JsonReader
from itertools import ifilter


from environment.kinetic_rate_laws import KineticFluxModel

CSV_DIALECT = csv.excel_tab
TRANSPORT_REACTIONS_FILE = os.path.join('environment', 'condition', 'look_up_tables', 'transport_reactions.tsv')
KINETIC_PARAMETERS_FILE = os.path.join('environment', 'condition', 'parameters', 'glt_family.tsv')
EXTERNAL_MOLECULES_FILE = os.path.join('environment', 'condition', 'environment_molecules.tsv')
WCM_SIMDATA_FILE = os.path.join('environment', 'condition', 'look_up_tables', 'wcm_sim_data.json')

# TODO -- reaction_ids and parameters should ultimately come from sim_data
REACTIONS = ["RXN0-5202", "TRANS-RXN-62B"]  # competing reactions

# TODO -- how to programatically add compartment to transporter? [i] should not be hardcoded
PARAMETERS = {
  "RXN0-5202": {
     "CYCA-MONOMER[i]": {
        "kcat_f": 1e4,
        "L-ALPHA-ALANINE[p]": 1e-3,
        "GLY[p]": 1e-3,
        "PROTON[p]": 1e-3
     }
  },
  "TRANS-RXN-62B": {
     "CYCA-MONOMER[i]": {
        "kcat_f": 1e2,
        "L-ALPHA-ALANINE[p]": 1e-3,
        "GLY[p]": 1e-3,
        "PROTON[p]": 1e-3
     }
  }
}

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

	# # Make kinetic_parameters in a nested format: {reaction_id: {transporter_id : {param_id: param_value}}}
	# kinetic_parameters = {}
	# with open(KINETIC_PARAMETERS_FILE, 'rU') as csvfile:
	# 	reader = JsonReader(
	# 		ifilter(lambda x: x.lstrip()[0] != "#", csvfile),  # Strip comments
	# 		dialect=CSV_DIALECT)
	# 	for row in reader:
	# 		reaction_id = row['reaction id']
	# 		transporter = row['transporter']
	# 		k_avg = float(row['k_avg'])
	# 		json_acceptable_max_conc = row['max_conc'].replace("'", "\"")
	# 		max_conc = json.loads(json_acceptable_max_conc)
	#
	# 		# Combine kinetics into dictionary
	# 		k_param = {'k_avg': k_avg}
	# 		k_param.update(max_conc)
	# 		transporter_kinetics = {transporter: k_param}
	#
	# 		# Add to kinetic_parameters dict
	# 		if reaction_id in kinetic_parameters:
	# 			kinetic_parameters[reaction_id].update(transporter_kinetics)
	# 		else:
	# 			kinetic_parameters[reaction_id] = transporter_kinetics

	# pass
	kinetic_parameters = PARAMETERS

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
	# TODO -- get these from wcm
	concentrations = {molecule_id: 1e-2 for molecule_id in molecule_ids}
	reaction_fluxes = kinetic_rate_laws.get_fluxes(concentrations)

	# kinetic_rate_laws.analyze_rate_laws

	import ipdb; ipdb.set_trace()




def initialize_state(self, set_concentrations={}):
	''' set all initial undefined molecular concentrations to their initial concentrations in the WCM'''
	concentrations = set_concentrations.copy()  # copy.deepcopy(CONDITIONS[0]['initial_concentrations'])

	time_index = int(len(self.wcm_sim_data['time']) / 2)  # get midpoint of timeseries
	cell_volume_fL = self.wcm_sim_data['volume'][time_index]  # [fL]
	cell_volume = cell_volume_fL / 1e15  # convert to L

	initial_concentrations = {}
	for molecule, series in self.wcm_sim_data.iteritems():
		if molecule not in ['time', 'cell_mass', 'volume']:
			# convert counts to molar concentrations
			initial_concentrations[molecule] = series[time_index] / self.avogadro / cell_volume  # [fM]

	# get all substrates in REACTIONS that are not yet set
	for rxn, specs in self.reactions.iteritems():
		substrates = specs['stoichiometry'].keys()
		transporters = specs['transporters']

		# loop through substrates
		for substrate in substrates:
			# if substrate is not in concentrations dict
			if substrate not in concentrations.keys() or concentrations[substrate] is None:
				concentrations[substrate] = initial_concentrations[substrate]

		# loop through transporters
		# this is different from the substrate loop in that it looks for the transporter name within the string
		for transporter in transporters:
			# if substrate is not in concentrations dict
			if transporter not in concentrations.keys() or concentrations[transporter] is None:
				transporter_id = [mol_id for mol_id in initial_concentrations.keys() if transporter in mol_id]

				# import ipdb; ipdb.set_trace()

				concentrations[transporter] = initial_concentrations[transporter_id[0]]

	return concentrations











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


# Run test
test_rate_laws()