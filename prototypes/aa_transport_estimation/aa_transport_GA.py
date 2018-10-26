
import matplotlib.pyplot as plt
import os, cPickle
import argparse
import random
import numpy as np
from scipy import constants
import datetime
import copy

from wholecell.io.tablereader import TableReader

# simDataFile = '/Users/eranagmon/code/wcEcoli/out/manual/condition_000002/kb/simData_Modified.cPickle'
# simOutDir = '/Users/eranagmon/code/wcEcoli/out/manual/condition_000002/000000/generation_000000/000000/simOut'

OUTDIR = os.path.join(
	os.path.split(__file__)[0],
	'out'
	)

# simulation parameters
TIME_TOTAL = 10.0 # seconds
TIME_STEP = 0.1 # seconds

# genetic algorithm parameters
POPULATION_SIZE = 100
MUTATION_VARIANCE = 0.5
MAX_GENERATIONS = 100

# set allowable parameter ranges
PARAM_RANGES = {
	'km': [1.0, 10.0],
	'kcat': [0.1, 2.0],
	}

# set estimation targets and penalities
CONCENTRATION_PENALTY = 0.0
FLUX_PENALTY = 0.0
MOLECULE_FLUX_PENALTY = 1.0

TARGETS = {
	'concentrations': {
		'LEU[c]': 701000.0,
		},
	'fluxes': {
		'ABC-35-RXN': 0.0,
		'TRANS-RXN-126B': 0.05, # millimol/sec (?)
		},
	'molecule_flux': {
		'LEU[c]' : 1.427, 
		}
	}

# set initial concentrations.
# those not set here will be taken from the WCM
INITIAL_CONCENTRATIONS = {
	'G6458-MONOMER': 0.0,
	}

# define the reactions. Parameters set to None will be assigned
# a random initial value within a range defined by PARAM_RANGES
INITIAL_REACTIONS = {
			'TRANS-RXN0-270': {'type': 'symport', # export
				'substrates': {
					'A1': 'LEU[c]',
					'B1': 'PROTON[p]',
					'A2': 'LEU[p]',
					'B2': 'PROTON[c]',
				},
				'transporter': ['G6984-MONOMER', 'B4141-MONOMER'], # TODO (eran) only using the first transporter
				'params': {
				   'kcat': None,
				   'km_A': None,
				   'km_B': None,
				}},
			'TRANS-RXN0-569-LEU': {'type': 'uniport', # export
				'substrates': {
				   'A1': 'LEU[c]',
				   'A2': 'LEU[p]',
				},
				'transporter': ['G6458-MONOMER'],
				'params': {
				   'kcat': None,
				   'km': None,
				}},
			'ABC-35-RXN': {'type': 'symport_reversible',
				'substrates': {
				   'A1': 'LEU[p]',
				   'A2': 'LEU[c]',
				   'B1': 'ATP[c]',
				   'B2': 'ADP[c]',
				},
				'transporter': ['ABC-15-CPLX', 'ABC-304-CPLX'], # TODO (eran) only using the first transporter
				'params': {
					'kcat_f': None,
					'kcat_r': None,
					'km_A1': None,
					'km_A2': None,
					'km_B1': None,
					'km_B2': None,
				}},
			'TRANS-RXN-126B': {'type': 'symport',
				'substrates': {
				   'A1': 'LEU[p]',
				   'A2': 'LEU[c]',
				   'B1': 'NA+[p]',
				   'B2': 'NA+[c]',
				},
				'transporter': ['BRNQ-MONOMER'],
				'params': {
				   'kcat': None,
				   'km_A': None,
				   'km_B': None,
				}},
			}

class aa_transport_estimation:


	def main(self, args):

		# initialize
		self.initialize_from_WCM(args.simout)
		self.concentrations = self.initialize_concentrations(args.simout)
		self.parameter_indices, parameter_values = self.initialize_parameters() # TODO pass ids separate from values? we only need ids once

		population = self.initialize_population()

		# genetic algorithm loop
		generations = 0
		top_fitness = []
		avg_fitness = []

		while generations < MAX_GENERATIONS:
			fitness = self.evaluate_fitness(population)
			top_fit = max(fitness.values())
			avg_fit = sum(fitness.values())/len(fitness.values())
			top_fitness.append(top_fit)
			avg_fitness.append(avg_fit)

			population = self.repopulate(population, fitness)

			generations += 1

			print('gen ' + str(generations) + ' top_fit: ' + str(top_fit))

		self.plot_evolution(top_fitness, avg_fitness)

		# run simulation and save output for best individual
		top_idx = fitness.values().index(max(fitness.values()))
		self.run_sim(args.simout, population[top_idx])


	def repopulate(self, population, fitness):

		new_population = {}

		# normalize fitness
		total = np.sum(fitness.values())
		# TODO -- don't update this in place
		fitness.update((idx, value/total) for idx, value in fitness.items())

		idx = 0
		# re-seed top individual
		top_idx = fitness.values().index(max(fitness.values()))
		new_population[idx] = population[top_idx]
		idx += 1

		while len(new_population) < POPULATION_SIZE:
			## Selection
			selection = 0
			sum = fitness[selection]
			rand = random.uniform(0,1)
			while sum < rand:
				selection += 1
				sum += fitness[selection]

			## Mutation
			genotype = population[selection]

			# gaussian distance
			magnitude = random.gauss(0, MUTATION_VARIANCE)

			# random unit vector
			direction = [random.gauss(0, 1) for i in range(len(genotype))]
			direction_mag = np.sum(x**2 for x in direction)**0.5
			vector = [magnitude * x / direction_mag for x in direction]

			# apply mutation
			new_population[idx] = [x + y for x, y in zip(genotype, vector)]

			# TODO -- enforce bounds

			idx += 1

		return new_population

	def evaluate_fitness(self, population):

		fitness = {}

		# get mean squared error
		for ind, params in population.iteritems():

			reaction_fluxes, molecule_fluxes = self.get_fluxes(params, self.concentrations)

			error = 0.0
			for molecule, target_conc in TARGETS['concentrations'].iteritems():
				error += CONCENTRATION_PENALTY * (self.concentrations[molecule] - target_conc) ** 2
			for rxn, target_flux in TARGETS['fluxes'].iteritems():
				error += FLUX_PENALTY * (reaction_fluxes[rxn] - target_flux) ** 2
			for molecule, target_flux in TARGETS['molecule_flux'].iteritems():
				error += MOLECULE_FLUX_PENALTY * (molecule_fluxes[molecule] - target_flux) ** 2

			fitness[ind] = 1 / (1 + error)

		return fitness


	## Initialization
	def initialize_from_WCM(self, simOutDir):
		# get mass, density, and compute volume
		mass = TableReader(os.path.join(simOutDir, "Mass"))
		self.cell_mass = mass.readColumn("cellMass")
		self.dry_cell_mass = mass.readColumn("dryMass")
		mass.close()
		self.density = 1100
		self.cell_volume = self.cell_mass / self.density

		# Coefficient to convert between flux (mol/g DCW/hr) basis and concentration (M) basis
		self.coefficient = self.dry_cell_mass[0] / self.cell_mass[0] * self.density * (TIME_STEP)

	def initialize_parameters(self):

		reactions = INITIAL_REACTIONS

		parameter_indices = {rxn : {} for rxn in reactions.keys()}
		parameter_values = []

		km_range = PARAM_RANGES['km']
		kcat_range = PARAM_RANGES['kcat']

		idx = 0
		# loop through all reactions
		for rxn, specs in reactions.iteritems():
			# loop through all of the reaction's parameters
			for param, value in specs['params'].iteritems():

				# fill unassigned parameter values
				if value is None and 'km' in param:
					parameter_indices[rxn][param] = idx
					parameter_values.append(random.uniform(km_range[0], km_range[1]))

				elif value is None and 'kcat' in param:
					parameter_indices[rxn][param] = idx
					parameter_values.append(random.uniform(kcat_range[0], kcat_range[1]))

				# if value already assigned
				elif 'km' in param:
					parameter_indices[rxn][param] = idx
					parameter_values.append(value)

				elif 'kcat' in param:
					parameter_indices[rxn][param] = idx
					parameter_values.append(value)

				idx += 1

		return parameter_indices, parameter_values

	def initialize_concentrations(self, simOutDir):
		''' set all initial undefined molecular concentrations to their initial concentrations in the WCM'''

		concentrations = copy.deepcopy(INITIAL_CONCENTRATIONS)
		bulkMolecules = TableReader(os.path.join(simOutDir, "BulkMolecules"))
		molecule_ids = bulkMolecules.readAttribute("objectNames")

		# get all substrates in self.reactions
		for rxn, specs in INITIAL_REACTIONS.iteritems():
			substrates = specs['substrates'].values()
			transporters = specs['transporter']

			# loop through substrates
			for substrate in substrates:
				# if substrate is not in concentrations dict
				if substrate not in concentrations.keys() or concentrations[substrate] is None:
					substrate_index = molecule_ids.index(substrate)
					substrate_count = bulkMolecules.readColumn("counts")[0, substrate_index]
					#convert to concentration
					concentrations[substrate] = substrate_count / self.cell_volume[0]

			# loop through transporters
			for transporter in transporters:
				# if substrate is not in concentrations dict
				if transporter not in concentrations.keys() or concentrations[transporter] is None:
					transporter_id = [id for id in molecule_ids if transporter in id]
					transporter_index = molecule_ids.index(transporter_id[0])
					transporter_count = bulkMolecules.readColumn("counts")[0, transporter_index]

					#convert to concentration
					concentrations[transporter] = transporter_count / self.cell_volume[0]

		bulkMolecules.close()

		return concentrations

	def initialize_population(self):
		''' fill self.population with {index: [parameters]} for index in POPULATION_SIZE'''

		population = {}

		for ind in xrange(POPULATION_SIZE):
			parameter_indices, parameter_values = self.initialize_parameters()
			population[ind] = parameter_values

		return population


	## Get flux for all reactions and molecules
	def get_fluxes(self, parameters, concentrations):

		reaction_fluxes = {}
		molecule_fluxes = {mol : 0.0 for mol in concentrations.keys()}

		# loop through all reactions, save reaction flux and molecule flux.
		for rxn, specs in INITIAL_REACTIONS.iteritems():

			params = {param : parameters[idx] for param, idx in self.parameter_indices[rxn].iteritems()}

			if specs['type'] is 'uniport':
				reaction_fluxes[rxn] = self.get_uniporter_flux(
					specs['substrates'],
					specs['transporter'][0],
					params,
					concentrations,
					)

				# molecule flux
				A1 = specs['substrates']['A1']
				A2 = specs['substrates']['A2']
				molecule_fluxes[A1] -= reaction_fluxes[rxn]
				molecule_fluxes[A2] += reaction_fluxes[rxn]

			# TODO -- add if specs['type'] is 'uniport_reversible'

			if specs['type'] is 'symport':
				reaction_fluxes[rxn] = self.get_symporter_flux(
					specs['substrates'],
					specs['transporter'][0],
					params,
					concentrations,
					)

				# molecule flux
				A1 = specs['substrates']['A1']
				A2 = specs['substrates']['A2']
				B1 = specs['substrates']['B1']
				B2 = specs['substrates']['B2']
				molecule_fluxes[A1] -= reaction_fluxes[rxn]
				molecule_fluxes[A2] += reaction_fluxes[rxn]
				molecule_fluxes[B1] -= reaction_fluxes[rxn]
				molecule_fluxes[B2] += reaction_fluxes[rxn]

			if specs['type'] is 'symport_reversible':
				reaction_fluxes[rxn] = self.get_symporter_reversible_flux(
					specs['substrates'],
					specs['transporter'][0],
					params,
					concentrations,
					)

				# molecule flux
				A1 = specs['substrates']['A1']
				A2 = specs['substrates']['A2']
				B1 = specs['substrates']['B1']
				B2 = specs['substrates']['B2']
				molecule_fluxes[A1] -= reaction_fluxes[rxn]
				molecule_fluxes[A2] += reaction_fluxes[rxn]
				molecule_fluxes[B1] -= reaction_fluxes[rxn]
				molecule_fluxes[B2] += reaction_fluxes[rxn]

		return reaction_fluxes, molecule_fluxes


	## Transporter models
	def get_symporter_flux(self, substrates_dict, transporter, parameters_dict, concentrations):

		# concentrations
		A1 = concentrations[substrates_dict['A1']]
		B1 = concentrations[substrates_dict['B1']]
		T = concentrations[transporter]

		# parameters
		km_A = parameters_dict['km_A']
		km_B = parameters_dict['km_B']
		kcat = parameters_dict['kcat']

		flux = T * (kcat * A1/km_A * B1/km_B) / ((1 + A1/km_A) * (1 + B1/km_B))

		return flux

	def get_symporter_reversible_flux(self, substrates_dict, transporter, parameters_dict, concentrations):

		# concentrations
		A1 = concentrations[substrates_dict['A1']]
		A2 = concentrations[substrates_dict['A2']]
		B1 = concentrations[substrates_dict['B1']]
		B2 = concentrations[substrates_dict['B2']]
		T = concentrations[transporter]

		# parameters
		km_A1 = parameters_dict['km_A1']
		km_A2 = parameters_dict['km_A2']
		km_B1 = parameters_dict['km_B1']
		km_B2 = parameters_dict['km_B2']
		kcat_f = parameters_dict['kcat_f']
		kcat_r = parameters_dict['kcat_r']

		flux =  T * (kcat_f * A1/km_A1 * B1/km_B1 - kcat_r * A2/km_A2 * B2/km_B2) / (
				(1 + A1/km_A1) * (1 + B1/km_B1) + (1 + A2/km_A2) * (1 + B2/km_B2) - 1
				)

		return flux

	def get_uniporter_flux(self, substrates_dict, transporter, parameters_dict, concentrations):

		# concentrations
		A1 = concentrations[substrates_dict['A1']]
		A2 = concentrations[substrates_dict['A2']]
		T = concentrations[transporter]

		# parameters
		km = parameters_dict['km']
		kcat = parameters_dict['kcat']

		flux = kcat * T * A1 / (A1 + km)

		return flux

	def get_uniporter_reversible_flux(self, substrates_dict, transporter, parameters_dict, concentrations):

		# concentrations
		A1 = concentrations[substrates_dict['A1']]
		A2 = concentrations[substrates_dict['A2']]
		T = concentrations[transporter]

		# parameters
		km_A1 = parameters_dict['km_A1']
		km_A2 = parameters_dict['km_A2']
		kcat_f = parameters_dict['kcat_f']
		kcat_r = parameters_dict['kcat_r']

		flux = T * (kcat_f * A1/km_A1 - kcat_r * A2/km_A2) / (1 + A1/km_A1 + A2/km_A2)

		return flux


	## Run simulation. plot time series of concentrations, flux
	def run_sim(self, simOutDir, parameters):
		time = 0

		# initialize concentrations
		self.concentrations = self.initialize_concentrations(simOutDir)

		# create dicts for saved timeseries, starting with initial concentrations and fluxes = 0
		saved_concentrations = {mol : [conc] for mol, conc in self.concentrations.iteritems()}

		reaction_fluxes, molecule_fluxes = self.get_fluxes(parameters, self.concentrations)
		saved_fluxes = {rxn : [0.0] for rxn, flux in reaction_fluxes.iteritems()}


		while time < TIME_TOTAL:
			reaction_fluxes, molecule_fluxes = self.get_fluxes(parameters, self.concentrations)

			for rxn, flux in reaction_fluxes.iteritems():
				exchange_flux = flux / self.coefficient
				substrates = INITIAL_REACTIONS[rxn]['substrates']

				if INITIAL_REACTIONS[rxn]['type'] is 'symport':
					self.concentrations[substrates['A1']] -= exchange_flux
					self.concentrations[substrates['B1']] -= exchange_flux
					self.concentrations[substrates['A2']] += exchange_flux
					self.concentrations[substrates['B2']] += exchange_flux

				if INITIAL_REACTIONS[rxn]['type'] is 'symport_reversible':
					self.concentrations[substrates['A1']] -= exchange_flux
					self.concentrations[substrates['B1']] -= exchange_flux
					self.concentrations[substrates['A2']] += exchange_flux
					self.concentrations[substrates['B2']] += exchange_flux

				if INITIAL_REACTIONS[rxn]['type'] is 'uniport':
					self.concentrations[substrates['A1']] -= exchange_flux
					self.concentrations[substrates['A2']] += exchange_flux

			time += TIME_STEP

			# save state
			for molecule, timeseries in saved_concentrations.iteritems():
				timeseries.append(self.concentrations[molecule])

			for rxn, timeseries in saved_fluxes.iteritems():
				timeseries.append(reaction_fluxes[rxn])

		self.plot_out(saved_concentrations, saved_fluxes, parameters)

	def plot_out(self, saved_concentrations, saved_fluxes, parameters):

		plt.figure(figsize=(8.5, 11))
		columns = 2
		rows = len(saved_concentrations.keys())

		idx = 1
		for molecule, timeseries in saved_concentrations.iteritems():
			plt.subplot(rows, columns, 2*idx-1)
			plt.plot(timeseries)
			plt.title(molecule)
			if idx < len(saved_concentrations.keys()):
				plt.tick_params(
					axis='x',  # changes apply to the x-axis
					which='both',  # both major and minor ticks are affected
					bottom=False,  # ticks along the bottom edge are off
					top=False,  # ticks along the top edge are off
					labelbottom=False)  # labels along the bottom edge are off
			else:
				plt.ylabel("Concentration")
				plt.xlabel("Time (s)")
			idx += 1

		idx = 1
		for reaction, timeseries in saved_fluxes.iteritems():
			plt.subplot(rows, columns, 2*idx)
			plt.plot(timeseries, 'r')
			plt.title(reaction)
			if idx < len(saved_fluxes.keys()):
				plt.tick_params(
					axis='x',  # changes apply to the x-axis
					which='both',  # both major and minor ticks are affected
					bottom=False,  # ticks along the bottom edge are off
					top=False,  # ticks along the top edge are off
					labelbottom=False)  # labels along the bottom edge are off
			else:
				plt.ylabel("Flux")
				plt.xlabel("Time (s)")
			idx += 1

		plt.subplot(int(rows/2), columns, int(rows/2) * columns)

		# save params in dictionary and add to figure as text
		reaction_params = {} #rxn : {} for rxn, params in self.parameter_indices.iteritems()}
		for rxn, params in self.parameter_indices.iteritems():
			reaction_params[rxn] = {param : parameters[idx] for param, idx in self.parameter_indices[rxn].iteritems()}

		plt.text(0.0, 1.0, reaction_params, wrap=True)
		plt.axis('off')

		plt.subplots_adjust(hspace=0.9, wspace=0.5)
		now = datetime.datetime.now()

		if not os.path.exists(OUTDIR):
			os.mkdir(OUTDIR)
		fig_name = ('sim_' + str(now.month) + '-' + str(now.day) + '_' + str(now.hour) + str(now.minute) + str(now.second))
		plt.savefig(os.path.join(OUTDIR,fig_name))


	def plot_evolution(self, top_fitness, avg_fitness):

		plt.plot(top_fitness)
		plt.plot(avg_fitness, 'r')
		plt.ylabel('Fitness')
		plt.xlabel('Generation')

		now = datetime.datetime.now()

		if not os.path.exists(OUTDIR):
			os.mkdir(OUTDIR)
		fig_name = ('GA_' + str(now.month) + '-' + str(now.day) + '_' + str(now.hour) + str(now.minute) + str(now.second))
		plt.savefig(os.path.join(OUTDIR,fig_name))


if __name__ == "__main__":
	parser = argparse.ArgumentParser(description='evolve parameters for transport')
	parser.add_argument('--simout', help='directory of sim out data', default='out/manual/condition_000002/000000/generation_000000/000000/simOut')
	args = parser.parse_args()
	aa_transport_estimation().main(args)
