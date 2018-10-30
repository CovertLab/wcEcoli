from __future__ import absolute_import, division, print_function

import os
import csv
import argparse
import random # TODO -- use numpy random
import datetime
import copy

import numpy as np
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA

from wholecell.io.tablereader import TableReader

## options
ENFORCE_BOUNDS = True
POPULATION_ANALYTICS = False
USE_WCM = True

SAVE_FITNESS_THRESHOLD = 0.95
PARAM_FILE = 'best_parameters.tsv'

PARAMOUTDIR = os.path.join(
	os.path.split(__file__)[0],
	'out/saved_parameters'
	)

OUTDIR = os.path.join(
	os.path.split(__file__)[0],
	'out/plot_out'
	)

# simulation parameters
TIME_TOTAL = 10.0 # seconds
TIME_STEP = 0.1 # seconds

# genetic algorithm parameters
COHORT_SIZE = 1
POPULATION_SIZE = 200
MUTATION_VARIANCE = 1.0
MAX_GENERATIONS = 200

# set allowable parameter ranges
PARAM_RANGES = {
	'km': [1e-10, 1e2],
	'kcat': [1e-2, 1e6],
	}

## from Parest:
# Concentration are bounded from 1 fM to 1 MM.  Both are more permissive than
# expected values.  A concentration of one molecule per E. coli cell is roughly
# 1 nM, while water, the most abundant species, has a concentration of about
# 50 M.  The range is consistent with the RESOLUTION.
#
# lower_KM = 1e-10
# upper_KM = 1e2
#
# lower_kcat = 1e-2  # gives average kcat of about 100 w/ upper kcat of 1e6
# upper_kcat = 1e6  # catalase is around 1e5 /s


# estimation targets and penalties
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
		},
	'parameters' : {
		'TRANS-RXN0-270': {
			'kcat': None,
			'km_A': None,
			'km_B': None,
		},
		'TRANS-RXN0-569-LEU': {
			'kcat': None,
			'km': None,
		},
		'ABC-35-RXN': {
			'kcat_f': None,
			'kcat_r': None,
			'km_A1': None,
			'km_A2': None,
			'km_B1': None,
			'km_B2': None,
		},
		'TRANS-RXN-126B': {
			'kcat': None,
			'km_A': None,
			'km_B': None,
		},
	},
	}

# initial concentrations.
# those not set here will be taken from the WCM, or from WCM_INITIAL_CONCS
INITIAL_CONCENTRATIONS = {
	'G6458-MONOMER': 0.0,
	}

# define the reactions. Parameters set to None will be assigned
# a random initial value within a range defined by PARAM_RANGES
REACTIONS = {
			'TRANS-RXN0-270': {'type': 'symport', # export
				'substrates': {
					'A1': 'LEU[c]',
					'B1': 'PROTON[p]',
					'A2': 'LEU[p]',
					'B2': 'PROTON[c]',
				},
				'transporter': ['G6984-MONOMER', 'B4141-MONOMER'], # TODO (eran) only using the first transporter
				},
			'TRANS-RXN0-569-LEU': {'type': 'uniport', # export
				'substrates': {
				   'A1': 'LEU[c]',
				   'A2': 'LEU[p]',
				},
				'transporter': ['G6458-MONOMER'],
				},
			'ABC-35-RXN': {'type': 'symport_reversible',
				'substrates': {
				   'A1': 'LEU[p]',
				   'A2': 'LEU[c]',
				   'B1': 'ATP[c]',
				   'B2': 'ADP[c]',
				},
				'transporter': ['ABC-15-CPLX', 'ABC-304-CPLX'], # TODO (eran) only using the first transporter
				},
			'TRANS-RXN-126B': {'type': 'symport',
				'substrates': {
				   'A1': 'LEU[p]',
				   'A2': 'LEU[c]',
				   'B1': 'NA+[p]',
				   'B2': 'NA+[c]',
				},
				'transporter': ['BRNQ-MONOMER'],
				},
			}

REACTION_PARAMS = {
	'uniport' : ['kcat', 'km'],
	'uniport_reversible' : ['kcat_f', 'kcat_r', 'km_A1', 'km_A2'],
	'symport' : ['kcat', 'km_A', 'km_B'],
	'symport_reversible' : ['kcat_f', 'kcat_r', 'km_A1', 'km_A2', 'km_B1', 'km_B2'],
	}

WCM_INITIAL_CONCS = {
	'LEU[c]': 699439.9961946806,
	'NA+[p]': 60086.898174455666,
	'ADP[c]': 336487.9122664919,
	'LEU[p]': 0.0,
	'BRNQ-MONOMER': 69.78251909824665,
	'ABC-304-CPLX': 68.27370787450077,
	'G6984-MONOMER': 7.1668533127929,
	'B4141-MONOMER': 18.86014029682342,
	'ATP[c]': 5768368.6289441595,
	'NA+[c]': 0.0,
	'G6458-MONOMER': 82.60741450008658,
	'PROTON[c]': 37.72028059364684,
	'ABC-15-CPLX': 69.40531629231019,
	'PROTON[p]': 0.0
	}


class TransportEstimation(object):

	def main(self, args):

		# initialize
		if USE_WCM:
			self.initialize_from_WCM(args.simout)
		else:
			self.coefficient = 33.14430980346108 # taken from WCM
			self.cell_volume = 2.6510937465518967

		# assign indices to all parameters in problem
		self.phenotype_function, self.parameter_indices = self.get_phenotype_function()

		# set concentrations of substrates and transporters
		self.concentrations = self.initialize_concentrations(args.simout)

		## Evolve populations
		# initial population
		population = self.initialize_population()

		# run genetic algorithm
		final_population, final_fitness, saved_fitness = self.evolve_population(population)

		## Visualization and Analysis
		# make cohort_id, for naming saved output in organized fashion
		# TODO -- what if there are files in outfolder that don't match pattern?
		now = datetime.datetime.now()
		time_stamp = (str(now.month) + str(now.day))

		cohort_nums = os.listdir(OUTDIR)
		if not cohort_nums:
			cohort_num = 1
		else:
			cohort_nums = [name.replace(name[0:name.find("__")+2], '') for name in cohort_nums]
			cohort_nums = [name.replace(name[name.find("."):], '') for name in cohort_nums]
			cohort_nums = [int(name) for name in cohort_nums]
			cohort_num = max(cohort_nums) + 1

		self.cohort_id = (time_stamp + '__' + str(cohort_num))


		# Parameter analysis
		self.parameter_analysis(final_population, final_fitness)

		# Plot evolution
		self.plot_evolution(saved_fitness)

		# Run simulation of the best individual and plot output
		top_index = final_fitness.values().index(max(final_fitness.values()))
		top_genotype = final_population[top_index]

		top_phenotype = self.get_phenotype(top_genotype)

		saved_concentrations, saved_fluxes = self.run_sim(args.simout, top_phenotype)
		self.plot_out(saved_concentrations, saved_fluxes, top_phenotype)


	## Genetic algorithm

	def evolve_population(self, population):

		fitness_max = 0.9999 # TODO -- make this a global

		generations = 0
		top_fit = 0
		saved_fitness = []

		# genetic algorithm loop
		while generations < MAX_GENERATIONS and top_fit < fitness_max:

			# evaluate fitness and repopulate
			fitness = self.evaluate_fitness(population)
			population = self.repopulate(population, fitness)

			generations += 1

			# save values
			saved_fitness.append(fitness.values())
			top_fit = max(fitness.values())

			print('gen ' + str(generations) + ' fittest: ' + str(top_fit))

		if top_fit >= fitness_max:
			print('Success!')
		else:
			print('Did not find solution :-(')

		return population, fitness, saved_fitness

	def repopulate(self, population, fitness):

		new_population = {}

		# normalize fitness
		total = np.sum(fitness.values())
		normalized_fitness = {index: value/total for index, value in fitness.items()}

		index = 0
		## Elitist selection: re-seed the top individual
		top_index = fitness.values().index(max(fitness.values()))
		new_population[index] = population[top_index]
		index += 1

		while len(new_population) < POPULATION_SIZE:
			## Selection
			# TODO -- use numpy multimodal instead, for fitness-proportionate selection
			selection_index = 0
			total = normalized_fitness[selection_index]
			rand = random.uniform(0,1)
			while total < rand:
				selection_index += 1
				total += normalized_fitness[selection_index]

			## Mutation
			# TODO -- vectorize these steps, apply mutations to the entire population at once.
			genotype = population[selection_index]

			# gaussian distance
			magnitude = random.gauss(0, MUTATION_VARIANCE)

			# random unit vector
			direction = [random.gauss(0, 1) for i in range(len(genotype))]
			direction_mag = np.sum(x**2 for x in direction)**0.5
			vector = [magnitude * x / direction_mag for x in direction]

			# apply mutation
			new_genotype = np.array([x + y for x, y in zip(genotype, vector)])

			# enforce bounds on genotype
			if ENFORCE_BOUNDS:
				#X clip at bounds
				# new_genotype[new_genotype >= 1.0] = 1.0
				# new_genotype[new_genotype <= 0.0] = 0.0

				# if parameter is not in range, initialize it randomly within range
				out_of_range = np.where(np.logical_or(new_genotype <= 0.0, new_genotype >= 1.0))
				new_genotype[out_of_range] = random.uniform(0.0, 1.0)

			new_population[index] = new_genotype

			index += 1

		return new_population

	def evaluate_fitness(self, population):

		fitness = {}

		# get mean squared error
		for individual, genotype in population.iteritems():

			phenotype = self.get_phenotype(genotype)

			reaction_fluxes, molecule_fluxes = self.get_fluxes(phenotype, self.concentrations)

			error = 0.0
			for molecule, target_conc in TARGETS['concentrations'].iteritems():
				error += CONCENTRATION_PENALTY * (self.concentrations[molecule] - target_conc) ** 2
			for rxn, target_flux in TARGETS['fluxes'].iteritems():
				error += FLUX_PENALTY * (reaction_fluxes[rxn] - target_flux) ** 2
			for molecule, target_flux in TARGETS['molecule_flux'].iteritems():
				error += MOLECULE_FLUX_PENALTY * (molecule_fluxes[molecule] - target_flux) ** 2

			fitness[individual] = 1 / (1 + error)

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

	def get_phenotype_function(self):
		''' create a list that maps each parameter to an index in the parameter array'''

		parameter_indices = {rxn: {} for rxn in REACTIONS.keys()}
		phenotype_function = []

		index = 0
		# loop through all reactions
		for rxn, specs in REACTIONS.iteritems():
			parameters = REACTION_PARAMS[specs['type']]
			for param in parameters:

				if 'km' in param:
					bounds = PARAM_RANGES['km']

					g_to_p = self.make_genotype_to_phenotype(bounds[0], bounds[1])
					p_to_g = self.make_phenotype_to_genotype(bounds[0], bounds[1])

				elif 'kcat' in param:
					bounds = PARAM_RANGES['kcat']

					g_to_p = self.make_genotype_to_phenotype(bounds[0], bounds[1])
					p_to_g = self.make_phenotype_to_genotype(bounds[0], bounds[1])

				map = {
					'bounds' : bounds,
					'geno_to_pheno' : g_to_p,
					'pheno_to_geno' : p_to_g,
					}

				phenotype_function.append(map)

				parameter_indices[rxn][param] = index
				index += 1

		return phenotype_function, parameter_indices

	def make_phenotype_to_genotype(self, min_phenotype, max_phenotype):
		a = min_phenotype
		b = max_phenotype / min_phenotype

		def p_to_g(x):
			return (np.log(x / a) / np.log(b))

		return p_to_g

	def make_genotype_to_phenotype(self, min_phenotype, max_phenotype):
		a = min_phenotype
		b = max_phenotype / min_phenotype

		def g_to_p(x):
			return (a * b ** (x))

		return g_to_p

	def get_phenotype(self, genotype):
		'''convert all genes to param values using geno_to_geno mapping'''

		phenotype = np.empty(len(self.phenotype_function))
		for index, gene in enumerate(genotype):
			phenotype[index] = self.phenotype_function[index]['geno_to_pheno'](gene)

		return phenotype

	def initialize_concentrations(self, simOutDir):
		''' set all initial undefined molecular concentrations to their initial concentrations in the WCM'''

		if USE_WCM:
			concentrations = copy.deepcopy(INITIAL_CONCENTRATIONS)
			bulkMolecules = TableReader(os.path.join(simOutDir, "BulkMolecules"))
			molecule_ids = bulkMolecules.readAttribute("objectNames")

			# get all substrates in self.reactions
			for rxn, specs in REACTIONS.iteritems():
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
				# Note -- this is different from the substrate loop in that it looks for the transporter name within the strin
				for transporter in transporters:
					# if substrate is not in concentrations dict
					if transporter not in concentrations.keys() or concentrations[transporter] is None:
						transporter_id = [mol_id for mol_id in molecule_ids if transporter in mol_id]
						transporter_index = molecule_ids.index(transporter_id[0])
						transporter_count = bulkMolecules.readColumn("counts")[0, transporter_index]

						#convert to concentration
						concentrations[transporter] = transporter_count / self.cell_volume[0]

		else:
			concentrations = copy.deepcopy(INITIAL_CONCENTRATIONS)

			# get all substrates in self.reactions
			for rxn, specs in REACTIONS.iteritems():
				substrates = specs['substrates'].values()
				transporters = specs['transporter']

				# loop through substrates
				for substrate in substrates:
					# if substrate is not in concentrations dict
					if substrate not in concentrations.keys() or concentrations[substrate] is None:
						#convert to concentration
						concentrations[substrate] = WCM_INITIAL_CONCS[substrate]

				# loop through transporters
				# Note -- this is different from the substrate loop in that it looks for the transporter name within the strin
				for transporter in transporters:
					# if substrate is not in concentrations dict
					if transporter not in concentrations.keys() or concentrations[transporter] is None:

						# retrieve from saved initial concentrations
						concentrations[transporter] = WCM_INITIAL_CONCS[transporter]

		return concentrations

	def initialize_population(self):
		''' fill self.population with {index: [parameters]} for index in POPULATION_SIZE'''

		population = {}
		for ind in xrange(POPULATION_SIZE):
			population[ind] = self.initialize_genotype()

		return population

	def initialize_genotype(self):

		genotype = np.random.uniform(0, 1, len(self.phenotype_function))

		# set parameters that have a target value
		for rxn, params in TARGETS['parameters'].iteritems():
			for param, pheno_target in params.iteritems():
				if pheno_target:
					index = self.parameter_indices[rxn][param]
					gene_value = self.phenotype_function[index]['pheno_to_geno'](pheno_target)
					genotype[index] = gene_value

		return genotype

	## Transporter models
	def get_fluxes(self, parameters, concentrations):

		# TODO -- pass in a phenotype as dictionary, with parameters in there.

		reaction_fluxes = {}
		molecule_fluxes = {mol : 0.0 for mol in concentrations.keys()}

		# loop through all reactions, save reaction flux and molecule flux.
		for rxn, specs in REACTIONS.iteritems():

			params = {param : parameters[index] for param, index in self.parameter_indices[rxn].iteritems()}

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

			elif specs['type'] is 'symport':
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

			elif specs['type'] is 'symport_reversible':
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

			else:
				raise Exception('Unknown reaction type: {}'.format(specs['type']))


		return reaction_fluxes, molecule_fluxes

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
				concentration_change = flux / self.coefficient
				substrates = REACTIONS[rxn]['substrates']

				if REACTIONS[rxn]['type'] is 'symport':
					self.concentrations[substrates['A1']] -= concentration_change
					self.concentrations[substrates['B1']] -= concentration_change
					self.concentrations[substrates['A2']] += concentration_change
					self.concentrations[substrates['B2']] += concentration_change

				if REACTIONS[rxn]['type'] is 'symport_reversible':
					self.concentrations[substrates['A1']] -= concentration_change
					self.concentrations[substrates['B1']] -= concentration_change
					self.concentrations[substrates['A2']] += concentration_change
					self.concentrations[substrates['B2']] += concentration_change

				if REACTIONS[rxn]['type'] is 'uniport':
					self.concentrations[substrates['A1']] -= concentration_change
					self.concentrations[substrates['A2']] += concentration_change

			time += TIME_STEP

			# save state
			for molecule, timeseries in saved_concentrations.iteritems():
				timeseries.append(self.concentrations[molecule])

			for rxn, timeseries in saved_fluxes.iteritems():
				timeseries.append(reaction_fluxes[rxn])

		return saved_concentrations, saved_fluxes

	def plot_out(self, saved_concentrations, saved_fluxes, parameters):

			plt.figure(figsize=(16, 16))
			columns = 4
			rows = max([len(saved_concentrations), len(saved_fluxes), len(parameters)])

			index = 1
			for molecule, timeseries in saved_concentrations.iteritems():
				plt.subplot(rows, columns, columns*index-(columns-1))
				plt.plot(timeseries)
				plt.title(molecule)
				if index < len(saved_concentrations):
					plt.tick_params(
						axis='x',  # changes apply to the x-axis
						which='both',  # both major and minor ticks are affected
						bottom=False,  # ticks along the bottom edge are off
						top=False,  # ticks along the top edge are off
						labelbottom=False)  # labels along the bottom edge are off
				else:
					plt.ylabel("Concentration")
					plt.xlabel("Time (s)")
				index += 1

			index = 1
			for reaction, timeseries in saved_fluxes.iteritems():
				plt.subplot(rows, columns, columns*index-(columns-2))
				plt.plot(timeseries, 'r')
				plt.title(reaction)
				if index < len(saved_fluxes):
					plt.tick_params(
						axis='x',  # changes apply to the x-axis
						which='both',  # both major and minor ticks are affected
						bottom=False,  # ticks along the bottom edge are off
						top=False,  # ticks along the top edge are off
						labelbottom=False)  # labels along the bottom edge are off
				else:
					plt.ylabel("Flux")
					plt.xlabel("Time (s)")
				index += 1

			index = 1
			km_range = PARAM_RANGES['km']
			kcat_range = PARAM_RANGES['kcat']
			for rxn, params in self.parameter_indices.iteritems():
				for param, param_idx in params.iteritems():

					plt.subplot(rows, columns, columns*index-(columns-3))

					param_value = parameters[param_idx]
					if 'km' in param:
						plt.axvline(x=km_range[0])
						plt.axvline(x=km_range[1])
					elif 'kcat' in param:
						plt.axvline(x=kcat_range[0])
						plt.axvline(x=kcat_range[1])

					plt.axhline(y=0.5)

					plt.plot(param_value, 0.5, 'bo', markersize=10)

					info = (rxn + ' -- ' + REACTIONS[rxn]['type'] + ': ' + param)
					plt.title(info)
					plt.ylim(0, 1)
					plt.xscale("log")
					# plt.yticks([])
					plt.tick_params(top=False, bottom=False, left=False, right=False, labelleft=False, labelbottom=True)
					# plt.axis('off')

					index += 1

			plt.subplots_adjust(hspace=2.0, wspace=0.5)

			if not os.path.exists(OUTDIR):
				os.mkdir(OUTDIR)
			fig_name = ('sim_' + self.cohort_id)
			plt.savefig(os.path.join(OUTDIR,fig_name))

			print('top simulation plot saved')


	## Analysis plots of parameters and evolution
	def parameter_analysis(self, final_population, final_fitness):

		# get indices of individualus with fitness higher than 0.95
		top_fit_indices = [index for index, value in enumerate(final_fitness.values()) if value >= SAVE_FITNESS_THRESHOLD]
		top_fit_parameters = [final_population[index] for index in top_fit_indices]

		# TODO -- convert genotype to phenotype

		# save top parameters to 'best_parameters' file
		if not os.path.exists(PARAMOUTDIR):
			os.mkdir(PARAMOUTDIR)

		with open(os.path.join(PARAMOUTDIR, PARAM_FILE), 'a') as tsv_file:
			writer = csv.writer(tsv_file)
			for parameter in top_fit_parameters:
				writer.writerow(parameter)
		tsv_file.close()

		# Plot parameter space
		if POPULATION_ANALYTICS:
			# gather all saved parameter values
			with open(os.path.join(PARAMOUTDIR,PARAM_FILE), 'r') as tsv_file:
				reader = csv.reader(tsv_file)
				best_parameters = list(reader)
			tsv_file.close()
			self.plot_parameters(best_parameters)

		print('parameters analyzed')

	def plot_parameters(self, best_parameters):

		best_param_array = np.array(best_parameters, np.float64)

		pca = PCA(n_components=2)
		pca.fit(best_param_array)
		best_params_reduced = pca.transform(best_param_array)


		# todo -- cluster analysis

		plt.figure(figsize=(8.5, 8.5))
		plt.scatter(best_params_reduced[:, 0], best_params_reduced[:, 1])
		plt.xlabel('PC1')
		plt.ylabel('PC2')

		if not os.path.exists(OUTDIR):
			os.mkdir(OUTDIR)
		fig_name = ('param_space_' + self.cohort_id)
		plt.savefig(os.path.join(OUTDIR,fig_name))

		print('parameters plot saved')

	def plot_evolution(self, saved_fitness):

		n_bins = 10
		max_gens_plot = 10
		n_saved_gens = len(saved_fitness)

		if n_saved_gens >= max_gens_plot:
			nth_gen = int(n_saved_gens/max_gens_plot)
		else:
			nth_gen = 1


		plot_gen = saved_fitness[0::nth_gen]
		gen_label = ['gen ' + str(index*nth_gen) for index in xrange(len(plot_gen))]

		top_fitness = [max(fit) for fit in saved_fitness]
		avg_fitness = [sum(fit) / len(fit) for fit in saved_fitness]

		plt.figure(figsize=(8.5, 11))

		# plot fitness over gens
		plt.subplot(2, 1, 1)
		plt.plot(top_fitness)
		plt.plot(avg_fitness, 'r')
		plt.ylabel('Fitness')
		plt.xlabel('Generation')

		plt.subplot(2, 1, 2)
		# plt.hist(plot_gen, bins=np.logspace(np.log10(0.01),np.log10(1.0), n_bins), label=gen_label)
		# plt.gca().set_xscale("log")
		plt.hist(plot_gen, bins=n_bins, label=gen_label)
		plt.legend(loc='upper right')
		plt.ylabel('Count')
		plt.xlabel('Fitness')

		plt.hist(saved_fitness[0], bins=10)

		plt.subplots_adjust(hspace=0.5)



		if not os.path.exists(OUTDIR):
			os.mkdir(OUTDIR)
		fig_name = ('GA_' + self.cohort_id)
		plt.savefig(os.path.join(OUTDIR,fig_name))

		print('evolution plot saved')


if __name__ == "__main__":
	parser = argparse.ArgumentParser(description='evolve parameters for transport')
	parser.add_argument('--simout', help='directory of sim out data', default='out/manual/condition_000002/000000/generation_000000/000000/simOut')
	args = parser.parse_args()
	TransportEstimation().main(args)
