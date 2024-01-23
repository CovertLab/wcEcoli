"""
The parca, aka parameter calculator.

TODO: establish a controlled language for function behaviors (i.e. create* set* fit*)
TODO: functionalize so that values are not both set and returned from some methods
"""

import binascii
import functools
import itertools
import os
import pickle
import time
import traceback
from typing import Callable, List

from stochastic_arrow import StochasticSystem
from cvxpy import Variable, Problem, Minimize, norm
import numpy as np
import scipy.optimize
import scipy.sparse

from reconstruction.ecoli.initialization import create_bulk_container
from reconstruction.ecoli.simulation_data import SimulationDataEcoli
from wholecell.containers.bulk_objects_container import BulkObjectsContainer
from wholecell.utils import filepath, parallelization, units
from wholecell.utils.fitting import normalize, masses_and_counts_for_homeostatic_target


# Fitting parameters
FITNESS_THRESHOLD = 1e-9
MAX_FITTING_ITERATIONS = 100
N_SEEDS = 10

# Parameters used in fitPromoterBoundProbability()
PROMOTER_PDIFF_THRESHOLD = 0.06  # Minimum difference between binding probabilities of a TF in conditions where TF is active and inactive
PROMOTER_REG_COEFF = 1e-3  # Optimization weight on how much probability should stay close to original values
PROMOTER_SCALING = 10  # Multiplied to all matrices for numerical stability
PROMOTER_NORM_TYPE = 1  # Matrix 1-norm
PROMOTER_MAX_ITERATIONS = 100
PROMOTER_CONVERGENCE_THRESHOLD = 1e-9
ECOS_0_TOLERANCE = 1e-10  # Tolerance to adjust solver output to 0

BASAL_EXPRESSION_CONDITION = "M9 Glucose minus AAs"

VERBOSE = 1

COUNTS_UNITS = units.dmol
VOLUME_UNITS = units.L
MASS_UNITS = units.g
TIME_UNITS = units.s

functions_run = []


def fitSimData_1(raw_data, **kwargs):
	"""
	Fits parameters necessary for the simulation based on the knowledge base

	Inputs:
		raw_data (KnowledgeBaseEcoli) - knowledge base consisting of the
			necessary raw data
		cpus (int) - number of processes to use (if > 1, use multiprocessing)
		debug (bool) - if True, fit only one arbitrarily-chosen transcription
			factor in order to speed up a debug cycle (should not be used for
			an actual simulation)
		save_intermediates (bool) - if True, save the state (sim_data and cell_specs)
			to disk in intermediates_directory after each Parca step
		intermediates_directory (str) - path to the directory to save intermediate
			sim_data and cell_specs files to
		load_intermediate (str) - the function name of the Parca step to load
			sim_data and cell_specs from; functions prior to and including this
			will be skipped but all following functions will run
		variable_elongation_transcription (bool) - enable variable elongation
			for transcription
		variable_elongation_translation (bool) - enable variable elongation for
			translation
		disable_ribosome_capacity_fitting (bool) - if True, ribosome expression
			is not fit to protein synthesis demands
		disable_rnapoly_capacity_fitting (bool) - if True, RNA polymerase
			expression is not fit to protein synthesis demands
	"""

	sim_data = SimulationDataEcoli()
	cell_specs = {}

	# Functions to modify sim_data and/or cell_specs
	# Functions defined below should be wrapped by @save_state to allow saving
	# and loading sim_data and cell_specs to skip certain functions while doing
	# development for faster testing and iteration of later functions that
	# might not need earlier functions to be rerun each time.
	sim_data, cell_specs = initialize(sim_data, cell_specs, raw_data=raw_data, **kwargs)
	sim_data, cell_specs = input_adjustments(sim_data, cell_specs, **kwargs)
	sim_data, cell_specs = basal_specs(sim_data, cell_specs, **kwargs)
	sim_data, cell_specs = tf_condition_specs(sim_data, cell_specs, **kwargs)
	sim_data, cell_specs = fit_condition(sim_data, cell_specs, **kwargs)
	sim_data, cell_specs = promoter_binding(sim_data, cell_specs, **kwargs)
	sim_data, cell_specs = adjust_promoters(sim_data, cell_specs, **kwargs)
	sim_data, cell_specs = set_conditions(sim_data, cell_specs, **kwargs)
	sim_data, cell_specs = final_adjustments(sim_data, cell_specs, **kwargs)

	if sim_data is None:
		raise ValueError('sim_data is not specified.  Check that the'
			f' load_intermediate function ({kwargs.get("load_intermediate")})'
			' is correct and matches a function to be run.')

	return sim_data

def save_state(func):
	"""
	Wrapper for functions called in fitSimData_1() to allow saving and loading
	of sim_data and cell_specs at different points in the parameter calculation
	pipeline.  This is useful for development in order to skip time intensive
	steps that are not required to recalculate in order to work with the desired
	stage of parameter calculation.

	This wrapper expects arguments in the kwargs passed into a wrapped function:
		save_intermediates (bool): if True, the state (sim_data and cell_specs)
			will be saved to disk in intermediates_directory
		intermediates_directory (str): path to the directory to save intermediate
			sim_data and cell_specs files to
		load_intermediate (str): the name of the function to load sim_data and
			cell_specs from, functions prior to and including this will be
			skipped but all following functions will run
	"""

	@functools.wraps(func)
	def wrapper(*args, **kwargs):
		func_name = func.__name__
		load_intermediate = kwargs.get('load_intermediate')
		intermediates_dir = kwargs.get('intermediates_directory', '')

		# Files to save to or load from
		sim_data_file = os.path.join(intermediates_dir, f'sim_data_{func_name}.cPickle')
		cell_specs_file = os.path.join(intermediates_dir, f'cell_specs_{func_name}.cPickle')

		# Run the wrapped function if the function to load is not specified or was already loaded
		if load_intermediate is None or load_intermediate in functions_run:
			start = time.time()
			sim_data, cell_specs = func(*args, **kwargs)
			end = time.time()
			print(f'Ran {func_name} in {end - start:.0f} s')
		# Load the saved results from the wrapped function if it is set to be loaded
		elif load_intermediate == func_name:
			if not os.path.exists(sim_data_file) or not os.path.exists(cell_specs_file):
				raise IOError(f'Could not find intermediate files ({sim_data_file}'
					f' or {cell_specs_file}) to load. Make sure to save intermediates'
					' before trying to load them.')
			with open(sim_data_file, 'rb') as f:
				sim_data = pickle.load(f)
			with open(cell_specs_file, 'rb') as f:
				cell_specs = pickle.load(f)
			print(f'Loaded sim_data and cell_specs for {func_name}')
		# Skip running or loading if a later function will be loaded
		else:
			print(f'Skipped {func_name}')
			sim_data = None
			cell_specs = {}

		# Save the current state of the parameter calculator after the function to disk
		if kwargs.get('save_intermediates', False) and intermediates_dir != '' and sim_data is not None:
			os.makedirs(intermediates_dir, exist_ok=True)
			with open(sim_data_file, 'wb') as f:
				pickle.dump(sim_data, f, protocol=pickle.HIGHEST_PROTOCOL)
			with open(cell_specs_file, 'wb') as f:
				pickle.dump(cell_specs, f, protocol=pickle.HIGHEST_PROTOCOL)
			print(f'Saved data for {func_name}')

		# Record which functions have been run to know if the loaded function has run
		functions_run.append(func_name)

		return sim_data, cell_specs
	return wrapper

@save_state
def initialize(sim_data, cell_specs, raw_data=None, **kwargs):
	sim_data.initialize(
		raw_data = raw_data,
		basal_expression_condition = BASAL_EXPRESSION_CONDITION,
		)

	return sim_data, cell_specs

@save_state
def input_adjustments(sim_data, cell_specs, debug=False, **kwargs):
	# Limit the number of conditions that are being fit so that execution time decreases
	if debug:
		print("Warning: Running the Parca in debug mode - not all conditions will be fit")
		key = list(sim_data.tf_to_active_inactive_conditions.keys())[0]
		sim_data.tf_to_active_inactive_conditions = {key: sim_data.tf_to_active_inactive_conditions[key]}

	# Make adjustments for metabolic enzymes
	setTranslationEfficiencies(sim_data)
	set_balanced_translation_efficiencies(sim_data)
	setRNAExpression(sim_data)
	setRNADegRates(sim_data)
	setProteinDegRates(sim_data)

	return sim_data, cell_specs

@save_state
def basal_specs(sim_data, cell_specs,
		disable_ribosome_capacity_fitting=False, disable_rnapoly_capacity_fitting=False,
		variable_elongation_transcription=True, variable_elongation_translation=False,
		**kwargs):

	cell_specs = buildBasalCellSpecifications(
		sim_data,
		variable_elongation_transcription,
		variable_elongation_translation,
		disable_ribosome_capacity_fitting,
		disable_rnapoly_capacity_fitting,
		)

	# Set expression based on ppGpp regulation from basal expression
	sim_data.process.transcription.set_ppgpp_expression(sim_data)
	# TODO (Travis): use ppGpp expression in condition fitting

	# Modify other properties
	# Compute Km's
	Km = setKmCooperativeEndoRNonLinearRNAdecay(sim_data, cell_specs["basal"]["bulkContainer"])
	n_transcribed_rnas = len(sim_data.process.transcription.rna_data)
	sim_data.process.transcription.rna_data['Km_endoRNase'] = Km[:n_transcribed_rnas]
	sim_data.process.transcription.mature_rna_data['Km_endoRNase'] = Km[n_transcribed_rnas:]

	## Calculate and set maintenance values
	# ----- Growth associated maintenance -----
	fitMaintenanceCosts(sim_data, cell_specs["basal"]["bulkContainer"])

	return sim_data, cell_specs

@save_state
def tf_condition_specs(sim_data, cell_specs, cpus=1,
		disable_ribosome_capacity_fitting=False, disable_rnapoly_capacity_fitting=False,
		variable_elongation_transcription=True, variable_elongation_translation=False,
		**kwargs):
	# Limit the number of CPUs before printing it to stdout.
	cpus = parallelization.cpus(cpus)

	# Apply updates to cell_specs from buildTfConditionCellSpecifications for each TF condition
	conditions = list(sorted(sim_data.tf_to_active_inactive_conditions))
	args = [
		(sim_data, tf, variable_elongation_transcription, variable_elongation_translation,
				disable_ribosome_capacity_fitting, disable_rnapoly_capacity_fitting)
		for tf in conditions]
	apply_updates(buildTfConditionCellSpecifications, args, conditions, cell_specs, cpus)

	for conditionKey in cell_specs:
		if conditionKey == "basal":
			continue

		sim_data.process.transcription.rna_expression[conditionKey] = cell_specs[conditionKey]["expression"]
		sim_data.process.transcription.rna_synth_prob[conditionKey] = cell_specs[conditionKey]["synthProb"]
		sim_data.process.transcription.cistron_expression[conditionKey] = cell_specs[conditionKey]["cistron_expression"]
		sim_data.process.transcription.fit_cistron_expression[conditionKey] = cell_specs[conditionKey]["fit_cistron_expression"]

	buildCombinedConditionCellSpecifications(
		sim_data,
		cell_specs,
		variable_elongation_transcription,
		variable_elongation_translation,
		disable_ribosome_capacity_fitting,
		disable_rnapoly_capacity_fitting)

	return sim_data, cell_specs

@save_state
def fit_condition(sim_data, cell_specs, cpus=1, **kwargs):
	# Apply updates from fitCondition to cell_specs for each fit condition
	conditions = list(sorted(cell_specs))
	args = [(sim_data, cell_specs[condition], condition)
		for condition in conditions]
	apply_updates(fitCondition, args, conditions, cell_specs, cpus)

	for condition_label in sorted(cell_specs):
		nutrients = sim_data.conditions[condition_label]["nutrients"]
		if nutrients not in sim_data.translation_supply_rate:
			sim_data.translation_supply_rate[nutrients] = cell_specs[condition_label]["translation_aa_supply"]

	return sim_data, cell_specs

@save_state
def promoter_binding(sim_data, cell_specs, **kwargs):
	if VERBOSE > 0:
		print('Fitting promoter binding')
	# noinspection PyTypeChecker
	fitPromoterBoundProbability(sim_data, cell_specs)

	return sim_data, cell_specs

@save_state
def adjust_promoters(sim_data, cell_specs, **kwargs):
	# noinspection PyTypeChecker
	fitLigandConcentrations(sim_data, cell_specs)
	calculateRnapRecruitment(sim_data, cell_specs)

	return sim_data, cell_specs

@save_state
def set_conditions(sim_data, cell_specs, **kwargs):
	sim_data.process.transcription.rnaSynthProbFraction = {}
	sim_data.process.transcription.rnapFractionActiveDict = {}
	sim_data.process.transcription.rnaSynthProbRProtein = {}
	sim_data.process.transcription.rnaSynthProbRnaPolymerase = {}
	sim_data.process.transcription.rnaPolymeraseElongationRateDict = {}
	sim_data.expectedDryMassIncreaseDict = {}
	sim_data.process.translation.ribosomeElongationRateDict = {}
	sim_data.process.translation.ribosomeFractionActiveDict = {}

	for condition_label in sorted(cell_specs):
		condition = sim_data.conditions[condition_label]
		nutrients = condition["nutrients"]

		if VERBOSE > 0:
			print("Updating mass in condition {}".format(condition_label))
		spec = cell_specs[condition_label]

		concDict = sim_data.process.metabolism.concentration_updates.concentrations_based_on_nutrients(media_id=nutrients)
		concDict.update(sim_data.mass.getBiomassAsConcentrations(sim_data.condition_to_doubling_time[condition_label]))

		avgCellDryMassInit, fitAvgSolublePoolMass = rescaleMassForSolubleMetabolites(
			sim_data, spec["bulkContainer"], concDict, sim_data.condition_to_doubling_time[condition_label]
			)

		if VERBOSE > 0:
			print('{} to {}'.format(spec["avgCellDryMassInit"], avgCellDryMassInit))

		spec["avgCellDryMassInit"] = avgCellDryMassInit
		spec["fitAvgSolublePoolMass"] = fitAvgSolublePoolMass

		mRnaSynthProb = sim_data.process.transcription.rna_synth_prob[condition_label][sim_data.process.transcription.rna_data['is_mRNA']].sum()
		tRnaSynthProb = sim_data.process.transcription.rna_synth_prob[condition_label][sim_data.process.transcription.rna_data['is_tRNA']].sum()
		rRnaSynthProb = sim_data.process.transcription.rna_synth_prob[condition_label][sim_data.process.transcription.rna_data['is_rRNA']].sum()

		if len(condition["perturbations"]) == 0:
			if nutrients not in sim_data.process.transcription.rnaSynthProbFraction:
				sim_data.process.transcription.rnaSynthProbFraction[nutrients] = {
					"mRna": mRnaSynthProb,
					"tRna": tRnaSynthProb,
					"rRna": rRnaSynthProb,
					}

			if nutrients not in sim_data.process.transcription.rnaSynthProbRProtein:
				prob = sim_data.process.transcription.rna_synth_prob[condition_label][sim_data.process.transcription.rna_data['includes_ribosomal_protein']]
				sim_data.process.transcription.rnaSynthProbRProtein[nutrients] = prob

			if nutrients not in sim_data.process.transcription.rnaSynthProbRnaPolymerase:
				prob = sim_data.process.transcription.rna_synth_prob[condition_label][sim_data.process.transcription.rna_data['includes_RNAP']]
				sim_data.process.transcription.rnaSynthProbRnaPolymerase[nutrients] = prob

			if nutrients not in sim_data.process.transcription.rnapFractionActiveDict:
				frac = sim_data.growth_rate_parameters.get_fraction_active_rnap(spec["doubling_time"])
				sim_data.process.transcription.rnapFractionActiveDict[nutrients] = frac

			if nutrients not in sim_data.process.transcription.rnaPolymeraseElongationRateDict:
				rate = sim_data.growth_rate_parameters.get_rnap_elongation_rate(spec["doubling_time"])
				sim_data.process.transcription.rnaPolymeraseElongationRateDict[nutrients] = rate

			if nutrients not in sim_data.expectedDryMassIncreaseDict:
				sim_data.expectedDryMassIncreaseDict[nutrients] = spec["avgCellDryMassInit"]

			if nutrients not in sim_data.process.translation.ribosomeElongationRateDict:
				rate = sim_data.growth_rate_parameters.get_ribosome_elongation_rate(spec["doubling_time"])
				sim_data.process.translation.ribosomeElongationRateDict[nutrients] = rate

			if nutrients not in sim_data.process.translation.ribosomeFractionActiveDict:
				frac = sim_data.growth_rate_parameters.get_fraction_active_ribosome(spec["doubling_time"])
				sim_data.process.translation.ribosomeFractionActiveDict[nutrients] = frac

	return sim_data, cell_specs

@save_state
def final_adjustments(sim_data, cell_specs, **kwargs):
	# Adjust expression for RNA attenuation
	sim_data.process.transcription.calculate_attenuation(sim_data, cell_specs)

	# Adjust ppGpp regulated expression after conditions have been fit for physiological constraints
	sim_data.process.transcription.adjust_polymerizing_ppgpp_expression(sim_data)
	sim_data.process.transcription.adjust_ppgpp_expression_for_tfs(sim_data)

	# Set supply constants for amino acids based on condition supply requirements
	average_basal_container = create_bulk_container(sim_data, n_seeds=5)
	average_with_aa_container = create_bulk_container(sim_data, condition='with_aa', n_seeds=5)
	sim_data.process.metabolism.set_phenomological_supply_constants(sim_data)
	sim_data.process.metabolism.set_mechanistic_supply_constants(sim_data, cell_specs,
		average_basal_container, average_with_aa_container)
	sim_data.process.metabolism.set_mechanistic_export_constants(sim_data, cell_specs,
		average_basal_container)
	sim_data.process.metabolism.set_mechanistic_uptake_constants(sim_data, cell_specs,
		average_with_aa_container)

	# Set ppGpp reaction parameters
	sim_data.process.transcription.set_ppgpp_kinetics_parameters(average_basal_container, sim_data.constants)

	return sim_data, cell_specs

def apply_updates(func, args, labels, dest, cpus):
	# type: (Callable[..., dict], List[tuple], List[str], dict, int) -> None
	"""
	Use multiprocessing (if cpus > 1) to apply args to a function to get
	dictionary updates for a destination dictionary.

	Args:
		func: function to call with args
		args: list of args to apply to func
		labels: label for each set of args for exception information
		dest: destination dictionary that will be updated with results
			from each function call
		cpus: number of cpus to use
	"""

	if cpus > 1:
		print("Starting {} Parca processes".format(cpus))

		# Apply args to func
		pool = parallelization.pool(cpus)
		results = {
			label: pool.apply_async(func, a)
			for label, a in zip(labels, args)
			}
		pool.close()
		pool.join()

		# Check results from function calls and update dest
		failed = []
		for label, result in results.items():
			if result.successful():
				dest.update(result.get())
			else:
				# noinspection PyBroadException
				try:
					result.get()
				except Exception as e:
					traceback.print_exc()
					failed.append(label)

		# Cleanup
		if failed:
			raise RuntimeError('Error(s) raised for {} while using multiple processes'
				.format(', '.join(failed)))
		pool = None
		print("End parallel processing")
	else:
		for a in args:
			dest.update(func(*a))

def buildBasalCellSpecifications(
		sim_data,
		variable_elongation_transcription=True,
		variable_elongation_translation=False,
		disable_ribosome_capacity_fitting=False,
		disable_rnapoly_capacity_fitting=False,
		):
	"""
	Creates cell specifications for the basal condition by fitting expression.
	Relies on expressionConverge() to set the expression and update masses.

	Inputs
	------
	- disable_ribosome_capacity_fitting (bool) - if True, ribosome expression
	is not fit
	- disable_rnapoly_capacity_fitting (bool) - if True, RNA polymerase
	expression is not fit

	Requires
	--------
	- Metabolite concentrations based on 'minimal' nutrients
	- 'basal' RNA expression
	- 'basal' doubling time

	Modifies
	--------
	- Average mass values of the cell
	- cistron expression
	- RNA expression and synthesis probabilities

	Returns
	--------
	- dict {'basal': dict} with the following keys in the dict from key 'basal':
		'concDict' {metabolite_name (str): concentration (float with units)} -
			dictionary of concentrations for each metabolite with a concentration
		'fit_cistron_expression' (array of floats) - hypothetical expression for
			each RNA cistron post-fit, total normalized to 1, if all
			transcription units were monocistronic
		'expression' (array of floats) - expression for each RNA, total normalized to 1
		'doubling_time' (float with units) - cell doubling time
		'synthProb' (array of floats) - synthesis probability for each RNA,
			total normalized to 1
		'avgCellDryMassInit' (float with units) - average initial cell dry mass
		'fitAvgSolubleTargetMolMass' (float with units) - the adjusted dry mass
			of the soluble fraction of a cell
		'bulkContainer' (BulkObjectsContainer object) - expected counts for
			bulk molecules based on expression

	Notes
	-----
	- TODO - sets sim_data attributes and returns values - change to only return values
	"""

	# Create dictionary for basal condition
	cell_specs = {}
	cell_specs["basal"] = {
		"concDict": sim_data.process.metabolism.concentration_updates.concentrations_based_on_nutrients(media_id="minimal"),
		"expression": sim_data.process.transcription.rna_expression["basal"].copy(),
		"doubling_time": sim_data.condition_to_doubling_time["basal"],
		}

	# Determine expression and synthesis probabilities
	expression, synthProb, fit_cistron_expression, avgCellDryMassInit, fitAvgSolubleTargetMolMass, bulkContainer, _ = expressionConverge(
		sim_data,
		cell_specs["basal"]["expression"],
		cell_specs["basal"]["concDict"],
		cell_specs["basal"]["doubling_time"],
		variable_elongation_transcription = variable_elongation_transcription,
		variable_elongation_translation = variable_elongation_translation,
		disable_ribosome_capacity_fitting = disable_ribosome_capacity_fitting,
		disable_rnapoly_capacity_fitting = disable_rnapoly_capacity_fitting
		)

	# Store calculated values
	cell_specs["basal"]["expression"] = expression
	cell_specs["basal"]["synthProb"] = synthProb
	cell_specs["basal"]["fit_cistron_expression"] = fit_cistron_expression
	cell_specs["basal"]["avgCellDryMassInit"] = avgCellDryMassInit
	cell_specs["basal"]["fitAvgSolubleTargetMolMass"] = fitAvgSolubleTargetMolMass
	cell_specs["basal"]["bulkContainer"] = bulkContainer

	# Modify sim_data mass
	sim_data.mass.avg_cell_dry_mass_init = avgCellDryMassInit
	sim_data.mass.avg_cell_dry_mass = sim_data.mass.avg_cell_dry_mass_init * sim_data.mass.avg_cell_to_initial_cell_conversion_factor
	sim_data.mass.avg_cell_water_mass_init = sim_data.mass.avg_cell_dry_mass_init / sim_data.mass.cell_dry_mass_fraction * sim_data.mass.cell_water_mass_fraction
	sim_data.mass.fitAvgSolubleTargetMolMass = fitAvgSolubleTargetMolMass

	# Modify sim_data expression
	sim_data.process.transcription.rna_expression["basal"][:] = cell_specs["basal"]["expression"]
	sim_data.process.transcription.rna_synth_prob["basal"][:] = cell_specs["basal"]["synthProb"]
	sim_data.process.transcription.fit_cistron_expression["basal"] = cell_specs["basal"]["fit_cistron_expression"]

	return cell_specs

def buildTfConditionCellSpecifications(
		sim_data,
		tf,
		variable_elongation_transcription=True,
		variable_elongation_translation=False,
		disable_ribosome_capacity_fitting=False,
		disable_rnapoly_capacity_fitting=False,
		):
	"""
	Creates cell specifications for a given transcription factor by
	fitting expression. Will set for the active and inactive TF condition.
	Relies on expressionConverge() to set the expression and masses.
	Uses fold change data relative to the 'basal' condition to determine
	expression for a given TF.

	Inputs
	------
	- tf (str) - label for the transcription factor to fit (eg. 'CPLX-125')
	- disable_ribosome_capacity_fitting (bool) - if True, ribosome expression
	is not fit
	- disable_rnapoly_capacity_fitting (bool) - if True, RNA polymerase
	expression is not fit

	Requires
	--------
	- Metabolite concentrations based on nutrients for the TF
	- Adjusted 'basal' cistron expression
	- Doubling time for the TF
	- Fold changes in expression for each gene given the TF

	Returns
	--------
	- dict {tf + '__active'/'__inactive': dict} with the following keys in each dict:
		'concDict' {metabolite_name (str): concentration (float with units)} -
			dictionary of concentrations for each metabolite with a concentration
		'expression' (array of floats) - expression for each RNA, total normalized to 1
		'doubling_time' (float with units) - cell doubling time
		'synthProb' (array of floats) - synthesis probability for each RNA,
			total normalized to 1
		'cistron_expression' (array of floats) - hypothetical expression for
			each RNA cistron, calculated from basal cistron expression levels
			and fold change data
		'fit_cistron_expression' (array of floats) - hypothetical expression for
			each RNA cistron post-fit, total normalized to 1, if all
			transcription units were monocistronic
		'avgCellDryMassInit' (float with units) - average initial cell dry mass
		'fitAvgSolubleTargetMolMass' (float with units) - the adjusted dry mass
			of the soluble fraction of a cell
		'bulkContainer' (BulkObjectsContainer object) - expected counts for
			bulk molecules based on expression
	"""

	cell_specs = {}
	for choice in ["__active", "__inactive"]:
		conditionKey = tf + choice
		conditionValue = sim_data.conditions[conditionKey]

		# Get expression for the condition based on fold changes over 'basal'
		# condition if the condition is not the same as 'basal'
		fcData = {}
		if choice == "__active" and conditionValue != sim_data.conditions["basal"]:
			fcData = sim_data.tf_to_fold_change[tf]
		if choice == "__inactive" and conditionValue != sim_data.conditions["basal"]:
			fcDataTmp = sim_data.tf_to_fold_change[tf].copy()
			for key, value in fcDataTmp.items():
				fcData[key] = 1. / value
		expression, cistron_expression = expressionFromConditionAndFoldChange(
			sim_data.process.transcription,
			conditionValue["perturbations"],
			fcData,
			)

		# Get metabolite concentrations for the condition
		concDict = sim_data.process.metabolism.concentration_updates.concentrations_based_on_nutrients(
			media_id=conditionValue["nutrients"]
			)
		concDict.update(sim_data.mass.getBiomassAsConcentrations(sim_data.condition_to_doubling_time[conditionKey]))

		# Create dictionary for the condition
		cell_specs[conditionKey] = {
			"concDict": concDict,
			"expression": expression,
			"doubling_time": sim_data.condition_to_doubling_time.get(
				conditionKey,
				sim_data.condition_to_doubling_time["basal"]
				)
			}

		# Determine expression and synthesis probabilities
		expression, synthProb, fit_cistron_expression, avgCellDryMassInit, fitAvgSolubleTargetMolMass, bulkContainer, concDict = expressionConverge(
			sim_data,
			cell_specs[conditionKey]["expression"],
			cell_specs[conditionKey]["concDict"],
			cell_specs[conditionKey]["doubling_time"],
			sim_data.process.transcription.rna_data['Km_endoRNase'],
			variable_elongation_transcription = variable_elongation_transcription,
			variable_elongation_translation = variable_elongation_translation,
			disable_ribosome_capacity_fitting = disable_ribosome_capacity_fitting,
			disable_rnapoly_capacity_fitting = disable_rnapoly_capacity_fitting
			)

		# Store calculated values
		cell_specs[conditionKey]["expression"] = expression
		cell_specs[conditionKey]["synthProb"] = synthProb
		cell_specs[conditionKey]["cistron_expression"] = cistron_expression
		cell_specs[conditionKey]["fit_cistron_expression"] = fit_cistron_expression
		cell_specs[conditionKey]["avgCellDryMassInit"] = avgCellDryMassInit
		cell_specs[conditionKey]["fitAvgSolubleTargetMolMass"] = fitAvgSolubleTargetMolMass
		cell_specs[conditionKey]["bulkContainer"] = bulkContainer

	return cell_specs

def buildCombinedConditionCellSpecifications(
		sim_data,
		cell_specs,
		variable_elongation_transcription=True,
		variable_elongation_translation=False,
		disable_ribosome_capacity_fitting=False,
		disable_rnapoly_capacity_fitting=False,
		):
	"""
	Creates cell specifications for sets of transcription factors being active.
	These sets include conditions like 'with_aa' or 'no_oxygen' where multiple
	transcription factors will be active at the same time.

	Inputs
	------
	- cell_specs {condition (str): dict} - information about each individual
	transcription factor condition
	- disable_ribosome_capacity_fitting (bool) - if True, ribosome expression
	is not fit
	- disable_rnapoly_capacity_fitting (bool) - if True, RNA polymerase
	expression is not fit

	Requires
	--------
	- Metabolite concentrations based on nutrients for the condition
	- Adjusted 'basal' RNA expression
	- Doubling time for the combined condition
	- Fold changes in expression for each gene given the TF

	Modifies
	--------
	- cell_specs dictionary for each combined condition
	- RNA expression and synthesis probabilities for each combined condition

	Notes
	-----
	- TODO - determine how to handle fold changes when multiple TFs change the
	same gene because multiplying both fold changes together might not be
	appropriate
	"""

	for conditionKey in sim_data.condition_active_tfs:
		# Skip adjustments if 'basal' condition
		if conditionKey == "basal":
			continue

		# Get expression from fold changes for each TF in the given condition
		fcData = {}
		conditionValue = sim_data.conditions[conditionKey]
		for tf in sim_data.condition_active_tfs[conditionKey]:
			for gene, fc in sim_data.tf_to_fold_change[tf].items():
				fcData[gene] = fcData.get(gene, 1) * fc
		for tf in sim_data.condition_inactive_tfs[conditionKey]:
			for gene, fc in sim_data.tf_to_fold_change[tf].items():
				fcData[gene] = fcData.get(gene, 1) / fc

		expression, cistron_expression = expressionFromConditionAndFoldChange(
			sim_data.process.transcription,
			conditionValue["perturbations"],
			fcData,
			)

		# Get metabolite concentrations for the condition
		concDict = sim_data.process.metabolism.concentration_updates.concentrations_based_on_nutrients(
			media_id=conditionValue["nutrients"]
			)
		concDict.update(sim_data.mass.getBiomassAsConcentrations(sim_data.condition_to_doubling_time[conditionKey]))

		# Create dictionary for the condition
		cell_specs[conditionKey] = {
			"concDict": concDict,
			"expression": expression,
			"doubling_time": sim_data.condition_to_doubling_time.get(
				conditionKey,
				sim_data.condition_to_doubling_time["basal"]
				)
			}

		# Determine expression and synthesis probabilities
		expression, synthProb, fit_cistron_expression, avgCellDryMassInit, fitAvgSolubleTargetMolMass, bulkContainer, concDict = expressionConverge(
			sim_data,
			cell_specs[conditionKey]["expression"],
			cell_specs[conditionKey]["concDict"],
			cell_specs[conditionKey]["doubling_time"],
			sim_data.process.transcription.rna_data['Km_endoRNase'],
			variable_elongation_transcription = variable_elongation_transcription,
			variable_elongation_translation = variable_elongation_translation,
			disable_ribosome_capacity_fitting = disable_ribosome_capacity_fitting,
			disable_rnapoly_capacity_fitting = disable_rnapoly_capacity_fitting)

		# Modify cell_specs for calculated values
		cell_specs[conditionKey]["expression"] = expression
		cell_specs[conditionKey]["synthProb"] = synthProb
		cell_specs[conditionKey]["cistron_expression"] = cistron_expression
		cell_specs[conditionKey]["fit_cistron_expression"] = fit_cistron_expression
		cell_specs[conditionKey]["avgCellDryMassInit"] = avgCellDryMassInit
		cell_specs[conditionKey]["fitAvgSolubleTargetMolMass"] = fitAvgSolubleTargetMolMass
		cell_specs[conditionKey]["bulkContainer"] = bulkContainer

		# Modify sim_data expression
		sim_data.process.transcription.rna_expression[conditionKey] = cell_specs[conditionKey]["expression"]
		sim_data.process.transcription.rna_synth_prob[conditionKey] = cell_specs[conditionKey]["synthProb"]
		sim_data.process.transcription.cistron_expression[conditionKey] = cell_specs[conditionKey]["cistron_expression"]
		sim_data.process.transcription.fit_cistron_expression[conditionKey] = cell_specs[conditionKey]["fit_cistron_expression"]

def expressionConverge(
		sim_data,
		expression,
		concDict,
		doubling_time,
		Km=None,
		variable_elongation_transcription=True,
		variable_elongation_translation=False,
		disable_ribosome_capacity_fitting=False,
		disable_rnapoly_capacity_fitting=False,
		):
	"""
	Iteratively fits synthesis probabilities for RNA. Calculates initial
	expression based on gene expression data and makes adjustments to match
	physiological constraints for ribosome and RNAP counts. Relies on
	fitExpression() to converge

	Inputs
	------
	- expression (array of floats) - expression for each RNA, normalized to 1
	- concDict {metabolite (str): concentration (float with units of mol/volume)} -
	dictionary for concentrations of each metabolite with location tag
	- doubling_time (float with units of time) - doubling time
	- Km (array of floats with units of mol/volume) - Km for each RNA associated
	with RNases
	- disable_ribosome_capacity_fitting (bool) - if True, ribosome expression
	is not fit
	- disable_rnapoly_capacity_fitting (bool) - if True, RNA polymerase
	expression is not fit

	Requires
	--------
	- MAX_FITTING_ITERATIONS (int) - number of iterations to adjust expression
	before an exception is raised
	- FITNESS_THRESHOLD (float) - acceptable change from one iteration to break
	the fitting loop

	Returns
	--------
	- expression (array of floats) - adjusted expression for each RNA,
	normalized to 1
	- synthProb (array of floats) - synthesis probability for each RNA which
	accounts for expression and degradation rate, normalized to 1
	- avgCellDryMassInit (float with units of mass) - expected initial dry cell mass
	- fitAvgSolubleTargetMolMass (float with units of mass) - the adjusted dry mass
	of the soluble fraction of a cell
	- bulkContainer (BulkObjectsContainer object) - expected counts for
	bulk molecules based on expression
	"""

	if VERBOSE > 0:
		print("Fitting RNA synthesis probabilities.")

	for iteration in range(MAX_FITTING_ITERATIONS):
		if VERBOSE > 1:
			print('Iteration: {}'.format(iteration))

		initialExpression = expression.copy()
		expression = setInitialRnaExpression(sim_data, expression, doubling_time)
		bulkContainer = createBulkContainer(sim_data, expression, doubling_time)
		avgCellDryMassInit, fitAvgSolubleTargetMolMass = rescaleMassForSolubleMetabolites(sim_data, bulkContainer, concDict, doubling_time)

		if not disable_rnapoly_capacity_fitting:
			setRNAPCountsConstrainedByPhysiology(
				sim_data,
				bulkContainer,
				doubling_time,
				avgCellDryMassInit,
				variable_elongation_transcription,
				Km)

		if not disable_ribosome_capacity_fitting:
			setRibosomeCountsConstrainedByPhysiology(
				sim_data,
				bulkContainer,
				doubling_time,
				variable_elongation_translation)

		# Normalize expression and write out changes
		expression, synthProb, fit_cistron_expression, cistron_expression_res = fitExpression(
			sim_data, bulkContainer, doubling_time, avgCellDryMassInit, Km)

		degreeOfFit = np.sqrt(np.mean(np.square(initialExpression - expression)))

		if VERBOSE > 1:
			print('degree of fit: {}'.format(degreeOfFit))
			print(f'Average cistron expression residuals: {np.linalg.norm(cistron_expression_res)}')

		if degreeOfFit < FITNESS_THRESHOLD:
			break

	else:
		raise Exception("Fitting did not converge")

	return expression, synthProb, fit_cistron_expression, avgCellDryMassInit, fitAvgSolubleTargetMolMass, bulkContainer, concDict

def fitCondition(sim_data, spec, condition):
	"""
	Takes a given condition and returns the predicted bulk average, bulk deviation,
	protein monomer average, protein monomer deviation, and amino acid supply to
	translation. This relies on calculateBulkDistributions and calculateTranslationSupply.

	Inputs
	------
	- condition (str) - condition to fit (eg 'CPLX0-7705__active')
	- spec {property (str): property values} - cell specifications for the given condition.
	This function uses the specs "expression", "concDict", "avgCellDryMassInit",
	and "doubling_time"

	Returns
	--------
	- A dictionary {condition (str): spec (dict)} with the updated spec dictionary
	with the following values updated:
		- bulkAverageContainer (BulkObjectsContainer object) - the mean of the bulk counts
		- bulkDeviationContainer (BulkObjectsContainer object) - the standard deviation of the bulk counts
		- proteinMonomerAverageContainer (BulkObjectsContainer object) - the mean of the protein monomer counts
		- proteinMonomerDeviationContainer (BulkObjectsContainer object) - the standard deviation of the protein
		monomer counts
		- translation_aa_supply (array with units of mol/(mass.time)) - the supply rates
		for each amino acid to translation
	"""

	if VERBOSE > 0:
		print("Fitting condition {}".format(condition))

	# Find bulk and protein distributions
	bulkAverageContainer, bulkDeviationContainer, proteinMonomerAverageContainer, proteinMonomerDeviationContainer = calculateBulkDistributions(
		sim_data,
		spec["expression"],
		spec["concDict"],
		spec["avgCellDryMassInit"],
		spec["doubling_time"],
		)
	spec["bulkAverageContainer"] = bulkAverageContainer
	spec["bulkDeviationContainer"] = bulkDeviationContainer
	spec["proteinMonomerAverageContainer"] = proteinMonomerAverageContainer
	spec["proteinMonomerDeviationContainer"] = proteinMonomerDeviationContainer

	# Find the supply rates of amino acids to translation given doubling time
	spec["translation_aa_supply"] = calculateTranslationSupply(
											sim_data,
											spec["doubling_time"],
											spec["proteinMonomerAverageContainer"],
											spec["avgCellDryMassInit"],
											)

	return {condition: spec}

def calculateTranslationSupply(sim_data, doubling_time, bulkContainer, avgCellDryMassInit):
	"""
	Returns the supply rates of all amino acids to translation given the desired
	doubling time. This creates a limit on the polypeptide elongation process,
	and thus on growth. The amino acid supply rate is found by calculating the
	concentration of amino acids per gram dry cell weight and multiplying by the
	loss to dilution given doubling time.

	Inputs
	------
	- doubling_time (float with units of time) - measured doubling times given the condition
	- bulkContainer (BulkObjectsContainer object) - a container that tracks the counts of all bulk molecules
	- avgCellDryMassInit (float with units of mass) - the average initial cell dry mass

	Notes
	-----
	- The supply of amino acids should not be based on a desired doubling time,
	but should come from a more mechanistic basis. This would allow simulations
	of environmental shifts in which the doubling time is unknown.
	"""

	aaCounts = sim_data.process.translation.monomer_data['aa_counts'] # the counts of each amino acid required for each protein
	proteinCounts = bulkContainer.counts(sim_data.process.translation.monomer_data["id"]) # the counts of all proteins
	nAvogadro = sim_data.constants.n_avogadro

	molAAPerGDCW = (
			units.sum(
				aaCounts * np.tile(proteinCounts.reshape(-1, 1), (1, 21)),
				axis = 0
			) * (
				(1 / (units.aa * nAvogadro)) *
				(1 / avgCellDryMassInit)
			)
		)

	# Calculate required amino acid supply to translation to counter dilution
	translation_aa_supply = molAAPerGDCW * np.log(2) / doubling_time

	return translation_aa_supply

# Sub-fitting functions

def setTranslationEfficiencies(sim_data):
	"""
	This function's goal is to set translation efficiencies for a subset of metabolic proteins.
	It first gathers the index of the proteins it wants to modify, then changes the monomer
	translation efficiencies based on the adjustment that is specified.
	These adjustments were made so that the simulation could run.

	Requires
	--------
	- For each protein that needs to be modified, it takes in an adjustment factor.

	Modifies
	--------
	- This function modifies, for a subset of proteins, their translational efficiencies in sim_data.
	It takes their current efficiency and multiplies them by the factor specified in adjustments.
	"""

	for protein in sim_data.adjustments.translation_efficiencies_adjustments:
		idx = np.where(sim_data.process.translation.monomer_data["id"] == protein)[0]
		sim_data.process.translation.translation_efficiencies_by_monomer[idx] *= sim_data.adjustments.translation_efficiencies_adjustments[protein]


def set_balanced_translation_efficiencies(sim_data):
	"""
	Sets the translation efficiencies of a group of proteins to be equal to the
	mean value of all proteins within the group.

	Requires
	--------
	- List of proteins that should have balanced translation efficiencies.

	Modifies
	--------
	- Translation efficiencies of proteins within each specified group.
	"""
	monomer_id_to_index = {
		monomer['id'][:-3]: i for (i, monomer) in enumerate(
			sim_data.process.translation.monomer_data)
		}

	for proteins in sim_data.adjustments.balanced_translation_efficiencies:
		protein_indexes = np.array([monomer_id_to_index[m] for m in proteins])
		mean_trl_eff = sim_data.process.translation.translation_efficiencies_by_monomer[protein_indexes].mean()
		sim_data.process.translation.translation_efficiencies_by_monomer[protein_indexes] = mean_trl_eff


def setRNAExpression(sim_data):
	"""
	This function's goal is to set expression levels for a subset of RNAs.
	It first gathers the index of the RNA's it wants to modify, then changes
	the expression levels of those RNAs, within sim_data, based on the
	specified adjustment factor. If the specified ID is an RNA cistron, the
	expression levels of all RNA molecules containing the cistron are adjusted.

	Requires
	--------
	- For each RNA that needs to be modified, it takes in an adjustment factor.

	Modifies
	--------
	- This function modifies the basal RNA expression levels set in sim_data,
	for the chosen RNAs. It takes their current basal expression and multiplies
	them by the factor specified in adjustments.
	- After updating the basal expression levels for the given RNAs, the
	function normalizes all the basal expression levels.
	"""
	cistron_ids = set(sim_data.process.transcription.cistron_data['id'])
	rna_id_to_index = {
		rna_id[:-3]: i for (i, rna_id)
		in enumerate(sim_data.process.transcription.rna_data['id'])}

	rna_index_to_adjustment = {}

	for mol_id, adj_factor in sim_data.adjustments.rna_expression_adjustments.items():
		if mol_id in cistron_ids:
			# Find indexes of all RNAs containing the cistron
			rna_indexes = sim_data.process.transcription.cistron_id_to_rna_indexes(mol_id)
		elif mol_id in rna_id_to_index:
			rna_indexes = rna_id_to_index[mol_id]
		else:
			raise ValueError(
				f'Molecule ID {mol_id} not found in list of cistrons or transcription units.')

		# If multiple adjustments are made for the same RNA, take the maximum
		# adjustment factor
		for rna_index in rna_indexes:
			rna_index_to_adjustment[rna_index] = max(
				rna_index_to_adjustment.get(rna_index, 0), adj_factor)

	# Multiply all degradation rates with the specified adjustment factor
	for rna_index, adj_factor in rna_index_to_adjustment.items():
		sim_data.process.transcription.rna_expression["basal"][rna_index] *= adj_factor

	sim_data.process.transcription.rna_expression["basal"] /= sim_data.process.transcription.rna_expression["basal"].sum()

def setRNADegRates(sim_data):
	"""
	This function's goal is to adjust the degradation rates for a subset of
	metabolic RNA's. It first gathers the index of the RNA's it wants to modify,
	then changes the degradation rates of those RNAs. If the specified ID is
	that of an RNA cistron, the degradation rates of all RNA molecules
	containing the cistron are adjusted. (Note: since RNA concentrations are
	assumed to be in equilibrium, increasing the degradation rate increases the
	synthesis rates of these RNAs)

	Requires
	--------
	- For each RNA that needs to be modified, it takes in an adjustment factor

	Modifies
	--------
	- This function modifies the RNA degradation rates for the chosen RNAs in
	sim_data. It takes their current degradation rate and multiplies them by the
	factor specified in adjustments.
	"""
	cistron_id_to_index = {
		cistron_id: i for (i, cistron_id)
		in enumerate(sim_data.process.transcription.cistron_data['id'])
		}
	rna_id_to_index = {
		rna_id[:-3]: i for (i, rna_id)
		in enumerate(sim_data.process.transcription.rna_data['id'])}

	rna_index_to_adjustment = {}

	for mol_id, adj_factor in sim_data.adjustments.rna_deg_rates_adjustments.items():
		if mol_id in cistron_id_to_index:
			# Multiply the cistron degradation rate with the specified
			# adjustment factor (Note: these rates are not actually used by the
			# simulation but are still adjusted for bookkeeping purposes)
			cistron_index = cistron_id_to_index[mol_id]
			sim_data.process.transcription.cistron_data.struct_array["deg_rate"][cistron_index] *= sim_data.adjustments.rna_deg_rates_adjustments[mol_id]

			# Find indexes of all RNAs containing the cistron
			rna_indexes = sim_data.process.transcription.cistron_id_to_rna_indexes(mol_id)
		elif mol_id in rna_id_to_index:
			rna_indexes = rna_id_to_index[mol_id]
		else:
			raise ValueError(
				f'Molecule ID {mol_id} not found in list of cistrons or transcription units.')

		# If multiple adjustments are made for the same RNA, take the maximum
		# adjustment factor
		for rna_index in rna_indexes:
			rna_index_to_adjustment[rna_index] = max(
				rna_index_to_adjustment.get(rna_index, 0), adj_factor)

	# Multiply all degradation rates with the specified adjustment factor
	for rna_index, adj_factor in rna_index_to_adjustment.items():
		sim_data.process.transcription.rna_data.struct_array["deg_rate"][rna_index] *= adj_factor

def setProteinDegRates(sim_data):
	"""
	This function's goal is to set the degradation rates for a subset of proteins.
	It first gathers the index of the proteins it wants to modify, then changes the degradation
	rates of those proteins. These adjustments were made so that the simulation could run.

	Requires
	--------
	- For each protein that needs to be modified it take in an adjustment factor.

	Modifies
	--------
	- This function modifies the protein degradation rates for the chosen proteins in sim_data.
	It takes their current degradation rate and multiplies them by the factor specified in adjustments.
	"""

	for protein in sim_data.adjustments.protein_deg_rates_adjustments:
		idx = np.where(sim_data.process.translation.monomer_data["id"] == protein)[0]
		sim_data.process.translation.monomer_data.struct_array['deg_rate'][idx] *= sim_data.adjustments.protein_deg_rates_adjustments[protein]


def rescaleMassForSolubleMetabolites(sim_data, bulkMolCntr, concDict, doubling_time):
	"""
	Adjust the cell's mass to accomodate target small molecule concentrations.

	Inputs
	------
	- bulkMolCntr (BulkObjectsContainer object) - a container that tracks the counts of all bulk molecules
	- concDict (dict) - a dictionary of metabolite ID (string) : concentration (unit'd number, dimensions of concentration) pairs
	- doubling_time (float with units of time) - measured doubling times given the condition

	Requires
	--------
	- Cell mass fraction data at a given doubling time.
	- Average cell density.
	- The conversion factor for transforming from the size of an average cell to the size of a cell
	  immediately following division.
	- Avogadro's number.
	- Concentrations of small molecules (including both dry mass components and water).

	Modifies
	--------
	- Adds small molecule counts to bulkMolCntr.

	Returns
	-------
	- newAvgCellDryMassInit, the adjusted dry mass of a cell immediately following division.
	- fitAvgSolubleTargetMolMass, the adjusted dry mass of the soluble fraction of a cell
	"""

	avgCellFractionMass = sim_data.mass.get_component_masses(doubling_time)

	non_small_molecule_initial_cell_mass = (
		avgCellFractionMass["proteinMass"]
		+ avgCellFractionMass["rnaMass"]
		+ avgCellFractionMass["dnaMass"]
		) / sim_data.mass.avg_cell_to_initial_cell_conversion_factor

	molar_units = units.mol / units.L

	targetMoleculeIds = sorted(concDict)
	targetMoleculeConcentrations = molar_units * np.array([
		concDict[key].asNumber(molar_units) for key in targetMoleculeIds
		]) # Have to strip and replace units to obtain the proper array data type

	assert np.all(targetMoleculeConcentrations.asNumber(molar_units) > 0), 'Homeostatic dFBA objective requires non-zero (positive) concentrations'

	molecular_weights = sim_data.getter.get_masses(targetMoleculeIds)

	massesToAdd, countsToAdd = masses_and_counts_for_homeostatic_target(
		non_small_molecule_initial_cell_mass,
		targetMoleculeConcentrations,
		molecular_weights,
		sim_data.constants.cell_density,
		sim_data.constants.n_avogadro
		)

	bulkMolCntr.countsIs(
		countsToAdd,
		targetMoleculeIds
		)

	# Increase avgCellDryMassInit to match these numbers & rescale mass fractions
	smallMoleculetargetMoleculesDryMass = units.hstack((
		massesToAdd[:targetMoleculeIds.index('WATER[c]')],
		massesToAdd[targetMoleculeIds.index('WATER[c]') + 1:]
		)) # remove water since it's not part of the dry mass

	newAvgCellDryMassInit = non_small_molecule_initial_cell_mass + units.sum(smallMoleculetargetMoleculesDryMass)
	fitAvgSolubleTargetMolMass = units.sum(smallMoleculetargetMoleculesDryMass) * sim_data.mass.avg_cell_to_initial_cell_conversion_factor

	return newAvgCellDryMassInit, fitAvgSolubleTargetMolMass

def setInitialRnaExpression(sim_data, expression, doubling_time):
	"""
	Creates a container that with the initial count and ID of each RNA,
	calculated based on the mass fraction, molecular weight, and expression
	distribution of each RNA. For rRNA the counts are set based on mass, while
	for tRNA and mRNA the counts are set based on mass and relative abundance.
	Relies on the math function totalCountFromMassesAndRatios.

	Requires
	--------
	- Needs information from the knowledge base about the mass fraction,
	molecular weight, and distribution of each RNA species.

	Inputs
	------
	- expression (array of floats) - expression for each RNA, normalized to 1
	- doubling_time (float with units of time) - doubling time for condition

	Returns
	--------
	- expression (array of floats) - contains the adjusted RNA expression,
	normalized to 1

	Notes
	-----
	- Now rnaData["synthProb"] does not match "expression"

	"""
	# Load from sim_data
	n_avogadro = sim_data.constants.n_avogadro
	transcription = sim_data.process.transcription
	cistron_data = transcription.cistron_data
	rna_data = transcription.rna_data
	get_average_copy_number = sim_data.process.replication.get_average_copy_number
	rna_mw = rna_data['mw']
	rna_rRNA_mw = rna_data['rRNA_mw']
	rna_tRNA_mw = rna_data['tRNA_mw']
	rna_coord = rna_data['replication_coordinate']

	# Mask arrays for each RNA type
	is_rRNA = rna_data['is_rRNA']
	is_tRNA = rna_data['is_tRNA']
	is_mRNA = rna_data['is_mRNA']

	# Get list of RNA IDs for each type and tRNA cistron IDs
	all_RNA_ids = rna_data['id']
	ids_rRNA = all_RNA_ids[is_rRNA]
	ids_mRNA = all_RNA_ids[is_mRNA]
	ids_tRNA = all_RNA_ids[is_tRNA]
	ids_tRNA_cistrons = cistron_data['id'][cistron_data['is_tRNA']]

	# Get mass fractions of each RNA type for this condition
	initial_rna_mass = (
		sim_data.mass.get_component_masses(doubling_time)['rnaMass']
		/ sim_data.mass.avg_cell_to_initial_cell_conversion_factor
		)
	ppgpp = sim_data.growth_rate_parameters.get_ppGpp_conc(doubling_time)
	rna_fractions = transcription.get_rna_fractions(ppgpp)
	total_mass_rRNA = initial_rna_mass * rna_fractions['rRNA']
	total_mass_tRNA = initial_rna_mass * rna_fractions['tRNA']
	total_mass_mRNA = initial_rna_mass * rna_fractions['mRNA']

	# Get molecular weights of each RNA. For rRNAs/tRNAs, we only account for
	# the masses of the mature rRNAs/tRNAs within each RNA, since these RNAs are
	# almost instantly processed to yield the mature RNAs.
	individual_masses_rRNA = rna_rRNA_mw[is_rRNA] / n_avogadro
	individual_masses_tRNA = rna_tRNA_mw[is_tRNA] / n_avogadro
	individual_masses_mRNA = rna_mw[is_mRNA] / n_avogadro

	# Set rRNA TU expression assuming equal per-copy transcription
	# probabilities. The combined expression levels of each rRNA TU are assumed
	# to be proportional to their expected average copy numbers, which are
	# dependent on the doubling time and the chromosomal position.
	tau = doubling_time.asNumber(units.min)
	coord_rRNA = rna_coord[is_rRNA]
	n_avg_copy_rRNA = get_average_copy_number(tau, coord_rRNA)
	distribution_rRNA = normalize(n_avg_copy_rRNA)

	total_count_rRNA = totalCountFromMassesAndRatios(
		total_mass_rRNA,
		individual_masses_rRNA,
		distribution_rRNA
		)
	counts_rRNA = total_count_rRNA * distribution_rRNA

	# Get the total mass of tRNAs that are expressed from rRNA TUs and subtract
	# this mass from the total tRNA mass
	tRNA_masses_in_each_rRNA = rna_tRNA_mw[is_rRNA] / n_avogadro
	total_mass_tRNA_in_rRNAs = units.dot(counts_rRNA, tRNA_masses_in_each_rRNA)
	total_mass_tRNA -= total_mass_tRNA_in_rRNAs

	# Get tRNA cistron distribution (see Dong 1996), while setting values for
	# cistrons that are expressed from rRNAs to zero
	tRNA_distribution = sim_data.mass.get_trna_distribution(doubling_time)
	tRNA_id_to_dist = {
		trna_id: dist for (trna_id, dist)
		in zip(tRNA_distribution['id'], tRNA_distribution['molar_ratio_to_16SrRNA'])
		}
	distribution_tRNA_cistrons = np.zeros(len(ids_tRNA_cistrons))
	for i, tRNA_id in enumerate(ids_tRNA_cistrons):
		distribution_tRNA_cistrons[i] = tRNA_id_to_dist[tRNA_id]

	tRNA_expressed_from_rRNA_mask = transcription.cistron_tu_mapping_matrix.dot(is_rRNA)[
		cistron_data['is_tRNA']].astype(bool)
	distribution_tRNA_cistrons[tRNA_expressed_from_rRNA_mask] = 0
	distribution_tRNA_cistrons = normalize(distribution_tRNA_cistrons)

	# Approximate distribution of tRNA-including transcripts from tRNA cistron
	# distribution by using NNLS
	distribution_tRNA_including_transcripts, _ = transcription.fit_trna_expression(
		distribution_tRNA_cistrons)

	# Get distribution of tRNA-including transcripts that are not rRNAs
	is_hybrid = rna_data['is_rRNA'][rna_data['includes_tRNA']]
	distribution_tRNA = distribution_tRNA_including_transcripts[~is_hybrid]
	distribution_tRNA = normalize(distribution_tRNA)

	# Assign tRNA counts based on this distribution
	total_count_tRNA = totalCountFromMassesAndRatios(
		total_mass_tRNA,
		individual_masses_tRNA,
		distribution_tRNA
		)
	counts_tRNA = total_count_tRNA * distribution_tRNA

	# Assign mRNA counts based on mass and relative abundances (microarrays)
	distribution_mRNA = normalize(expression[is_mRNA])
	total_count_mRNA = totalCountFromMassesAndRatios(
		total_mass_mRNA,
		individual_masses_mRNA,
		distribution_mRNA
		)
	counts_mRNA = total_count_mRNA * distribution_mRNA

	# Set expression counts in container
	rna_expression_container = BulkObjectsContainer(
		all_RNA_ids, dtype=np.float64)
	rna_expression_container.countsIs(counts_rRNA, ids_rRNA)
	rna_expression_container.countsIs(counts_tRNA, ids_tRNA)
	rna_expression_container.countsIs(counts_mRNA, ids_mRNA)

	expression = normalize(rna_expression_container.counts())

	return expression


def totalCountIdDistributionProtein(sim_data, expression, doubling_time):
	"""
	Calculates the total counts of proteins from the relative expression of RNA,
	individual protein mass, and total protein mass. Relies on the math functions
	netLossRateFromDilutionAndDegradationProtein, proteinDistributionFrommRNA,
	totalCountFromMassesAndRatios.

	Inputs
	------
	- expression (array of floats) - relative frequency distribution of RNA expression
	- doubling_time (float with units of time) - measured doubling time given the condition

	Returns
	--------
	- total_count_protein (float) - total number of proteins
	- ids_protein (array of str) - name of each protein with location tag
	- distribution_protein (array of floats) - distribution for each protein,
	normalized to 1
	"""
	ids_protein = sim_data.process.translation.monomer_data["id"]
	total_mass_protein = sim_data.mass.get_component_masses(doubling_time)["proteinMass"] / sim_data.mass.avg_cell_to_initial_cell_conversion_factor
	individual_masses_protein = sim_data.process.translation.monomer_data["mw"] / sim_data.constants.n_avogadro

	mRNA_cistron_expression = sim_data.process.transcription.cistron_tu_mapping_matrix.dot(
		expression)[sim_data.process.transcription.cistron_data['is_mRNA']]
	distribution_transcripts_by_protein = normalize(
		sim_data.relation.monomer_to_mRNA_cistron_mapping().dot(mRNA_cistron_expression)
		)

	translation_efficiencies_by_protein = normalize(sim_data.process.translation.translation_efficiencies_by_monomer)
	degradationRates = sim_data.process.translation.monomer_data['deg_rate']

	# Find the net protein loss
	netLossRate_protein = netLossRateFromDilutionAndDegradationProtein(doubling_time, degradationRates)

	# Find the protein distribution
	distribution_protein = proteinDistributionFrommRNA(
		distribution_transcripts_by_protein,
		translation_efficiencies_by_protein,
		netLossRate_protein
		)

	# Find total protein counts
	total_count_protein = totalCountFromMassesAndRatios(
		total_mass_protein,
		individual_masses_protein,
		distribution_protein
		)

	return total_count_protein, ids_protein, distribution_protein

def totalCountIdDistributionRNA(sim_data, expression, doubling_time):
	"""
	Calculates the total counts of RNA from their relative expression,
	individual mass, and total RNA mass. Relies on the math function
	totalCountFromMassesAndRatios.

	Inputs
	------
	- expression (array of floats) - relative frequency distribution of RNA
		expression
	- doubling_time (float with units of time) - measured doubling time given
		the condition

	Returns
	--------
	- total_count_RNA (float) - total number of RNAs
	- ids_rnas (array of str) - name of each RNA with location tag
	- distribution_RNA (array of floats) - distribution for each RNA,
		normalized to 1
	"""
	transcription = sim_data.process.transcription
	ids_rnas = transcription.rna_data["id"]
	total_mass_RNA = (
		sim_data.mass.get_component_masses(doubling_time)["rnaMass"] /
		sim_data.mass.avg_cell_to_initial_cell_conversion_factor
		)
	mws = transcription.rna_data["mw"]
	# Use only the rRNA/tRNA mass for rRNA/tRNA transcription units
	is_rRNA = transcription.rna_data['is_rRNA']
	is_tRNA = transcription.rna_data['is_tRNA']
	mws[is_rRNA] = transcription.rna_data['rRNA_mw'][is_rRNA]
	mws[is_tRNA] = transcription.rna_data['tRNA_mw'][is_tRNA]
	individual_masses_RNA = mws / sim_data.constants.n_avogadro

	distribution_RNA = normalize(expression)

	total_count_RNA = totalCountFromMassesAndRatios(
		total_mass_RNA,
		individual_masses_RNA,
		distribution_RNA
		)

	return total_count_RNA, ids_rnas, distribution_RNA

def createBulkContainer(sim_data, expression, doubling_time):
	"""
	Creates a container that tracks the counts of all bulk molecules. Relies on
	totalCountIdDistributionRNA and totalCountIdDistributionProtein to set the
	counts and IDs of all RNAs and proteins.

	Inputs
	------
	- expression (array of floats) - relative frequency distribution of RNA expression
	- doubling_time (float with units of time) - measured doubling time given the condition

	Returns
	-------
	- bulkContainer (BulkObjectsContainer object) - a wrapper around a NumPy
	array that tracks the counts of bulk molecules
	"""

	total_count_RNA, ids_rnas, distribution_RNA = totalCountIdDistributionRNA(sim_data, expression, doubling_time)
	total_count_protein, ids_protein, distribution_protein = totalCountIdDistributionProtein(sim_data, expression, doubling_time)

	ids_molecules = sim_data.internal_state.bulk_molecules.bulk_data["id"]

	# Construct bulk container
	bulkContainer = BulkObjectsContainer(ids_molecules, dtype = np.float64)

	# Assign RNA counts based on mass and expression distribution
	counts_RNA = total_count_RNA * distribution_RNA
	bulkContainer.countsIs(counts_RNA, ids_rnas)

	# Assign protein counts based on mass and mRNA counts
	counts_protein = total_count_protein * distribution_protein
	bulkContainer.countsIs(counts_protein, ids_protein)

	return bulkContainer

def setRibosomeCountsConstrainedByPhysiology(
		sim_data,
		bulkContainer,
		doubling_time,
		variable_elongation_translation):
	"""
	Set counts of ribosomal protein subunits based on three constraints:
	(1) Expected protein distribution doubles in one cell cycle
	(2) Measured rRNA mass fractions
	(3) Expected ribosomal protein subunit counts based on RNA expression data

	Inputs
	------
	bulkContainer (BulkObjectsContainer object) - counts of bulk molecules
	doubling_time (float with units of time) - doubling time given the condition
	variable_elongation_translation (bool) - whether there is variable elongation for translation

	Modifies
	--------
	- counts of ribosomal protein subunits in bulkContainer
	"""

	active_fraction = sim_data.growth_rate_parameters.get_fraction_active_ribosome(doubling_time)

	# Get IDs and stoichiometry of ribosome subunits
	ribosome_30S_subunits = sim_data.process.complexation.get_monomers(
		sim_data.molecule_ids.s30_full_complex)['subunitIds']
	ribosome_50S_subunits = sim_data.process.complexation.get_monomers(
		sim_data.molecule_ids.s50_full_complex)['subunitIds']
	ribosome_30S_stoich = sim_data.process.complexation.get_monomers(
		sim_data.molecule_ids.s30_full_complex)['subunitStoich']
	ribosome_50S_stoich = sim_data.process.complexation.get_monomers(
		sim_data.molecule_ids.s50_full_complex)['subunitStoich']

	# Remove rRNA subunits from each array
	monomer_ids = set(sim_data.process.translation.monomer_data['id'])
	def remove_rRNA(subunit_ids, subunit_stoich):
		is_protein = np.array([
			(subunit_id in monomer_ids) for subunit_id in subunit_ids])
		return (subunit_ids[is_protein], subunit_stoich[is_protein])

	ribosome_30S_subunits, ribosome_30S_stoich = remove_rRNA(
		ribosome_30S_subunits, ribosome_30S_stoich)
	ribosome_50S_subunits, ribosome_50S_stoich = remove_rRNA(
		ribosome_50S_subunits, ribosome_50S_stoich)

	# -- CONSTRAINT 1: Expected protein distribution doubling -- #
	## Calculate minimium number of 30S and 50S subunits required in order to double our expected
	## protein distribution in one cell cycle
	proteinLengths = units.sum(sim_data.process.translation.monomer_data['aa_counts'], axis = 1)
	proteinDegradationRates = sim_data.process.translation.monomer_data['deg_rate']
	proteinCounts = bulkContainer.counts(sim_data.process.translation.monomer_data["id"])

	netLossRate_protein = netLossRateFromDilutionAndDegradationProtein(
		doubling_time,
		proteinDegradationRates,
		)

	elongation_rates = sim_data.process.translation.make_elongation_rates(
		None,
		sim_data.growth_rate_parameters.get_ribosome_elongation_rate(doubling_time).asNumber(units.aa / units.s),
		1,
		variable_elongation_translation)

	nRibosomesNeeded = np.ceil(
		calculateMinPolymerizingEnzymeByProductDistribution(
			proteinLengths,
			elongation_rates,
			netLossRate_protein,
			proteinCounts).asNumber(units.aa / units.s) / active_fraction
		)

	# Minimum number of ribosomes needed
	constraint1_ribosome30SCounts = (
		nRibosomesNeeded * ribosome_30S_stoich
		)

	constraint1_ribosome50SCounts = (
		nRibosomesNeeded * ribosome_50S_stoich
		)


	# -- CONSTRAINT 2: Measured rRNA mass fraction -- #
	# Get rRNA counts
	rna_data = sim_data.process.transcription.rna_data
	rRNA_tu_counts = bulkContainer.counts(rna_data["id"][rna_data['is_rRNA']])
	rRNA_cistron_counts = sim_data.process.transcription.rRNA_cistron_tu_mapping_matrix.dot(
		rRNA_tu_counts)
	rRNA_cistron_indexes = np.where(
		sim_data.process.transcription.cistron_data['is_rRNA'])[0]
	rRNA_23S_counts = rRNA_cistron_counts[
		sim_data.process.transcription.cistron_data['is_23S_rRNA'][rRNA_cistron_indexes]]
	rRNA_16S_counts = rRNA_cistron_counts[
		sim_data.process.transcription.cistron_data['is_16S_rRNA'][rRNA_cistron_indexes]]
	rRNA_5S_counts = rRNA_cistron_counts[
		sim_data.process.transcription.cistron_data['is_5S_rRNA'][rRNA_cistron_indexes]]

	## 16S rRNA is in the 30S subunit
	massFracPredicted_30SCount = rRNA_16S_counts.sum()
	## 23S and 5S rRNA are in the 50S subunit
	massFracPredicted_50SCount = min(rRNA_23S_counts.sum(), rRNA_5S_counts.sum())

	constraint2_ribosome30SCounts = massFracPredicted_30SCount * ribosome_30S_stoich
	constraint2_ribosome50SCounts = massFracPredicted_50SCount * ribosome_50S_stoich

	# -- CONSTRAINT 3: Expected ribosomal subunit counts based expression
	## Calculate fundamental ribosomal subunit count distribution based on RNA expression data
	## Already calculated and stored in bulkContainer
	ribosome30SCounts = bulkContainer.counts(ribosome_30S_subunits)
	ribosome50SCounts = bulkContainer.counts(ribosome_50S_subunits)

	# -- SET RIBOSOME FUNDAMENTAL SUBUNIT COUNTS TO MAXIMUM CONSTRAINT -- #
	constraint_names = np.array(["Insufficient to double protein counts", "Too small for mass fraction", "Current level OK"])
	rib30lims = np.array([nRibosomesNeeded, massFracPredicted_30SCount, (ribosome30SCounts / ribosome_30S_stoich).min()])
	rib50lims = np.array([nRibosomesNeeded, massFracPredicted_50SCount, (ribosome50SCounts / ribosome_50S_stoich).min()])
	if VERBOSE > 1:
		print('30S limit: {}'.format(constraint_names[np.where(rib30lims.max() == rib30lims)[0]][-1]))
		print('30S actual count: {}'.format((ribosome30SCounts / ribosome_30S_stoich).min()))
		print('30S count set to: {}'.format(rib30lims[np.where(rib30lims.max() == rib30lims)[0]][-1]))
		print('50S limit: {}'.format(constraint_names[np.where(rib50lims.max() == rib50lims)[0]][-1]))
		print('50S actual count: {}'.format((ribosome50SCounts / ribosome_50S_stoich).min()))
		print('50S count set to: {}'.format(rib50lims[np.where(rib50lims.max() == rib50lims)[0]][-1]))

	bulkContainer.countsIs(
		np.fmax(np.fmax(ribosome30SCounts, constraint1_ribosome30SCounts), constraint2_ribosome30SCounts),
		ribosome_30S_subunits
		)

	bulkContainer.countsIs(
		np.fmax(np.fmax(ribosome50SCounts, constraint1_ribosome50SCounts), constraint2_ribosome50SCounts),
		ribosome_50S_subunits
		)

def setRNAPCountsConstrainedByPhysiology(
		sim_data,
		bulkContainer,
		doubling_time,
		avgCellDryMassInit,
		variable_elongation_transcription,
		Km=None):
	"""
	Set counts of RNA polymerase based on two constraints:
	(1) Number of RNAP subunits required to maintain steady state of mRNAs
	(2) Expected RNAP subunit counts based on (mRNA) distribution recorded in
		bulkContainer

	Inputs
	------
	- bulkContainer (BulkObjectsContainer object) - counts of bulk molecules
	- doubling_time (float with units of time) - doubling time given the condition
	- avgCellDryMassInit (float with units of mass) - expected initial dry cell mass
	- Km (array of floats with units of mol/volume) - Km for each RNA associated
	with RNases

	Modifies
	--------
	- bulkContainer (BulkObjectsContainer object) - the counts of RNA polymerase
	subunits are set according to Constraint 1

	Notes
	-----
	- Constraint 2 is not being used -- see final line of this function.
	"""

	# -- CONSTRAINT 1: Expected RNA distribution doubling -- #
	rnaLengths = units.sum(sim_data.process.transcription.rna_data['counts_ACGU'], axis = 1)

	rnaLossRate = None

	if Km is None:
		# RNA loss rate is in units of counts/time, and computed by summing the
		# contributions of degradation and dilution.
		rnaLossRate = netLossRateFromDilutionAndDegradationRNALinear(
			doubling_time,
			sim_data.process.transcription.rna_data['deg_rate'],
			bulkContainer.counts(sim_data.process.transcription.rna_data['id'])
		)
	else:
		# Get constants to compute countsToMolar factor
		cellDensity = sim_data.constants.cell_density
		cellVolume = avgCellDryMassInit / cellDensity / sim_data.mass.cell_dry_mass_fraction
		countsToMolar = 1 / (sim_data.constants.n_avogadro * cellVolume)

		# Gompute input arguments for netLossRateFromDilutionAndDegradationRNA()
		rnaConc = countsToMolar * bulkContainer.counts(sim_data.process.transcription.rna_data['id'])
		endoRNaseConc = countsToMolar * bulkContainer.counts(sim_data.process.rna_decay.endoRNase_ids)
		kcatEndoRNase = sim_data.process.rna_decay.kcats
		totalEndoRnaseCapacity = units.sum(endoRNaseConc * kcatEndoRNase)

		# RNA loss rate is in units of counts/time, and computed by accounting
		# for the competitive inhibition of RNase by other RNA targets.
		rnaLossRate = netLossRateFromDilutionAndDegradationRNA(
			doubling_time,
			(1 / countsToMolar) * totalEndoRnaseCapacity,
			Km,
			rnaConc,
			countsToMolar,
			)

	# Compute number of RNA polymerases required to maintain steady state of mRNA
	elongation_rates = sim_data.process.transcription.make_elongation_rates(
		None,
		sim_data.growth_rate_parameters.get_rnap_elongation_rate(doubling_time).asNumber(units.nt / units.s),
		1,
		variable_elongation_transcription)

	nActiveRnapNeeded = calculateMinPolymerizingEnzymeByProductDistributionRNA(
		rnaLengths,
		elongation_rates,
		rnaLossRate).asNumber(units.nt / units.s)

	nRnapsNeeded = np.ceil(
		nActiveRnapNeeded / sim_data.growth_rate_parameters.get_fraction_active_rnap(doubling_time))

	# Convert nRnapsNeeded to the number of RNA polymerase subunits required
	rnapIds = sim_data.process.complexation.get_monomers(sim_data.molecule_ids.full_RNAP)['subunitIds']
	rnapStoich = sim_data.process.complexation.get_monomers(sim_data.molecule_ids.full_RNAP)['subunitStoich']
	minRnapSubunitCounts = nRnapsNeeded * rnapStoich

	# -- CONSTRAINT 2: Expected RNAP subunit counts based on distribution -- #
	rnapCounts = bulkContainer.counts(rnapIds)

	## -- SET RNAP COUNTS TO MAXIMUM CONSTRAINTS -- #
	constraint_names = np.array(["Current level OK", "Insufficient to double RNA distribution"])
	rnapLims = np.array([(rnapCounts / rnapStoich).min(), (minRnapSubunitCounts / rnapStoich).min()])
	if VERBOSE > 1:
		print('rnap limit: {}'.format(constraint_names[np.where(rnapLims.max() == rnapLims)[0]][0]))
		print('rnap actual count: {}'.format((rnapCounts / rnapStoich).min()))
		print('rnap counts set to: {}'.format(rnapLims[np.where(rnapLims.max() == rnapLims)[0]][0]))

	if np.any(minRnapSubunitCounts < 0):
		raise ValueError('RNAP protein counts must be positive.')

	bulkContainer.countsIs(minRnapSubunitCounts, rnapIds)

def fitExpression(sim_data, bulkContainer, doubling_time, avgCellDryMassInit, Km=None):
	"""
	Determines expression and synthesis probabilities for RNA molecules to fit
	protein levels and RNA degradation rates. Assumes a steady state analysis
	where the RNA synthesis probability will be the same as the degradation rate.
	If no Km is given, then RNA degradation is assumed to be linear otherwise
	degradation is calculated based on saturation with RNases.

	Inputs
	------
	- bulkContainer (BulkObjectsContainer object) - expected counts for
	bulk molecules based on expression
	- doubling_time (float with units of time) - doubling time
	- avgCellDryMassInit (float with units of mass) - expected initial dry cell mass
	- Km (array of floats with units of mol/volume) - Km for each RNA associated
	with RNases

	Modifies
	--------
	- bulkContainer counts of RNA and proteins

	Returns
	--------
	- expression (array of floats) - adjusted expression for each RNA,
	normalized to 1
	- synth_prob (array of floats) - synthesis probability for each RNA which
	accounts for expression and degradation rate, normalized to 1
	- fit_cistron_expression (array of floats) - target expression levels of
	each cistron (gene) used to calculate RNA expression levels
	- cistron_expression_res (array of floats) - the residuals of the NNLS
	problem solved to calculate RNA expression levels

	Notes
	-----
	- TODO - sets bulkContainer counts and returns values - change to only return values
	"""
	# Load required parameters
	transcription = sim_data.process.transcription
	translation = sim_data.process.translation
	translation_efficiencies_by_protein = normalize(
		translation.translation_efficiencies_by_monomer)
	degradation_rates_protein = translation.monomer_data['deg_rate']
	net_loss_rate_protein = netLossRateFromDilutionAndDegradationProtein(
		doubling_time, degradation_rates_protein)
	avg_cell_fraction_mass = sim_data.mass.get_component_masses(doubling_time)
	total_mass_RNA = avg_cell_fraction_mass["rnaMass"] / sim_data.mass.avg_cell_to_initial_cell_conversion_factor
	cistron_tu_mapping_matrix = transcription.cistron_tu_mapping_matrix

	# Calculate current expression fraction of mRNA transcription units
	view_RNA = bulkContainer.countsView(transcription.rna_data["id"])
	rna_expression_container = BulkObjectsContainer(
		list(transcription.rna_data["id"]), dtype = np.dtype("float64"))
	rna_expression_container.countsIs(normalize(view_RNA.counts()))

	mRNA_tu_expression_view = rna_expression_container.countsView(
		transcription.rna_data["id"][transcription.rna_data['is_mRNA']])
	mRNA_tu_expression_frac = np.sum(mRNA_tu_expression_view.counts())

	# Calculate current expression levels of each cistron given the RNA
	# expression levels
	fit_cistron_expression = normalize(
		cistron_tu_mapping_matrix.dot(view_RNA.counts()))
	mRNA_cistron_expression_frac = fit_cistron_expression[
		transcription.cistron_data['is_mRNA']].sum()

	# Calculate required mRNA expression from monomer counts
	counts_protein = bulkContainer.counts(
		translation.monomer_data["id"])
	mRNA_cistron_distribution_per_protein = mRNADistributionFromProtein(
		normalize(counts_protein),
		translation_efficiencies_by_protein,
		net_loss_rate_protein)

	mRNA_cistron_distribution = normalize(
		sim_data.relation.monomer_to_mRNA_cistron_mapping().T.dot(
			mRNA_cistron_distribution_per_protein)
		)

	# Replace mRNA cistron expression with values calculated from monomer counts
	fit_cistron_expression[transcription.cistron_data['is_mRNA']] = (
		mRNA_cistron_expression_frac * mRNA_cistron_distribution)

	# Use least squares to calculate expression of transcription units required
	# to generate the given cistron expression levels and the residuals for
	# the expression of each cistron
	fit_tu_expression, cistron_expression_res = transcription.fit_rna_expression(
		fit_cistron_expression)
	fit_mRNA_tu_expression = fit_tu_expression[transcription.rna_data['is_mRNA']]

	mRNA_tu_expression_view.countsIs(
		mRNA_tu_expression_frac * normalize(fit_mRNA_tu_expression))
	expression = normalize(rna_expression_container.counts())

	# Set number of RNAs based on expression we just set
	mws = transcription.rna_data["mw"]

	# Use only the rRNA/tRNA mass for rRNA/tRNA transcription units
	is_rRNA = transcription.rna_data['is_rRNA']
	is_tRNA = transcription.rna_data['is_tRNA']
	mws[is_rRNA] = transcription.rna_data['rRNA_mw'][is_rRNA]
	mws[is_tRNA] = transcription.rna_data['tRNA_mw'][is_tRNA]

	n_rnas = totalCountFromMassesAndRatios(
		total_mass_RNA, mws / sim_data.constants.n_avogadro, expression)
	view_RNA.countsIs(n_rnas * expression)

	if Km is None:
		rnaLossRate = netLossRateFromDilutionAndDegradationRNALinear(
			doubling_time,
			transcription.rna_data['deg_rate'],
			bulkContainer.counts(transcription.rna_data['id'])
		)
	else:
		# Get constants to compute countsToMolar factor
		cellDensity = sim_data.constants.cell_density
		dryMassFraction = sim_data.mass.cell_dry_mass_fraction
		cellVolume = avgCellDryMassInit / cellDensity / dryMassFraction
		countsToMolar = 1 / (sim_data.constants.n_avogadro * cellVolume)

		endoRNaseConc = countsToMolar * bulkContainer.counts(sim_data.process.rna_decay.endoRNase_ids)
		kcatEndoRNase = sim_data.process.rna_decay.kcats
		totalEndoRnaseCapacity = units.sum(endoRNaseConc * kcatEndoRNase)

		rnaLossRate = netLossRateFromDilutionAndDegradationRNA(
			doubling_time,
			(1 / countsToMolar) * totalEndoRnaseCapacity,
			Km,
			countsToMolar * view_RNA.counts(),
			countsToMolar,
		)

	synth_prob = normalize(rnaLossRate.asNumber(1 / units.min))

	return expression, synth_prob, fit_cistron_expression, cistron_expression_res

def fitMaintenanceCosts(sim_data, bulkContainer):
	"""
	Fits the growth-associated maintenance (GAM) cost associated with metabolism.

	The energetic costs associated with growth have been estimated utilizing flux-balance analysis
	and are used with FBA to obtain accurate growth predictions.  In the whole-cell model, some of
	these costs are explicitly associated with the energetic costs of translation, a biomass
	assembly process.  Consequently we must estimate the amount of energy utilized by translation
	per unit of biomass (i.e. dry mass) produced, and subtract that quantity from reported GAM to
	acquire the modified GAM that we use in the metabolic submodel.

	Requires
	--------
	- amino acid counts associated with protein monomers
	- average initial dry mass
	- energetic (GTP) cost of translation (per amino acid polymerized)
	- observed growth-associated maintenance (GAM)
	In dimensions of ATP or ATP equivalents consumed per biomass

	Modifies
	--------
	- the "dark" ATP, i.e. the modified GAM

	Notes
	-----
	As more non-metabolic submodels account for energetic costs, this function should be extended
	to subtract those costs off the observed GAM.

	There also exists, in contrast, non-growth-associated-maintenance (NGAM), which is relative to
	total biomass rather than the biomass accumulation rate.  As the name would imply, this
	accounts for the energetic costs of maintaining the existing biomass.  It is also accounted for
	in the metabolic submodel.

	TODO (John): Rewrite as a true function.
	"""


	aaCounts = sim_data.process.translation.monomer_data['aa_counts']
	proteinCounts = bulkContainer.counts(sim_data.process.translation.monomer_data["id"])
	nAvogadro = sim_data.constants.n_avogadro
	avgCellDryMassInit = sim_data.mass.avg_cell_dry_mass_init
	gtpPerTranslation = sim_data.constants.gtp_per_translation
	atp_per_charge = 2  # ATP -> AMP is explicitly used in charging reactions so can remove from GAM

	# GTPs used for translation (recycled, not incorporated into biomass)
	aaMmolPerGDCW = (
			units.sum(
				aaCounts * np.tile(proteinCounts.reshape(-1, 1), (1, 21)),
				axis = 0
			) * (
				(1 / (units.aa * nAvogadro)) *
				(1 / avgCellDryMassInit)
			)
		)

	aasUsedOverCellCycle = units.sum(aaMmolPerGDCW)
	explicit_mmol_maintenance_per_gdcw = (atp_per_charge + gtpPerTranslation) * aasUsedOverCellCycle

	darkATP = ( # This has everything we can't account for
		sim_data.constants.growth_associated_maintenance -
		explicit_mmol_maintenance_per_gdcw
		)

	# We do not want to create energy with growth by having a negative darkATP
	# value. GAM measurements have some error so it's possible explicit
	# accounting could be more accurate or the GAM value used is too low which
	# would lead to a negative value. Easy fix is setting darkATP = 0 if this
	# error is raised.
	if darkATP.asNumber() < 0:
		raise ValueError('GAM has been adjusted too low. Explicit energy accounting should not exceed GAM.'
			' Consider setting darkATP to 0 if energy corrections are accurate.')

	sim_data.constants.darkATP = darkATP

def calculateBulkDistributions(sim_data, expression, concDict, avgCellDryMassInit, doubling_time):
	"""
	Finds a distribution of copy numbers for macromolecules. While RNA and protein
	expression can be approximated using well-described statistical	distributions,
	complexes require absolute copy numbers. To get these distributions, this
	function instantiates many cells with a reduced set of molecules, forms complexes,
	and iterates through equilibrium and two-component system processes until
	metabolite counts reach a steady-state. It then computes the resulting
	statistical distributions.

	Requires
	--------
	- N_SEEDS (int) - the number of instantiated cells

	Inputs
	------
	- expression (array of floats) - expression for each RNA, normalized to 1
	- concDict {metabolite (str): concentration (float with units of mol/volume)} -
	dictionary for concentrations of each metabolite with location tag
	- avgCellDryMassInit (float with units of mass) - initial dry cell mass
	- doubling_time (float with units of time) - doubling time for condition

	Returns
	--------
	- bulkAverageContainer (BulkObjectsContainer object) - the mean of the bulk counts
	- bulkDeviationContainer (BulkObjectsContainer object) - the standard deviation of the bulk counts
	- proteinMonomerAverageContainer (BulkObjectsContainer object) - the mean of the protein monomer counts
	- proteinMonomerDeviationContainer (BulkObjectsContainer object) - the standard deviation of the protein
	"""

	# Ids
	totalCount_RNA, ids_rnas, distribution_RNA = totalCountIdDistributionRNA(sim_data, expression, doubling_time)
	totalCount_protein, ids_protein, distribution_protein = totalCountIdDistributionProtein(sim_data, expression, doubling_time)
	ids_complex = sim_data.process.complexation.molecule_names
	ids_equilibrium = sim_data.process.equilibrium.molecule_names
	ids_twoComponentSystem = sim_data.process.two_component_system.molecule_names
	ids_metabolites = sorted(concDict)
	conc_metabolites = (units.mol / units.L) * np.array([concDict[key].asNumber(units.mol / units.L) for key in ids_metabolites])
	allMoleculesIDs = sorted(
		set(ids_rnas) | set(ids_protein) | set(ids_complex) | set(ids_equilibrium) | set(ids_twoComponentSystem) | set(ids_metabolites)
		)

	# Data for complexation
	complexationStoichMatrix = sim_data.process.complexation.stoich_matrix().astype(np.int64, order ="F")
	# Data for equilibrium binding
	# equilibriumDerivatives = sim_data.process.equilibrium.derivatives
	# equilibriumDerivativesJacobian = sim_data.process.equilibrium.derivativesJacobian

	# Data for metabolites
	cellDensity = sim_data.constants.cell_density
	cellVolume = avgCellDryMassInit / cellDensity / sim_data.mass.cell_dry_mass_fraction

	# Construct bulk container

	# We want to know something about the distribution of the copy numbers of
	# macromolecules in the cell.  While RNA and protein expression can be
	# approximated using well-described statistical distributions, we need
	# absolute copy numbers to form complexes.  To get a distribution, we must
	# instantiate many cells, form complexes, and finally compute the
	# statistics we will use in the fitting operations.

	bulkContainer = BulkObjectsContainer(sim_data.internal_state.bulk_molecules.bulk_data['id'])
	rnaView = bulkContainer.countsView(ids_rnas)
	proteinView = bulkContainer.countsView(ids_protein)
	complexationMoleculesView = bulkContainer.countsView(ids_complex)
	equilibriumMoleculesView = bulkContainer.countsView(ids_equilibrium)
	twoComponentSystemMoleculesView = bulkContainer.countsView(ids_twoComponentSystem)
	metabolitesView = bulkContainer.countsView(ids_metabolites)
	allMoleculesView = bulkContainer.countsView(allMoleculesIDs)

	allMoleculeCounts = np.empty((N_SEEDS, allMoleculesView.counts().size), np.int64)
	proteinMonomerCounts = np.empty((N_SEEDS, proteinView.counts().size), np.int64)

	if VERBOSE > 1:
		print("Bulk distribution seed:")

	# Instantiate cells to find average copy numbers of macromolecules
	for seed in range(N_SEEDS):
		if VERBOSE > 1:
			print('seed = {}'.format(seed))

		allMoleculesView.countsIs(0)

		rnaView.countsIs(totalCount_RNA * distribution_RNA)

		proteinView.countsIs(totalCount_protein * distribution_protein)
		proteinMonomerCounts[seed, :] = proteinView.counts()
		complexationMoleculeCounts = complexationMoleculesView.counts()

		# Form complexes
		time_step = 2**31 # don't stop until all complexes are formed.
		complexation_rates = sim_data.process.complexation.rates
		system = StochasticSystem(complexationStoichMatrix.T, random_seed=seed)
		complexation_result = system.evolve(
			time_step, complexationMoleculeCounts, complexation_rates)

		updatedCompMoleculeCounts = complexation_result['outcome']
		complexationMoleculesView.countsIs(updatedCompMoleculeCounts)

		metDiffs = np.inf * np.ones_like(metabolitesView.counts())
		nIters = 0

		# Iterate processes until metabolites converge to a steady-state
		while np.linalg.norm(metDiffs, np.inf) > 1:
			random_state = np.random.RandomState(seed)
			metCounts = conc_metabolites * cellVolume * sim_data.constants.n_avogadro
			metCounts.normalize()
			metCounts.checkNoUnit()
			metabolitesView.countsIs(
				metCounts.asNumber().round()
				)

			# Find reaction fluxes from equilibrium process
			# Do not use jit to avoid compiling time (especially when running
			# in parallel since sim_data needs to be pickled and reconstructed
			# each time)
			rxnFluxes, _ = sim_data.process.equilibrium.fluxes_and_molecules_to_SS(
				equilibriumMoleculesView.counts(),
				cellVolume.asNumber(units.L),
				sim_data.constants.n_avogadro.asNumber(1 / units.mol),
				random_state, jit=False,
				)
			equilibriumMoleculesView.countsInc(
				np.dot(sim_data.process.equilibrium.stoich_matrix().astype(np.int64), rxnFluxes)
				)
			assert np.all(equilibriumMoleculesView.counts() >= 0)

			# Find changes from two component system
			_, moleculeCountChanges = sim_data.process.two_component_system.molecules_to_ss(
				twoComponentSystemMoleculesView.counts(),
				cellVolume.asNumber(units.L),
				sim_data.constants.n_avogadro.asNumber(1 / units.mmol)
				)

			twoComponentSystemMoleculesView.countsInc(moleculeCountChanges)

			metDiffs = metabolitesView.counts() - metCounts.asNumber().round()

			nIters += 1
			if nIters > 100:
				raise Exception("Equilibrium reactions are not converging!")

		allMoleculeCounts[seed, :] = allMoleculesView.counts()

	# Update counts in bulk objects container
	bulkAverageContainer = BulkObjectsContainer(sim_data.internal_state.bulk_molecules.bulk_data['id'], np.float64)
	bulkDeviationContainer = BulkObjectsContainer(sim_data.internal_state.bulk_molecules.bulk_data['id'], np.float64)
	proteinMonomerAverageContainer = BulkObjectsContainer(sim_data.process.translation.monomer_data["id"], np.float64)
	proteinMonomerDeviationContainer = BulkObjectsContainer(sim_data.process.translation.monomer_data["id"], np.float64)

	bulkAverageContainer.countsIs(allMoleculeCounts.mean(0), allMoleculesIDs)
	bulkDeviationContainer.countsIs(allMoleculeCounts.std(0), allMoleculesIDs)
	proteinMonomerAverageContainer.countsIs(proteinMonomerCounts.mean(0), sim_data.process.translation.monomer_data["id"])
	proteinMonomerDeviationContainer.countsIs(proteinMonomerCounts.std(0), sim_data.process.translation.monomer_data["id"])

	return bulkAverageContainer, bulkDeviationContainer, proteinMonomerAverageContainer, proteinMonomerDeviationContainer


# Math functions

def totalCountFromMassesAndRatios(totalMass, individualMasses, distribution):
	"""
	Function to determine the expected total counts for a group of molecules
	in order to achieve a total mass with a given distribution of individual
	molecules.

	Math:
		Total mass = dot(mass, count)

		Fraction of i:
		f = count / Total counts

		Substituting:
		Total mass = dot(mass, f * Total counts)
		Total mass = Total counts * dot(mass, f)

		Total counts = Total mass / dot(mass, f)

	Requires
	--------
	- totalMass (float with mass units): total mass of the group of molecules
	- individualMasses (array of floats with mass units): mass for individual
	molecules in the group
	- distribution (array of floats): distribution of individual molecules,
	normalized to 1

	Returns
	--------
	- counts (float): total counts (does not need to be a whole number)
	"""

	assert np.allclose(np.sum(distribution), 1)
	counts = 1 / units.dot(individualMasses, distribution) * totalMass
	return units.strip_empty_units(counts)

def proteinDistributionFrommRNA(distribution_mRNA, translation_efficiencies, netLossRate):
	"""
	dP_i / dt = k * M_i * e_i - P_i * Loss_i

	At steady state:
	P_i = k * M_i * e_i / Loss_i

	Fraction of mRNA for ith gene is defined as:
	f_i = M_i / M_total

	Substituting in:
	P_i = k * f_i * e_i * M_total / Loss_i

	Normalizing P_i by summing over all i cancels out k and M_total
	assuming constant translation rate.

	Inputs
	------
	- distribution_mRNA (array of floats) - distribution for each mRNA,
	normalized to 1
	- translation_efficiencies (array of floats) - translational efficiency for each mRNA,
	normalized to 1
	- netLossRate (array of floats with units of 1/time) - rate of loss for each protein

	Returns
	-------
	- array of floats for the distribution of each protein, normalized to 1
	"""

	assert np.allclose(np.sum(distribution_mRNA), 1)
	assert np.allclose(np.sum(translation_efficiencies), 1)
	distributionUnnormed = 1 / netLossRate * distribution_mRNA * translation_efficiencies
	distributionNormed = distributionUnnormed / units.sum(distributionUnnormed)
	distributionNormed.normalize()
	distributionNormed.checkNoUnit()

	return distributionNormed.asNumber()

def mRNADistributionFromProtein(distribution_protein, translation_efficiencies, netLossRate):
	"""
	dP_i / dt = k * M_i * e_i - P_i * Loss_i

	At steady state:
	M_i = Loss_i * P_i / (k * e_i)

	Fraction of protein for ith gene is defined as:
	f_i = P_i / P_total

	Substituting in:
	M_i = Loss_i * f_i * P_total / (k * e_i)

	Normalizing M_i by summing over all i cancels out k and P_total
	assuming a constant translation rate.

	Inputs
	------
	- distribution_protein (array of floats) - distribution for each protein,
	normalized to 1
	- translation_efficiencies (array of floats) - translational efficiency for each mRNA,
	normalized to 1
	- netLossRate (array of floats with units of 1/time) - rate of loss for each protein

	Returns
	-------
	- array of floats for the distribution of each mRNA, normalized to 1
	"""

	assert np.allclose(np.sum(distribution_protein), 1)
	distributionUnnormed = netLossRate * distribution_protein / translation_efficiencies
	distributionNormed = distributionUnnormed / units.sum(distributionUnnormed)
	distributionNormed.normalize()
	distributionNormed.checkNoUnit()

	return distributionNormed.asNumber()

def calculateMinPolymerizingEnzymeByProductDistribution(
		productLengths, elongationRates, netLossRate, productCounts):
	"""
	Compute the number of ribosomes required to maintain steady state.

	dP/dt = production rate - loss rate
	dP/dt = e_r * (1/L) * R - (k_loss * P)

	At steady state: dP/dt = 0
	R = sum over i ((L_i / e_r) * k_loss_i * P_i)

	Multiplying both sides by volume gives an equation in terms of counts.

	P = protein concentration
	e_r = polypeptide elongation rate per ribosome
	L = protein length
	R = ribosome concentration
	k_loss = net protein loss rate
	i = ith protein

	Inputs
	------
	- productLengths (array of ints with units of amino_acids) - L, protein lengths
	- elongationRates (array of ints with units of amino_acid/time) e_r, polypeptide elongation rate
	- netLossRate (array of floats with units of 1/time) - k_loss, protein loss rate
	- productCounts (array of floats) - P, protein counts

	Returns
	--------
	- float with dimensionless units for the number of ribosomes required to
	maintain steady state
	"""

	nPolymerizingEnzymeNeeded = units.sum(
		productLengths / elongationRates
		* netLossRate
		* productCounts)
	return nPolymerizingEnzymeNeeded

def calculateMinPolymerizingEnzymeByProductDistributionRNA(productLengths, elongationRates, netLossRate):
	"""
	Compute the number of RNA polymerases required to maintain steady state of mRNA.

	dR/dt = production rate - loss rate
	dR/dt = e_r * (1/L) * RNAp - k_loss

	At steady state: dR/dt = 0
	RNAp = sum over i ((L_i / e_r) * k_loss_i)

	Multiplying both sides by volume gives an equation in terms of counts.

	R = mRNA transcript concentration
	e_r = transcript elongation rate per RNAp
	L = transcript length
	RNAp = RNAp concentration
	k_loss = net transcript loss rate (unit: concentration / time)
	i = ith transcript

	Inputs
	------
	- productLengths (array of ints with units of nucleotides) - L, transcript lengths
	- elongationRates (array of ints with units of nucleotide/time) - e_r, transcript elongation rate
	- netLossRate (array of floats with units of 1/time) - k_loss, transcript loss rate

	Returns
	-------
	- float with dimensionless units for the number of RNA polymerases required to
	maintain steady state
	"""

	nPolymerizingEnzymeNeeded = units.sum(
		productLengths / elongationRates
		* netLossRate)
	return nPolymerizingEnzymeNeeded

def netLossRateFromDilutionAndDegradationProtein(doublingTime, degradationRates):
	"""
	Compute total loss rate (summed contributions of degradation and dilution).

	Inputs
	------
	- doublingTime (float with units of time) - doubling time of the cell
	- degradationRates (array of floats with units of 1/time) - protein degradation rate

	Returns
	--------
	- array of floats with units of 1/time for the total loss rate for each protein
	"""

	return np.log(2) / doublingTime + degradationRates

def netLossRateFromDilutionAndDegradationRNA(doublingTime, totalEndoRnaseCountsCapacity, Km, rnaConc, countsToMolar):
	"""
	Compute total loss rate (summed impact of degradation and dilution).
	Returns the loss rate in units of (counts/time) in preparation for use in
	the steady state analysis in fitExpression() and
	setRNAPCountsConstrainedByPhysiology()
	(see calculateMinPolymerizingEnzymeByProductDistributionRNA()).

	Derived from steady state analysis of Michaelis-Menten enzyme kinetics with
	competitive inhibition: for a given RNA, all other RNAs compete for RNase.

	V_i = k_cat * [ES_i]
	v_i = k_cat * [E]0 * ([S_i]/Km_i) / (1 + sum over j genes([S_j] / Km_j))

	Inputs
	------
	- doublingTime (float with units of time) - doubling time of the cell
	- totalEndoRnaseCountsCapacity (float with units of 1/time) total kinetic
	capacity of all RNases in the cell
	- Km (array of floats with units of mol/volume) - Michaelis-Menten constant
	for each RNA
	- rnaConc (array of floats with units of mol/volume) - concentration for each RNA
	- countsToMolar (float with units of mol/volume) - conversion between counts and molar

	Returns
	--------
	- array of floats with units of 1/time for the total loss rate for each RNA
	"""

	fracSaturated = rnaConc / Km / (1 + units.sum(rnaConc / Km))
	rnaCounts = (1 / countsToMolar) * rnaConc

	return (np.log(2) / doublingTime) * rnaCounts + (totalEndoRnaseCountsCapacity * fracSaturated)

def netLossRateFromDilutionAndDegradationRNALinear(doublingTime, degradationRates, rnaCounts):
	"""
	Compute total loss rate (summed contributions of degradation and dilution).
	Returns the loss rate in units of (counts/time) in preparation for use in
	the steady state analysis in fitExpression() and
	setRNAPCountsConstrainedByPhysiology()
	(see calculateMinPolymerizingEnzymeByProductDistributionRNA()).

	Requires
	--------
	- doublingTime (float with units of time) - doubling time of the cell
	- degradationRates (array of floats with units of 1/time) - degradation rate
	for each RNA
	- rnaCounts (array of floats) - counts for each RNA

	Returns
	--------
	- array of floats with units of 1/time for the total loss rate for each RNA
	"""

	return (np.log(2) / doublingTime + degradationRates) * rnaCounts


def expressionFromConditionAndFoldChange(transcription, condPerturbations, tfFCs):
	"""
	Adjusts expression of RNA based on fold changes from basal for a given
	condition. Since fold changes are reported for individual RNA cistrons, the
	changes are applied to the basal expression levels of each cistron and the
	resulting vector is mapped back to RNA expression through nonnegative least
	squares. For genotype perturbations, the expression of all RNAs that include
	the given cistron are set to the given value.

	Inputs
	------
	- transcription: Instance of the Transcription class from
		reconstruction.ecoli.dataclasses.process.transcription
	- condPerturbations {cistron ID (str): fold change (float)} -
		dictionary of fold changes for cistrons based on the given condition
	- tfFCs {cistron ID (str): fold change (float)} -
		dictionary of fold changes for cistrons based on transcription factors
		in the given condition

	Returns
	--------
	- expression (array of floats) - adjusted expression for each RNA,
	normalized to 1

	Notes
	-----
	- TODO (Travis) - Might not properly handle if an RNA is adjusted from both a
	perturbation and a transcription factor, currently RNA self regulation is not
	included in tfFCs
	"""
	cistron_ids = transcription.cistron_data['id']
	cistron_expression = transcription.fit_cistron_expression['basal'].copy()

	# Gather indices and fold changes for each cistron that will be adjusted
	cistron_id_to_index = {
		cistron_id: i for (i, cistron_id) in enumerate(cistron_ids)
		}
	cistron_indexes = []
	cistron_fcs = []

	# Compile indexes and fold changes of each cistron
	for cistron_id, fc_value in tfFCs.items():
		if cistron_id in condPerturbations:
			continue
		cistron_indexes.append(cistron_id_to_index[cistron_id])
		cistron_fcs.append(fc_value)

	def apply_fcs_to_expression(expression, indexes, fcs):
		"""
		Applys the fold-change values to an expression vector while keeping the
		sum of expression values at one.

		Args:
			expression (np.ndarray of floats): Original expression vector of
				cistrons or RNAs
			indexes (List of floats): Indexes of cistrons/RNAs that the
				fold-changes should be applied to
			fcs (List of floats): Fold-changes of cistron/RNA expression
		"""
		fcs = [fc for (idx, fc) in
			sorted(zip(indexes, fcs), key=lambda pair: pair[0])]
		indexes = [idx for (idx, fc) in
			sorted(zip(indexes, fcs), key=lambda pair: pair[0])]

		# Adjust expression based on fold change and normalize
		indexes_bool = np.zeros(len(expression), dtype=bool)
		indexes_bool[indexes] = 1
		fcs = np.array(fcs)
		scaleTheRestBy = (1. - (expression[indexes] * fcs).sum()) / (
			1. - (expression[indexes]).sum())
		expression[indexes_bool] *= fcs
		expression[~indexes_bool] *= scaleTheRestBy

		return expression

	cistron_expression = apply_fcs_to_expression(
		cistron_expression, cistron_indexes, cistron_fcs)

	# Use NNLS to map new cistron expression to RNA expression
	expression, _ = transcription.fit_rna_expression(cistron_expression)
	expression = normalize(expression)

	# Apply genotype perturbations to all RNAs that contain each cistron
	rna_indexes = []
	rna_fcs = []
	cistron_perturbation_indexes = []
	cistron_perturbation_values = []

	for cistron_id, perturbation_value in condPerturbations.items():
		rna_indexes_with_cistron = transcription.cistron_id_to_rna_indexes(cistron_id)
		rna_indexes.extend(rna_indexes_with_cistron)
		rna_fcs.extend([perturbation_value] * len(rna_indexes_with_cistron))
		cistron_perturbation_indexes.append(cistron_id_to_index[cistron_id])
		cistron_perturbation_values.append(perturbation_value)

	expression = apply_fcs_to_expression(expression, rna_indexes, rna_fcs)
	# Also apply perturbations to cistrons for bookkeeping purposes
	cistron_expression = apply_fcs_to_expression(
		cistron_expression, cistron_perturbation_indexes, cistron_perturbation_values)

	return expression, cistron_expression

def fitPromoterBoundProbability(sim_data, cell_specs):
	r"""
	Calculates the probabilities (P) that each transcription factor will bind
	to its target RNA. This function initially calculates these probabilities
	from the bulk average counts of the TFs and ligands calculated from
	previous steps. Then, values of parameters alpha and r in the equation
	below are fit such that the computed RNA synthesis probabilities converge
	to the measured RNA synthesis probabilities.

	v_{synth, j} = \alpha_j + \sum_{i} P_{T,i}*r_{ij}

	Due to constraints applied in the optimization, both v and P need to
	be shifted from their initial values.

	Requires
	--------
	- Bulk average counts of transcription factors and associated ligands
	for each condition (in cell_specs)

	Inputs
	------
	- cell_specs {condition (str): dict} - information about each condition

	Modifies
	--------
	- Probabilities of TFs binding to their promoters
	- RNA synthesis probabilities
	- cell_specs['basal']['r_vector']: Fit parameters on how the recruitment of
	a TF affects the expression of a gene. High (positive) values of r indicate
	that the TF binding increases the probability that the gene is expressed.
	- cell_specs['basal']['r_columns']: mapping of column name to index in r

	Notes
	--------
	See supplementary materials on transcription regulation for details on
	the parameters being fit.
	"""

	def build_vector_k(sim_data, cell_specs):
		"""
		Construct vector k that contains existing fit transcription
		probabilities of RNAs in each relevant condition, normalized by the
		average copy number of the gene encoding the RNA while the cell grows
		in that condition.

		Returns
		--------
		- k: List of RNA synthesis probabilities for each RNA and condition,
		normalized by gene copy number.
		- kInfo: List of dictionaries that hold information on values of k -
		kInfo[i]["condition"] and kInfo[i]["idx"] hold what condition and RNA
		index the probability k[i] refers to, respectively.
		"""

		k, kInfo = [], []

		for idx, (rnaId, rnaCoordinate) in enumerate(
				zip(sim_data.process.transcription.rna_data["id"],
				sim_data.process.transcription.rna_data['replication_coordinate'])):

			# Get list of TFs that regulate this RNA
			tfs = sim_data.relation.rna_id_to_regulating_tfs.get(rnaId, [])
			conditions = ["basal"]
			tfsWithData = []

			# Take only those TFs with active/inactive conditions data
			# TODO (Gwanggyu): cache this list of conditions for each RNA
			for tf in tfs:
				if tf not in sorted(sim_data.tf_to_active_inactive_conditions):
					continue

				# Add conditions for selected TFs
				conditions.append(tf + "__active")
				conditions.append(tf + "__inactive")
				tfsWithData.append(tf)

			for condition in conditions:
				# Skip basal conditions, unless the RNA is not regulated by any TFs
				if len(tfsWithData) > 0 and condition == "basal":
					continue

				# Get specific doubling time for this condition
				tau = cell_specs[condition]["doubling_time"].asNumber(units.min)

				# Calculate average copy number of gene for this condition
				n_avg_copy = sim_data.process.replication.get_average_copy_number(tau, rnaCoordinate)

				# Compute synthesis probability per gene copy
				prob_per_copy = sim_data.process.transcription.rna_synth_prob[condition][idx] / n_avg_copy

				# Gather RNA synthesis probabilities for each RNA per condition
				k.append(prob_per_copy)
				kInfo.append({"condition": condition, "idx": idx})

		k = np.array(k)

		return k, kInfo

	def build_matrix_G(sim_data, pPromoterBound):
		"""
		Construct matrix G that contains probabilities of pPromoterBound as
		elements. Each row of the matrix is named "[RNA]__[condition]", where
		there are two conditions [active/inactive] for each TF that regulates
		the expression of the given RNA. For RNAs that are not regulated by
		any TFs, a single row named "[RNA]__basal" represents the RNA. Each
		column is named "[RNA]__[TF]", for each TF that regulates the
		expression of the given RNA. Each element is set to the value in
		pPromoterBound that corresponds to the condition given by the row,
		and the TF given by the column. For each RNA, there is an additional
		column named "[RNA]__alpha", and all elements in this column that
		corresponds to the rows for the RNA are set to 1.

		Inputs
		------
		- pPromoterBound: Probabilities that a given TF is bound to its
		promoter in a given condition, calculated from bulk average
		concentrations of the TF and its associated ligands.

		Returns
		--------
		- G: Matrix of values in pPromoterBound, rearranged based on each RNA
		- row_name_to_index: Dict[str, int] of row names of G to row index
		- col_name_to_index: Dict[str, int] of column names of G to column index
		"""

		gI, gJ, gV = [], [], []
		row_name_to_index, col_name_to_index = {}, {}

		for idx, rnaId in enumerate(sim_data.process.transcription.rna_data["id"]):
			rnaIdNoLoc = rnaId[:-3]  # Remove compartment ID from RNA ID

			# Get list of TFs that regulate this RNA
			tfs = sim_data.relation.rna_id_to_regulating_tfs.get(rnaId, [])
			conditions = ["basal"]
			tfsWithData = []

			# Take only those TFs with active/inactive conditions data
			for tf in tfs:
				if tf not in sorted(sim_data.tf_to_active_inactive_conditions):
					continue

				# Add conditions for selected TFs
				conditions.append(tf + "__active")
				conditions.append(tf + "__inactive")
				tfsWithData.append(tf)

			for condition in conditions:
				# Skip basal conditions, unless the RNA is not regulated by any TFs
				if len(tfsWithData) > 0 and condition == "basal":
					continue

				# Add row for each condition specific to each RNA
				row_name = rnaIdNoLoc + "__" + condition
				row_name_to_index[row_name] = len(row_name_to_index)

				for tf in tfsWithData:
					# Add column for each TF that regulates each RNA
					col_name = rnaIdNoLoc + "__" + tf

					# TODO (Gwanggyu): Are these checks necessary?
					if col_name not in col_name_to_index:
						col_name_to_index[col_name] = len(col_name_to_index)

					gI.append(row_name_to_index[row_name])
					gJ.append(col_name_to_index[col_name])
					gV.append(pPromoterBound[condition][tf])  # Probability that TF is bound in given condition

				# Add alpha column for each RNA
				col_name = rnaIdNoLoc + "__alpha"

				if col_name not in col_name_to_index:
					col_name_to_index[col_name] = len(col_name_to_index)

				gI.append(row_name_to_index[row_name])
				gJ.append(col_name_to_index[col_name])
				gV.append(1.)

		gI, gJ, gV = np.array(gI), np.array(gJ), np.array(gV)
		G = np.zeros((len(row_name_to_index), len(col_name_to_index)), np.float64)
		G[gI, gJ] = gV

		return G, row_name_to_index, col_name_to_index

	def build_matrix_Z(sim_data, col_name_to_index):
		"""
		Construct matrix Z that connects all possible TF combinations with
		each TF. Each row of the matrix corresponds to an RNA-(TF combination)
		pair, and each column corresponds to an RNA-TF pair, with an additional
		RNA-alpha column for each RNA (identical to matrix G). Matrix values
		are set to one if the TF specified by the column is "active" in the
		combination specified by the row or if the column is an RNA-alpha
		column, and zero otherwise.

		Inputs
		------
		- col_name_to_index: Dict[str, int] of column names of G to column index

		Returns
		--------
		- Z: Matrix of zeros and ones, specifying which TFs in the columns
		correspond to combinations in the rows.
		"""

		zI, zJ, zV = [], [], []
		row_idx = 0

		for rna_id in sim_data.process.transcription.rna_data["id"]:
			rna_id_no_loc = rna_id[:-3]  # Remove compartment ID from RNA ID

			# Get list of TFs that regulate this RNA
			tfs = sim_data.relation.rna_id_to_regulating_tfs[rna_id]
			tfs_with_data = []

			# Get column index of the RNA's alpha column
			col_idxs = [col_name_to_index[rna_id_no_loc + "__alpha"]]

			# Take only those TFs with active/inactive conditions data
			for tf in tfs:
				if tf not in sim_data.tf_to_active_inactive_conditions:
					continue

				tfs_with_data.append(tf)

				# Get column index of the RNA-TF pair
				col_idxs.append(col_name_to_index[rna_id_no_loc + "__" + tf])

			n_tfs = len(tfs_with_data)

			# For all possible combinations of TFs
			for n_combinations in range(n_tfs + 1):
				for combination in itertools.combinations(range(1, n_tfs + 1), n_combinations):
					# Always include alpha column
					zI.append(row_idx)
					zJ.append(col_idxs[0])
					zV.append(1)

					# Set matrix value to one if the TF specified by the column is
					# present in the combination of TFs specified by the row
					for col_idx in combination:
						zI.append(row_idx)
						zJ.append(col_idxs[col_idx])
						zV.append(1)

					row_idx += 1

		# Build matrix Z
		zI, zJ, zV = np.array(zI), np.array(zJ), np.array(zV)
		Z = np.zeros((zI.max() + 1, zJ.max() + 1), np.float64)
		Z[zI, zJ] = zV

		return Z

	def build_matrix_T(sim_data, col_name_to_index):
		"""
		Construct matrix T that specifies the direction of regulation for each
		RNA-TF pair.

		Inputs
		------
		- col_name_to_index: Dict[str, int] of column names of G to column index

		Returns
		--------
		- T: Diagonal matrix. Diagonal value is +1 if the direction of
		regulation by the TF-RNA pair specified by the row is positive, -1 if
		this is negative, and 0 if the row is an RNA_alpha row.
		"""

		tI, tJ, tV = [], [], []
		row_idx = 0

		for rnaId in sim_data.process.transcription.rna_data["id"]:
			rnaIdNoLoc = rnaId[:-3]  # Remove compartment ID from RNA ID

			# Get list of TFs that regulate this RNA
			tfs = sim_data.relation.rna_id_to_regulating_tfs[rnaId]
			tfsWithData = []

			# Get list of constituent cistron IDs
			constituent_cistron_ids = [
				sim_data.process.transcription.cistron_data['id'][i]
				for i in sim_data.process.transcription.rna_id_to_cistron_indexes(rnaId)
				]

			# Take only those TFs with active/inactive conditions data
			for tf in tfs:
				if tf not in sim_data.tf_to_active_inactive_conditions:
					continue

				tfsWithData.append(tf)

			for tf in tfsWithData:
				# Calculate the consensus regulation direction of the TF on the
				# cistrons that constitute this transcription unit
				directions = np.array([
					sim_data.tf_to_direction[tf].get(cistron_id, 0)
					for cistron_id in constituent_cistron_ids
					])
				consensus_direction = -1 + 2*(directions.sum() >= 0)

				# Add row for TF and find column for TF in col_name_to_index
				col_name = rnaIdNoLoc + "__" + tf

				# Set matrix value to regulation direction (+1 or -1)
				tI.append(row_idx)
				tJ.append(col_name_to_index[col_name])
				tV.append(consensus_direction)
				row_idx += 1

			# Add RNA_alpha rows and columns, and set matrix value to zero
			col_name = rnaIdNoLoc + "__alpha"

			tI.append(row_idx)
			tJ.append(col_name_to_index[col_name])
			tV.append(0)
			row_idx += 1

		tI, tJ, tV = np.array(tI), np.array(tJ), np.array(tV)
		T = np.zeros((tI.max() + 1, tJ.max() + 1), np.float64)
		T[tI, tJ] = tV

		return T

	def build_matrix_H(sim_data, col_name_to_index, pPromoterBound, r, fixedTFs, cell_specs):
		r"""
		Construct matrix H that contains values of vector r as elements.
		Each row of the matrix is named "[RNA]__[condition]", where
		there are two conditions [active/inactive] for each TF that regulates
		the expression of the given RNA. For RNAs that are not regulated by
		any TFs, a single row named "[RNA]__basal" represents the RNA. Each
		column is named "[TF]__[condition]", for each TF that regulates the
		expression of the given RNA, and each condition given for the RNA.
		Each element is set to the optimized value in r that corresponds to
		the RNA given by the row, and the TF given by the column. For each RNA,
		there is an additional column named "[RNA]__alpha", and all elements in
		this column are set to the value of parameter alpha for the RNA
		optimized in r.

		Inputs
		------
		- col_name_to_index: Dict[str, int] of column names of G to column index
		- pPromoterBound: Probabilities that a given TF is bound to its
		promoter in a given condition, calculated from bulk average
		concentrations of the TF and its associated ligands.
		- r: Optimized values of \Delta r (the effect of TF on RNAP
		recruitment) and alpha (basal recruitment of RNAP)
		- fixedTFs: List of TFs whose activities do not change with the
		nutrient conditions

		Returns
		--------
		- H: Matrix of values in optimized r, rearranged for each RNA and
		condition
		- pInit: Vector of values in pPromoterBound, rearranged to the ordering
		of the columns of H
		- pAlphaIdxs: Indexes of columns that correspond to alpha's in H and pInit
		- pNotAlphaIdxs: Indexes of columns that correspond to r's in H and pInit
		- fixedTFIdxs: Indexes of columns that correspond to fixed TFs in H and pInit
		- pPromoterBoundIdxs: Dictionary of indexes to pInit.
		- H_col_name_to_index: Dict[str, int] of column names of H to column index
		"""

		rDict = dict([(col_name, value) for col_name, value in zip(col_name_to_index, r)])

		pPromoterBoundIdxs = dict([(condition, {}) for condition in pPromoterBound])
		hI, hJ, hV, pInitI, pInitV = [], [], [], [], []
		H_row_name_to_index, H_col_name_to_index = {}, {}

		for idx, rnaId in enumerate(sim_data.process.transcription.rna_data["id"]):
			rnaIdNoLoc = rnaId[:-3]  # Remove compartment ID from RNA ID

			tfs = sim_data.relation.rna_id_to_regulating_tfs[rnaId]
			conditions = ["basal"]
			tfsWithData = []

			# Take only those TFs with active/inactive conditions data
			for tf in tfs:
				if tf not in sorted(sim_data.tf_to_active_inactive_conditions):
					continue

				# Add conditions for selected TFs
				conditions.append(tf + "__active")
				conditions.append(tf + "__inactive")
				tfsWithData.append(tf)

			for condition in conditions:
				# Skip basal conditions, unless the RNA is not regulated by any TFs
				if len(tfsWithData) > 0 and condition == "basal":
					continue

				# Add row for each condition specific to each RNA
				row_name = rnaIdNoLoc + "__" + condition
				H_row_name_to_index[row_name] = len(H_row_name_to_index)

				for tf in tfsWithData:
					# Add column for each TF and condition
					col_name = tf + "__" + condition

					if col_name not in H_col_name_to_index:
						H_col_name_to_index[col_name] = len(H_col_name_to_index)

					hI.append(H_row_name_to_index[row_name])
					hJ.append(H_col_name_to_index[col_name])

					# Handle the case of the TF being knocked out (admittedly not the cleanest solution)
					if cell_specs[condition]["bulkAverageContainer"].count(tf + "[c]") == 0:
						hV.append(0)  # TF is knocked out in the given condition
					else:
						hV.append(rDict[rnaIdNoLoc + "__" + tf])  # Optimized r value for TF-RNA pair

					# Rearrange values in pPromoterBound in the same order
					# given by the columns of H
					pInitI.append(H_col_name_to_index[col_name])
					pInitV.append(pPromoterBound[condition][tf])
					pPromoterBoundIdxs[condition][tf] = H_col_name_to_index[col_name]

				# Add alpha column for each RNA
				col_name = rnaIdNoLoc + "__alpha"

				if col_name not in H_col_name_to_index:
					H_col_name_to_index[col_name] = len(H_col_name_to_index)

				# Add optimized value of alpha in r to H
				hI.append(H_row_name_to_index[row_name])
				hJ.append(H_col_name_to_index[col_name])
				hV.append(rDict[col_name])

				# Set corresponding value in pInit to one
				pInitI.append(H_col_name_to_index[col_name])
				pInitV.append(1.)

		# Save indices to update promoter binding for active TFs in combined conditions
		for condition, tfs in sim_data.condition_active_tfs.items():
			for tf in tfs:
				col_name = f'{tf}__{tf}__active'
				pPromoterBoundIdxs[condition][tf] = H_col_name_to_index[col_name]

		# Save indices to update promoter binding for inactive TFs in combined conditions
		for condition, tfs in sim_data.condition_inactive_tfs.items():
			for tf in tfs:
				col_name = f'{tf}__{tf}__inactive'
				pPromoterBoundIdxs[condition][tf] = H_col_name_to_index[col_name]

		# Build vector pInit and matrix H
		pInit = np.zeros(len(set(pInitI)))
		pInit[pInitI] = pInitV

		hI, hJ, hV = np.array(hI), np.array(hJ), np.array(hV)
		Hshape = (hI.max() + 1, hJ.max() + 1)
		H = np.zeros(Hshape, np.float64)
		H[hI, hJ] = hV

		# Get indexes of alpha and non-alpha columns in pInit and H
		pAlphaIdxs = np.array([idx for col_name, idx in H_col_name_to_index.items() if col_name.endswith("__alpha")])
		pNotAlphaIdxs = np.array([idx for col_name, idx in H_col_name_to_index.items() if not col_name.endswith("__alpha")])

		# Get indexes of columns that correspond to fixed TFs
		fixedTFIdxs = []
		for col_name, idx in H_col_name_to_index.items():
			secondElem = col_name.split("__")[1]

			if secondElem in fixedTFs:
				fixedTFIdxs.append(idx)

		fixedTFIdxs = np.array(fixedTFIdxs, dtype=int)

		return H, pInit, pAlphaIdxs, pNotAlphaIdxs, fixedTFIdxs, pPromoterBoundIdxs, H_col_name_to_index

	def build_matrix_pdiff(sim_data, H_col_name_to_index):
		"""
		Construct matrix Pdiff that specifies the indexes of corresponding
		TFs and conditions.

		Inputs
		------
		- H_col_name_to_index: Dict[str, int] of column names of H to column index

		Returns
		--------
		- Pdiff: Matrix with [TF] as rows and [TF]_[condition] as columns.
		Matrix value is set to 1 when the TF of the column matches with the TF
		of the row, and the condition is TF__active. Matrix value is set to -1
		when the TF of the column matches with the TF of the row, and the
		condition is TF__inactive.
		"""

		PdiffI, PdiffJ, PdiffV = [], [], []

		for rowIdx, tf in enumerate(sorted(sim_data.tf_to_active_inactive_conditions)):
			# For each TF, find condition [TF]__[TF]__active and set element to 1
			condition = tf + "__active"
			col_name = tf + "__" + condition
			PdiffI.append(rowIdx)
			PdiffJ.append(H_col_name_to_index[col_name])
			PdiffV.append(1)

			# Find condition [TF]__[TF]__inactive and set element to -1
			condition = tf + "__inactive"
			col_name = tf + "__" + condition
			PdiffI.append(rowIdx)
			PdiffJ.append(H_col_name_to_index[col_name])
			PdiffV.append(-1)

		# Build matrix Pdiff
		PdiffI, PdiffJ, PdiffV = np.array(PdiffI), np.array(PdiffJ), np.array(PdiffV)
		Pdiffshape = (PdiffI.max() + 1, len(H_col_name_to_index))
		Pdiff = np.zeros(Pdiffshape, np.float64)
		Pdiff[PdiffI, PdiffJ] = PdiffV

		return Pdiff

	def fromArray(p, pPromoterBound, pPromoterBoundIdxs):
		"""
		Updates values in pPromoterBound with fit probabilities.

		Inputs
		------
		- p: Vector of probabilities optimized in the current step.
		- pPromoterBoundIdxs: Dictionary of indexes to p

		Modifies
		--------
		Values in pPromoterBound - probabilities that each transcription factor
		is bound to its promoter for each growth condition
		"""

		for condition in sorted(pPromoterBoundIdxs):
			for tf in sorted(pPromoterBoundIdxs[condition]):
				pPromoterBound[condition][tf] = p[pPromoterBoundIdxs[condition][tf]]

	def updateSynthProb(sim_data, cell_specs, kInfo, k):
		"""
		Updates RNA synthesis probabilities with fit values of P and R. The
		expected average copy number of genes for the condition are multiplied
		to the per-copy probabilities.

		Inputs
		------
		- kInfo: List of dictionaries that hold information on values of k -
		kInfo[i]["condition"] and kInfo[i]["idx"] hold what condition and RNA
		index the probability k[i] refers to, respectively.
		- k: RNA synthesis probabilities computed from fit P and R.

		Modifies
		--------
		- RNA synthesis probabilities

		Notes
		--------
		These values are used to calculate pre-set probabilities for
		pools of mRNA, tRNA, and rRNA, and the genes that encode RNAP and
		ribosome subunits.
		"""

		# Get replication coordinates of each RNA
		replication_coordinate = sim_data.process.transcription.rna_data['replication_coordinate']

		# Update sim_data values with fit values
		for D, k_value in zip(kInfo, k):
			condition = D["condition"]
			rna_idx = D["idx"]

			# Get coordinate of RNA
			rnaCoordinate = replication_coordinate[rna_idx]

			# Get specific doubling time for this condition
			tau = cell_specs[condition]["doubling_time"].asNumber(units.min)

			# Calculate average copy number of gene for this condition
			n_avg_copy = sim_data.process.replication.get_average_copy_number(tau, rnaCoordinate)

			# Multiply copy number by k_value to get synthesis probability
			# (if k_value is negative clip at zero)
			sim_data.process.transcription.rna_synth_prob[condition][rna_idx] = max(0, k_value) * n_avg_copy

		# Normalize values such that probabilities for each condition sum to one
		for condition in sim_data.process.transcription.rna_synth_prob:
			assert np.all(sim_data.process.transcription.rna_synth_prob[condition] >= 0)
			sim_data.process.transcription.rna_synth_prob[condition] /= sim_data.process.transcription.rna_synth_prob[condition].sum()

	# Initialize pPromoterBound using mean TF and ligand concentrations
	pPromoterBound = calculatePromoterBoundProbability(sim_data, cell_specs)
	pInit0 = None
	lastNorm = np.inf

	fixedTFs = []
	for tf in sim_data.tf_to_active_inactive_conditions:
		if sim_data.process.transcription_regulation.tf_to_tf_type[tf] == "2CS":
			fixedTFs.append(tf)
		if (sim_data.process.transcription_regulation.tf_to_tf_type[tf] == "1CS"
				and sim_data.tf_to_active_inactive_conditions[tf]["active nutrients"] ==
				sim_data.tf_to_active_inactive_conditions[tf]["inactive nutrients"]):
			fixedTFs.append(tf)

	# Build vector of existing fit transcription probabilities
	k, kInfo = build_vector_k(sim_data, cell_specs)

	# Repeat for a fixed maximum number of iterations
	for i in range(PROMOTER_MAX_ITERATIONS):
		# Build matrices used in optimizing R
		G, G_row_name_to_index, G_col_name_to_index = build_matrix_G(sim_data, pPromoterBound)
		Z = build_matrix_Z(sim_data, G_col_name_to_index)
		T = build_matrix_T(sim_data, G_col_name_to_index)

		# Optimize R such that transcription initiation probabilities computed
		# from existing values of P in matrix G are close to fit values.
		R = Variable(G.shape[1])  # Vector of r's and alpha's

		# Objective: minimize difference between k and G @ R
		objective_r = Minimize(
			norm(G @ (PROMOTER_SCALING * R) - PROMOTER_SCALING * k, PROMOTER_NORM_TYPE))

		# Optimization constraints
		# 1) 0 <= Z @ R <= 1 : Assuming P = 1 for all TFs, all possible
		# combinations of TFs should yield a valid transcription probability
		# value between zero and one.
		# 2) T @ R >= 0 : Values of r for positive regulation should be positive,
		# and values of r for negative regulation should be negative.
		constraint_r = [
			0 <= Z @ (PROMOTER_SCALING * R), Z @ (PROMOTER_SCALING * R) <= PROMOTER_SCALING,
			T @ (PROMOTER_SCALING * R) >= 0]

		# Solve optimization problem
		prob_r = Problem(objective_r, constraint_r)
		prob_r.solve(solver='ECOS', max_iters=1000)

		if prob_r.status == 'optimal_inaccurate':
			raise RuntimeError('Solver found an optimum that is inaccurate.'
				' Try increasing max_iters or adjusting tolerances.')
		elif prob_r.status != 'optimal':
			raise RuntimeError('Solver could not find optimal value')

		# Get optimal value of R
		r = np.array(R.value).reshape(-1)
		r[np.abs(r) < ECOS_0_TOLERANCE] = 0  # Adjust to 0 for small values from solver tolerance

		# Use optimal value of R to construct matrix H and vector Pdiff
		H, pInit, pAlphaIdxs, pNotAlphaIdxs, fixedTFIdxs, pPromoterBoundIdxs, H_col_name_to_index = build_matrix_H(
			sim_data, G_col_name_to_index, pPromoterBound, r, fixedTFs, cell_specs)
		pdiff = build_matrix_pdiff(sim_data, H_col_name_to_index)

		# On first iteration, save the value of the initial p
		if i == 0:
			pInit0 = pInit.copy()

		# Optimize P such that the transcription probabilities computed from
		# current values of R in matrix H are close to fit values.
		P = Variable(H.shape[1])

		# Construct a boolean vector that marks column indexes of H
		# corresponding to alpha's and fixed TFs
		D = np.zeros(H.shape[1])
		D[pAlphaIdxs] = 1
		D[fixedTFIdxs] = 1

		# Mask initial p with boolean vector constructed above
		Drhs = pInit0.copy()
		Drhs[D != 1] = 0

		# Objective: minimize difference between k (fit RNAP initiation
		# probabilities) and H @ P (computed initiation probabilities) while
		# also minimizing deviation of P from the original value calculated
		# from mean TF and ligand concentrations
		objective_p = Minimize(
			norm(H @ (PROMOTER_SCALING * P) - PROMOTER_SCALING * k, PROMOTER_NORM_TYPE)
			+ PROMOTER_REG_COEFF * norm(P - pInit0, PROMOTER_NORM_TYPE))

		# Constraints
		# 1) 0 <= P <= 1 : All DNA-bound probabilities should be between zero
		# and one.
		# 2) D @ P == Drhs : Values of P that correspond to alpha's and fixed TFs
		# should not change.
		# 3) pdiff @ P >= PROMOTER_PDIFF_THRESHOLD : There must be at least a
		# certain difference between binding probabilities of a TF in conditions
		# TF__active and TF__inactive
		constraint_p = [
			0 <= PROMOTER_SCALING * P, PROMOTER_SCALING * P <= PROMOTER_SCALING,
			np.diag(D) @ (PROMOTER_SCALING * P) == PROMOTER_SCALING * Drhs,
			pdiff @ (PROMOTER_SCALING * P) >= PROMOTER_SCALING * PROMOTER_PDIFF_THRESHOLD,
			]

		# Solve optimization problem
		prob_p = Problem(objective_p, constraint_p)
		prob_p.solve(solver='ECOS')

		if prob_p.status == 'optimal_inaccurate':
			raise RuntimeError('Solver found an optimum that is inaccurate.'
				' Try increasing max_iters or adjusting tolerances.')
		elif prob_p.status != 'optimal':
			raise RuntimeError('Solver could not find optimal value')

		# Get optimal value of P
		p = np.array(P.value).reshape(-1)

		# Adjust for solver tolerance over bounds to get proper probabilities
		p[p < ECOS_0_TOLERANCE] = 0
		p[p > (1 - ECOS_0_TOLERANCE)] = 1

		# Update pPromoterBound with fit p
		fromArray(p, pPromoterBound, pPromoterBoundIdxs)

		# Break from loop if parameters have converged
		if np.abs(np.linalg.norm(np.dot(H, p) - k, PROMOTER_NORM_TYPE) - lastNorm) < PROMOTER_CONVERGENCE_THRESHOLD:
			break
		else:
			lastNorm = np.linalg.norm(np.dot(H, p) - k, PROMOTER_NORM_TYPE)

	# Update sim_data with fit bound probabilities and RNAP initiation
	# probabilities computed from these bound probabilities
	sim_data.pPromoterBound = pPromoterBound
	updateSynthProb(sim_data, cell_specs, kInfo, np.dot(H, p))

	cell_specs['basal']['r_vector'] = r
	cell_specs['basal']['r_columns'] = G_col_name_to_index

def fitLigandConcentrations(sim_data, cell_specs):
	"""
	Using the fit values of pPromoterBound, updates the set concentrations of
	ligand metabolites and the kd's of the ligand-TF binding reactions.

	Requires
	--------
	- Fitted pPromoterBound: probabilities that a TF will bind to its promoter,
	fit by function fitPromoterBoundProbability().

	Inputs
	------
	- cell_specs {condition (str): dict} - information about each condition

	Modifies
	--------
	- Set concentrations of metabolites that are ligands in 1CS
	- kd's of equilibrium reactions in 1CS
	"""
	cellDensity = sim_data.constants.cell_density
	pPromoterBound = sim_data.pPromoterBound

	for tf in sorted(sim_data.tf_to_active_inactive_conditions):
		# Skip TFs that are not 1CS or are linked to genotypic perturbations
		if sim_data.process.transcription_regulation.tf_to_tf_type[tf] != "1CS":
			continue
		if (len(sim_data.tf_to_active_inactive_conditions[tf]["active genotype perturbations"]) > 0
				or len(sim_data.tf_to_active_inactive_conditions[tf]["inactive genotype perturbations"]) > 0):
			continue

		activeKey = tf + "__active"
		inactiveKey = tf + "__inactive"

		# Determine if metabolite-bound form of the TF is the active form
		boundId = sim_data.process.transcription_regulation.active_to_bound[tf]
		negativeSignal = (tf != boundId)  # True if unbound form is the active TF

		# Calculate kd of bound TF
		fwdRate = sim_data.process.equilibrium.get_fwd_rate(boundId + "[c]")
		revRate = sim_data.process.equilibrium.get_rev_rate(boundId + "[c]")
		kd = revRate/fwdRate

		# Get the metabolite that binds to the TF and its stoich coefficient
		metabolite = sim_data.process.equilibrium.get_metabolite(boundId + "[c]")
		metaboliteCoeff = sim_data.process.equilibrium.get_metabolite_coeff(boundId + "[c]")

		# Calculate the concentrations of the metabolite under conditions where
		# TF is active and inactive
		activeCellVolume = (cell_specs[activeKey]["avgCellDryMassInit"] /
			cellDensity / sim_data.mass.cell_dry_mass_fraction)
		activeCountsToMolar = 1/(sim_data.constants.n_avogadro * activeCellVolume)
		activeSignalConc = ((activeCountsToMolar * cell_specs[activeKey]["bulkAverageContainer"].count(metabolite))
			.asNumber(units.mol/units.L))
		inactiveCellVolume = (cell_specs[inactiveKey]["avgCellDryMassInit"] /
			cellDensity / sim_data.mass.cell_dry_mass_fraction)
		inactiveCountsToMolar = 1/(sim_data.constants.n_avogadro * inactiveCellVolume)
		inactiveSignalConc = ((inactiveCountsToMolar * cell_specs[inactiveKey]["bulkAverageContainer"].count(metabolite))
			.asNumber(units.mol/units.L))

		# Update kd with fitted values of P and the bulk average concentrations
		# of the metabolite, and use this fitted kd to recalculate the set
		# amounts of the metabolite in metabolism
		p_active = pPromoterBound[activeKey][tf]
		p_inactive = pPromoterBound[inactiveKey][tf]

		if negativeSignal:
			if p_inactive == 0:
				raise ValueError('Inf ligand concentration from p_inactive = 0.'
					' Check results from fitPromoterBoundProbability and Kd values.')
			if 1 - p_active < 1e-9:
				kdNew = kd  # Concentration of metabolite-bound TF is negligible
			else:
				kdNew = ((activeSignalConc**metaboliteCoeff) * p_active/(1 - p_active))**(1/metaboliteCoeff)

			# Reset metabolite concentration with fitted P and kd
			sim_data.process.metabolism.concentration_updates.molecule_set_amounts[metabolite] = (
				(kdNew**metaboliteCoeff*(1 - p_inactive)/p_inactive)**(1./metaboliteCoeff)*(units.mol/units.L))

		else:
			if p_active == 1:
				raise ValueError('Inf ligand concentration from p_active = 1.'
					' Check results from fitPromoterBoundProbability and Kd values.')
			if p_inactive < 1e-9:
				kdNew = kd  # Concentration of metabolite-bound TF is negligible
			else:
				kdNew = ((inactiveSignalConc**metaboliteCoeff) * (1 - p_inactive)/p_inactive)**(1/metaboliteCoeff)

			# Reset metabolite concentration with fitted P and kd
			sim_data.process.metabolism.concentration_updates.molecule_set_amounts[metabolite] = (
				(kdNew**metaboliteCoeff*p_active/(1 - p_active))**(1./metaboliteCoeff)*(units.mol/units.L))

		# Fit reverse rate in line with fitted kd
		sim_data.process.equilibrium.set_rev_rate(boundId + "[c]", kdNew * fwdRate)


def calculatePromoterBoundProbability(sim_data, cell_specs):
	"""
	Calculate the probability that a transcription factor is bound to its
	associated promoter for all simulated growth conditions. The bulk
	average concentrations calculated for TFs and their ligands are used to
	compute the probabilities based on the type (0CS, 1CS, 2CS) of the TF.

	Requires
	--------
	- Bulk average counts of transcription factors and associated ligands
	for each condition (in cell_specs)

	Returns
	--------
	- pPromoterBound: Probability that a transcription factor is bound to
	its promoter, per growth condition and TF. Each probability is indexed by
	pPromoterBound[condition][TF].
	"""

	pPromoterBound = {}  # Initialize return value
	cellDensity = sim_data.constants.cell_density
	init_to_average = sim_data.mass.avg_cell_to_initial_cell_conversion_factor

	# Matrix to determine number of promoters each TF can bind to in a given condition
	tf_idx = {tf: i for i, tf in enumerate(sim_data.tf_to_active_inactive_conditions)}
	cistron_id_to_tu_indexes = {
		cistron_id: sim_data.process.transcription.cistron_id_to_rna_indexes(cistron_id)
		for cistron_id in sim_data.process.transcription.cistron_data['id']}
	regulation_i = []
	regulation_j = []
	regulation_v = []
	for tf, cistrons in sim_data.tf_to_fold_change.items():
		if tf not in tf_idx:
			continue

		for cistron in cistrons:
			for tu_index in cistron_id_to_tu_indexes[cistron]:
				regulation_i.append(tf_idx[tf])
				regulation_j.append(tu_index)
				regulation_v.append(1)

	regulation = scipy.sparse.csr_matrix(
		(regulation_v, (regulation_i, regulation_j)),
		shape=(len(tf_idx), len(sim_data.process.transcription.rna_data)))
	rna_coords = sim_data.process.transcription.rna_data['replication_coordinate']

	for conditionKey in sorted(cell_specs):
		pPromoterBound[conditionKey] = {}
		tau = sim_data.condition_to_doubling_time[conditionKey].asNumber(units.min)
		n_avg_copy = sim_data.process.replication.get_average_copy_number(tau, rna_coords)
		n_promoter_targets = regulation.dot(n_avg_copy)

		cellVolume = cell_specs[conditionKey]["avgCellDryMassInit"]/cellDensity/sim_data.mass.cell_dry_mass_fraction
		countsToMolar = 1/(sim_data.constants.n_avogadro * cellVolume)

		for tf in sorted(sim_data.tf_to_active_inactive_conditions):
			tfType = sim_data.process.transcription_regulation.tf_to_tf_type[tf]
			tf_counts = cell_specs[conditionKey]["bulkAverageContainer"].count(tf + "[c]")
			tf_targets = n_promoter_targets[tf_idx[tf]]
			limited_tf_counts = min(1, tf_counts * init_to_average / tf_targets)
			if tfType == "0CS":
				pPromoterBound[conditionKey][tf] = limited_tf_counts  # If TF exists, the promoter is always bound to the TF

			elif tfType == "1CS":
				boundId = sim_data.process.transcription_regulation.active_to_bound[tf]  # ID of TF bound to ligand
				kd = sim_data.process.equilibrium.get_rev_rate(boundId + "[c]") / sim_data.process.equilibrium.get_fwd_rate(boundId + "[c]")

				signal = sim_data.process.equilibrium.get_metabolite(boundId + "[c]")  # ID of ligand that binds to TF
				signalCoeff = sim_data.process.equilibrium.get_metabolite_coeff(boundId + "[c]")  # Stoichiometric coefficient of ligand

				# Get bulk average concentrations of ligand and TF
				signalConc = (countsToMolar*cell_specs[conditionKey]["bulkAverageContainer"].count(signal)).asNumber(units.mol/units.L)
				tfConc = (countsToMolar*tf_counts).asNumber(units.mol/units.L)

				# If TF is active in its bound state
				if tf == boundId:
					if tfConc > 0:
						pPromoterBound[conditionKey][tf] = limited_tf_counts * sim_data.process.transcription_regulation.p_promoter_bound_SKd(signalConc, kd, signalCoeff)
					else:
						pPromoterBound[conditionKey][tf] = 0.

				# If TF is active in its unbound state
				else:
					if tfConc > 0:
						pPromoterBound[conditionKey][tf] = 1. - limited_tf_counts * sim_data.process.transcription_regulation.p_promoter_bound_SKd(signalConc, kd, signalCoeff)
					else:
						pPromoterBound[conditionKey][tf] = 0.

			elif tfType == "2CS":
				# Get bulk average concentrations of active and inactive TF
				activeTfConc = (countsToMolar*tf_counts).asNumber(units.mol/units.L)
				inactiveTf = sim_data.process.two_component_system.active_to_inactive_tf[tf + "[c]"]
				inactiveTfConc = (countsToMolar*cell_specs[conditionKey]["bulkAverageContainer"].count(inactiveTf)).asNumber(units.mol/units.L)

				if activeTfConc == 0 and inactiveTfConc == 0:
					pPromoterBound[conditionKey][tf] = 0.
				else:
					pPromoterBound[conditionKey][tf] = limited_tf_counts * activeTfConc/(activeTfConc + inactiveTfConc)

	# Check for any inconsistencies that could lead to feasbility issues when fitting
	for condition in pPromoterBound:
		if 'inactive' in condition:
			tf = condition.split('__')[0]
			active_p = pPromoterBound[f'{tf}__active'][tf]
			inactive_p = pPromoterBound[f'{tf}__inactive'][tf]

			if inactive_p >= active_p:
				print('Warning: active condition does not have higher binding'
					f' probability than inactive condition for {tf}'
					f' ({active_p:.3f} vs {inactive_p:.3f}).')

	return pPromoterBound


def calculateRnapRecruitment(sim_data, cell_specs):
	"""
	Constructs the basal_prob vector and delta_prob matrix from values of r.
	The basal_prob vector holds the basal transcription probabilities of each
	transcription unit. The delta_prob matrix holds the differences in
	transcription probabilities when transcription factors bind to the
	promoters of each transcription unit. Both values are stored in sim_data.

	Requires
	--------
	- cell_specs['basal']:
		- ['r_vector']: Fit parameters on how the recruitment of a TF affects the expression
		of a gene. High (positive) values of r indicate that the TF binding
		increases the probability that the gene is expressed.
		- ['r_columns']: mapping of column name to index in r

	Modifies
	--------
	- Rescales values in basal_prob such that all values are positive
	- Adds basal_prob and delta_prob arrays to sim_data
	"""

	r = cell_specs['basal']['r_vector']
	col_names_to_index = cell_specs['basal']['r_columns']

	# Get list of transcription units and TF IDs
	transcription = sim_data.process.transcription
	transcription_regulation = sim_data.process.transcription_regulation
	all_TUs = transcription.rna_data["id"]
	all_tfs = transcription_regulation.tf_ids

	# Initialize basal_prob vector and delta_prob sparse matrix
	basal_prob = np.zeros(len(all_TUs))
	deltaI, deltaJ, deltaV = [], [], []

	for rna_idx, rnaId in enumerate(all_TUs):
		rnaIdNoLoc = rnaId[:-3]  # Remove compartment ID from RNA ID

		# Take only those TFs with active/inactive conditions data
		for tf in sim_data.relation.rna_id_to_regulating_tfs.get(rnaId, []):
			if tf not in sorted(sim_data.tf_to_active_inactive_conditions):
				continue

			colName = rnaIdNoLoc + "__" + tf

			# Set element in delta to value in r that corresponds to the
			# transcription unit of the row, and the TF of the column
			deltaI.append(rna_idx)
			deltaJ.append(all_tfs.index(tf))
			deltaV.append(r[col_names_to_index[colName]])

		# Add alpha column for each RNA
		colName = rnaIdNoLoc + "__alpha"

		# Set element in basal_prob to the transcription unit's value for alpha
		basal_prob[rna_idx] = r[col_names_to_index[colName]]

	# Convert to arrays
	deltaI, deltaJ, deltaV = np.array(deltaI), np.array(deltaJ), np.array(deltaV)
	delta_shape = (len(all_TUs), len(all_tfs))

	# Adjust any negative basal probabilities to 0
	basal_prob[basal_prob < 0] = 0

	# Add basal_prob vector and delta_prob matrix to sim_data
	transcription_regulation.basal_prob = basal_prob
	transcription_regulation.delta_prob = {
		"deltaI": deltaI,
		"deltaJ": deltaJ,
		"deltaV": deltaV,
		"shape": delta_shape,
		}


def crc32(*arrays: np.ndarray, initial: int = 0) -> int:
	"""Return a CRC32 checksum of the given ndarrays."""
	def crc_next(initial: int, array: np.ndarray) -> int:
		shape = str(array.shape).encode()
		values = array.tobytes()
		return binascii.crc32(values, binascii.crc32(shape, initial))

	return functools.reduce(crc_next, arrays, initial)


def setKmCooperativeEndoRNonLinearRNAdecay(sim_data, bulkContainer):
	"""
	Fits the affinities (Michaelis-Menten constants) for RNAs binding to
	endoRNAses.

	EndoRNAses perform the first step of RNA decay by cleaving a whole RNA
	somewhere inside its extent.  This results in RNA fragments, which are then
	digested into monomers by exoRNAses. To model endoRNAse activity, we need to
	determine an affinity (Michaelis-Menten constant) for each RNA that is
	consistent with experimentally observed half-lives.  The Michaelis-Menten
	constants must be determined simultaneously, as the RNAs must compete for
	the active site of the endoRNAse. (See the RnaDegradation Process class for
	more information about the dynamical model.) The parameters are estimated
	using a root solver (scipy.optimize.fsolve).  (See the
	sim_data.process.rna_decay.kmLossFunction method for more information about
	the optimization problem.)

	Requires
	--------
	- cell density, dry mass fraction, and average initial dry mass
		Used to calculate the cell volume, which in turn is used to calculate
		concentrations.
	- observed RNA degradation rates (half-lives)
	- endoRNAse counts
	- endoRNAse catalytic rate constants
	- RNA counts
	- boolean options that enable sensitivity analyses (see Notes below)

	Modifies
	--------
	- Michaelis-Menten constants for first-order decay (initially set to zeros)
	- Several optimization-related values
		Sensitivity analyses (optional, see Notes below)
		Terminal values for optimization-related functions

	Returns
	-------
	- enoRNAse Km values, in units of M

	Notes
	-----
	If certain options are set, a sensitivity analysis will be performed using a
	range of metaparameters. Outputs will be cached and utilized instead of
	running the optimization if possible.
	The function that generates the optimization functions is defined under
	sim_data but has no dependency on sim_data, and therefore could be moved
	here or elsewhere. (TODO)

	TODO (John): Refactor as a pure function.
	TODO (John): Why is this function called 'cooperative'?  It seems to instead
		assume and model competitive binding.
	TODO (John): Determine what part (if any) of the 'linear' parameter fitting
		should be retained.
	"""

	def arrays_differ(a: np.ndarray, b: np.ndarray) -> bool:
		return a.shape != b.shape or not np.allclose(a, b, equal_nan=True)

	cellDensity = sim_data.constants.cell_density
	cellVolume = sim_data.mass.avg_cell_dry_mass_init / cellDensity / sim_data.mass.cell_dry_mass_fraction
	countsToMolar = 1 / (sim_data.constants.n_avogadro * cellVolume)

	degradable_rna_ids = np.concatenate((
		sim_data.process.transcription.rna_data['id'],
		sim_data.process.transcription.mature_rna_data['id']
	))
	degradation_rates = (1 / units.s) * np.concatenate((
		sim_data.process.transcription.rna_data['deg_rate'].asNumber(1/units.s),
		sim_data.process.transcription.mature_rna_data['deg_rate'].asNumber(1/units.s)
		))
	endoRNaseConc = countsToMolar * bulkContainer.counts(
		sim_data.process.rna_decay.endoRNase_ids)
	kcatEndoRNase = sim_data.process.rna_decay.kcats
	totalEndoRnaseCapacity = units.sum(endoRNaseConc * kcatEndoRNase)

	endoRnaseRnaIds = sim_data.molecule_groups.endoRNase_rnas
	isEndoRnase = np.array([
		(x in endoRnaseRnaIds) for x in degradable_rna_ids])

	rna_counts = bulkContainer.counts(degradable_rna_ids)
	rna_conc = countsToMolar * rna_counts
	Km_counts = ((1 / degradation_rates * totalEndoRnaseCapacity) - rna_conc).asNumber()
	sim_data.process.rna_decay.Km_first_order_decay = Km_counts

	# Residuals can be written as follows: Res = f(Km) = 0, then Km = g(Km)
	# Compute derivative g(Km) in counts:
	KmQuadratic = 1 / np.power((1 / countsToMolar * Km_counts).asNumber(), 2)
	denominator = np.power(np.sum(rna_counts / (1 / countsToMolar * Km_counts).asNumber()),2)
	numerator = (1 / countsToMolar * totalEndoRnaseCapacity).asNumber() * (denominator - (rna_counts / (1 / countsToMolar * Km_counts).asNumber()))
	gDerivative = np.abs(KmQuadratic * (1 - (numerator / denominator)))
	if VERBOSE: print("Max derivative (counts) = %f" % max(gDerivative))

	# Compute derivative g(Km) in concentrations:
	KmQuadratic = 1 / np.power(Km_counts, 2)
	denominator = np.power(np.sum(rna_conc.asNumber() / Km_counts),2)
	numerator = totalEndoRnaseCapacity.asNumber() * (denominator - (rna_conc.asNumber() / Km_counts))
	gDerivative = np.abs(KmQuadratic * (1 - (numerator / denominator)))
	if VERBOSE: print("Max derivative (concentration) = %f" % max(gDerivative))

	# Sensitivity analysis: alpha (regularization term)
	Alphas = []
	if sim_data.constants.sensitivity_analysis_alpha:
		Alphas = [0.0001, 0.001, 0.01, 0.1, 1, 10]

	total_endo_rnase_capacity_mol_l_s = totalEndoRnaseCapacity.asNumber(units.mol / units.L / units.s)
	rna_conc_mol_l = rna_conc.asNumber(units.mol / units.L)
	degradation_rates_s = degradation_rates.asNumber(1 / units.s)

	for alpha in Alphas:
		if VERBOSE: print('Alpha = %f' % alpha)

		LossFunction, Rneg, R, LossFunctionP, R_aux, L_aux, Lp_aux, Jacob, Jacob_aux = sim_data.process.rna_decay.km_loss_function(
				total_endo_rnase_capacity_mol_l_s,
				rna_conc_mol_l,
				degradation_rates_s,
				isEndoRnase,
				alpha)
		Km_cooperative_model = scipy.optimize.fsolve(LossFunction, Km_counts, fprime = LossFunctionP)
		sim_data.process.rna_decay.sensitivity_analysis_alpha_residual[alpha] = np.sum(np.abs(R_aux(Km_cooperative_model)))
		sim_data.process.rna_decay.sensitivity_analysis_alpha_regulari_neg[alpha] = np.sum(np.abs(Rneg(Km_cooperative_model)))

	alpha = 0.5

	# Sensitivity analysis: kcatEndoRNase
	kcatEndo = []
	if sim_data.constants.sensitivity_analysis_kcat_endo:
		kcatEndo = [0.0001, 0.001, 0.01, 0.1, 1, 10]

	for kcat in kcatEndo:
		if VERBOSE: print('Kcat = %f' % kcat)

		totalEndoRNcap = units.sum(endoRNaseConc * kcat)
		LossFunction, Rneg, R, LossFunctionP, R_aux, L_aux, Lp_aux, Jacob, Jacob_aux = sim_data.process.rna_decay.km_loss_function(
				totalEndoRNcap.asNumber(units.mol / units.L),
				rna_conc_mol_l,
				degradation_rates_s,
				isEndoRnase,
				alpha)
		KmcountsIni = (( totalEndoRNcap / degradation_rates.asNumber() ) - rna_conc).asNumber()
		Km_cooperative_model = scipy.optimize.fsolve(LossFunction, KmcountsIni, fprime = LossFunctionP)
		sim_data.process.rna_decay.sensitivity_analysis_kcat[kcat] = Km_cooperative_model
		sim_data.process.rna_decay.sensitivity_analysis_kcat_res_ini[kcat] = np.sum(np.abs(R_aux(Km_counts)))
		sim_data.process.rna_decay.sensitivity_analysis_kcat_res_opt[kcat] = np.sum(np.abs(R_aux(Km_cooperative_model)))

	# Loss function, and derivative
	LossFunction, Rneg, R, LossFunctionP, R_aux, L_aux, Lp_aux, Jacob, Jacob_aux = sim_data.process.rna_decay.km_loss_function(
			total_endo_rnase_capacity_mol_l_s,
			rna_conc_mol_l,
			degradation_rates_s,
			isEndoRnase,
			alpha)

	# The checksum in the filename picks independent caches for distinct cases
	# such as different Parca options or Parca code in different git branches.
	# `make clean` will delete the cache files.
	needToUpdate = ''
	cache_dir = filepath.makedirs(filepath.ROOT_PATH, "cache")
	checksum = crc32(Km_counts, isEndoRnase, np.array(alpha))
	km_filepath = os.path.join(cache_dir, f'parca-km-{checksum}.cPickle')

	if os.path.exists(km_filepath):
		with open(km_filepath, "rb") as f:
			Km_cache = pickle.load(f)

		# KmCooperativeModel fits a set of Km values to give the expected degradation rates.
		# It takes 1.5 - 3 minutes to recompute.
		# R_aux calculates the difference of the degradation rate based on these
		# Km values and the expected rate so this sum seems like a good test of
		# whether the cache fits current input data, but cross-check additional
		# inputs to avoid Issue #996.
		Km_cooperative_model = Km_cache['Km_cooperative_model']
		if (Km_counts.shape != Km_cooperative_model.shape
				or np.sum(np.abs(R_aux(Km_cooperative_model))) > 1e-15
				or arrays_differ(Km_cache['total_endo_rnase_capacity_mol_l_s'], total_endo_rnase_capacity_mol_l_s)
				or arrays_differ(Km_cache['rna_conc_mol_l'], rna_conc_mol_l)
				or arrays_differ(Km_cache['degradation_rates_s'], degradation_rates_s)):
			needToUpdate = 'recompute'
	else:
		needToUpdate = 'compute'

	if needToUpdate:
		if VERBOSE: print(f'Running non-linear optimization to {needToUpdate} {km_filepath}')
		Km_cooperative_model = scipy.optimize.fsolve(LossFunction, Km_counts, fprime = LossFunctionP)
		Km_cache = dict(
			Km_cooperative_model=Km_cooperative_model,
			total_endo_rnase_capacity_mol_l_s=total_endo_rnase_capacity_mol_l_s,
			rna_conc_mol_l=rna_conc_mol_l,
			degradation_rates_s=degradation_rates_s)

		with open(km_filepath, "wb") as f:
			pickle.dump(Km_cache, f, protocol=pickle.HIGHEST_PROTOCOL)
	else:
		if VERBOSE:
			print("Not running non-linear optimization--using cached result {}".format(km_filepath))

	if VERBOSE > 1:
		print("Loss function (Km inital) = %f" % np.sum(np.abs(LossFunction(Km_counts))))
		print("Loss function (optimized Km) = %f" % np.sum(np.abs(LossFunction(Km_cooperative_model))))
		print("Negative km ratio = %f" % np.sum(np.abs(Rneg(Km_cooperative_model))))
		print("Residuals (Km initial) = %f" % np.sum(np.abs(R(Km_counts))))
		print("Residuals optimized = %f" % np.sum(np.abs(R(Km_cooperative_model))))
		print("EndoR residuals (Km initial) = %f" % np.sum(np.abs(isEndoRnase * R(Km_counts))))
		print("EndoR residuals optimized = %f" % np.sum(np.abs(isEndoRnase * R(Km_cooperative_model))))
		print("Residuals (scaled by Kdeg * RNAcounts) Km initial = %f" % np.sum(np.abs(R_aux(Km_counts))))
		print("Residuals (scaled by Kdeg * RNAcounts) optimized = %f" % np.sum(np.abs(R_aux(Km_cooperative_model))))

	# Evaluate Jacobian around solutions (Kmcounts and KmCooperativeModel)
	JacobDiag = np.diag(Jacob(Km_cooperative_model))
	Jacob_auxDiag = np.diag(Jacob_aux(Km_cooperative_model))

	# Compute convergence of non-linear optimization: g'(Km)
	Gkm = np.abs(1. - JacobDiag)
	Gkm_aux = np.abs(1. - Jacob_auxDiag)
	sim_data.process.rna_decay.Km_convergence = Gkm_aux

	# Convergence is guaranteed if g'(Km) <= K < 1
	if VERBOSE:
		print("Convergence (Jacobian) = %.0f%% (<K> = %.5f)" % (len(Gkm[Gkm < 1.]) / float(len(Gkm)) * 100., np.mean(Gkm)))
		print("Convergence (Jacobian_aux) = %.0f%% (<K> = %.5f)" % (len(Gkm_aux[Gkm_aux < 1.]) / float(len(Gkm_aux)) * 100., np.mean(Gkm_aux[Gkm_aux < 1.])))

	# Save statistics KM optimization
	sim_data.process.rna_decay.stats_fit['LossKm'] = np.sum(np.abs(LossFunction(Km_counts)))
	sim_data.process.rna_decay.stats_fit['LossKmOpt'] = np.sum(np.abs(LossFunction(Km_cooperative_model)))
	sim_data.process.rna_decay.stats_fit['RnegKmOpt'] = np.sum(np.abs(Rneg(Km_cooperative_model)))
	sim_data.process.rna_decay.stats_fit['ResKm'] = np.sum(np.abs(R(Km_counts)))
	sim_data.process.rna_decay.stats_fit['ResKmOpt'] = np.sum(np.abs(R(Km_cooperative_model)))
	sim_data.process.rna_decay.stats_fit['ResEndoRNKm'] = np.sum(np.abs(isEndoRnase * R(Km_counts)))
	sim_data.process.rna_decay.stats_fit['ResEndoRNKmOpt'] = np.sum(np.abs(isEndoRnase * R(Km_cooperative_model)))
	sim_data.process.rna_decay.stats_fit['ResScaledKm'] = np.sum(np.abs(R_aux(Km_counts)))
	sim_data.process.rna_decay.stats_fit['ResScaledKmOpt'] = np.sum(np.abs(R_aux(Km_cooperative_model)))

	return units.mol / units.L * Km_cooperative_model
