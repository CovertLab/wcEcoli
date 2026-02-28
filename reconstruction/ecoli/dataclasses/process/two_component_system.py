"""
Two component systems.

Note: Ligand binding to histidine kinases is modeled by equilibrium.

TODOs:
moleculesToNextTimeStep()
	Consider relocating (since it's useful for both the parca and simulation)
"""

import numpy as np
import scipy
import scipy.integrate
import re

import sympy as sp

from wholecell.utils import build_ode
from wholecell.utils import data
from wholecell.utils import units


# Alternative methods to try (in order of priority) when solving ODEs to the next time step
IVP_METHODS = ['LSODA', 'BDF']

# TEMP ADDITION:
class TwoComponentSystemError(Exception):
	pass

class MoleculeNotFoundError(TwoComponentSystemError):
	pass
# END TEMP ADDITION

class TwoComponentSystem(object):
	def __init__(self, raw_data, sim_data):
		# Store two component system raw data for use in analysis
		sim_data.molecule_groups.twoComponentSystems = raw_data.two_component_systems

		# Build the abstractions needed for two component systems
		molecules = []  # list of all molecules involved in two component system
		moleculeTypes = []  # the type of each molecule (metabolite, protein monomer, protein complex, etc.)

		ratesFwd = []  # rate of reaction fwd
		ratesRev = []  # rate of reaction reverse (most/all are 0 in flat file)
		rxnIds = []  # ID tied to each rxn equation

		stoichMatrixI = []  # Molecule indices
		stoichMatrixJ = []  # Reaction indices
		stoichMatrixV = []  # Stoichometric coefficients

		stoichMatrixMass = []  # molecular mass of molecules in stoichMatrixI

		independentMolecules = []  # list of all specific independent molecule names
		independent_molecule_indexes = []  # index of each of the independent molecules
		independentToDependentMolecules = {}  # holds the phosphorylated version of the independent molecules

		activeToInactiveTF = {}  # convention: active TF is the DNA-binding form (active form is phosphorylated version of RR)

		# TEMP:
		common_molecules = []

		# Build template reactions
		signalingTemplate = {
			1: ["POS-LIGAND-BOUND-HK-PHOSPHORYLATION_RXN",
				"POS-LIGAND-BOUND-HK-PHOSPHOTRANSFER_RXN",
				"POS-RR-DEPHOSPHORYLATION_RXN",
				"POS-HK-PHOSPHORYLATION_RXN",
				"POS-HK-PHOSPHOTRANSFER_RXN",
				],
			-1: ["NEG-HK-PHOSPHORYLATION_RXN",
				"NEG-HK-PHOSPHOTRANSFER_RXN",
				"NEG-RR-DEPHOSPHORYLATION_RXN",
				],
			}

		reactionTemplate = {}
		for reactionIndex, reaction in enumerate(raw_data.two_component_system_templates):
			reactionTemplate[str(reaction["id"])] = reaction


		# Build stoichiometry matrix
		for systemIndex, system in enumerate(raw_data.two_component_systems):
			for reaction in signalingTemplate[system["orientation"]]:
				reactionName = self.get_reaction_name(reaction, system["molecules"])

				if reactionName not in rxnIds:
					rxnIds.append(reactionName)
					ratesFwd.append(reactionTemplate[reaction]["forward_rate"])
					ratesRev.append(reactionTemplate[reaction]["reverse_rate"])
					reactionIndex = len(rxnIds) - 1

				else:
					reactionIndex = rxnIds.index(reactionName)

				for molecule in reactionTemplate[reaction]["stoichiometry"]:

					# Build name for system molecules
					if molecule["molecule"] in system["molecules"]:
						moleculeName = "{}[{}]".format(
							system["molecules"][molecule["molecule"]],
							sim_data.getter.get_compartment(system["molecules"][molecule["molecule"]])[0]
							)

					# Build name for common molecules (ATP, ADP, PI, WATER, PROTON)
					else:
						moleculeName = "{}[{}]".format(
							molecule["molecule"],
							molecule["location"]
							)
						if moleculeName not in common_molecules:
							common_molecules.append(moleculeName)

					if moleculeName not in molecules:
						molecules.append(moleculeName)
						moleculeIndex = len(molecules) - 1
						moleculeTypes.append(str(molecule["molecule"]))

					else:
						moleculeIndex = molecules.index(moleculeName)

					coefficient = molecule["coeff"]

					# Store indices for the row and column, and molecule coefficient for building the stoichiometry matrix
					stoichMatrixI.append(moleculeIndex)
					stoichMatrixJ.append(reactionIndex)
					stoichMatrixV.append(coefficient)

					# Build matrix with linearly independent rows based on network orientation
					if str(molecule["molecule"]) in ["HK", "RR", "ATP"] and moleculeName not in independentMolecules:
						independentMolecules.append(moleculeName)
						independent_molecule_indexes.append(moleculeIndex)

						# Map linearly independent molecules (rows) to their dependents (phosphorylated forms of histidine kinases and response regulators)
						if str(molecule["molecule"]) != "ATP":
							independentToDependentMolecules[moleculeName] = "{}[{}]".format(
								system["molecules"]["PHOSPHO-" + str(molecule["molecule"])],
								molecule["location"]
								)

						# Map active transcription factors (phosphorylated response regulators) to their inactive forms (unphosphorylated response regulators)
						if str(molecule["molecule"]) == "RR":
							activeTF = "{}[{}]".format(
								system["molecules"]["PHOSPHO-RR"],
								molecule["location"]
								)
							activeToInactiveTF[activeTF] = moleculeName

					# Account for ligand-bound histidine kinases for positively oriented networks
					if system["orientation"] == 1:
						if str(molecule["molecule"]) == "HK-LIGAND" and moleculeName not in independentMolecules:
							independentMolecules.append(moleculeName)
							independent_molecule_indexes.append(moleculeIndex)

							# Map the linearly independent ligand-bound histidine kinases to their dependents (phosphorylated forms of ligand-bound histidine kinases)
							independentToDependentMolecules[moleculeName] = "{}[{}]".format(
								system["molecules"]["PHOSPHO-" + str(molecule["molecule"])],
								molecule["location"]
								)

					# Find molecular mass
					molecularMass = sim_data.getter.get_mass(moleculeName).asNumber(units.g / units.mol)
					stoichMatrixMass.append(molecularMass)

		# TODO(jerry): Move most of the rest to a subroutine for __init__ and __setstate__?
		self._stoichMatrixI = np.array(stoichMatrixI)  # array of molecule indices
		self._stoichMatrixJ = np.array(stoichMatrixJ)  # array of reaction indices
		self._stoichMatrixV = np.array(stoichMatrixV)  # arrary of stoichometric coefficients

		self.molecule_names = np.array(molecules, dtype='U')
		self.molecule_types = np.array(moleculeTypes, dtype='U')
		self.rxn_ids = rxnIds
		self.rates_fwd = np.array(ratesFwd)
		self.rates_rev = np.array(ratesRev)

		self.independent_molecules = np.array(independentMolecules, dtype='U')
		self.independent_molecule_indexes = np.array(independent_molecule_indexes)
		self.independent_to_dependent_molecules = independentToDependentMolecules

		self.independent_molecules_atp_index = np.where(self.independent_molecules == "ATP[c]")[0][0]

		# TODO (mia): figure out why raw_data.modified_proteins is used here to
		#  generate this dictionary , as it incorperates molecules with
		#  compartment tag combinations that are not "valid" according to
		#  is_valid_molecule(), thus, causing stoich_matrix_monomers()
		#  built in this listener to not include certain molecules and not match
		#  the stoich_matrix() built using the template above well.
		self.complex_to_monomer = self._buildComplexToMonomer(raw_data.modified_proteins, self.molecule_names)
		self.molecules_to_parent_complexes_dict = {}
		self.molecules_to_all_downstream_complexes_dict = {}


		# Mass balance matrix
		self._stoich_matrix_mass = np.array(stoichMatrixMass)
		self.balance_matrix = self.stoich_matrix() * self.mass_matrix()

		# Find the mass balance of each equation in the balanceMatrix
		massBalanceArray = self.mass_balance()

		# The stoichometric matrix should balance out to numerical zero.
		assert np.max(np.absolute(massBalanceArray)) < 1e-9

		# Map active TF to inactive TF
		self.active_to_inactive_tf = activeToInactiveTF

		# Build matrices
		self._populate_derivative_and_jacobian()
		self.dependency_matrix = self._make_dependency_matrix()

		# Molecules that are required to produce ATP with the independent stoich matrix
		self.atp_reaction_reactant_mask = self.dependency_matrix[:, self.independent_molecules_atp_index] < 0

		# TODO (mia): test out altering stoich_matrix_monomers to be built the
		#  same way as it is in the complexation and equilibrium process files.
		# Generate the stoich matrix for complexes to all consituent monomers:
		stoich_matrix_monomers = self.stoich_matrix_monomers()

		# Generate dictionary mapping molecules to the direct parent complexes they form:
		for subunit in self.molecule_names:
			# find the matrix index where this subunit is as a molecule:
			subunit_index = np.where(self.molecule_names == subunit)

			# Find the indicies of self._stoich_matrix_J where the value will
			# correspond to the index of the correct reaction in self.ids_reactions
			reaction_indicies = np.where(
				(self._stoichMatrixI == subunit_index) &
				(self._stoichMatrixV < 0))[0]

			# For each reaction index, find the complex(es) that is(are) formed:
			parent_complexes = {}
			for reaction_idx in reaction_indicies:
				# find the value of the reaction index (which will be the index
				# that corresponds to the index of the reaction in self.ids_reactions):
				rxn_idx = self._stoichMatrixJ[reaction_idx]

				# Initialize data structures to hold complex information
				complex_information = []
				stoich = {}
				stoich_known = {}
				complex_type = {}
				reaction_name = {}

				# Find the indexes of the complexes formed in this reaction:
				complex_indices = np.where(
					(self._stoichMatrixJ == rxn_idx) &
					(self._stoichMatrixV > 0))[0]

				for complex_index in complex_indices:
					# Check if the complex is in the complex_to_monomer dictionary
					# (i.e. if it is a complex with known subunit composition):
					complex_name = self.molecule_names[self._stoichMatrixI[complex_index]]
					if complex_name not in self.complex_to_monomer.keys():
						# If the product is not a complex, skip it.
						continue
					else:
						# Find the number of unique subunits in this complex reaction:
						cplx_idx = list(self.complex_to_monomer.keys()).index(complex_name)
						components = np.where(stoich_matrix_monomers[:, cplx_idx] < 0)[0]
						num_unique_subunits = len(components)
						# NOTE: the complex_to_monomer dictionary output might
						# not be aligned with the stoich_matrix() currently.
						# Ideally, it would be best to confirm the modified_proteins
						# should even be used ot generate complex_to_monomer.

						if num_unique_subunits > 1:
							cplx_type = 'heterogeneous'
						else:
							cplx_type = 'homogeneous'

						# Add complex information to lists
						complex_name = self.molecule_names[self._stoichMatrixI[complex_index]]
						reaction_name['reaction_id'] = self.rxn_ids[rxn_idx]
						stoich['stoichiometry'] = self._stoichMatrixV[reaction_idx]
						stoich_known['stoich_unknown'] = 'not applicable'
						complex_type['complex_type'] = cplx_type
						complex_information.append(reaction_name)
						complex_information.append(stoich)
						complex_information.append(stoich_known)
						complex_information.append(complex_type)

						# Append the complex name and stoich as a dictionary entry
						parent_complexes[complex_name] = complex_information

			self.molecules_to_parent_complexes_dict[subunit] = parent_complexes


		# Make a dictionary mapping molecules to all downstream complexes they form
		# (both directly and indirectly via another complex):
		for subunit in self.molecule_names:
			# find the matrix index where this subunit is as a molecule:
			subunit_index = np.where(self.molecule_names == subunit)

			# Find the indices of complexes containing the subunit:
			complex_indices = np.where(stoich_matrix_monomers[subunit_index, :] < 0)[0]

			downstream_complexes = {}
			for complex_idx in complex_indices:
				# Find the complex's name:
				complex_name = list(self.complex_to_monomer.keys())[complex_idx]

				# Obtain the index of the complex within self.molecule_names:
				cplx_idx = np.where(self.molecule_names == complex_name)

				# Use the stoichMatrix() to find the reaction index:
				reaction_indices = np.where(
					(self._stoichMatrixI == cplx_idx) & (self._stoichMatrixV > 0))[0]

				# Obtain the value that corresponds to the index of the reaction in self.ids_reactions:
				reaction_idx = self._stoichMatrixJ[reaction_indices]

				# Initialize data structures to hold complex information
				downstream_complex_information = []
				stoich = {}
				stoich_known = {}
				complex_type = {}
				reaction_name = {}

				# Find the number of unique subunits in this complex reaction:
				unique_subunits_in_complex = np.where(
					(self._stoichMatrixJ == reaction_idx) &
					(self._stoichMatrixV < 0))[0]
				num_unique_subunits = len(unique_subunits_in_complex)
				if num_unique_subunits > 1:
					cplx_type = 'heterogeneous'
				elif num_unique_subunits == 1:
					cplx_type = 'homogeneous'
				else:
					cplx_type = 'unknown'

				# Add complex information to lists:
				reaction_name['reaction_id'] = self.rxn_ids[reaction_idx[0]]
				stoich['stoichiometry'] = stoich_matrix_monomers[subunit_index, complex_idx]
				stoich_known['stoich_unknown'] = 'not applicable'
				complex_type['complex_type'] = cplx_type
				downstream_complex_information.append(reaction_name)
				downstream_complex_information.append(stoich)
				downstream_complex_information.append(stoich_known)
				downstream_complex_information.append(complex_type)

				# Append the complex name and stoich as a dictionary entry
				downstream_complexes[complex_name] = downstream_complex_information

			self.molecules_to_all_downstream_complexes_dict[subunit] = downstream_complexes

		# TEST making a different SMM
		self._stoichMatrix = self.stoich_matrix()
		self.modified_molecules = self.make_modified_molecule_names(sim_data)
		SMM = self.stoich_matrix_monomers_TEMP(sim_data)



	def __getstate__(self):
		"""Return the state to pickle, omitting derived attributes that
		__setstate__() will recompute, esp. the ode_derivatives
		that don't pickle.
		"""
		return data.dissoc_strict(self.__dict__, (
			'symbolic_rates', 'symbolic_rates_jacobian',
			'derivatives_parca_symbolic', 'derivatives_parca_jacobian_symbolic',
			'_rates', '_rates_jacobian',
			'derivatives_parca', 'derivatives_parca_jacobian',
			'dependency_matrix', '_stoich_matrix'))

	def __setstate__(self, state):
		"""Restore instance attributes, recomputing some of them."""
		self.__dict__.update(state)
		self._populate_derivative_and_jacobian()
		self.dependency_matrix = self._make_dependency_matrix()

	def _buildComplexToMonomer(self, modifiedFormsMonomers, tcsMolecules):
		'''
		Maps each complex to a dictionary that maps each subunit of the complex to its stoichiometry
		'''
		D = {}
		for row in modifiedFormsMonomers:
			# tags on the molecule compartment found in tcsMolecules
			molecule_and_location = f"{row['id']}[{row['compartment']}]"
			# filter by molecules actively being used in the simulation:
			if molecule_and_location in tcsMolecules:
				D[molecule_and_location] = {}
				for subunit in row["subunits"]:
					# TODO (mia): NOTE: I think some of the hardcoded subunit compartment tags are incorrect. consider using getter.get_compartment()
					# We only care about mapping to protein monomers for now
					# and PI[c] stoichiometry is off for some complexes so we
					# can skip it for now (see #975)
					if subunit['monomer'] == 'PI[c]':
						continue
					D[molecule_and_location][str(subunit["monomer"])] = float(subunit["stoichiometry"])

		return D


	def stoich_matrix(self):
		'''
		Builds stoichiometry matrix
		Rows: molecules
		Columns: reactions
		Values: reaction stoichiometry
		'''
		shape = (self._stoichMatrixI.max()+1, self._stoichMatrixJ.max()+1)
		out = np.zeros(shape, np.float64)
		out[self._stoichMatrixI, self._stoichMatrixJ] = self._stoichMatrixV
		return out


	def mass_matrix(self):
		'''
		Builds stoichiometry mass matrix
		Rows: molecules
		Columns: reactions
		Values: molecular mass
		'''
		shape = (self._stoichMatrixI.max()+1, self._stoichMatrixJ.max()+1)
		out = np.zeros(shape, np.float64)
		out[self._stoichMatrixI, self._stoichMatrixJ] = self._stoich_matrix_mass
		return out


	def mass_balance(self):
		'''
		Sum along the columns of the massBalance matrix to check for reaction mass balance
		'''
		return np.sum(self.balance_matrix, axis=0)


	def stoich_matrix_monomers(self):
		'''
		Builds stoichiometry matrix for monomers (complex subunits)
		Rows: molecules (complexes and monomers)
		Columns: complexes
		Values: monomer stoichiometry
		'''
		# TODO (mia): figure out if this can be replaced by finding a way to
		#  get a list of the complex_IDs similar to how they are generated in
		#  the other complexation process files, as the current method relies
		#  on complex_to_monomer, which is not built using the same reactions
		#  as the stoich_matrix() matrix is.
		ids_complexes = self.complex_to_monomer.keys()
		stoichMatrixMonomersI = []
		stoichMatrixMonomersJ = []
		stoichMatrixMonomersV = []
		for colIdx, id_complex in enumerate(ids_complexes):
			D = self.get_monomers(id_complex)

			rowIdx = self.molecule_names.tolist().index(id_complex)
			stoichMatrixMonomersI.append(rowIdx)
			stoichMatrixMonomersJ.append(colIdx)
			stoichMatrixMonomersV.append(1.)

			for subunitId, subunitStoich in zip(D["subunitIds"], D["subunitStoich"]):
				# NOTE (mia): since the complex_to_monomer dictionary does not
				# match how molecule_names is generated, many subunits will not
				# be appended here.
				if subunitId in self.molecule_names.tolist():
					rowIdx = self.molecule_names.tolist().index(subunitId)
					stoichMatrixMonomersI.append(rowIdx)
					stoichMatrixMonomersJ.append(colIdx)
					stoichMatrixMonomersV.append(-1. * subunitStoich)

		stoichMatrixMonomersI = np.array(stoichMatrixMonomersI)
		stoichMatrixMonomersJ = np.array(stoichMatrixMonomersJ)
		stoichMatrixMonomersV = np.array(stoichMatrixMonomersV)
		shape = (stoichMatrixMonomersI.max() + 1, stoichMatrixMonomersJ.max() + 1)

		out = np.zeros(shape, np.float64)
		out[stoichMatrixMonomersI, stoichMatrixMonomersJ] = stoichMatrixMonomersV
		return out


	def _populate_derivative_and_jacobian(self):
		'''Compile callable functions for computing the derivative and the Jacobian.'''
		self._make_derivative()
		self._make_derivative_parca()

		self._rates = build_ode.derivatives(self.symbolic_rates)
		self._rates_jacobian = build_ode.derivatives_jacobian(self.symbolic_rates_jacobian)
		self._stoich_matrix = self.stoich_matrix()  # Matrix is small and can be cached for derivatives

		# WORKAROUND: Avoid Numba LoweringError JIT-compiling these functions by selecting the
		#   non-JIT versions (at index 0):
		self.derivatives_parca = build_ode.derivatives(self.derivatives_parca_symbolic)[0]
		self.derivatives_parca_jacobian = build_ode.derivatives_jacobian(self.derivatives_parca_jacobian_symbolic)[0]


	def _make_y_dy(self):
		S = self.stoich_matrix()

		yStrings = ["y[%d]" % x for x in range(S.shape[0])]
		y = sp.symbols(yStrings)

		rates = []
		for colIdx in range(S.shape[1]):
			negIdxs = np.where(S[:, colIdx] < 0)[0]
			posIdxs = np.where(S[:, colIdx] > 0)[0]

			reactantFlux = self.rates_fwd[colIdx]
			for negIdx in negIdxs:
				stoich = -S[negIdx, colIdx]
				if stoich == 1:
					reactantFlux *= y[negIdx]
				else:
					reactantFlux *= y[negIdx]**stoich

			productFlux = self.rates_rev[colIdx]
			for posIdx in posIdxs:
				stoich = S[posIdx, colIdx]
				if stoich == 1:
					productFlux *= y[posIdx]
				else:
					productFlux *= y[posIdx]**stoich

			rates.append(reactantFlux - productFlux)

		return y, rates


	def _make_derivative(self):
		'''
		Creates symbolic representation of the ordinary differential equations
		and the Jacobian. Used during simulations.
		'''
		y, rates = self._make_y_dy()

		rates = sp.Matrix(rates)
		J = rates.jacobian(y)

		self.symbolic_rates = rates
		self.symbolic_rates_jacobian = J


	def _make_derivative_parca(self):
		'''
		Creates symbolic representation of the ordinary differential equations
		and the Jacobian assuming ATP, ADP, Pi, water and protons are at
		steady state. Used in the parca.
		'''
		y, rates = self._make_y_dy()
		dy = self.stoich_matrix().dot(rates)

		# Metabolism will keep these molecules at steady state
		constantMolecules = ["ATP[c]", "ADP[c]", "Pi[c]", "WATER[c]", "PROTON[c]"]
		for molecule in constantMolecules:
			moleculeIdx = np.where(self.molecule_names == molecule)[0][0]
			dy[moleculeIdx] = sp.S.Zero

		dy = sp.Matrix(dy)
		J = dy.jacobian(y)

		self.derivatives_parca_jacobian_symbolic = J
		self.derivatives_parca_symbolic = dy


	def molecules_to_next_time_step(self, moleculeCounts, cellVolume,
			nAvogadro, timeStepSec, random_state, method="LSODA",
			min_time_step=None, jit=True, methods_tried=None):
		"""
		Calculates the changes in the counts of molecules in the next timestep
		by solving an initial value ODE problem.

		Args:
			moleculeCounts (1d ndarray, ints): current counts of molecules
				involved in the ODE
			cellVolume (float): current volume of cell
			nAvogadro (float): Avogadro's number
			timeStepSec (float): current length of timestep in seconds
			random_state (RandomState object): process random state
			method (str): name of the ODE method to use
			min_time_step (int): if not None, timeStepSec will be scaled down until
				it is below min_time_step if negative counts are encountered
			jit (bool): if True, use the jit compiled version of derivatives
				functions
			methods_tried (Optional[Set[str]]): methods for the solver that have
				already been tried

		Returns:
			moleculesNeeded (1d ndarray, ints): counts of molecules that need
				to be consumed
			allMoleculesChanges (1d ndarray, ints): expected changes in
				molecule counts after timestep
		"""
		y_init = moleculeCounts / (cellVolume * nAvogadro)

		# In this version of SciPy, solve_ivp does not support args so need to
		# select the derivatives functions to use. Could be simplified to single
		# functions that take a jit argument from solve_ivp in the future.
		if jit:
			derivatives = self.derivatives_jit
			derivatives_jacobian = self.derivatives_jacobian_jit
		else:
			derivatives = self.derivatives
			derivatives_jacobian = self.derivatives_jacobian

		sol = scipy.integrate.solve_ivp(
			derivatives, [0, timeStepSec], y_init,
			method=method, t_eval=[0, timeStepSec], atol=1e-8,
			jac=derivatives_jacobian
			)
		y = sol.y.T

		# Handle negative counts by attempting to solve again with different options
		if np.any(y[-1, :] * (cellVolume * nAvogadro) <= -1e-3):
			if min_time_step and timeStepSec > min_time_step:
				# Call method again with a shorter time step until min_time_step is reached
				return self.molecules_to_next_time_step(
					moleculeCounts, cellVolume, nAvogadro, timeStepSec/2, random_state,
					method=method, min_time_step=min_time_step, jit=jit)

			# Try with different method for better stability
			if methods_tried is None:
				methods_tried = set()
			methods_tried.add(method)
			for new_method in IVP_METHODS:
				# Skip methods that have already been tried
				if new_method in methods_tried:
					continue

				print(f'Warning: switching to {new_method} method in TCS')
				return self.molecules_to_next_time_step(
					moleculeCounts, cellVolume, nAvogadro, timeStepSec, random_state,
					method=new_method, min_time_step=min_time_step, jit=jit,
					methods_tried=methods_tried)
			else:
				raise Exception(
					"Solution to ODE for two-component systems has negative values."
					)

		y[y < 0] = 0
		yMolecules = y * (cellVolume * nAvogadro)
		dYMolecules = yMolecules[-1, :] - yMolecules[0, :]

		independentMoleculesCounts = np.round(dYMolecules[self.independent_molecule_indexes])

		max_atp_rxns = moleculeCounts[self.atp_reaction_reactant_mask].min()
		# To ensure that we have non-negative counts of phosphate, we must
		# have the following (which can be seen from the dependency matrix)
		independentMoleculesCounts[self.independent_molecules_atp_index] = np.fmin(
			independentMoleculesCounts[:self.independent_molecules_atp_index].sum()
			+ independentMoleculesCounts[(self.independent_molecules_atp_index + 1):].sum(),
			max_atp_rxns
			)

		# Calculate changes in molecule counts for all molecules
		allMoleculesChanges = self.dependency_matrix.dot(independentMoleculesCounts)

		# Calculate molecules needed by assuming other molecules that would produce necessary
		# phosphate won't be allocated
		negative = independentMoleculesCounts.copy()
		negative[negative > 0] = 0
		negative[self.independent_molecules_atp_index] = (
			negative[:self.independent_molecules_atp_index].sum()
			+ negative[(self.independent_molecules_atp_index + 1):].sum()
			)
		moleculesNeeded = self.dependency_matrix.dot(-negative).clip(min=0)
		positive = independentMoleculesCounts.copy()
		positive[positive < 0] = 0
		moleculesNeeded += self.dependency_matrix.dot(-positive).clip(min=0)

		# Adjust molecules to prevent using more than allocated
		iteration = 0
		final_molecules = allMoleculesChanges + moleculeCounts
		while np.any(final_molecules < 0):
			stoich = self.stoich_matrix()
			mol_idx = np.where(final_molecules < 0)[0][0]
			rxns = stoich[mol_idx, :] < 0  # reactions that consume the molecule that has been depleted

			# Get products of reactions to turn back into reactants
			consuming_stoich = stoich[:, rxns]
			consuming_stoich[consuming_stoich < 0] = 0  # exclude molecules that are also consumed in these reactions
			consuming_stoich[consuming_stoich.sum(axis=1) > 1] = 0  # exclude molecules that are shared between reactions

			# Weight possible reactions by how different the rounded solution is to the integrated solution
			rxn_propensity = (allMoleculesChanges - dYMolecules).dot(consuming_stoich)
			rxn_propensity[rxn_propensity < 0] = 0
			rxn_propensity /= rxn_propensity.sum()

			# Sample from propensities to find reaction to reverse
			rxn = np.where(random_state.multinomial(1, rxn_propensity))[0][0]
			allMoleculesChanges -= stoich[:, rxns][:, rxn]
			final_molecules = allMoleculesChanges + moleculeCounts

			# Prevent possibility of infinite loop - should never need to reduce each reaction more than once
			iteration += 1
			if iteration > stoich.shape[1]:
				raise ValueError('Could not get positive molecule counts for {} in two_component_system'
					.format(self.molecule_names[mol_idx]))

		return moleculesNeeded, allMoleculesChanges


	def molecules_to_ss(self, moleculeCounts, cellVolume, nAvogadro, timeStepSec=1e20):
		"""
		Calculates the changes in the counts of molecules as the system
		reaches steady state

		Args:
			moleculeCounts: current counts of molecules involved in the ODE
			cellVolume: current volume of cell
			nAvogadro: Avogadro's number
			timeStepSec: current length of timestep (set to large number)

		Returns:
			moleculesNeeded: counts of molecules that need to be consumed
			allMoleculesChanges: expected changes in molecule counts after
			timestep
		"""
		# TODO (Gwanggyu): This function should probably get merged with the
		# 	function above.
		y_init = moleculeCounts / (cellVolume * nAvogadro)

		y = scipy.integrate.odeint(
			self.derivatives_parca, y_init,
			t=[0, timeStepSec], Dfun=self.derivatives_parca_jacobian
		)

		if np.any(y[-1, :] * (cellVolume * nAvogadro) <= -1):
			raise Exception(
				"Solution to ODE for two-component systems has negative values."
				)

		y[y < 0] = 0
		yMolecules = y * (cellVolume * nAvogadro)
		dYMolecules = yMolecules[-1, :] - yMolecules[0, :]

		independentMoleculesCounts = np.array(
			[np.round(dYMolecules[x]) for x in self.independent_molecule_indexes]
			)

		# To ensure that we have non-negative counts of phosphate, we must
		# have the following (which can be seen from the dependency matrix)
		independentMoleculesCounts[self.independent_molecules_atp_index] = (
			independentMoleculesCounts[:self.independent_molecules_atp_index].sum()
			+ independentMoleculesCounts[(self.independent_molecules_atp_index + 1):].sum()
			)

		# Calculate changes in molecule counts for all molecules
		allMoleculesChanges = np.dot(
			self.dependency_matrix, independentMoleculesCounts)

		moleculesNeeded = np.negative(allMoleculesChanges).clip(min=0)

		return moleculesNeeded, allMoleculesChanges


	def get_monomers(self, cplxId):
		'''
		Returns subunits for a complex (or any ID passed).
		If the ID passed is already a monomer returns the
		monomer ID again with a stoichiometric coefficient
		of zero.
		'''

		info = self.complex_to_monomer
		if cplxId in info:
			out = {
				'subunitIds': list(info[cplxId].keys()),
				'subunitStoich': list(info[cplxId].values())}
		else:
			out = {'subunitIds': cplxId, 'subunitStoich': 1}
		return out


	def get_reaction_name(self, templateName, systemMolecules):
		'''
		Returns reaction name for a particular system.
		'''
		startIndex = 0
		reactionName = templateName
		for endIndex in [x.start() for x in re.finditer("-", templateName)]:
			if templateName[startIndex:endIndex] in systemMolecules:
				reactionName = reactionName.replace(templateName[startIndex:endIndex], str(systemMolecules[templateName[startIndex:endIndex]]))
			startIndex = endIndex + 1

		return reactionName


	def _make_dependency_matrix(self):
		'''
		Builds matrix mapping linearly independent molecules (ATP, histidine kinases,
		response regulators, and ligand-bound histidine kinases for positively oriented
		networks) to their dependents.
		'''
		molecule_name_to_index = {
			name: i for (i, name) in enumerate(self.molecule_names)}

		dependencyMatrixI = []
		dependencyMatrixJ = []
		dependencyMatrixV = []
		dependencyMatrixATPJ = -1

		for independentMoleculeIndex, independentMoleculeId in enumerate(self.independent_molecule_indexes):
			dependencyMatrixI.append(independentMoleculeId)
			dependencyMatrixJ.append(independentMoleculeIndex)
			dependencyMatrixV.append(1)

			if self.molecule_names[independentMoleculeId] == "ATP[c]":
				dependencyMatrixATPJ = independentMoleculeIndex
			else:
				dependentMoleculeId = molecule_name_to_index[
					self.independent_to_dependent_molecules[self.molecule_names[independentMoleculeId]]]
				dependencyMatrixI.append(dependentMoleculeId)
				dependencyMatrixJ.append(independentMoleculeIndex)
				dependencyMatrixV.append(-1)

		# ATP dependents: ADP, PI, WATER, PROTON)
		for ATPdependent in ["ADP[c]", "Pi[c]", "WATER[c]", "PROTON[c]"]:
			dependencyMatrixI.append(molecule_name_to_index[ATPdependent])
			dependencyMatrixJ.append(dependencyMatrixATPJ)
			if ATPdependent == "WATER[c]":
				dependencyMatrixV.append(1)
			else:
				dependencyMatrixV.append(-1)

		for col in np.arange(self.independent_molecule_indexes.size):
			if col == dependencyMatrixATPJ:
				continue
			else:
				dependencyMatrixI.append(molecule_name_to_index["Pi[c]"])
				dependencyMatrixJ.append(col)
				dependencyMatrixV.append(1)

				dependencyMatrixI.append(molecule_name_to_index["WATER[c]"])
				dependencyMatrixJ.append(col)
				dependencyMatrixV.append(-1)

		shape = (np.max(dependencyMatrixI) +1, np.max(dependencyMatrixJ) +1)

		out = np.zeros(shape, np.float64)

		out[dependencyMatrixI, dependencyMatrixJ] = dependencyMatrixV
		return out

	def derivatives(self, t, y):
		"""
		Calculate derivatives from stoichiometry and rates with argument order
		for solve_ivp.
		"""
		return self._stoich_matrix.dot(self._rates[0](y, t))

	def derivatives_jacobian(self, t, y):
		"""
		Calculate the jacobian of derivatives from stoichiometry and rates
		with argument order for solve_ivp.
		"""
		return self._stoich_matrix.dot(self._rates_jacobian[0](y, t))

	def derivatives_jit(self, t, y):
		"""
		Calculate derivatives from stoichiometry and rates with argument order
		for solve_ivp.
		"""
		return self._stoich_matrix.dot(self._rates[1](y, t))

	def derivatives_jacobian_jit(self, t, y):
		"""
		Calculate the jacobian of derivatives from stoichiometry and rates
		with argument order for solve_ivp.
		"""
		return self._stoich_matrix.dot(self._rates_jacobian[1](y, t))

	# TODO: add temporary functions for generating a different version of the stoich matrix:

	def make_modified_molecule_names(self, sim_data):

		# orginal monomers:
		tcs_molecules = self.molecule_names.tolist()
		# get the new molecules from the complex_to_monomer dictionary:
		new_molecules = []
		for complex in self.complex_to_monomer.keys():
			# index to the individual subunit dictionaries within:
			subunits = self.get_monomers_TEMP(sim_data, complex)
			for subunit in subunits['subunitIds']:
				if subunit not in tcs_molecules and subunit not in new_molecules:
					new_molecules.append(subunit)
		return np.array(tcs_molecules + new_molecules)


	def get_monomers_TEMP(self, sim_data, cplxId):
		'''
		Returns subunits for a complex (or any ID passed).
		If the ID passed is already a monomer returns the
		monomer ID again with a stoichiometric coefficient
		of zero.
		'''

		info = self.complex_to_monomer
		if cplxId in info:
			subunits = []
			subunit_stoich = []
			for subunitId, subunitStoich in zip(info[cplxId].keys(), info[cplxId].values()):
				# Check that the compartment tags are correct and valid:
				if sim_data.getter.is_valid_molecule(subunitId):
					subunits.append(subunitId)
					subunit_stoich.append(subunitStoich)
				else:
					subunit_wo_tag = subunitId.split("[")[0]
					compartment_tag = sim_data.getter.get_compartment(subunit_wo_tag)
					new_subunitId = f"{subunit_wo_tag}[{compartment_tag[0]}]"
					subunits.append(new_subunitId)
					subunit_stoich.append(subunitStoich)

			out = {
				'subunitIds': list(subunits),
				'subunitStoich': list(subunit_stoich)}
		else:
			out = {'subunitIds': cplxId, 'subunitStoich': 1}
		return out

	def stoich_matrix_monomers_TEMP(self, sim_data):
		"""
		Builds a stoichiometric matrix where each column is a reaction that
		forms a complex directly from its constituent monomers. Since some
		reactions from the raw data are complexation reactions of complexes,
		this is different from the stoichiometric matrix generated by
		stoichMatrix().
		"""
		modified_molecule_names = self.make_modified_molecule_names(sim_data)
		ids_complexes = self.complex_to_monomer.keys()
		mol_names = list(modified_molecule_names)
		stoichMatrixMonomersI = []
		stoichMatrixMonomersJ = []
		stoichMatrixMonomersV = []

		for colIdx, id_complex in enumerate(ids_complexes):
			D = self.get_monomers_TEMP(sim_data, id_complex)
			rowIdx = mol_names.index(id_complex)
			stoichMatrixMonomersI.append(rowIdx)
			stoichMatrixMonomersJ.append(colIdx)
			stoichMatrixMonomersV.append(1.)

			# for some reason, some of these are not coming out as lists, so this handles that:
			subunitIds = D["subunitIds"] if isinstance(D["subunitIds"], list) else [
				D["subunitIds"]]
			subunitStoichs = D["subunitStoich"] if isinstance(D["subunitStoich"], list) else [
				D["subunitStoich"]]
			for subunitId, subunitStoich in zip(subunitIds, subunitStoichs):

				if subunitId in mol_names:
					rowIdx = mol_names.index(subunitId)
					stoichMatrixMonomersI.append(rowIdx)
					stoichMatrixMonomersJ.append(colIdx)
					stoichMatrixMonomersV.append(-1. * subunitStoich)

		stoichMatrixMonomersI = np.array(stoichMatrixMonomersI)
		stoichMatrixMonomersJ = np.array(stoichMatrixMonomersJ)
		stoichMatrixMonomersV = np.array(stoichMatrixMonomersV)
		shape = (stoichMatrixMonomersI.max() + 1, stoichMatrixMonomersJ.max() + 1)

		out = np.zeros(shape, np.float64)
		out[stoichMatrixMonomersI, stoichMatrixMonomersJ] = stoichMatrixMonomersV
		return out



