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

		# Build dictionary mapping complexes to their base monomer subunit
		# composition from the modified_proteins table in the flat file (since
		# its not possible to get the base monomers from the raw TCS tsv files):
		self.complex_to_monomer = self._buildComplexToMonomer(
			sim_data, raw_data.modified_proteins, self.molecule_names)

		# Build list of molecules that include those from the original molecule
		# list and molecules from the modified_proteins.tsv table that are not
		# in the oiginal molecule list:
		self.modified_molecules = self.make_modified_molecule_list()

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

	def _buildComplexToMonomer(self, sim_data, modifiedFormsMonomers, tcsMolecules):
		'''
		Maps each complex to a dictionary that maps each subunit of the complex
		to its stoichiometry.

		This function also handles correcting compartment tags for subunits
		where the compartment tag listed in the modified_proteins.tsv table
		is inconsistent with the molecule's tag saved in the bulk molecule data.
		'''
		D = {}
		for row in modifiedFormsMonomers:
			# tags on the molecule compartment found in tcsMolecules
			molecule_and_location = f"{row['id']}[{row['compartment']}]"
			# filter by molecules actively being used in the simulation:
			if molecule_and_location in tcsMolecules:
				D[molecule_and_location] = {}
				for subunit in row["subunits"]:
					# We only care about mapping to protein monomers for now
					# and PI[c] stoichiometry is off for some complexes so we
					# can skip it for now (see #975)
					if subunit['monomer'] == 'Pi[c]' or subunit['monomer'] == 'PI[c]':
						continue
					# Since the compartment tag for some subunits may be off,
					# check if the subunit tag needs to be corrected before
					# adding it to the dictionary:
					if sim_data.getter.is_valid_molecule(subunit["monomer"]):
						D[molecule_and_location][str(subunit["monomer"])] = float(
							subunit["stoichiometry"])
					else:
						# Correct the compartment tag:
						subunit_wo_tag = subunit["monomer"].split("[")[0]
						compartment_tag = sim_data.getter.get_compartment(subunit_wo_tag)
						new_subunit_ID = f"{subunit_wo_tag}[{compartment_tag[0]}]"
						D[molecule_and_location][str(new_subunit_ID)] = float(
							subunit["stoichiometry"])

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
		Builds stoichiometry matrix for complexes to their base monomers subunits
		Rows: modified molecules (complexes and monomers, including base monomer
		subunits that were not included in the original molecule list extracted
		from two_component_systems.tsv, but were included in the modified_proteins.tsv table)
		Columns: complexes
		Values: base monomer stoichiometry (not including non-monomer subunits,
		like metabolites, ATP, etc.)
		'''
		ids_complexes = self.complex_to_monomer.keys()
		molecule_names = list(self.modified_molecules)
		stoichMatrixMonomersI = []
		stoichMatrixMonomersJ = []
		stoichMatrixMonomersV = []
		for colIdx, id_complex in enumerate(ids_complexes):
			D = self.get_monomers(id_complex)
			rowIdx = molecule_names.index(id_complex)
			stoichMatrixMonomersI.append(rowIdx)
			stoichMatrixMonomersJ.append(colIdx)
			stoichMatrixMonomersV.append(1.)

			for subunitId, subunitStoich in zip(D["subunitIds"], D["subunitStoich"]):
				if subunitId in molecule_names:
					rowIdx = molecule_names.index(subunitId)
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
		of one.
		'''
		info = self.complex_to_monomer
		if cplxId in info:
			out = {
				'subunitIds': list(info[cplxId].keys()),
				'subunitStoich': list(info[cplxId].values())}
		else:
			# Return stoich of 1 for monomers passed through:
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

	def make_modified_molecule_list(self):
		"""
		Since the raw modified_proteins table contains proteins that are either

		1. not in the original molecule pool provided in two_component_systems.tsv
		(PHOQ-MONOMER, PHOR-MONOMER, ARCB-MONOMER, and NARX-MONOMER), or,
		2. in the molecule pool but have different compartment tags from what
		the simulation expects (DUCS-MONOMER[i] should be DCUS-MONOMER[c]),

		the list of molecules that the stoich_matrix_monomers() function pulls
		from needs to be manually expanded to include and correct for these cases.

		This function generates a modified molecule list that includes the original
		molecules from the TCS molecule pool (built in __init__()) and the
		new proteins from the modified_proteins.tsv file that were not included
		in the original list (or were, but needed corrected compartment tags).

		# TODO (mia): eventually go through and fix the compartment tag for
		DCUS-MONOMER from [c] to [i] to avoid needing to check compartment
		validity and to be more consistent with true biology.
		"""
		# Orginal monomers (extracted from two_component_systems.tsv):
		tcs_molecules = self.molecule_names.tolist()
		# Obtain the new molecules from the complex_to_monomer dictionary
		# (generated from modified_proteins.tsv):
		new_molecules = []
		for complex in self.complex_to_monomer.keys():
			# Index to the individual subunit dictionaries for the complex:
			subunits = self.get_monomers(complex)
			for subunit in subunits['subunitIds']:
				# Append subunits not in the original list to the new list:
				if subunit not in tcs_molecules and subunit not in new_molecules:
					new_molecules.append(subunit)
		return np.array(tcs_molecules + new_molecules)

	def _view_matrix_with_row_and_col_names(self, rows, cols, matrix):
		"""
		Returns stoichiometry matrix as DataFrame with row and column labels.

		NOTE:
			for self.stoich_matrix(): rows=self.molecule_names, cols=self.ids_reactions
			for self.stoich_matrix_monomers(): rows=self.modified_molecules, cols=self.complex_to_monomer.keys()
		"""
		import pandas as pd
		return pd.DataFrame(matrix, index=rows, columns=cols)



