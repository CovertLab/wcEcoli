"""
SimulationData for metabolism process

TODO:
- improved estimate of ILE/LEU abundance or some external data point
- implement L1-norm minimization for AA concentrations
- find concentration for PI[c]
- add (d)NTP byproduct concentrations

@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 03/06/2015
"""

from __future__ import absolute_import, division, print_function

from copy import copy
import re
from typing import Any, Dict, List, Optional, Set

import numpy as np
import sympy as sp

from reconstruction.ecoli.knowledge_base_raw import KnowledgeBaseEcoli
from wholecell.utils import units

PPI_CONCENTRATION = 0.5e-3  # M, multiple sources
ILE_LEU_CONCENTRATION = 3.03e-4  # M, Bennett et al. 2009
ILE_FRACTION = 0.360  # the fraction of iso/leucine that is isoleucine; computed from our monomer data
ECOLI_PH = 7.2
MICROMOLAR_UNITS = units.umol / units.L
METABOLITE_CONCENTRATION_UNITS = units.mol / units.L

USE_ALL_CONSTRAINTS = False # False will remove defined constraints from objective

REVERSE_TAG = ' (reverse)'
REVERSE_REACTION_ID = '{{}}{}'.format(REVERSE_TAG)
ENZYME_REACTION_ID = '{}__{}'

# threshold (units.mmol / units.L) separates concentrations that are import constrained with
# max flux = 0 from unconstrained molecules.
IMPORT_CONSTRAINT_THRESHOLD =  1e-5

VERBOSE = False


class Metabolism(object):
	""" Metabolism """

	def __init__(self, raw_data, sim_data):
		# set solver and kinetic objective weight
		self.solver = "glpk-linear"
		if "linear" in self.solver:
			self.kinetic_objective_weight = sim_data.constants.metabolismKineticObjectiveWeightLinear
		else:
			self.kinetic_objective_weight = sim_data.constants.metabolismKineticObjectiveWeightQuadratic

		self.boundary = Boundary(raw_data, sim_data)

		# make a list of transport reactions
		transport_reactions_raw = raw_data.transport_reactions
		self.transport_reactions = []
		for transport_reaction in transport_reactions_raw:
			reaction_id = transport_reaction.get('reaction id')
			self.transport_reactions.append(reaction_id)

		self._buildBiomass(raw_data, sim_data)
		self._buildMetabolism(raw_data, sim_data)
		self._build_ppgpp_reactions(raw_data, sim_data)

	def _buildBiomass(self, raw_data, sim_data):
		wildtypeIDs = set(entry["molecule id"] for entry in raw_data.biomass)

		# Load the biomass function flat file as a dict
		self.biomassFunction = {entry['molecule id']:entry['coefficient'] for entry in raw_data.biomass}

		self.previousBiomassMeans = {entry['molecule id']:entry['mean flux'] for entry in raw_data.previousBiomassFluxes}
		self.previousBiomassLog10Means = {entry['molecule id']:entry['mean log10 flux'] for entry in raw_data.previousBiomassFluxes}
		self.previousBiomassStds = {entry['molecule id']:entry['standard deviation'] for entry in raw_data.previousBiomassFluxes}

		# Create vector of metabolite target concentrations

		# Since the data only covers certain metabolites, we need to rationally
		# expand the dataset to include the other molecules in the biomass
		# function.

		# First, load in metabolites that do have concentrations, then assign
		# compartments according to those given in the biomass objective.  Or,
		# if there is no compartment, assign it to the cytoplasm.

		metaboliteIDs = []
		metaboliteConcentrations = []
		metaboliteConcentrationData = dict(
			(m["Metabolite"], m["Concentration"].asNumber(METABOLITE_CONCENTRATION_UNITS))
			for m in raw_data.metaboliteConcentrations)

		wildtypeIDtoCompartment = {
			wildtypeID[:-3] : wildtypeID[-3:]
			for wildtypeID in wildtypeIDs
			} # this assumes biomass reaction components only exist in a single compartment

		for metaboliteID, concentration in metaboliteConcentrationData.viewitems():
			if metaboliteID in wildtypeIDtoCompartment:
				metaboliteIDs.append(
					metaboliteID + wildtypeIDtoCompartment[metaboliteID]
					)

			else:
				metaboliteIDs.append(
					metaboliteID + "[c]"
					)

			metaboliteConcentrations.append(concentration)

		# ILE/LEU: split reported concentration according to their relative abundances
		ileRelative = ILE_FRACTION
		leuRelative = 1 - ileRelative

		metaboliteIDs.append("ILE[c]")
		metaboliteConcentrations.append(ileRelative * ILE_LEU_CONCENTRATION)

		metaboliteIDs.append("LEU[c]")
		metaboliteConcentrations.append(leuRelative * ILE_LEU_CONCENTRATION)

		# CYS/SEC/GLY: concentration based on other amino acids
		aaConcentrations = []

		for aaIndex, aaID in enumerate(sim_data.amino_acid_1_to_3_ordered.values()):
			if aaID in metaboliteIDs:
				metIndex = metaboliteIDs.index(aaID)
				aaConcentrations.append(metaboliteConcentrations[metIndex])

		aaSmallestConc = min(aaConcentrations)

		metaboliteIDs.append("GLY[c]")
		metaboliteConcentrations.append(
			metaboliteConcentrations[metaboliteIDs.index("L-ALPHA-ALANINE[c]")]
			)

		metaboliteIDs.append("CYS[c]")
		metaboliteConcentrations.append(aaSmallestConc)

		metaboliteIDs.append("L-SELENOCYSTEINE[c]")
		metaboliteConcentrations.append(aaSmallestConc)

		# DGTP: set to smallest of all other DNTP concentrations
		dntpConcentrations = []

		for dntpIndex, dntpID in enumerate(sim_data.moleculeGroups.dNtpIds):
			if dntpID in metaboliteIDs:
				metIndex = metaboliteIDs.index(dntpID)
				dntpConcentrations.append(metaboliteConcentrations[metIndex])

		dntpSmallestConc = min(dntpConcentrations)

		metaboliteIDs.append("DGTP[c]")
		metaboliteConcentrations.append(dntpSmallestConc)

		# H: from reported pH
		hydrogenConcentration = 10**(-ECOLI_PH)

		metaboliteIDs.append("PROTON[c]")
		metaboliteConcentrations.append(hydrogenConcentration)

		# PPI: multiple sources report 0.5 mM
		metaboliteIDs.append("PPI[c]")
		metaboliteConcentrations.append(PPI_CONCENTRATION)

		metaboliteIDs.append("PI[c]")
		metaboliteConcentrations.append(PPI_CONCENTRATION)

		# Add byproducts with no annotated concentration to force recycling
		metaboliteIDs.append("UMP[c]")
		metaboliteConcentrations.append(2.40e-5)

		# include metabolites that are part of biomass
		for key, value in sim_data.mass.getBiomassAsConcentrations(sim_data.doubling_time).iteritems():
			metaboliteIDs.append(key)
			metaboliteConcentrations.append(value.asNumber(METABOLITE_CONCENTRATION_UNITS))

		# save concentrations as class variables
		self.concentrationUpdates = ConcentrationUpdates(dict(zip(
			metaboliteIDs,
			METABOLITE_CONCENTRATION_UNITS * np.array(metaboliteConcentrations)
			)),
			raw_data.equilibriumReactions,
			self.boundary.exchange_data_dict,
		)
		self.concDict = self.concentrationUpdates.concentrationsBasedOnNutrients("minimal")
		self.nutrientsToInternalConc = {}
		self.nutrientsToInternalConc["minimal"] = self.concDict.copy()

	def _buildMetabolism(self, raw_data, sim_data):
		"""
		Build the matrices/vectors for metabolism (FBA)
		Reads in and stores reaction and kinetic constraint information
		"""

		(reactionStoich, reversibleReactions, catalysts
			) = self.extract_reactions(raw_data, sim_data)

		# Load kinetic reaction constraints from raw_data
		constraints = self.extract_kinetic_constraints(raw_data, sim_data,
			stoich=reactionStoich, catalysts=catalysts)

		# Make modifications from kinetics data
		constraints, reactionStoich, catalysts, reversibleReactions = self._replace_enzyme_reactions(
			constraints, reactionStoich, catalysts, reversibleReactions)
		# TODO: check saturation with concentrations
		# TODO: lambdify saturation


		# TODO (cleanup): remove start
		# Initialize variables to store kinetic constraint information
		constraintDict = {}
		constraintIdList = []
		enzymeIdList = []
		kineticsSubstratesList = []
		reactionsToConstraintsDict = {}

		for constraint in raw_data.enzymeKinetics:
			if constraint["rateEquationType"] == "custom":
				continue

			constraintId = constraint["reactionID"].encode("utf-8") + "__" + constraint["enzymeIDs"].encode("utf-8") + "__%f" % (constraint["kcat"].asNumber(1 / units.s))
			assert len(constraint["Concentration Substrates"]) == len(constraint["kM"]) + len(constraint["kI"]), "Concentration Substrates are wrong length"
			assert constraintId not in constraintIdList, "constraintId already exists"

			# Get compartment for enzyme
			enzymeId = (constraint["enzymeIDs"] + "[" + sim_data.getter.getLocation([constraint["enzymeIDs"]])[0][0] + "]").encode("utf-8")
			constraint["enzymeIDs"] = enzymeId
			assert enzymeId in catalysts[constraint["reactionID"]], "%s is not a catalyst for %s according to FBA reconstruction" % (enzymeId, constraint["reactionID"])

			# Get compartments for Concentration Substrates
			concentrationSubstrates = []
			for substrate in constraint["Concentration Substrates"]:
				# In current implementation, anything with a concentration exists in the cytosol
				substrateWithCompartment = substrate.encode("utf-8") + "[c]"
				if substrateWithCompartment not in self.concDict:
					raise Exception, "Don't have concentration for %s" % substrateWithCompartment
				concentrationSubstrates.append(substrateWithCompartment)
				kineticsSubstratesList.append(substrateWithCompartment)
			constraint["Concentration Substrates"] = concentrationSubstrates

			# Get compartments for substrates
			substrates = []
			for substrate in constraint["substrateIDs"]:
				# In current implementation, anything with a concentration exists in the cytosol
				substrateFound = False
				for rxnSubstrate in reactionStoich[constraint["reactionID"]]:
					if rxnSubstrate.startswith(substrate + "["):
						substrateFound = True
						substrates.append(rxnSubstrate.encode("utf-8"))
				if not substrateFound:
					raise Exception, "Could not find compartment for substrate %s" % substrate
			constraint["substrateIDs"] = substrates

			# Adjust kcat for temperature
			temperature = constraint["Temp"]

			# If temperature not reported, assume 25 C
			if type(temperature) == str:
				temperature = 25
			constraint["kcatAdjusted"] = 2**((37. - temperature) / 10.) * constraint["kcat"]

			# Fix reactionID based on directionality
			stoichiometry = reactionStoich[constraint["reactionID"]]
			stoichVals = []
			for substrate in constraint["substrateIDs"]:
				stoichVals.append(stoichiometry[substrate])
			stoichVals = np.array(stoichVals)

			if np.all(stoichVals < 0):
				forward = True
			elif np.all(stoichVals > 0):
				forward = False
			else:
				raise Exception, "Have data for some reactants and some products (this is an inconsistency)"

			constraint["reactionID"] = constraint["reactionID"].encode("utf-8")
			if forward == False:
				constraint["reactionID"] = REVERSE_REACTION_ID.format(constraint["reactionID"])

			# Get rid of constraints for reverse reactions that the FBA reconstruction says should not exist
			# (i.e., if the FBA reconstruction says the reaction is irreversible but we have a constraint on
			#  the reverse reaction, drop the constraint)
			if constraint["reactionID"] not in reactionStoich:
				continue

			enzymeIdList.append(enzymeId)
			constraintIdList.append(constraintId)
			constraintDict[constraintId] = constraint
			if constraint["reactionID"] not in reactionsToConstraintsDict:
				reactionsToConstraintsDict[constraint["reactionID"]] = []
			reactionsToConstraintsDict[constraint["reactionID"]].append(constraintId)

		constraintIdList = sorted(constraintIdList)
		constrainedReactionList = sorted(reactionsToConstraintsDict)
		kineticsSubstratesList = sorted(set(kineticsSubstratesList))
		enzymeIdList = sorted(set(enzymeIdList))

		# split out reactions that are kinetically constrained and that have more than one enzyme that catalyzes the reaction
		for rxn in constrainedReactionList:
			rxn_catalysts = catalysts[rxn]
			if len(rxn_catalysts) > 1:
				for catalyst in rxn_catalysts:
					# create new reaction name with enzyme appended to the end
					if rxn.endswith(REVERSE_TAG):
						newReaction = REVERSE_REACTION_ID.format("%s__%s" % (rxn[:-len(REVERSE_TAG)], catalyst[:-3]))
					else:
						newReaction = "%s__%s" % (rxn, catalyst[:-3])

					# add the new reaction to appropriate lists and dicts
					if rxn in reversibleReactions:
						reversibleReactions.append(newReaction)

					reactionStoich[newReaction] = copy(reactionStoich[rxn])
					catalysts[newReaction] = [catalyst]
					for constraint in reactionsToConstraintsDict[rxn]:
						if constraintDict[constraint]["enzymeIDs"] == catalyst:
							constraintDict[constraint]["reactionID"] = newReaction
							if newReaction not in reactionsToConstraintsDict:
								reactionsToConstraintsDict[newReaction] = []
							reactionsToConstraintsDict[newReaction].append(constraint)

				# remove old reaction name
				reactionStoich.pop(rxn)
				catalysts.pop(rxn)
				reactionsToConstraintsDict.pop(rxn)
				if rxn in reversibleReactions:
					reversibleReactions.pop(reversibleReactions.index(rxn))

		constrainedReactionList = sorted(reactionsToConstraintsDict)
		catalystsList = sorted({c for rxn_catalysts in catalysts.values() for c in rxn_catalysts})
		reactionCatalystsList = sorted(catalysts)

		# TODO (cleanup): remove end

		# Create catalysis matrix (to be used in the simulation)
		catalysisMatrixI = []
		catalysisMatrixJ = []
		catalysisMatrixV = []

		for row, reaction in enumerate(reactionCatalystsList):
			for catalyst in catalysts[reaction]:
				col = catalystsList.index(catalyst)
				catalysisMatrixI.append(row)
				catalysisMatrixJ.append(col)
				catalysisMatrixV.append(1)

		catalysisMatrixI = np.array(catalysisMatrixI)
		catalysisMatrixJ = np.array(catalysisMatrixJ)
		catalysisMatrixV = np.array(catalysisMatrixV)

#		shape = (catalysisMatrixI.max() + 1, catalysisMatrixJ.max() + 1)
#		catalysisMatrix = np.zeros(shape, np.float64)
#		catalysisMatrix[catalysisMatrixI, catalysisMatrixJ] = catalysisMatrixV

		# Create constraint to reaction matrix (to be used in the simulation)
		constraintToReactionMatrixI = []
		constraintToReactionMatrixJ = []
		constraintToReactionMatrixV = []

		for row, reaction in enumerate(constrainedReactionList):
			for constraintId in reactionsToConstraintsDict[reaction]:
				col = constraintIdList.index(constraintId)
				constraintToReactionMatrixI.append(row)
				constraintToReactionMatrixJ.append(col)
				constraintToReactionMatrixV.append(1)

		constraintToReactionMatrixI = np.array(constraintToReactionMatrixI)
		constraintToReactionMatrixJ = np.array(constraintToReactionMatrixJ)
		constraintToReactionMatrixV = np.array(constraintToReactionMatrixV)

#		shape = (constraintToReactionMatrixI.max() + 1, constraintToReactionMatrixJ.max() + 1)
#		constraintToReactionMatrix = np.zeros(shape, np.float64)
#		constraintToReactionMatrix[constraintToReactionMatrixI, constraintToReactionMatrixJ] = constraintToReactionMatrixV

		# Use Sympy to create vector function that returns all kinetic constraints
		kineticsSubstrates = sp.symbols(["kineticsSubstrates[%d]" % idx for idx in xrange(len(kineticsSubstratesList))])
		enzymes = sp.symbols(["enzymes[%d]" % idx for idx in xrange(len(enzymeIdList))])
		constraints = [sp.symbol.S.Zero] * len(constraintIdList)
		constraintIsKcatOnly = np.zeros(len(constraintIdList))

		for constraintIdx, constraintId in enumerate(constraintIdList):
			constraint = constraintDict[constraintId]
			if constraint["rateEquationType"] == "custom":
				raise Exception

			enzymeIdx = enzymeIdList.index(constraint["enzymeIDs"])
			constraints[constraintIdx] = constraint["kcatAdjusted"].asNumber(1 / units.s) * enzymes[enzymeIdx]
			if len(constraint["kM"]) == 0 and len(constraint["kI"]) == 0:
				constraintIsKcatOnly[constraintIdx] = 1

			concSubstratePos = 0
			for kM in constraint["kM"]:
				kineticsSubstrateIdx = kineticsSubstratesList.index(constraint["Concentration Substrates"][concSubstratePos])
				S = kineticsSubstrates[kineticsSubstrateIdx]
				constraints[constraintIdx] *= (S / (S + kM))
				concSubstratePos += 1

			for kI in constraint["kI"]:
				kineticsSubstrateIdx = kineticsSubstratesList.index(constraint["Concentration Substrates"][concSubstratePos])
				I = kineticsSubstrates[kineticsSubstrateIdx]
				constraints[constraintIdx] *= (kI / (kI + I))
				concSubstratePos += 1

		self._kineticConstraints = str(constraints)
		self._compiledConstraints = None

		# Properties for FBA reconstruction
		self.reactionStoich = reactionStoich
		self.maintenanceReaction = {"ATP[c]": -1, "WATER[c]": -1, "ADP[c]": +1, "PI[c]": +1, "PROTON[c]": +1,}
		self.reversibleReactions = reversibleReactions

		# Properties for catalysis matrix (to set hard bounds)
		self.reactionCatalysts = catalysts
		self.catalystsList = catalystsList
		self.reactionCatalystsList = reactionCatalystsList
		self.catalysisMatrixI = catalysisMatrixI
		self.catalysisMatrixJ = catalysisMatrixJ
		self.catalysisMatrixV = catalysisMatrixV

		# Properties for setting flux targets
		self.constraintIdList = constraintIdList
		self.constrainedReactionList = constrainedReactionList
		self.constraintToReactionMatrixI = constraintToReactionMatrixI
		self.constraintToReactionMatrixJ = constraintToReactionMatrixJ
		self.constraintToReactionMatrixV = constraintToReactionMatrixV
		self.enzymeIdList = enzymeIdList
		self.kineticsSubstratesList = kineticsSubstratesList
		self.constraintDict = constraintDict
		self.reactionsToConstraintsDict = reactionsToConstraintsDict
		self.constraintIsKcatOnly = constraintIsKcatOnly
		self.useAllConstraints = USE_ALL_CONSTRAINTS
		self.constraintsToDisable = [rxn["disabled reaction"] for rxn in raw_data.disabledKineticReactions]

	def _build_ppgpp_reactions(self, raw_data, sim_data):
		'''
		Creates structures for ppGpp reactions for use in polypeptide_elongation.

		Attributes set:
			ppgpp_synthesis_reaction (str): reaction ID for ppGpp synthesis
				(catalyzed by RelA and SpoT)
			ppgpp_degradation_reaction (str): reaction ID for ppGpp degradation
				(catalyzed by SpoT)
			ppgpp_reaction_names (list[str]): names of reaction involved in ppGpp
			ppgpp_reaction_metabolites (list[str]): names of metabolites in
				ppGpp reactions
			ppgpp_reaction_stoich (array[int]): 2D array with metabolites on rows
				and reactions on columns containing the stoichiometric coefficient
		'''

		self.ppgpp_synthesis_reaction = 'GDPPYPHOSKIN-RXN'
		self.ppgpp_degradation_reaction = 'PPGPPSYN-RXN'

		self.ppgpp_reaction_names = [
			self.ppgpp_synthesis_reaction,
			self.ppgpp_degradation_reaction,
			]

		self.ppgpp_reaction_metabolites = []

		# Indices (i: metabolite, j: reaction) and values (v: stoichiometry)
		# for sparse reaction matrix
		metabolite_indices = {}
		new_index = 0
		rxn_i = []
		rxn_j = []
		rxn_v = []

		# Record sparse indices in the matrix
		for j, rxn in enumerate(self.ppgpp_reaction_names):
			for met, stoich in self.reactionStoich[rxn].items():
				idx = metabolite_indices.get(met, new_index)

				if idx == new_index:
					metabolite_indices[met] = new_index
					self.ppgpp_reaction_metabolites.append(met)
					new_index += 1

				rxn_i.append(idx)
				rxn_j.append(j)
				rxn_v.append(stoich)

		# Assemble matrix based on indices
		# new_index is number of metabolites, j+1 is number of reactions
		self.ppgpp_reaction_stoich = np.zeros((new_index, j+1), dtype=np.int32)
		self.ppgpp_reaction_stoich[rxn_i, rxn_j] = rxn_v

	def getKineticConstraints(self, enzymes, substrates):
		'''
		Allows for dynamic code generation for kinetic constraint calculation
		for use in Metabolism process. Inputs should be unitless but the order
		of magnitude should match the kinetics parameters (umol/L/s).

		If trying to pickle sim_data object after function has been called,
		_compiledConstraints might not be able to be pickled.  See
		__getstate__(), __setstate__() comments on PR 111 to address.

		Returns np.array of floats of the kinetic constraint target for each
		reaction with kinetic parameters
		Inputs:
			enzymes (np.array of floats) - concentrations of enzymes associated
				with kinetics constraints
			substrates (np.array of floats) - concentrations of substrates
				associated with kinetics constraints
		'''

		if self._compiledConstraints is None:
			self._compiledConstraints = eval('lambda enzymes, kineticsSubstrates: np.array(%s)\n'
				% self._kineticConstraints, {'np': np}, {}
				)

		return self._compiledConstraints(enzymes, substrates)

	def exchangeConstraints(self, exchangeIDs, coefficient, targetUnits, currentNutrients, exchange_data, concModificationsBasedOnCondition = None):
		"""
		Called during Metabolism process
		Returns the homeostatic objective concentrations based on the current nutrients
		Returns levels for external molecules available to exchange based on the current nutrients
		"""

		unconstrained_exchange_molecules = exchange_data["importUnconstrainedExchangeMolecules"]
		constrained_exchange_molecules = exchange_data["importConstrainedExchangeMolecules"]

		concDict = self.concentrationUpdates.concentrationsBasedOnNutrients(currentNutrients)
		if concModificationsBasedOnCondition is not None:
			concDict.update(concModificationsBasedOnCondition)

		conversion = targetUnits.asUnit(self.concentrationUpdates.units)
		newObjective = dict((key, (val / conversion).asNumber()) for key, val in concDict.iteritems())

		externalMoleculeLevels = np.zeros(len(exchangeIDs), np.float64)

		for index, moleculeID in enumerate(exchangeIDs):
			if moleculeID in unconstrained_exchange_molecules:
				externalMoleculeLevels[index] = np.inf
			elif moleculeID in constrained_exchange_molecules:
				externalMoleculeLevels[index] = (
					constrained_exchange_molecules[moleculeID] * coefficient
					).asNumber(targetUnits)
			else:
				externalMoleculeLevels[index] = 0.

		return externalMoleculeLevels, newObjective

	def set_supply_constants(self, sim_data):
		"""
		Sets constants to determine amino acid supply during translation.

		Args:
			sim_data (SimulationData object)

		Sets class attributes:
			KI_aa_synthesis (ndarray[float]): KI for each AA for synthesis
				portion of supply (in units of METABOLITE_CONCENTRATION_UNITS)
			KM_aa_export (ndarray[float]): KM for each AA for export portion
				of supply (in units of METABOLITE_CONCENTRATION_UNITS)
			fraction_supply_rate (float): fraction of AA supply that comes from
				a base synthesis rate
			fraction_import_rate (ndarray[float]): fraction of AA supply that
				comes from AA import if nutrients are present
			base_aa_conc (ndarray[float]): expected AA conc in basal condition
				(in units of METABOLITE_CONCENTRATION_UNITS)

		Assumptions:
			- Each internal amino acid concentration in 'minimal_plus_amino_acids'
			media is not lower than in 'minimal' media

		TODO (Travis):
			Base on measured KI and KM values.
			Add impact from synthesis enzymes and transporters.
			Better handling of concentration assumption
		"""

		aa_ids = sim_data.moleculeGroups.aaIDs
		conc = self.concentrationUpdates.concentrationsBasedOnNutrients

		aa_conc_basal = np.array([
			conc('minimal')[aa].asNumber(METABOLITE_CONCENTRATION_UNITS)
			for aa in aa_ids])
		aa_conc_aa_media = np.array([
			conc('minimal_plus_amino_acids')[aa].asNumber(METABOLITE_CONCENTRATION_UNITS)
			for aa in aa_ids])

		# Lower concentrations might produce strange rates (excess supply or
		# negative import when present externally) and constants so raise
		# to double check the implementation
		if not np.all(aa_conc_basal <= aa_conc_aa_media):
			aas = np.array(aa_ids)[np.where(aa_conc_basal > aa_conc_aa_media)]
			raise ValueError('Check that amino acid concentrations should be lower in amino acid media for {}'.format(aas))

		f_inhibited = sim_data.constants.fraction_supply_inhibited
		f_exported = sim_data.constants.fraction_supply_exported

		# Assumed units of METABOLITE_CONCENTRATION_UNITS for KI and KM
		self.KI_aa_synthesis = f_inhibited * aa_conc_basal / (1 - f_inhibited)
		self.KM_aa_export = (1 / f_exported - 1) * aa_conc_aa_media
		self.fraction_supply_rate = 1 - f_inhibited + aa_conc_basal / (self.KM_aa_export + aa_conc_basal)
		self.fraction_import_rate = 1 - (self.fraction_supply_rate + 1 / (1 + aa_conc_aa_media / self.KI_aa_synthesis) - f_exported)

	def aa_supply_scaling(self, aa_conc, aa_present):
		"""
		Called during polypeptide_elongation process
		Determine amino acid supply rate scaling based on current amino acid
		concentrations.

		Args:
			aa_conc (ndarray[float] with mol / volume units): internal
				concentration for each amino acid
			aa_present (ndarray[bool]): whether each amino acid is in the
				external environment or not

		Returns:
			ndarray[float]: scaling for the supply of each amino acid with
				higher supply rate if >1, lower supply rate if <1
		"""

		aa_conc = aa_conc.asNumber(METABOLITE_CONCENTRATION_UNITS)

		aa_supply = self.fraction_supply_rate
		aa_import = aa_present * self.fraction_import_rate
		aa_synthesis = 1 / (1 + aa_conc / self.KI_aa_synthesis)
		aa_export = aa_conc / (self.KM_aa_export + aa_conc)
		supply_scaling = aa_supply + aa_import + aa_synthesis - aa_export

		return supply_scaling

	@staticmethod
	def extract_reactions(raw_data, sim_data):
		# type: (KnowledgeBaseEcoli, SimulationDataEcoli) -> (Dict[str, Dict[str, int]], List[str], Dict[str, List[str]])
		"""
		Args:
			raw_data: knowledge base data
			sim_data: simulation data

		Returns:
			reactionStoich: {reaction ID: {metabolite ID with location tag: stoichiometry}}
				stoichiometry of metabolites for each reaction
			reversibleReactions: reaction IDs for reactions that have a reverse
				complement, does not have reverse tag
			reactionCatalysts: {reaction ID: enzyme IDs with location tag}
				enzyme catalysts for each reaction with known catalysts, likely
				a subset of reactions in stoich
		"""

		# Initialize variables to store reaction information
		reactionStoich = {}
		reversibleReactions = []
		reactionCatalysts = {}

		# Load and parse reaction information from raw_data
		for reaction in raw_data.reactions:
			reactionID = reaction["reaction id"]
			stoich = reaction["stoichiometry"]
			reversible = reaction["is reversible"]

			if len(stoich) <= 1:
				raise Exception("Invalid biochemical reaction: {}, {}".format(reactionID, stoich))

			reactionStoich[reactionID] = stoich

			catalystsForThisRxn = []
			for catalyst in reaction["catalyzed by"]:
				try:
					catalystWithLoc = (catalyst + "[" + sim_data.getter.getLocation([catalyst])[0][0] + "]").encode("utf-8")
					catalystsForThisRxn.append(catalystWithLoc)
				# If we don't have the catalyst in our reconstruction, drop it
				except KeyError:
					if VERBOSE:
						print('Skipping catalyst {} for {} since it is not in the model'
							.format(catalyst, reactionID))

			if len(catalystsForThisRxn) > 0:
				reactionCatalysts[reactionID] = catalystsForThisRxn

			# Add the reverse reaction
			if reversible:
				reverseReactionID = REVERSE_REACTION_ID.format(reactionID)
				reactionStoich[reverseReactionID] = {
					moleculeID:-stoichCoeff
					for moleculeID, stoichCoeff in reactionStoich[reactionID].viewitems()
					}

				reversibleReactions.append(reactionID)
				if len(catalystsForThisRxn) > 0:
					reactionCatalysts[reverseReactionID] = list(reactionCatalysts[reactionID])

		return reactionStoich, reversibleReactions, reactionCatalysts

	@staticmethod
	def match_reaction(stoich, catalysts, rxn, enz, mets, direction=None):
		# type: (Dict[str, Dict[str, int]], Dict[str, List[str]], str, str, List[str], Optional[str]) -> Optional[str]
		"""
		Matches a given reaction (rxn) to reactions that exist in stoich given
		that enz is known to catalyze the reaction and mets are reactants in
		the reaction.  Can perform a fuzzy reaction match since rxn just needs
		to be part of the actual reaction name to match specific instances of a
		reaction (eg. rxn="ALCOHOL-DEHYDROG-GENERIC-RXN" can match
		"ALCOHOL-DEHYDROG-GENERIC-RXN-ETOH/NAD//ACETALD/NADH/PROTON.30.").

		Args:
			stoich: {reaction ID: {metabolite ID with location tag: stoichiometry}}
				stoichiometry of metabolites for each reaction
			catalysts: {reaction ID: enzyme IDs with location tag}
				enzyme catalysts for each reaction with known catalysts,
				likely a subset of reactions in stoich
			rxn: reaction ID from kinetics to match to existing reactions
			enz: enzyme ID with location tag
			mets: metabolite IDs with no location tag from kinetics
			direction: reaction directionality, 'forward' or 'reverse' or None

		Returns:
			rxn: matched reaction ID to reaction in stoich with reverse tag
				if in the reverse direction, returns None if no match
		"""

		# Mapping to handle instances of metabolite classes in kinetics
		# Keys: specific molecules in kinetics file
		# Values: class of molecules in reactions file that contain the key
		class_mets = {
			'RED-THIOREDOXIN-MONOMER': 'Red-Thioredoxin',
			'RED-THIOREDOXIN2-MONOMER': 'Red-Thioredoxin',
			'RED-GLUTAREDOXIN': 'Red-Glutaredoxins',
			'GRXB-MONOMER': 'Red-Glutaredoxins',
			'GRXC-MONOMER': 'Red-Glutaredoxins',
			'OX-FLAVODOXIN1': 'Oxidized-flavodoxins',
			'OX-FLAVODOXIN2': 'Oxidized-flavodoxins',
		}

		# Match full reaction name from partial reaction in kinetics. Must
		# also match metabolites since there can be multiple reaction instances.
		match = False
		if rxn not in stoich:
			for long_rxn, long_mets in stoich.items():
				if rxn in long_rxn and not long_rxn.endswith(REVERSE_TAG):
					match = True
					stripped_enzs = {e[:-3] for e in catalysts.get(long_rxn, [])}
					stripped_mets = {m[:-3] for m in long_mets}
					if (np.all([class_mets.get(m, m) in stripped_mets for m in mets])
							and enz in stripped_enzs):
						# TODO: check if other reactions match instead of breaking on first
						rxn = long_rxn
						break
			else:
				if VERBOSE:
					if match:
						print('Partial reaction match: {} {} {} {} {}'.format(
							rxn, enz, stripped_enzs, mets, stripped_mets))
					else:
						print('No reaction match: {}'.format(rxn))
				return None

		# Determine direction of kinetic reaction from annotation or
		# metabolite stoichiometry.
		reverse_rxn = REVERSE_REACTION_ID.format(rxn)
		reverse_rxn_exists = reverse_rxn in stoich
		if direction:
			reverse = direction == 'reverse'
		else:
			s = {k[:-3]: v for k, v in stoich.get(rxn, {}).items()}
			direction = np.unique(np.sign([
				s.get(class_mets.get(m, m), 0) for m in mets]))
			if len(direction) == 0 and not reverse_rxn_exists:
				reverse = False
			elif len(direction) != 1 or direction[0] == 0:
				if VERBOSE:
					print('Conflicting directionality: {} {} {}'.format(
						rxn, mets, direction))
				return None
			else:
				reverse = direction[0] > 0

		# Verify a reverse reaction exists in the model
		if reverse:
			if reverse_rxn_exists:
				rxn = reverse_rxn
			else:
				if VERBOSE:
					print('No reverse reaction: {} {}'.format(rxn, mets))
				return None

		return rxn

	@staticmethod
	def temperature_adjusted_kcat(kcat, temp):
		# type: (units.Unum, float) -> np.ndarray[float]
		"""
		Args:
			kcat: enzyme turnover number(s) (1 / time)
			temp: temperature of measurement, defaults to 25 if ''

		Returns:
			temperature adjusted kcat values, in units of 1/s
		"""

		if temp == '':
			temp = 25
		return 2**((37. - temp) / 10.) * kcat.asNumber(1 / units.s)

	@staticmethod
	def _construct_default_saturation_equation(mets, kms, kis):
		# type: (List[str], List[float], List[float]) -> str
		"""
		Args:
			mets: metabolite IDs without location tag for KM and KI
				parameters ordered to match order of kms then kis
			kms: KM parameters associated with mets
			kis: KI parameters associated with mets

		Returns:
			saturation equation with metabolites to replace delimited
				by double quote (eg. "metabolite")
		"""

		# Check input dimensions
		n_params = len(kms) + len(kis)
		if n_params == 0:
			return '1'
		if n_params != len(mets):
			if VERBOSE:
				print('Saturation parameter mismatch: {} {} {}'.format(mets, kms, kis))
			return '1'

		terms = []
		# Add KM terms
		for m, k in zip(mets, kms):
			terms.append('1+{}/"{}"'.format(k, m))
		# Add KI terms
		for m, k in zip(mets[len(kms):], kis):
			terms.append('1+"{}"/{}'.format(m, k))

		# Enclose groupings if being multiplied together
		if len(terms) > 1:
			terms[0] = '(' + terms[0]
			terms[-1] += ')'

		return '1/({})'.format(')*('.join(terms))

	@staticmethod
	def _extract_custom_constraint(constraint):
		# type: (Dict[str, Any]) -> (Optional[np.ndarray[float]], List[str])
		"""
		Args:
			constraint: values defining a kinetic constraint with key:
				'customRateEquation' (str): mathematical representation of
					rate, must contain 'kcat*E'
				'customParameterVariables' (Dict[str, str]): mapping of
					variable names in the rate equation to metabolite IDs
					without location tags, must contain key 'E' (enzyme)
				'customParameterConstants' (List[str]): constant strings
					in the rate equation that correspond to values, must
					contain 'kcat'
				'customParameterConstantValues' (List[float]): values for
					each of the constant strings
				'Temp' (float or ''): temperature of measurement

		Returns:
			kcats: temperature adjusted kcat value, in units of 1/s
			saturation: saturation equation with metabolites to replace
				delimited by double quote (eg. "metabolite")
		"""

		equation = constraint['customRateEquation']
		variables = constraint['customParameterVariables']
		constant_keys = constraint['customParameterConstants']
		constant_values = constraint['customParameterConstantValues']
		temp = constraint['Temp']

		# Need to have these in the constraint
		kcat_str = 'kcat'
		enzyme_str = 'E'
		capacity_str = '{}*{}'.format(kcat_str, enzyme_str)

		# Make sure kcat exists
		if kcat_str not in constant_keys:
			if VERBOSE:
				print('Missing {} in custom constants: {}'.format(
					kcat_str, constant_keys))
			return None, []

		custom_kcat = 1 / units.s * np.array([constant_values[constant_keys.index(kcat_str)]])
		kcats = Metabolism.temperature_adjusted_kcat(custom_kcat, temp)

		# Make sure equation can be parsed, otherwise just return kcat
		if enzyme_str not in variables:
			if VERBOSE:
				print('Missing enzyme key ({}) in custom variables: {}'.format(
					enzyme_str, variables))
			return kcats, []
		if not capacity_str in equation:
			if VERBOSE:
				print('Expected to find {} in custom equation: {}'.format(
					capacity_str, equation))
			return kcats, []
		if len(constant_keys) != len(constant_values):
			if VERBOSE:
				print('Mismatch between constants: {} {}'.format(
					constant_keys, constant_values))
			return kcats, []

		# Substitute values into custom equations
		## Remove capacity to get only saturation
		new_equation = equation.replace(capacity_str, '1')

		## Tokenize equation to terms and symbols
		parsed_variables = re.findall('\w*', new_equation)[:-1]  # Remove trailing empty match
		parsed_symbols = re.findall('\W', new_equation)
		tokenized_equation = np.array(parsed_variables)
		tokenized_equation[tokenized_equation == ''] = parsed_symbols
		if ''.join(tokenized_equation) != new_equation:
			if VERBOSE:
				print('Error parsing custom equation: {}'.format(equation))
			return kcats, []

		## Replace terms with known constant values or sim molecule IDs
		custom_subs = {k: str(v) for k, v in zip(constant_keys, constant_values)}
		custom_subs.update({k: '"{}"'.format(v) for k, v in variables.items()})

		# Reconstruct saturation equation with replacements
		saturation = [''.join([custom_subs.get(token, token) for token in tokenized_equation])]

		return kcats, saturation

	@staticmethod
	def extract_kinetic_constraints(raw_data, sim_data, stoich=None, catalysts=None):
		# type: (KnowledgeBaseEcoli, SimulationDataEcoli, Optional[Dict[str, Dict[str, int]]], Optional[Dict[str, List[str]]]) -> Dict[(str, str), Dict[str, List[Any]]]
		"""
		Load and parse kinetic constraint information from raw_data

		Args:
			raw_data: knowledge base data
			sim_data: simulation data
			stoich: {reaction ID: {metabolite ID with location tag: stoichiometry}}
				stoichiometry of metabolites for each reaction, if None, data
				is loaded from raw_data and sim_data
			catalysts: {reaction ID: enzyme IDs with location tag}
				enzyme catalysts for each reaction with known catalysts, likely
				a subset of reactions in stoich, if None, data is loaded from
				raw_data and sim_data

		Returns:
			constraints: valid kinetic constraints for each reaction/enzyme pair
				{(reaction ID, enzyme with location tag): {
					'kcat': kcat values (List[float]),
					'saturation': saturation equations (List[str])
				}}

		TODO:
			use function in reflect script
		"""

		# Load data for optional args if needed
		if stoich is None or catalysts is None:
			loaded_stoich, _, loaded_catalysts = Metabolism.extract_reactions(raw_data, sim_data)

			if stoich is None:
				stoich = loaded_stoich
			if catalysts is None:
				catalysts = loaded_catalysts

		constraints = {}
		for constraint in raw_data.metabolism_kinetics:
			rxn = constraint['reactionID']
			enzyme = constraint['enzymeID']
			metabolites = constraint['substrateIDs']
			direction = constraint['direction']
			kms = list(constraint['kM'].asNumber(MICROMOLAR_UNITS))
			kis = list(constraint['kI'].asNumber(MICROMOLAR_UNITS))
			n_reactants = len(metabolites) - len(kis)
			matched_rxn = Metabolism.match_reaction(stoich, catalysts, rxn, enzyme,
				metabolites[:n_reactants], direction)
			if matched_rxn is None:
				continue

			# Ensure enzyme catalyzes reaction in model
			enzymes_tag_conversion = {e[:-3]: e for e in catalysts.get(matched_rxn, [])}
			if enzyme not in enzymes_tag_conversion:
				if VERBOSE:
					print('{} does not catalyze {}'.format(enzyme, matched_rxn))
				continue
			else:
				enzyme = enzymes_tag_conversion[enzyme]

			# Update metabolites with a location tag from the reaction
			# First look in reactants but some products can inhibit
			reactant_tags = {k[:-3]: k for k, v in stoich[matched_rxn].items() if v < 0}
			product_tags = {k[:-3]: k for k, v in stoich[matched_rxn].items() if v > 0}
			mets_with_tag = [
				reactant_tags.get(met, product_tags.get(met, None))
				for met in metabolites
				if met in reactant_tags or met in product_tags
			]
			if len(mets_with_tag) != len(metabolites):
				# Warn if verbose but no continue since we can still use kcat
				if VERBOSE:
					print('Could not match all metabolites: {} {}'.format(
						metabolites, mets_with_tag))

			# Extract kcat and saturation parameters
			if constraint['rateEquationType'] == 'custom':
				kcats, saturation = Metabolism._extract_custom_constraint(constraint)
				if kcats is None:
					continue
			else:
				kcats = Metabolism.temperature_adjusted_kcat(constraint['kcat'], constraint['Temp'])
				if len(kcats) > 1:
					if len(kcats) != len(kms) or len(kms) != len(mets_with_tag):
						if VERBOSE:
							print('Could not align kcats and kms: {} {} {} {}'.format(
								rxn, kcats, kms, mets_with_tag))
						continue

					saturation = [
						Metabolism._construct_default_saturation_equation([m], [km], [])
						for m, km in zip(mets_with_tag, kms)
					]
				else:
					saturation = [
						Metabolism._construct_default_saturation_equation(mets_with_tag, kms, kis)
					]

			# Add new kcats and saturation terms for the enzymatic reaction
			key = (matched_rxn, enzyme)
			entries = constraints.get(key, {})
			entries['kcat'] = entries.get('kcat', []) + list(kcats)
			entries['saturation'] = entries.get('saturation', []) + saturation
			constraints[key] = entries

		return constraints

	@staticmethod
	def _replace_enzyme_reactions(constraints, stoich, rxn_catalysts, reversible_rxns):
		# type: (Dict[(str, str), Dict[str, List[Any]]], Dict[str, Dict[str, int]], Dict[str, List[str]], List[str]) -> (Dict[str, Any], Dict[str, Dict[str, int]], Dict[str, List[str]], List[str])
		"""
		Modifies reaction IDs in data structures to duplicate reactions with
		kinetic constraints and multiple enzymes.

		Args:
			constraints: valid kinetic constraints for each reaction/enzyme pair
				{(reaction ID, enzyme with location tag): {
					'kcat': kcat values (List[float]),
					'saturation': saturation equations (List[str])
				}}
			stoich: {reaction ID: {metabolite ID with location tag: stoichiometry}}
				stoichiometry of metabolites for each reaction, if None, data
				is loaded from raw_data and sim_data
			rxn_catalysts: {reaction ID: enzyme IDs with location tag}
				enzyme catalysts for each reaction with known catalysts, likely
				a subset of reactions in stoich, if None, data is loaded from
				raw_data and sim_data
			reversible_rxns: reaction IDs for reactions that have a reverse
				complement, does not have reverse tag

		Returns:
			new_constraints: valid kinetic constraints for each reaction
				{reaction ID: {
					'enzyme': enzyme catalyst (str),
					'kcat': kcat values (List[float]),
					'saturation': saturation equations (List[str])
				}}
			stoich: {reaction ID: {metabolite ID with location tag: stoichiometry}}
				stoichiometry of metabolites for each reaction with updated
				reactions for enzyme catalyzed kinetic reactions
			rxn_catalysts: {reaction ID: enzyme IDs with location tag}
				enzyme catalysts for each reaction with known catalysts, likely
				a subset of reactions in stoich with updated
				reactions for enzyme catalyzed kinetic reactions
			reversible_rxns: reaction IDs for reactions that have a reverse
				complement with updated reactions for enzyme catalyzed kinetic
				reactions, does not have reverse tag
		"""

		new_constraints = {}

		n_catalysts = {rxn: len(catalysts) for rxn, catalysts in rxn_catalysts.items()}

		# Split out reactions that are kinetically constrained and that have
		# more than one enzyme that catalyzes the reaction
		for (rxn, enzyme), constraint in constraints.items():
			if n_catalysts[rxn] > 1:
				# Create new reaction name with enzyme appended to the end
				if rxn.endswith(REVERSE_TAG):
					new_rxn = REVERSE_REACTION_ID.format(ENZYME_REACTION_ID.format(rxn[:-len(REVERSE_TAG)], enzyme[:-3]))
				else:
					new_rxn = ENZYME_REACTION_ID.format(rxn, enzyme[:-3])

				# Add the new reaction to appropriate lists and dicts
				stoich[new_rxn] = copy(stoich[rxn])
				rxn_catalysts[new_rxn] = [enzyme]
				if rxn in reversible_rxns:
					reversible_rxns.append(new_rxn)

				# Remove enzyme from old reaction and remove old reaction if no
				# more enzyme catalysts
				rxn_catalysts[rxn].pop(rxn_catalysts[rxn].index(enzyme))
				if len(rxn_catalysts[rxn]) == 0:
					stoich.pop(rxn)
					rxn_catalysts.pop(rxn)
					if rxn in reversible_rxns:
						reversible_rxns.pop(reversible_rxns.index(rxn))
			else:
				new_rxn = rxn

			new_constraints[new_rxn] = dict(constraints[(rxn, enzyme)], enzyme=enzyme)

		return new_constraints, stoich, rxn_catalysts, reversible_rxns



# Class used to update metabolite concentrations based on the current nutrient conditions
class ConcentrationUpdates(object):
	def __init__(self, concDict, equilibriumReactions, exchange_data_dict):
		self.units = units.getUnit(concDict.values()[0])
		self.defaultConcentrationsDict = dict((key, concDict[key].asNumber(self.units)) for key in concDict)
		self.exchange_fluxes = self._exchange_flux_present(exchange_data_dict)

		# factor of internal amino acid increase if amino acids present in nutrients
		self.moleculeScaleFactors = {
			"L-ALPHA-ALANINE[c]": 2.,
			"ARG[c]": 2.,
			"ASN[c]": 2.,
			"L-ASPARTATE[c]": 2.,
			"CYS[c]": 2.,
			"GLT[c]": 1.1,
			"GLN[c]": 2.,
			"GLY[c]": 2.,
			"HIS[c]": 2.,
			"ILE[c]": 2.,
			"LEU[c]": 2.,
			"LYS[c]": 2.,
			"MET[c]": 2.,
			"PHE[c]": 2.,
			"PRO[c]": 2.,
			"SER[c]": 2.,
			"THR[c]": 2.,
			"TRP[c]": 2.,
			"TYR[c]": 2.,
			"L-SELENOCYSTEINE[c]": 2.,
			"VAL[c]": 2.,
		}

		self.moleculeSetAmounts = self._addMoleculeAmounts(equilibriumReactions, self.defaultConcentrationsDict)

	# return adjustments to concDict based on nutrient conditions
	def concentrationsBasedOnNutrients(self, media_id=None):
		concentrationsDict = self.defaultConcentrationsDict.copy()

		metaboliteTargetIds = sorted(concentrationsDict.keys())
		concentrations = self.units * np.array([concentrationsDict[k] for k in metaboliteTargetIds])
		concDict = dict(zip(metaboliteTargetIds, concentrations))

		if media_id is not None:
			exchanges = self.exchange_fluxes[media_id]
			for moleculeName, setAmount in self.moleculeSetAmounts.iteritems():
				if ((moleculeName in exchanges and (moleculeName[:-3] + "[c]" not in self.moleculeScaleFactors or moleculeName == "L-SELENOCYSTEINE[c]"))
						or (moleculeName in self.moleculeScaleFactors and moleculeName[:-3] + "[p]" in exchanges)):
					concDict[moleculeName] = setAmount

		return concDict

	def _exchange_flux_present(self, exchange_data):
		# type: (Dict[str, Any]) -> Dict[str, Set[str]]
		"""
		Caches the presence of exchanges in each media condition based on
		exchange_data to set concentrations in concentrationsBasedOnNutrients().

		Args:
			exchange_data: dictionary of exchange data for all media conditions with keys:
				importUnconstrainedExchangeMolecules (dict[str, set[str]]): for each media ID key,
					exchange molecules (with location tag) that do not have an upper bound on their flux
				importConstrainedExchangeMolecules (dict[str, dict[str, float with mol/mass/time units]]):
					for each media ID key, constrained molecules (with location tag)
					with upper bound flux constraints

		Returns:
			sets of molecules IDs (with location tags) that can be imported for each
			media ID
		"""

		exchange_fluxes = {}

		all_unconstrained = exchange_data['importUnconstrainedExchangeMolecules']
		all_constrained = exchange_data['importConstrainedExchangeMolecules']

		for media in all_unconstrained:
			fluxes = set(all_unconstrained[media])
			fluxes.update([molecule
				for molecule, conc in all_constrained[media].items()
				if conc.asNumber() > 0])
			exchange_fluxes[media] = fluxes

		return exchange_fluxes

	def _addMoleculeAmounts(self, equilibriumReactions, concDict):
		moleculeSetAmounts = {}
		for reaction in equilibriumReactions:
			# We only want to do this for species with standard Michaelis-Menten kinetics initially
			if len(reaction["stoichiometry"]) != 3:
				continue

			moleculeName = [x["molecule"].encode("utf-8") for x in reaction["stoichiometry"] if x["type"] == "metabolite"][0]
			amountToSet = 1e-4
			moleculeSetAmounts[moleculeName + "[p]"] = amountToSet * self.units
			moleculeSetAmounts[moleculeName + "[c]"] = amountToSet * self.units

		for moleculeName, scaleFactor in self.moleculeScaleFactors.iteritems():
			moleculeSetAmounts[moleculeName] = scaleFactor * concDict[moleculeName] * self.units
		return moleculeSetAmounts


class Boundary(object):
	'''
	Boundary provides an interface between metabolism and the environment.
	This class builds exchange_data_dict, with keys for the five categories
	of external molecules that set up the FBA problem space (see _getExchangeDataDict for a description).
	'''
	def __init__(self, raw_data, sim_data):

		self.import_constraint_threshold = IMPORT_CONSTRAINT_THRESHOLD
		self.env_to_exchange_map = sim_data.external_state.environment.env_to_exchange_map

		self.all_external_exchange_molecules = self._getAllExternalExchangeMolecules(raw_data)
		self.secretion_exchange_molecules = self._getSecretionExchangeMolecules(raw_data)
		self.exchange_data_dict = self._getExchangeDataDict(sim_data)

	def _getAllExternalExchangeMolecules(self, raw_data):
		'''
		Returns:
			list[str]: all external exchange molecules
		'''
		externalExchangeData = []
		# initiate all molecules with 0 concentrations
		for row in raw_data.condition.environment_molecules:
			externalExchangeData.append(row["molecule id"] + row["exchange molecule location"])

		return externalExchangeData

	def _getSecretionExchangeMolecules(self, raw_data):
		'''
		Returns:
			set[str]: all secretion exchange molecules
		'''
		secretionExchangeMolecules = []
		for secretion in raw_data.secretions:
			if secretion["lower bound"] and secretion["upper bound"]:
				# "non-growth associated maintenance", not included in our metabolic model
				continue
			else:
				secretionExchangeMolecules.append(secretion["molecule id"])

		return set(secretionExchangeMolecules)

	def _getExchangeDataDict(self, sim_data):
		'''
		Returns:
			dict[str, Any]: keys are the five exchange_data variables with the following keys:
				externalExchangeMolecules (dict[str, set[str]]): for each media ID key,
					all exchange molecules (with location tag), includes both import
					and secretion exchanged molecules
				importExchangeMolecules (dict[str, set[str]]): for each media ID key,
					molecules (with location tag) that can be imported from the
					environment into the cell
				importConstrainedExchangeMolecules (dict[str, dict[str, float with mol/mass/time units]]):
					for each media ID key, constrained molecules (with location tag)
					with upper bound flux constraints
				importUnconstrainedExchangeMolecules (dict[str, set[str]]): for each media ID key,
					exchange molecules (with location tag) that do not have an upper bound on their flux
				secretionExchangeMolecules (set[str]): molecules (with location tag)
					that can be secreted by the cell into the environment
		'''

		saved_media = sim_data.external_state.environment.saved_media

		externalExchangeMolecules = {}
		importExchangeMolecules = {}
		importConstrainedExchangeMolecules = {}
		importUnconstrainedExchangeMolecules = {}
		secretionExchangeMolecules = self.secretion_exchange_molecules

		for environment_name, molecules in saved_media.iteritems():
			exchange_data = self.exchangeDataFromConcentrations(molecules)

			externalExchangeMolecules[environment_name] = exchange_data['externalExchangeMolecules']
			importExchangeMolecules[environment_name] = exchange_data['importExchangeMolecules']
			importConstrainedExchangeMolecules[environment_name] = exchange_data['importConstrainedExchangeMolecules']
			importUnconstrainedExchangeMolecules[environment_name] = exchange_data['importUnconstrainedExchangeMolecules']

		return {
			"externalExchangeMolecules": externalExchangeMolecules,
			"importExchangeMolecules": importExchangeMolecules,
			"importConstrainedExchangeMolecules": importConstrainedExchangeMolecules,
			"importUnconstrainedExchangeMolecules": importUnconstrainedExchangeMolecules,
			"secretionExchangeMolecules": secretionExchangeMolecules,
		}

	def exchangeDataFromConcentrations(self, molecules):
		# type: (Dict[str, float]) -> Dict[str, Any]
		'''
		Update importExchangeMolecules for FBA based on current nutrient concentrations.
		This provides a simple type of transport to accommodate changing nutrient
		concentrations in the environment. Transport is modeled as a binary switch:
		When there is a high concentrations of environment nutrients, transporters
		are unconstrained and nutrients are transported as needed by metabolism.
		When concentrations fall below the threshold, that nutrient's transport
		is constrained to max flux of 0.

		Args:
			molecules: external molecules (no location tag) with external concentration,
				concentration can be inf

		Returns dict with the following keys:
			externalExchangeMolecules (set[str]): all exchange molecules (with
				location tag), includes both import and secretion exchanged molecules
			importExchangeMolecules (set[str]): molecules (with location tag) that
				can be imported from the environment into the cell
			importConstrainedExchangeMolecules (dict[str, float with mol/mass/time units]):
				constrained molecules (with location tag) with upper bound flux constraints
			importUnconstrainedExchangeMolecules (set[str]): exchange molecules
				(with location tag) that do not have an upper bound on their flux
			secretionExchangeMolecules (set[str]): molecules (with location tag)
				that can be secreted by the cell into the environment
		'''

		externalExchangeMolecules = set()
		importExchangeMolecules = set()
		secretionExchangeMolecules = self.secretion_exchange_molecules

		glc_id = 'GLC[p]'
		oxygen_id = 'OXYGEN-MOLECULE[p]'

		exchange_molecules = {self.env_to_exchange_map[mol]: conc for mol, conc in molecules.iteritems()}

		# Unconstrained uptake if greater than import threshold
		importUnconstrainedExchangeMolecules = {molecule_id
			for molecule_id, concentration in exchange_molecules.items()
			if concentration >= self.import_constraint_threshold}
		importExchangeMolecules.update(importUnconstrainedExchangeMolecules)
		externalExchangeMolecules.update(importUnconstrainedExchangeMolecules)

		# Constrain molecules below import threshold at 0
		importConstrainedExchangeMolecules = {molecule_id: 0. * (units.mmol / units.g / units.h)
			for molecule_id, concentration in exchange_molecules.items()
			if concentration < self.import_constraint_threshold}

		# Limit glucose uptake if present depending on the presence of oxygen
		if glc_id in importUnconstrainedExchangeMolecules:
			if oxygen_id in importUnconstrainedExchangeMolecules:
				importConstrainedExchangeMolecules[glc_id] = 20. * (units.mmol / units.g / units.h)
			else:
				importConstrainedExchangeMolecules[glc_id] = 100. * (units.mmol / units.g / units.h)
			importUnconstrainedExchangeMolecules.remove(glc_id)

		externalExchangeMolecules.update(secretionExchangeMolecules)

		return {
			"externalExchangeMolecules": externalExchangeMolecules,
			"importExchangeMolecules": importExchangeMolecules,
			"importConstrainedExchangeMolecules": importConstrainedExchangeMolecules,
			"importUnconstrainedExchangeMolecules": importUnconstrainedExchangeMolecules,
			"secretionExchangeMolecules": secretionExchangeMolecules,
		}

	def exchangeDataFromMedia(self, media_label):
		'''
		Returns:
			dict: exchange_data for a media_label saved in exchange_data_dict.
		'''

		externalExchangeMolecules = self.exchange_data_dict['externalExchangeMolecules'][media_label]
		importExchangeMolecules = self.exchange_data_dict['importExchangeMolecules'][media_label]
		importConstrainedExchangeMolecules = self.exchange_data_dict['importConstrainedExchangeMolecules'][media_label]
		importUnconstrainedExchangeMolecules = self.exchange_data_dict['importUnconstrainedExchangeMolecules'][media_label]
		secretionExchangeMolecules = self.exchange_data_dict['secretionExchangeMolecules']

		return {
			"externalExchangeMolecules": externalExchangeMolecules,
			"importExchangeMolecules": importExchangeMolecules,
			"importConstrainedExchangeMolecules": importConstrainedExchangeMolecules,
			"importUnconstrainedExchangeMolecules": importUnconstrainedExchangeMolecules,
			"secretionExchangeMolecules": secretionExchangeMolecules,
		}

	def getImportConstraints(self, exchange_data):
		'''
		Returns:
			import_constraint (list[bool]): the indices of all importConstrainedExchangeMolecules
				in self.all_external_exchange_molecules are true, the rest as false.
			import_exchange (list[bool]): the indices of all importExchangeMolecules
				in self.all_external_exchange_molecules are true, the rest as false.
		'''

		# molecules from all_external_exchange_molecules set to 'true' if they are current importExchangeMolecules.
		import_exchange = [
			molecule_id in exchange_data['importExchangeMolecules']
			for molecule_id in self.all_external_exchange_molecules
			]

		# molecules from all_external_exchange_molecules set to 'true' if they are current importConstrainedExchangeMolecules.
		import_constraint = [
			molecule_id in exchange_data['importConstrainedExchangeMolecules']
			for molecule_id in self.all_external_exchange_molecules
			]

		return import_exchange, import_constraint
