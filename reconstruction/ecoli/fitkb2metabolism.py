
from __future__ import division
import os
import collections

import numpy as np
import tables
import cvxopt

from wholecell.utils import units
from wholecell.utils.modular_fba import FluxBalanceAnalysis

COUNTS_UNITS = units.mmol
VOLUME_UNITS = units.L
MASS_UNITS = units.g

def fitKb_2_metabolism(kb, simOutDir, bulkAverageContainer, bulkDeviationContainer):

	# Load the simulation output

	## Effective biomass reaction
	with tables.open_file(os.path.join(simOutDir, "ConcentrationChange.hdf")) as h5file:
		time = h5file.root.ConcentrationChange.col("time")
		timeStep = h5file.root.ConcentrationChange.col("timeStep")

		# NOTE: units are M/s
		concentrationChange = h5file.root.ConcentrationChange.col("concentrationChange")

		names = h5file.root.names
		biomassMoleculeIDs = np.array(names.moleculeIDs.read())

	## Find the most extreme concentration flux, after removing the first few time steps

	# TODO: intelligent filtering - use the correlation coefficient?

	concentrationChange = concentrationChange[3:, :] # NOTE: magic number

	concentrationChangeMostPositive = concentrationChange.max(0)
	concentrationChangeMostNegative = concentrationChange.min(0)

	effectiveBiomassReaction = concentrationChangeMostPositive.copy()

	negativeIsMostExtreme = (np.abs(concentrationChangeMostNegative)
		> concentrationChangeMostPositive)

	effectiveBiomassReaction[negativeIsMostExtreme] = concentrationChangeMostNegative[negativeIsMostExtreme]

	biomassReaction = dict(zip(biomassMoleculeIDs, effectiveBiomassReaction*1000)) # conversion to mM/s
	biomassReaction.pop('PPGPP[c]')# TODO: Make this work with metabolism
	reactionRates = {reactionID:1 for reactionID in kb.metabolismReactionEnzymes.viewkeys()}

	fba = FluxBalanceAnalysis(
		kb.metabolismReactionStoich,
		kb.metabolismExternalExchangeMolecules,
		biomassReaction,
		reversibleReactions = kb.metabolismReversibleReactions,
		reactionEnzymes = kb.metabolismReactionEnzymes.copy(), # TODO: copy in class
		reactionRates = reactionRates,
		)

	# Set constraints
	## External molecules
	externalMoleculeIDs = fba.externalMoleculeIDs()

	initWaterMass = kb.avgCellWaterMassInit
	initDryMass = kb.avgCellDryMassInit

	initCellMass = (
		initWaterMass
		+ initDryMass
		)

	coefficient = initDryMass / initCellMass * kb.cellDensity * (1 * units.s)

	externalMoleculeLevels = kb.metabolismExchangeConstraints(
		externalMoleculeIDs,
		coefficient,
		COUNTS_UNITS / VOLUME_UNITS
		)

	fba.externalMoleculeLevelsIs(externalMoleculeLevels)

	## Set enzymes unlimited
	fba.enzymeLevelsIs(np.inf)

	## Constrain FBA to 100% efficacy
	fba.maxReactionFluxIs(fba._standardObjectiveReactionName, 1)

	## Solve and assert feasibility
	fba.run()

	assert abs(fba.objectiveReactionFlux() - 1) < 1e-10, "vMax fitting is infeasible"

	enzymeUsage = fba.enzymeUsage()

	TOLERANCE = 2

	min_vMax = enzymeUsage[enzymeUsage > 1e-16].min()

	enzyme_vMax = dict(zip(
		fba.enzymeIDs(), np.fmax(enzymeUsage, min_vMax) * TOLERANCE
		))

	kb.metabolismReactionMaxRates = { # NOTE: adding attribute to KB
		reactionID:enzyme_vMax[enzymeID]
		for reactionID, enzymeID in kb.metabolismReactionEnzymes.viewitems()
		}


def _expressionFitting():
	## Build the enzyme-fitting problem

	# TODO: write a class for setting up LP problems

	values = []
	rowIndexes = []
	colIndexes = []

	lowerValues = []
	lowerIndexes = []

	upperValues = []
	upperIndexes = []

	objValues = []
	objIndexes = []

	rowNames = []
	colNames = []

	### Set up reverse reactions

	dt = 1 * units.s

	reactionStoich = kb.metabolismReactionStoich.copy()
	reactionEnzymes = kb.metabolismReactionEnzymes.copy()
	reactionRates = kb.metabolismReactionRates(dt)

	for reactionID in kb.metabolismReversibleReactions:
		reverseReactionID = "{} (reverse)".format(reactionID)
		assert reverseReactionID not in reactionStoich.viewkeys()
		reactionStoich[reverseReactionID] = {
			moleculeID:-coeff
			for moleculeID, coeff in reactionStoich[reactionID].viewitems()
			}

		if reactionID in reactionEnzymes.viewkeys():
			reactionEnzymes[reverseReactionID] = reactionEnzymes[reactionID]

		if reactionID in reactionRates.viewkeys():
			reactionRates[reverseReactionID] = reactionRates[reactionID]

	### Set up metabolites and biochemical reactions

	for reactionID, stoich in reactionStoich.viewitems():
		assert reactionID not in colNames
		reactionIndex = len(colNames)
		colNames.append(reactionID)

		for moleculeID, coeff in stoich.viewitems():
			try:
				moleculeIndex = rowNames.index(moleculeID)

			except ValueError:
				moleculeIndex = len(rowNames)
				rowNames.append(moleculeID)

			rowIndexes.append(moleculeIndex)
			colIndexes.append(reactionIndex)
			values.append(coeff)

	### Set up exchange reactions

	initWaterMass = kb.avgCellWaterMassInit
	initDryMass = kb.avgCellDryMassInit

	initCellMass = initWaterMass + initDryMass

	initCellVolume = initCellMass / kb.cellDensity

	coefficient = initDryMass / initCellVolume * dt

	exchangeConstraints = kb.metabolismExchangeConstraints(
		kb.metabolismExternalExchangeMolecules,
		coefficient,
		units.mmol / units.L
		)

	for moleculeID, constraint in zip(kb.metabolismExternalExchangeMolecules, exchangeConstraints):
		exchangeID = "{} exchange".format(moleculeID)

		assert exchangeID not in colNames
		exchangeIndex = len(colNames)
		colNames.append(exchangeID)

		moleculeIndex = rowNames.index(moleculeID)

		rowIndexes.append(moleculeIndex)
		colIndexes.append(exchangeIndex)
		values.append(-1.0)

		lowerIndexes.append(exchangeIndex)
		lowerValues.append(-min(constraint, 1e6))

	### Set up biomass reaction

	biomassID = "biomass reaction"
	assert biomassID not in colNames
	biomassIndex = len(colNames)
	colNames.append(biomassID)

	effectiveBiomassReaction *= 10**3 # change to mmol

	for moleculeID, coeff in zip(biomassMoleculeIDs, effectiveBiomassReaction):
		moleculeIndex = rowNames.index(moleculeID)

		rowIndexes.append(moleculeIndex)
		colIndexes.append(biomassIndex)
		values.append(-coeff)

		lowerIndexes.append(biomassIndex)
		lowerValues.append(+1) # must be capable of producing 100% of the biomass in a step

	### Set up enzyme usage

	enzymeRatesAll = collections.defaultdict(set)

	for reactionID, enzymeID in reactionEnzymes.viewitems():
		reactionRate = reactionRates[reactionID]

		enzymeRatesAll[enzymeID].add(reactionRate)

	enzymeIDs = enzymeRatesAll.keys()
	perEnzymeRates = {
		enzymeID:max(enzymeRates)
		for enzymeID, enzymeRates in enzymeRatesAll.viewitems()
		}

	minimalEnzymeCounts = np.fmax(
		bulkAverageContainer.counts(enzymeIDs) - 2 * bulkDeviationContainer.counts(enzymeIDs),
		0
		)

	enzymeConc = (
		1 / kb.nAvogadro / initCellVolume * minimalEnzymeCounts
		).asNumber(units.mmol / units.L)

	fullEnzymeRates = {
		enzymeID:perEnzymeRates[enzymeID] * enzymeConc[index]
		for index, enzymeID in enumerate(enzymeIDs)
		}

	for enzymeID, rateConstraint in fullEnzymeRates.viewitems():
		assert enzymeID not in rowNames
		enzymeIndex = len(rowNames)
		rowNames.append(enzymeID)

		constraintID = "{} constraint".format(enzymeID)
		assert constraintID not in colNames
		constraintIndex = len(colNames)
		colNames.append(constraintID)

		excessID = "{} excess capacity".format(enzymeID)
		assert excessID not in colNames
		excessIndex = len(colNames)
		colNames.append(excessID)

		rowIndexes.append(enzymeIndex)
		colIndexes.append(constraintIndex)
		values.append(+1.0)

		upperIndexes.append(constraintIndex)
		upperValues.append(rateConstraint)

		rowIndexes.append(enzymeIndex)
		colIndexes.append(excessIndex)
		values.append(+1.0)

		objIndexes.append(excessIndex)
		objValues.append(+1.0) # TODO: weighting

	for reactionID, enzymeID in reactionEnzymes.viewitems():
		if reactionID not in reactionRates.viewkeys():
			raise Exception("This code was not designed to handle enzymatic constraints without annotated rates.")

		reactionIndex = colNames.index(reactionID)
		enzymeIndex = rowNames.index(enzymeID)

		rowIndexes.append(enzymeIndex)
		colIndexes.append(reactionIndex)
		values.append(-1)

	import cvxopt
	import cvxopt.solvers

	nRows = max(rowIndexes) + 1
	nCols = max(colIndexes) + 1

	assert len(values) == len(rowIndexes) == len(colIndexes)

	A = cvxopt.spmatrix(values, rowIndexes, colIndexes)

	b = cvxopt.matrix(np.zeros(nRows, np.float64))

	assert len(objIndexes) == len(objValues)

	objectiveFunction = np.zeros(nCols, np.float64)
	objectiveFunction[objIndexes] = objValues

	f = cvxopt.matrix(objectiveFunction)

	G = cvxopt.matrix(np.concatenate(
		[np.identity(nCols, np.float64), -np.identity(nCols, np.float64)]
		))

	assert len(upperIndexes) == len(upperValues)

	upperBound = np.empty(nCols, np.float64)
	upperBound.fill(1e6)
	upperBound[upperIndexes] = upperValues

	assert len(lowerIndexes) == len(lowerValues)

	lowerBound = np.empty(nCols, np.float64)
	lowerBound.fill(0)
	lowerBound[lowerIndexes] = lowerValues

	h = cvxopt.matrix(np.concatenate([upperBound, -lowerBound], axis = 0))

	oldOptions = cvxopt.solvers.options.copy()

	cvxopt.solvers.options["LPX_K_MSGLEV"] = 0

	solution = cvxopt.solvers.lp(f, G, h, A, b, solver = "glpk")

	cvxopt.solvers.options.update(oldOptions)

	import ipdb; ipdb.set_trace()
