#!/usr/bin/env python

"""
EnzymeKinetics

Takes in enzyme kinetics data on initialization, and returns dicts of rate estimates when passed
metabolite and enzyme concentrations at runtime.

@author: Morgan Paull
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 2/8/2016
"""

import numpy as np

from wholecell.utils import units
import re
from Equation import Expression

COUNTS_UNITS = units.umol
TIME_UNITS = units.s
VOLUME_UNITS = units.L

class EnzymeKinetics(object):
	"""
	EnzymeKinetics

	Returns rate estimates from kinetic equation information stored in reactionRateInfo.
	"""

	def __init__(self, reactionRateInfo, kcatsOnly=False, noCustoms=False, moreThanKcat=False):

		# Set default reaction rate limit, to which reactions are set absent other information
		self.defaultRate = (COUNTS_UNITS / TIME_UNITS / VOLUME_UNITS) * np.inf

		# Load rate functions from enzymeKinetics.tsv flat file
		self.reactionRateInfo = reactionRateInfo

		## Filter the reactions as specified
		# Exclude any rows with more than a kcat
		if kcatsOnly:
			reactionRateInfoNew = {}
			for constraintID, reactionInfo in self.reactionRateInfo.iteritems():
				# Kcat-only reactions will have no kMs, kIs, or custom equations
				if len(reactionInfo["kM"]) or len(reactionInfo["kI"]) or reactionInfo["customRateEquation"]:
					continue
				reactionRateInfoNew[constraintID] = reactionInfo
			self.reactionRateInfo = reactionRateInfoNew

		# Exclude any custom equation rows
		if noCustoms:
			reactionRateInfoNew = {}
			for constraintID, reactionInfo in self.reactionRateInfo.iteritems():
				if reactionInfo["customRateEquation"] == None:
					reactionRateInfoNew[constraintID] = reactionInfo
			self.reactionRateInfo = reactionRateInfoNew

		# Throw out any kcat-only reactions
		if moreThanKcat:
			reactionRateInfoNew = {}
			for constraintID, reactionInfo in self.reactionRateInfo.iteritems():
				if len(reactionInfo["kM"]) or len(reactionInfo["kI"]) or reactionInfo["customRateEquation"]:
					reactionRateInfoNew[constraintID] = reactionInfo
			self.reactionRateInfo = reactionRateInfoNew

		self.allConstraintIDs = self.reactionRateInfo.keys()

		self.allReactionIDs = [x["reactionID"] for x in self.reactionRateInfo.values()]

		self.inputsChecked = False

	def checkKnownSubstratesAndEnzymes(self, metaboliteConcentrationsDict, enzymeConcentrationsDict, removeUnknowns=False):
		knownConstraints = {}
		unusableConstraints = {}
		unknownSubstrates = set()
		unknownEnzymes = set()
		unknownCustomVars = set()


		for constraintID, reactionInfo in self.reactionRateInfo.iteritems():
			keepReaction = True
			reactionType = reactionInfo["rateEquationType"]
			if reactionType == "standard":
				
				# Check if the substrates are known
				for substrateID in reactionInfo["substrateIDs"]:
					if substrateID not in metaboliteConcentrationsDict:
						unknownSubstrates.add(substrateID)
						unusableConstraints[constraintID] = reactionInfo
						keepReaction = False

				# Check if the enzymes are known
				for enzymeID in reactionInfo["enzymeIDs"]:
					if enzymeID not in enzymeConcentrationsDict:
						unknownEnzymes.add(enzymeID)
						unusableConstraints[constraintID] = reactionInfo
						keepReaction = False


			elif reactionType == "custom":

				for variable in reactionInfo["customParameterVariables"].values():
					if variable not in metaboliteConcentrationsDict:
						notSubstrate = True
					if variable not in enzymeConcentrationsDict:
						notEnzyme = True

					if notSubstrate and notEnzyme:
						unknownCustomVars.add(variable)
						unusableConstraints[constraintID] = reactionInfo
						keepReaction = False
			else:
				# Reaction type is unknown
				raise Exception("Reaction type '%s' is unknown. Must be either 'standard' or 'custom'." % (reactionType))

			# Keep the reaction only if both substrates and enzymes are known
			if keepReaction:
				knownConstraints[constraintID] = reactionInfo

		if removeUnknowns:
			self.reactionRateInfo = knownConstraints
			self.inputsChecked = True

		unknownVals = {"unknownSubstrates":unknownSubstrates, "unknownEnzymes":unknownEnzymes, "unknownCustomVars":unknownCustomVars}

		return knownConstraints, unusableConstraints, unknownVals


	def reactionRate(self, reactionInfo, metaboliteConcentrationsDict, enzymeConcentrationsDict):

		if reactionInfo["rateEquationType"] == "standard":
			return self.reactionStandard(reactionInfo, metaboliteConcentrationsDict, enzymeConcentrationsDict)
		elif reactionInfo["rateEquationType"] == "custom":
			return self.reactionCustom(reactionInfo, metaboliteConcentrationsDict, enzymeConcentrationsDict)
		else:
			raise NameError("rateEquationType %s not recognized! Must be either 'standard' or 'custom'." % (reactionInfo["rateEquationType"]))

	def reactionStandard(self, reactionInfo, metaboliteConcentrationsDict, enzymeConcentrationsDict):
		# Find the enzymes needed for this rate
		enzymeConc = enzymeConcentrationsDict[reactionInfo["enzymeIDs"][0]]
		kMs = reactionInfo["kM"]
		kIs = reactionInfo["kI"]

		rate = np.amax(reactionInfo["kcat"])*enzymeConc

		idx = 0
		for kM in reactionInfo["kM"]:
			substrateConc = metaboliteConcentrationsDict[reactionInfo["substrateIDs"][idx]]
			rate *= (substrateConc / (float(substrateConc) + kM))
			idx+=1

		for kI in reactionInfo["kI"]:
			try:
				substrateConc = metaboliteConcentrationsDict[reactionInfo["substrateIDs"][idx]]
				rate *= 1.0 / (1.0 + (substrateConc / kI))
				idx+=1
			except:
				import ipdb; ipdb.set_trace()

		return (COUNTS_UNITS / TIME_UNITS / VOLUME_UNITS) * rate

	def reactionCustom(self, reactionInfo, metaboliteConcentrationsDict, enzymeConcentrationsDict):
		enzymeConcArray = [enzymeConcentrationsDict[reactionInfo["enzymeIDs"][0]]]

		customParameters = reactionInfo["customParameters"]
		customParameterVariables = reactionInfo["customParameterVariables"]
		customParameterConstants = reactionInfo["customParameterConstantValues"]
		equationString = reactionInfo["customRateEquation"]
		parameterDefinitionArray = reactionInfo["customParameters"]

		parametersArray = customParameterConstants
		for customParameter in customParameters[len(customParameterConstants):]:
			variable = customParameterVariables[customParameter]
			if variable in enzymeConcentrationsDict:
				parametersArray.append(enzymeConcentrationsDict[variable])
			elif variable in metaboliteConcentrationsDict:
				parametersArray.append(metaboliteConcentrationsDict[variable])

		assert (len(parametersArray) == len(parameterDefinitionArray))

		customRateLaw = Expression(equationString, parameterDefinitionArray)

		return (COUNTS_UNITS / TIME_UNITS / VOLUME_UNITS) * customRateLaw(*parametersArray)


	def allConstraintsDict(self, metaboliteConcentrationsDict, enzymeConcentrationsDict):
		constraintsDict = {}
		for constraintID, reactionInfo in self.reactionRateInfo.iteritems():
			constraintsDict[constraintID] = self.reactionRate(reactionInfo, metaboliteConcentrationsDict, enzymeConcentrationsDict)

		return constraintsDict

	def allReactionsDict(self, metaboliteConcentrationsDict, enzymeConcentrationsDict):
		"""
		Create a dict of dicts from reactionID to constraintIDs for that reaction, to rates for each constraintID.
		"""
		reactionsDict = {}
		for constraintID, reactionInfo in self.reactionRateInfo.iteritems():
			reactionID = reactionInfo["reactionID"]
			if reactionID not in reactionsDict:
				reactionsDict[reactionID] = {}
			reactionsDict[reactionID][constraintID] = self.reactionRate(reactionInfo, metaboliteConcentrationsDict, enzymeConcentrationsDict)

		return reactionsDict