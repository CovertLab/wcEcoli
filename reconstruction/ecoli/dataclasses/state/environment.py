"""
SimulationData for environment molecules state

@author: Eran Agmon
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 05/17/2018
"""

from __future__ import division

import numpy as np

from wholecell.utils import units
from wholecell.utils.unit_struct_array import UnitStructArray

from reconstruction.ecoli.dataclasses.state.stateFunctions import addToStateEnvironment

VERBOSE = True

class Environment(object):
	""" EnvironmentMolecules """

	def __init__(self, raw_data, sim_data):
		environmentData = np.zeros(
			0,
			dtype = [
				("id", "a50"),
				# TODO (Eran) -- mass might not be needed here
				("mass", "{}f8".format(len(sim_data.molecular_weight_order))),
				]
			)

		# Add units to values
		field_units = {
			"id"		:	None,
			"mass"				:	units.g / units.mol,
			}

		self.environmentData = UnitStructArray(environmentData, field_units)

	def addToEnvironmentState(self, ids, masses):
		self.environmentData = addToStateEnvironment(self.environmentData, ids, masses)


	# TODO (ERAN) -- the functions below were brought in from simulation_data.py, once integrated in state they should be removed from there
	def addConditionData(self, raw_data):
		self.conditionToDoublingTime = dict([(x["condition"].encode("utf-8"), x["doubling time"]) for x in raw_data.condition.condition_defs])
		self.conditionActiveTfs = dict([(x["condition"].encode("utf-8"), x["active TFs"]) for x in raw_data.condition.condition_defs])

		abbrToActiveId = dict([(x["TF"].encode("utf-8"), x["activeId"].encode("utf-8").split(", ")) for x in raw_data.tfIds if len(x["activeId"]) > 0])

		geneIdToRnaId = dict([(x["id"].encode("utf-8"), x["rnaId"].encode("utf-8")) for x in raw_data.genes])

		abbrToRnaId = dict(
			[(x["symbol"].encode("utf-8"), x["rnaId"].encode("utf-8")) for x in raw_data.genes] +
			[(x["name"].encode("utf-8"), geneIdToRnaId[x["geneId"].encode("utf-8")]) for x in raw_data.translationEfficiency if x["geneId"] != "#N/A"]
			)

		self.tfToFC = {}
		self.tfToDirection = {}
		notFound = []
		for row in raw_data.foldChanges:
			tf = abbrToActiveId[row["TF"].encode("utf-8")][0]
			try:
				target = abbrToRnaId[row["Target"].encode("utf-8")]
			except KeyError:
				notFound.append(row["Target"].encode("utf-8"))
				continue
			if tf not in self.tfToFC:
				self.tfToFC[tf] = {}
				self.tfToDirection[tf] = {}
			FC = row["F_avg"]
			if row["Regulation_direct"] < 0:
				FC *= -1.
				self.tfToDirection[tf][target] = -1
			else:
				self.tfToDirection[tf][target] = 1
			FC = 2**FC
			self.tfToFC[tf][target] = FC

		if VERBOSE:
			print "The following target genes listed in foldChanges.tsv have no corresponding entry in genes.tsv:"
			for item in notFound:
				print item


		self.tfToActiveInactiveConds = {}
		for row in raw_data.condition.tf_condition:
			tf = row["active TF"].encode("utf-8")
			activeGenotype = row["active genotype perturbations"]
			activeNutrients = row["active nutrients"].encode("utf-8")
			inactiveGenotype = row["inactive genotype perturbations"]
			inactiveNutrients = row["inactive nutrients"].encode("utf-8")

			if tf not in self.tfToActiveInactiveConds:
				self.tfToActiveInactiveConds[tf] = {}
			else:
				print "Warning: overwriting TF fold change conditions for %s" % tf

			self.tfToActiveInactiveConds[tf]["active genotype perturbations"] = activeGenotype
			self.tfToActiveInactiveConds[tf]["active nutrients"] = activeNutrients
			self.tfToActiveInactiveConds[tf]["inactive genotype perturbations"] = inactiveGenotype
			self.tfToActiveInactiveConds[tf]["inactive nutrients"] = inactiveNutrients

		self.conditions = {}
		for row in raw_data.condition.condition_defs:
			condition = row["condition"].encode("utf-8")
			self.conditions[condition] = {}
			self.conditions[condition]["nutrients"] = row["nutrients"].encode("utf-8")
			self.conditions[condition]["perturbations"] = row["genotype perturbations"]

		for tf in sorted(self.tfToActiveInactiveConds):
			activeCondition = tf + "__active"
			inactiveCondition = tf + "__inactive"

			if activeCondition in self.conditionActiveTfs:
				del self.conditionActiveTfs[activeCondition]
			if inactiveCondition in self.conditionActiveTfs:
				del self.conditionActiveTfs[inactiveCondition]

			self.conditions[activeCondition] = {}
			self.conditions[inactiveCondition] = {}
			self.conditions[activeCondition]["nutrients"] = self.tfToActiveInactiveConds[tf]["active nutrients"]
			self.conditions[inactiveCondition]["nutrients"] = self.tfToActiveInactiveConds[tf]["inactive nutrients"]
			self.conditions[activeCondition]["perturbations"] = self.tfToActiveInactiveConds[tf]["active genotype perturbations"]
			self.conditions[inactiveCondition]["perturbations"] = self.tfToActiveInactiveConds[tf]["inactive genotype perturbations"]

		self.nutrientsTimeSeries = {}
		for label in dir(raw_data.condition.timeseries):
			if label.startswith("__"):
				continue

			# TODO (Eran) -- what is collections and is it needed?
			# self.nutrientsTimeSeries[label] = collections.deque()
			# timeseries = getattr(raw_data.condition.timeseries, label)
			# for row in timeseries:
			# 	self.nutrientsTimeSeries[label].append((
			# 		row["time"].asNumber(units.s),
			# 		row["nutrients"].encode("utf-8")
			# 		))

		self.nutrientToDoublingTime = {}
		for condition in self.conditionToDoublingTime:
			if len(self.conditions[condition]["perturbations"]) > 0:
				continue
			nutrientLabel = self.conditions[condition]["nutrients"]
			if nutrientLabel in self.nutrientToDoublingTime and self.conditionToDoublingTime[condition] != self.nutrientToDoublingTime[nutrientLabel]:
				raise Exception, "Multiple doubling times correspond to the same media conditions"
			self.nutrientToDoublingTime[nutrientLabel] = self.conditionToDoublingTime[condition]


	def getNutrientData(self, raw_data):

		externalExchangeMolecules = {}
		importExchangeMolecules = {}
		secretionExchangeMolecules = set()
		importConstrainedExchangeMolecules = {}
		importUnconstrainedExchangeMolecules = {}
		nutrientsList = [(x, getattr(raw_data.condition.nutrient, x)) for x in dir(raw_data.condition.nutrient) if not x.startswith("__")]
		for nutrientsName, nutrients in nutrientsList:
			externalExchangeMolecules[nutrientsName] = set()
			importExchangeMolecules[nutrientsName] = set()
			importConstrainedExchangeMolecules[nutrientsName] = {}
			importUnconstrainedExchangeMolecules[nutrientsName] = []
			for nutrient in nutrients:
				if not np.isnan(nutrient["lower bound"].asNumber()) and not np.isnan(nutrient["upper bound"].asNumber()):
					continue
				elif not np.isnan(nutrient["upper bound"].asNumber()):
					importConstrainedExchangeMolecules[nutrientsName][nutrient["molecule id"]] = nutrient["upper bound"]
					externalExchangeMolecules[nutrientsName].add(nutrient["molecule id"])
					importExchangeMolecules[nutrientsName].add(nutrient["molecule id"])
				else:
					importUnconstrainedExchangeMolecules[nutrientsName].append(nutrient["molecule id"])
					externalExchangeMolecules[nutrientsName].add(nutrient["molecule id"])
					importExchangeMolecules[nutrientsName].add(nutrient["molecule id"])

			for secretion in raw_data.secretions:
				if secretion["lower bound"] and secretion["upper bound"]:
					# "non-growth associated maintenance", not included in our metabolic model
					continue

				else:
					externalExchangeMolecules[nutrientsName].add(secretion["molecule id"])
					secretionExchangeMolecules.add(secretion["molecule id"])

			externalExchangeMolecules[nutrientsName] = sorted(externalExchangeMolecules[nutrientsName])
			importExchangeMolecules[nutrientsName] = sorted(importExchangeMolecules[nutrientsName])
		secretionExchangeMolecules = sorted(secretionExchangeMolecules)

		return {
			"externalExchangeMolecules": externalExchangeMolecules,
			"importExchangeMolecules": importExchangeMolecules,
			"importConstrainedExchangeMolecules": importConstrainedExchangeMolecules,
			"importUnconstrainedExchangeMolecules": importUnconstrainedExchangeMolecules,
			"secretionExchangeMolecules": secretionExchangeMolecules,
		}