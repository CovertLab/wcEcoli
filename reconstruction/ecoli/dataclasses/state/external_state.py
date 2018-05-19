"""
SimulationData state associated data

@author: Nick Ruggero
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 02/12/2015
"""

from __future__ import division

from wholecell.utils import units
from wholecell.utils.unit_struct_array import UnitStructArray

from reconstruction.ecoli.dataclasses.state.bulkMolecules import BulkMolecules
from reconstruction.ecoli.dataclasses.state.uniqueMolecules import UniqueMolecules
from reconstruction.ecoli.dataclasses.state.environment import Environment

from reconstruction.ecoli.dataclasses.state import stateFunctions as sf

import re
import numpy as np

# import collections #TODO (Eran) -- check what this is, remove
# VERBOSE = True

class ExternalState(object):
	""" External State """

	def __init__(self, raw_data, sim_data):

		self.environment = Environment(raw_data, sim_data)

		# self._buildEnvironmentMolecules(raw_data, sim_data)
		self._buildEnvironment(raw_data, sim_data)


	def _buildEnvironment(self, raw_data, sim_data):

		#condition data is brought into environment
		# self.environment.addConditionData(raw_data)
		self.environment.nutrientData = self._getNutrientData(raw_data)
		self.environment.condition = "basal"
		self.environment.nutrientsTimeSeriesLabel = "000000_basal"



	# def _buildEnvironmentMolecules(self, raw_data, sim_data):

		# # TODO (Eran) get nutrients instead of metabolites
		# # Set metabolites
		# nutrientIds = sf.createIdsWithCompartments(raw_data.metabolites)
		# nutrientMasses = units.g / units.mol * sf.createMetaboliteMassesByCompartments(raw_data.metabolites, 7, 11)
		#
		# self.environment.addToEnvironmentState(nutrientIds, nutrientMasses)
		#
		# # TODO (Eran) -- get environmental water
		# # Set water
		# waterIds = sf.createIdsWithCompartments(raw_data.water)
		# waterMasses = units.g / units.mol * sf.createMetaboliteMassesByCompartments(raw_data.water, 8, 11)
		#
		# self.environment.addToEnvironmentState(waterIds, waterMasses)



	# TODO (ERAN) -- the function below was brought in from simulation_data.py, once integrated in state they should be removed from there
	def _getNutrientData(self, raw_data):

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






