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

class State(object):
	""" State """

	def __init__(self, raw_data, sim_data):

		self.bulkMolecules = BulkMolecules(raw_data, sim_data)
		self.uniqueMolecules = UniqueMolecules(raw_data, sim_data)
		self.environment = Environment(raw_data, sim_data)

		self._buildBulkMolecules(raw_data, sim_data)
		self._buildUniqueMolecules(raw_data, sim_data)
		self._buildCompartments(raw_data, sim_data)
		# self._buildEnvironmentMolecules(raw_data, sim_data)
		self._buildEnvironment(raw_data, sim_data)


	def _buildBulkMolecules(self, raw_data, sim_data):

		# Set metabolites
		metaboliteIds = sf.createIdsWithCompartments(raw_data.metabolites)
		metaboliteMasses = units.g / units.mol * sf.createMetaboliteMassesByCompartments(raw_data.metabolites, 7, 11)

		self.bulkMolecules.addToBulkState(metaboliteIds, metaboliteMasses)

		# Set water
		waterIds = sf.createIdsWithCompartments(raw_data.water)
		waterMasses = units.g / units.mol * sf.createMetaboliteMassesByCompartments(raw_data.water, 8, 11)

		self.bulkMolecules.addToBulkState(waterIds, waterMasses)

		# Set RNA
		rnaIds = sf.createIdsWithCompartments(raw_data.rnas)
		rnaMasses = units.g / units.mol * sf.createMassesByCompartments(raw_data.rnas)

		self.bulkMolecules.addToBulkState(rnaIds, rnaMasses)

		# Set proteins
		proteinIds = sf.createIdsWithCompartments(raw_data.proteins)
		proteinMasses = units.g / units.mol * sf.createMassesByCompartments(raw_data.proteins)

		self.bulkMolecules.addToBulkState(proteinIds, proteinMasses)

		# Set complexes
		complexIds = sf.createIdsWithCompartments(raw_data.proteinComplexes)
		complexMasses = units.g / units.mol * sf.createMassesByCompartments(raw_data.proteinComplexes)

		self.bulkMolecules.addToBulkState(complexIds, complexMasses)

		# Set modified forms
		modifiedFormIds = sf.createIdsWithCompartments(raw_data.modifiedForms)
		modifiedFormMasses = units.g / units.mol * sf.createModifiedFormMassesByCompartments(raw_data.modifiedForms)

		self.bulkMolecules.addToBulkState(modifiedFormIds, modifiedFormMasses)

		# Set chromosome
		chromosomeIds = sf.createIdsWithCompartments(raw_data.chromosome)
		chromosomeMasses = units.g / units.mol * sf.createMassesByCompartments(raw_data.chromosome)

		self.bulkMolecules.addToBulkState(chromosomeIds, chromosomeMasses)

		# Set fragments
		test = []
		for x in raw_data.polymerized:
			if x['is_ntp']:
				if not x['is_end']:
					temp = x
					temp['id'] = x['id'].replace('Polymerized','Fragment')
					test.append(temp)
		fragmentsIds = sf.createIdsWithCompartments(test)
		fragmentsMasses = units.g / units.mol * sf.createMassesByCompartments(test)

		self.bulkMolecules.addToBulkState(fragmentsIds, fragmentsMasses)

	def _buildUniqueMolecules(self, raw_data, sim_data):
		# Add active RNA polymerase
		rnaPolyComplexMass = self.bulkMolecules.bulkData["mass"][self.bulkMolecules.bulkData["id"] == "APORNAP-CPLX[c]"]
		rnaPolyAttributes = {
				'rnaIndex' : 'i8',
				'transcriptLength' : 'i8'
				}
		self.uniqueMolecules.addToUniqueState('activeRnaPoly', rnaPolyAttributes, rnaPolyComplexMass)

		# Add active ribosome
		# TODO: This is a bad hack that works because in the fitter
		# I have forced expression to be these subunits only
		ribosome30SMass = self.bulkMolecules.bulkData["mass"][
		self.bulkMolecules.bulkData["id"] == sim_data.moleculeGroups.s30_fullComplex[0]
			]
		ribosome50SMass = self.bulkMolecules.bulkData["mass"][
		self.bulkMolecules.bulkData["id"] == sim_data.moleculeGroups.s50_fullComplex[0]
			]
		ribosomeMass = ribosome30SMass + ribosome50SMass
		ribosomeAttributes = {
				'proteinIndex' : 'i8',
				'peptideLength': 'i8'
				}
		self.uniqueMolecules.addToUniqueState('activeRibosome', ribosomeAttributes, ribosomeMass)

		# Add active DNA polymerase
		dnaPolyMass = units.g / units.mol * np.zeros_like(rnaPolyComplexMass) # NOTE: dnaPolymerases currently have no mass
		dnaPolymeraseAttributes = {
				'sequenceIdx' : 'i8',
				'sequenceLength' : 'i8',
				'replicationRound' : 'i8',
				'chromosomeIndex' : 'i8',
				}
		self.uniqueMolecules.addToUniqueState('dnaPolymerase', dnaPolymeraseAttributes, dnaPolyMass)

		# Origin of replication
		originMass = units.g / units.mol * np.zeros_like(rnaPolyComplexMass) # NOTE: origins currently have no mass
		originAttributes = {}
		self.uniqueMolecules.addToUniqueState('originOfReplication', originAttributes, originMass)

		# Full chromosome
		fullChromosomeMass = units.g / units.mol * np.zeros_like(rnaPolyComplexMass) # NOTE: origins currently have no mass
		fullChromosomeAttributes = {"division_time" : "f8"}
		self.uniqueMolecules.addToUniqueState('fullChromosome', fullChromosomeAttributes, fullChromosomeMass)

	def _buildCompartments(self, raw_data, sim_data):
		compartmentData = np.empty(len(raw_data.compartments),
			dtype = [('id','a20'),('compartmentAbbreviation', 'a1')])

		compartmentData['id'] = [x['id'] for x in raw_data.compartments]
		compartmentData['compartmentAbbreviation'] = [x['abbrev'] for x in raw_data.compartments]
		self.compartments = compartmentData

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

	def _buildEnvironment(self, raw_data, sim_data):

		# TODO (ERAN) -- condition data should all be brought into self.state.environment
		self.environment.addConditionData(raw_data)
		self.environment.nutrientData = self._getNutrientData(raw_data)
		self.environment.condition = "basal"
		self.environment.nutrientsTimeSeriesLabel = "000000_basal"

		import ipdb;
		ipdb.set_trace()

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






