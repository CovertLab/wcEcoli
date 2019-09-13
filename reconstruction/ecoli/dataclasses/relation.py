"""
SimulationData relation functions

@author: Nick Ruggero
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 03/10/2015
"""

from __future__ import division

import re
import numpy as np

# Unit imports
from wholecell.utils import units
from wholecell.utils.unit_struct_array import UnitStructArray

class Relation(object):
	""" Relation """

	def __init__(self, raw_data, sim_data):
		self._buildRnaIndexToMonomerMapping(raw_data, sim_data)
		self._buildMonomerIndexToRnaMapping(raw_data, sim_data)
		#self._buildRnaIndexToGeneMapping(raw_data, sim_data)
	
	def _buildRnaIndexToMonomerMapping(self, raw_data, sim_data):
		'''
		Input:
		sim_data.process.transcription.rnaData and sim_data.process.monomerData.
		Output:
		np.ndarray, where each element index corresponds to the index of a monomer in monomerData, containing a list
		of the indices where the RNA's that can be used to make that monomer are located within rnaData.

		Example:
		self.rnaIndexToMonomerMapping[1454] = [2099, 2100]
		The monomer is: EG12197-MONOMER
		The RNA's are: "EG12197_RNA" and "EG12197_EG12144_RNA"
		'''

		rnaData_id_index = {}
		for idx, row in enumerate(sim_data.process.transcription.rnaData):
			rnaData_id_index[row['id']] = idx

		self.rnaIndexToMonomerMapping_new = []
		for protein_row in sim_data.process.translation.monomerData:
			set_indices = []
			for rna_id in protein_row['rnaSet']:
				set_indices.append(rnaData_id_index[rna_id])
			self.rnaIndexToMonomerMapping.append(set_indices)
		self.rnaIndexToMonomerMapping_new = np.array(self.rnaIndexToMonomerMapping)

		self.rnaIndexToMonomerMapping = np.array(
			[np.where(x == sim_data.process.transcription.rnaData["id"])[0][0] 
			for x in sim_data.process.translation.monomerData["rnaId"]])

	def _buildMonomerIndexToRnaMapping(self, raw_data, sim_data):
		'''
		Input:
		sim_data.process.transcription.rnaData and sim_data.process.monomerData.

		Output:

		'''
		monomerData_id_index = {}
		for idx, row in enumerate(sim_data.process.translation.monomerData):
			monomerData_id_index[row['id']] = idx

		self.monomerIndexToRnaMapping_all = []
		for rna_row in sim_data.process.transcription.rnaData:
			set_indices = []
			for monomer_id in rna_row['monomerSet']:
				set_indices.append(monomerData_id_index[monomer_id])
			self.monomerIndexToRnaMapping_all.append(set_indices)
		#remove all empty lists.
		self.monomerIndexToRnaMapping_new = [x for x in self.monomerIndexToRnaMapping_all if x != []]
		self.monomerIndexToRnaMapping_new = np.array(self.monomerIndexToRnaMapping)

		#new and old mapping does not match here!
		#new mapping maps directly from rnas to proteins. mimics the functioning of _buildRnaIndexToMonomerMapping
		#old mapping takes out all non-monomers from the indexing, so the indexing does not directly match 
		#the rnas and proteins indexing.

		self.monomerIndexToRnaMapping = np.array([np.where(x == sim_data.process.translation.monomerData["rnaId"])[0][0] for x in sim_data.process.transcription.rnaData["id"] if len(np.where(x == sim_data.process.translation.monomerData["rnaId"])[0])])
		#self.monomerIndexToRnaMapping_old_no_if = np.array([np.where(x == sim_data.process.translation.monomerData["rnaId"])[0][0] for x in sim_data.process.transcription.rnaData["id"]])

		#import pdb; pdb.set_trace()
	#def _buildRnaIndexToGeneMapping(self, raw_data, sim_data):
	#	self.rnaIndexToGeneMapping = np.array([np.where(x + "[c]" == sim_data.process.transcription.rnaData["id"])[0][0] for x in sim_data.process.replication.geneData["rnaId"]])