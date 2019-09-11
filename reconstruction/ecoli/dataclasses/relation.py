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
		Rewriting this from the original to pull the relationships from operon_rnas.tsv
		rather than from rnas.tsv and proteins.tsv

		Will go through all protein monomers and will spit out all the relevant RNA indices within a list.
		Go into RNA
		'''
		
		#returns a dictionary where the keys are the RNA_ids and the values are the indices
		#where that ID appears in operon_rnas.tsv
		rnaData_id_index = {}
		for idx, row in enumerate(sim_data.process.transcription.rnaData):
			rnaData_id_index[row['id']] = idx

		
		import pdb; pdb.set_trace()
		#{rnaData_id_index[row['id']] = idx for idx, row in enumerate(sim_data.process.transcription.rnaData)}
		self.rnaIndexToMonomerMapping = []
		for protein_row in sim_data.process.translation.monomerData:
			set_indices = []
			for rna_id in protein_row['rnaSet']:
				set_indices.append(rnaData_id_index[rna_id])
			rnaData_id_index.append(set_indices)
		

		#self.rnaIndexToMonomerMapping = np.array([np.where(x == sim_data.process.transcription.rnaData["id"])[0][0] for x in sim_data.process.translation.monomerData["rnaId"]])

	def _buildMonomerIndexToRnaMapping(self, raw_data, sim_data):
		self.monomerIndexToRnaMapping = np.array([
			np.where(x == sim_data.process.translation.monomerData["rnaId"])[0][0] 
			for x in sim_data.process.transcription.rnaData["id"] 
			if len(np.where(x == sim_data.process.translation.monomerData["rnaId"])[0])])

	#def _buildRnaIndexToGeneMapping(self, raw_data, sim_data):
	#	self.rnaIndexToGeneMapping = np.array([np.where(x + "[c]" == sim_data.process.transcription.rnaData["id"])[0][0] for x in sim_data.process.replication.geneData["rnaId"]])