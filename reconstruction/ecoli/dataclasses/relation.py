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
	'''
	def make_individual_rna_ids(operon_rna_id):
			#move this to the top if we keep it
			SPLIT_DELIMITER = '_'
			#This solution will work for now, but will need to look for type and add the appropriate suffix
			#There has to be an easier way, maybe store the first gene in the operon_rnas file.
			rna_ids = re.split(SPLIT_DELIMITER, operon_rna_id)
			split_rna_ids = []
			for rna_id in rna_ids:
				split_rna_ids.append(rna_id + '_RNA')
			return split_rna_ids
	'''
	def _buildRnaIndexToMonomerMapping(self, raw_data, sim_data):
		#import pdb; pdb.set_trace()
		self.rnaIndexToMonomerMapping = np.array([np.where(x == sim_data.process.transcription.rnaData["id"])[0][0] for x in sim_data.process.translation.monomerData["rnaId"]])

	def _buildMonomerIndexToRnaMapping(self, raw_data, sim_data):
		self.monomerIndexToRnaMapping = np.array([np.where(x == sim_data.process.translation.monomerData["rnaId"])[0][0] for x in sim_data.process.transcription.rnaData["id"] if len(np.where(x == sim_data.process.translation.monomerData["rnaId"])[0])])

	#def _buildRnaIndexToGeneMapping(self, raw_data, sim_data):
	#	self.rnaIndexToGeneMapping = np.array([np.where(x + "[c]" == sim_data.process.transcription.rnaData["id"])[0][0] for x in sim_data.process.replication.geneData["rnaId"]])