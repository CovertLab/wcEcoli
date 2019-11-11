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
		self.is_mrna = np.where(sim_data.process.transcription.rnaData['isMRna'])[0]
		self.mrna_data = sim_data.process.transcription.rnaData[self.is_mrna]

		self._buildRnaIndexToMonomerMapping(raw_data, sim_data)
		self._buildMonomerIndexToRnaMapping(raw_data, sim_data)
		#self._check_if_monomer_in_multiple_tus(tu_dict, monomer)

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
			self.rnaIndexToMonomerMapping_new.append(set_indices)
		self.rnaIndexToMonomerMapping_new = np.array(self.rnaIndexToMonomerMapping_new)
		'''
		self.rnaIndexToMonomerMapping = np.array(
			[np.where(x == sim_data.process.transcription.rnaData["id"])[0][0] 
			for x in sim_data.process.translation.monomerData["rnaId"]])
		'''
	def check_if_monomer_in_multiple_tus(self, tu_dict, monomer):
		all_overlapping_tus = []
		for key, value in tu_dict.items():
			overlapping_tu = []
			if monomer in key:
				print monomer + " is in the tu " + key
				overlapping_tu.append(key)
			if overlapping_tu:
				all_overlapping_tus.append(overlapping_tu)
			#if monomer in key:
				#print key, monomer
		if all_overlapping_tus:
			print "Overlapping TUS are" 
			for tu in all_overlapping_tus:
				print tu
		return


	def _buildMonomerIndexToRnaMapping(self, raw_data, sim_data):
		'''
		Input:
		sim_data.process.transcription.rnaData and sim_data.process.monomerData.

		Output:

		'''

		mrna = sim_data.process.transcription.rnaData[
			np.where(sim_data.process.transcription.rnaData['isMRna'])[0]]
		mrna_tu_count = len(self.mrna_data)
		monomer_count = len(sim_data.process.translation.monomerData)

		# This matrix is a transform from monomer space to mrna space.
		# It can be applied by writing:
		#     np.matmul(monomer_vector, self.monomerToRnaTransform)
		# which returns a vector of len(mrna)
		self.monomerToMrnaTransform = np.zeros((monomer_count, mrna_tu_count))

		#Pull in transcription unit data and convert to a dictionary for quicker access.
		#tu_counts = raw_data.transcription_units

		

		tu_dict = {tu['tu_id'] : tu['tu_count'] for tu in raw_data.transcription_units}

		

		tu_count_update = {}
		for key, value in tu_dict.items():
			if '_' in key:
				polycistronic_tu_count = value
				monocistronic_tu_in_poly = str(key).split('_')
				#first draft of this assumed the polycistronic mRNA would be the only one..
				#but in the full set of data there might be a ton of polys that encode for the 
				#same mRNA. Need to look for all the overlaps.
				print "Working on TU " + key 
				#look for the union in a list.
				fraction_count_poly = []
				#work on this once the simpler version is working.
				#look if the mono is actually present as a TU, if its not, move on.
				for mono in monocistronic_tu_in_poly:
					print "Looking for monomer " + mono
					self.check_if_monomer_in_multiple_tus(tu_dict, mono)
					
					#check if components of the TU are present in other TUS.
					
					if mono in tu_dict:
						print "Mono " + mono + " in tu_dict"
						fraction_count_mono = tu_dict[mono] / (value + tu_dict[mono])
						fraction_count_poly.append(value / (value + tu_dict[mono]))
						tu_count_update.update({mono:fraction_count_mono})
						tu_count_update.update({str(key):fraction_count_poly})
				# put something to check if monomers are present in multiple tus
				#if key not in tu_count_update:
					#tu_count_update.update({key: value})

		import pdb; pdb.set_trace()
			#print key
		y = []
		for x in self.mrna_data['monomerSet']:
			if len(x) > 1:
				y.append(x)

		for gene in y:
			#import pdb; pdb.set_trace()
			overlapping_tus = set(gene).intersection(list(self.mrna_data['monomerSet'])) 

			
		

				


		monomer_id_indexes = {}
		for idx, row in enumerate(sim_data.process.translation.monomerData):
			monomer_id_indexes[row['id']] = idx

		for mrna_index, mrna_data in enumerate(self.mrna_data):
			monomer_set = mrna_data['monomerSet']

			# TODO(Ryan): this number should be based on data
			value = 1.0 / len(monomer_set)

			for monomer_id in monomer_set:
				monomer_index = monomer_id_indexes[monomer_id]
				
				self.monomerToMrnaTransform[monomer_index][mrna_index] = value

		# TODO(Ryan): check to see if this is a legit way to find the inverse mapping?
		self.mrnaToMonomerTransform = self.monomerToMrnaTransform.T


		# ---------------------
		# This part is not used right now, mapping that Mialy added.
		# this mapping is a vector like the old one
		self.monomerIndexToRnaMapping_all = []
		for rna_row in sim_data.process.transcription.rnaData:
			set_indices = []
			for monomer_id in rna_row['monomerSet']:
				set_indices.append(monomer_id_indexes[monomer_id])
			self.monomerIndexToRnaMapping_all.append(set_indices)

		#remove all empty lists.
		self.monomerIndexToRnaMapping_new = [x for x in self.monomerIndexToRnaMapping_all if x != []]
		self.monomerIndexToRnaMapping_new = np.array(self.monomerIndexToRnaMapping_new)
		
		#new and old mapping does not match here!
		#new mapping maps directly from rnas to proteins. mimics the functioning of _buildRnaIndexToMonomerMapping
		#old mapping takes out all non-monomers from the indexing, so the indexing does not directly match 
		#the rnas and proteins indexing.

		self.monomerIndexToRnaMapping = np.array([np.where(x == sim_data.process.translation.monomerData["rnaId"])[0][0] for x in sim_data.process.transcription.rnaData["id"] if len(np.where(x == sim_data.process.translation.monomerData["rnaId"])[0])])
		#self.monomerIndexToRnaMapping_old_no_if = np.array([np.where(x == sim_data.process.translation.monomerData["rnaId"])[0][0] for x in sim_data.process.transcription.rnaData["id"]])

		#import pdb; pdb.set_trace()
	#def _buildRnaIndexToGeneMapping(self, raw_data, sim_data):
	#	self.rnaIndexToGeneMapping = np.array([np.where(x + "[c]" == sim_data.process.transcription.rnaData["id"])[0][0] for x in sim_data.process.replication.geneData["rnaId"]])
