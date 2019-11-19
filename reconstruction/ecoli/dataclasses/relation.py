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
	def find_overlapping_tus(self, tu_dict):
		'''
		Input: Dictionary of all the transcription units
		Output: A dictionary where the keys, are the monocistrons of every multi-gene 
		transcription unit, and the values are a list of all the TUs that monocistron belongs to.

		This function also adds the locations of the mRNAs which are assumed
		to all be in the cytoplasm.
		'''
		overlapping_tu_dict = {}
		polycistronic_tus = []

		for tu in tu_dict.keys():
			# Find all polycistronic RNA's - not exclusive to mRNA's
			if '_' in tu:
				polycistronic_tus.append([tu, str(tu).split('_')])
		for polycistron in polycistronic_tus:
			for monocistron in polycistron[1]:
				for tu in tu_dict.keys():
					#Check all the TU's that a monocistron is present in.
					if monocistron in tu:
						monocistron_loc = monocistron + '_RNA[c]'
						tu_loc = tu + '_RNA[c]'
						overlapping_tu_dict.setdefault(monocistron_loc, [])
						# Want to check if a key is already present as a value to
						# prevent adding it multiple times.
						if tu_loc  not in overlapping_tu_dict[monocistron_loc]:
							overlapping_tu_dict[monocistron_loc].append(tu_loc)

		return overlapping_tu_dict, polycistronic_tus


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
		#self.monomerToMrnaTransform = np.zeros((monomer_count, mrna_tu_count))

		tu_dict = {tu['tu_id']: tu['tu_count'] for tu in raw_data.transcription_units}
		overlapping_tu_dict, polycistronic_tus = self.find_overlapping_tus(tu_dict)
		
		
		#Note this is actually matching mono-cistrons with polycistrons, with the
		#expectation that proteins monomers are still linked to a single monocistrons. 
		tu_fraction_dict = {}
		for key, value in overlapping_tu_dict.items():
			tu_fraction_dict.setdefault(key, [])
			count_sum = 0
			#first calculate the sum counts of the tus that belong to each monomer.
			for v in value:
				#Need to remove the location tag from v in order to find it in the tu_dict.
				count_sum += tu_dict[v[:-7]]
			for v in value:
				tu_fraction_dict[key].append(tu_dict[v[:-7]] / count_sum)

		overlapping_tu_fractions = {}
		for monomer_monocistron, overlapping_tus in overlapping_tu_dict.items():
			for mrna_monocistron, tu_fraction in tu_fraction_dict.items():
				if monomer_monocistron == mrna_monocistron:
					overlapping_tu_fractions[monomer_monocistron] = {fraction: tu_fraction[index] 
						for index, fraction in enumerate(overlapping_tus)}
		




		

		# Now need to fit these values into the monomer/mRNA matrix
		# Note sim_data contains the locations, so will probably have to add
		# the locations to the TU's -- all should be in the cytoplasm

		#find row first: this is from monomer data
		#find rna's associated with monomer from rnaSet within monomer data
		#find column from mrnaData.
		#make sure length of rnaSet and the fractions associated with the monomer are correct.

		# Need to go back and do 1:1 mapping
		self.monomer_to_mrna_transform = np.zeros((monomer_count, mrna_tu_count))
		self.mrna_to_monomer_transform = np.zeros((mrna_tu_count, monomer_count))

		monomer_data = sim_data.process.translation.monomerData

		#create a dictionary of mRNA IDs attached to their indexes.
		mrna_id_index_dict = {mrna_id: index for index, mrna_id in enumerate(mrna['id'])}
		monomer_id_index_dict = {rnaId: index for index, rnaId in enumerate(monomer_data['rnaId'])}

		#this is a bit of a hack since for some reason the rnaIds for 
		#the lac genes arent assigned a location. Need to figure out why this
		#is happening, what other locations are available?

		'''
		The following works to assign a 1:1 ratio for single transcripts identified
		in the monomer data to single transcripts in the mRNA data.


		TO DO: Put this in its own function.
		'''

		for monocistron, moncistron_index in monomer_id_index_dict.items():
			for tu, tu_index in mrna_id_index_dict.items():
				if monocistron == tu:
					self.monomer_to_mrna_transform[moncistron_index, tu_index] = 1
					self.mrna_to_monomer_transform[tu_index ,moncistron_index] = 1
				elif monocistron[:-3] + '[None]' == tu:
					self.monomer_to_mrna_transform[moncistron_index, tu_index] = 1
					self.mrna_to_monomer_transform[tu_index, moncistron_index] = 1		

		'''
		TO DO: Put this in its own function.

		Note the above code would not do anything for the lac genes since they
		are not available as monocistrons. Need to capture all poly cistronic mRNAs.

		The following code should add the polycistrons that are linked to each monomer.
		'''

		for monocistron, transcription_units in overlapping_tu_fractions.items():
			#get the index of the monocistron according to the monomer indexing.
			if monocistron in monomer_id_index_dict:
				monocistron_index = monomer_id_index_dict[monocistron]
			elif monocistron[:-3] + '[None]' in monomer_id_index_dict:
				monocistron_index = monomer_id_index_dict[monocistron[:-3] + '[None]']
			#Get the index of the transcription_units associated to the mono-cistron according to the mRNA indexing.
			print "Monocistron is " + monocistron, monocistron_index
			for tu in transcription_units:
				if tu in mrna_id_index_dict:
					tu_index = mrna_id_index_dict[tu]
					self.monomer_to_mrna_transform[monocistron_index, tu_index] = overlapping_tu_fractions[monocistron][tu]
					self.mrna_to_monomer_transform[tu_index, monocistron_index] = 1
				print "Tu is " + tu, tu_index
		

		'''
		This is a quick check that can be implemented once a bunch of transcription units are added.
		Basically, the sum of each row should ideally = 1. For one row currently the value is 
		0.9999999999999999, likely a rounding issue. 

		To Do: Put this in its own function. 

		'''
		for index, row in enumerate(self.monomer_to_mrna_transform):
			if sum(row) < 1:
				print "Row " + str(index) + " is less than 1: " + str(sum(row))
			elif sum(row) > 1:
				print "Row " + str(index) + " is greater than 1: " + str(sum(row))


		import pdb; pdb.set_trace()



		'''
		#This is Ryans	
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
	 	# Mialy's note: This will not work.
		self.mrnaToMonomerTransform = self.monomerToMrnaTransform.T

		'''
		

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
