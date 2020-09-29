"""
SimulationData relation functions
"""

from __future__ import absolute_import, division, print_function

import numpy as np


class Relation(object):
	""" Relation """

	def __init__(self, raw_data, sim_data):
		# pull mrna and monomer data
		self.is_mrna = np.where(sim_data.process.transcription.rna_data['is_mRNA'])[0]
		self.mrna = sim_data.process.transcription.rna_data[self.is_mrna]
		self.monomer = sim_data.process.translation.monomer_data

		self._buildRnaIndexToMonomerMapping(raw_data, sim_data)
		self._buildMonomerIndexToRnaMapping(raw_data, sim_data)

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

		rna_data_id_index = {}
		for idx, row in enumerate(sim_data.process.transcription.rna_data):
			rna_data_id_index[row['id']] = idx
		self.rnaIndexToMonomerMapping_new = []
		for protein_row in sim_data.process.translation.monomer_data:
			set_indices = []
			for rna_id in protein_row['rna_set']:
				set_indices.append(rna_data_id_index[rna_id])
			self.rnaIndexToMonomerMapping_new.append(set_indices)
		self.rnaIndexToMonomerMapping_new = np.array(self.rnaIndexToMonomerMapping_new)
		'''
		self.rnaIndexToMonomerMapping = np.array(
			[np.where(x == sim_data.process.transcription.rnaData["id"])[0][0] 
			for x in sim_data.process.translation.monomerData["rnaId"]])
		'''
	def find_overlapping_tus(self, tu_dict):
		'''
		Input: 
			tu_dict: Dictionary of all the transcription units. Keys are the tu_ids. 
				Values are the counts of each TU.
		Output: 
			overlapping_tu_dict: Contains all monocistronic mRNAs found within overlapping
			transcription units (TUS) as keys (str), and all the TUS they are a 
			member of as a value stored in a list.
			This function also adds the locations of the mRNAs which are assumed
			to all be in the cytoplasm - Though this does not always seem to be 
			the case and is dealt with in a later function. 
			Example:
			{
			'EG10526_RNA[c]': [u'EG10527_EG10526_EG10524_RNA[c]', u'EG10526_EG10524_RNA[c]'], 
			'EG10527_RNA[c]': [u'EG10527_EG10526_EG10524_RNA[c]'], 
			'EG12197_RNA[c]': [u'EG12197_RNA[c]', u'EG12197_EG12144_RNA[c]'], 
			'EG12144_RNA[c]': [u'EG12144_RNA[c]', u'EG12197_EG12144_RNA[c]'], 
			'EG10524_RNA[c]': [u'EG10527_EG10526_EG10524_RNA[c]', u'EG10526_EG10524_RNA[c]']
			}
		TODO: 
			-Find out why the locations of the RNA's is not always set to 
			be in the cytoplasm.
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
					if monocistron in tu.split('_'):
						monocistron_loc = monocistron + '_RNA[c]'
						tu_loc = tu + '_RNA[c]'
						overlapping_tu_dict.setdefault(monocistron_loc, [])
						# Want to check if a key is already present as a value to
						# prevent adding it multiple times.
						if tu_loc  not in overlapping_tu_dict[monocistron_loc]:
							overlapping_tu_dict[monocistron_loc].append(tu_loc)
		return overlapping_tu_dict

	def find_overlapping_tu_fractions(self, overlapping_tu_dict, tu_dict):
		'''
		Input:
			overlapping_tu_dict created in the function: find_overlapping_tus
			tu_dict: created in the setup portion of _buildMonomerIndexToRnaMapping
			You can find descriptions of both of these dictionaries in: 
				find_overlapping_tus
		Output:
			overlapping_tu_fractions: a nested dictionary taking the structure of the
			overlapping_tu_dict but attaches counts to each tu. See the following 
			example.
			{
			'EG10526_RNA[c]': {u'EG10526_EG10524_RNA[c]': 0.5436849253059814, 
							   u'EG10527_EG10526_EG10524_RNA[c]': 0.4563150746940185}, 
			'EG10527_RNA[c]': {u'EG10527_EG10526_EG10524_RNA[c]': 1.0}, 
			'EG12197_RNA[c]': {u'EG12197_RNA[c]': 0.43015561634652816, 
			 				   u'EG12197_EG12144_RNA[c]': 0.5698443836534718}, 
			'EG12144_RNA[c]': {u'EG12144_RNA[c]': 0.19687419413435814, 
			 				   u'EG12197_EG12144_RNA[c]': 0.8031258058656418}, 
			'EG10524_RNA[c]': {u'EG10526_EG10524_RNA[c]': 0.5436849253059814, 
			 				   u'EG10527_EG10526_EG10524_RNA[c]': 0.4563150746940185}
			 }

		Note:
			This is matching rna_ids gathered from the monomer data 
			with TU_ids from the TU data. It is not actually matching the monomer
			ID with a TU ID. This is bc each monomer is still attached
			to a single rna_id (even if its never actually transcribed as a mono-
			cistronic unit). This makes it easier to set up this relationship.

			Assumption: if the total counts to a TU are zero, 
			just assign an equal fraction for any overlapping TU's
		'''
		tu_fraction_dict = {}
		for key, value in overlapping_tu_dict.items():
			tu_fraction_dict.setdefault(key, [])
			count_sum = 0
			#first calculate the sum counts of the tus that belong to each monomer.
			for v in value:
				#Need to remove the location tag from v in order to find it in the tu_dict.
				count_sum += tu_dict[v[:-7]]
			for v in value:
				#if the total counts to a TU are zero, just assign an equal fraction for any overlapping TU's
				if count_sum == 0:
					tu_fraction_dict[key].append(1/len(value))
				else:
					tu_fraction_dict[key].append(tu_dict[v[:-7]] / count_sum)

		overlapping_tu_fractions = {}
		for monomer_monocistron, overlapping_tus in overlapping_tu_dict.items():
			for mrna_monocistron, tu_fraction in tu_fraction_dict.items():
				if monomer_monocistron == mrna_monocistron:
					overlapping_tu_fractions[monomer_monocistron] = {
						fraction: tu_fraction[index] 
						for index, fraction in enumerate(overlapping_tus)}

		return overlapping_tu_fractions

	def create_transformation_matrices(self, monomer_id_index_dict, 
		mrna_id_index_dict, overlapping_tu_fractions):
		'''
		Input:
			monomer_id_index_dict: Dictionary of rna_id (key) and its associated
			index from the monomer data.
			mrna_id_index_dict: Dictionary of rna_id (key) and its associated
			index from the transcription unit data.
			overlapping_tu_fractions: Described in find_overlapping_tu_fractions
		Output:
			self.monomer_to_mrna_transform:
				Dimensions:
					rows: length of monomers
					cols: length of the mRNA transcription units
				Mapping between monomers and transcription units that can be 
				used to make the monomers. The mapping is 1:1 if a monomer can
				only be made by a single TU. The mapping is a fraction if a 
				monomer can be made by >1 TU. The rows should always equal to 1.

				There is a check that can be run to ensure this is true:
					check_monomer_to_mrna_transformation_matrix
			self.mrna_to_monomer_transform:
				Dimensions: 
					rows = length of the mRNA transcription units
					cols = length of monomers
				Boolean mapping between the monomers and all the transcription
				units that can be used to make the monomer.
		'''
		# Initialize the matrices with a 1:1 mapping for all monocistronic mRNA's
		# from monomer indexing, to the corresponding monocistronic transcription
		# unit.

		mrna_tu_count = len(self.mrna)
		monomer_count = len(self.monomer)
	
		monomer_to_mrna_transform = np.zeros((monomer_count, mrna_tu_count))
		mrna_to_monomer_transform = np.zeros((mrna_tu_count, monomer_count))

		for monocistron, moncistron_index in monomer_id_index_dict.items():
			for tu, tu_index in mrna_id_index_dict.items():
				if monocistron == tu or monocistron[:-3] + '[None]' == tu:
					monomer_to_mrna_transform[moncistron_index, tu_index] = 1
					mrna_to_monomer_transform[tu_index ,moncistron_index] = 1
		# -For monomer_to_mrna_transform, add the fractions of each overlapping 
		# transcription unit that can be expected to represent each monomer.
		# -For mrna_to_monomer_transform, using the overlapping tu data, add mapping
		# between all transcription units that can be used to make a particular
		# monomer. 
		# 
		for monocistron, transcription_units in overlapping_tu_fractions.items():
			#get the index of the monocistron according to the monomer indexing.
			if monocistron in monomer_id_index_dict:
				monocistron_index = monomer_id_index_dict[monocistron]
			elif monocistron[:-3] + '[None]' in monomer_id_index_dict:
				monocistron_index = monomer_id_index_dict[monocistron[:-3] + '[None]']
			#Get the index of the transcription_units associated to the 
			#mono-cistron according to the mRNA indexing.
			for tu in transcription_units:
				if tu in mrna_id_index_dict:
					fraction = overlapping_tu_fractions[monocistron][tu]
					tu_index = mrna_id_index_dict[tu]
					monomer_to_mrna_transform[
						monocistron_index, tu_index] = fraction
					mrna_to_monomer_transform[tu_index, monocistron_index] = 1
		return 	monomer_to_mrna_transform, mrna_to_monomer_transform

	def check_monomer_to_mrna_transformation_matrix(self, monomer_to_mrna_transform):
		'''
		Input:
			monomer_to_mrna_transform: Described in create_transformation_matrices
		Output:
			Will print to screen all the rows in self.monomer_to_mrna_transform 
			that are greater than or less than 1.
			This is because this is a matrix of the fraction of that a monomer will
			be coded by different transcription units. If it is over or under 1, 
			there is likely an error in how the matrix is being constructed.
			If the number is extremely close to 1 (i.e. 0.9999999), it's likely a 
			rounding 'error'. Do not know yet if this will negatively impact anything.
		'''
		for index, row in enumerate(self.monomer_to_mrna_transform):
				if sum(row) < 1:
					print("Row " + str(index) + " is less than 1: " + str(sum(row)))
				elif sum(row) > 1:
					print("Row " + str(index) + " is greater than 1: " + str(sum(row)))
		return

	def _buildMonomerIndexToRnaMapping(self, raw_data, sim_data):
		'''
		Input:
			raw_data
			sim_data
		Output:
		self.mrna_to_monomer_transform:
			Description in create_transformation_matrices
		self.monomer_to_mrna_transform:
			Description in create_transformation_matrices
		Note:
			If check_matrix = True: 
			Will print to screen all the rows in self.monomer_to_mrna_transform 
			that are greater than or less than 1.
			This is because this is a matrix of the fraction of that a monomer will
			be coded by different transcription units. If it is over or under 1, 
			there is likely an error in how the matrix is being constructed.
			If the number is extremely close to 1 (i.e. 0.9999999), it's likely a 
			rounding 'error'. Do not know yet if this will negatively impact anything.
		'''
		# Read above for description of this boolean.
		check_matrix = True

		# Create dictionaries of mRNA IDs attached to their indexes from both the transcription 
		# unit and monomer space.
		mrna_id_index_dict = {mrna_id: index 
			for index, mrna_id in enumerate(self.mrna['id'])}
		monomer_id_index_dict = {rna_id: index 
			for index, rna_id in enumerate(self.monomer['rna_id'])}
		# tu_dict: Contains all the TU's as keys (str), and their counts as values
		tu_dict = {tu['tu_id']: tu['tu_count'] 
			for tu in raw_data.transcription_units}
		overlapping_tu_dict = self.find_overlapping_tus(tu_dict)
	
		overlapping_tu_fractions = self.find_overlapping_tu_fractions(
			overlapping_tu_dict, tu_dict)

		self.monomer_to_mrna_transform, self.mrna_to_monomer_transform = self.create_transformation_matrices(
			monomer_id_index_dict, mrna_id_index_dict, overlapping_tu_fractions)

		if check_matrix:
			self.check_monomer_to_mrna_transformation_matrix(
				self.monomer_to_mrna_transform)
