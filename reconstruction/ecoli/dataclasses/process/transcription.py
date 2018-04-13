"""
SimulationData for transcription process

@author: Nick Ruggero
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 03/06/2015

TODO - solve optimization for degRates and expression for TU instead of average
"""

from __future__ import division

import numpy as np

from wholecell.utils import units
from wholecell.utils.unit_struct_array import UnitStructArray
from wholecell.utils.polymerize import polymerize

#RNA_SEQ_ANALYSIS = "seal_rpkm"
RNA_SEQ_ANALYSIS = "rsem_tpm"

class Transcription(object):
	""" Transcription """

	def __init__(self, raw_data, sim_data):
		self._buildRnaData(raw_data, sim_data)
		self._buildTuData(raw_data, sim_data)
		self._buildTranscription(raw_data, sim_data)

	def _buildRnaData(self, raw_data, sim_data):
		assert all([len(rna['location']) == 1 for rna in raw_data.rnas])
		rnaIds = ['{}[{}]'.format(rna['id'], rna['location'][0]) for rna in raw_data.rnas if len(rna['location']) == 1]
		rnaDegRates = np.log(2) / np.array([rna['halfLife'] for rna in raw_data.rnas]) # TODO: units
		rnaLens = np.array([len(rna['seq']) for rna in raw_data.rnas])
		ntCounts = np.array([
			(rna['seq'].count('A'), rna['seq'].count('C'),
				rna['seq'].count('G'), rna['seq'].count('U'))
			for rna in raw_data.rnas
			])

		# Load expression from RNA-seq data
		expression = []
		for rna in raw_data.rnas:
			arb_exp = [x[sim_data.basal_expression_condition] for x in eval("raw_data.rna_seq_data.rnaseq_{}_mean".format(RNA_SEQ_ANALYSIS)) if x['Gene'] == rna['geneId']]
			if len(arb_exp):
				expression.append(arb_exp[0])
			elif rna['type'] == 'mRNA' or rna['type'] == 'miscRNA':
				raise Exception('No RNA-seq data found for {}'.format(rna['id']))
			elif rna['type'] == 'rRNA' or rna['type'] == 'tRNA':
				expression.append(0.)
			else:
				raise Exception('Unknown RNA {}'.format(rna['id']))

		expression = np.array(expression)
		synthProb = expression * (
			np.log(2) / sim_data.doubling_time.asNumber(units.s)
			+ rnaDegRates
			)

		synthProb /= synthProb.sum()

		KcatEndoRNase = 0.001
		EstimateEndoRNases = 5000

		Km = (KcatEndoRNase * EstimateEndoRNases / rnaDegRates) - expression

		mws = np.array([rna['mw'] for rna in raw_data.rnas]).sum(axis = 1)

		geneIds = np.array([rna['geneId'] for rna in raw_data.rnas])

		size = len(rnaIds)

		is23S = np.zeros(size, dtype = np.bool)
		is16S = np.zeros(size, dtype = np.bool)
		is5S = np.zeros(size, dtype = np.bool)

		for rnaIndex, rna in enumerate(raw_data.rnas):
			if rna["type"] == "rRNA" and rna["id"].startswith("RRL"):
				is23S[rnaIndex] = True

			if rna["type"] == "rRNA" and rna["id"].startswith("RRS"):
				is16S[rnaIndex] = True

			if rna["type"] == "rRNA" and rna["id"].startswith("RRF"):
				is5S[rnaIndex] = True

		sequences = [rna['seq'] for rna in raw_data.rnas]

		maxSequenceLength = max(len(sequence) for sequence in sequences)

		monomerIds = [x['monomerId'] for x in raw_data.rnas]

		rnaData = np.zeros(
			size,
			dtype = [
				('id', 'a50'),
				('degRate', 'f8'),
				('length', 'i8'),
				('countsACGU', '4i8'),
				('mw', 'f8'),
				('isMRna', 'bool'),
				('isMiscRna', 'bool'),
				('isRRna', 'bool'),
				('isTRna', 'bool'),
				('isRRna23S', 'bool'),
				('isRRna16S', 'bool'),
				('isRRna5S', 'bool'),
				('isRProtein', 'bool'),
				('isRnap',	'bool'),
				('sequence', 'a{}'.format(maxSequenceLength)),
				('geneId', 'a50'),
				('KmEndoRNase', 'f8'),
				]
			)

		rnaData['id'] = rnaIds
		rnaData['degRate'] = rnaDegRates
		rnaData['length'] = rnaLens
		rnaData['countsACGU'] = ntCounts
		rnaData['mw'] = mws
		rnaData['isMRna'] = [rna["type"] == "mRNA" for rna in raw_data.rnas]
		rnaData['isMiscRna'] = [rna["type"] == "miscRNA" for rna in raw_data.rnas]
		rnaData['isRRna'] = [rna["type"] == "rRNA" for rna in raw_data.rnas]
		rnaData['isTRna'] = [rna["type"] == "tRNA" for rna in raw_data.rnas]
		rnaData['isRProtein'] = ["{}[c]".format(x) in sim_data.moleculeGroups.rProteins for x in monomerIds]
		rnaData['isRnap'] = ["{}[c]".format(x) in sim_data.moleculeGroups.rnapIds for x in monomerIds]
		rnaData['isRRna23S'] = is23S
		rnaData['isRRna16S'] = is16S
		rnaData['isRRna5S'] = is5S
		rnaData['sequence'] = sequences
		rnaData['geneId'] = geneIds
		rnaData['KmEndoRNase'] = Km

		field_units = {
			'id'			:	None,
			'degRate'		:	1 / units.s,
			'length'		:	units.nt,
			'countsACGU'	:	units.nt,
			'mw'			:	units.g / units.mol,
			'isMRna'		:	None,
			'isMiscRna'		:	None,
			'isRRna'		:	None,
			'isTRna'		:	None,
			'isRRna23S'		:	None,
			'isRRna16S'		:	None,
			'isRRna5S'		:	None,
			'isRProtein'	:	None,
			'isRnap'		:	None,
			'sequence'		:   None,
			'geneId'		:	None,
			'KmEndoRNase'	:	units.mol / units.L,
			}

		self.rnaExpression = {}
		self.rnaSynthProb = {}

		self.rnaExpression["basal"] = expression / expression.sum()
		self.rnaSynthProb["basal"] = synthProb / synthProb.sum()

		self.rnaData = UnitStructArray(rnaData, field_units)

	def _buildTuData(self, raw_data, sim_data):
		'''
		Loads information about transcription units from raw_data and builds
		other necessary information related to transcription units

		Creates:
			self.rnaExpression (dict {condition (str): expression (np.array)}) -
				dictionary containing expression of each transcript
				(in 'basal' condition)
			self.rnaSynthProb (dict {condition (str): synthProb (np.array)}) -
				dictionary containing synthesis probabilities of each
				transcript (in 'basal' condition)
			self.tuData (UnitStructArray) - contains TU relevant info
			self.tuMappingMatrix (# RNA by # TU np.array) - matrix with entry
				of 1 if an RNA is in an associated TU
		'''

		location = '[c]'  # TODO - store location as part of flat file?
		rnaIndices = {rna['id']: idx for idx, rna in enumerate(self.rnaData)}
		tuIds = ['{}{}'.format(tu['id'], location) for tu in raw_data.transcription_units]
		tuDegRates = np.array([self._geometricMean([
			self.rnaData[rnaIndices[rnaId + location]]['degRate']
			for rnaId in tu['rnas']]) for tu in raw_data.transcription_units])
		tuLens = np.array([tu['length'] for tu in raw_data.transcription_units])
		# TODO - save ntCounts in flat file and load from there
		ntCounts = np.array([
			(tu['seq'].count('A'), tu['seq'].count('C'),
				tu['seq'].count('G'), tu['seq'].count('U'))
			for tu in raw_data.transcription_units
			])

		expression = np.array([self._geometricMean([
			self.rnaExpression['basal'][rnaIndices[rnaId + location]]
			for rnaId in tu['rnas']]) for tu in raw_data.transcription_units])
		synthProb = expression * (
			np.log(2) / sim_data.doubling_time.asNumber(units.s)
			+ tuDegRates
			)
		synthProb /= synthProb.sum()

		KcatEndoRNase = 0.001
		EstimateEndoRNases = 5000

		Km = (KcatEndoRNase * EstimateEndoRNases / tuDegRates) - expression

		mws = np.array([tu['mw'] for tu in raw_data.transcription_units])

		sequences = [tu['seq'] for tu in raw_data.transcription_units]
		maxSequenceLength = max(len(sequence) for sequence in sequences)

		tuData = np.zeros(
			len(tuIds),
			dtype = [
				('id', 'a50'),
				('degRate', 'f8'),
				('length', 'i8'),
				('countsACGU', '4i8'),
				('mw', '{}f8'.format(len(sim_data.molecular_weight_keys))),
				('isMRna', 'bool'),
				('isMiscRna', 'bool'),
				('isRRna', 'bool'),
				('isTRna', 'bool'),
				('isRRna23S', 'bool'),
				('isRRna16S', 'bool'),
				('isRRna5S', 'bool'),
				('isRProtein', 'bool'),
				('isRnap',	'bool'),
				('sequence', 'a{}'.format(maxSequenceLength)),
				('KmEndoRNase', 'f8'),
				('isProcessed', 'bool'),
				]
			)

		tuData['id'] = tuIds
		tuData['degRate'] = tuDegRates
		tuData['length'] = tuLens
		tuData['countsACGU'] = ntCounts
		tuData['mw'] = mws
		# TODO - functionalize array creation or place in a single loop?
		tuData['isMRna'] = [np.any([self.rnaData[rnaIndices[rnaId + location]]['isMRna'] for rnaId in tu['rnas']]) for tu in raw_data.transcription_units]
		tuData['isMiscRna'] = [np.any([self.rnaData[rnaIndices[rnaId + location]]['isMiscRna'] for rnaId in tu['rnas']]) for tu in raw_data.transcription_units]
		tuData['isRRna'] = [np.any([self.rnaData[rnaIndices[rnaId + location]]['isRRna'] for rnaId in tu['rnas']]) for tu in raw_data.transcription_units]
		tuData['isTRna'] = [np.any([self.rnaData[rnaIndices[rnaId + location]]['isTRna'] for rnaId in tu['rnas']]) for tu in raw_data.transcription_units]
		tuData['isRProtein'] = [np.any([self.rnaData[rnaIndices[rnaId + location]]['isRProtein'] for rnaId in tu['rnas']]) for tu in raw_data.transcription_units]
		tuData['isRnap'] = [np.any([self.rnaData[rnaIndices[rnaId + location]]['isRnap'] for rnaId in tu['rnas']]) for tu in raw_data.transcription_units]
		tuData['isRRna23S'] = [np.any([self.rnaData[rnaIndices[rnaId + location]]['isRRna23S'] for rnaId in tu['rnas']]) for tu in raw_data.transcription_units]
		tuData['isRRna16S'] = [np.any([self.rnaData[rnaIndices[rnaId + location]]['isRRna16S'] for rnaId in tu['rnas']]) for tu in raw_data.transcription_units]
		tuData['isRRna5S'] = [np.any([self.rnaData[rnaIndices[rnaId + location]]['isRRna5S'] for rnaId in tu['rnas']]) for tu in raw_data.transcription_units]
		tuData['sequence'] = sequences
		tuData['KmEndoRNase'] = Km
		tuData['isProcessed'] = [tu['processed'] for tu in raw_data.transcription_units]

		field_units = {
			'id'			:	None,
			'degRate'		:	1 / units.s,
			'length'		:	units.nt,
			'countsACGU'	:	units.nt,
			'mw'			:	units.g / units.mol,
			'isMRna'		:	None,
			'isMiscRna'	:	None,
			'isRRna'		:	None,
			'isTRna'		:	None,
			'isRRna23S'	:	None,
			'isRRna16S'	:	None,
			'isRRna5S'		:	None,
			'isRProtein'	:	None,
			'isRnap'		:	None,
			'sequence'		:   None,
			'KmEndoRNase'	:	units.mol / units.L,
			'isProcessed'	:	None,
			}

		self.rnaExpression = {}
		self.rnaSynthProb = {}

		self.rnaExpression["basal"] = expression / expression.sum()
		self.rnaSynthProb["basal"] = synthProb / synthProb.sum()

		self.tuData = UnitStructArray(tuData, field_units)

		tuMappingMatrix = np.zeros((len(self.rnaData), len(self.tuData)))
		for idx, tu in enumerate(raw_data.transcription_units):
			for rna in tu['rnas']:
				tuMappingMatrix[rnaIndices[rna + location], idx] = 1

		self.tuMappingMatrix = tuMappingMatrix

	def _geometricMean(self, values):
		'''
		Returns a float of the geometric mean for a list of values

		Inputs:
			values (list or np.array) - list of values to calculate mean of
		'''

		return np.exp(np.mean(np.log(values)))

	def _buildTranscription(self, raw_data, sim_data):
		sequences = self.tuData["sequence"] # TODO: consider removing sequences

		maxLen = np.int64(
			self.tuData["length"].asNumber().max()
			+ 2 * sim_data.growthRateParameters.rnaPolymeraseElongationRate.asNumber(units.nt / units.s) # hardcode!
			)

		self.transcriptionSequences = np.empty((sequences.shape[0], maxLen), np.int8)
		self.transcriptionSequences.fill(polymerize.PAD_VALUE)

		ntMapping = {ntpId:i for i, ntpId in enumerate(["A", "C", "G", "U"])}

		for i, sequence in enumerate(sequences):
			for j, letter in enumerate(sequence):
				self.transcriptionSequences[i, j] = ntMapping[letter]

		self.transcriptionMonomerWeights = (
			(
				sim_data.getter.getMass(sim_data.moleculeGroups.ntpIds)
				- sim_data.getter.getMass(["PPI[c]"])
				)
			/ raw_data.constants['nAvogadro']
			).asNumber(units.fg)

		self.transcriptionEndWeight = (sim_data.getter.getMass(["PPI[c]"]) / raw_data.constants['nAvogadro']).asNumber(units.fg)

	def getRnaInTranscriptionUnit(self, transcriptionUnitId):
		'''
		Gets the RNA IDs of the RNA in a transcription unit (must include
		compartment id with transcriptionUnit).	Will return for multiple
		transcription units if passed an iterable

		Returns np array of np arrays of RNA ID str in transcriptionUnitId
		Inputs:
			transcriptionUnitId (str or iterable of str) - ID of transcription
				unit(s) to get the RNA ID that are part of that TU
		'''

		if isinstance(transcriptionUnitId, basestring):
			transcriptionUnitId = [transcriptionUnitId]

		rnaIds = []
		for tu in transcriptionUnitId:
			idx = np.where(self.tuData['id'] == tu)[0]
			rnaIds.append(self.rnaData['id'][np.where(self.tuMappingMatrix[:, idx])[0]])

		return np.array(rnaIds)

	def getTranscriptionUnitWithRna(self, rnaId):
		'''
		Gets the TU IDs that the RNA is a part of (must include compartment
		id with rnaId). Will return for multiple RNAs if passed an iterable

		Returns np array of np arrays of transcription unit ID str with rnaId
		Inputs:
			rnaId (str or iterable of str) - ID of RNA(s) to get the
				transcription unit IDs that contain that RNA
		'''

		if isinstance(rnaId, basestring):
			rnaId = [rnaId]

		tuIds = []
		for rna in rnaId:
			idx = np.where(self.rnaData['id'] == rna)[0]
			tuIds.append(self.tuData['id'][np.where(self.tuMappingMatrix[idx, :])[1]])

		return np.array(tuIds)