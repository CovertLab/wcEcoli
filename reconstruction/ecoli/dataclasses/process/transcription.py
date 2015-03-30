"""
SimulationData for transcription process

@author: Nick Ruggero
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 03/06/2015
"""

from __future__ import division

from wholecell.utils import units
from wholecell.utils.unit_struct_array import UnitStructArray
import numpy as np

class Transcription(object):
	""" Transcription """

	def __init__(self, raw_data, sim_data):
		self._buildRnaData(raw_data, sim_data)
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
		expression = np.array([rna['expression'] for rna in raw_data.rnas])
		synthProb = expression * (
			np.log(2) / raw_data.parameters['cellCycleLen'].asNumber(units.s)
			+ rnaDegRates
			)
		
		synthProb /= synthProb.sum()

		expression = [x['expression'] for x in raw_data.rnas]

		##################################
		### Vary gene copy number code ###
		synthProbPopAvg = expression * (
			np.log(2) / sim_data.constants.cellCycleLen.asNumber(units.s)
			+ rnaDegRates
			)

		synthProbPopAvg /= synthProbPopAvg.sum()

		isRRna = np.array([rna["type"] == "rRNA" for rna in raw_data.rnas])
		proteinIds = np.array([x['id'] for x in raw_data.proteins])
		proteinRnaIds = np.array([x['rnaId'] for x in raw_data.proteins])
		rProteins50 = proteinRnaIds[np.array([np.where(proteinIds == x[:-3])[0][0] for x in sim_data.moleculeGroups.s50_proteins if len(np.where(proteinIds == x[:-3])[0])])]
		rProteins30 = proteinRnaIds[np.array([np.where(proteinIds == x[:-3])[0][0] for x in sim_data.moleculeGroups.s30_proteins if len(np.where(proteinIds == x[:-3])[0])])]
		#rProteinsIdxs = np.array([np.where(rnaIds==k)[0][0] for k in rProteins if len(np.where(rnaIds==k)[0])])
		#rProteinsBool = np.zeros(len(rnaIds),dtype=bool)
		#rProteinsBool[rProteinsIdxs]=np.ones(len(rProteinsIdxs),dtype=bool)
		rProteins50Bool = np.in1d(rnaIds,rProteins50)
		rProteins30Bool = np.in1d(rnaIds,rProteins30)

		#Ribosome associated RNAs (rRNAs and mRNAs for rProteins)
		highlyRegulated = (isRRna | rProteins50Bool | rProteins30Bool)

		fractionSynthProbHighlyRegulated = np.sum(synthProbPopAvg[highlyRegulated])
		fractionSynthProbNotHighlyRegulated = 1 - fractionSynthProbHighlyRegulated

		geneCoordinates = np.array([gene['coordinate'] for rna in raw_data.rnas for gene in raw_data.genes if rna['geneId'] == gene['id']])
		genomeLength = len(raw_data.genome_sequence)
		geneEndCoordinates = np.array([(x['coordinate'] + x['length']) % genomeLength if x['direction'] == '+' else (x['coordinate'] - x['length']) % genomeLength for x in raw_data.genes])
		

		minDistFromOriC = np.minimum(np.abs(sim_data.constants.oriCCenter.asNumber()-geneEndCoordinates-genomeLength),
						np.abs(geneEndCoordinates-sim_data.constants.oriCCenter.asNumber()))
		ageReplicated = minDistFromOriC / sim_data.constants.dnaPolymeraseElongationRate.asNumber()

		synthProb = np.zeros(len(synthProbPopAvg))
		synthProb[highlyRegulated] = synthProbPopAvg[highlyRegulated]
		synthProb[~highlyRegulated] = synthProbPopAvg[~highlyRegulated] / (2 * np.exp(-np.log(2)*ageReplicated[~highlyRegulated]/sim_data.constants.cellCycleLen.asNumber()))
		synthProb[~highlyRegulated] = synthProb[~highlyRegulated] / np.sum(synthProb[~highlyRegulated]) * fractionSynthProbNotHighlyRegulated

		#synthProb = synthProbPopAvg / (2 * np.exp(-np.log(2)*ageReplicated/sim_data.constants.cellCycleLen.asNumber()))
		#synthProb /= synthProb.sum()
		
		synthProbTimeAvg = np.zeros(len(synthProbPopAvg))
		synthProbTimeAvg[highlyRegulated] = synthProb[highlyRegulated]
		synthProbTimeAvg[~highlyRegulated] = synthProb[~highlyRegulated] * timeAvgProbAsFuncofCopyAge(ageReplicated[~highlyRegulated])
		synthProbTimeAvg[~highlyRegulated] = synthProbPopAvg[~highlyRegulated] / np.sum(synthProbPopAvg[~highlyRegulated]) * fractionSynthProbNotHighlyRegulated

		#synthProbTimeAvg = synthProb * timeAvgProbAsFuncofCopyAge(ageReplicated)
		#synthProbTimeAvg /= synthProbTimeAvg.sum()

		### Vary gene copy number code ###
		##################################

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

		# TODO: Add units
		rnaData = np.zeros(
			size,
			dtype = [
				('id', 'a50'),
				('synthProb', 'f8'),
				('synthProbTimeAvg', 'f8'),
				('ageReplicated', 'f8'),
				('expression', 'float64'),
				('expressionInitial',	'float64'),
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
				('isHighlyRegulated', 'bool'),
				('sequence', 'a{}'.format(maxSequenceLength)),
				('geneId', 'a50')
				]
			)

		rnaData['id'] = rnaIds
		rnaData['synthProb'] = synthProb
		rnaData['synthProbTimeAvg'] = synthProbTimeAvg
		rnaData['ageReplicated'] = ageReplicated
		rnaData['expression'] = expression
		rnaData['expressionInitial']	= np.ones(len(rnaIds)) #This will be set after fitting
		rnaData['degRate'] = rnaDegRates
		rnaData['length'] = rnaLens
		rnaData['countsACGU'] = ntCounts
		rnaData['mw'] = mws
		rnaData['isMRna'] = [rna["type"] == "mRNA" for rna in raw_data.rnas]
		rnaData['isMiscRna'] = [rna["type"] == "miscRNA" for rna in raw_data.rnas]
		rnaData['isRRna'] = [rna["type"] == "rRNA" for rna in raw_data.rnas]
		rnaData['isTRna'] = [rna["type"] == "tRNA" for rna in raw_data.rnas]
		rnaData['isRRna23S'] = is23S
		rnaData['isRRna16S'] = is16S
		rnaData['isRRna5S'] = is5S
		rnaData['isHighlyRegulated'] = highlyRegulated
		rnaData['sequence'] = sequences
		rnaData['geneId'] = geneIds

		field_units = {
			'id'		:	None,
			'synthProb' :	None,
			'synthProbTimeAvg' :	None,
			'ageReplicated' :	units.s,
			'expression':	None,
			'expressionInitial':	None,
			'degRate'	:	1 / units.s,
			'length'	:	units.nt,
			'countsACGU':	units.nt,
			'mw'		:	units.g / units.mol,
			'isMRna'	:	None,
			'isMiscRna'	:	None,
			'isRRna'	:	None,
			'isTRna'	:	None,
			'isRRna23S'	:	None,
			'isRRna16S'	:	None,
			'isRRna5S'	:	None,
			'isHighlyRegulated'	:	None,
			'sequence'  :   None,
			'geneId'	:	None,
			}


		self.rnaData = UnitStructArray(rnaData, field_units)
		#self.getTrnaAbundanceData = getTrnaAbundanceAtGrowthRate

	def _buildTranscription(self, raw_data, sim_data):
		from wholecell.utils.polymerize import PAD_VALUE

		sequences = self.rnaData["sequence"] # TODO: consider removing sequences

		maxLen = np.int64(
			self.rnaData["length"].asNumber().max()
			+ raw_data.parameters['rnaPolymeraseElongationRate'].asNumber(units.nt / units.s)
			)

		self.transcriptionSequences = np.empty((sequences.shape[0], maxLen), np.int8)
		self.transcriptionSequences.fill(PAD_VALUE)

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

def timeAvgProbAsFuncofCopyAge(ageReplicatedInSeconds):
	ageReplicatedInMinutes = ageReplicatedInSeconds/60
	return -8.022e-7*ageReplicatedInMinutes**3+1.179e-4*ageReplicatedInMinutes**2-1.434e-2*ageReplicatedInMinutes+1.215
