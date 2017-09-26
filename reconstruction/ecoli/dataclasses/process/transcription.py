"""
SimulationData for transcription process

@author: Nick Ruggero
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 03/06/2015

TODO: add mapping of tRNA to charged tRNA if allowing more than one modified form of tRNA and separate mappings for tRNA and charged tRNA to AA
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
		self._buildTranscription(raw_data, sim_data)
		self._buildChargedTrna(raw_data, sim_data)

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
				raise Exception('Unknonw RNA {}'.format(rna['id']))

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

		# TODO: Add units
		rnaData = np.zeros(
			size,
			dtype = [
				('id', 'a50'),
				# ('synthProb', 'f8'),
				# ('expression', 'float64'),
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
		# rnaData["synthProb"] = synthProb
		# rnaData["expression"] = expression
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
			# 'synthProb' 	:	None,
			# 'expression'	:	None,
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
		#self.getTrnaAbundanceData = getTrnaAbundanceAtGrowthRate

	def _buildTranscription(self, raw_data, sim_data):
		sequences = self.rnaData["sequence"] # TODO: consider removing sequences

		maxLen = np.int64(
			self.rnaData["length"].asNumber().max()
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

	def _buildChargedTrna(self, raw_data, sim_data):
		# create list of charged tRNAs
		trnaNames = self.rnaData['id'][self.rnaData['isTRna']]
		chargedTrnas = [x['modifiedForms'] for x in raw_data.rnas if x['id'] + '[c]' in trnaNames]
		filteredCharged = []
		for chargedArray in chargedTrnas:
			for charged in chargedArray:
				if 'FMET' in charged or 'modified' in charged:
					continue

				assert('c' in sim_data.getter.getLocation([charged])[0])
				filteredCharged += [charged + '[c]']

		self.chargedTrnaNames = filteredCharged
		assert(len(self.chargedTrnaNames) == len(trnaNames))

		# create mapping of each tRNA/charged tRNA to associated AA
		trnaDict = {
			'RNA0-300[c]': 'VAL',
			'RNA0-301[c]': 'LYS',
			'RNA0-302[c]': 'LYS',
			'RNA0-303[c]': 'LYS',
			'RNA0-304[c]': 'ASN',
			'RNA0-305[c]': 'ILE',
			'RNA0-306[c]': 'MET',
			}
		aaNames = sim_data.moleculeGroups.aaIDs
		self.aaFromTrna = np.zeros((len(aaNames), len(trnaNames)))
		for trna in trnaNames:
			aa = trna[:3].upper()
			if aa == 'ALA':
				aa = 'L-ALPHA-ALANINE'
			elif aa == 'ASP':
				aa = 'L-ASPARTATE'
			elif aa == 'SEL':
				aa = 'L-SELENOCYSTEINE'
			elif aa == 'RNA':
				aa = trnaDict[trna]

			assert('c' in sim_data.getter.getLocation([aa])[0])
			aa += '[c]'
			if aa in aaNames:
				aaIdx = aaNames.index(aa)
				trnaIdx = np.where(trnaNames == trna)[0]
				self.aaFromTrna[aaIdx, trnaIdx] = 1

		# array for stoichiometry matrix
		molecules = []

		stoichMatrixI = []
		stoichMatrixJ = []
		stoichMatrixV = []

		stoichMatrixMass = []

		# remove reactions from modificationReactions that don't have both an uncharged and charged tRNA
		deleteReactions = []
		for reactionIndex, reaction in enumerate(raw_data.modificationReactions):
			noChargedTrnaInReaction = True
			noTrnaInReaction = True
			for mol in [molecule['molecule'] + '[' + molecule['location'] + ']' for molecule in reaction['stoichiometry']]:
				if mol in self.chargedTrnaNames:
					noChargedTrnaInReaction = False

				if mol in trnaNames:
					noTrnaInReaction = False

			if noChargedTrnaInReaction or noTrnaInReaction:
				deleteReactions.append(reactionIndex)

		for reactionIndex in deleteReactions[::-1]:
			del raw_data.modificationReactions[reactionIndex]

		# create stoichiometry matrix for charging reactions
		for reaction in raw_data.modificationReactions:
			assert reaction['process'] == 'rna'
			assert reaction['dir'] == 1

			trna = None
			for mol in [molecule['molecule'] + '[' + molecule['location'] + ']' for molecule in reaction['stoichiometry']]:
				if mol in trnaNames:
					trna = mol
					break

			if trna is None:
				continue
			trnaIndex = np.where(trnaNames == trna)[0][0]

			for molecule in reaction['stoichiometry']:
				if molecule['type'] == 'metabolite':
					moleculeName = '{}[{}]'.format(
						molecule['molecule'].upper(),
						molecule['location']
						)
				else:
					moleculeName = '{}[{}]'.format(
						molecule['molecule'],
						molecule['location']
						)

				if moleculeName not in molecules:
					molecules.append(moleculeName)
					moleculeIndex = len(molecules) - 1
				else:
					moleculeIndex = molecules.index(moleculeName)

				coefficient = molecule['coeff']

				assert coefficient % 1 == 0

				stoichMatrixI.append(moleculeIndex)
				stoichMatrixJ.append(trnaIndex)
				stoichMatrixV.append(coefficient)

		self._stoichMatrixI = np.array(stoichMatrixI)
		self._stoichMatrixJ = np.array(stoichMatrixJ)
		self._stoichMatrixV = np.array(stoichMatrixV)

		self.chargingMolecules = molecules

	# returns matrix with rows of metabolites for each tRNA charging reaction on the column
	def chargingStoichMatrix(self):
		shape = (self._stoichMatrixI.max()+1, self._stoichMatrixJ.max()+1)

		out = np.zeros(shape, np.float64)

		out[self._stoichMatrixI, self._stoichMatrixJ] = self._stoichMatrixV

		return out
