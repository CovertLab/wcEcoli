#!/usr/bin/env python

"""
PolypeptideProcessing

Sub-model for processing nascent polypeptides into their mature form.

TODO:

@author: Kalli Kappel
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 10/14/14
"""

from __future__ import division

from itertools import izip

import numpy as np

import wholecell.processes.process
from wholecell.utils.polymerize import buildSequences, polymerize, computeMassIncrease, PAD_VALUE
from wholecell.utils.random import stochasticRound
from wholecell.utils import units

class PolypeptideProcessing(wholecell.processes.process.Process):
	""" PolypeptideProcessing """

	_name = "PolypeptideProcessing"

	# Constructor
	def __init__(self):
		# Parameters
		self.proteinSequences = None
		self.proteinSequencesNascent = None
		self.h2oWeight = None
		self.aaWeightsIncorporated = None

		# Views
		self.bulkMonomers = None
		self.bulkMonomersNascent = None
		self.aas = None
		self.h2o = None

		super(PolypeptideProcessing, self).__init__()


	# Construct object graph
	def initialize(self, sim, kb):
		super(PolypeptideProcessing, self).initialize(sim, kb)

		# Load parameters

		proteinIds = kb.proteinData['id']
		proteinNascentIds = kb.proteinNascentData['id']

		self.proteinLengths = kb.proteinData["length"].asNumber()
		self.proteinLengthsNascent = kb.proteinNascentData["length"].asNumber()

		#self.proteinSequences = kb.translationSequences
		self.proteinSequencesNascent = kb.translationSequencesNascent

		self.aaWeightsIncorporated = kb.translationMonomerWeights

		# Views

		self.bulkMonomers = self.bulkMoleculesView(proteinIds)
		self.bulkMonomersNascent = self.bulkMoleculesView(proteinNascentIds)

		self.aas = self.bulkMoleculesView(kb.aaIDs)
		self.h2o = self.bulkMoleculeView('H2O[c]')

		self.pi = self.bulkMoleculeView("PI[c]")
		self.h   = self.bulkMoleculeView("H[c]")

		self.numberAAs = len(kb.aaIDs)
		self.pdf = self.bulkMoleculeView("EG11440-MONOMER[c]")
		self.pdfKcat = kb.pdfKcat
		self.map = self.bulkMoleculeView("EG10570-MONOMER[c]")
		self.mapKcat = kb.mapKcat

	def calculateRequest(self):

		self.pdf.requestAll()
		self.map.requestAll()
		self.bulkMonomersNascent.requestAll()


	# Calculate temporal evolution
	def evolveState(self):
		if sum(self.bulkMonomersNascent.counts()) == 0:
			return

		pdfCapacity = self.pdf.count() * self.pdfKcat
		mapCapacity = self.map.count() * self.mapKcat

		countToBeDeformylated = sum(self.bulkMonomersNascent.counts()[np.where(self.bulkMonomersNascent.counts())])
		metCleavedIndices = np.where(self.proteinLengthsNascent-self.proteinLengths)[0]

		if len(metCleavedIndices) > 0:
			countMetCleaved = sum(self.bulkMonomersNascent.counts()[metCleavedIndices])
			import ipdb; ipdb.set_trace()
			aasVector = np.zeros(self.numberAAs)
			aasVector[12] = countMetCleaved
			self.aas.countsInc(aasVector)
		else: countMetCleaved = 0

		# Logic needed if we take enzyme capactiy into consideration
		#if pdfCapacity > countToBeDeformylated and mapCapacity > countMetCleaved:
		#	self.aas.countsInc(aasVector)
		#	self.bulkMonomers.countsInc(self.bulkMonomersNascent.counts())
		#	self.bulkMonomersNascent.countsDec(self.bulkMonomersNascent.counts())

		self.bulkMonomers.countsInc(self.bulkMonomersNascent.counts())
		self.bulkMonomersNascent.countsDec(self.bulkMonomersNascent.counts())
