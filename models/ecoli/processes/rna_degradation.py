#!/usr/bin/env python

"""
RnaDegradation
<<<<<<< HEAD

RNA degradation sub-model. Encodes molecular simulation of RNA degradation as a Poisson process

TODO:
- handle complexes

@author: Derek Macklin
=======
RNA degradation sub-model. 
Mathematical formulation:
dr/dt = kb - kcatEndoRNase * EndoRNase * r / (Km + r) = r * ln(2)/tau
	where	r = RNA counts
			kb = RNAP synthesis rate 
			tau = doubling time
			kcatEndoRNase = enzymatic activity for EndoRNases
			kd = RNA degradation rates 
			Km = Michaelis-Menten constants fitted to recapitulate first-order RNA decay ( kd * r = kcatEndoRNase * EndoRNase * r / (Km + r) )
This sub-model encodes molecular simulation of RNA degradation as two main steps guided by RNases, "endonucleolytic cleavage" and "exonucleolytic digestion":
1. Compute total counts of RNA to be degraded (D) and total capacity for endo-cleavage (C) at each time point
2. Evaluate C and D. If C > D, then define a fraction of active endoRNases 
3. Dissect RNA degraded into different species (mRNA, tRNA, and rRNA) by accounting endoRNases specificity
4. Update RNA fragments (assumption: fragments are represented as a pull of nucleotides) because of endonucleolytic cleavage
5. Compute total capacity of exoRNases and determine fraction of nucleotides that can be diggested
6. Update pull of metabolites (H and H2O) because of exonucleolytic digestion
@author: Javier Carrera
>>>>>>> master
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 4/2/2013
"""

from __future__ import division

import numpy as np

import wholecell.processes.process
from wholecell.utils.constants import REQUEST_PRIORITY_DEGRADATION
from wholecell.utils import units

class RnaDegradation(wholecell.processes.process.Process):
	""" RnaDegradation """

	_name = "RnaDegradation"

	# Constructor
	def __init__(self):
		# Parameters
		self.rnaLens = None			# RNA lengths
		self.rnaDegRates = None		# RNA degradation rates (1/s)
		self.rnaDegSMat = None		# RNA degradation stoichiometry matrix [metabolite x rna]

		# Views
		self.metabolites = None
		self.rnas = None
		self.rnase = None

		super(RnaDegradation, self).__init__()

	# Construct object graph
	def initialize(self, sim, kb):
		super(RnaDegradation, self).initialize(sim, kb)

		metaboliteIds = ["AMP[c]", "CMP[c]", "GMP[c]", "UMP[c]",
			"WATER[c]", "PPI[c]", "PROTON[c]"]

		nmpIdxs = np.arange(0, 4)
		h2oIdx = metaboliteIds.index('WATER[c]')
		ppiIdx = metaboliteIds.index('PPI[c]')
		hIdx = metaboliteIds.index('PROTON[c]')

		rnaIds = kb.process.transcription.rnaData['id']

		# Rna
		self.rnaDegRates = kb.process.transcription.rnaData['degRate'].asNumber(1 / units.s) * self.timeStepSec
		self.rnaLens = kb.process.transcription.rnaData['length'].asNumber()

		# Build stoichiometric matrix
		# TODO: account for NTP on 5' end
		self.rnaDegSMat = np.zeros((len(metaboliteIds), len(rnaIds)), np.int64)
		self.rnaDegSMat[nmpIdxs, :] = units.transpose(kb.process.transcription.rnaData['countsACGU']).asNumber()
		# self.rnaDegSMat[h2oIdx, :]  = -(self.rnaLens - 1)
		self.rnaDegSMat[h2oIdx, :]  = -self.rnaLens # using one additional water to hydrolyze PPI on 5' end
		self.rnaDegSMat[ppiIdx, :]    =  1
		self.rnaDegSMat[hIdx, :] = self.rnaLens

		# Views
		self.metabolites = self.bulkMoleculesView(metaboliteIds)
		self.rnas = self.bulkMoleculesView(rnaIds)
		self.rnase = self.bulkMoleculeView('EG11259-MONOMER[c]')

		self.bulkMoleculesRequestPriorityIs(REQUEST_PRIORITY_DEGRADATION)

		self.Km = kb.process.transcription.rnaData["KmEndoRNase"]


	# Calculate temporal evolution

	def calculateRequest(self):
		nRNAsToDegrade = np.fmin(
			self.randomState.poisson(self.rnaDegRates * self.rnas.total()),
			self.rnas.total()
			)

		# nReactions = np.dot(self.rnaLens, nRNAsToDegrade)

		# self.h2o.requestIs(nReactions)
		self.rnas.requestIs(nRNAsToDegrade)
		self.rnase.requestAll()

		metaboliteUsage = np.fmax(
			-np.dot(self.rnaDegSMat, nRNAsToDegrade),
			0
			)

		self.metabolites.requestIs(metaboliteUsage)


	def evolveState(self):
		# Check if RNAse R expressed
		if self.rnase.count() == 0:
			return

		# Degrade RNA
		self.metabolites.countsInc(np.dot(
			self.rnaDegSMat,
			self.rnas.counts()
			))

		self.rnas.countsIs(0)
