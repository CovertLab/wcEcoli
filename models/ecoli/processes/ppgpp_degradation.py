#!/usr/bin/env python

"""
ppGppDegradation

ppGpp degradation submodel

@author: Nick Ruggero
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 11/04/2014 - the birthday process
"""

from __future__ import division

import numpy as np

import wholecell.processes.process
from wholecell.utils import units

class ppGppDegradation(wholecell.processes.process.Process):
	""" ppGppDegradation """

	_name = "ppGppDegradation"

	# Constructor
	def __init__(self):
		# Parameters
		self.ppGpp = None

		super(ppGppDegradation, self).__init__()


	# Construct object graph
	def initialize(self, sim, kb):
		super(ppGppDegradation, self).initialize(sim, kb)

		self.k_m = 800. * units.mmol / units.L
		self.k_cat = 3. * 1 / units.s
		self.v_max = self.k_cat * 58

		self.nAvogadro = kb.nAvogadro
		self.cellDensity = kb.cellDensity

		# Views
		self.ppGpp = self.bulkMoleculeView("PPGPP[c]")

	def calculateRequest(self):
		cellMass = (self.readFromListener("Mass", "cellMass") * units.fg)
		cellVolume = cellMass / self.cellDensity
		ppGpp_conc = (1 / self.nAvogadro) * (1 / cellVolume) * self.ppGpp.total()[0]
		maxTurnover = (self.v_max * (ppGpp_conc / (ppGpp_conc + self.k_m))).asNumber(1 / units.s)  * self.timeStepSec
		maxTurnover = np.floor(maxTurnover)
		
		print 'ppGpp requested for degradation: {}'.format(maxTurnover)
		#import ipdb; ipdb.set_trace()
		self.ppGpp.requestIs(maxTurnover)

	# Calculate temporal evolution
	def evolveState(self):
		cellMass = (self.readFromListener("Mass", "cellMass") * units.fg)
		cellVolume = cellMass / self.cellDensity
		ppGpp_conc = (1 / self.nAvogadro) * (1 / cellVolume) * self.ppGpp.total()[0]
		maxTurnover = (self.v_max * (ppGpp_conc / (ppGpp_conc + self.k_m))).asNumber(1 / units.s)  * self.timeStepSec
		maxTurnover = np.floor(maxTurnover)

		if maxTurnover > self.ppGpp.count():
			print 'Allocated less ppGpp than requested to degradation process!'

		self.ppGpp.countDec(maxTurnover)