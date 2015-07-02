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

		self.k_m = 800. * units.umol / units.L
		self.k_cat = 3. * 1 / units.s
		self.v_max = self.k_cat * 10000

		self.nAvogadro = kb.constants.nAvogadro
		self.cellDensity = kb.constants.cellDensity

		self.ppGpp_base_conc = kb.process.metabolism.metabolitePoolConcentrations[kb.process.metabolism.metabolitePoolIDs.index("GUANOSINE-5DP-3DP[c]")]

		# Views
		self.ppGpp = self.bulkMoleculeView("GUANOSINE-5DP-3DP[c]")
		self.gdp = self.bulkMoleculeView("GDP[c]")
		self.h2o = self.bulkMoleculeView("WATER[c]")
		self.ppi = self.bulkMoleculeView("PPI[c]")

	def calculateRequest(self):
		cellMass = (self.readFromListener("Mass", "cellMass") * units.fg)
		cellVolume = cellMass / self.cellDensity
		ppGpp_conc = (1 / self.nAvogadro) * (1 / cellVolume) * self.ppGpp.total()[0]
		if ppGpp_conc.asNumber(units.mol / units.L) > self.ppGpp_base_conc.asNumber(units.mol / units.L):
			spoT_saturation = (ppGpp_conc / (ppGpp_conc + self.k_m)).normalize()
			spoT_saturation.checkNoUnit()
			spoT_saturation = spoT_saturation.asNumber()
			maxTurnover = (self.v_max * spoT_saturation).asNumber(1 / units.s)  * self.timeStepSec
			maxTurnover = np.floor(maxTurnover)
		else:
			maxTurnover = 0.
			spoT_saturation = (ppGpp_conc / (ppGpp_conc + self.k_m)).normalize()
			spoT_saturation.checkNoUnit()
			spoT_saturation = spoT_saturation.asNumber()

		self.ppGpp.requestIs(maxTurnover)
		self.h2o.requestIs(maxTurnover)

	# Calculate temporal evolution
	def evolveState(self):
		cellMass = (self.readFromListener("Mass", "cellMass") * units.fg)
		cellVolume = cellMass / self.cellDensity
		ppGpp_conc = (1 / self.nAvogadro) * (1 / cellVolume) * self.ppGpp.total()[0]
		if ppGpp_conc.asNumber(units.mol / units.L) > self.ppGpp_base_conc.asNumber(units.mol / units.L):
			spoT_saturation = (ppGpp_conc / (ppGpp_conc + self.k_m)).normalize()
			spoT_saturation.checkNoUnit()
			spoT_saturation = spoT_saturation.asNumber()
			maxTurnover = (self.v_max * spoT_saturation).asNumber(1 / units.s)  * self.timeStepSec
			maxTurnover = np.floor(maxTurnover)
		else:
			maxTurnover = 0.
			spoT_saturation = (ppGpp_conc / (ppGpp_conc + self.k_m)).normalize()
			spoT_saturation.checkNoUnit()
			spoT_saturation = spoT_saturation.asNumber()

		if maxTurnover > self.ppGpp.count():
			print 'Allocated less ppGpp than requested to degradation process!'

		self.ppGpp.countDec(maxTurnover)
		self.h2o.countDec(maxTurnover)

		self.gdp.countInc(maxTurnover)
		self.ppi.countInc(maxTurnover)

		self.writeToListener("GrowthRateControl", "spoT_saturation", spoT_saturation)
