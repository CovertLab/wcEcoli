# PURPOSE
# modulate transcription via binding of AT complex to DNA
# cleave mRNAs - 16S

from __future__ import division

import numpy as np
import scipy.sparse

import wholecell.processes.process
from wholecell.utils import units

from itertools import izip

class ToxinAntitoxin(wholecell.processes.process.Process):

	_name = 'ToxinAntitoxin'

	def __init__(self):
		super(ToxinAntitoxin, self).__init__()

	def initialize(self, sim, sim_data):
		super(TranscriptInitiation, self).initialize(sim, sim_data)

		# get views of things you need in here Toxin, antitoxin, rna etc views = state
		# pull any info you need out of sim data
		self.mazE = self.bulkMoleculeView('EG10571-MONOMER[c]')
		self.mazEdimer = self.bulkMoleculeView('CPLX0-841[c]')
		self.mazF = self.bulkMoleculeView('EG11249-MONOMER[c]')
		self.mazFdimer = self.bulkMoleculeView('CPLX0-1241[c]')
		self.mazEFcomplex = self.bulkMoleculeView('CPLX0-1242[c]')
		# free promoter?
		# 16S rRNA - in the subunit or in the full ribosome
		self.ribosome30S = self.bulkMoleculeView(sim_data.moleculeIds.s30_fullComplex)
		self.ribosome70S = self.uniqueMoleculesView('activeRibosome')




	def calculateRequest(self):
		# declare what is wanted by the process based on views

		# how many 16S rRNAs will be cleaved? free toxin * rate * free mRNA?

		# how many operators will be bound? free operator * free complex * rate
		pass

	def evolveState(self):
		# take what is allocated and operate on it

		# complexation reactions?

		# binding of complex to promoter

		# cleaving of 16S ACA by mazF

		pass


