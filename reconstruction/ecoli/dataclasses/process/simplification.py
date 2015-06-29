
from __future__ import division

import numpy as np

class Simplification(object):
	def __init__(self, raw_data, sim_data):
		# Build the abstractions needed for simplification

		pRevs = []

		# Remove complexes that are currently not simulated
		FORBIDDEN_MOLECULES = {
			"modified-charged-selC-tRNA", # molecule does not exist
			}

		deleteReactions = []
		for reactionIndex, reaction in enumerate(raw_data.complexationReactions):
			for molecule in reaction["stoichiometry"]:
				if molecule["molecule"] in FORBIDDEN_MOLECULES:
					deleteReactions.append(reactionIndex)
					break

		for reactionIndex in deleteReactions[::-1]:
			del raw_data.complexationReactions[reactionIndex]

		for reaction in raw_data.complexationReactions:
			assert reaction["process"] == "complexation"
			assert reaction["dir"] == 1

			pRevs.append(reaction["p_rev"])


		self.pRevs = np.array(pRevs)
