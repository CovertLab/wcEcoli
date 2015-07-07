
from __future__ import division

import numpy as np

from reconstruction.ecoli.dataclasses.process.complexation import Complexation

class Simplification(object):
	def __init__(self, raw_data, sim_data):
		# Build the abstractions needed for simplification

		complexation = Complexation(raw_data, sim_data)

		pDissoc = []

		moleculeIds = []

		stoichMatrixI = []
		stoichMatrixJ = []
		stoichMatrixV = []

		# Remove complexes that are currently not simulated
		FORBIDDEN_MOLECULES = {
			"modified-charged-selC-tRNA", # molecule does not exist
			}

		complexIds = []
		for proteinComplex in raw_data.proteinComplexes:
			if proteinComplex["id"] in FORBIDDEN_MOLECULES:
				continue
			for loc in proteinComplex["location"]:
				complexIds.append(
					"%s[%s]" % (proteinComplex["id"].encode("utf8"), loc.encode("utf8"))
					)
			pDissoc.append(proteinComplex["p_dissoc"])

		for complexIdx, complexId in enumerate(complexIds):
			composition = complexation.getMonomers(complexId)
			subunitIds = composition["subunitIds"].tolist()
			stoichiometry = composition["subunitStoich"]

			for subunitId in subunitIds:
				if subunitId not in moleculeIds:
					moleculeIds.append(subunitId)
			if complexId in moleculeIds:
				raise Exception, "Duplicate complex entered: %s" % complexId
			moleculeIds.append(complexId)

			for coeff, subunitId in zip(stoichiometry, subunitIds):
				stoichMatrixI.append(moleculeIds.index(subunitId))
				stoichMatrixJ.append(complexIdx)
				stoichMatrixV.append(coeff)
			stoichMatrixI.append(moleculeIds.index(complexId))
			stoichMatrixJ.append(complexIdx)
			stoichMatrixV.append(-1)

		self.moleculeIds = moleculeIds
		self.complexIds = complexIds

		shape = (np.array(stoichMatrixI).max() + 1, np.array(stoichMatrixJ).max() + 1)
		self.stoichMatrix = np.zeros(shape, np.float64)
		self.stoichMatrix[stoichMatrixI, stoichMatrixJ] = stoichMatrixV

		self.pDissoc = np.array(pDissoc)