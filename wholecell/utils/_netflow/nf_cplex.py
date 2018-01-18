
from __future__ import division

from collections import defaultdict

import cplex
import numpy as np
from scipy.sparse import coo_matrix

from ._base import NetworkFlowProblemBase

'''
TODO:
- change naming convention (set, get)
- handle QP vs LP
'''

class NetworkFlowCPLEX(NetworkFlowProblemBase):
	_lowerBoundDefault = 0
	_upperBoundDefault = cplex.infinity

	def __init__(self):
		self._model = cplex.Cplex()

		self._model.set_log_stream(None)
		self._model.set_error_stream(None)
		self._model.set_warning_stream(None)
		self._model.set_results_stream(None)

		self._flows = {}
		self._lb = {}
		self._ub = {}
		self._objective = {}
		self._materialCoeffs = defaultdict(list)

		self._linear = True

		self._eqConstBuilt = False
		self._solved = False

	def _getVar(self, flow):
		if flow in self._flows:
			idx = self._flows[flow]
		else:
			self._model.variables.add(obj=[0])
			idx = len(self._flows)
			self._lb[flow] = self._lowerBoundDefault
			self._ub[flow] = self._upperBoundDefault
			self._model.variables.set_lower_bounds(idx, self._lb[flow])
			self._model.variables.set_upper_bounds(idx, self._ub[flow])
			self._flows[flow] = idx

		return idx

	def flowMaterialCoeffIs(self, flow, material, coefficient):
		if self._eqConstBuilt:
			if material not in self._materialIdxLookup:
				raise Exception("Invalid material")
			if flow not in self._flows:
				raise Exception("Invalid flow")
			materialIdx = self._materialIdxLookup[material]
			flowIdx = self._flows[flow]
			coeffs, flowIdxs = zip(*self._materialCoeffs[material])
			coeffs = list(coeffs)
			flowLoc = flowIdxs.index(flowIdx)
			coeffs[flowLoc] = coefficient
			self._materialCoeffs[material] = zip(coeffs, flowIdxs)

			self._model.linear_constraints.set_coefficients(zip([materialIdx], [flowIdx], [coefficient]))
		else:
			idx = self._getVar(flow)
			self._materialCoeffs[material].append((coefficient, idx))

		self._solved = False

	def flowLowerBoundIs(self, flow, lowerBound):
		idx = self._getVar(flow)
		self._lb[flow] = lowerBound
		self._model.variables.set_lower_bounds(idx, self._lb[flow])

		self._solved = False

	def flowLowerBound(self, flow):
		return self._lb[flow]

	def flowUpperBoundIs(self, flow, upperBound):
		idx = self._getVar(flow)
		self._ub[flow] = upperBound
		self._model.variables.set_upper_bounds(idx, self._ub[flow])

		self._solved = False

	def flowUpperBound(self, flow):
		return self._ub[flow]

	def flowObjectiveCoeffIs(self, flow, coefficient):
		idx = self._getVar(flow)
		self._objective[flow] = coefficient
		if self._linear:
			self._model.objective.set_linear(
				idx,
				coefficient
				)
		else:
			self._model.objective.set_quadratic_coefficients(
				idx,
				idx,
				coefficient
				)

		self._solved = False

	def flowRates(self, flows):
		if isinstance(flows, basestring):
			flows = (flows,)

		self._solve()

		return np.array(
			[self._model.solution.get_values(self._flows[flow]) if flow in self._flows else None for flow in flows]
			)

	# TODO: implement
	def rowDualValues(self, materials):
		# reduced cost of row
		if not self._eqConstBuilt:
			raise Exception("Equality constraints not yet built. Finish construction of the problem before accessing dual values.")

		self._solve()

		return np.array(
			[self._model.solution.get_dual_values(self._materialIdxLookup[material]) if material in self._materialIdxLookup else None for material in materials]
			)
		return

	# TODO: implement
	def columnDualValues(self, fluxNames):
		# reduced cost of col
		if not self._eqConstBuilt:
			raise Exception("Equality constraints not yet built. Finish construction of the problem before accessing dual values.")

		self._solve()

		return np.array(
			[self._model.solution.get_reduced_costs(self._flows[fluxName]) if fluxName in self._flows else None for fluxName in fluxNames]
			)
		return

	def objectiveValue(self):
		return self._model.solution.get_objective_value()

	def getSMatrix(self):
		if not self._eqConstBuilt:
			raise Exception("Equality constraints not yet built. Finish construction of the problem before accessing S matrix.")
		A = np.zeros((len(self._materialCoeffs), len(self._flows)))
		self._materialIdxLookup = {}
		for materialIdx, (material, pairs) in enumerate(sorted(self._materialCoeffs.viewitems())):
			self._materialIdxLookup[material] = materialIdx
			for pair in pairs:
				A[materialIdx, pair[1]] = pair[0]
		return A

	def getFlowNames(self):
		if not self._eqConstBuilt:
			raise Exception("Equality constraints not yet built. Finish construction of the problem before accessing flow names.")
		return sorted(self._flows, key=self._flows.__getitem__)

	def getMaterialNames(self):
		if not self._eqConstBuilt:
			raise Exception("Equality constraints not yet built. Finish construction of the problem before accessing material names.")
		return sorted(self._materialIdxLookup, key=self._materialIdxLookup.__getitem__)

	def getUpperBounds(self):
		return self._ub.copy()

	def getLowerBounds(self):
		return self._lb.copy()

	def getObjective(self):
		return self._objective.copy()

	def buildEqConst(self):
		if self._eqConstBuilt:
			raise Exception("Equality constraints already built.")

		nMaterials = len(self._materialCoeffs)
		nFlows = len(self._flows)
		A = np.zeros((nMaterials, nFlows))

		# avoid creating duplicate constraints
		self._materialIdxLookup = {}
		for materialIdx, (material, pairs) in enumerate(sorted(self._materialCoeffs.viewitems())):
			self._materialIdxLookup[material] = materialIdx
			for pair in pairs:
				A[materialIdx, pair[1]] = pair[0]
		A_coo = coo_matrix(A)
		row = A_coo.row.astype(np.int32).astype(np.int32)
		col = A_coo.col.astype(np.int32).astype(np.int32)
		data = A_coo.data.astype(np.float64).astype(np.float64)

		# set solver constraints
		self._model.linear_constraints.add(rhs=np.zeros(nMaterials), senses='E'*nMaterials)
		self._model.linear_constraints.set_coefficients(zip(row.tolist(), col.tolist(), data.tolist()))

		self._eqConstBuilt = True

	def _solve(self):
		if self._solved:
			return

		if self._maximize:
			self._model.objective.set_sense(self._model.objective.sense.maximize)
		else:
			self._model.objective.set_sense(self._model.objective.sense.minimize)
		
		self._model.solve()

		self._solved = True
