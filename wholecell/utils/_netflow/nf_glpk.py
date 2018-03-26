# TODO(Jerry): Review ecyglpki source for more prechecks to avoid GLPK exits.
# TODO(Jerry): Add more comments from the GLPK docs.
# TODO(Jerry): Reconsider which exceptions to raise.

from __future__ import division

from collections import defaultdict
from enum import Enum
from math import isinf
import numpy as np
from scipy.sparse import coo_matrix
import swiglpk as glp

from ._base import NetworkFlowProblemBase

class MessageLevel(Enum):
	OFF = glp.GLP_MSG_OFF  # no output
	ERR = glp.GLP_MSG_ERR  # error and warning messages only (default)
	ON  = glp.GLP_MSG_ON   # normal output
	ALL = glp.GLP_MSG_ALL  # full output (including informational messages)
	DBG = glp.GLP_MSG_DBG

class SimplexMethod(Enum):
	PRIMAL = glp.GLP_PRIMAL  # two-phase primal simplex
	DUALP = glp.GLP_DUALP    # two-phase dual simplex
	DUAL = glp.GLP_DUAL      # two-phase dual simplex; fallback to the primal simplex

SOLUTION_STATUS_TO_STRING = {
    glp.GLP_OPT:	'GLP_OPT: optimal',  # solution is optimal
    glp.GLP_FEAS:	'GLP_FEAS: feasible',  # solution is feasible
    glp.GLP_INFEAS:	'GLP_INFEAS: infeasible',  # solution is infeasible
    glp.GLP_NOFEAS:	'GLP_NOFEAS: no feasible',  # problem has no feasible solution
    glp.GLP_UNBND:	'GLP_UNBND: unbounded',  # problem has no unbounded solution
    glp.GLP_UNDEF:	'GLP_UNDEF: undefined',  # solution is undefined
}

_MAXED_OUT = ("GLP_EOBJ{}L: Dual simplex: The objective function being"
			  + " maximized reached its {} limit and continues {}")

SIMPLEX_RETURN_CODE_TO_STRING = {
	0: '',  # successfully solved; not necessarily an optimal solution
	glp.GLP_EBADB:  "GLP_EBADB: Basis is invalid, number of basic variables != number of rows",
	glp.GLP_ESING:  "GLP_ESING: Basis matrix is singular",
	glp.GLP_ECOND:  "GLP_ECOND: Basis matrix is ill-conditioned, it's condition number is too large",
	glp.GLP_EBOUND: "GLP_EBOUND: Some double-bounded variables have incorrect bounds",
	glp.GLP_EFAIL:  "GLP_EFAIL: Solver failure",
	glp.GLP_EOBJLL: _MAXED_OUT.format('L', 'lower', 'decreasing'),
	glp.GLP_EOBJUL: _MAXED_OUT.format('U', 'upper', 'increasing'),
	glp.GLP_EITLIM: "GLP_EITLIM: Iteration limit exceeded",
	glp.GLP_ETMLIM: "GLP_ETMLIM: Time limit exceeded",
	glp.GLP_ENOPFS: "GLP_ENOPFS: Presolver: Problem has no primal feasible solution",
	glp.GLP_ENODFS: "GLP_ENODFS: Presolver: Problem has no dual feasible solution",
}


def _ndarrayToIntArray(array):
	# assert array.dtype == np.int64
	# return glp.intArray_frompointer(array.data)  # TODO: Should work, but from int32 or int64?
	ia = glp.intArray(len(array))  # TODO: Add an element here so the caller doesn't have to?
	for i in xrange(1, len(array)):
		ia[i] = int(array[i])
	return ia

def _ndarrayToDoubleArray(array):
	# assert array.dtype == np.float64
	# return glp.doubleArray_frompointer(array.data)  # TODO: Should work...
	fa = glp.doubleArray(len(array))  # TODO: Add an element here so the caller doesn't have to?
	for i in xrange(1, len(array)):
		fa[i] = float(array[i])
	return fa


class NetworkFlowGLPK(NetworkFlowProblemBase):
	_lowerBoundDefault = 0
	_upperBoundDefault = np.inf

	def __init__(self):
		self._lp = glp.glp_create_prob()
		self._smcp = glp.glp_smcp()
		glp.glp_init_smcp(self._smcp)
		self._smcp.msg_lev = glp.GLP_MSG_ERR
		self._n_vars = 0
		self._n_eq_constraints = 0

		self._flows = {}
		self._lb = {}
		self._ub = {}
		self._objective = {}
		self._materialCoeffs = defaultdict(list)
		self._materialIdxLookup = {}

		self._eqConstBuilt = False
		self._solved = False

	def __del__(self):
		glp.glp_delete_prob(self._lp)


	@property
	def message_level(self):
		"""The message level for terminal output, as an enum value."""
		return MessageLevel(self._smcp.msg_lev)

	@message_level.setter
	def message_level(self, message_level):
		"""Set the message for terminal output using an enum value."""
		self._smcp.msg_lev = message_level.value

	@property
	def simplex_method(self):
		"""The Simplex method option, as an enum value."""
		return SimplexMethod(self._smcp.meth)

	@simplex_method.setter
	def simplex_method(self, simplex_method):
		"""Set the Simplex method option using an enum value."""
		self._smcp.meth = simplex_method.value

	@property
	def simplex_iteration_limit(self):
		"""The Simplex iteration limit, an integer."""
		return self._smcp.it_lim

	@simplex_iteration_limit.setter
	def simplex_iteration_limit(self, limit):
		"""Set the Simplex iteration limit."""
		self._smcp.it_lim = int(limit)

	@property
	def primal_feasible_tolerance(self):
		"""Tolerance used to check if the basic solution is primal feasible."""
		return self._smcp.tol_bnd

	@primal_feasible_tolerance.setter
	def primal_feasible_tolerance(self, tolerance):
		"""Tolerance used to check if the basic solution is primal feasible.
		(Do not change this parameter without detailed understanding its purpose.)
		Default: 1e-7.
		"""
		self._smcp.tol_bnd = float(tolerance)

	@property
	def status_code(self):
		"""The generic status code for the current basic solution."""
		return glp.glp_get_status(self._lp)

	@property
	def status_string(self):
		"""Return the generic status message for the current basic solution."""
		return SOLUTION_STATUS_TO_STRING.get(
			self.status_code, "GLP_?: UNKNOWN SOLUTION STATUS CODE")

	def columnPrimalValue(self, j):
		"""Return the primal value of the structural variable for j-th column."""
		return glp.glp_get_col_prim(self._lp, j)

	def columnDualValue(self, j):
		"""Return the dual value (i.e. reduced cost) of the structural variable
		for the j-th column.
		"""
		return glp.glp_get_col_dual(self._lp, j)

	def rowPrimalValue(self, i):
		"""Return the primal value of the structural variable for i-th row."""
		return glp.glp_get_row_prim(self._lp, i)

	def rowDualValue(self, i):
		"""Return the dual value (i.e. reduced cost) of the structural variable
		for the i-th row.
		"""
		return glp.glp_get_row_dual(self._lp, i)


	def _add_rows(self, n_rows):
		glp.glp_add_rows(self._lp, n_rows)
		self._n_eq_constraints += 1

	def _add_cols(self, n_cols):
		glp.glp_add_cols(self._lp, n_cols)
		self._n_vars += 1

	def _set_col_bounds(self, index, lower, upper):
		if isinf(lower) and isinf(upper):
			typ = glp.GLP_FR  # free (unbounded) variable
		elif lower == upper:
			typ = glp.GLP_FX  # fixed variable
		elif not isinf(lower) and not isinf(upper):
			typ = glp.GLP_DB  # double-bounded variable
		elif isinf(upper):
			typ = glp.GLP_LO  # variable with lower bound
		else:
			typ = glp.GLP_UP  # variable with upper bound

		glp.glp_set_col_bnds(self._lp, index, typ, lower, upper)

	def _set_mat_row(self, row, length, ja_array, ar_array):
		if ja_array.dtype != np.int32:
			raise ValueError("Expected ja_array of dtype np.int32")
		if ar_array.dtype != np.float64:
			raise ValueError("Expected ar_array of dtype np.float64")
		if ja_array.shape != ar_array.shape:
			raise ValueError("Array sizes must match")

		ja = _ndarrayToIntArray(ja_array)
		ar = _ndarrayToDoubleArray(ar_array)
		glp.glp_set_mat_row(self._lp, row, length, ja, ar)

	def _add_eq_constrs(self, ia_array, ja_array, ar_array):
		if ia_array.dtype != np.int32:
			raise ValueError("Expected ia_array of dtype np.int32")
		if ja_array.dtype != np.int32:
			raise ValueError("Expected ja_array of dtype np.int32")
		if ar_array.dtype != np.float64:
			raise ValueError("Expected ar_array of dtype np.float64")
		if ia_array.shape != ja_array.shape:
			raise ValueError("Array sizes must match")
		if ja_array.shape != ar_array.shape:
			raise ValueError("Array sizes must match")

		ia = _ndarrayToIntArray(ia_array)
		ja = _ndarrayToIntArray(ja_array)
		ar = _ndarrayToDoubleArray(ar_array)
		for row in xrange(1, self._n_eq_constraints + 1):
			glp.glp_set_row_bnds(self._lp, row, glp.GLP_FX, 0.0, 0.0)
		n_elems = ia_array.shape[0]
		glp.glp_load_matrix(self._lp, n_elems - 1, ia, ja, ar)

	def _getVar(self, flow):
		if flow in self._flows:
			idx = self._flows[flow]
		else:
			self._add_cols(1)
			idx = len(self._flows)
			self._lb[flow] = self._lowerBoundDefault
			self._ub[flow] = self._upperBoundDefault
			self._set_col_bounds(
				1 + idx,			# GLPK does 1-indexing
				self._lb[flow],
				self._ub[flow],
				)

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

			rowIdx = np.int32(materialIdx + 1)
			colIdxs = np.hstack((np.int32(-1), 1 + np.array(flowIdxs, dtype = np.int32))
				).astype(np.int32, copy = False)
			data = np.hstack((np.nan, np.array(coeffs, dtype = np.float64))
				).astype(np.float64, copy = False)
			self._set_mat_row(rowIdx, np.int32(len(colIdxs) - 1), colIdxs, data)
		else:
			idx = self._getVar(flow)
			self._materialCoeffs[material].append((coefficient, idx))

		self._solved = False


	def flowLowerBoundIs(self, flow, lowerBound):
		idx = self._getVar(flow)
		self._lb[flow] = lowerBound
		self._set_col_bounds(
			1 + idx,				# GLPK does 1 indexing
			self._lb[flow],
			self._ub[flow],
			)

		self._solved = False


	def flowLowerBound(self, flow):
		return self._lb[flow]


	def flowUpperBoundIs(self, flow, upperBound):
		idx = self._getVar(flow)
		self._ub[flow] = upperBound
		self._set_col_bounds(
			1 + idx,				# GLPK does 1 indexing
			self._lb[flow],
			self._ub[flow],
			)

		self._solved = False


	def flowUpperBound(self, flow):
		return self._ub[flow]


	def flowObjectiveCoeffIs(self, flow, coefficient):
		idx = self._getVar(flow)
		self._objective[flow] = coefficient
		glp.glp_set_obj_coef(
			self._lp,
			1 + idx,				# GLPK does 1 indexing
			coefficient
			)

		self._solved = False


	def flowRates(self, flows):
		if isinstance(flows, basestring):
			flows = (flows,)

		self._solve()

		return np.array(
			[glp.glp_get_col_prim(self._lp, 1 + self._flows[flow])
			 if flow in self._flows else None
			 for flow in flows]
			)

	def rowDualValues(self, materials):
		if not self._eqConstBuilt:
			raise Exception("Equality constraints not yet built. Finish construction of the problem before accessing dual values.")

		self._solve()

		return np.array(
			[glp.glp_get_row_dual(self._lp, 1 + self._materialIdxLookup[material])
			 if material in self._materialIdxLookup else None
			 for material in materials]
			)

	def columnDualValues(self, fluxNames):
		if not self._eqConstBuilt:
			raise Exception("Equality constraints not yet built. Finish construction of the problem before accessing dual values.")

		self._solve()

		return np.array(
			[glp.glp_get_col_dual(self._lp, 1 + self._flows[fluxName])
			 if fluxName in self._flows else None
			 for fluxName in fluxNames]
			)

	def objectiveValue(self):
		"""The current value of the objective function."""
		return glp.glp_get_obj_val(self._lp)


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

		self._add_rows(len(self._materialCoeffs))
		A = np.zeros((len(self._materialCoeffs), len(self._flows)))
		# avoid creating duplicate constraints
		self._materialIdxLookup = {}
		for materialIdx, (material, pairs) in enumerate(sorted(self._materialCoeffs.viewitems())):
			self._materialIdxLookup[material] = materialIdx
			for pair in pairs:
				A[materialIdx, pair[1]] = pair[0]
		A_coo = coo_matrix(A)
		row = np.hstack((np.int32(-1), 1 + A_coo.row.astype(np.int32))
			).astype(np.int32, copy = False)
		col = np.hstack((np.int32(-1), 1 + A_coo.col.astype(np.int32))
			).astype(np.int32, copy = False)
		data = np.hstack((np.nan, A_coo.data.astype(np.float64))
			).astype(np.float64, copy = False)
		self._add_eq_constrs(row, col, data)

		self._eqConstBuilt = True


	def _solve(self):
		if self._solved:
			return

		if self._maximize:
			glp.glp_set_obj_dir(self._lp, glp.GLP_MAX)
		else:
			glp.glp_set_obj_dir(self._lp, glp.GLP_MIN)

		result = glp.glp_simplex(self._lp, self._smcp)

		if result != 0:
			raise RuntimeError(SIMPLEX_RETURN_CODE_TO_STRING.get(
				result, "GLP_?: UNKNOWN SOLVER RETURN VALUE"))
		if self.status_code != glp.GLP_OPT:
			raise RuntimeError(self.status_string)

		self._solved = True
