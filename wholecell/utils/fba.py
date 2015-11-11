
from __future__ import division

from collections import defaultdict

import numpy as np
import cvxpy as cvx

class FluxBalanceAnalysis(object):
	def __init__(self, reaction_stoich, external_molecules, objective,
			reversible_reactions):

		# Parse reactions

		self._reaction_rates = {}

		metabolite_constraints = defaultdict(dict)

		for reaction_id, stoich in reaction_stoich.viewitems():
			rate = cvx.Variable(name = reaction_id)

			self._reaction_rates[reaction_id] = rate

			for met_id, coeff in stoich.viewitems():
				metabolite_constraints[met_id][rate] = coeff

		# Parse exchanges

		self._external_exchange_rates = {}
		self._external_exchange_limits = {}

		for met_id in external_molecules:
			rate = cvx.Variable(name = '{} external exchange flux'.format(met_id))

			self._external_exchange_rates[met_id] = rate

			metabolite_constraints[met_id][rate] = -1 # points OUT

			limit = cvx.Parameter(name = '{} external exchange limit'.format(met_id), value = 0)

			self._external_exchange_limits[rate] = limit

		# Parse objective, internal exchange

		objective_components = {}
		self._internal_exchange_rates = {}
		self._internal_exchange_limits = {}
		self._output = {}

		for met_id, conc in objective.viewitems():
			rate_out = cvx.Variable(name = '{} output flux'.format(met_id))
			objective_components[rate_out] = conc
			metabolite_constraints[met_id][rate_out] = -1 # points OUT

			self._output[met_id] = rate_out

			rate_in = cvx.Variable(name = '{} internal exchange flux'.format(met_id))
			self._internal_exchange_rates[met_id] = rate_in
			metabolite_constraints[met_id][rate_in] = -1 # points OUT

			limit = cvx.Parameter(name = '{} internal exchange limit'.format(met_id), value = 0)
			self._internal_exchange_limits[rate_in] = limit

		self._objective = sum(
			cvx.abs(1 - rate/conc)
			# cvx.square(1 - rate/conc)
			for rate, conc in objective_components.viewitems()
			)

		# Create constraints
		constraints = []

		constraints += [
			(sum(coeff * rate for rate, coeff in const.viewitems()) == 0)
			for const in metabolite_constraints.viewvalues()
			]

		# reaction reversibility (TODO: handle better)
		constraints += [
			(rate >= 0)
			for reaction_id, rate in self._reaction_rates.viewitems()
			if reaction_id not in reversible_reactions
			]

		# exchange constraints
		constraints += [
			(rate >= -limit)
			for rate, limit in self._external_exchange_limits.viewitems()
			]

		constraints += [
			(rate == -limit) # inputs are forced in
			for rate, limit in self._internal_exchange_limits.viewitems()
			]

		constraints += [
			(rate >= 0)
			for rate in objective_components.viewkeys()
			]

		# Create problem

		self._problem = cvx.Problem(
			cvx.Minimize(self._objective),
			constraints
			)

		self._solved = False


	def external_limit_is(self, met_id, value):
		# TODO: disentangle this
		variable = self._external_exchange_rates[met_id]
		parameter = self._external_exchange_limits[variable]

		parameter.value = value

		self._solved = False


	def internal_limit_is(self, met_id, value):
		# TODO: see above
		variable = self._internal_exchange_rates[met_id]
		parameter = self._internal_exchange_limits[variable]

		parameter.value = value

		self._solved = False


	def output(self, met_id):
		self._solve()

		return self._output[met_id].value


	def _solve(self):
		if self._solved:
			return


		self._problem.solve(solver = cvx.GUROBI, verbose = True)


		if self._problem.status == cvx.OPTIMAL:
			self._solved = True

		else:
			raise Exception('Solver failed with status "{}"'.format(self._problem.status))

