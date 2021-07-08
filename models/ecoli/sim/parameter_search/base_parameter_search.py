"""

"""

import os
import pickle
from typing import Any, Dict, Iterable

import numpy as np

from wholecell.io.tablereader import TableReader
from wholecell.sim.simulation import DEFAULT_SIMULATION_KWARGS


DEFAULT_CLI_KWARGS = {
	'init_sims': 1,
	'generations': 1,
	}


class RawParameter():
	def __init__(self, attr: str, id_columns: Dict[str, Any], columns: Iterable[str], name: str):
		self._attr = attr
		self._id = id_columns
		self._columns = columns
		self._name = name
		self._cached_index = None

	def get_row(self, raw_data):
		def find_row(attr):
			for i, row in enumerate(attr):
				for col, val in self._id.items():
					if row[col] != val:
						break
				else:
					self._cached_index = i
					return row

			raise RuntimeError('Could not find a row matching the given ID columns.')

		attr = getattr(raw_data, self._attr)

		if self._cached_index is not None:
			row = attr[self._cached_index]
			for col, val in self._id.items():
				if row[col] != val:
					row = find_row(attr)
					break
		else:
			row = find_row(attr)

		return row

	def get_param(self, raw_data):
		row = self.get_row(raw_data)
		return np.mean([row[col] for col in self._columns])

	def set_param(self, raw_data, value):
		row = self.get_row(raw_data)
		for col in self._columns:
			row[col] = value

	def __str__(self):
		return self._name


class BaseParameterSearch():
	parca_args = {}
	# TODO: handle raw and sim params the same - create a class for SimParameter and combine attributes below
	_raw_params = ()
	_sim_params = ()
	_init_raw_params = {}
	_init_sim_params = {}
	sims_to_run = ()

	def __init__(self):
		self.variant_name = self.__class__.__name__
		self.n_parameters = len(self._raw_params) + len(self._sim_params)
		self.raw_params = {p: None for p in self._raw_params}
		self.sim_params = {p: None for p in self._sim_params}
		self.initialized = False

	def get_objective(self, sim_out_dirs, sim_data_files):
		raise NotImplementedError('Need to implement in a subclass.')

	def initialize(self, raw_data_file, sim_data_file, iteration):
		# If no raw params, raw_data will not be saved with each iteration so
		# this would cause issues loading from a specific iteration
		if self.raw_params:
			with open(raw_data_file, 'rb') as f:
				raw_data = pickle.load(f)
			for param in self.raw_params:
				value = self.get_attr(raw_data, param)
				if iteration == 0:
					self.raw_params[param] = self._init_raw_params.get(str(param), value)
				else:
					self.raw_params[param] = value

		with open(sim_data_file, 'rb') as f:
			sim_data = pickle.load(f)
		for param in self.sim_params:
			value = self.get_attr(sim_data, param)
			if iteration == 0:
				self.sim_params[param] = self._init_sim_params.get(param, value)
			else:
				self.sim_params[param] = value

		self.initialized = True

	def get_sim_params(self, sim_dir, variants):
		all_params = []
		for variant in variants:
			for index, sim_params in enumerate(self.sims_to_run):
				params = DEFAULT_SIMULATION_KWARGS.copy()
				params.update(DEFAULT_CLI_KWARGS)
				params.update(sim_params)

				params['variant type'] = self.variant_name
				params['variant'] = variant
				params['sim dir'] = sim_dir
				params['index'] = index

				all_params.append(params)

		return all_params

	def get_attr(self, obj, attr, default=None):
		if isinstance(attr, RawParameter):
			return attr.get_param(obj)
		else:
			attrs = attr.split('.')
			for a in attrs:
				if hasattr(obj, a):
					obj = getattr(obj, a)
				else:
					return default

			return obj

	def set_attr(self, obj, attr, val):
		if isinstance(attr, RawParameter):
			attr.set_param(obj, val)
		else:
			attrs = attr.split('.')
			for a in attrs[:-1]:
				obj = getattr(obj, a)
			setattr(obj, attrs[-1], val)

	def print_update(self):
		def print_params(params, label):
			if params:
				print(f'{label} parameters:')
				for p, val in params.items():
					print(f'\t{p}: {val}')

		print_params(self.raw_params, 'Raw data')
		print_params(self.sim_params, 'Sim data')

	def read_column(self, out_dir, table, column):
		return TableReader(os.path.join(out_dir, table)).readColumn(column)

	def read_attribute(self, out_dir, table, attr):
		return TableReader(os.path.join(out_dir, table)).readAttribute(attr)
