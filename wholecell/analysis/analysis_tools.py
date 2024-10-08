"""
Analysis script toolbox functions
"""

from __future__ import annotations

import os
from typing import Callable, Iterator, List, Optional, Sequence, Tuple, Union

import numpy as np

from wholecell.io.tablereader import TableReader
from wholecell.utils import filepath

LOW_RES_DIR = 'low_res_plots'
SVG_DIR = 'svg_plots'
HTML_DIR = 'html_plots'
LOW_RES_DPI = 120


def exportFigure(plt, plotOutDir, plotOutFileName, metadata=None, transparent=False,
        dpi=LOW_RES_DPI, extension=None):

	if metadata is not None and "analysis_type" in metadata:
		analysis_type = metadata["analysis_type"]

		if analysis_type == 'single':
			# Format metadata signature for single gen figure
			metadata_signature = "_".join([
				str(metadata["time"])[:13],
				str(metadata["variant_function"]),
				str(metadata["variant_index"]),
				"Seed",
				str(metadata["seed"]),
				"Gen",
				str(metadata["gen"]) + '/' + str(int(metadata["total_gens"]) - 1),
				"Githash",
				str(metadata["git_hash"])[:10],
				"Desc",
				str(metadata["description"])
				])

		elif analysis_type == 'multigen':
			# Format metadata signature for multi gen figure
			metadata_signature = "_".join([
				str(metadata["time"][:13]),
				str(metadata["variant_function"]),
				str(metadata["variant_index"]),
				"Seed",
				str(metadata["seed"]),
				str(metadata["total_gens"]),
				"gens",
				"Githash",
				str(metadata["git_hash"])[:10],
				"Desc",
				str(metadata["description"])
				])

		elif analysis_type == 'cohort':
			# Format metadata signature for cohort figure
			metadata_signature = "_".join([
				str(metadata["time"][:13]),
				str(metadata["variant_function"]),
				str(metadata["variant_index"]),
				str(metadata["total_gens"]),
				"gens",
				"Githash",
				str(metadata["git_hash"])[:10],
				"Desc",
				str(metadata["description"])
				])

		elif analysis_type == 'variant':
			# Format metadata signature for variant figure
			metadata_signature = "_".join([
				str(metadata["time"][:13]),
				str(metadata["total_variants"]),
				"variants",
				str(metadata["total_gens"]),
				"gens",
				"Githash",
				str(metadata["git_hash"])[:10],
				"Desc",
				str(metadata["description"])
				])

		elif analysis_type == 'parca':
			# Format metadata signature for parca figure
			metadata_signature = "_".join([
				str(metadata["time"][:13]),
				"Githash",
				str(metadata["git_hash"])[:10],
				"Desc",
				str(metadata["description"])
				])

		elif analysis_type == 'comparison':
			# Format metadata signature for a comparison figure
			metadata_signature = "_".join([
				str(metadata["time"][:13]),
				str(metadata["total_variants"]),
				"variants",
				str(metadata["total_gens"]),
				"gens",
				"Githash",
				str(metadata["git_hash"])[:10],
				"Desc",
				str(metadata["description"])
				])

		else:
			raise ValueError('Unknown analysis_type {}'.format(analysis_type))

		# Add metadata signature to the bottom of the plot
		# Don't accidentally trigger $TeX formatting$.
		metadata_signature = metadata_signature.replace('$', '')
		plt.figtext(0,0, metadata_signature, size=8)

	# Make folders for holding alternate types of images
	filepath.makedirs(plotOutDir, LOW_RES_DIR)
	filepath.makedirs(plotOutDir, SVG_DIR)

	# Save images
	if extension:
		# Only save one type in main analysis directory if extension is given
		plt.savefig(os.path.join(plotOutDir, plotOutFileName + extension), dpi=dpi, transparent=transparent)
	else:
		# Save all image types
		plt.savefig(os.path.join(plotOutDir, plotOutFileName + '.pdf'), transparent=transparent)
		plt.savefig(os.path.join(plotOutDir, SVG_DIR, plotOutFileName + '.svg'), transparent=transparent)
		plt.savefig(os.path.join(plotOutDir, LOW_RES_DIR, plotOutFileName + '.png'), dpi=dpi, transparent=transparent)

def _check_bulk_inputs(mol_names: Union[Tuple[Sequence[str], ...], Sequence[str]]) -> Tuple[Sequence[str], ...]:
	"""
	Use to check and adjust mol_names inputs for functions that read bulk
	molecules to get consistent argument handling in both functions.
	"""

	# Wrap an array in a tuple to ensure correct dimensions
	if not isinstance(mol_names, tuple):
		mol_names = (mol_names,)

	# Check for string instead of array since it will cause mol_indices lookup to fail
	for names in mol_names:
		if isinstance(names, (bytes, str)):
			raise Exception('mol_names tuple must contain arrays not strings like {!r}'.format(names))

	return mol_names

def _remove_first(remove_first: bool):
	"""Slice to remove the first time point entry if remove_first is True"""
	return slice(1, None) if remove_first else slice(None)

def read_bulk_molecule_counts(sim_out_dir, mol_names):
	# type: (str, Union[Tuple[Sequence[str], ...], Sequence[str]]) -> Iterator[np.ndarray]
	'''
	Reads a subset of molecule counts from BulkMolecules using the indexing method
	of readColumn. Should only be called once per simulation being analyzed with
	all molecules of interest.

	Args:
		sim_out_dir: path to the directory with simulation output data
		mol_names: tuple of list-likes of strings or a list-like of strings
			that name molecules to read the counts for. A single array will be
			converted to a tuple for processing.

	Returns:
		generator of ndarray: int counts with all time points on the first dimension
		and each molecule of interest on the second dimension. The number of
		generated arrays will be separated based on the input dimensions of mol_names
		(ie if mol_names is a tuple of two arrays, two arrays will be generated).

	Example use cases:
		names1 = ['ATP[c]', 'AMP[c]']
		names2 = ['WATER[c]']

		# Read one set of molecules
		(counts1,) = read_bulk_molecule_counts(sim_out_dir, names1)

		# Read two or more sets of molecules
		(counts1, counts2) = read_bulk_molecule_counts(sim_out_dir, (names1, names2))

	TODO: generalize to any TableReader, not just BulkMolecules, if readColumn method
	is used for those tables.
	'''

	mol_names = _check_bulk_inputs(mol_names)

	bulk_reader = TableReader(os.path.join(sim_out_dir, 'BulkMolecules'))

	bulk_molecule_names = bulk_reader.readAttribute("objectNames")
	mol_indices = {mol: i for i, mol in enumerate(bulk_molecule_names)}

	lengths = [len(names) for names in mol_names]
	indices = np.hstack([[mol_indices[mol] for mol in names] for names in mol_names])
	bulk_counts = bulk_reader.readColumn('counts', indices, squeeze=False)

	start_slice = 0
	for length in lengths:
		counts = bulk_counts[:, start_slice:start_slice + length].squeeze()
		start_slice += length
		yield counts

def read_stacked_bulk_molecules(
		cell_paths: np.ndarray,
		mol_names: Union[Tuple[Sequence[str], ...], Sequence[str]],
		remove_first: bool = False,
		ignore_exception: bool = False,
		) -> List[np.ndarray]:
	"""
	Reads bulk molecule counts from multiple cells and assembles each group
	into a single array.

	Args:
		cell_paths: paths to all cells to read data from (directories should
			contain a simOut/ subdirectory), typically the return from
			AnalysisPaths.get_cells()
		mol_names: tuple of list-likes of strings or a list-like of strings
			that name molecules to read the counts for. A single array will be
			converted to a tuple for processing.
		remove_first: if True, removes the first column of data from each cell
			which might be set to a default value in some cases
		ignore_exception: if True, ignores any exceptions encountered while reading

	Returns:
		stacked data (n time points) if single molecule or
			(n time points, m molecules) if multiple molecules for each group
			in mol_names
	"""

	mol_names = _check_bulk_inputs(mol_names)

	data = [[] for _ in mol_names]  # type: List[List[np.ndarray]]
	for sim_dir in cell_paths:
		sim_out_dir = os.path.join(sim_dir, 'simOut')
		try:
			for i, counts in enumerate(read_bulk_molecule_counts(sim_out_dir, mol_names)):
				data[i].append(counts[_remove_first(remove_first)])
		except Exception as e:
			if ignore_exception:
				print(f'Ignored exception in read_stacked_bulk_molecules for {sim_out_dir}: {e!r}')
				continue
			else:
				raise

	# Use vstack for 2D or hstack for 1D to get proper dimension alignments
	return [np.vstack(d) if len(d[0].shape) > 1 else np.hstack(d) for d in data]

def read_stacked_columns(cell_paths: np.ndarray, table: str, column: str,
		remove_first: bool = False, ignore_exception: bool = False,
		fun: Optional[Callable] = None) -> np.ndarray:
	"""
	Reads column data from multiple cells and assembles into a single array.

	Args:
		cell_paths: paths to all cells to read data from (directories should
			contain a simOut/ subdirectory), typically the return from
			AnalysisPaths.get_cells()
		table: name of the table to read data from
		column: name of the column to read data from
		remove_first: if True, removes the first column of data from each cell
			which might be set to a default value in some cases
		ignore_exception: if True, ignores any exceptions encountered while reading
		fun: function to apply to data in each generation (eg. np.mean will
			return and array with the mean value for each generation instead
			of each time point)

	Returns:
		stacked data (n time points, m subcolumns)

	TODO:
		add common processing options:
		- mean per cell cycle
		- normalize to a time point
	"""

	if fun is None:
		fun = lambda x: x

	data = []
	for sim_dir in cell_paths:
		sim_out_dir = os.path.join(sim_dir, 'simOut')
		try:
			reader = TableReader(os.path.join(sim_out_dir, table))
			data.append(fun(reader.readColumn(column, squeeze=False)[_remove_first(remove_first)]))
		except Exception as e:
			if ignore_exception:
				print(f'Ignored exception in read_stacked_columns for {sim_out_dir}: {e!r}')
				continue
			else:
				raise

	return np.vstack(data) if data else np.array([])

def stacked_cell_identification(cell_paths: np.ndarray, table: str, column: str,
		remove_first: bool = False, ignore_exception: bool = False,
		) -> np.ndarray:
	"""
	Reads column data from multiple cells and assembles into a single array.

	Args:
		cell_paths: paths to all cells to read data from (directories should
			contain a simOut/ subdirectory), typically the return from
			AnalysisPaths.get_cells()
		table: name of the table to read data from
		column: name of the column to read data from
		remove_first: if True, removes the first column of data from each cell
			which might be set to a default value in some cases
		ignore_exception: if True, ignores any exceptions encountered while reading

	Returns:
		stacked data (n time points, m subcolumns) with a unique numerical
			identifier for each cell in cell_paths

	"""

	data = []
	counter = 0
	for sim_dir in cell_paths:
		sim_out_dir = os.path.join(sim_dir, 'simOut')
		try:
			reader = TableReader(os.path.join(sim_out_dir, table))
			column_data = reader.readColumn(column, squeeze=False)
			# [_remove_first(remove_first)]
			data.append(np.ones_like(column_data)*counter)

		except Exception as e:
			if ignore_exception:
				print(f'Ignored exception in stacked_cell_identification for'
					  f' {sim_out_dir}: {e!r}')
				continue
			else:
				raise

		counter += 1

	return np.vstack(data) if data else np.array([])

def stacked_cell_threshold_mask(cell_paths: np.ndarray, table: str, column:
str,
		threshold_value: float, remove_first: bool = False,
		ignore_exception: bool = False, fun: Optional[Callable] = None) -> np.ndarray:
	"""
	Returns single boolean array to indicate whether each value in the column
	data (from multiple cells) occurs before the first
	value >= threshold_value in that simulation (seed).

	Args:
		cell_paths: paths to all cells to read data from (directories should
			contain a simOut/ subdirectory), typically the return from
			AnalysisPaths.get_cells()
		table: name of the table to read data from
		column: name of the column to read data from
		threshold_value: values after the first value >= threshold_value in
			each cell will correspond to false
		remove_first: if True, removes the first column of data from each cell
			which might be set to a default value in some cases
		ignore_exception: if True, ignores any exceptions encountered while
			reading
		fun: function to apply to data in each generation (eg. np.mean will
			return and array with the mean value for each generation instead
			of each time point)

	Returns:
		stacked data (n time points, m subcolumns) with a boolean mask for each
		 	cell in cell_paths

	TODO: extend to be compatible with multiple values per cell

	"""

	if fun is None:
		fun = lambda x: x

	boolean_data = []
	"""
	For each simulation, first_threshold_not_seen will keep track of whether we
	have already seen a value >= threshold_value
	"""
	first_threshold_not_seen = True #
	for sim_dir in cell_paths:
		sim_out_dir = os.path.join(sim_dir, 'simOut')
		try:
			gen_index = int(os.path.basename(os.path.dirname(sim_dir))[-6:])
			if gen_index == 0:
				first_threshold_not_seen = True # reset flag for new simulation
			reader = TableReader(os.path.join(sim_out_dir, table))
			column_data = fun(reader.readColumn(column, squeeze=False)
							  [_remove_first(remove_first)])
			assert len(column_data) == 1, \
				"stacked_cell_threshold_mask cannot yet accomodate multiple " \
				"values per generation"
			if (first_threshold_not_seen == True) and \
					(column_data >= threshold_value):
				first_threshold_not_seen = False
			boolean_data.append(bool(first_threshold_not_seen))

		except Exception as e:
			if ignore_exception:
				print(f'Ignored exception in stacked_cell_threshold_mask for'
					  f' {sim_out_dir}: {e!r}')
				continue
			else:
				raise

	return np.vstack(boolean_data) if boolean_data else np.array([])