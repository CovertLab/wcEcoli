'''
AnalysisPaths: object for easily accessing file paths to simulations based on
variants, seeds, and generation.
'''

from __future__ import annotations

from os import listdir
from os.path import dirname, isdir, isfile, join
import pickle
import re
from itertools import chain
from typing import cast, Iterable, List, Optional, Tuple, Union

import numpy as np

from reconstruction.ecoli.simulation_data import SimulationDataEcoli
from validation.ecoli.validation_data import ValidationDataEcoli
from wholecell.utils import constants
from wholecell.utils import filepath as fp


class AnalysisPaths(object):
	'''
	Object for easily accessing file paths to simulations based on
	variants, seeds, and generation. Within a specified variant you
	can then access all seeds and/or generations.

	Example:
		from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
		ap = AnalysisPaths(simOutDir, variant_plot=True)
		ap.get_cells(variant = [0,3], seed = [0], generation = [1,2])

	Above example should return all paths corresponding to variants
	0 and 3, seed 0, and generations 1 and 2. If a field is left blank
	it is assumed that all values for that field are desired. If all
	fields are left blank all cells will be returned.

	* For an all_variant_plot, out_dir must be a top level simulation output
	  dir. It runs over all variant directories, including those with the high
	  order variant digit (1xxxxx) created for the Parca "--operons=both"
	  secondary simData case. These variant analysis classes must:
	  * not assume that variant indexes are contiguous;
	  * use the high order variant index digit to select the correct simData and
	    variantData for each variant.
	  AnalysisPaths has support methods for all_variant_plot classes.
	* For a variant_plot, out_dir must be a top level simulation output dir.
	  This plot type skips the secondary simData case for compatibility with
	  existing variant analysis classes.
	* For a cohort_plot, out_dir must be a variant output dir.
	* For a multi_gen_plot, out_dir must be a seed output dir.
	'''

	def __init__(self, out_dir, *,
				 all_variant_plot: bool = False, variant_plot: bool = False,
				 multi_gen_plot: bool = False, cohort_plot: bool = False) -> None:
		assert (all_variant_plot + variant_plot + multi_gen_plot + cohort_plot
				) == 1, "Must specify exactly one plot type!"

		generation_dirs = []  # type: List[str]
		if variant_plot or all_variant_plot:
			# Find suitable variant directories in the given simulation output dir
			variant_pattern = re.compile(r'.+_0\d{5}' if variant_plot else r'.+_\d{6}')
			all_dirs = listdir(out_dir)
			variant_out_dirs = []
			for directory in all_dirs:
				if variant_pattern.fullmatch(directory):
					variant_out_dirs.append(join(out_dir, directory))

			if len(variant_out_dirs) == 0:
				raise Exception("Variant out_dir doesn't contain variants!")

			# Get all seed directories in each variant directory
			seed_out_dirs = []
			for variant_dir in variant_out_dirs:
				all_dirs = listdir(variant_dir)
				for directory in all_dirs:
					if re.match(r'^\d{6}$', directory):
						seed_out_dirs.append(join(variant_dir, directory))

			# Get all generation files for each seed
			generation_dirs = []
			for seed_dir in seed_out_dirs:
				generation_dirs.extend(chain.from_iterable(self._get_generations(seed_dir)))

		elif cohort_plot:
			# Find all seed directories in the given variant directory
			seed_out_dirs = []
			all_dirs = listdir(out_dir)
			for directory in all_dirs:
				if re.match(r'^\d{6}$', directory):
					seed_out_dirs.append(join(out_dir, directory))

			# Get all generation files for each seed
			generation_dirs = []
			for seed_dir in seed_out_dirs:
				generation_dirs.extend(chain.from_iterable(self._get_generations(seed_dir)))


		elif multi_gen_plot:
			# Find all generation directories in the given seed directory
			generation_dirs = list(chain.from_iterable(self._get_generations(out_dir)))

		self._path_data = np.zeros(len(generation_dirs), dtype=[
			("path", "U500"),
			("variant", "i8"),
			("seed", "i8"),
			("generation", "i8"),
			("variantkb", "U500")
			])

		generations = []
		seeds = []
		variants = []
		variant_kb = []
		for filePath in generation_dirs:
			# Find generation
			matches = re.findall(r'generation_\d{6}', filePath)
			if len(matches) > 1:
				raise Exception("Expected only one match for generation!")
			generations.append(int(matches[0][-6:]))

			gen_subdir_index = filePath.rfind('generation_')

			# Extract the seed index
			# Assumes: 6-digit seed index SSSSSS in 'SSSSSS/generation_...'
			seeds.append(int(filePath[gen_subdir_index - 7 : gen_subdir_index - 1]))

			# Extract the variant index
			# Assumes: 6-digit variant index VVVVVV in 'VARIANT-TYPE_VVVVVV/SSSSSS/generation_...'
			variants.append(int(filePath[gen_subdir_index - 14 : gen_subdir_index - 8]))

			# Find the variant kb pickle
			variant_kb.append(
				join(filePath[: gen_subdir_index - 8],
					 constants.VKB_DIR, constants.SERIALIZED_SIM_DATA_MODIFIED))

		self._path_data["path"] = generation_dirs
		self._path_data["variant"] = variants
		self._path_data["seed"] = seeds
		self._path_data["generation"] = generations
		self._path_data["variantkb"] = variant_kb

		self.variants = sorted(set(variants))
		operons = set()
		base_variants = set()
		for variant_index in self.variants:
			operon, base_variant = fp.split_variant_index(variant_index)
			operons.add(operon)
			base_variants.add(base_variant)
		self.operons = sorted(operons)  # unique operon indexes
		self.base_variants = sorted(base_variants)  # unique base variant indexes

		self.n_generation = len(set(generations))
		self.n_variant = len(self.variants)  # CAUTION: variant indexes are discontiguous when operons=both
		self.n_seed = len(set(seeds))

	def get_cells(self, variant = None, seed = None, generation = None):
		# type: (Optional[Iterable[Union[int, str]]], Optional[Iterable[int]], Optional[Iterable[int]]) -> np.ndarray
		"""Returns file paths for all the simulated cells matching the given
		indexes.
		"""
		if variant is None:
			variantBool = np.ones(self._path_data.shape)
		else:
			variantBool = self._set_match("variant", variant)

		if seed is None:
			seedBool = np.ones(self._path_data.shape)
		else:
			seedBool = self._set_match("seed", seed)

		if generation is None:
			generationBool = np.ones(self._path_data.shape)
		else:
			generationBool = self._set_match("generation", generation)

		return self._path_data['path'][np.logical_and.reduce((variantBool, seedBool, generationBool))]

	def get_variant_kb(self, variant):
		# type: (Union[int, str]) -> str
		kb_path = np.unique(self._path_data['variantkb'][np.where(self._path_data["variant"] == variant)])
		assert kb_path.size == 1
		return kb_path[0]

	def get_variants(self):
		# type: () -> List[Union[int, str]]
		"""Return all the variants (like self.variants)."""
		return sorted(np.unique(self._path_data["variant"]))

	def get_cell_variant(self, path: str) -> int:
		"""Return the variant index for the given get_cells() path."""
		return self._path_data['variant'][self._path_index(path)]

	def get_cell_seed(self, path: str) -> int:
		"""Return the seed for the given get_cells() path."""
		return self._path_data['seed'][self._path_index(path)]

	def get_cell_generation(self, path: str) -> int:
		"""Return the generation number for the given get_cells() path."""
		return self._path_data['generation'][self._path_index(path)]

	def get_cell_variant_kb(self, path: str) -> str:
		"""Return the variant kb path for the given get_cells() path."""
		return self._path_data['variantkb'][self._path_index(path)]

	def _path_index(self, path: str) -> int:
		"""Return the index into _path_data of the given get_cells() path."""
		indexes = np.where(self._path_data['path'] == path)[0]
		assert indexes.size == 1
		return indexes[0]

	def read_sim_data_files(self, simDataFile: str) -> Tuple[
			SimulationDataEcoli, SimulationDataEcoli]:
		"""Return the primary sim_data object read from simDataFile and the
		secondary sim_data object read from the adjacent kb-poly/ directory,
		or the primary one if there isn't a secondary.
		"""
		with open(simDataFile, 'rb') as f:
			sim_data1 = sim_data2 = pickle.load(f)

		sim_path = dirname(dirname(simDataFile))
		secondary_sim_data_file = join(
			sim_path, constants.PKB_DIR, constants.SERIALIZED_SIM_DATA_FILENAME)
		if isfile(secondary_sim_data_file):
			with open(secondary_sim_data_file, 'rb') as f:
				sim_data2 = pickle.load(f)

		return sim_data1, sim_data2

	def read_validation_data_files(self, validationDataFile: str) -> Tuple[
			ValidationDataEcoli, ValidationDataEcoli]:
		"""Return the primary validation_data object read from validationDataFile
		and the secondary validation_data object read from the adjacent kb-poly/
		directory, or the primary one if there isn't a secondary.
		"""
		with open(validationDataFile, 'rb') as f:
			validation_data1 = validation_data2 = pickle.load(f)

		sim_path = dirname(dirname(validationDataFile))
		secondary_validation_data_file = join(
			sim_path, constants.PKB_DIR, constants.SERIALIZED_VALIDATION_DATA)
		if isfile(secondary_validation_data_file):
			with open(secondary_validation_data_file, 'rb') as f:
				validation_data2 = pickle.load(f)

		return validation_data1, validation_data2

	def _get_generations(self, directory):
		# type: (str) -> List[List[str]]
		"""Get a sorted list of the directory's generation paths, each as a list
		of daughter cell paths.
		ASSUMES: directory contains "generation_000000" thru "generation_GGGGGG".
		"""
		generation_files = [
			join(directory, f) for f in listdir(directory)
			if isdir(join(directory, f)) and "generation" in f]  # type: List[str]
		generations = [[] for _ in generation_files]  # type: List[List[str]]
		for gen_file in generation_files:
			generations[int(gen_file[gen_file.rfind('_') + 1:])] = self._get_individuals(gen_file)
		return generations

	def _get_individuals(self, directory):
		# type: (str) -> List[str]
		"""Get a sorted list of the directory's daughter cell paths, each of
		them a place for simOut/ and plotOut/ subdirs.
		ASSUMES: directory is a generation directory like "generation_000001",
		each containing numbered daughter subdirs "000000" thru "DDDDDD".
		"""
		individual_files = [
			join(directory, f) for f in listdir(directory)
			if isdir(join(directory, f))]  # type: List[str]
		individuals = [''] * len(individual_files)  # type: List[str]
		for ind_file in individual_files:
			individuals[int(ind_file[ind_file.rfind('/') + 1:])] = ind_file
		return individuals

	def _set_match(self, field, value):
		# type: (str, Iterable[Union[int, str]]) -> np.ndarray
		union = np.zeros(self._path_data[field].size)
		for x in value:
			union = cast(np.ndarray, np.logical_or(self._path_data[field] == x, union))
		return union
