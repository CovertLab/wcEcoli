"""
Abstract base class for a Firetask that runs a category of analysis plots.

If the `DEBUG_GC` environment variable is true, enable memory leak detection.
"""

import abc
from collections import OrderedDict
import datetime
import importlib
import os
import sys
import time
import traceback
from typing import Dict, List

from fireworks import FiretaskBase
import matplotlib as mpl
from PIL import Image

from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from wholecell.utils import data
from wholecell.utils import parallelization
import wholecell.utils.filepath as fp


# Used to set the backend to Agg before pyplot imports in other scripts.
# Other configuration settings can be added to the file as well.
mpl.rc_file(fp.MATPLOTLIBRC_FILE)

SUB_DIRECTORIES = {'.png': 'low_res_plots'}

CORE_TAG = 'CORE'
DEFAULT_TAG = 'DEFAULT'
VARIANT_TAG = 'VARIANT'

class AnalysisBase(FiretaskBase):
	"""Base class for analysis plot Firetasks.

	Each subclass should set the usual Firetask class variables _fw_name,
	required_params, and optional_params; also

		TAGS = the dictionary that maps tag names to lists of analysis plot
			module filenames (in this category) to run.

			Use all uppercase for tag names so they don't conflict with module
			filenames.

			The tag 'DEFAULT' will run 'CORE' and 'VARIANT' when the argument
			list is empty.

			The tag 'CORE' lists the plots that are most commonly used for
			analysis.

			The tag 'VARIANT' lists any analysis plots specific to the variant
			that was simulated.

			The tag 'ACTIVE' lists all active plots in this category. The
			nightly build should run 'ACTIVE'.

		MODULE_PATH = the module pathname for loading this subclass's analysis
			plots.

	Optional params include plot, output_filename_prefix, cpus.
	"""

	analysis_path_options = {}  # type: Dict[str, bool]

	@abc.abstractmethod
	def plotter_args(self, module_filename):
		"""(Abstract) Return a tuple of arguments to pass to the analysis plot
		class' `main()` method.
		"""
		raise NotImplementedError

	def expand_plot_names(self, plot_names, name_dict):
		'''Recursively expand TAGS and plot names that lack the '.py' suffix,
		adding the names to name_dict. name_dict is an OrderedDict doing the
		job of an OrderedSet class.
		'''
		for name in plot_names:
			if name == DEFAULT_TAG:
				self.expand_plot_names([CORE_TAG, VARIANT_TAG], name_dict)
			elif name == VARIANT_TAG:
				variant = self['metadata'].get('variant', '').upper()
				if variant in self.TAGS:
					self.expand_plot_names(self.TAGS[variant], name_dict)
			elif name in self.TAGS:
				self.expand_plot_names(self.TAGS[name], name_dict)
			elif name.endswith('.py'):
				name_dict[name] = True
			elif name.isupper():
				print(f'Ignoring an unknown analysis TAG: "{name}"')
			else:
				name_dict[name + '.py'] = True

	def list_plot_files(self, plot_names):
		'''List the plot module files (within self.MODULE_PATH) named by the
		given list of plot names, doing these transformations:

			* Default to the 'DEFAULT' tag
			* Expand all TAGS as defined by self.TAGS.
			* Append '.py' to filenames as needed.
			* Deduplicate entries but preserve the order.
		'''
		if not plot_names:
			plot_names = [DEFAULT_TAG]
		name_dict = OrderedDict()
		self.expand_plot_names(plot_names, name_dict)
		return list(name_dict.keys())

	def compile_images(self, file_list, extension='.png'):
		# type: (List[str], str) -> None
		"""
		Compiles images into a single file.

		Args:
			file_list: list of python analysis files that were run
			extension: plot filename extension, expected sub directory
				should be in SUB_DIRECTORIES

		Notes:
			- Only handles .png extension but other python packages should be able
			to handle vector based outputs (.pdf or .svg)
			- Will not handle all analysis plots produced if saved under a name
			other than the default filename
		"""

		# Identify images to compile
		sub_dir = SUB_DIRECTORIES.get(extension, '')
		output_dir = os.path.join(self['output_plots_directory'], sub_dir)
		output_bases = [self.plotter_args(f)[2] for f in file_list]

		# Load images and data
		image_files = [os.path.join(output_dir, base + extension) for base in output_bases]
		images = [Image.open(f) for f in image_files if os.path.exists(f)]
		if not images:
			return
		widths, heights = list(zip(*(i.size for i in images)))

		# Create and save compiled image
		compiled_image = Image.new('RGB', (max(widths), sum(heights)), (255, 255, 255))
		pos = 0
		for image in images:
			compiled_image.paste(image, (0, pos))
			pos += image.size[1]
		compiled_image.save(os.path.join(output_dir, 'compiled' + extension))

		for image in images:
			image.close()
		compiled_image.close()

	def run_task(self, fw_spec):
		start_real_sec = time.monotonic()
		print("\n{}: --- Starting {} ---".format(
			time.ctime(), type(self).__name__))

		plot_names = self.get("plot", [])
		fileList = self.list_plot_files(plot_names)

		self['output_filename_prefix'] = self.get('output_filename_prefix', '')
		self['metadata'] = data.expand_keyed_env_vars(self['metadata'])

		cpus = parallelization.cpus(self.get("cpus", 1))
		pool = None
		results = {}

		if cpus > 1:
			pool = parallelization.pool(cpus)

		# Set analysis paths from args
		input_dir = self.plotter_args('')[0]
		variant_paths = self.get('variant_paths')
		seed_paths = self.get('seed_paths')
		generation_paths = self.get('generation_paths')
		only_successful = self.get('only_successful', False)
		analysis_paths = AnalysisPaths(input_dir, **self.analysis_path_options)
		analysis_paths.update_cells(variant=variant_paths, seed=seed_paths,
			generation=generation_paths, only_successful=only_successful)

		exceptionFileList = []
		for f in fileList:
			try:
				mod = importlib.import_module(self.MODULE_PATH + '.' + f[:-3])
			except ModuleNotFoundError as e:
				print(f'Ignoring an unknown analysis class: {e}')
				continue
			except ImportError:
				traceback.print_exc()
				exceptionFileList.append(f)
				continue

			args = self.plotter_args(f)

			if pool:
				results[f] = pool.apply_async(run_plot, args=(mod.Plot, args, f, analysis_paths))
			else:
				print("{}: Running {}".format(time.ctime(), f))
				# noinspection PyBroadException
				try:
					mod.Plot.main(*args, analysis_paths=analysis_paths)
				except Exception:
					traceback.print_exc()
					exceptionFileList.append(f)

		if pool:
			pool.close()
			pool.join()
			for f, result in results.items():
				if not result.successful():
					exceptionFileList.append(f)

		if self.get('compile', False):
			print('{}: Compiling images'.format(time.ctime()))
			self.compile_images(fileList)

		elapsed_real_sec = time.monotonic() - start_real_sec

		duration = datetime.timedelta(seconds=elapsed_real_sec)
		if exceptionFileList:
			print('Completed analysis in {} with an exception in:'.format(duration))
			for f in exceptionFileList:
				print('\t{}'.format(f))
			raise RuntimeError('Error in analysis plot(s): {}'.format(', '.join(exceptionFileList)))
		else:
			print('Completed analysis in {}'.format(duration))


def run_plot(plot_class, args, name, analysis_paths):
	"""Run the given plot class in a Pool worker.
	Since this Firetask is running multiple plot classes in parallel, ask them
	to use just 1 CPU core each.
	"""
	try:
		print("{}: Running {}".format(time.ctime(), name))
		plot_class.main(*args, cpus=1, analysis_paths=analysis_paths)
	except KeyboardInterrupt:
		sys.exit(1)
	except Exception as e:
		traceback.print_exc()
		raise Exception(e)  # TODO: Return e so the caller can print it?
