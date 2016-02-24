"""
AnalysisVariantTask

Analyzes across variants. Has access to all cells in the entire simulation run.

@author: Morgan Paull
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 1/06/2015
"""

import cPickle
import time
import os

from fireworks import FireTaskBase, explicit_serialize
import models.ecoli.analysis.variant
import importlib
import multiprocessing as mp

from wholecell.analysis.analysis_firetask_utils import run_specific_order, run_function

@explicit_serialize
class AnalysisVariantTask(FireTaskBase):

	_fw_name = "AnalysisVariantTask"
	required_params = [
		"input_directory",
		"input_validation_data",
		"output_plots_directory",
		"metadata",
		]

	def run_task(self, fw_spec):

		startTime = time.time()
		print "%s: Running variant analysis" % time.ctime(startTime)

		directory = os.path.dirname(models.ecoli.analysis.variant.__file__)

		# Run analysis scripts in order of modification, most recently edited first
		fileList = os.listdir(directory)
		fileList.sort(key=lambda x: os.stat(os.path.join(directory, x)).st_mtime, reverse=True)

		# Run files in runlast.txt after others
		if "runlast.txt" in fileList:
			fileList = run_specific_order(directory, fileList, "runlast.txt", position=-1, verbose=False)
		# Run files in runfirst.txt before others
		if "runfirst.txt" in fileList:
			fileList = run_specific_order(directory, fileList, "runfirst.txt", position=0, verbose=False)


		if "WC_ANALYZE_FAST" in os.environ:
			pool = mp.Pool(processes = 8)

		for f in fileList:
			if f.endswith(".pyc") or f == "__init__.py":
				continue

			mod = importlib.import_module("models.ecoli.analysis.variant." + f[:-3])
			args = (
				self["input_directory"],
				self["output_plots_directory"],
				f[:-3],
				self['input_validation_data'],
				self["metadata"]
				)

			if "WC_ANALYZE_FAST" in os.environ:
				pool.apply_async(run_function, args = (mod.main, args, f))
			else:
				print "%s: Running %s" % (time.ctime(), f)
				mod.main(*args)

		if "WC_ANALYZE_FAST" in os.environ:
			pool.close()
			pool.join()
		timeTotal = time.time() - startTime
		print "Completed variant analysis in %s" % (time.strftime("%H:%M:%S", time.gmtime(timeTotal)))
