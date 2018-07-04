"""
AnalysisCohortTask

Analyzes all cells, all seeds, all generations.

@author: Morgan Paull
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 11/30/2015
"""

import time
import os
import traceback

from fireworks import FireTaskBase, explicit_serialize
import models.ecoli.analysis.cohort
import importlib
import multiprocessing as mp

@explicit_serialize
class AnalysisCohortTask(FireTaskBase):

	_fw_name = "AnalysisCohortTask"
	required_params = [
		"input_variant_directory",
		"input_sim_data",
		"input_validation_data",
		"output_plots_directory",
		"metadata",
	]
	optional_params = [
		"plots_to_run",  # absent or empty => run all active analysis plots
	]

	def run_task(self, fw_spec):

		startTime = time.time()
		print "\n%s: Running cohort analysis" % time.ctime(startTime)

		fileList = self.get("plots_to_run", [])
		if not fileList:
			fileList = models.ecoli.analysis.cohort.ACTIVE

		if "WC_ANALYZE_FAST" in os.environ:
			pool = mp.Pool(processes = 8)
			results = {}

		exceptionFileList = []
		for f in fileList:
			mod = importlib.import_module("models.ecoli.analysis.cohort." + f[:-3])
			args = (
				self["input_variant_directory"],
				self["output_plots_directory"],
				f[:-3],
				self["input_sim_data"],
				self['input_validation_data'],
				self["metadata"],
				)

			if "WC_ANALYZE_FAST" in os.environ:
				results[f] = pool.apply_async(run_plot, args=(mod.Plot, args, f))
			else:
				print "%s: Running %s" % (time.ctime(), f)
				try:
					mod.Plot.main(*args)
				except Exception:
					traceback.print_exc()
					exceptionFileList += [f]

		if "WC_ANALYZE_FAST" in os.environ:
			pool.close()
			pool.join()
			for f, result in results.items():
				if not result.successful():
					exceptionFileList += [f]

		timeTotal = time.time() - startTime

		if exceptionFileList:
			print "Completed cohort analysis in %s with an exception in:" % (time.strftime("%H:%M:%S", time.gmtime(timeTotal)))
			for file in exceptionFileList:
				print "\t%s" % file
			raise Exception("Error in cohort analysis")
		else:
			print "Completed cohort analysis in %s" % (time.strftime("%H:%M:%S", time.gmtime(timeTotal)))

def run_plot(plot_class, args, name):
	try:
		print "%s: Running %s" % (time.ctime(), name)
		plot_class.main(*args)
	except KeyboardInterrupt:
		import sys; sys.exit(1)
	except Exception as e:
		traceback.print_exc()
		raise Exception(e)
