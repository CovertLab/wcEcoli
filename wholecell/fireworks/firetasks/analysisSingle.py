import cPickle
import time
import os

from fireworks import FireTaskBase, explicit_serialize
import models.ecoli.analysis.single
import importlib
import multiprocessing as mp

@explicit_serialize
class AnalysisSingleTask(FireTaskBase):

	_fw_name = "AnalysisSingleTask"
	required_params = [
		"input_results_directory",
		"input_sim_data",
		"input_validation_data",
		"output_plots_directory",
		"metadata",
		]

	def run_task(self, fw_spec):
		return

def run_function(f, args, name):
	try:
		print "%s: Running %s" % (time.ctime(), name)
		f(*args)
	except KeyboardInterrupt:
		import sys; sys.exit(1)
