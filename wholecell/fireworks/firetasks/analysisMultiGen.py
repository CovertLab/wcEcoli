import cPickle
import time
import os

from fireworks import FireTaskBase, explicit_serialize
import models.ecoli.analysis.multigen
import importlib
import multiprocessing as mp

@explicit_serialize
class AnalysisMultiGenTask(FireTaskBase):

	_fw_name = "AnalysisMultiGenTask"
	required_params = [
		"input_seed_directory",
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