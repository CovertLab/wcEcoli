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
		return

def run_function(f, args, name):
	try:
		print "%s: Running %s" % (time.ctime(), name)
		f(*args)
	except KeyboardInterrupt:
		import sys; sys.exit(1)