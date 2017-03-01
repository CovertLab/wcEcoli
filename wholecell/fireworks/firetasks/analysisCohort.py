"""
AnalysisCohortTask

Analyzes all cells, all seeds, all generations.

@author: Morgan Paull
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 11/30/2015
"""

import cPickle
import time
import os

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

	def run_task(self, fw_spec):
		return

def run_function(f, args, name):
	try:
		print "%s: Running %s" % (time.ctime(), name)
		f(*args)
	except KeyboardInterrupt:
		import sys; sys.exit(1)
