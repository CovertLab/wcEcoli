"""
Run all multigen analysis plots for a given sim.
Set PYTHONPATH when running this.

Arguments:
	location	(string, default) run_name, a subdirectory of "out/". Defaults
				to the highest timestamped subdirectory name, else to "manual".

Set the environment variable WC_ANALYZE_FAST to run multiple analysis scripts
in parallel processes.
"""

from __future__ import absolute_import
from __future__ import division

import importlib
import multiprocessing as mp
import datetime
import time
import os
import sys

import models.ecoli.analysis.multigen
import wholecell
from wholecell.utils import filepath

# On Sherlock, "out" is usually a symlink to "/scratch/users/$USER/wcEcoli_out/".
ROOT_DIR = os.path.dirname(os.path.dirname(wholecell.__file__))
BASE_DIR = os.path.join(ROOT_DIR, "out")
DIR = "000000"

if len(sys.argv) < 2:
	path = os.path.join(BASE_DIR, "manual")
	for folder in sorted(os.listdir(BASE_DIR), reverse=True):
		if os.path.isdir(folder) and folder[0].isdigit():
			path = os.path.join(BASE_DIR, folder)
			break
else:
	arg = sys.argv[1]
	if arg.startswith("out/"):
		arg = arg[4:]
	path = os.path.join(BASE_DIR, arg)

variant = "basal"
for folder in os.listdir(path):
	if os.path.isdir(folder) and "_" in folder:
		variant = folder

# variant = "wildtype_000000"
# variant = "condition_000000"
# variant = "geneKnockdown_000030"

if variant == "basal":
	inputSeedDirectory = os.path.join(path)  # TODO
	outputDir = filepath.makedirs(path, variant, "plotOut")
	simData = os.path.join(path, "kb/simData_Fit_1.cPickle")  # TODO
	validationData = os.path.join(path, "kb/validationData.cPickle")  # TODO
else:
	inputSeedDirectory = os.path.join(path, variant, DIR)
	outputDir = filepath.makedirs(path, variant, DIR, "plotOut")
	simData = os.path.join(path, variant, "kb/simData_Modified.cPickle")
	validationData = os.path.join(path, "kb/validationData.cPickle")

metadata = ""


start_sec = time.clock()
print "%s: Running multigen simulation analysis for %s" % (time.ctime(), inputSeedDirectory)

directory = os.path.dirname(models.ecoli.analysis.multigen.__file__)

# Run analysis scripts in order of modification, most recently edited first
fileList = os.listdir(directory)
fileList.sort(key=lambda x: os.stat(os.path.join(directory, x)).st_mtime, reverse=True)

def run_function(analysis_function, args, name):
	analysis_function(*args)

pool = None
if "WC_ANALYZE_FAST" in os.environ:
	pool = mp.Pool(processes=8)

for f in fileList:
	if not f.endswith(".py") or f == "__init__.py":
		continue

	mod = importlib.import_module("models.ecoli.analysis.multigen." + f[:-3])
	args = (
		inputSeedDirectory,
		outputDir,
		f[:-3],
		simData,
		validationData,
		metadata,
		)

	if pool:
		pool.apply_async(run_function, args=(mod.main, args, f))
	else:
		print "%s: Running %s" % (time.ctime(), f)
		mod.main(*args)

if pool:
	pool.close()
	pool.join()

end_sec = time.clock()
elapsed = datetime.timedelta(seconds=end_sec - start_sec)
print "Completed multiple generation analysis in {}h {}m {}s total".format(
	*str(elapsed).split(':'))
