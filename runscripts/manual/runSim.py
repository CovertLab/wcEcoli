"""
Run a simulation.
Only runs first gen at the moment.
Set PYTHONPATH when running this.

Arguments:
	run_name	(string, "manual") outputDir subdirectory of "out/".
"""

from __future__ import absolute_import
from __future__ import division

import datetime
import time
import sys
import os

from models.ecoli.sim.simulation import EcoliSimulation
import wholecell

options = {}
# options['lengthSec'] = 30

ROOT_DIR = os.path.dirname(os.path.dirname(wholecell.__file__))
BASE_DIR = os.path.join(ROOT_DIR, "out")

if len(sys.argv) > 1:
	run_name = sys.argv[1]
else:
	run_name = "manual"

location = os.path.join(BASE_DIR, run_name)

# TODO: Prefer <location>/<variant>/kb/simData_Modified.cPickle
# TODO: What about <location>/kb/simData_Most_Fit.cPickle?
simDataFile = os.path.join(location, "kb", "simData_Fit_1.cPickle")
if not os.path.exists(simDataFile):
	raise IOError("Missing '{}'.  Run the fitter.".format(simDataFile))

options["simDataLocation"] = simDataFile
options["outputDir"] = os.path.join(location, "sim")

start_sec = time.clock()

print "{}: Running simulation".format(time.ctime())
print "Options: {}\n".format(options)

sim = EcoliSimulation(**options)
sim.run()

print "{}: Done running simulation".format(time.ctime())

end_sec = time.clock()
elapsed = datetime.timedelta(seconds=end_sec - start_sec)
print "Ran a simulation in {}h {}m {}s total".format(*str(elapsed).split(':'))
