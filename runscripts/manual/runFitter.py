"""
Run the Fitter. The output goes into out/manual/.
Set PYTHONPATH when running this.
"""

from __future__ import absolute_import
from __future__ import division

import cPickle
import datetime
import time
import os

from reconstruction.ecoli.knowledge_base_raw import KnowledgeBaseEcoli
from reconstruction.ecoli.fit_sim_data_1 import fitSimData_1
import wholecell
from wholecell.utils import filepath

DEBUG = False
ROOT_DIR = os.path.dirname(os.path.dirname(wholecell.__file__))
print "ROOT_DIR = {}".format(ROOT_DIR)

location = filepath.makedirs(ROOT_DIR, "out", "manual", "kb")
rawDataFile = os.path.join(location, "raw_data.cPickle")
simDataFile = os.path.join(location, "simData_Fit_1.cPickle")
print "simDataFile = {}".format(simDataFile)

start_sec = time.clock()

print "{}: Loading raw".format(time.ctime())
raw_data = KnowledgeBaseEcoli()

print "{}: Fitting".format(time.ctime())
sim_data = fitSimData_1(raw_data, debug=DEBUG)
print "{}: Done fitting".format(time.ctime())

with open(rawDataFile, "wb") as f:
    cPickle.dump(raw_data, f)
with open(simDataFile, "wb") as f:
    cPickle.dump(sim_data, f)
print "{}: Done writing parameter files".format(time.ctime())

end_sec = time.clock()
elapsed = datetime.timedelta(seconds=end_sec - start_sec)
print "Ran the Fitter in {}h {}m {}s total".format(*str(elapsed).split(':'))
