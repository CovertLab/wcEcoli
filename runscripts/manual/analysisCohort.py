# runs all cohort analysis plots for a given sim with data from one location to another
# only need parent directory in wcEcoli_out as argument

import os
import sys
import cPickle

from wholecell.fireworks.firetasks.analysisCohort import AnalysisCohortTask

BASE_DIR = '/scratch/users/thorst/wcEcoli_out/'

if len(sys.argv) < 2:
	for folder in sorted(os.listdir(BASE_DIR), reverse = True):
		if folder[0].isdigit():
			path = os.path.join(BASE_DIR, folder)
			break
else:
	arg = sys.argv[1]
	if arg.startswith('out/'):
		arg = arg[4:]
	path = os.path.join(BASE_DIR, arg)

for folder in os.listdir(path):
	if '_' in folder:
		variant = folder

#variant = 'condition_000001'

inputVariantDirectory = os.path.join(path, variant)
outputDir = os.path.join(path, variant, 'plotOut')
simData = os.path.join(path, variant, 'kb/simData_Modified.cPickle')
validationData = os.path.join(path, 'kb/validationData.cPickle')
with open(os.path.join(path, 'metadata', 'metadata.cPickle')) as f:
	metadata = cPickle.load(f)
metadata['analysis_type'] = 'cohort'
metadata['variant_function'] = variant
metadata['variant_index'] = None

print('Cohort analysis from {}'.format(inputVariantDirectory))
task = AnalysisCohortTask(input_variant_directory=inputVariantDirectory, input_sim_data=simData, input_validation_data=validationData, output_plots_directory=outputDir, metadata=metadata)
task.run_task(None)