# runs all single analysis plots for a given sim
# if no folder is given it runs for most recent in wcEcoli_out otherwise picks first arg passed

import os
import sys
import cPickle

from wholecell.fireworks.firetasks.analysisSingle import AnalysisSingleTask

BASE_DIR = '/scratch/users/thorst/wcEcoli_out/'
SEED = '000000'
GEN = 'generation_000000'
DAUGHTER = '000000'
DIRS = os.path.join(SEED, GEN, DAUGHTER)

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

# variant = 'condition_000000'
# variant = 'lambdaWeight_000009'

resultsDir = os.path.join(path, variant, DIRS, 'simOut')
outputDir = os.path.join(path, variant, DIRS, 'plotOut')
simData = os.path.join(path, variant, 'kb/simData_Modified.cPickle')
validationData = os.path.join(path, 'kb/validationData.cPickle')
with open(os.path.join(path, 'metadata', 'metadata.cPickle')) as f:
	metadata = cPickle.load(f)
metadata['analysis_type'] = 'single'
metadata['variant_function'] = variant
metadata['variant_index'] = None
metadata['seed'] = SEED
metadata['gen'] = GEN

print('Single analysis from {}'.format(resultsDir))
task = AnalysisSingleTask(input_results_directory=resultsDir, input_sim_data=simData, input_validation_data=validationData, output_plots_directory=outputDir, metadata=metadata)
task.run_task(None)
