# runs all variant analysis plots for a given sim with data from one location to another
# only need parent directory in wcEcoli_out as argument

import os
import sys
import cPickle

from wholecell.fireworks.firetasks.analysisVariant import AnalysisVariantTask

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

inputDirectory = path
outputDir = os.path.join(path, 'plotOut')
validationData = os.path.join(path, 'kb/validationData.cPickle')
with open(os.path.join(path, 'metadata', 'metadata.cPickle')) as f:
	metadata = cPickle.load(f)
metadata['analysis_type'] = 'variant'
metadata['total_variants'] = None

print('Variant analysis from {}'.format(inputDirectory))
task = AnalysisVariantTask(input_directory=inputDirectory, input_validation_data=validationData, output_plots_directory=outputDir, metadata=metadata)
task.run_task(None)
