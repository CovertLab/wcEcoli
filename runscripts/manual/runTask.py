import sys

from wholecell.fireworks.firetasks import InitRawDataTask
from wholecell.fireworks.firetasks import InitRawValidationDataTask
from wholecell.fireworks.firetasks import InitValidationDataTask
from wholecell.fireworks.firetasks import FitSimDataTask
from wholecell.fireworks.firetasks import VariantSimDataTask
from wholecell.fireworks.firetasks import SimulationTask
from wholecell.fireworks.firetasks import SimulationDaughterTask
from wholecell.fireworks.firetasks import AnalysisVariantTask
from wholecell.fireworks.firetasks import AnalysisCohortTask
from wholecell.fireworks.firetasks import AnalysisSingleTask
from wholecell.fireworks.firetasks import AnalysisMultiGenTask
from wholecell.fireworks.firetasks import BuildCausalityNetworkTask

def trim_leading_dashes(s):
	return s[2:]

if __name__ == '__main__':
	shell_args = sys.argv[1:]
	args = {}
	for arg in range(0, len(shell_args), 2):
		key = trim_leading_dashes(shell_args[arg])
		args[key] = shell_args[arg + 1]
	
	tasks = {
		'init_raw_data': InitRawDataTask,
		'init_raw_validation_data': InitRawValidationDataTask,
		'init_validation_data': InitValidationDataTask,
		'fit_sim_data': FitSimDataTask,
		'variant_sim_data': VariantSimDataTask,
		'simulation': SimulationTask,
		'simulation_daughter': SimulationDaughterTask,
		'analysis_variant': AnalysisVariantTask,
		'analysis_cohort': AnalysisCohortTask,
		'analysis_single': AnalysisSingleTask,
		'analysis_multigen': AnalysisMultiGenTask,
		'build_causality_network': BuildCausalityNetworkTask}

	task = tasks[args['task']]()
	for key, value in args.iteritems():
		setattr(task, key, value)

	task.run_task({})
