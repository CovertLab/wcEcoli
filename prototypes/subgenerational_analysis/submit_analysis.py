import json
import sys
def analysis_command(num_gens):
	inputs = {}
	for generation in range(num_gens):
		generation_key = format(generation, '06')
		inputs['counts_{}'.format(generation_key)] = '/wcEcoli/out/counts/wildtype_000000/{{seed}}/generation_'+generation_key+'/000000/simOut/BulkMolecules/counts'
		inputs['counts_attr_{}'.format(generation_key)] = '/wcEcoli/out/counts/wildtype_000000/{{seed}}/generation_'+generation_key+'/000000/simOut/BulkMolecules/attributes.json'
		inputs['main_{}'.format(generation_key)] = '/wcEcoli/out/counts/wildtype_000000/{{seed}}/generation_'+generation_key+'/000000/simOut/Main/time'
		inputs['main_attr_{}'.format(generation_key)] = '/wcEcoli/out/counts/wildtype_000000/{{seed}}/generation_'+generation_key+'/000000/simOut/Main/attributes.json'

	command = {
		'name': 'pull_rna_protein_counts_cloud',
		'vars': {
			'seed': '000000'},
		'command': [
			'python', 
			'-u', 
			'prototypes/subgenerational_analysis/pull_rna_protein_counts_cloud.py', 
			'out/counts/wildtype_000000/{{seed}}'],
		'image': 'gcr.io/allen-discovery-center-mcovert/mialydefelice-wcm-code:latest',
		'outputs': {
			'output': '/wcEcoli/out/counts/wildtype_000000/count_out/{{seed}}_multi_gen_rna_protein_counts_ids.tsv'},
		'inputs': inputs}
	return command

def analysis_step(seed, num_gens):
	seed_key = format(seed, '06')
	inputs = {}
	for generation in range(num_gens):
		generation_key = format(generation, '06')
		inputs['counts_{}'.format(generation_key)] = 'sisyphus:data/mialydefelice/20190812.122845__Test_32gen_100seeds_basal/wildtype_000000/'+seed_key+'/generation_'+generation_key+'/000000/simOut/BulkMolecules/counts'
		inputs['counts_attr_{}'.format(generation_key)] = 'sisyphus:data/mialydefelice/20190812.122845__Test_32gen_100seeds_basal/wildtype_000000/'+seed_key+'/generation_'+generation_key+'/000000/simOut/BulkMolecules/attributes.json'
		inputs['main_{}'.format(generation_key)] = 'sisyphus:data/mialydefelice/20190812.122845__Test_32gen_100seeds_basal/wildtype_000000/'+seed_key+'/generation_'+generation_key+'/000000/simOut/Main/time'
		inputs['main_attr_{}'.format(generation_key)] = 'sisyphus:data/mialydefelice/20190812.122845__Test_32gen_100seeds_basal/wildtype_000000/'+seed_key+'/generation_'+generation_key+'/000000/simOut/Main/attributes.json'
	step = {
		'name': 'pull_rna_protein_counts_cloud_{}'.format(seed_key),
		'command': 'pull_rna_protein_counts_cloud',
		'vars': {
			'seed': seed_key},
		'outputs': {
			'output': 'sisyphus:data/mialydefelice/20190812.122845__Test_32gen_100seeds_basal/out/{}_multi_gen_rna_protein_counts_ids.tsv'.format(seed_key)},
		'inputs': inputs}
	return step

if __name__ == '__main__':
	seeds = [0]
	num_gens = int(sys.argv[1])
	commands =  [analysis_command(num_gens)]
	steps = [analysis_step(seed, num_gens) for seed in seeds]
	json.dump(commands, open('out/analysis-commands.json', 'w'), sort_keys=True, indent=4, separators=(',', ': '))
	json.dump(steps, open('out/analysis-steps.json', 'w'), sort_keys=True, indent=4, separators=(',', ': '))
