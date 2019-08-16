import json
import sys

SEEDS = [0, 1, 2, 3, 4, 5, 6, 7, 8, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 29, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 50, 51, 52, 53, 54, 55, 56, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 75, 76, 77, 78, 79, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 95, 96, 97, 98, 99, 128, 300, 301, 302, 303, 304, 305, 306, 307, 308, 309]

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
	num_gens = int(sys.argv[1])
	commands =  [analysis_command(num_gens)]
	steps = [analysis_step(seed, num_gens) for seed in SEEDS]
	json.dump(commands, open('out/analysis-commands.json', 'w'), sort_keys=True, indent=4, separators=(',', ': '))
	json.dump(steps, open('out/analysis-steps.json', 'w'), sort_keys=True, indent=4, separators=(',', ': '))
