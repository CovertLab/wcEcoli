"""
This script is intended to simply pull RNA and protein counts from a simulation output.

The data can be subsampled or not.

@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 12/21/2018
"""

from __future__ import absolute_import

import os
import cPickle
from datetime import datetime
import numpy as np
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import pandas as pd
import sys
from datetime import date
from functools import partial
from reconstruction import spreadsheets
import csv
from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from wholecell.io.tablereader import TableReader

#Using memory_profiler for optimization purposes
#from memory_profiler import profile

startTime = datetime.now()

#subsampling is currenty not avaialable.
#subsample = False
#subsample_degree = 100 #subsamples every 100 timepoints
print_runtime = True #Prints to screen the amount of time it took to run the analysis

DIALECT = "excel-tab"
JsonReader = partial(spreadsheets.JsonReader, dialect = DIALECT)
JsonWriter = partial(spreadsheets.JsonWriter, dialect = DIALECT)

def make_output_dirs(seed):
	"""
	Make output directory, by assuming the cwd is the base wcEcoli folder.
	"""
	output_dir = os.path.join('out', 'counts', 'wildtype_000000', 'count_out', seed, 'complex')
	print(output_dir)
	if not os.path.exists(output_dir):
		os.makedirs(output_dir)
	return output_dir

def parse_tsv(tsv_file):
	
#Takes in a tsv file, and creates a list of lists of the rows 
#contained within the TSV.
	
	tsv_list = []
	with open(tsv_file) as tsvfile:
		reader = JsonReader(row for row in tsvfile if not row.startswith('#'))
		fieldnames = reader.fieldnames
		for row in reader:
			tsv_list.append(row)

	return tsv_list, fieldnames

def load_simulation_data(sim_data_file):
	f = open(sim_data_file, 'rb')
	print(f)
	sim_data = cPickle.load(f)
	f.close()
	return sim_data

def import_proteins_rnas_tsv():
	'''
	Assumes the current directory is wcEcoli. 
	need to structure outputs at np.ndarray
	'''
	FLAT_DIR = os.path.join('reconstruction', 'ecoli', 'flat')
	COMPLEX_FILE = os.path.join(FLAT_DIR, 'complexationReactions.tsv')
	complex_info, fieldnames = parse_tsv(COMPLEX_FILE)

	return complex_info

def extract_protein_rna_ids(complex_data):
	'''
	The original code that pulls the protein and rna ids from sim_data
	puts it into an array. I am mimicing that action here. 

	The point of pulling in data directly from the flat files, instead 
	of from sim_data is that since we only need this data, we dont want to 
	pay money to download the sim_data file from the cloud.
	'''
	complex_ids = np.asarray([x['stoichiometry'][0]['molecule'] + '[' + str(x['stoichiometry'][0]['location']) + ']' for x in complex_data], dtype='S50')

	return complex_ids

#@profile
def gather_time_info(sim_out_dir, TableReader, time, time_eachGen, gen_num):
	time += TableReader(os.path.join(sim_out_dir, "Main")).readColumn("time").tolist()
	time_eachGen.append(TableReader(os.path.join(sim_out_dir, "Main")).readColumn("time").tolist()[0])
	num_timesteps = len(TableReader(os.path.join(sim_out_dir, "Main")).readColumn("time").tolist())
	gen_num = gen_num + [i+1]* num_timesteps
	return time, time_eachGen, num_timesteps, gen_num

def initialize_ids(bulkMolecules, complex_ids):
	molecule_ids = bulkMolecules.readAttribute("objectNames")
	complex_indices = np.array([molecule_ids.index(x) for x in complex_ids])
	return molecule_ids, complex_indices
#@profile
def extract_rna_protein_counts(bulkMolecules, complex_indices):
	complex_counts = bulkMolecules.readColumn("counts")[:, complex_indices]
	return complex_counts
#@profile
def save_counts_per_gen(i, complex_counts, output_dir):
	file_name_per_gen_complex = format(i, '02') + '_gen_data_complex.tsv'
	with open(os.path.join(output_dir, file_name_per_gen_complex), 'wb') as fp:
		np.savetxt(fp, complex_counts, '%s','\t')
	return

def save_seed_level_files(filenames, data_to_save, output_dir):
	for i in range(0, len(filenames)):
		with open(os.path.join(output_dir, filenames[i]), 'wb') as fp:
			np.savetxt(fp, data_to_save[i], '%s','\t')
	return
#assumes sys.argv[1] goes straight to the seed folder.
data_path = sys.argv[1]
seed_name = data_path.split('/')[-1]
output_dir = make_output_dirs(seed_name)

complex_info = import_proteins_rnas_tsv()
complex_ids = extract_protein_rna_ids(complex_info)


#get all cells
ap = AnalysisPaths(data_path, multi_gen_plot = True)
all_dir = ap.get_cells()

# For each generation
non_zero_sum_rna_counts_all_gens = []
non_zero_sum_protein_counts_all_gens = []
time = []
time_eachGen = []
gen_num = []

for i, sim_dir in enumerate(all_dir):
	print 'Start of ' + str(i) 
	print datetime.now() - startTime
	sim_out_dir = os.path.join(sim_dir, "simOut")
	print(sim_out_dir)
	time, time_eachGen, num_timesteps, gen_num = gather_time_info(sim_out_dir, TableReader, time, time_eachGen, gen_num)
	bulkMolecules = TableReader(os.path.join(sim_out_dir, "BulkMolecules"))
	if i == 0:
		molecule_ids, complex_indices = initialize_ids (bulkMolecules, complex_ids)
	complex_counts = extract_rna_protein_counts(bulkMolecules, complex_indices)
	save_counts_per_gen(i, complex_counts, output_dir)
	bulkMolecules.close()
	print 'End of ' + str(i) 
	print datetime.now() - startTime


'''
TO DO:rewrite subsampling to move it into the main function!
This was originally written to take in the full concatenated data set. Now need to break it up per gen.
'''
'''
if subsample:
	rna_counts = rna_counts[0::subsample_degree].copy()
	protein_counts = protein_counts[0::subsample_degree].copy()
	time = np.array(time)[0::subsample_degree].copy()
	gen_num = np.array(gen_num)[0::subsample_degree].copy()
	save_file_name = seed_name + '_multi_gen_rna_protein_counts_ids_subsampled.tsv'
else:
	save_file_name = seed_name + '_multi_gen_rna_protein_counts_ids.tsv'
'''	


#save the remainder of collected information:

save_file_name_ids_complex = 'ids_complex.tsv'
save_file_name_time = 'time_info.tsv'
save_file_name_gen = 'gen_info.tsv'

seed_level_filenames = [save_file_name_ids_complex, save_file_name_time, save_file_name_gen]
seed_level_data = [complex_ids, np.array(time), gen_num]
save_seed_level_files(seed_level_filenames, seed_level_data, output_dir)

if print_runtime:
	print datetime.now() - startTime 