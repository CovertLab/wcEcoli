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

startTime = datetime.now()

subsample = False
subsample_degree = 100 #subsamples every 100 timepoints
print_runtime = True #Prints to screen the amount of time it took to run the analysis

DIALECT = "excel-tab"
JsonReader = partial(spreadsheets.JsonReader, dialect = DIALECT)
JsonWriter = partial(spreadsheets.JsonWriter, dialect = DIALECT)

def make_output_dirs():
	"""
	Make output directory, by assuming the cwd is the base wcEcoli folder.
	"""
	output_dir = os.path.join('out', 'counts', 'wildtype_000000', 'count_out')
	
	if not os.path.exists(output_dir):
		os.makedirs(output_dir)
	
	return output_dir

def parse_tsv(tsv_file):
	
#Takes in a tsv file, and creates a list of lists of the rows 
#contained within the TSV.
	
	tsv_list = []
	with open(tsv_file) as tsvfile:
		reader = JsonReader(tsvfile)
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
	RNA_FILE = os.path.join(FLAT_DIR, 'rnas.tsv')
	PROTEIN_FILE = os.path.join(FLAT_DIR, 'proteins.tsv')
	protein_info, fieldnames = parse_tsv(PROTEIN_FILE)
	rna_info, fieldnames = parse_tsv(RNA_FILE)

	return protein_info, rna_info

def extract_protein_rna_ids(protein_data, rna_data):
	'''
	The original code that pulls the protein and rna ids from sim_data
	puts it into an array. I am mimicing that action here. 

	The point of pulling in data directly from the flat files, instead 
	of from sim_data is that since we only need this data, we dont want to 
	pay money to download the sim_data file from the cloud.
	'''
	protein_ids = np.asarray([x['id'] + '[' + str(x['location'][0]) + ']' for x in protein_data], dtype = 'S50')
	rna_ids = np.asarray([x['id'] + '[' + str(x['location'][0]) + ']' for x in rna_data], dtype='S50')

	return protein_ids, rna_ids


#assumes sys.argv[1] goes straight to the seed folder.
data_path = sys.argv[1]
output_dir = make_output_dirs()
seed_name = data_path.split('/')[-1]

protein_info, rna_info = import_proteins_rnas_tsv()
protein_ids, rna_ids = extract_protein_rna_ids(protein_info, rna_info)
#
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
	sim_out_dir = os.path.join(sim_dir, "simOut")
	
	time += TableReader(os.path.join(sim_out_dir, "Main")).readColumn("time").tolist()
	time_eachGen.append(TableReader(os.path.join(sim_out_dir, "Main")).readColumn("time").tolist()[0])
	num_timesteps = len(TableReader(os.path.join(sim_out_dir, "Main")).readColumn("time").tolist())
	gen_num = gen_num + [i+1]* num_timesteps

	# Read counts of transcripts
	bulkMolecules = TableReader(os.path.join(sim_out_dir, "BulkMolecules"))
	
	if i == 0:
		molecule_ids = bulkMolecules.readAttribute("objectNames")
		#import pdb; pdb.set_trace()
		rna_indices = np.array([molecule_ids.index(x) for x in rna_ids])
		protein_indices = np.array([molecule_ids.index(x) for x in protein_ids])
	try:
		rna_counts
	except NameError:
		rna_counts = None

	if rna_counts is None:
		rna_counts = bulkMolecules.readColumn("counts")[:, rna_indices]
		protein_counts = bulkMolecules.readColumn("counts")[:, protein_indices]
	else:
		rna_counts = np.vstack((rna_counts, bulkMolecules.readColumn("counts")[:, rna_indices]))
		protein_counts = np.vstack((protein_counts, bulkMolecules.readColumn("counts")[:, protein_indices]))
	bulkMolecules.close()

if subsample:
	rna_counts = rna_counts[0::subsample_degree].copy()
	protein_counts = protein_counts[0::subsample_degree].copy()
	time = np.array(time)[0::subsample_degree].copy()
	gen_num = np.array(gen_num)[0::subsample_degree].copy()
	save_file_name = seed_name + '_multi_gen_rna_protein_counts_ids_subsampled.tsv'
else:
	save_file_name = seed_name + '_multi_gen_rna_protein_counts_ids.tsv'

rna_ids_counts = np.vstack((rna_ids, rna_counts))
protein_ids_counts = np.vstack((protein_ids, protein_counts))
time_labeled = np.hstack((np.array('time'), np.array(time)))
gen_labeled = np.hstack((np.array('gen'), gen_num))
#import ipdb; ipdb.set_trace()
time_gen_stack = np.vstack((gen_labeled, time_labeled)).transpose()
combined_rna_protein_counts = np.hstack((time_gen_stack, protein_ids_counts, rna_ids_counts))

with open(os.path.join(output_dir, save_file_name), 'wb') as fp:
	np.savetxt(fp, combined_rna_protein_counts, '%s','\t')

if print_runtime:
	print datetime.now() - startTime 