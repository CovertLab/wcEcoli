"""
This script is intended to simply pull RNA and protein counts from a simulation output.

The data can be subsampled or not.

@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 12/21/2018
"""

from __future__ import absolute_import

import os
import cPickle
import numpy as np
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import pandas as pd
import sys
from datetime import date
from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from wholecell.io.tablereader import TableReader

subsample = False
subsample_degree = 100 #subsamples every 100 timepoints
'''
When data is stored on the cloud, need to still save outputs locally.
Take in this boolean then choose where to make the outputs folder
accordingly.
'''
data_on_cloud = True

def get_todays_date():
	"""
	Want to grab today's date so that data in output folders is not 
	overwritten (just in case its needed later)
	"""
	today = date.today()
	date_str = str(today.month) + str(today.day) + str(today.year)[2:]
	return date_str

def load_simulation_data(sim_data_file):
	f = open(sim_data_file, 'rb')
	sim_data = cPickle.load(f)
	f.close()
	return sim_data

def make_output_dirs(sim_data_path, data_on_cloud):
	"""
	If data is being store on Google Cloud: data_on_cloud = True
	Store outputs in the sub_gen_data sets on todays date.

	If data is being stored Locally: data_on_cloud = False
	Store outputs in E. coli output directory, near where simulation data 
	is being kept.

	Make the output directory if it does not already exist
	"""
	output_file_name = get_todays_date()
	if data_on_cloud:
		base_output_dir = os.path.join(os.path.dirname(os.getcwd()), 
			'data_sets', 'sub_gen_data_sets', output_file_name)
	else:
		base_output_dir = os.path.join(sim_data_path, 
			'subgenerational_analysis', output_file_name)
	output_dir = os.path.join(base_output_dir, 'outputs')
	#make output directory
	if not os.path.exists(output_dir):
		os.makedirs(output_dir)
	return output_dir


def get_input_paths(argument):
	# Check for input file.
	if len(argument) < 2:
		raise Exception, "Oops you didn't provide a path to the data!"
	else:
		sim_data_path = argument[-1]
		#check that parent path is in fact a directory:
		if not os.path.isdir(sim_data_path):
			raise Exception, "input_data_parent_dir does not currently " \
				"exist as a directory"
		sim_data_file = os.path.join(sim_data_path, 'wildtype_000000', 'kb', 'simData_Modified.cPickle')
		if not os.path.isfile(sim_data_file):
			raise Exception, "There is something wrong with the data path " \
				"provided. Make sure you are pointing to the main sim output " \
				"directory."
		sim_data_path_all_sims = os.path.join(sim_data_path, 'wildtype_000000')
		#need to get paths for all the sims:
		current_sub_dir = os.listdir(sim_data_path_all_sims)
		sim_data_folders = [direct for direct in current_sub_dir if direct.isdigit()]
		sim_data_path_gen = []
		for sim_data_folder in sim_data_folders:
			sim_data_path_gen.append(os.path.join(sim_data_path_all_sims, sim_data_folder))
	return sim_data_path, sim_data_path_gen, sim_data_file 

'''
def get_input_paths(argument): 
		# Check for input file.
		if len(argument) < 2:
			raise Exception, "Oops you didn't provide a path to the data!"
		else:
			sim_data_path = argument[-1]
			sim_data_path_gen = os.path.join(sim_data_path, 'wildtype_000000', '000000')
			#check that parent path is in fact a directory:
			if not os.path.isdir(sim_data_path):
				raise Exception, "input_data_parent_dir does not currently " \
					"exist as a directory"
			
		return sim_data_path, sim_data_path_gen, sim_data_file 
'''
simulation_data_path, simulation_data_path_gen, simulation_data_file = get_input_paths(sys.argv)
sim_data = load_simulation_data(simulation_data_file)
output_dir = make_output_dirs(simulation_data_path, data_on_cloud)

#get all cells
for data_path in simulation_data_path_gen:
	seed_name = data_path.split('/')[-1]
	#import pdb; pdb.set_trace()
	ap = AnalysisPaths(data_path, multi_gen_plot = True)
	all_dir = ap.get_cells()


	rna_ids = sim_data.process.transcription.rnaData["id"]
	gene_ids = sim_data.process.transcription.rnaData["geneId"]
	protein_ids = sim_data.process.translation.monomerData["id"]

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