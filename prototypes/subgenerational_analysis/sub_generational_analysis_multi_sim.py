"""
This analysis script is intended as a scratch pad to pull out sub-generational gene expression info 
for a set of specified genes, as well as for further data exploration.

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
#from wholecell.containers.bulk_objects_container import BulkObjectsContainer
#from wholecell.analysis.analysis_tools import exportFigure
#from models.ecoli.analysis import multigenAnalysisPlot
#from wholecell.utils.sparkline import whitePadSparklineAxis
'''
#TODO:

change all names to having underscore from all camel case.

'''





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

def make_output_dirs(sim_data_path):
	"""
	Store outputs in E. coli output directory, near where simulation data 
	is being kept.

	Make the output directory if it does not already exist
	"""
	output_file_name = get_todays_date()
	base_output_dir = os.path.join(os.path.dirname(sim_data_path), 
		'subgenerational_analysis', output_file_name)
	output_dir = os.path.join(base_output_dir, 'outputs')
	if not os.path.exists(output_dir):
		os.makedirs(output_dir)
	return output_dir

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
		sim_data_file = os.path.join(sim_data_path, 'wildtype_000000', 'kb', 
			'simData_Modified.cPickle')
		if not os.path.isfile(sim_data_file):
			raise Exception, "There is something wrong with the data path " \
				"provided. Make sure you are pointing to the main sim output " \
				"directory."
	return sim_data_path, sim_data_path_gen, sim_data_file 

simulation_data_path, simulation_data_path_gen, simulation_data_file = get_input_paths(sys.argv)
sim_data = load_simulation_data(simulation_data_file)
output_dir = make_output_dirs(simulation_data_path)

#get all cells
ap = AnalysisPaths(simulation_data_path_gen, multi_gen_plot = True)
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
	#I think i initially called this _all cause its different than how heejo classified protein ids.
	bulkMolecules.close()

'''
	# Sum counts over timesteps
	sum_rna_counts = rna_counts.sum(axis = 0)
	sum_protein_counts = protein_counts_all.sum(axis = 0)

	# Flag where the sum is nonzero (True if nonzero, False if zero)
	non_zero_sum_rna_counts = sum_rna_counts != 0
	non_zero_sum_rna_counts_all_gens.append(non_zero_sum_rna_counts)
	non_zero_sum_protein_counts = sum_protein_counts != 0
	non_zero_sum_protein_counts_all_gens.append(non_zero_sum_protein_counts)

'''
#rna_ids_counts = np.vstack((rna_ids, rna_counts))
#protein_ids_counts = np.vstack((protein_ids, protein_counts))
#time_labeled = np.hstack((np.array('time'), np.array(time)))
#gen_labeled = np.hstack((np.array('gen'), gen_num))




# The file with all the data is too large to work with easily.
# Am just taking every 100 timepoints to make it easier to work with.

rna_ids_counts_subsampled = np.vstack((rna_ids, rna_counts[0::100].copy()))
protein_ids_counts_subsampled = np.vstack((protein_ids, protein_counts[0::100].copy()))
time_labeled_subsampled = np.hstack((np.array('time'), np.array(time)[0::100].copy()))
gen_labeled_subsampled = np.hstack((np.array('gen'), np.array(gen_num)[0::100].copy()))
#import ipdb; ipdb.set_trace()
time_gen_stack = np.vstack((gen_labeled_subsampled, time_labeled_subsampled)).transpose()

combined_rna_protein_counts = np.hstack((time_gen_stack, protein_ids_counts_subsampled, rna_ids_counts_subsampled))
'''
with open(os.path.join(output_dir, "multi_gen_rnacounts_ids.tsv"), 'wb') as fp:
	np.savetxt(fp, rna_ids_counts, '%s','\t')
with open(os.path.join(output_dir, "multi_gen_proteincounts_ids.tsv"), 'wb') as fp:
	np.savetxt(fp, protein_ids_counts, '%s','\t')
'''
with open(os.path.join(output_dir, "multi_gen_rna_protein_counts_ids.tsv"), 'wb') as fp:
	np.savetxt(fp, combined_rna_protein_counts, '%s','\t')

'''
# Average (mean) over generations
non_zero_sum_rna_counts_all_gens = np.array(non_zero_sum_rna_counts_all_gens)
avg_rna_counts = non_zero_sum_rna_counts_all_gens.mean(axis = 0)

# Average (mean) monomers over generations
non_zero_sum_protein_counts_all_gens = np.array(non_zero_sum_protein_counts_all_gens)
avg_protein_counts = non_zero_sum_protein_counts_all_gens.mean(axis = 0)

# Identify subgenerationally transcribed genes
subgen_rna_indices = np.where(np.logical_and(avg_rna_counts != 0., avg_rna_counts != 1.))[0]
subgen_rna_ids = rna_ids[subgen_rna_indices]
subgen_monomer_indices = [np.where(sim_data.process.translation.monomerData["rnaId"] == x)[0][0] for x in subgen_rna_ids]
subgen_monomer_ids = sim_data.process.translation.monomerData["id"][subgen_monomer_indices]

# Identify subgenerationally translated genes (not directly tied to subgenerational transcription)
subgen_monomer_indices_all = np.where(np.logical_and(avg_protein_counts != 0., avg_protein_counts != 1.))[0]
subgen_monomer_ids_all = protein_ids[subgen_monomer_indices_all]


with open(os.path.join(output_dir, "subgenerational_rna_ids.tsv"), 'wb') as fp:
	np.savetxt(fp, subgen_rna_ids, '%s','\t')
with open(os.path.join(output_dir, "subgenerational_protein_ids.tsv"), 'wb') as fp:
	np.savetxt(fp, subgen_monomer_ids_all, '%s','\t')

# plot mRNA and 
'''