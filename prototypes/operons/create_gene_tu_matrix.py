'''
The goal of this script is to take in a file containing all the possible TU's (this is most like what we would expect from the sequencing data) and create a gene to TU matrix. 
For now extract the TU list from the rna_operons.tsv file. Later this will be provided. 
Extract the gene list from rnas.tsv.
'''

import numpy as np
import os
import re
from reconstruction import spreadsheets
from functools import partial
DIALECT = "excel-tab"
JsonReader = partial(spreadsheets.JsonReader, dialect = DIALECT)

def parse_tsv(tsv_file):
	'''
	Takes in a tsv file, and creates a list of dicts of the rows 
	contained within the TSV.
	'''
	tsv_list = []
	with open(tsv_file) as tsvfile:
		reader = JsonReader(tsvfile)
		for row in reader:
			tsv_list.append(row)
	return tsv_list

def find_tu_indices(tu_info):
	'''
	Input:
	A list of dicts containing all the information found in the TU_FILE.
	Do:
	Find all the TU's (based on the presence of the '_' dividng them.
	Retreive the 'index', based on count down the file.
	Return:
	A list containin the indices of all the TUs.
	Preserves order.
	TODO: Add another column to the operon_rnas.tsv file to assign a type as
	a TU.
	'''
	#count = 0
	#tu_ind = []
	return [
		count 
		for count, row in enumerate(tu_info) 
		if SPLIT_DELIMITER in row['geneId']]
	'''
	for row in tu_info:
		gene_id = row['geneId']
		if SPLIT_DELIMITER in gene_id:
			tu_ind.append(count)
		count += 1
	return tu_ind
	'''

def create_gene_to_tu_matrix(rna_info, tu_info):
	'''
	Input:
	Parsed rna and tu data files. 
	Do:
	Parse tsv files then create a gene to tu matrix mapping genes in each TU
	to its partner RNA.
	Return:
	Sparse numpy matrix, mapping TU to rna's. 0 = no mapping; 1 = mapping.

	TODO: Add a list containing the geneIds in a TU, can pull this in to
	genes_in_tu.
	'''

	num_rnas = len(rna_info)
	num_tus = len(tu_info)

	gene_to_tu_matrix = np.zeros((num_rnas, num_tus))
	tu_index = find_tu_indices(tu_info)

	rnas_gene_order = [row['geneId'] for row in rna_info]

	reverse_index = {
		row['geneId']: gene_index 
		for gene_index, row in enumerate(rna_info)}

	for index in tu_index:
		genes_in_tu = re.split(SPLIT_DELIMITER, tu_info[index]['geneId'])
		for gene in genes_in_tu:
			gene_index = reverse_index[gene]
			gene_to_tu_matrix[gene_index, index] = 1
	import ipdb; ipdb.set_trace()
	'''
	for gene_index, row in enumerate(rna_info):
		gene_id = row['geneId']
		for index in tu_index:
			genes_in_tu = re.split(SPLIT_DELIMITER, tu_info[index]['geneId'])
			for gene in genes_in_tu:
				if gene == gene_id:
					gene_to_tu_matrix[gene_index, index] = 1
	'''
	return gene_to_tu_matrix, rnas_gene_order

def create_rnaseq_count_vector(rnas_gene_order):
	'''
	The counts vector is not in the same order as the Gene_TU_Matrix.
	Need to reoder and pull out count information. 

	Gathers information needed based on the condition the model is 
	being run in.
	'''
	rna_seq_data_all_cond = parse_tsv(RNA_SEQ_FILE)
	#rna_seq_gene_id_list = []
	rna_seq_counts_vector = []
	#have this nexted for loop structure bc, I am trying to preserve the 
	#gene order from rnas.tsv so that the matrix and the vector match up.

	rna_seq_data_index = {
		row['Gene']: row[CONDITION] 
		for row in rna_seq_data_all_cond}
		
	rna_seq_counts_vector = [
		rna_seq_data_index[gene] 
		for gene in rnas_gene_order]
	
	'''
	for gene in rnas_gene_order:
		rna_seq_counts_vector.append(rna_seq_data_index[gene])
		
		for row in rna_seq_data_all_cond:
			if row['Gene'] == gene:
				#rna_seq_gene_id_list.append(row['Gene'])
				rna_seq_counts_vector.append(row[CONDITION])
	'''	

	return rna_seq_counts_vector


#parent_directory = os.path.dirname(os.path.dirname(os.getcwd()))
FLAT_DIR = os.path.join('reconstruction', 'ecoli', 'flat')
RNA_FILE = os.path.join(FLAT_DIR, 'rnas.tsv')
TU_FILE = os.path.join(FLAT_DIR, 'operon_rnas_2.tsv')
RNA_SEQ_FILE = os.path.join(FLAT_DIR, 'rna_seq_data', 'rnaseq_rsem_tpm_mean.tsv')
CONDITION = 'M9 Glucose minus AAs'
SPLIT_DELIMITER = '_'

rna_info = parse_tsv(RNA_FILE)
tu_info = parse_tsv(TU_FILE)
gene_tu_matrix, rnas_gene_order = create_gene_to_tu_matrix(rna_info, tu_info)
rna_seq_ids, rna_seq_counts_vector = create_rnaseq_count_vector(rnas_gene_order)
#returns 0's for anything that is not in a TU.
#do i need to make a new vector contining the rna seq counts with tu_counts_soln?
tu_counts_vector = np.linalg.lstsq(gene_tu_matrix, rna_seq_counts_vector)[0]


