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
import csv
DIALECT = "excel-tab"
JsonReader = partial(spreadsheets.JsonReader, dialect = DIALECT)
JsonWriter = partial(spreadsheets.JsonWriter, dialect = DIALECT)

FLAT_DIR = os.path.join('reconstruction', 'ecoli', 'flat')
RNA_FILE = os.path.join(FLAT_DIR, 'rnas.tsv')
TU_FILE = os.path.join(FLAT_DIR, 'operon_rnas.tsv')
RNA_SEQ_FILE = os.path.join(FLAT_DIR, 'rna_seq_data', 'rnaseq_rsem_tpm_mean.tsv')
CONDITION = 'M9 Glucose minus AAs'
SPLIT_DELIMITER = '_'

output_tu_counts = os.path.join(FLAT_DIR, "transcription_units.tsv")
output_gene_tu_matrix = os.path.join(FLAT_DIR, "gene_to_tu_matrix.tsv")

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

def map_genes_to_tu(gene, index, reverse_index, gene_to_tu_matrix):

	return gene_to_tu_matrix


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

	rnas_gene_order = [row['geneId'] for row in rna_info]
	#import pdb; pdb.set_trace()
	reverse_index = {
		row['geneId']: gene_index 
		for gene_index, row in enumerate(rna_info)}


	for index, tu in enumerate(tu_info):
		if not SPLIT_DELIMITER in tu['geneId']:
			gene = tu['geneId']
			gene_index = reverse_index[gene]
			gene_to_tu_matrix[gene_index, index] = 1
		else:
			genes_in_tu = re.split(SPLIT_DELIMITER, tu['geneId'])
			for gene in genes_in_tu:
				gene_index = reverse_index[gene]
				gene_to_tu_matrix[gene_index, index] = 1

	return gene_to_tu_matrix, rnas_gene_order

def create_rnaseq_count_vector(rnas_gene_order):
	'''
	The counts vector is not in the same order as the Gene_TU_Matrix.
	Need to reoder and pull out count information. 

	Gathers information needed based on the condition the model is 
	being run in.
	'''
	rna_seq_data_all_cond = parse_tsv(RNA_SEQ_FILE)

	rna_seq_data_index = {
		row['Gene']: row[CONDITION] 
		for row in rna_seq_data_all_cond}

	rna_seq_counts_vector = [
		rna_seq_data_index[gene] 
		for gene in rnas_gene_order]

	return rna_seq_counts_vector

def create_tu_counts_vector(gene_tu_matrix, rna_seq_counts_vector, tu_info):
	tu_counts_vector = np.linalg.lstsq(gene_tu_matrix, rna_seq_counts_vector)[0]

	#import pdb; pdb.set_trace()	
	
	tu_gene_order = [row['geneId'] for row in tu_info]
	tu_genes_counts = []

	for i in range(0, len(tu_gene_order)):
		tu_genes_counts.append({'tu_id': tu_gene_order[i], 'tu_count': tu_counts_vector[i]})

	return tu_genes_counts


def calculate_tu_counts_vector():
	rna_info = parse_tsv(RNA_FILE)
	tu_info = parse_tsv(TU_FILE)
	gene_tu_matrix, rnas_gene_order = create_gene_to_tu_matrix(rna_info, tu_info)
	rna_seq_counts_vector = create_rnaseq_count_vector(rnas_gene_order)
	tu_counts = create_tu_counts_vector(gene_tu_matrix, rna_seq_counts_vector, tu_info)
	fieldnames = ['tu_id', 'tu_count']



	with open(output_tu_counts, "w") as f:
		writer = JsonWriter(f, fieldnames)
		writer.writeheader()
		for tu_count in tu_counts:
			writer.writerow(tu_count)

	with open(output_gene_tu_matrix, "w") as f:
		writer = csv.writer(f, delimiter=' ')
		for row in gene_tu_matrix:
			writer.writerow(row)

	#return tu_counts


if __name__ == "__main__":
	calculate_tu_counts_vector()
