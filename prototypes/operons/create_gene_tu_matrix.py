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
	Takes in a tsv file, and creates a list of lists of the rows 
	contained within the TSV.
	'''
	tsv_list = []
	with open(tsv_file) as tsvfile:
		reader = JsonReader(tsvfile)
		#reader.next()
		for row in reader:
			tsv_list.append(row)
	return tsv_list


#parent_directory = os.path.dirname(os.path.dirname(os.getcwd()))
flat_file_dirctory = os.path.join('reconstruction', 'ecoli', 'flat')
rnas_file = os.path.join(flat_file_dirctory, 'rnas.tsv')
tus_file = os.path.join(flat_file_dirctory, 'operon_rnas.tsv')

#rna_info contains all the rna info directly pulled from rnas.tsv as a list of lists.
rna_info = parse_tsv(rnas_file)
tu_info = parse_tsv(tus_file)


#Number of RNA's
num_rnas = len(rna_info)
num_tus = len(tu_info)

#Create empty matrix of the intended size.
gene_tu_matrix = np.zeros((num_rnas, num_tus))

#Now go through the tu_info and find all the multi gene tus
#This is probably not a good way to do this.

multi_gene_tu = []
#looking for length of molecule_id instead of string parsing gene_id, bc
# i didnt feel like doing that for now. But def something to change later.

#actually go and do just gene_id only look if geneId contains an underscore
#if it does break up the gene_id to invidual ids to map. 

count = 0
tu_index = []
for row in tu_info:
	
	gene_id = row['geneId']
	split_delimiter = '_'
	if split_delimiter in gene_id:
		#genes_in_tu = re.split(split_delimiter, gene_id)
		tu_index.append(count)
		#multi_gene_tu.append(genes_in_tu)
	count += 1

# Now we have the locations of the TUs. Need to find individual genes in 
# rnas.tsv.
#LOTS of nested shit.


count = 0
for row in rna_info:
	# Get the gene_id of the RNAs
	gene_id = row['geneId']
	# Go through the tus
	gene_index = count
	for index in tu_index:
		genes_in_tu = re.split(split_delimiter, tu_info[tu_index[0]]['geneId'])

		# Try to find if the current RNA id matches an id with a TU.
		for gene in genes_in_tu:
			if gene == gene_id:
				gene_tu_matrix[gene_index, tu_index] = 1
				








import ipdb; ipdb.set_trace()
#with open 