'''
The goal of this script is to take in a file containing all the possible TU's (this is most like what we would expect from the sequencing data) and create a gene to TU matrix. 
For now extract the TU list from the rna_operons.tsv file. Later this will be provided. 
Extract the gene list from rnas.tsv.
'''

import csv
import os

def parse_tsv(tsv_file):
	'''
	Takes in a tsv file, and creates a list of lists of the rows 
	contained within the TSV.
	'''
	tsv_list = []
	with open(tsv_file) as tsvfile:
		reader = csv.reader(tsvfile, delimiter = '\t')
		reader.next()
		for row in reader:
			tsv_list.append(row)
	return tsv_list
#this gives directory into wcEcoli. Particular to the current storage location 
# of this folder within prototypes.
parent_directory = os.path.dirname(os.path.dirname(os.getcwd()))
flat_file_dirctory = os.path.join(parent_directory, 'reconstruction', 'ecoli', 'flat')
rnas_file = os.path.join(flat_file_dirctory, 'rnas.tsv')
tus_file = os.path.join(flat_file_dirctory, 'operon_rnas.tsv')

#rna_info contains all the rna info directly pulled from rnas.tsv as a list of lists.
rna_info = parse_tsv(rnas_file)
tu_info = parse_tsv(tus_file)


#Number of RNA's
num_rnas = len(rna_info)
num_tus = len(tu_info)




import ipdb; ipdb.set_trace()
#with open 