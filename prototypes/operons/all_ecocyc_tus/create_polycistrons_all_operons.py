# Import packages
import itertools
import os
from operator import itemgetter
import csv
import numpy as np
import warnings

from functools import partial
from reconstruction import spreadsheets
from reconstruction.spreadsheets import tsv_reader

'''
NOTE:
This is only extracting polycistrons that form mRNAs.


TODO: right now im extracting waaaayyy more information than is needed. Will need
to condense it down to just what is needed.

Use tsvreader instead of json reader.

- output a readme file with relevant information : what genes are not encoded...etc
- write output with misordered genes
- Hand edit operons

When was the ecocyc file pulled?
- Get new TU's and coordinates file.

- Can we use a different synonyms file? Do we need the synonyms file?


'''

# --- Functions

def parseStringList(string):
	newStrArray = []
	splitStr = string.split('//')
	for st in splitStr:
		newStrArray.append(st.strip())
	return newStrArray

def gather_tu_information(tu_data, synonyms_data):
	'''
	Assemble data - record the object, for each gene and the operon direction.
	Note I was able to get this TU data a while ago, the ids and syns are new...
	Some 'genes' from before are now no longer genes so they dont show up in the new list
	so need to do a check to see that object ids are even found first before adding them to our list.	
	'''
	tu_information = []
	for row in tu_data:
		object_ids = parseStringList(row[1])
		if object_ids != []:
			tu_information.append([object_ids, row[4]])
	return tu_information

def rna_id_to_gene_id(gene_data_path):
	'''
	Returns a dictionary linking geneIds to RNA ids
	'''
	rna_to_gene_id_dict = {}
	with open(gene_data_path) as tsvfile:
		reader = JsonReader(tsvfile)
		fieldnames = reader.fieldnames
		for row in reader:
			rna_to_gene_id_dict[row['rna_id']] = row['id']
	return rna_to_gene_id_dict

def find_rna_ids(file_path, rna_to_gene_id_dict):
	mrnas_list = []
	all_rnas = []
	with open(file_path) as tsvfile:
		reader = JsonReader(tsvfile)
		fieldnames = reader.fieldnames
		for row in reader:
			try:
				#all_rnas.append(row['id']['geneId'])
				all_rnas.append(rna_to_gene_id_dict[row['id']])
			except:
				pass 
			if row['type'] == 'mRNA':
				mrnas_list.append(rna_to_gene_id_dict[row['id']])
	return all_rnas, mrnas_list

def returnNotMatches(a, b):
    return [[x for x in a if x not in b], [x for x in b if x not in a]]

def gather_gene_coordinates(gene_data_path):
	gene_coordinates = {}
	with open(gene_data_path) as tsvfile:
		reader = JsonReader(tsvfile)
		fieldnames = reader.fieldnames
		for row in reader:
			gene_coordinates[row['id']] = {}
			gene_coordinates[row['id']]['coordinate'] = row['coordinate']
			gene_coordinates[row['id']]['direction'] = row['direction']
	return gene_coordinates

def sort_strands(gene_coordinates):
	directions = ('+', '-')
	sorted_dictionary = {}
	sorted_strands = [[], []]
	reverse_bool = False
	for idx, direct in enumerate(directions):
		if direct == '-':
			reverse_bool = True
		#{key:value['coordinate'] for key, value in gene_coordinates.items() if value['direction'] == '+'}
		sorted_dictionary[direct] = {}
		sorted_dictionary[direct] = {key:value['coordinate'] for key, value in gene_coordinates.items() if value['direction'] == direct}
		sorted_strands[idx] = sorted(sorted_dictionary[direct], key=sorted_dictionary[direct].__getitem__, reverse=reverse_bool)
	return sorted_strands

def create_neighbor_dictionary(gene_coordinates):
	sorted_strands = sort_strands(gene_coordinates)
	neighbor_dictionary = {}
	for strand in sorted_strands:
		for idx, gene in enumerate(strand):
			previous_index = idx
			next_index = idx
			if next_index == len(strand)-1:
				next_index = -1
			neighbor_dictionary[gene] = {}
			neighbor_dictionary[gene]['previous_gene'] = strand[previous_index-1]
			neighbor_dictionary[gene]['next_gene'] = strand[next_index+1]
	return neighbor_dictionary

def check_gene_order(sorted_tu, neighbor_dictionary):
	'''
	The purpose of this function is to check that the nighboring genes to the target one match.
	'''
	tu_len = len(sorted_tu)
	mismatch_check = []

	def make_dict(tu_name, gene, direction, neighbor, ecocyc_tu_neighbor):
		check_dict = {}
		check_dict['Transcription Unit'] = tu_name
		check_dict['Current Gene'] = gene
		check_dict['Direction'] = direction
		check_dict['Real Neighbor'] = neighbor
		check_dict['Ecocyc TU Neighbor'] = ecocyc_tu_neighbor
		return check_dict

	def check_neighbor(idx, gene, sorted_tu, direction):
		if direction == 'previous_gene':
			dir_index = idx-1
		else:
			dir_index = idx+1
		check = []
		if sorted_tu[dir_index] != neighbor_dictionary[gene][direction]:
			check.append(make_dict('_'.join(sorted_tu), gene, direction, neighbor_dictionary[gene][direction], sorted_tu[dir_index]))
		return check

	directions = ['next_gene', 'previous_gene']
	for idx, gene in enumerate(sorted_tu):
		if idx == 0:
			mismatch_check.append(check_neighbor(idx, gene, sorted_tu, directions[0]))
		elif idx == tu_len - 1:
			mismatch_check.append(check_neighbor(idx, gene, sorted_tu, directions[1]))
		else:
			for direct in directions:
				mismatch_check.append(check_neighbor(idx,gene, sorted_tu, direct))
	mismatch_check = [x for x in mismatch_check if x]
	return mismatch_check

def sort_tu(row, gene_coordinates, mrna_ids, neighbor_dictionary):
	'''
	Returns:
	List of genes in sorted order.

	Sorts genes using a series of 'checks'
	'''
	# Prespecify a warning message:
	mismatch_warning_message = "For this operon ({}), the directions do not match"

	# initialize variables:
	row_dup_removed = set(row[0])
	tu_coordinates = []
	misordered_operon = []

	# Check 1:
	# Look to see if the gene has coordinates within our version of the genome. If not flag it.
	check_1 = 0
	for gene in row_dup_removed:
		try:
			tu_coordinates.append((gene, gene_coordinates[gene]['coordinate'], gene_coordinates[gene]['direction']))
		except KeyError:
			check_1 += 1

	# Check 2:
	# Check if the genes in the TU are mRNAs
	check_2 = 0
	for gene in row_dup_removed:
		if gene not in mrna_ids:
			check_2 +=1

	# Detect if the directions of the genes in the opeorn match
	if len(set([tu[2] for tu in tu_coordinates])) > 1:
		warnings.warn(mismatch_warning_message.format(tu[0]))

	# If any of the checks fail, then dont recored that TU
	if check_1 > 0 or check_2 > 0:
		sorted_tu = []
	# If the checks pass then sort the TU and record, reverse order if gene is on - strand
	elif check_1 == 0 and check_2 == 0:
		sorted_tu = [gene[0] for gene in sorted(tu_coordinates, key=itemgetter(1), reverse=False)]
		if row[1] == '-':
			sorted_tu.reverse()
		mismatch_check = check_gene_order(sorted_tu, neighbor_dictionary)
		if mismatch_check:
			misordered_operon.append(mismatch_check)
	return sorted_tu, misordered_operon

def find_polycistrons():
	return

def write_output_file(data, header, output_file):
	with open(output_file, "w") as f:
		writer = JsonWriter(f, header)
		writer.writeheader()
		for row in data:
			writer.writerow(row)
	return

def write_ouput_file_omit_genes(data, header, omission_list, output_file):
	with open(output_file, "w") as f:
		writer = JsonWriter(f, header)
		writer.writeheader()
		for row in data:
			count = 0
			for gene in row['transcription_units']:
				if gene in omission_list:
					count+=1
			if count == 0:
				writer.writerow(row)
	return

DIALECT = "excel-tab"
JsonReader = partial(spreadsheets.JsonReader, dialect = DIALECT)
JsonWriter = partial(spreadsheets.JsonWriter, dialect = DIALECT)


# --- Base file path to data in prototypes folder
data_loc = os.path.join(os.getcwd(), 'prototypes', 'operons', 'all_ecocyc_tus')

# --- Upload TU structure data
tu_path = os.path.join(data_loc, 'new_frameID_tus_coordinates_directions.txt')
tu_data = list(csv.reader(open(tu_path, 'r'), delimiter='\t'))[1:]

# --- Upload gene data
gene_data_path = os.path.join('reconstruction', 'ecoli', 'flat', 'genes.tsv')

# --- Upload synonyms data
synonyms_path = os.path.join(data_loc, '070320_gene_ids_syns.txt')
synonyms_data = list(csv.reader(open(synonyms_path, 'r'), delimiter = '\t'))

# --- Upload RNA data
rna_file_path = os.path.join('reconstruction', 'ecoli', 'flat', 'rnas.tsv')

# --- Specify output path
polycistrons_output_file = os.path.join(data_loc, 'all_polycistrons.tsv')
monocistrons_output_file = os.path.join(data_loc, 'monocistrons_to_keep.tsv')

# --- Get rid of duplicate rows in tu_data
tu_data.sort()
tu_data = list(tu_data for tu_data,_ in itertools.groupby(tu_data))

tu_information = gather_tu_information(tu_data, synonyms_data)

# how many unique genes are coded into operons
genes_in_operons = set([j for i in tu_information for j in i[0]])

rna_to_gene_id_dict = rna_id_to_gene_id(gene_data_path)
# find ids of all mRNAs
rna_ids, mrna_ids = find_rna_ids(rna_file_path, rna_to_gene_id_dict)

# find the genes that are not mrnas (doing the negative to save space)
#non_mrnas = [gene for gene in genes_in_operons if gene not in mrna_ids]

# Create polycistrons file, for each operon greater than one gene, delete all monomers.
polycistrons = [[row[0], row[1]] for row in tu_information if len(row[0]) > 1]
monocistrons = [row[0] for row in tu_information if len(row[0]) == 1]

gene_coordinates = gather_gene_coordinates(gene_data_path)
neighbor_dictionary = create_neighbor_dictionary(gene_coordinates)

polycistron_dictionary = find_polycistrons()


polycistrons_file = []
header = ['transcription_units', 'monomers_to_remove']
poly_dict_list = []
misordered_operons = []
for row in polycistrons:
	rowdict = {}
	sorted_tu, misordered_operon = sort_tu(row, gene_coordinates, mrna_ids, neighbor_dictionary)
	if sorted_tu != []:
		rowdict['transcription_units'] = sorted_tu
		rowdict['monomers_to_remove'] = [gene for gene in sorted_tu if [gene] not in monocistrons]
		#make sure not to repeat rows
		if rowdict not in poly_dict_list:
			poly_dict_list.append(rowdict)
	if misordered_operon:
		misordered_operons.append(misordered_operon)

#write_output_file(poly_dict_list, header, polycistrons_output_file)

functional_genes_file_path = os.path.join('prototypes/troubleshooting_mass_loss', 'functional_genes.tsv')
type_to_exclued = 'Metabolism'
def create_exclusion_file(functional_genes_file_path, type_to_exclued):
	'''
	with open(functional_genes_file_path, 'r') as j:
		file_reader = csv.reader(j, delimiter='\t')
		next(file_reader) # skip header
		next(file_reader)
		data = []
		for row in file_reader:
			data.append(row)
		breakpoint()
	'''
	exclusion_reader = csv.reader(open(functional_genes_file_path, 'r'), delimiter='\t')
	exclusion_list = list(exclusion_reader)[2:]
	exclusion_rnas = [row[1] for row in exclusion_list if row[3] == 'Transcription' or row[3] == 'RNA Decay']
	#exclusion_reader.close()
	return exclusion_rnas

exclusion_list = create_exclusion_file(functional_genes_file_path, type_to_exclued)
excluded = [rna[:-4] for rna in exclusion_list]

omission_output_file_path = os.path.join(data_loc, 'all_polycistrons_except_transcription_rna_decay_fg.tsv')
write_ouput_file_omit_genes(poly_dict_list, header, excluded, omission_output_file_path)

# use this portion in the case that you want to exclude certain RNAs from being in operons at all
'''
high_expressors_file_path = os.path.join('out', 'master_troubleshooting', 'top_expressers_master_rna_data.tsv')
high_expressers_list = list(csv.reader(open(high_expressors_file_path, 'r'), delimiter = '\t'))[0]
high_expressers = [rna[:-4] for rna in high_expressers_list]

omission_output_file_path = os.path.join(data_loc, 'polycistrons_no_master_rna_high_expressers.tsv')
write_ouput_file_omit_genes(poly_dict_list, header, high_expressers, omission_output_file_path)
'''



