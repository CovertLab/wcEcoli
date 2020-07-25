# Import packages
import itertools
import os
import csv
import numpy as np

from functools import partial
from reconstruction import spreadsheets

'''
NOTE:
This is only extracting polycistrons that form mRNAs.


TODO: right now im extracting waaaayyy more information than is needed. Will need
to condense it down to just what is needed.

- output a readme file with relevant information : what genes are not encoded...etc

When was the ecocyc file pulled?
- Get new TU's and coordinates file.
'''

# --- Functions

def parseStringList(string):
	newStrArray = []
	splitStr = string.split('//')
	for st in splitStr:
		newStrArray.append(st.strip())
	return newStrArray

def get_ids(genes, synonyms_data):
	object_ids = []
	for gene in genes:
		for row in synonyms_data:
			if gene in row[1] or gene in row[2]:
				object_ids.append(row[0])
	return object_ids
def gather_tu_information(tu_data, synonyms_data):
	'''
	Assemble data - record the object, for each gene and the operon direction.
	Note I was able to get this TU data a while ago, the ids and syns are new...
	Some 'genes' from before are now no longer genes so they dont show up in the new list
	so need to do a check to see that object ids are even found first before adding them to our list.	
	'''
	tu_information = []
	for row in tu_data:
		object_ids = get_ids(parseStringList(row[0]), synonyms_data)
		if object_ids != []:
			tu_information.append([object_ids, row[1]])
	return tu_information

def find_rna_ids(file_path):
	mrnas_list = []
	all_rnas = []
	with open(file_path) as tsvfile:
		reader = JsonReader(tsvfile)
		fieldnames = reader.fieldnames
		for row in reader:
			try:
				all_rnas.append(row['geneId'])
			except:
				pass 
			if row['type'] == 'mRNA':
				mrnas_list.append(row['geneId'])
	return all_rnas, mrnas_list

def returnNotMatches(a, b):
    return [[x for x in a if x not in b], [x for x in b if x not in a]]


DIALECT = "excel-tab"
JsonReader = partial(spreadsheets.JsonReader, dialect = DIALECT)
JsonWriter = partial(spreadsheets.JsonWriter, dialect = DIALECT)


# --- Base file path to data in prototypes folder
data_loc = os.path.join(os.getcwd(), 'prototypes', 'operons', 'all_ecocyc_tus')

# --- Upload TU structure data
tu_path = os.path.join(data_loc, 'tus_coordinates_directions.txt')
tu_data = list(csv.reader(open(tu_path, 'r'), delimiter='\t'))[1:]

# --- Upload gene data

gene_data_path = os.path.join('reconstruction', 'ecoli', 'flat', 'genes.tsv')

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

def gather_gene_directions(gene_coordinates):
	postive_strand_gene_dict = {key:value['coordinate'] for key, value in gene_coordinates.items() if value['direction'] == '+'}
	negative_strand_gene_dict = {key:value['coordinate'] for key, value in gene_coordinates.items() if value['direction'] == '-'}

	positive_strand_sorted = sorted(postive_strand_gene_dict, key=postive_strand_gene_dict.__getitem__)
	negative_strand_sorted = sorted(negative_strand_gene_dict, key=negative_strand_gene_dict.__getitem__, reverse=True)

	import pdb; pdb.set_trace
	return positive_strand_sorted, negative_strand_sorted

def sort_tu(row, gene_coordinates, mrna_ids):
	row_dup_removed = list(dict.fromkeys(row[0]))
	tu_coordinates = []
	check_1 = 0
	for gene in row_dup_removed:
		try:
			tu_coordinates.append(gene_coordinates[gene])
		except KeyError:
			check_1 += 1
	check_2 = 0
	for gene in row_dup_removed:
		if gene not in mrna_ids:
			check_2 +=1
	if check_1 > 0 or check_2 > 0:
		sorted_tu = []
	elif check_1 == 0 and check_2 == 0:
		sorted_tu = [x for (y,x) in sorted(zip(np.argsort(tu_coordinates),row_dup_removed), key=lambda pair: pair[0])]
		if row[1] == '-':
			sorted_tu.reverse()
	return sorted_tu

def write_output_file(data, header, output_file):
	with open(output_file, "w") as f:
		writer = JsonWriter(f, header)
		writer.writeheader()
		for row in data:
			writer.writerow(row)
	return


# --- Upload synonyms data
synonyms_path = os.path.join(data_loc, '070320_gene_ids_syns.txt')
synonyms_data = list(csv.reader(open(synonyms_path, 'r'), delimiter = '\t'))

# --- Upload RNA data
rna_file_path = os.path.join('reconstruction', 'ecoli', 'flat', 'rnas.tsv')

# --- Specify output path
output_file = os.path.join(data_loc, 'all_polycistrons.tsv')

# --- Get rid of duplicate rows in tu_data
tu_data.sort()
tu_data = list(tu_data for tu_data,_ in itertools.groupby(tu_data))

tu_information = gather_tu_information(tu_data, synonyms_data)

# how many unique genes are coded into operons
genes_in_operons = set([j for i in tu_information for j in i[0]])

# find ids of all mRNAs
rna_ids, mrna_ids = find_rna_ids(rna_file_path)

# find the genes that are not mrnas (doing the negative to save space)
#non_mrnas = [gene for gene in genes_in_operons if gene not in mrna_ids]

# Create polycistrons file, for each operon greater than one gene, delete all monomers.
polycistrons = [[row[0], row[1]] for row in tu_information if len(row[0]) > 1]


# TODO: Remove repeat genes in operon.

# Reorder polycistrons:
'''
for pc in polycistrons:
	pc_coordinate = []
	for gene in pc[0]:
		pc_coordinate.append(gene_coordinates[gene])
	sorted_index = np.argsort(pc_coordinate)

# Some operons repeat but just in a different order. Keep only one order version.

for row in polycistrons:
	for check_row in polycistrons:
		if row[0] != check_row[0] and set(row[0]).difference(set(check_row[0])) == set([]) and set(check_row[0]).difference(set(row[0])) == set([]):
			polycistrons.remove(row[0])
'''
monomers_to_include = [row[0] for row in polycistrons if len(row[0]) == 1 and row[0] not in non_mrnas]
gene_coordinates = gather_gene_coordinates(gene_data_path)


polycistrons_file = []
header = ['transcription_units', 'monomers_to_remove']
poly_dict_list = []
for row in polycistrons:
	rowdict = {}
	sorted_tu = sort_tu(row, gene_coordinates, mrna_ids)
	if sorted_tu != []:
		rowdict['transcription_units'] = sorted_tu
		rowdict['monomers_to_remove'] = [gene for gene in sorted_tu if gene not in monomers_to_include]
		if rowdict not in poly_dict_list:
			poly_dict_list.append(rowdict)
# remove repeat rows in list.

# check for bi-directional operon genes

genes_ordered_pos, genes_ordered_neg = gather_gene_directions(gene_coordinates)

test_operon = genes_ordered_pos[0:4]

import pdb; pdb.set_trace()

write_output_file(poly_dict_list, header, output_file)


# - Collect some extra stats so we know more about the genes being included
# - Or excluded from the model.

pc_genes_excluded = set(list(itertools.chain.from_iterable(all_polycistrons_with_misc_mrnas)))
pc_genes_included = set(list(itertools.chain.from_iterable(find_all_polycistrons)))

def gather_tu_stats():
	return

# are any ids not present in our lists?
what_gene_is_this = [gene for gene in genes_in_operons if gene not in all_rnas]
# which mrnas are not not in a TU?
mRNAs_not_in_ecocyc_data = returnNotMatches(mrnas_list, genes_in_operons)

'''
find_all_polycistrons = []
all_polycistrons_with_misc_mrnas = []
all_genes_in_misc_operons = []
all_monomers_that_are_staying = []
for row in polycistrons:

		all_monomers_that_are_staying.append([gene for gene in row if gene in monomers_to_include])
		rowdict['monomers_to_remove'] = monomers_to_remove
		polycistrons_file.append(rowdict)
'''

# if __name__ == "__main__":

