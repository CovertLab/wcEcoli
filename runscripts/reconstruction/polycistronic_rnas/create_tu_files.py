import argparse
from collections import Counter
import csv
from functools import partial
import itertools
import numpy as np
import os
import re
from scipy.optimize import nnls
import time
import warnings

from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

from functools import partial
from reconstruction import spreadsheets

'''
Purpose:
This file autogenerates several files used to allow polycistronic transcription
units in the model (operon_rnas.tsv - functionally replaces rnas.tsv; gene_tu_matrix.tsv;
and transcription_units.tsv). 


Inputs:
- polycistronic_mrnas_in_the_model:
	Supplies the polycistronic mRNAs to add as a list of geneIds, 
	in the column "transcription_units". For each transcription unit added 
	we need to specify if the monocistrons in that operon should be removed.
	This is because often the mono-cistrons in a transcription unit are not 
	individually transcribed. Add these as a list in the column called 
	"monomers_to_remove".
	Note, as this file gets more complicated will need to make sure that these 
	mono-cistrons are consistent across all transcription units that they
	are a part of.
- rnas.tsv: 
	Supplies data about the individual mono-cistrons. operon_rnas.tsv mimics
	the layout/structure of this file.
- genes.tsv:
	Supplies data about the genes in the model. This file is used to gather
	the direction and coordinate of the genes in the pcic TUs for the
	purpose of finding the sequence of the pcic TUs including the 
	intergenic regions.
- flattened_sequence.fasta:
	A single line version of the sequence.fasta file so that its easier to 
	search and only has to be generated once outside the model instead of 
	within. 
Assumptions:
1. The half-life matches the first gene in the TU.
2. The location matches the first gene in the TU or the most prevalent one.
3. The miroarray value matches the first gene in the TU.

Returns:
	operon_rnas.tsv
		- Contains all the monomers and polycistronic mRNAs used in the model.
		- This is a flat file that mimics the structure of rnas.tsv (with
		an additional column 'monomerSet', that allows mapping of mRNAs
		back to their componenet monomer - this allows more than one
		monomer to map to a mRNA). 

		-Values are defined as follows:
		'halfLife' = half life of first gene in TU, taken from rnas.tsv
		'name' =  string of all names in TU separated by '_'
		'seq' = sequence of all genes in TU including all intergenic regions
		'type' = str indiciated the type of RNA
		'modifiedForms' = list containing modified form each RNA within
			a TU can form.
		'monomerId' = string of all monomerIds in TU separated by '_'
		'comments' = "Transcription unit created within script, for 
			individual RNA comments look at rnas.tsv for that RNA"
		'mw' = mw of the TU sequence
		'location' = location of the first gene in a 
		'ntCount' = count_ntps_rna(tu_info[pc_gene_id]['seq'])
		'id' = pc_gene_id + '_RNA'
		'geneId' = pc_gene_id
		'monomerSet' = pc_monomer_id_list
		'microarray expression' = tu_genes_info[pc_gene_id][rna]['microarray expression']
	gene_to_tu_matrix.tsv
		- A boolean matrix, rows = genes, cols = tus, which serves as a maping of
		genes to tus.
	transcription_units.tsv
		- Counts of each transcription unit included in the model. Functionally replaces
		rna seq counts in the model, at least for the basal case. See todos for future
		upgrades to this file
Note:
- I have only functionally checked mRNAs; if adding tRNA or rRNA operons please 
	double check that all naming conventions hold correctly - updates
	may need to be made throughout the model.
- Sequence included for polycistronic TUs includes the intergenic regions so will
	not sum to the length, width and nt counts of the component mono-cistronic
	mRNAs.
- A comment line is added to the top of the file to indicate that it generated
	from this script, and to record the date it was compiled.
TODO:
- UPDATE PROTEINS.TSV
- Take polycistrons file as an input argument
- Allow for certain RNA type mixing (rRNA, tRNA) - handle this type differently
- Allow mass to determine type and location for importing mass fractions in the Parca
- Right now I am just making the transcription_units file for a single condition. In the 
future we will want it for all conditions. Will need to take care in the parca that the 
correct condiditon data is being used.
- Dont allow for deletion of a gene without incorporation in a TU.

-Import helper functions from another file.
- Fix parse tsv2.
- Allow protiens.tsv to overwrite itself, check if rnaset column already exists
- Make parameters consistent across functions.
- check that directions of genes in a TU match
-fix furure warning

'''
def parse_args():
	'''
	return error message if argument is not given
	'''
	parser = argparse.ArgumentParser()
	arser = argparse.ArgumentParser()
	parser.add_argument('filename')
	args = parser.parse_args()
	return args.filename

DIALECT = "excel-tab"
JsonReader = partial(spreadsheets.JsonReader, dialect = DIALECT)
JsonWriter = partial(spreadsheets.JsonWriter, dialect = DIALECT)

# File paths to all necessary flat files.
FLAT_DIR = os.path.join('reconstruction', 'ecoli', 'flat')
RNA_FILE = os.path.join(FLAT_DIR, 'rnas.tsv')
GENES_FILE = os.path.join(FLAT_DIR, 'genes.tsv')
#POLY_CISTRON_FILE = os.path.join(FLAT_DIR, 'polycistronic_mrnas_in_model.tsv')
POLY_CISTRON_FILE = parse_args()
GENOME_SEQUENCE_FILE = os.path.join(FLAT_DIR, 'flattened_sequence.fasta')
km_file = os.path.join('fixtures', 'endo_km', 'km.cPickle')
RNA_SEQ_FILE = os.path.join(FLAT_DIR, 'rna_seq_data', 'rnaseq_rsem_tpm_mean.tsv')
PROTEIN_FILE = os.path.join(FLAT_DIR, 'proteins_old.tsv')

# output files
TU_FILE = os.path.join(FLAT_DIR, 'operon_rnas.tsv')
output_tu_counts = os.path.join(FLAT_DIR, "transcription_units.tsv")
output_gene_tu_matrix = os.path.join(FLAT_DIR, "gene_to_tu_matrix.tsv")
output_proteins = os.path.join(FLAT_DIR, "proteins.tsv")

CONDITION = 'M9 Glucose minus AAs'
SPLIT_DELIMITER = '_'

'''
parser = argparse.ArgumentParser()
parser.add_argument('filename')
args = parser.parse_args()
with open(args.filename) as file:
	print(file.readline())
import pdb; pdb.set_trace()
#with open(args.filename) as file:
'''


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

def read_file_skipping_comments(tsvfile, tsv_list, tempfile):
	num_lines = 0
	for curline in tsvfile:
		if not curline.startswith('#'):
			tempfile.write(curline)
	tempfile.close() #need to close file for IOS memory reasons.
	with open('temp.csv') as tempfile:
		reader = JsonReader(tempfile)
		fieldnames = reader.fieldnames
		for row in reader:
			tsv_list.append(row)
	if os.path.exists("temp.csv"):
		os.remove("temp.csv")
	else:
		print("temp.csv does not exist")
	return tsv_list, fieldnames

def parse_tsv_2(tsv_file):
	'''
	Takes in a tsv file, and creates a list of dicts of the rows 
	contained within the TSV.
	'''
	tsv_list = []
	with open(tsv_file, 'r') as tsvfile, open('temp.csv', 'w') as tempfile:
		# read in first line in tsvfile
		if not tsvfile.readline().startswith('#'):
			tsv_list=[]
			reader = JsonReader(tsvfile)
			fieldnames = reader.fieldnames
			for row in reader:
				tsv_list.append(row)
		else:
			tsv_list, fieldnames = read_file_skipping_comments(tsvfile, tsv_list, tempfile)
	return tsv_list, fieldnames

def find_tu_type(tu_genes_info):
	'''
	Purpose:
	Need to assign an rna type to the pcic RNA (e.g. rRNA,
	mRNA...)
	
	If all the RNA types in a TU are the same, the type will be assigned
	to that type.

	If they are not the same, the most prevalent one will be assigned
	for the whole TU, or if all types are equally represented, the 
	type from the first RNA in the TU will be assigned.

	Inputs:
	Dictionary of transcription unit info for target transcription unit.

	Returns:
	RNA type for the transcription unit.
	'''
	mismatch_output = """Types of RNA's in transcription unit {} dont 
	match, please double check that your transcription unit is correct. 
	Might have to do additional checks. The naming for this transcription 
	unit will follow the most prevalant RNA type, or the first RNA in the TU 
	if there is an even number of RNAs."""
	#tu_types = [tu_genes_info[gene]['type'].lower() for gene in tu_genes_info]
	tu_types = [tu_genes_info[gene]['type'] for gene in tu_genes_info]
	if len(set(tu_types)) > 1:
		warnings.warn(mismatch_output.format('_'.join(tu_genes_info.keys())))
	tu_type = max(tu_types, key=Counter(tu_types).get)
	return tu_type

def find_tu_location(tu_genes_info):
	'''
	Mimics structure and output of find_tu_type but comparing locations
	of the TU instead.

	Return: 
		List containing a unicode str of the location.
	'''
	mismatch_output = """Locations of RNA's in transcription unit {} dont 
	match, please double check that your transcription unit is correct. 
	Might have to do additional checks. The naming for this transcription 
	unit will follow the most prevalant RNA location, or the first RNA in the TU 
	if there is an even number of RNAs."""
	tu_locations = [tu_genes_info[gene]['location'] for gene in tu_genes_info]
	loc_elements = [loc[0] for loc in tu_locations]
	if len(set(loc_elements)) > 1:
		warnings.warn(mismatch_output.format('_'.join(tu_genes_info.keys())))
	location = [max(loc_elements, key=Counter(loc_elements).get)]
	return location

def load_sequence(sequence):
	'''
	Input:
	- Path to genomic sequence fasta file, who is formatted as a single line.
	Returns:
	- List containing a string of the genomic sequence.
	'''
	with open(sequence, "r") as f:
		genome_sequence = f.readlines()[1:]
	return genome_sequence

def get_tu_sequence(tu_genes_info, first_gene, last_gene):
	'''
	Purpose:
	Obtain the RNA sequence for the entire transcription unit.
	Cannot simply concatenate component individual RNA sequences since
	they would not include intergenic regions.
	
	Inputs:
	- Genomic sequence. 
	- Dictionary containing information for a transcription unit (Dict 
	contains direction and coordinate info).
	- First and last gene in the transcription unit.

	Return:
	-The rna sequence for the transcription unit.
	
	Note:
	This was developed with only implementing mRNAs for now. If 
	incorporating untranslated RNAs (like rRNAs or tRNAs) in the future
	would have to add a process for RNA processing.
	'''
	direction = tu_genes_info[first_gene]['direction'].encode('utf-8')
	rna_sequence = []
	if direction is '+':
		first_gene_coordinate = tu_genes_info[first_gene]['chromosome_coordinate']
		last_gene_coordinate = tu_genes_info[last_gene]['chromosome_coordinate'] + tu_genes_info[last_gene]['length']
		sequence = Seq(GENOMIC_SEQUENCE[0][first_gene_coordinate:last_gene_coordinate], IUPAC.unambiguous_dna)
		rna_sequence = str(sequence.transcribe())
	elif direction is '-':
		first_gene_coordinate = tu_genes_info[first_gene]['chromosome_coordinate'] +1
		last_gene_coordinate = tu_genes_info[last_gene]['chromosome_coordinate'] - tu_genes_info[last_gene]['length'] + 1
		sequence = Seq(GENOMIC_SEQUENCE[0][last_gene_coordinate:first_gene_coordinate], IUPAC.unambiguous_dna)
		rna_sequence = str(sequence.reverse_complement().transcribe())
	else:
		print("For some reason your gene hasnt been assigned a direction")
		
	return rna_sequence

def calculate_rna_biomass(sequence):
	'''
	Purpose: Calculate the RNA biomass of a particular RNA sequence.
	Input: 
	- RNA sequence as a str.
	Return: 
	- Mass of the RNA sequence as a float.
	'''
	rna_masses = {
		"A": 503.15,
		"C": 479.124,
		"G": 519.149,
		"U": 480.108,
	}

	ppi_mass = 174.949

	base_order = ['A', 'C', 'G', 'U']
	ntp_order = {base + 'TP': index for index, base in enumerate(base_order)}
	counts = {base: (rna_masses[base] - ppi_mass) * sequence.count(base) for base in base_order}
	
	return [0.0, 0.0, 0.0, np.sum(counts.values()) + ppi_mass, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
	
def count_ntps_rna(sequence):
	'''
	Input: 
	- RNA sequence.
	Return: 
	- Counts of the nucleotidees in the sequence
	'''

	return [sequence.count('A'), sequence.count('C'),
			sequence.count('G'), sequence.count('U')]

def gather_tu_genes_info():
	'''
	Purpose: Gather data for each rna within each TU.

	Input:
	- PC_INFO: A list of dictionaries, where each row is a new 
	polycistronic TU that needs to be added to the model. 
		- 'transcription_units': list of genes in the TU.
		- 'monomers_to_remove': list of monocistronic genes from the
			TU to remove from the model.
	Returns:
	-A nested dictionary of all the TUs to be added to the model. 
	First level - and information for all the genes within that
	TU. This data will used in a later function to define 
	all the features of a specific TU. 
	'''
	tu_genes_info = {}
	for pc in PC_INFO:
		pc_gene_id = '_'.join(pc['transcription_units'])
		tu_genes_info[pc_gene_id] = {}
		for rna in pc['transcription_units']:
			tu_genes_info[pc_gene_id][rna] = {}
			for rna_row in RNA_INFO:
				if rna_row['geneId'] == rna:
					tu_genes_info[pc_gene_id][rna]['halfLife'] = rna_row['halfLife']
					tu_genes_info[pc_gene_id][rna]['microarray expression'] = rna_row['microarray expression']
					tu_genes_info[pc_gene_id][rna]['name'] = rna_row['name'][0:-6]
					tu_genes_info[pc_gene_id][rna]['type'] = rna_row['type']
					tu_genes_info[pc_gene_id][rna]['modifiedForms'] = rna_row['modifiedForms']
					tu_genes_info[pc_gene_id][rna]['monomerId'] = rna_row['monomerId']
					tu_genes_info[pc_gene_id][rna]['location'] = rna_row['location']
			for gene_row in GENE_INFO:
				if gene_row['id'] == rna:
					tu_genes_info[pc_gene_id][rna]['chromosome_coordinate'] = gene_row['coordinate']
					tu_genes_info[pc_gene_id][rna]['direction'] = gene_row['direction']
					tu_genes_info[pc_gene_id][rna]['length'] = gene_row['length']
	return tu_genes_info

def gather_tu_info(tu_genes_info):
	'''
	Purpose:
	Need to gather all the info for each TU we need to add it to the 
	model. Format closely follows rnas.tsv with modifications outlined
	at the top of this file.
	Inputs
	- PC_INFO: A list of dictionaries, where each row is a new 
	polycistronic TU that needs to be added to the model. 
		- 'transcription_units': list of genes in the TU.
		- 'monomers_to_remove': list of monocistronic genes from the
			TU to remove from the model.
	- tu_genes_info:
		Dictionary created within gather_tu_genes_info, nested dictionary
		containing info for each gene within each TU.
	Returns:
	- A nested dictionary. First level keys are the gene ids for 
	each polycistron being added to the model, nested keys are all the 
	columns that will be needed within the operon_rnas.tsv file. 
	For more detail on each key, please see the comment at the top
	of this file. Naming mimics rnas.tsv.
	'''
	tu_info = {}
	for pc in PC_INFO:
		# Get the gene ID of the polycistron
		pc_gene_id = '_'.join(pc['transcription_units'])

		# Find the monomerIds, names and modifiedForms of all the RNAs in the TU
		pc_monomer_id_list = []
		pc_name_list = []
		pc_modified_forms_list = []
		for rna in pc['transcription_units']:
			pc_monomer_id_list.append(tu_genes_info[pc_gene_id][rna]['monomerId'])
			pc_name_list.append(tu_genes_info[pc_gene_id][rna]['name'])
			pc_modified_forms_list.append(tu_genes_info[pc_gene_id][rna]['modifiedForms'])
		
		# Find the first and last genes in the TU
		first_gene = pc['transcription_units'][0]
		last_gene = pc['transcription_units'][-1]

		# Construct dictionary of transcription unit info for each tu.
		tu_info[pc_gene_id] = {}
		tu_info[pc_gene_id]['halfLife'] = tu_genes_info[pc_gene_id][first_gene]['halfLife']
		tu_info[pc_gene_id]['name'] =  ':'.join(pc_name_list)
		tu_info[pc_gene_id]['seq'] = get_tu_sequence(tu_genes_info[pc_gene_id], 
			first_gene, last_gene)
		tu_info[pc_gene_id]['type'] = find_tu_type(tu_genes_info[pc_gene_id])
		tu_info[pc_gene_id]['modifiedForms'] = []
		tu_info[pc_gene_id]['monomerId'] = '_'.join(pc_monomer_id_list)
		tu_info[pc_gene_id]['comments'] = """Transcription unit created within script, for individual RNA comments look at rnas.tsv for that RNA"""
		tu_info[pc_gene_id]['mw'] = calculate_rna_biomass(tu_info[pc_gene_id]['seq'])
		tu_info[pc_gene_id]['location'] = find_tu_location(tu_genes_info[pc_gene_id])
		tu_info[pc_gene_id]['ntCount'] = count_ntps_rna(tu_info[pc_gene_id]['seq'])
		if not tu_info[pc_gene_id]['type'] == 'mRNA':
			tu_info[pc_gene_id]['id'] = pc_gene_id + '_' + tu_info[pc_gene_id]['type'].upper()
		else:
			tu_info[pc_gene_id]['id'] = pc_gene_id + '_RNA'
		tu_info[pc_gene_id]['geneId'] = pc_gene_id
		tu_info[pc_gene_id]['monomerSet'] = pc_monomer_id_list
		tu_info[pc_gene_id]['microarray expression'] = tu_genes_info[pc_gene_id][first_gene]['microarray expression']
	return tu_info


def find_monomers_to_remove():
	'''
	Purpose: Find all the momnomers that should not be included in
	the final operon_rnas.tsv file.
	Input:
	- Global variable, of the list of user input data on polycistrons
	to include and monomers to remove.
	Return:
	- Set of unique monomers to not include in final tsv.
	'''
	return set(itertools.chain(*[pc['monomers_to_remove'] 
		for pc in PC_INFO]))

def update_rna_info(fieldnames):
	'''
	Purpose:
	Add an additional key for each RNA within RNA_INFO specificing the 
	monomers assigned to that mRNA. 
	Input:
	- RNA_INFO
	Return:
	- List of Dictionaries for all the RNA_INFO but this time including
	a key for monomerSet which contains a list of all the monomers assigned
	to a particular mRNA, list will be empty if there are no monomers to
	assign to an RNA.
	'''
	for rna_row in RNA_INFO:
		if not rna_row['monomerId']:
			rna_row['monomerSet'] = []
		else:
			rna_row['monomerSet'] = [rna_row['monomerId']]
	# add fieldname for 'monomersets'
	fieldnames.append('monomerSet')
	return RNA_INFO, fieldnames

def write_output_file(tu_info, fieldnames, monomers_to_remove):
	'''
	Construct a tsv file that mimics the structure and formatting
	of rnas.tsv.
	'''
	# Create file with JSONWriter
	with open(TU_FILE, "w") as f:
		# Add comment line with file creation info
		f.write('# Generated by {} on {} \n'.format(__file__, time.ctime()))

		# Write file with JsonWriter
		writer = JsonWriter(f, fieldnames)
		writer.writeheader()
		for rna_row in RNA_INFO:
			# only write monomers that we are keeping.
			if rna_row['geneId'] not in monomers_to_remove:
				writer.writerow(rna_row)
		for pc_data in tu_info:
			writer.writerow(tu_info[pc_data])
	return

def make_operon_rnas_file():
	global PC_INFO, RNA_INFO, GENE_INFO, GENOMIC_SEQUENCE
	RNA_INFO, fieldnames = parse_tsv(RNA_FILE)
	PC_INFO, pc_fieldnames = parse_tsv(POLY_CISTRON_FILE)
	GENE_INFO, g_fieldnames = parse_tsv(GENES_FILE)
	GENOMIC_SEQUENCE = load_sequence(GENOME_SEQUENCE_FILE)

	tu_genes_info = gather_tu_genes_info()
	tu_info = gather_tu_info(tu_genes_info)

	monomers_to_remove = find_monomers_to_remove()
	RNA_INFO, fieldnames = update_rna_info(fieldnames)

	write_output_file(tu_info, fieldnames, monomers_to_remove)


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
	#import pdb; pdb.set_trace()
	rna_seq_data_index = {
		row['Gene']: row[CONDITION] 
		for row in rna_seq_data_all_cond[0]}

	rna_seq_counts_vector = [
		rna_seq_data_index[gene] 
		for gene in rnas_gene_order]

	return rna_seq_counts_vector

def create_tu_counts_vector(gene_tu_matrix, rna_seq_counts_vector, tu_info):
	'''
	Purpose:
		Calculate transcription unit counts. Will functionally replace how rna seq
		counts were used in the model. If no TUs are added it will exactly match
		rna-seq counts. Using scipy.optimize..nnls solver to give a non-negative least squares fit.
	Inputs:
		gene_tu_matrix: Sparse numpy matrix, mapping TU to rna's. 0 = no mapping; 1 = mapping.
		rna_seq_counts_vector: A vector of RNA seq counts for the condition we are interested in.
		tu_info: list of dictionaries containing all the information we have for all the
		transcription units we want to include in the model.
	Outpus:
		A list of dictionaries. Each dictionary contains the 'tu_id' and the 'tu_count' for
		each TU in the model.
	'''
	tu_counts_vector = nnls(gene_tu_matrix, rna_seq_counts_vector)[0]
	tu_gene_order = [row['geneId'] for row in tu_info]
	tu_genes_counts = []

	for i in range(0, len(tu_gene_order)):
		tu_genes_counts.append({'tu_id': tu_gene_order[i], 'tu_count': tu_counts_vector[i]})

	return tu_genes_counts


def make_transcription_units_file():
	'''
	Purpose: 
		Run all the functions necessary to get the transcription counts vector 
		output as transcription_units.tsv.
	'''
	tu_info, fieldnames_rna = parse_tsv_2(TU_FILE)
	gene_tu_matrix, rnas_gene_order = create_gene_to_tu_matrix(RNA_INFO, tu_info)
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

def make_new_proteins_file(output_file):
	'''
	TODO:
	Check if the key monomer set exists.
	'''

	#import pdb; pdb.set_trace()
	protein_info, protein_fieldnames = parse_tsv(PROTEIN_FILE)
	rna_info, rna_fieldnames = parse_tsv_2(TU_FILE)

	#Go through monomerSet line by line. Find the matching monomers within
	#those lists then find the corresponding monomer in proteins.tsv.
	#Add the id from operon_rnas to the rnaSet list

	protein_index = {}
	for protein_row in protein_info:
		protein_row['rnaSet'] = []
		protein_index[protein_row['id']] = protein_row

	
	for rna_row in rna_info:
		for monomer in rna_row['monomerSet']:
			protein_row = protein_index[monomer]
			protein_row['rnaSet'].append(rna_row['id'])	

	
	#add fieldname for 'monomersets'
	protein_fieldnames.append('rnaSet')

	with open(output_file, "w") as f:
		writer = JsonWriter(f, protein_fieldnames)
		writer.writeheader()
		for protein_row in protein_info:
			writer.writerow(protein_row)


def remove_kms_file(km_file):
	'''
	Purpose:
	The parca checks if the km's have already been calculated, if they have then
	it does not calculate them again and just pulls them from the km.cpickle file.
	However, if the operon_rnas file gets modified the kms need to be recalculated.
	This function will remove the .cpickle file and force the parca to re-calculate the
	Kms, this way we can be sure we are working with a correct dataset.
	'''

	if os.path.exists(km_file):
		os.remove(km_file)
	return

if __name__ == "__main__":
	make_operon_rnas_file()
	make_transcription_units_file()
	remove_kms_file(km_file)
	make_new_proteins_file(output_proteins)
