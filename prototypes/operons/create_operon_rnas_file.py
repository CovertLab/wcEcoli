from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from collections import Counter
import csv
import itertools
import numpy as np
import os
import time

from functools import partial
from reconstruction import spreadsheets

'''
Purpose:
This file autogenerates the operon_rnas.tsv file. This file contains the 
the monocistronic and pcic mRNAs used within the whole cell model. Functionally
replaces rnas.tsv within the model.

How are the flat files used:
- polycistronic_mrnas_in_the_model:
	Supplies the pcic mRNAs to add as a list of geneIds, 
	in the column "transcription_units". For each transcription unit added 
	we need to specify if the monocistrons in that operon should be removed.
	This is because often the mono-cistrons in a transcription unit are not 
	individually transcribed. 
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
2. The location matches the first gene in the TU.
3. The miroarray value matches the first gene in the TU.

Returns:
	operon_rnas.tsv
		-This is a flat file that mimics the structure of rnas.tsv (with
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



Note:
- I have only functionally checked mRNAs; if adding tRNA or rRNA operons please 
	double check that all naming conventions hold correctly.
- Sequence included for polycistronic TUs includes the intergenic regions so will
	not sum to the length, width and nt counts of the component mono-cistronic
	mRNAs.
- A comment line is added to the top of the file to indicate that it generated
	from this script, and to record the date it was compiled.

To Do:

'''

WRITE_WITH_DICTWRITER = True

DIALECT = "excel-tab"
JsonReader = partial(spreadsheets.JsonReader, dialect = DIALECT)
JsonWriter = partial(spreadsheets.JsonWriter, dialect = DIALECT)

# File paths to all necessary flat files.
PROTOTYPES_DIR = os.path.join('prototypes', 'operons')
FLAT_DIR = os.path.join('reconstruction', 'ecoli', 'flat')
RNA_FILE = os.path.join(FLAT_DIR, "rnas.tsv")
GENES_FILE = os.path.join(FLAT_DIR, "genes.tsv")
POLY_CISTRON_FILE = os.path.join(PROTOTYPES_DIR, 'polycistronic_mrnas_in_model.tsv')
GENOME_SEQUENCE_FILE = os.path.join(FLAT_DIR, 'flattened_sequence.fasta')

#saving to a new file for now so that all manually input TUs in 
#operon_rnas.tsv are not overwritten.
OUTPUT_FILE = os.path.join(FLAT_DIR, "operon_rnas_testing.tsv")

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


	tu_types = [tu_genes_info[gene]['type'].lower() for gene in tu_genes_info]
	 
	if len(set(tu_types)) > 1:
		print(mismatch_output.format('_'.join(tu_genes_info.keys())))
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
		print(mismatch_output.format('_'.join(tu_genes_info.keys())))

	location = [max(loc_elements, key=Counter(loc_elements).get)]
	return location

def load_sequence(sequence):
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
	The rna sequence for the transcription unit.
	
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
	Input: RNA sequence
	Return: Mass of the RNA sequence
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
	Input: RNA sequence.
	Return: Counts of the nucleotidees in the sequence
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
		A nested dictionary of all the TUs to be added to the model 
		- first level - and information for all the genes within that
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
		tu_info[pc_gene_id]['comments'] = """Transcription unit created within script,
		"for individual RNA comments look at rnas.tsv for that RNA"""
		tu_info[pc_gene_id]['mw'] = calculate_rna_biomass(tu_info[pc_gene_id]['seq'])
		tu_info[pc_gene_id]['location'] = find_tu_location(tu_genes_info[pc_gene_id])
		tu_info[pc_gene_id]['ntCount'] = count_ntps_rna(tu_info[pc_gene_id]['seq'])
		if not tu_info[pc_gene_id]['type'] == 'mrna':
			tu_info[pc_gene_id]['id'] = pc_gene_id + '_' + tu_info[pc_gene_id]['type'].upper()
		else:
			tu_info[pc_gene_id]['id'] = pc_gene_id + '_RNA'
		tu_info[pc_gene_id]['geneId'] = pc_gene_id
		tu_info[pc_gene_id]['monomerSet'] = pc_monomer_id_list
		tu_info[pc_gene_id]['microarray expression'] = tu_genes_info[pc_gene_id][first_gene]['microarray expression']
	return tu_info


def find_monomers_to_remove():
	return set(itertools.chain(*[pc['monomers_to_remove'] 
		for pc in PC_INFO]))

def update_rna_info(fieldnames):
	#look for instances where rnas do not have an assigned monomerId
	#replace with a empty list, else put monomerId(s) into a list.
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
	Write either with dict writer or regular writer. The choice is for
	troubleshooting purposes only. 
	Dict write currently outputs in the correct format, but cannot
	accomidate a commented line, whereas the normal csv writer can. 
	'''
	if WRITE_WITH_DICTWRITER:
		# Create file with JSONWriter
		# Formatting is correct, but doesnt add top line.
		with open(OUTPUT_FILE, "w") as f:
			writer = JsonWriter(f, fieldnames)
			writer.writeheader()
			for rna_row in RNA_INFO:
				# only write monomers that we are keeping.
				if rna_row['geneId'] not in monomers_to_remove:
					writer.writerow(rna_row)
			for pc_data in tu_info:
				writer.writerow(tu_info[pc_data])
	else:
		# Create file with standard csv writer -
		# Appropriately adds header but all formatting is off.
		with open(OUTPUT_FILE, "w") as f:
			writer = csv.writer(f, delimiter='\t', quotechar = "'", 
				quoting = csv.QUOTE_MINIMAL, lineterminator="\n")
			writer.writerow(['# Generated by {} on {}'.format(__file__, time.ctime())])
			writer.writerow(RNA_INFO[0].keys())
			for rna_row in RNA_INFO:
				# only write monomers that we are keeping.
				if rna_row['geneId'] not in monomers_to_remove:
					writer.writerow(rna_row.values())
			for pc_data in tu_info:
				writer.writerow(tu_info[pc_data].values())
	return
	

def make_collection():
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

if __name__ == "__main__":
	make_collection()
