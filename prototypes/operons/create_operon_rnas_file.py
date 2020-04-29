from collections import Counter
import itertools
import os
import csv
import time
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
import numpy as np
from functools import partial
from reconstruction import spreadsheets


'''
Purpose:
This file autogenerates the operon_rnas.tsv file. This file contains the 
the monocistronic and polycistronic mRNAs used within the whole cell model.

How are the flat files used:
- polycistronic_mrnas_in_the_model:
	Supplies the polycistronic mRNAs to add as a list of geneIds, 
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
	the direction and coordinate of the genes in the polycistronic TUs for the
	purpose of finding the sequence of the polycistronic TUs including the 
	intergenic regions.
- flattened_sequence.fasta:
	A single line version of the sequence.fasta file so that its easier to 
	search and only has to be generated once outside the model instead of 
	within. 
Assumptions:
1. The half-life matches the first gene in the TU.
2. The location matches the first gene in the TU.
3. The miroarray value matches the first gene in the TU.

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

DIALECT = "excel-tab"
JsonReader = partial(spreadsheets.JsonReader, dialect = DIALECT)
JsonWriter = partial(spreadsheets.JsonWriter, dialect = DIALECT)


PROTOTYPES_DIR = os.path.join('prototypes', 'operons')
FLAT_DIR = os.path.join('reconstruction', 'ecoli', 'flat')
RNA_FILE = os.path.join(FLAT_DIR, "rnas.tsv")
GENES_FILE = os.path.join(FLAT_DIR, "genes.tsv")
POLY_CISTRON_FILE = os.path.join(PROTOTYPES_DIR, 'polycistronic_mrnas_in_model.tsv')
#saving to a new file for now so that all manually input TUs in 
#operon_rnas.tsv are not overwritten.
output_file = os.path.join(FLAT_DIR, "operon_rnas_testing.tsv")
GENOME_SEQUENCE_FILE = os.path.join(FLAT_DIR, 'flattened_sequence.fasta')

'''

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

def find_tu_type(transcription_unit_info):

	'''
	Purpose:
	Need to assign an rna type to the polycistronic RNA (e.g. rRNA,
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

	tu_types = []
	for gene in transcription_unit_info:
		tu_types.append(transcription_unit_info[gene]['type'].lower())
	if len(set(tu_types)) > 1:
		print(mismatch_output.format(tu))
	tu_type = max(tu_types, key=Counter(tu_types).get)	
	return tu_type

def load_sequence(sequence):
	with open(sequence, "r") as f:
		genome_sequence = f.readlines()[1:]
	return genome_sequence

def get_tu_sequence(gene_sequence, tu_info, first_gene, last_gene):
	'''
	Purpose:

	Inputs:


	'''
	direction = tu_info[first_gene]['direction'].encode('utf-8')
	rna_sequence = []
	if direction is '+':
		first_gene_coordinate = tu_info[first_gene]['chromosome_coordinate']
		last_gene_coordinate = tu_info[last_gene]['chromosome_coordinate'] + tu_info[last_gene]['length']
		sequence = Seq(gene_sequence[0][first_gene_coordinate:last_gene_coordinate], IUPAC.unambiguous_dna)
		rna_sequence = str(sequence.transcribe())
		#convert to rna
	elif direction is '-':
		first_gene_coordinate = tu_info[first_gene]['chromosome_coordinate'] +1
		last_gene_coordinate = tu_info[last_gene]['chromosome_coordinate'] - tu_info[last_gene]['length'] + 1
		sequence = Seq(gene_sequence[0][last_gene_coordinate:first_gene_coordinate], IUPAC.unambiguous_dna)
		rna_sequence = str(sequence.reverse_complement().transcribe())
	else:
		print("For some reason your gene hasnt been assigned a direction")
		
	return rna_sequence

def calculate_rna_biomass(sequence):
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
	
def counts_ntps(sequence):

	return [sequence.count('A'), sequence.count('C'),
			sequence.count('G'), sequence.count('U')]

def make_collection():
	#pull in rna_info as a list
	rna_info, fieldnames = parse_tsv(RNA_FILE)
	polycistron_info, pc_fieldnames = parse_tsv(POLY_CISTRON_FILE)
	genes_info, g_fieldnames = parse_tsv(GENES_FILE)
	gene_sequence = load_sequence(GENOME_SEQUENCE_FILE)




	# Gather data for each rna within each TU.
	transcription_unit_info = {}
	for polycistron in polycistron_info:
		polycistron_gene_id = '_'.join(polycistron['transcription_units'])
		transcription_unit_info[polycistron_gene_id] = {}
		for rna in polycistron['transcription_units']:
			transcription_unit_info[polycistron_gene_id][rna] = {}
			for rna_row in rna_info:
				if rna_row['geneId'] == rna:
					transcription_unit_info[polycistron_gene_id][rna]['halfLife'] = rna_row['halfLife']
					transcription_unit_info[polycistron_gene_id][rna]['microarray expression'] = rna_row['microarray expression']
					transcription_unit_info[polycistron_gene_id][rna]['name'] = rna_row['name'][0:-6]
					transcription_unit_info[polycistron_gene_id][rna]['type'] = rna_row['type']
					transcription_unit_info[polycistron_gene_id][rna]['modifiedForms'] = rna_row['modifiedForms']
					transcription_unit_info[polycistron_gene_id][rna]['monomerId'] = rna_row['monomerId']
					transcription_unit_info[polycistron_gene_id][rna]['location'] = rna_row['location']
			for gene_row in genes_info:
				if gene_row['id'] == rna:
					transcription_unit_info[polycistron_gene_id][rna]['chromosome_coordinate'] = gene_row['coordinate']
					transcription_unit_info[polycistron_gene_id][rna]['direction'] = gene_row['direction']
					transcription_unit_info[polycistron_gene_id][rna]['length'] = gene_row['length']
	
	# Check that RNA types match:

	#check_rna_types_match(transcription_unit_info)
	# Go through the indivuidual gene data to construct the multigene tu data
	combined_tu_info = {}
	for polycistron in polycistron_info:
		polycistron_gene_id = '_'.join(polycistron['transcription_units'])
		first_gene = polycistron['transcription_units'][0]
		last_gene = polycistron['transcription_units'][-1]
		polycistron_monomer_id_list = []
		polycistron_name_list = []
		combined_tu_info[polycistron_gene_id] = {}
		for rna in polycistron['transcription_units']:
			polycistron_monomer_id_list.append(transcription_unit_info[polycistron_gene_id][rna]['monomerId'])
			polycistron_name_list.append(transcription_unit_info[polycistron_gene_id][rna]['name'])
			#polycistron_id_list.append(transcription_unit_info[rna]['id'])
		first_gene_id = polycistron['transcription_units'][0]
		combined_tu_info[polycistron_gene_id]['halfLife'] = transcription_unit_info[polycistron_gene_id][rna]['halfLife']
		combined_tu_info[polycistron_gene_id]['name'] =  ':'.join(polycistron_name_list)
		
		combined_tu_info[polycistron_gene_id]['seq'] = get_tu_sequence(gene_sequence, transcription_unit_info[polycistron_gene_id], 
			first_gene, last_gene)
		
		combined_tu_info[polycistron_gene_id]['type'] = find_tu_type(transcription_unit_info[polycistron_gene_id])
		combined_tu_info[polycistron_gene_id]['modifiedForms'] = []
		combined_tu_info[polycistron_gene_id]['monomerId'] = '_'.join(polycistron_monomer_id_list)
		combined_tu_info[polycistron_gene_id]['comments'] = "Transcription unit created within script, for individual RNA comments look at rnas.tsv for that RNA"
		combined_tu_info[polycistron_gene_id]['mw'] = calculate_rna_biomass(combined_tu_info[polycistron_gene_id]['seq'])
		combined_tu_info[polycistron_gene_id]['location'] = transcription_unit_info[polycistron_gene_id][rna]['location']
		combined_tu_info[polycistron_gene_id]['ntCount'] = counts_ntps(combined_tu_info[polycistron_gene_id]['seq'])
		combined_tu_info[polycistron_gene_id]['id'] = polycistron_gene_id + '_RNA'
		combined_tu_info[polycistron_gene_id]['geneId'] = polycistron_gene_id
		combined_tu_info[polycistron_gene_id]['monomerSet'] = polycistron_monomer_id_list
		combined_tu_info[polycistron_gene_id]['microarray expression'] = transcription_unit_info[polycistron_gene_id][rna]['microarray expression']
	

	#Creating a new column for monomerSets. Which are the monomers for given
	#mRNA transcripts. Right now we are manually adding in multi-gene TUs, 
	#so Im not going to worry about anything except taking the current rnas.tsv
	#and creating an output that has a new column putting the monomerID into a list.
	#for a multigene TU, will have two IDs within this list.

	#look for instances where rnas do not have an assigned monomerId
	#replace with a empty list, else put monomerId(s) into a list.
	monomers_to_remove = []
	for pc_info in polycistron_info:
		monomers_to_remove.append(pc_info['monomers_to_remove'])
	monomers_to_remove = set(itertools.chain(*monomers_to_remove))


	
	for rna_row in rna_info:
		if not rna_row['monomerId']:
			rna_row['monomerSet'] = []
		else:
			rna_row['monomerSet'] = [rna_row['monomerId']]

	# add fieldname for 'monomersets'
	fieldnames.append('monomerSet')

	for name in fieldnames:
		if name not in rna_info[0].keys():
			print("The fieldname {} is not in rna_info".format(name))
	for field in rna_info[0].keys():
		if field not in fieldnames:
			print("The field {} is not in fieldnames".format(field))
	'''
	# Create file with standard csv writer -
	# Appropriately adds header but all formatting is off.
	with open(output_file, "w") as f:
		writer = csv.writer(f, delimiter='\t', quotechar = "'", quoting = csv.QUOTE_MINIMAL, lineterminator="\n")
		writer.writerow(['# Generated by {} on {}'.format(__file__, time.ctime())])
		writer.writerow(rna_info[0].keys())
		for rna_row in rna_info:
			# only write monomers that we are keeping.
			if rna_row['geneId'] not in monomers_to_remove:
				writer.writerow(rna_row.values())
		for pc_data in combined_tu_info:
			writer.writerow(combined_tu_info[pc_data].values())
	'''
	# Create file with JSONWriter
	# Formatting is correct, but doesnt add top line.
	with open(output_file, "w") as f:
		writer = JsonWriter(f, fieldnames)
		writer.writeheader()
		for rna_row in rna_info:
			# only write monomers that we are keeping.
			if rna_row['geneId'] not in monomers_to_remove:
				writer.writerow(rna_row)
		for pc_data in combined_tu_info:
			writer.writerow(combined_tu_info[pc_data])


if __name__ == "__main__":
	make_collection()
