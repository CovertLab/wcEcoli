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
from reconstruction.spreadsheets import tsv_reader

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
		'seq' = sequence of all genes in TU including all intergenic regions
		'type' = str indiciated the type of RNA
		'modifiedForms' = list containing modified form each RNA within
			a TU can form.
		'monomerId' = string of all monomerIds in TU separated by '_'
		'comments' = "Transcription unit created within script, for 
			individual RNA comments look at rnas.tsv for that RNA"
		'mw' = mw of the TU sequence
		'location' = location of the first gene in a 
		'id' = pc_gene_id + '_RNA'
		'geneId' = pc_gene_id
		'monomerSet' = pc_monomer_id_list
		'gene_starts_stops' = list of lists recording the start and stop position of each gene in the TU, 
			to be used later within initial_conditions
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
- Allow for certain RNA type mixing (rRNA, tRNA) - handle this type differently
- Allow mass to determine type and location for importing mass fractions in the Parca
- Dont allow for deletion of a gene without incorporation in a TU.

-Import helper functions from another file.
- Make parameters consistent across functions.
- check that directions of genes in a TU match
-fix furure warning
- Combine with file to make the polycistrons file, add flag so user can stil make their own.

'''


def parse_args():
	"""Return error message if argument is not given"""
	parser = argparse.ArgumentParser()
	arser = argparse.ArgumentParser()
	parser.add_argument('filename')
	args = parser.parse_args()
	return args.filename

def parse_tsv(tsv_file):
	"""
	Args:
		tsv_file: flat file to read

	Returns:
		list of lists of the rows contained in the tsv, column names
	"""
	with tsv_reader(tsv_file) as reader:
		tsv_fieldnames = reader.fieldnames
		tsv_list = list(reader)
		return tsv_list, tsv_fieldnames

def load_sequence(sequence):
	"""
	Args:
		sequence: Path to genomic sequence fasta file, who is formatted as a single line.

	Returns:
		List containing a string of the genomic sequence.
	"""
	with open(sequence, "r") as f:
		genome_sequence = f.readlines()[1:]
	return genome_sequence


DIALECT = "excel-tab"
JsonReader = partial(spreadsheets.JsonReader, dialect = DIALECT)
JsonWriter = partial(spreadsheets.JsonWriter, dialect = DIALECT)

# File paths to all necessary flat files.
FLAT_DIR = os.path.join('reconstruction', 'ecoli', 'flat')
RNA_FILE = os.path.join(FLAT_DIR, 'rnas.tsv')
RNA_HALF_LIVES_FILE = os.path.join(FLAT_DIR, 'rna_half_lives.tsv')
GENES_FILE = os.path.join(FLAT_DIR, 'genes.tsv')
#POLY_CISTRON_FILE = os.path.join(FLAT_DIR, 'polycistronic_mrnas_in_model.tsv')
POLY_CISTRON_FILE = parse_args()
GENOME_SEQUENCE_FILE = os.path.join(FLAT_DIR, 'flattened_sequence.fasta')
km_file = os.path.join('fixtures', 'endo_km', 'km3.cPickle')
RNA_SEQ_FILE = os.path.join(FLAT_DIR, 'rna_seq_data', 'rnaseq_rsem_tpm_mean.tsv')
PROTEIN_FILE = os.path.join(FLAT_DIR, 'proteins_old.tsv')
TF_COND_FILE = os.path.join(FLAT_DIR, 'condition','tf_condition_old.tsv')

# output files
TU_FILE = os.path.join(FLAT_DIR, 'operon_rnas.tsv')
TU_HALF_LIVES_FILE = os.path.join(FLAT_DIR, 'operon_rnas_half_lives.tsv')
output_tu_counts = os.path.join(FLAT_DIR, "transcription_units.tsv")
output_gene_tu_matrix = os.path.join(FLAT_DIR, "gene_to_tu_matrix.tsv")
output_proteins = os.path.join(FLAT_DIR, "proteins.tsv")
output_tf_conditions = os.path.join(FLAT_DIR, 'condition', 'tf_condition.tsv')

CONDITION = 'M9 Glucose minus AAs'
SPLIT_DELIMITER = '_'

# Read flat file data
RNA_INFO, fieldnames = parse_tsv(RNA_FILE)
RNA_HALF_LIVES_INFO, rna_hl_fieldnames = parse_tsv(RNA_HALF_LIVES_FILE)
PC_INFO, pc_fieldnames = parse_tsv(POLY_CISTRON_FILE)
GENE_INFO, g_fieldnames = parse_tsv(GENES_FILE)
PROTEIN_INFO, protein_fieldnames = parse_tsv(PROTEIN_FILE)
TF_INFO, TF_FIELDNAMES = parse_tsv(TF_COND_FILE)
GENOMIC_SEQUENCE = load_sequence(GENOME_SEQUENCE_FILE)
rna_seq_data_all_cond = parse_tsv(RNA_SEQ_FILE)


def find_tu_type(tu_genes_info):
	"""
	Purpose:
	Need to assign an rna type to the pcic RNA (e.g. rRNA,
	mRNA...)

	If all the RNA types in a TU are the same, the type will be assigned
	to that type.

	If they are not the same, the most prevalent one will be assigned
	for the whole TU, or if all types are equally represented, the
	type from the first RNA in the TU will be assigned.

	Args:
		tu_genes_info: Dictionary of transcription unit info for target transcription unit.

	Returns:
		RNA type for the transcription unit.
	"""
	mismatch_output = """Types of RNA's in transcription unit {} dont 
	match, please double check that your transcription unit is correct. 
	Might have to do additional checks. The naming for this transcription 
	unit will follow the most prevalant RNA type, or the first RNA in the TU 
	if there is an even number of RNAs."""

	tu_types = [tu_genes_info[gene]['type'] for gene in tu_genes_info]
	if len(set(tu_types)) > 1:
		warnings.warn(mismatch_output.format('_'.join(tu_genes_info.keys())))
	tu_type = max(tu_types, key=Counter(tu_types).get)
	return tu_type

def find_tu_location(tu_genes_info):
	"""
	Mimics structure and output of find_tu_type but comparing locations
	of the TU instead.
	Args:
		tu_genes_info: Dictionary of transcription unit info for target transcription unit.

	Returns:
		List containing a unicode str of the location.
	"""
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

def get_tu_sequence(tu_genes_info, first_gene, last_gene):
	"""
	Purpose:
	Obtain the RNA sequence for the entire transcription unit.
	Cannot simply concatenate component individual RNA sequences since
	they would not include intergenic regions.

	Args:
		tu_genes_info: Dictionary containing information for a transcription unit (Dict
		contains direction and coordinate info).
		first_gene: First gene in the transcription unit
		last_gene: Last gene in the transcription unit

	Returns:
		The rna sequence for the transcription unit.

	Note:
	This was developed with only implementing mRNAs for now. If
	incorporating untranslated RNAs (like rRNAs or tRNAs) in the future
	would have to add a process for RNA processing.
	"""
	direction = tu_genes_info[first_gene]['direction']
	rna_sequence = []
	if direction == '+':
		first_gene_coordinate = tu_genes_info[first_gene]['chromosome_coordinate']
		last_gene_coordinate = tu_genes_info[last_gene]['chromosome_coordinate'] + tu_genes_info[last_gene]['length']
		sequence = Seq(GENOMIC_SEQUENCE[0][first_gene_coordinate:last_gene_coordinate], IUPAC.unambiguous_dna)
		rna_sequence = str(sequence.transcribe())
	elif direction == '-':
		first_gene_coordinate = tu_genes_info[first_gene]['chromosome_coordinate'] +1
		last_gene_coordinate = tu_genes_info[last_gene]['chromosome_coordinate'] - tu_genes_info[last_gene]['length'] + 1
		sequence = Seq(GENOMIC_SEQUENCE[0][last_gene_coordinate:first_gene_coordinate], IUPAC.unambiguous_dna)
		rna_sequence = str(sequence.reverse_complement().transcribe())
	else:
		print("For some reason your gene hasnt been assigned a direction")

	return rna_sequence

def calculate_rna_biomass(sequence):
	"""
	Purpose: Calculate the RNA biomass of a particular RNA sequence.

	Args:
		sequence: RNA sequence as a str.

	Returns:
		Mass of the RNA sequence as a float.
	"""
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

	return [0.0, 0.0, 0.0, np.sum(list(counts.values())) + ppi_mass, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

def count_ntps_rna(sequence):
	"""
	Args:
		sequence: RNA sequence.

	Returns:
		Counts of the nucleotidees in the sequence
	"""
	return [sequence.count('A'), sequence.count('C'),
			sequence.count('G'), sequence.count('U')]

def gather_tu_genes_info():
	"""
	Purpose: Gather data for each rna within each TU.
	PC_INFO (imported at the top-level):
		A list of dictionaries, where each row is a new polycistronic TU
		that needs to be added to the model.
		- 'transcription_units': list of genes in the TU.
		- 'monomers_to_remove': list of monocistronic genes from the
			TU to remove from the model.

	Returns:
		A nested dictionary of all the TUs to be added to the model.
		First level - and information for all the genes within that
		TU. This data will used in a later function to define all the
		features of a specific TU.
	"""
	tu_genes_info = {}
	for pc in PC_INFO:
		pc_gene_id = '_'.join(pc['transcription_units'])
		tu_genes_info[pc_gene_id] = {}
		for rna in pc['transcription_units']:
			tu_genes_info[pc_gene_id][rna] = {}
			rnaId = rna + "_RNA"  # by default, rnaId is set to geneId_RNA (for half lives file)
			for rna_row in RNA_INFO:
				if rna_row['geneId'] == rna:
					rnaId = rna_row['id']  # rnaId updated based on rnas.tsv, if suffix is not _RNA
					tu_genes_info[pc_gene_id][rna]['type'] = rna_row['type']
					tu_genes_info[pc_gene_id][rna]['modifiedForms'] = rna_row['modifiedForms']
					tu_genes_info[pc_gene_id][rna]['monomerId'] = rna_row['monomerId']
					tu_genes_info[pc_gene_id][rna]['location'] = rna_row['location']
			for rna_hl_row in RNA_HALF_LIVES_INFO:  # create separate half lives file
				if rna_hl_row['id'] == rnaId:
					tu_genes_info[pc_gene_id][rna]['halfLife'] = rna_hl_row['half_life']
			for gene_row in GENE_INFO:
				if gene_row['id'] == rna:
					tu_genes_info[pc_gene_id][rna]['chromosome_coordinate'] = gene_row['coordinate']
					tu_genes_info[pc_gene_id][rna]['direction'] = gene_row['direction']
					tu_genes_info[pc_gene_id][rna]['length'] = gene_row['length']
	return tu_genes_info

def find_gene_starts_stops(pc_gene_id, pc_gene_info):
	# getting the genes in the PC this way to preserve the proper order
	genes_in_pc = pc_gene_id.split('_')
	gene_starts_stops = []
	for idx, gene in enumerate(genes_in_pc):
		if idx == 0:
			start = 0
			stop = pc_gene_info[gene]['length'] - 1
			gene_starts_stops.append([start, stop])
		else:
			if pc_gene_info[gene]['direction'] == '+':
				start = pc_gene_info[gene]['chromosome_coordinate'] - pc_gene_info[genes_in_pc[idx-1]]['chromosome_coordinate'] - 1
			elif pc_gene_info[gene]['direction'] == '-':
				start = pc_gene_info[genes_in_pc[idx-1]]['chromosome_coordinate'] - pc_gene_info[gene]['chromosome_coordinate'] - 1
			stop = start + pc_gene_info[gene]['length']
			gene_starts_stops.append([start, stop])
	return gene_starts_stops

def gather_tu_info(tu_genes_info):
	"""
	Purpose:
	Need to gather all the info for each TU we need to add it to the
	model. Format closely follows rnas.tsv with modifications outlined
	at the top of this file.

	Args:
		tu_genes_info: Dictionary created within gather_tu_genes_info, nested
		dictionary containing info for each gene within each TU.

	Returns:
		- A nested dictionary. First level keys are the gene ids for each
		polycistron being added to the model, nested keys are all the columns
		that will be needed within the operon_rnas.tsv file. For more detail on each
		key, please see the comment at the top of this file. Naming mimics rnas.tsv.
		- A nested dictionary specifically for RNA half lives. First level keys
		are the gene ids for each polycistron being added to the model, nested keys
		are RNA ids and half lives. Naming mimics rna_half_lives.tsv.
	"""
	tu_info = {}
	tu_half_life_info = {}

	for pc in PC_INFO:
		# Get the gene ID of the polycistron
		pc_gene_id = '_'.join(pc['transcription_units'])

		# Find the monomerIds, names and modifiedForms of all the RNAs in the TU
		pc_monomer_id_list = []
		pc_name_list = []
		pc_modified_forms_list = []
		for rna in pc['transcription_units']:
			pc_monomer_id_list.append(tu_genes_info[pc_gene_id][rna]['monomerId'])
			# pc_name_list.append(tu_genes_info[pc_gene_id][rna]['name'])
			pc_modified_forms_list.append(tu_genes_info[pc_gene_id][rna]['modifiedForms'])

		# Find the first and last genes in the TU
		first_gene = pc['transcription_units'][0]
		last_gene = pc['transcription_units'][-1]

		# Construct dictionary of transcription unit info for each tu.
		tu_info[pc_gene_id] = {}
		tu_info[pc_gene_id]['seq'] = get_tu_sequence(tu_genes_info[pc_gene_id],
			first_gene, last_gene)
		tu_info[pc_gene_id]['type'] = find_tu_type(tu_genes_info[pc_gene_id])
		tu_info[pc_gene_id]['modifiedForms'] = []
		tu_info[pc_gene_id]['monomerId'] = '{}{}'.format(
											'_'.join([x.replace('-MONOMER', '')
											for x in pc_monomer_id_list]), '-MONOMER')
		tu_info[pc_gene_id]['comments'] = """Transcription unit created within script, for individual RNA comments look at rnas.tsv for that RNA"""
		tu_info[pc_gene_id]['mw'] = calculate_rna_biomass(tu_info[pc_gene_id]['seq'])
		tu_info[pc_gene_id]['location'] = find_tu_location(tu_genes_info[pc_gene_id])
		if not tu_info[pc_gene_id]['type'] == 'mRNA':
			tu_info[pc_gene_id]['id'] = pc_gene_id + '_' + tu_info[pc_gene_id]['type'].upper()
		else:
			tu_info[pc_gene_id]['id'] = pc_gene_id + '_RNA'
		tu_info[pc_gene_id]['geneId'] = pc_gene_id
		tu_info[pc_gene_id]['monomerSet'] = pc_monomer_id_list
		tu_info[pc_gene_id]['gene_starts_stops'] = find_gene_starts_stops(pc_gene_id, tu_genes_info[pc_gene_id])

		tu_half_life_info[pc_gene_id] = {}
		tu_half_life_info[pc_gene_id]['id'] = pc_gene_id + '_RNA'
		try:
			tu_half_life_info[pc_gene_id]['half_life'] = tu_genes_info[pc_gene_id][first_gene]['halfLife']
		except:
			#for now dont include a half life if the first gene does not have one, and let it be solved for in the parca
			del tu_half_life_info[pc_gene_id]
	return tu_info, tu_half_life_info


def find_monomers_to_remove():
	"""
	Purpose: Find all the momnomers that should not be included in
	the final operon_rnas.tsv file.
	Returns: Set of unique monomers to not include in final tsv.
	"""
	return set(itertools.chain(*[pc['monomers_to_remove']
		for pc in PC_INFO]))

def update_rna_info():
	"""
	Purpose:
	Add an additional key for each RNA within RNA_INFO specificing the
	monomers assigned to that mRNA.
	RNA_INFO imported at the top level.

	Returns:
		List of Dictionaries for all the RNA_INFO but this time including a key
		for monomerSet which contains a list of all the monomers assigned to a
		particular mRNA, list will be empty if there are no monomers to assign to an RNA.
	"""
	for rna_row in RNA_INFO:
		if not rna_row['monomerId'] or rna_row['monomerId'] == "null":
			rna_row['monomerSet'] = []
		else:
			rna_row['monomerSet'] = [rna_row['monomerId']]
		for gene in GENE_INFO:
			if rna_row['id'] == gene['rnaId']:
				rna_row['gene_starts_stops'] = [[0, gene['length']-1]]
	# add fieldname for 'monomersets'
	fieldnames.append('monomerSet')
	fieldnames.append('gene_starts_stops')

def write_output_file(tu_info, tu_half_life_info, monomers_to_remove):
	"""
	Construct a tsv file that mimics the structure and formatting
	of rnas.tsv.
	"""
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

	# Write separate operon_rnas_half_lives file
	with open(TU_HALF_LIVES_FILE, "w") as f:
		f.write('# Generated by {} on {} \n'.format(__file__, time.ctime()))

		writer = JsonWriter(f, rna_hl_fieldnames)
		writer.writeheader()

		for rna_hl_row in RNA_HALF_LIVES_INFO:
			# need to check if monomer id is contained (not ==) in RNA id, since monomer id does not end with '_RNA'
			if not any(monomer in rna_hl_row['id'] for monomer in monomers_to_remove):
				writer.writerow(rna_hl_row)
		for pc_data in tu_info:
			try:
				writer.writerow(tu_half_life_info[pc_data])
			except:
				pass
	return

def make_operon_rnas_file():
	tu_genes_info = gather_tu_genes_info()
	tu_info, tu_half_life_info = gather_tu_info(tu_genes_info)

	monomers_to_remove = find_monomers_to_remove()
	update_rna_info()

	write_output_file(tu_info, tu_half_life_info, monomers_to_remove)


def create_gene_to_tu_matrix(rna_info, tu_info):
	"""
	Purpose: Parse tsv files then create a gene to tu matrix mapping genes
	in each TU to its partner RNA.

	Args:
		rna_info: Parsed rna data file.
		tu_info: Parsed tu data file.

	Returns:
		Sparse numpy matrix, mapping TU to rna's. 0 = no mapping; 1 = mapping.

	TODO: Add a list containing the geneIds in a TU, can pull this in to genes_in_tu.
	"""

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
	"""
	The counts vector is not in the same order as the Gene_TU_Matrix.
	Need to reoder and pull out count information.

	Gathers information needed based on the condition the model is
	being run in.
	"""
	#import pdb; pdb.set_trace()
	rna_seq_data_index = {
		row['Gene']: row[CONDITION]
		for row in rna_seq_data_all_cond[0]}

	rna_seq_counts_vector = [
		rna_seq_data_index[gene]
		for gene in rnas_gene_order]

	return rna_seq_counts_vector

def create_tu_counts_vector(gene_tu_matrix, rna_seq_counts_vector, tu_info):
	"""
	Purpose:
		Calculate transcription unit counts. Will functionally replace how rna seq
		counts were used in the model. If no TUs are added it will exactly match
		rna-seq counts. Using scipy.optimize..nnls solver to give a non-negative least squares fit.

	Args:
		gene_tu_matrix: Sparse numpy matrix, mapping TU to rna's. 0 = no mapping; 1 = mapping.
		rna_seq_counts_vector: A vector of RNA seq counts for the condition we are interested in.
		tu_info: list of dictionaries containing all the information we have for all the
		transcription units we want to include in the model.

	Returns:
		A list of dictionaries. Each dictionary contains the 'tu_id' and the 'tu_count' for
		each TU in the model.
	"""
	tu_counts_vector = nnls(gene_tu_matrix, rna_seq_counts_vector)[0]
	tu_gene_order = [row['geneId'] for row in tu_info]
	tu_genes_counts = []

	for i in range(0, len(tu_gene_order)):
		tu_genes_counts.append({'tu_id': tu_gene_order[i], 'tu_count': tu_counts_vector[i]})

	return tu_genes_counts


def make_transcription_units_file():
	"""
	Run all the functions necessary to get the transcription counts vector
	output as transcription_units.tsv.
	"""
	tu_info, fieldnames_rna = parse_tsv(TU_FILE)
	gene_tu_matrix, rnas_gene_order = create_gene_to_tu_matrix(RNA_INFO, tu_info)
	rna_seq_counts_vector = create_rnaseq_count_vector(rnas_gene_order)
	tu_counts = create_tu_counts_vector(gene_tu_matrix, rna_seq_counts_vector, tu_info)
	tu_fieldnames = ['tu_id', 'tu_count']

	with open(output_tu_counts, "w") as f:
		writer = JsonWriter(f, tu_fieldnames)
		writer.writeheader()
		for tu_count in tu_counts:
			writer.writerow(tu_count)

	with open(output_gene_tu_matrix, "w") as f:
		writer = csv.writer(f, delimiter=' ')
		for row in gene_tu_matrix:
			writer.writerow(row)

def make_new_proteins_file(output_file):
	rna_info, rna_fieldnames = parse_tsv(TU_FILE)

	#Go through monomerSet line by line. Find the matching monomers within
	#those lists then find the corresponding monomer in proteins.tsv.
	#Add the id from operon_rnas to the rnaSet list

	protein_index = {}
	for protein_row in PROTEIN_INFO:
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
		for protein_row in PROTEIN_INFO:
			writer.writerow(protein_row)

	#breakpoint()

def remove_kms_file(km_file):
	"""
	Purpose:
	The parca checks if the km's have already been calculated, if they have then
	it does not calculate them again and just pulls them from the km.cpickle file.
	However, if the operon_rnas file gets modified the kms need to be recalculated.
	This function will remove the .cpickle file and force the parca to re-calculate the
	Kms, this way we can be sure we are working with a correct dataset.
	"""
	if os.path.exists(km_file):
		os.remove(km_file)
	return

def make_new_tf_conditions_file(output_file):
	protein_info, _ = parse_tsv(output_proteins)

	rnaIdToRnaSet = {x["rnaId"]: x["rnaSet"] for x in protein_info}
	fields_to_modify = ["active genotype perturbations", "inactive genotype perturbations"]

	for row in TF_INFO:
		for field in fields_to_modify:
			if len(row[field]) > 0:
				genotype = {}
				for key, val in row[field].items():
					for new_key in rnaIdToRnaSet[key.strip("[c]")]:
						genotype[new_key + "[c]"] = val
				row[field] = genotype
	with open(output_file, "w") as f:
		writer = JsonWriter(f, TF_FIELDNAMES)
		writer.writeheader()
		for tf_row in TF_INFO:
			writer.writerow(tf_row)


if __name__ == "__main__":
	make_operon_rnas_file()
	make_transcription_units_file()
	remove_kms_file(km_file)
	make_new_proteins_file(output_proteins)
	make_new_tf_conditions_file(output_tf_conditions)

