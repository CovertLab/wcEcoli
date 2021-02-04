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
from wholecell.utils.fast_nonnegative_least_squares import fast_nnls

'''
Purpose:
This file autogenerates several files used to allow polycistronic transcription
units in the model (operon_rnas.tsv - functionally replaces rnas.tsv; gene_tu_matrix.tsv;
and transcription_units.tsv). 


Inputs:
- polycistronic_mrnas_in_the_model:
	Supplies the polycistronic mRNAs to add as a list of gene_ids, 
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
		an additional column 'monomer_set', that allows mapping of mRNAs
		back to their componenet monomer - this allows more than one
		monomer to map to a mRNA). 

		-Values are defined as follows:
		'seq' = sequence of all genes in TU including all intergenic regions
		'type' = str indiciated the type of RNA
		'modifiedForms' = list containing modified form each RNA within
			a TU can form.
		#'monomer_id' = string of all monomer_ids in TU separated by '_'
		'comments' = "Transcription unit created within script, for 
			individual RNA comments look at rnas.tsv for that RNA"
		'mw' = mw of the TU sequence
		'location' = location of the first gene in a 
		'id' = pc_gene_id + '_RNA'
		'gene_id' = pc_gene_id
		'monomer_set' = pc_monomer_id_list
		'gene_starts_stops' = list of lists recording the start and stop position of each gene in the TU, 
			to be used later within initial_conditions
	transcription_units.tsv
		- Counts of each transcription unit included in the model. Functionally replaces
		rna seq counts in the model, at least for the basal case. See todos for future
		upgrades to this file
	operon_half_.tsv
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

- is rnas_genes_order used?
- update all documentation for changes that needed to happen due to changes in the flat file
- Delete functions that are no longer necessary.
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

def strip_units(parsed_file_output):
	# hack for stripping units off any field that contains units.
	# hacky way to look for units too
	parsed_tsv, parsed_fieldnames = parsed_file_output
	fields_to_strip_units = [idx
		for idx, field in enumerate(parsed_fieldnames)
			if '(' in field
		]
	for field_idx in fields_to_strip_units:
		unitless_fieldname = parsed_fieldnames[field_idx].split(' ')[0]
		for row in parsed_tsv:
			row[unitless_fieldname] = int(re.sub("[^0-9]", "", str(row[unitless_fieldname])))
			row[parsed_fieldnames[field_idx]] = row.pop(unitless_fieldname)
		#unit_keys = parsed_fieldnames[field_idx].split(' ')[0]
	return parsed_tsv, parsed_fieldnames


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
PROTEIN_FILE = os.path.join(FLAT_DIR, 'proteins.tsv')
PROTEIN_OLD_FILE = os.path.join(FLAT_DIR, 'proteins_old.tsv')
TF_COND_FILE = os.path.join(FLAT_DIR, 'condition','tf_condition_old.tsv')

# output files
TU_FILE = os.path.join(FLAT_DIR, 'operon_rnas.tsv')
TU_HALF_LIVES_FILE = os.path.join(FLAT_DIR, 'operon_rnas_half_lives.tsv')
output_tu_counts = os.path.join(FLAT_DIR, "transcription_units.tsv")
output_proteins = os.path.join(FLAT_DIR, "proteins.tsv")
output_tf_conditions = os.path.join(FLAT_DIR, 'condition', 'tf_condition.tsv')
CONDITION = 'M9 Glucose minus AAs'
SPLIT_DELIMITER = '_'

# Read flat file data
RNA_INFO, fieldnames = parse_tsv(RNA_FILE)
RNA_HALF_LIVES_INFO, rna_hl_fieldnames = strip_units(parse_tsv(RNA_HALF_LIVES_FILE))
RNA_HALF_LIVES_DICTIONARY = {
	rna['id'].replace('_RNA', '').replace('-tRNA','').replace('-RNA',''): rna['half_life (units.s)'] 
	for rna in RNA_HALF_LIVES_INFO}
PC_INFO, pc_fieldnames = parse_tsv(POLY_CISTRON_FILE)
GENE_INFO, g_fieldnames = parse_tsv(GENES_FILE)
PROTEIN_INFO, protein_fieldnames = parse_tsv(PROTEIN_FILE)
PROTEIN_OLD_INFO, protein_old_fieldnames = parse_tsv(PROTEIN_OLD_FILE)
TF_INFO, TF_FIELDNAMES = parse_tsv(TF_COND_FILE)
GENOMIC_SEQUENCE = load_sequence(GENOME_SEQUENCE_FILE)
rna_seq_data_all_cond = parse_tsv(RNA_SEQ_FILE)

'''
How to calculate Half-Life of a polycistronic mRNA:
1. For a multi-gene operon, take the value of the first gene. If the first gene does not 
have an experimentally determined half-life, let the average be used.
2. For a multi-gene operon, take the value of the first gene. If the first gene does not 
have an experimentally determined half-life, take the value of the second gene, and so on.
If no genes have half life, allow to use the average.
3. Take the average of all available half-lives.

TODO: Add half life calc value as an argument
'''
HALF_LIFE_CALC = 3


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
	rna_to_gene_map = {gene['rna_id']: gene['id'] for gene in GENE_INFO}
	for pc in PC_INFO:
		pc_gene_id = '_'.join(pc['transcription_units'])
		first_gene = pc['transcription_units'][0] + '_RNA'
		tu_genes_info[pc_gene_id] = {}
		for rna in pc['transcription_units']:
			tu_genes_info[pc_gene_id][rna] = {}
			rna_id = rna + "_RNA"  # by default, rnaId is set to gene_id_RNA (for half lives file)
			# make dictionary mapping rna_id to gene_id
			for rna_row in RNA_INFO:
				if rna_to_gene_map[rna_row['id']] == rna:
					rnaId = rna_row['id']  # rnaId updated based on rnas.tsv, if suffix is not _RNA
					tu_genes_info[pc_gene_id][rna]['type'] = rna_row['type']
					tu_genes_info[pc_gene_id][rna]['modified_forms'] = rna_row['modified_forms']
					tu_genes_info[pc_gene_id][rna]['monomer_id'] = rna_row['monomer_id']
			# Remove this region?
			'''
			if rna_id in RNA_HALF_LIVES_DICTIONARY:
				tu_genes_info[pc_gene_id][rna]['halfLife'] = RNA_HALF_LIVES_DICTIONARY[rna_id]
			'''
	return tu_genes_info

def calculate_half_life(pc_gene_id, pc):
	'''
	How to calculate Half-Life of a polycistronic mRNA:
	1. For a multi-gene operon, take the value of the first gene. If the first gene does not 
	have an experimentally determined half-life, let the average be used.
	2. For a multi-gene operon, take the value of the first gene. If the first gene does not 
	have an experimentally determined half-life, take the value of the second gene, and so on.
	If no genes have half life, allow to use the average.
	3. Take the average of all available half-lives.
	'''
	if HALF_LIFE_CALC == 1:
		first_gene = pc['transcription_units'][0]
		try:
			half_life = RNA_HALF_LIVES_DICTIONARY[first_gene]
			# Update Dictionary to have new half life info
			RNA_HALF_LIVES_DICTIONARY[pc_gene_id] = half_life
		except:
			#for now dont include a half life if the first gene does not have one, and let it be solved for in the parca
			half_life = False
	if HALF_LIFE_CALC == 2:
		half_lives = []
		for gene in pc['transcription_units']:
			try:
				half_lives.append(RNA_HALF_LIVES_DICTIONARY[gene])
			except:
				pass
		if half_lives:
			half_life = half_lives[0]
			RNA_HALF_LIVES_DICTIONARY[pc_gene_id] = half_life
		else: 
			half_life = False
	if HALF_LIFE_CALC == 3:
		half_lives = []
		for gene in pc['transcription_units']:
			try:
				half_lives.append(RNA_HALF_LIVES_DICTIONARY[gene])
			except:
				pass
		if half_lives:
			half_life = np.mean(half_lives)
			RNA_HALF_LIVES_DICTIONARY[pc_gene_id] = half_life
		else: 
			half_life = False
	if HALF_LIFE_CALC == 4:
		half_lives = []
		for gene in pc['transcription_units']:
			try:
				half_lives.append(RNA_HALF_LIVES_DICTIONARY[gene])
			except:
				pass
		if half_lives:
			half_life = np.max(half_lives)
			RNA_HALF_LIVES_DICTIONARY[pc_gene_id] = half_life
		else: 
			half_life = False
	if HALF_LIFE_CALC == 5:
		half_lives = []
		for gene in pc['transcription_units']:
			try:
				half_lives.append(RNA_HALF_LIVES_DICTIONARY[gene])
			except:
				half_lives.append(344.50836653386455)
		if half_lives:
			half_life = np.mean(half_lives)
			RNA_HALF_LIVES_DICTIONARY[pc_gene_id] = half_life
		else: 
			half_life = False
	return half_life

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

		# Find the monomer_ids, names and modified_forms of all the RNAs in the TU
		pc_monomer_id_list = []
		pc_name_list = []
		pc_modified_forms_list = []
		for rna in pc['transcription_units']:
			pc_monomer_id_list.append(tu_genes_info[pc_gene_id][rna]['monomer_id'])
			# pc_name_list.append(tu_genes_info[pc_gene_id][rna]['name'])
			pc_modified_forms_list.append(tu_genes_info[pc_gene_id][rna]['modified_forms'])

		# Find the first and last genes in the TU
		first_gene = pc['transcription_units'][0]
		last_gene = pc['transcription_units'][-1]

		# Construct dictionary of transcription unit info for each tu.
		tu_info[pc_gene_id] = {}
		#tu_info[pc_gene_id]['seq'] = get_tu_sequence(tu_genes_info[pc_gene_id],
		#	first_gene, last_gene)
		tu_info[pc_gene_id]['type'] = find_tu_type(tu_genes_info[pc_gene_id])
		tu_info[pc_gene_id]['modified_forms'] = []
		#tu_info[pc_gene_id]['monomer_id'] = '{}{}'.format(
											#'_'.join([x.replace('-MONOMER', '')
											#for x in pc_monomer_id_list]), '-MONOMER')
		#tu_info[pc_gene_id]['comments'] = """Transcription unit created within script, for individual RNA comments look at rnas.tsv for that RNA"""
		#tu_info[pc_gene_id]['mw'] = calculate_rna_biomass(tu_info[pc_gene_id]['seq'])
		#tu_info[pc_gene_id]['location'] = find_tu_location(tu_genes_info[pc_gene_id])
		if not tu_info[pc_gene_id]['type'] == 'mRNA':
			tu_info[pc_gene_id]['id'] = pc_gene_id + '_' + tu_info[pc_gene_id]['type'].upper()
		else:
			tu_info[pc_gene_id]['id'] = pc_gene_id + '_RNA'
		tu_info[pc_gene_id]['gene_set'] = pc['transcription_units']
		tu_info[pc_gene_id]['monomer_set'] = pc_monomer_id_list
		#tu_info[pc_gene_id]['gene_starts_stops'] = find_gene_starts_stops(pc_gene_id, tu_genes_info[pc_gene_id])

		# Half-Life Info
		# Note: Will assign average half life if first gene does not have an associated half life.
		tu_half_life_info[pc_gene_id] = {}
		tu_half_life_info[pc_gene_id]['id'] = pc_gene_id + '_RNA'
		half_life = calculate_half_life(pc_gene_id, pc)
		# might not need this part, just save the new dictionary thats being created.
		if half_life:
			tu_half_life_info[pc_gene_id]['half_life (units.s)'] = half_life
		else:
			del tu_half_life_info[pc_gene_id]

		#for now dont include a half life if the first gene does not have one, and let it be solved for in the parca
			
		'''
		try:
			tu_half_life_info[pc_gene_id]['half_life (units.s)'] = tu_genes_info[pc_gene_id][first_gene]['halfLife']
		except:
			#for now dont include a half life if the first gene does not have one, and let it be solved for in the parca
			del tu_half_life_info[pc_gene_id]
		'''
	return tu_info, tu_half_life_info

def find_monomers_to_remove():
	"""
	Purpose: Find all the monoomers that should not be included in
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
		for monomer_set which contains a list of all the monomers assigned to a
		particular mRNA, list will be empty if there are no monomers to assign to an RNA.
	"""
	rna_to_gene_map = {gene['rna_id']: gene['id'] for gene in GENE_INFO}

	for rna_row in RNA_INFO:
		if not rna_row['monomer_id'] or rna_row['monomer_id'] == "null":
			rna_row['monomer_set'] = []
		else:
			rna_row['monomer_set'] = [rna_row['monomer_id']]
		rna_row.pop('monomer_id')
		rna_row['gene_set'] = [rna_to_gene_map[rna_row['id']]]
		#for gene in GENE_INFO:
			#if rna_row['id'] == gene['rna_id']:
				#rna_row['gene_starts_stops'] = [[0, gene['length']-1]]
	# add fieldname for 'monomer_set'
	fieldnames.append('gene_set')
	fieldnames.append('monomer_set')
	fieldnames.remove('monomer_id')

	#fieldnames.append('gene_starts_stops')

def write_output_file(tu_info, tu_half_life_info, monomers_to_remove):
	"""
	Construct a tsv file that mimics the structure and formatting
	of rnas.tsv.
	"""
	rna_to_gene_map = {gene['rna_id']: gene['id'] for gene in GENE_INFO}
	# Create file with JSONWriter
	with open(TU_FILE, "w") as f:
		# Add comment line with file creation info
		f.write('# Generated by {} on {} \n'.format(__file__, time.ctime()))

		# Write file with JsonWriter
		writer = JsonWriter(f, fieldnames)
		writer.writeheader()
		for rna_row in RNA_INFO:
			# only write monomers that we are keeping.
			if rna_to_gene_map[rna_row['id']] not in monomers_to_remove:
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
			if not any(monomer + '_RNA' == rna_hl_row['id'] for monomer in monomers_to_remove):
				writer.writerow(rna_hl_row)
		for pc_data in tu_info:
			#breakpoint()
			try:
				writer.writerow(tu_half_life_info[pc_data])
			except:
				# Note: Need to do this try:except since sometimes the operon is not given a half life if the first gene 
				# is not originally assigned a half life.
				#breakpoint()
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

	TODO: Add a list containing the gene_ids in a TU, can pull this in to genes_in_tu.
	"""

	num_rnas = len(rna_info)
	num_tus = len(tu_info)

	gene_to_tu_matrix = np.zeros((num_rnas, num_tus))
	rnas_gene_order = [row['gene_set'][0] for row in rna_info]
	reverse_index = {
		row['gene_set'][0]: gene_index
		for gene_index, row in enumerate(rna_info)}

	for index, tu in enumerate(tu_info):
		if len(tu['gene_set']) == 1:
			#if tu['gene_set'][0] == "G0-8903":
				#breakpoint()
			gene = tu['gene_set'][0]
			gene_index = reverse_index[gene]
			gene_to_tu_matrix[gene_index, index] = 1
		else:
			for gene in tu['gene_set']:
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
	rna_seq_data_index = {
		row['Gene']: row[CONDITION]
		for row in rna_seq_data_all_cond[0]}

	rna_seq_counts_vector = [
		rna_seq_data_index[gene]
		for gene in rnas_gene_order]
	return np.array(rna_seq_counts_vector)

def create_tu_counts_vector(gene_tu_matrix, rna_seq_counts_vector, tu_info):
	"""
	Purpose:
		Calculate transcription unit counts. Will functionally replace how rna seq
		counts were used in the model. If no TUs are added it will exactly match
		rna-seq counts. Using scipy.optimize.nnls solver to give a non-negative least squares fit.

	Args:
		gene_tu_matrix: Sparse numpy matrix, mapping TU to rna's. 0 = no mapping; 1 = mapping.
		rna_seq_counts_vector: A vector of RNA seq counts for the condition we are interested in.
		tu_info: list of dictionaries containing all the information we have for all the
		transcription units we want to include in the model.

	Returns:
		A list of dictionaries. Each dictionary contains the 'tu_id' and the 'tu_count' for
		each TU in the model.

	TODO: enumerate for the i in range stuff below
	"""
	tu_counts_vector, rnorm = fast_nnls(gene_tu_matrix, rna_seq_counts_vector)

	tu_gene_order = [('_').join(row['gene_set']) for row in tu_info]
	tu_genes_counts = []

	for i, val in enumerate(tu_gene_order):
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
	tu_counts_vector = np.array([x['tu_count'] for x in tu_counts])

	def get_monocistronic_gene_indexes(gene_tu_matrix):
		indexes = []

		for i in np.where(gene_tu_matrix.sum(axis=1) == 1)[0]:
			TU_index = np.where(gene_tu_matrix[i, :] == 1)[0]
			if gene_tu_matrix[:, TU_index].sum() == 1:
				indexes.append(i)

		return np.array(indexes)

	import matplotlib.pyplot as plt
	mc_gene_indexes = get_monocistronic_gene_indexes(gene_tu_matrix)
	rna_counts_new = gene_tu_matrix.dot(tu_counts_vector)
	rna_counts_old = rna_seq_counts_vector

	plt.figure(figsize=(12, 12))
	plt.scatter(rna_counts_old[~mc_gene_indexes], rna_counts_new[~mc_gene_indexes], s=3, c='b', label='poly')
	plt.scatter(rna_counts_old[mc_gene_indexes], rna_counts_new[mc_gene_indexes], s=3, c='r', label='mono')
	plt.legend()
	plt.xlabel('Original counts')
	plt.ylabel('Counts calculated through NNLS')
	plt.savefig('nnlq_rna_counts.pdf')


	tu_fieldnames = ['tu_id', 'tu_count']

	with open(output_tu_counts, "w") as f:
		writer = JsonWriter(f, tu_fieldnames)
		writer.writeheader()
		for tu_count in tu_counts:
			writer.writerow(tu_count)

def make_new_proteins_file(output_file):
	rna_info, rna_fieldnames = parse_tsv(TU_FILE)

	#Go through monomer_set line by line. Find the matching monomers within
	#those lists then find the corresponding monomer in proteins.tsv.
	#Add the id from operon_rnas to the rna_set list

	protein_index = {}
	for protein_row in PROTEIN_INFO:
		protein_row['rna_set'] = []
		#protein_row['gene_id'] = []
		protein_index[protein_row['id']] = protein_row
		for old_protein in PROTEIN_OLD_INFO:
			if old_protein['id'] == protein_row['id']:
				protein_row['gene_id'] = old_protein['geneId']

	for rna_row in rna_info:
		for monomer in rna_row['monomer_set']:
			protein_row = protein_index[monomer]
			protein_row['rna_set'].append(rna_row['id'])
	'''
	for rna_row in rna_info:
		for monomer in rna_row['monomer_set']:
			protein_row = protein_index[monomer]
			#list_of_rnas_in_set = rna_row['id'].split('_')[0:-1]
			protein_row['rna_set'] = [rna + '_RNA' for rna in rna_row['id'].split('_')[0:-1]]
			#protein_row['rna_set'].append(rna_row['id'].split('_')[0:-1])
	for protein in PROTEIN_INFO:
		if type(protein['rna_set'][0]) == list:
			protein['rna_set'] = protein['rna_set'][0]
	'''
	#add fieldname for 'rna_set' if it doesnt already exist.
	if 'rna_set' and 'gene_id' not in protein_fieldnames:
		protein_fieldnames.append('rna_set')
		protein_fieldnames.append('gene_id')

	with open(output_file, "w") as f:
		writer = JsonWriter(f, protein_fieldnames)
		writer.writeheader()
		for protein_row in PROTEIN_INFO:
			writer.writerow(protein_row)


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
	tu_info, _ = parse_tsv(TU_FILE)
	rna_id_to_rna_set = {}
	for monocistronic_rna in RNA_INFO:
		rna_id = monocistronic_rna['id']
		for protein_row in protein_info:
			for rna_in_set in protein_row['rna_set']:
				mono_in_poly = rna_in_set.replace('_RNA', '').replace('-tRNA','').replace('-RNA','').split('_')
				monocis_rnaid = monocistronic_rna['id'].replace('_RNA', '').replace('-tRNA','').replace('-RNA','')
				if len(mono_in_poly) > 1:
					for cistron in mono_in_poly:
						if monocis_rnaid == cistron:
							if rna_id in rna_id_to_rna_set and rna_in_set not in rna_id_to_rna_set[rna_id]:
								rna_id_to_rna_set[rna_id].append(rna_in_set)
							elif rna_id not in rna_id_to_rna_set:
								rna_id_to_rna_set[rna_id] = [rna_in_set]
				elif len(mono_in_poly) == 1:
					if monocis_rnaid == mono_in_poly[0]:
						if rna_id in rna_id_to_rna_set and rna_in_set not in rna_id_to_rna_set[rna_id]:
							rna_id_to_rna_set[rna_id].append(rna_in_set)
						elif rna_id not in rna_id_to_rna_set:
							rna_id_to_rna_set[rna_id] = [rna_in_set]

	fields_to_modify = ["active genotype perturbations", "inactive genotype perturbations"]
	
	for row in TF_INFO:
		for field in fields_to_modify:
			if len(row[field]) > 0:
				genotype = {}
				for key, val in row[field].items():
					#if key == 'EG10440_RNA[c]':
						#breakpoint()
					for new_key in rna_id_to_rna_set[key.strip("[c]")]:
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

