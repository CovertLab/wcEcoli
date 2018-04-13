'''
Script to generate flat file for TU information
Inputs (from transcription_units flat directory:
	geneTUMatrix.pkl - matrix of all genes to the transcription unit they belong
	geneTUMatrix_Rows.pkl - row labels for the matrix (genes)
	geneTUMatrix_Columns.pkl - column labels for the matrix (list of genes in TU)

TODO - better annotation
TODO - explain choices
TODO - explain future considerations
TODO - document functions
'''

import json
import cPickle
import numpy as np
import os
import re

from reconstruction.ecoli.knowledge_base_raw import KnowledgeBaseEcoli


FLAT_DIR = os.path.join('reconstruction', 'ecoli', 'flat')
TU_FLAT_OUTPUT_FILE = os.path.join(FLAT_DIR, 'transcription_units.tsv')
TU_ID_FORMAT = 'TU{:05d}'

TU_DIR = os.path.join(FLAT_DIR, 'transcription_units')
MATRIX_FILE = os.path.join(TU_DIR, 'geneTUMatrix.pkl')
ROWS_FILE = os.path.join(TU_DIR, 'geneTUMatrix_Rows.pkl')
COLS_FILE = os.path.join(TU_DIR, 'geneTUMatrix_Columns.pkl')
JSON_OUTPUT_FILE = os.path.join(TU_DIR, 'tu.json')
MANUAL_DIRECTIONALITY_FILE = os.path.join(TU_DIR, '041018_TU_Directionality.csv')

OTHER_RNA_IDX = 10  # molecular weight index for other RNA
N_MW_TYPES = 11  # number of molecular weight groups that are tracked

def load_tu_matrix(raw_data, rerun=True):
	if rerun:
		# read cPickle data from Mialy
		with open(MATRIX_FILE) as f:
			mat = np.array(cPickle.load(f))
		with open(ROWS_FILE) as f:
			rows = cPickle.load(f)
		with open(COLS_FILE) as f:
			cols = cPickle.load(f)

		# map ecocyc ID to wcm ID for row and col annotation
		ec_to_wc_gene_map =  {gene['id']: gene['rnaId'] for gene in raw_data.genes}
		wcm_rows = [ec_to_wc_gene_map[gene] for gene in rows]
		wcm_cols = [sorted(set([ec_to_wc_gene_map[gene] for gene in col if gene in ec_to_wc_gene_map])) for col in cols]

		# remove phantom genes in data but not included in model
		n_genes_in_tu = np.sum(mat, axis=0)
		wcm_cols = np.array(wcm_cols)[n_genes_in_tu > 0].tolist()
		mat = mat[:, n_genes_in_tu > 0]

		assert np.all(np.sum(mat, axis=0) == np.array([len(col) for col in wcm_cols]))

		# represent matrix more compactly
		mati, matj = np.where(mat)

		# create new entries mapping genes that do not appear in data to their own operon
		for gene_idx in np.where(np.sum(mat, axis=1) == 0)[0]:
			mati = np.hstack((mati, gene_idx))
			matj = np.hstack((matj, len(wcm_cols)))
			wcm_cols.append([wcm_rows[gene_idx]])

		# sort rows by name to match order in sim_data for synth prob
		sorted_indices = {j: i for i, j in enumerate(np.argsort(wcm_rows))}
		mati_sorted = np.array([sorted_indices[idx] for idx in mati])

		# reconstruct matrix with new columns
		dims = (len(wcm_rows), len(wcm_cols))
		expanded_mat = np.zeros(dims)
		expanded_mat[mati_sorted, matj] = 1

		# create dict with all information
		tu = {}
		tu['mati'] = mati_sorted.tolist()
		tu['matj'] = matj.tolist()
		tu['rows'] = sorted(wcm_rows)
		tu['cols'] = wcm_cols
		tu['dims'] = dims
	else:
		with open(JSON_OUTPUT_FILE) as f:
			tu = json.load(f)

	return tu

def save_tu_matrix(tu):
	# save json object
	with open(JSON_OUTPUT_FILE, 'w') as f:
		json.dump(tu, f)

def create_tu_matrix(tu):
	# reconstruct matrix with new columns
	mat = np.zeros(tu['dims'])
	mat[tu['mati'], tu['matj']] = 1

	return mat


# Load raw data
raw_data = KnowledgeBaseEcoli()

# Load genome data
genome = raw_data.genome_sequence
reverse_complement = genome.reverse_complement()

# Load masses of interest
ntps = ['ATP', 'CTP', 'GTP', 'UTP']
small_molecule_mws = {mol['id']: mol['mw7.2'] for mol in raw_data.metabolites}
end_weight = small_molecule_mws["PPI"]
ntp_weights = np.array([small_molecule_mws[ntp] - end_weight for ntp in ntps])
letter_to_index = {nt[0]: i for i, nt in enumerate(ntps)}
rna_mws = {rna['id']: rna['mw'] for rna in raw_data.rnas}

# Create transcription unit matrix (genes x tu)
tu = load_tu_matrix(raw_data)
save_tu_matrix(tu)
tu_mat = create_tu_matrix(tu)
rows = tu['rows']
cols = tu['cols']

# Load RNA data from raw_data
gene_dict = {gene['rnaId']: gene for gene in raw_data.genes}
rna_seqs = {rna['id']: rna['seq'] for rna in raw_data.rnas}
rna_in_complexes = np.hstack([
	np.array([stoich['molecule']
	for stoich in rxn['stoichiometry'] if stoich['type'] == 'rna'])
	for rxn in raw_data.complexationReactions
	])

# Pull out information about genes in TUs
single_genes = set([tu[0] for tu in cols if len(tu) == 1])
n_tu_per_gene = {gene: count for gene, count in zip(rows, np.sum(tu_mat, axis=1))}

# Write to csv for manual curation of some directionality
import csv
writer = csv.writer(open('opposite_direction_genes.tsv', 'w'), delimiter='\t')

# Create a TU for each column in the matrix
count = 0
with open(TU_FLAT_OUTPUT_FILE, 'w') as f:
	f.write('"length"\t"id"\t"seq"\t"coordinate"\t"direction"\t"rnas"\t"mw"\t"processed"\n')
	for mat_col, col_name in zip(tu_mat.T, cols):
		starts = []
		ends = []
		dirs = []
		ids = []
		types = []
		rna_ids = []
		seqs = []

		# Get gene info for each transcription unit
		for gene_idx in np.where(mat_col)[0]:
			gene_id = rows[gene_idx]

			assert rows[gene_idx] in col_name

			dir = str(gene_dict[gene_id]['direction'])
			start = int(gene_dict[gene_id]['coordinate'])
			if dir == '-':
				start += 1
			length = int(gene_dict[gene_id]['length'])
			if dir == '+':
				end = start + length
			else:
				end = start - length

			dirs += [dir]
			starts += [start]
			ends += [end]
			ids += [str(gene_dict[gene_id]['id'])]
			rna_ids += [str(gene_dict[gene_id]['rnaId'])]
			types += [str(gene_dict[gene_id]['type'])]
			seqs += [str(gene_dict[gene_id]['seq'])]

		sorted_idx = np.argsort(starts)

		dirs = np.array(dirs)[sorted_idx]
		starts = np.array(starts)[sorted_idx]
		ends = np.array(ends)[sorted_idx]
		ids = np.array(ids)[sorted_idx]
		rna_ids = np.array(rna_ids)[sorted_idx]
		types = np.array(types)[sorted_idx]

		# Determine directionality of TU
		if len(dirs) > 1:
			ignored_misc = types == 'miscRNA'

			count_fwd = np.sum(dirs[~ignored_misc] == '+')
			count_rev = np.sum(dirs[~ignored_misc] == '-')

			skip = False
			if count_fwd == count_rev:
				# Ignore TUs that can represented as single gene TUs
				if len(dirs) == 2:
					if ((rna_ids[0] in single_genes or n_tu_per_gene[rna_ids[0]] == 1
							or (n_tu_per_gene[rna_ids[0]] == 2 and rna_ids[0] in single_genes))
							and (rna_ids[1] in single_genes or n_tu_per_gene[rna_ids[1]] == 1
							or (n_tu_per_gene[rna_ids[1]] == 2 and rna_ids[1] in single_genes))):
						skip = True


				# TODO - remove skipped TUs
				# TODO - handle manual curation from file
				# TODO - verify skipped TUs don't leave a gene without a TU
				if not skip:
					writer.writerow([ids, dirs, types, starts])
					# print col_name

				# flagged = False
				# internal = False
				# start_dir = ''
				# for d in sorted_dir:
				# 	if start_dir == '':
				# 		start_dir = d
				# 	elif d != start_dir:
				# 		flagged = True
				# 	elif flagged and d == start_dir:
				# 		internal = True
				#
				# if flagged:
			elif count_fwd > count_rev:
				tu_dir = '+'
			else:
				tu_dir = '-'
		else:
			tu_dir = dirs[0]

		# Determine mask for which entries match TU direction
		dir_mask = dirs == tu_dir

		# Determine which proteins can be translated from TU based on directionality
		tu_genes = rna_ids[dir_mask]

		# Determine sequence of TU
		if tu_dir == '+':
			tu_start = np.min(np.hstack((starts[dir_mask], ends[~dir_mask])))
			tu_end = np.max(np.hstack((ends[dir_mask], starts[~dir_mask])))
			tu_seq = re.sub('T', 'U', str(genome[tu_start:tu_end]))
		else:
			tu_start = np.max(np.hstack((starts[dir_mask], ends[~dir_mask])))
			tu_end = np.min(np.hstack((ends[dir_mask], starts[~dir_mask])))
			tu_seq = re.sub('T', 'U', str(reverse_complement[-tu_start:-tu_end]))

		tu_length = len(tu_seq)

		# Determine molecule weight of TU
		## Total MW of TU sequence
		ntp_counts = np.array([
			tu_seq.count('A'), tu_seq.count('C'),
			tu_seq.count('G'), tu_seq.count('U')
			])
		tu_mw = np.zeros(N_MW_TYPES)
		tu_mw[OTHER_RNA_IDX] = np.dot(ntp_counts, ntp_weights) + end_weight

		## MW for each gene and its type that is part of TU
		component_types = np.zeros((len(tu_seq), N_MW_TYPES))
		for rna in tu_genes:
			rna_mw = np.array(rna_mws[rna])

			seq = rna_seqs[rna]
			start = tu_seq.index(seq)
			end = start + len(seq)
			component_types[start:end, np.where(rna_mws[rna])[0]] += 1

		## Account for overlapping genes in a TU and split MW evenly between
		## different types for shared nucleotides
		sequence_mw = np.array([ntp_weights[letter_to_index[ntp]] for ntp in tu_seq])
		sequence_mw[-1] += end_weight
		shared_idx = (np.sum(component_types, axis=1) > 1)
		component_types[shared_idx, :] = (component_types[shared_idx, :]
			/ np.sum(component_types[shared_idx, :], axis=1).reshape(-1, 1)
			)
		component_mw = np.dot(sequence_mw, component_types)
		tu_mw += component_mw
		tu_mw[OTHER_RNA_IDX] -= np.sum(component_mw)

		## Adjust for precision errors
		tu_mw[np.where(np.abs(tu_mw) < 1e-8)[0]] = 0

		# Determine if transcription unit needs to be processed to individual RNA
		# Ribosomal RNA, tRNA and RNA in complexes will need to be processed
		processed = 'false'
		if ('rRNA' in types or 'tRNA' in types
				or np.any([rna in rna_in_complexes for rna in tu_genes])):
			processed = 'true'

		# Set TU ID
		tu_id = TU_ID_FORMAT.format(count)
		count += 1

		# Write info for TU
		f.write('{}\t"{}"\t"{}"\t{}\t"{}"\t{}\t{}\t{}\n'.format(
			tu_length, tu_id, tu_seq, tu_start, tu_dir,
			json.dumps(tu_genes.tolist()), json.dumps(tu_mw.tolist()), processed
		))

import ipdb; ipdb.set_trace()
