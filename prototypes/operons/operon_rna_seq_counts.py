import argparse
import csv
from reconstruction import spreadsheets
from functools import partial
import os
import numpy as np

DIALECT = "excel-tab"
JsonReader = partial(spreadsheets.JsonReader, dialect = DIALECT)
JsonWriter = partial(spreadsheets.JsonWriter, dialect = DIALECT)

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

def rna_seq_disparity(polycistron_file_path):

    rna_seq_path =  os.path.join( "reconstruction", "ecoli", "flat", "rna_seq_data","rnaseq_rsem_tpm_mean.tsv")

    rna_seq_data, fieldnames = parse_tsv(rna_seq_path)

    rna_seq_dict = {}
    for row in rna_seq_data:
        rna_seq_dict[row['Gene']] = row['M9 Glucose minus AAs']



    rna_list, fieldnames = parse_tsv(polycistron_file_path)

    operon_rna_seq = np.empty(shape=(0,3), dtype=object)
    for row in rna_list:
        operon = '_'.join(row['transcription_units'])
        gene_list = []
        for gene in row["transcription_units"]:
            gene_list.append(rna_seq_dict[gene])

        # expression_diff = (max(gene_list) - min(gene_list))/np.mean(gene_list)
        gl_nomax = [x for x in gene_list if x != max(gene_list)]
        expression_diff = (max(gene_list) - np.median(gl_nomax))/np.median(gl_nomax)

        new_row = [operon, gene_list, expression_diff]

        operon_rna_seq = np.vstack((operon_rna_seq, new_row))

    operon_rna_seq = operon_rna_seq[np.array(operon_rna_seq[:,2], dtype=float).argsort(),:]


    with open('out/polycistron_rna_seq_diff.tsv', 'wt') as f:
        tsv_writer = csv.writer(f, delimiter='\t')
        tsv_writer.writerow(['operon', 'rna abundances per gene', 'normalized rna count diff'])
        for row in operon_rna_seq:
            tsv_writer.writerow(row)



# -----------------------------------------------------------------------------------------------------------------
if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('polycistron file', type=str, help='polycistornic_mrnas_in_model.tsv')
    args = vars(parser.parse_args())
    rna_seq_disparity(args['polycistron file'])

