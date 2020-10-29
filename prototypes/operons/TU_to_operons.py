from reconstruction.spreadsheets import tsv_reader
from collections import defaultdict
import copy


def TU_to_operons(operon_rnas):
    with tsv_reader(operon_rnas) as reader:
        tsv_data = list(reader)

    # make dictionaries to map rnaid <> geneid
    # make list of polycistronic rnas
    geneid_to_rna_names = defaultdict(list)
    rna_name_to_geneids = defaultdict(list)
    polycistron_list = []

    for rna in tsv_data:
        if rna['id'].count('_') > 1:
            polycistron_list.append(rna['id'])
        rna_name_to_geneids[rna['id']] = rna['gene_set']
        for gene in rna['gene_set']:
            geneid_to_rna_names[gene].append(rna['id'])

    # sort by operon length
    polycistron_list= sorted(polycistron_list, key=lambda x: x.count("_"), reverse=True)

    # gather all rnas that cover a set of genes
    operons = defaultdict(list)
    genes_in_operons = []

    while len(polycistron_list) > 0:
        for pc in polycistron_list:
            print(pc)
            for gene in rna_name_to_geneids[pc]:
                
                for rna in geneid_to_rna_names[gene]:
                    if rna not in operons[pc]:
                        operons[pc].append(rna)
                    if rna in polycistron_list:
                        polycistron_list.remove(rna)

    # get all genes contained in an operon
    operon_to_genes = {}
    for operon_rna in operons:
        gene_list = []
        for rna in operons[operon_rna]:
            gene_list = gene_list + rna_name_to_geneids[rna]

        operon_to_genes[operon_rna] = set(gene_list)









    return operons, operon_to_genes


def count_disparity(operons, operon_to_genes):
    rna_seq_file = 'reconstruction/ecoli/flat/rna_seq_data/rnaseq_rsem_tpm_mean.tsv'
    with tsv_reader(rna_seq_file) as reader:
        tsv_data = list(reader)

    gene_to_count = {}
    for row in tsv_data:
        gene_to_count[row['Gene']] = row['M9 Glucose minus AAs']

    operon_count_disparity = {}

    for operon_rna in operons.keys():
        operon_count_disparity[operon_rna] = {x : gene_to_count[x] for x in operon_to_genes[operon_rna]}
        TUs = sorted(operons[operon_rna], key=lambda x: x.count('_'), reverse=True)

        for rna in TUs:
            genes_in_rna = rna.split('_')[0:-1]
            min_count = min([gene_to_count[x] for x in genes_in_rna])

            for gene in genes_in_rna:
                try:
                    operon_count_disparity[operon_rna][gene] -= min_count
                except:
                    import ipdb; ipdb.set_trace()

    import ipdb; ipdb.set_trace()



# -----------------------------------------------------------------------------------------------------------------
if __name__ == "__main__":
    operons, operon_to_genes = TU_to_operons('reconstruction/ecoli/flat/operon_rnas.tsv')
    count_disparity(operons, operon_to_genes)














