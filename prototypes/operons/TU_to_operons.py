from reconstruction.spreadsheets import tsv_reader
from collections import defaultdict


def TU_to_operons(operon_rnas):
    with tsv_reader(operon_rnas) as reader:
        tsv_fieldnames = reader.fieldnames
        tsv_data = list(reader)


    geneid_to_rna_names = defaultdict(list)
    rna_name_to_geneids = defaultdict(list)
    polycistron_list = []

    for rna in tsv_data:
        if rna['id'].count('_') > 1:
            polycistron_list.append(rna['id'])
        rna_name_to_geneids[rna['id']] = rna['gene_set']
        for gene in rna['gene_set']:
            geneid_to_rna_names[gene].append(rna['id'])


    polycistron_list = sorted(polycistron_list, key=lambda x: x.count("_"), reverse=True)
    rna_names_list = sorted(list(rna_name_to_geneids.keys()), key=len, reverse=True)

    operons = defaultdict(list)
    while len(polycistron_list) > 0:
        for pc in polycistron_list:
            for gene in rna_name_to_geneids[pc]:
                for rna in geneid_to_rna_names[gene]:
                    if rna not in operons[pc]:
                        operons[pc].append(rna)
                    if rna in polycistron_list:
                        polycistron_list.remove(rna)


    



