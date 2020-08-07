import csv
import os
from collections import defaultdict
from prototypes.operons.get_operon_rna_list import pc_mc_conversion_dicts

def inspect_functional_genes(pc_rnas_path):
    # check if functional_genes.tsv exists
    if not os.path.isfile('runscripts/reflect/functional_genes.tsv'):
        print("functional_genes.tsv does not exist. create this file using runscripts/reflect/model_inspection")
        return


    with open('runscripts/reflect/functional_genes.tsv') as f:

        tsvreader = csv.reader(f, delimiter='\t')

        process_dict = defaultdict(list)

        lcount = 0
        for line in tsvreader:

            if lcount > 1:
                process_dict[line[3]].append(line[1])

            lcount += 1

    # some of the annotations seem redundant, so let's concatenate them
    process_dict['Transcription'].append(process_dict['Transcription Regulation'])
    process_dict['Metabolism'].append(process_dict['Metabolism, RNA Decay'])
    process_dict['Metabolism'].append(process_dict['Transcription Regulation, Metabolism'])
    process_dict['RNA Decay'].append(process_dict['Metabolism, RNA Decay'])
    process_dict['Transcription'].append(process_dict['Transcription Regulation, Metabolism'])


    # ------------------------------------------------------------------------------------------------------------
    # PLOTS
    # plot the RNA and the protein levels for all the functional genes

    keys_to_plot = ['Metabolism', 'RNA Decay', 'Transcription', ]

    # RNA PLOT
    pc_to_mc_dict, mc_to_pc_dict = pc_mc_conversion_dicts(pc_rnas_path)

    for key in keys_to_plot:

