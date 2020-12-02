from wholecell.io.tablereader import TableReader
from reconstruction.spreadsheets import JsonReader
import os
import glob
import copy
from collections import defaultdict
import csv
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


def get_functional_genes():
    # generate this file by running runscripts/reflect/model_inspection.py
    with open('runscripts/reflect/functional_genes.tsv') as f:

        tsvreader = csv.reader(f, delimiter='\t')

        process_dict = defaultdict(list)

        lcount = 0
        for line in tsvreader:

            if lcount > 1:
                process_dict[line[3]].append(line[2])

            lcount += 1

    # some of the annotations seem redundant, so let's concatenate them
    process_dict['Transcription'] = process_dict['Transcription'] + process_dict['Transcription Regulation']
    process_dict['Metabolism'] = process_dict['Metabolism'] + process_dict['Metabolism, RNA Decay']
    process_dict['Metabolism'] = process_dict['Metabolism'] + process_dict['Transcription Regulation, Metabolism']
    process_dict['RNA Decay'] = process_dict['RNA Decay'] + process_dict['Metabolism, RNA Decay']
    process_dict['Transcription'] = process_dict['Transcription'] + process_dict['Transcription Regulation, Metabolism']

    return process_dict


def cat_gen_data(gen_path):
    gen_data = {}
    gen_data['monomer_data'] = {}
    gen_data['rna_data'] = {}
    gen_data['monomer_data']['gens'] = {}
    gen_data['rna_data']['gens'] = {}
    for idx, simdir in enumerate(sorted(glob.glob(gen_path + '/generation*'))):

        simOutDir = simdir + '/000000/simOut'

        monomerCounts = TableReader(os.path.join(simOutDir, "MonomerCounts"))
        ids = monomerCounts.readAttribute("monomerIds")
        counts = monomerCounts.readColumn("monomerCounts")
        gen_data['monomer_data']['gens'][idx] = counts


        mRNA_counts_reader = TableReader(os.path.join(simOutDir, 'mRNACounts'))
        mRNA_names = mRNA_counts_reader.readAttribute('mRNA_ids')
        mRNACounts = mRNA_counts_reader.readColumn('mRNA_counts')
        gen_data['rna_data']['gens'][idx] = mRNACounts


    gen_data['rna_data']['names'] = [name[:-3] for name in mRNA_names]
    gen_data['monomer_data']['names'] = [name[:-3] for name in ids]

    return gen_data

def operon_id_to_prot_rna_mapping():

    with open('reconstruction/ecoli/flat/proteins.tsv') as f:

        tsvreader = csv.reader(f, delimiter='\t')

        id_to_prot = {}
        lcount = 0
        for line in tsvreader:

            if lcount > 0:
                id_to_prot[line[5]] = line[0]
            elif line[0] != 'id' or line[5] != "gene_id":
                print("proteins.tsv headers do not match expected values")
                break

            lcount += 1

    with open('reconstruction/ecoli/flat/operon_rnas.tsv') as f:

        tsvreader = csv.reader(f, delimiter='\t')

        id_to_rna = defaultdict(list)
        lcount = 0
        for line in tsvreader:

            if lcount > 1:
                for gene in line[3]:
                    id_to_rna[gene].append(line[0])



            lcount += 1
    return id_to_prot, id_to_rna



def plot_operon_master(operon_data, master_data, genes, fname):
    cmap = [plt.cm.tab20(i) for i in np.linspace(0, 1, 15)]
    min_gens = min([len(master_data['monomer_data']['gens'].keys()) , len(operon_data['monomer_data']['gens'].keys())])
    with PdfPages('out/operon_vs_master_multigen_prots' + fname + '.pdf') as pdf:
        for gene in genes:
            op_t = 0
            m_t = 0
            # import ipdb; ipdb.set_trace()
            operon_idx = operon_data['monomer_data']['names'].index(gene)
            master_idx = master_data['monomer_data']['names'].index(gene)

            fig = plt.figure(dpi=300)
            for g in range(0, min_gens):
                # import ipdb; ipdb.set_trace()
                op_t_end = operon_data['monomer_data']['gens'][g].shape[0]
                m_t_end = master_data['monomer_data']['gens'][g].shape[0]

                plt.plot(range(op_t, op_t + op_t_end), operon_data['monomer_data']['gens'][g][:,operon_idx], color=cmap[g])
                plt.plot(range(m_t, m_t + m_t_end), master_data['monomer_data']['gens'][g][:, master_idx], color=cmap[g], linestyle=':')

                op_t += op_t_end
                m_t += m_t_end

            plt.title(gene)
            pdf.savefig(fig)
            plt.close()


# operon_gen_data = cat_gen_data('/Users/taryn/GoogleDrive/code/wcEcoli/out/allTU/wildtype_000000/000002')
# operon_gen_data = cat_gen_data('/Users/taryn/GoogleDrive/code/wcEcoli/out/no_operon_32g/wildtype_000000/000000')
operon_gen_data = cat_gen_data('/Users/taryn/GoogleDrive/code/wcEcoli/out/no_func_genes/wildtype_000000/000000/')

master_gen_data = cat_gen_data('/Users/taryn/GoogleDrive/code/wcEcoli_master/wcEcoli/out/master1019/wildtype_000000/000000')

func_gene_dict = get_functional_genes()
plot_operon_master(operon_gen_data, master_gen_data, func_gene_dict['Transcription'], '_transcription')
plot_operon_master(operon_gen_data, master_gen_data, func_gene_dict['Translation'], '_tranlsation')
