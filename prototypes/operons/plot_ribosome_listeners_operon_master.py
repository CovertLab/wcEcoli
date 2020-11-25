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


def cat_gen_data(gen_path):
    gen_data = {}
    for idx, simdir in enumerate(sorted(glob.glob(gen_path + '/generation*'))):

        try:
            gen_data[idx] = {}

            simOutDir = simdir + '/000000/simOut'

            ribosome_data = TableReader(os.path.join(simOutDir, "RibosomeData"))
            col_names = ribosome_data.columnNames()
            if 'Icon\r' in col_names:
                col_names.remove('Icon\r')
            for attribute in col_names:
                gen_data[idx][attribute] = ribosome_data.readColumn(attribute)

        except:
            import ipdb; ipdb.set_trace()

    return gen_data


def plot_operon_master(operon_data, master_data):
    cmap = [plt.cm.tab20(i) for i in np.linspace(0, 1, 20)]
    min_gens = min([len(master_data.keys()) , len(operon_data.keys())])
    plot_list = ['actualElongations',
                 'didTerminate',
                 'probTranslationPerTranscript',
                 'rrn23S_init_prob',
                 'terminationLoss',
                 'rrn5S_produced',
                 'didInitialize',
                 'n_ribosomes_on_partial_mRNA_per_transcript',
                 'effectiveElongationRate',
                 'n_ribosomes_per_transcript',
                 'rrn16S_produced',
                 'total_rna_init',
                 'translationSupply',
                 'aaCounts',
                 # 'aaCountInSequence',
                 'rrn5S_init_prob',
                 'processElongationRate',
                 'numTrpATerminated',
                 'rrn23S_produced']

    with PdfPages('out/operon_vs_master_ribosome_data_multigen.pdf') as pdf:
        for p in plot_list:
            op_t = 0
            m_t = 0


            print(p)
            print(operon_data[0][p].shape)
            fig = plt.figure(dpi=75)
            for g in range(0, min_gens):
                # import ipdb; ipdb.set_trace()
                op_t_end = operon_data[g][p].shape[0]
                m_t_end = master_data[g][p].shape[0]

                plt.plot(range(op_t, op_t + op_t_end), operon_data[g][p], color=cmap[g])
                plt.plot(range(m_t, m_t + m_t_end), master_data[g][p], color='k', linestyle=':', alpha=0.5)

                op_t += op_t_end
                m_t += m_t_end

            plt.title(p)
            pdf.savefig(fig)
            plt.close()



operon_gen_data = cat_gen_data('/Users/taryn/GoogleDrive/code/wcEcoli/out/allTU_1124/wildtype_000000/000000')
master_gen_data = cat_gen_data('/Users/taryn/GoogleDrive/code/wcEcoli_master/wcEcoli/out/master1019/wildtype_000000/000000')
plot_operon_master(operon_gen_data, master_gen_data)