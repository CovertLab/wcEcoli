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

            listener_data = TableReader(os.path.join(simOutDir, "GrowthLimits"))
            col_names = listener_data.columnNames()
            if 'Icon\r' in col_names:
                col_names.remove('Icon\r')
            for attribute in col_names:
                gen_data[idx][attribute] = listener_data.readColumn(attribute)

        except:
            import ipdb; ipdb.set_trace()

    return gen_data


def plot_operon_master(operon_data, master_data):
    cmap = [plt.cm.tab20(i) for i in np.linspace(0, 1, 15)]
    min_gens = min([len(master_data.keys()) , len(operon_data.keys())])
    plot_list = ['aaAllocated',
                 # 'simulationStep',
                 'ntpAllocated',
                 'aaPoolSize',
                 'aaRequestSize',
                 'spot_deg',
                 'net_charged',
                 'rela_syn',
                 'activeRibosomeAllocated',
                 'aasUsed',
                 'ntpRequestSize',
                 # 'time',
                 'ntpPoolSize',
                 'ntpUsed',
                 'spot_syn',
                 'fraction_trna_charged']

    with PdfPages('out/rend_seq_ribosome_only/operon_vs_master_growth_limits_multigen.pdf') as pdf:
        for p in plot_list:
            op_t = 0
            m_t = 0

            print(p)
            print(operon_data[0][p].shape)

            fig = plt.figure(dpi=150)
            for g in range(0, min_gens):
                # import ipdb; ipdb.set_trace()
                op_t_end = operon_data[g][p].shape[0]
                m_t_end = master_data[g][p].shape[0]
                op_t_vec = range(op_t, op_t + op_t_end)
                m_t_vec = range(m_t, m_t + m_t_end)

                if len(operon_data[g][p].shape) > 1:

                    for col in range(0,operon_data[g][p].shape[1]):
                        plt.plot(op_t_vec[0::25], operon_data[g][p][0::25,col], color=cmap[g])
                    for col in range(0, master_data[g][p].shape[1]):
                        plt.plot(m_t_vec[0::25], master_data[g][p][0::25,col], color='k', linestyle=':', alpha=0.5)

                else:
                    plt.plot(op_t_vec[0::25], operon_data[g][p][0::25], color=cmap[g])
                    plt.plot(m_t_vec[0::25], master_data[g][p][0::25], color='k', linestyle=':', alpha=0.5)


                op_t += op_t_end
                m_t += m_t_end

            plt.title(p)
            pdf.savefig(fig)
            plt.close()



operon_gen_data = cat_gen_data('/Users/mialydefelice/Documents/code_repositories/wcEcoli_2/wcEcoli/out/rend_seq_ribosome_only/wildtype_000000/000000')
master_gen_data = cat_gen_data('/Users/mialydefelice/Documents/code_repositories/wcEcoli_2/wcEcoli/out/master_troubleshooting/wildtype_000000/000000')
plot_operon_master(operon_gen_data, master_gen_data)