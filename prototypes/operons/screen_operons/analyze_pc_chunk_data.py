
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



def cat_gen_data(gen_path, listener):
    gen_data = {}
    for idx, simdir in enumerate(sorted(glob.glob(gen_path + '/generation*'))):

        gen_data[idx] = {}

        simOutDir = simdir + '/000000/simOut'
        try:
            listener_data = TableReader(os.path.join(simOutDir, listener))
            col_names = listener_data.columnNames()
            if 'Icon\r' in col_names:
                col_names.remove('Icon\r')
            for attribute in col_names:
                gen_data[idx][attribute] = listener_data.readColumn(attribute)
        except:
            import ipdb; ipdb.set_trace()
    return gen_data



#---------------------------------------------------------------------------------
folder_list = sorted(glob.glob('out/screen_50_tus_*'))
listener_column_dict = {'Mass': ['cellMass', 'dnaMass', 'mRnaMass', 'proteinMass', 'rRnaMass'],
                        'RibosomeData': ['actualElongations', 'didInitialize', 'effectiveElongationRate', 'total_rna_init'],
                        'RnapData': ['didTerminate', 'didInitialize', 'actualElongations'],
                        'RnaDegradationListener': ['FractionActiveEndoRNases', 'fragmentBasesDigested', 'nucleotidesFromDegradation']}

# get all the listener data for all sims and gens
listener_data = {}
for idx, folder in enumerate(folder_list):
    sim_name = folder.split('/')[-1]

    data_path = os.path.join(folder, "wildtype_000000/000000")
    listener_data[sim_name] = {}

    for listener in listener_column_dict.keys():
        listener_data[sim_name][listener] = cat_gen_data(data_path, listener)

master_data = {}
master_data_path = '/Users/taryn/GoogleDrive/code/wcEcoli_master/wcEcoli/out/master1019/wildtype_000000/000000'
for listener in listener_column_dict.keys():
    master_data[listener] = cat_gen_data(master_data_path, listener)

# make a plot for each listener
cmap = [plt.cm.tab20(i) for i in np.linspace(0, 1, 15)]
line_list = [None] * len(listener_data.keys())
with PdfPages('out/screen_50_tus.pdf') as pdf:
    for listener in listener_column_dict.keys():
        for column in listener_column_dict[listener]:
            fig = plt.figure(dpi=150)

            for c, sim_name in enumerate(listener_data.keys()):
                label = sim_name
                for gen in listener_data[sim_name][listener].keys():
                    t_vec = listener_data[sim_name][listener][gen]['time']
                    if gen > 0:
                        label = '_nolegend_'

                    plt.plot(t_vec[0::25], listener_data[sim_name][listener][gen][column][0::25], color=cmap[c], label=label)

            label='master'
            for g in range(0,gen+1):
                if g > 0:
                    label = '_nolegend_'
                t_vec = master_data[listener][g]['time']

                plt.plot(t_vec[0::25], master_data[listener][g][column][0::25], color='k', label=label, alpha=0.6, linewidth=2)

            plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.05), ncol=int(np.ceil(len(listener_data.keys())/2)), fontsize=5)
            print(listener + ': ' + column)
            plt.title(listener + ': ' + column)
            pdf.savefig(fig)
            plt.close()

