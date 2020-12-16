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
import argparse



def cat_gen_data(gen_path, listeners):
    gen_data = {}
    for listener in listeners:
        gen_data[listener] = {}

        for idx, simdir in enumerate(sorted(glob.glob(gen_path + '/generation*'))):

            gen_data[listener][idx] = {}

            simOutDir = simdir + '/000000/simOut'

            listener_data = TableReader(os.path.join(simOutDir, listener))
            col_names = listener_data.columnNames()
            if 'Icon\r' in col_names:
                col_names.remove('Icon\r')
            for attribute in col_names:
                gen_data[listener][idx][attribute] = listener_data.readColumn(attribute)

    return gen_data


def plot_listeners(sim1_data, sim2_data, plot_dict, out_file):
    firstkey = list(plot_dict.keys())[0]
    min_gens = min([len(sim2_data[firstkey].keys()) , len(sim1_data[firstkey].keys())])
    max_gens = max([len(sim2_data[firstkey].keys()) , len(sim1_data[firstkey].keys())])
    cmap = [plt.cm.tab20(i) for i in np.linspace(0, 1, max_gens)]

    with PdfPages(out_file) as pdf:
        for k in plot_dict:
            for p in plot_dict[k]:

                fig = plt.figure(dpi=75)
                for g in range(0, min_gens):

                    plt.plot(sim1_data[k][g]['time'][0::25], sim1_data[k][g][p][0::25], color=cmap[g])
                    plt.plot(sim2_data[k][g]['time'][0::25], sim2_data[k][g][p][0::25], color='k', linestyle=':', alpha=0.5)

                plt.title(k + ': ' + p)
                pdf.savefig(fig)
                plt.close()

def run_listeners(sim1, sim2, outdir, plotoptions):
    listener_dict = {}
    listener_dict['minimal'] ={'Mass': ['cellMass', 'dnaMass', 'mRnaMass', 'proteinMass', 'rRnaMass'],
                        'RibosomeData': ['actualElongations', 'didInitialize', 'effectiveElongationRate', 'total_rna_init'],
                        'RnapData': ['didTerminate', 'didInitialize', 'actualElongations'],
                        'RnaDegradationListener': ['FractionActiveEndoRNases', 'fragmentBasesDigested', 'nucleotidesFromDegradation']}

    listener_dict['all'] = {'RnapData' : [
                # 'active_rnap_coordinates',
                #  'active_rnap_domain_indexes',
                 'didTerminate',
                 # 'rnaInitEvent',
                 'n_removed_ribosomes',
                 # 'active_rnap_n_bound_ribosomes',
                 # 'headon_collision_coordinates',
                 # 'n_codirectional_collisions',
                 'terminationLoss',
                 'didInitialize',
                 # 'codirectional_collision_coordinates',
                 'actualElongations',
                 # 'simulationStep',
                 # 'time',
                 'n_total_collisions',
                 # 'active_rnap_unique_indexes',
                 'n_headon_collisions'],

                'GrowthLimits': ['aaAllocated',
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
                 'fraction_trna_charged'],

                'RnaSynthProb':[
                # 'gene_copy_number',
                #  'bound_TF_coordinates',
                 # 'simulationStep',
                 # 'rnaSynthProb',
                 # 'bound_TF_domains',
                 'nPromoterBound',
                 'nActualBound',
                 # 'bound_TF_indexes',
                 # 'n_bound_TF_per_TU',
                 # 'time',
                 'pPromoterBound'],

                'RibosomeData': ['actualElongations',
                 'didTerminate',
                 # 'probTranslationPerTranscript',
                 'rrn23S_init_prob',
                 'terminationLoss',
                 'rrn5S_produced',
                 'didInitialize',
                 # 'n_ribosomes_on_partial_mRNA_per_transcript',
                 'effectiveElongationRate',
                 # 'n_ribosomes_per_transcript',
                 'rrn16S_produced',
                 'total_rna_init',
                 'translationSupply',
                 'aaCounts',
                 # 'aaCountInSequence',
                 'rrn5S_init_prob',
                 'processElongationRate',
                 'numTrpATerminated',
                 'rrn23S_produced']}

    sim1_data = cat_gen_data(sim1, [*listener_dict[plotoptions]])
    sim2_data = cat_gen_data(sim2, [*listener_dict[plotoptions]])

    out_file = os.path.join(outdir, 'listener_plots.pdf')
    plot_listeners(sim1_data, sim2_data, listener_dict[plotoptions], out_file)


def make_sim_paths(dirname, seednum):
    seed = '000000'
    str_seed = str(seednum)
    zero_pad = '0'*(len(seed)-len(str_seed))
    seed = zero_pad + str_seed

    seed_path = 'wildtype_000000/' + seed +'/'

    if dirname.count('/') < 1: # assume it is a sim name not a directory and is int out/name
        simdir = os.path.join('out', dirname, seed_path)
    else:
        simdir = os.path.join(dirname, seed_path)

    return simdir

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('sim1Dir', type=str, help='directory containing sim (str)')
    parser.add_argument('sim2Dir', type=str, help='directory containing master sim (str)')
    parser.add_argument('-s1', default=0, help='seed number for first sim, default is 0')
    parser.add_argument('-s2', default=0, help='seed number for second sim, default is 0')
    parser.add_argument('-plots', default='minimal', help='minimal or all, default=minimal')


    args = vars(parser.parse_args())

    sim1dir = make_sim_paths(args['sim1Dir'], args['s1'])
    sim2dir = make_sim_paths(args['sim2Dir'], args['s2'])


    # save plots in directory of first sim
    outdir = os.path.join(sim1dir, 'plotOut')
    if not os.path.isdir(outdir):
        os.mkdir(outdir)

    run_listeners(sim1dir, sim2dir, outdir, args['plots'])