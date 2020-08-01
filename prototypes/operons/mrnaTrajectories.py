from __future__ import absolute_import, division, print_function

import os
import argparse
from matplotlib import pyplot as plt
from scipy.stats import pearsonr
import numpy as np

from wholecell.io.tablereader import TableReader
from wholecell.analysis.analysis_tools import exportFigure
from prototypes.operons.get_operon_rna_list import get_operon_rna_list


def mrnaTrajectories(sim1Name, sim1Dir, sim2Name, sim2Dir, polycistron_file, plotOutDir):

    # GET CANDIDATE mRNAs
    # TAKE 25 mRNAs WITH GREATEST EXPRESSION DISPARITY BETWEEN POLY AND MONO CISTRON SIMS
    mRNA_names = {}
    mRNA_counts = {}
    for simName, simDir in zip([sim1Name, sim2Name], [sim1Dir, sim2Dir]):
        mRNA_counts_reader = TableReader(os.path.join(simDir, 'mRNACounts'))
        mRNA_names[simName] = mRNA_counts_reader.readAttribute('mRNA_ids')
        mRNA_counts[simName] = mRNA_counts_reader.readColumn('mRNA_counts')

    pc_dict = get_operon_rna_list(polycistron_file)

    # filter out subgen genes based on sim2 monocistron expression levels
    gene_list = []
    expression_disparity = []
    for polycis in pc_dict:
        polycis_idx = mRNA_names[sim1Name].index(polycis)
        if type(polycis_idx) == tuple:
            import ipdb; ipdb.set_trace()
        t_expressed = []
        monocis_expression = np.array([]).reshape(0,mRNA_counts[simName].shape[0])

        for monocis in pc_dict[polycis]:
            monocis_idx = mRNA_names[sim2Name].index(monocis)
            if type(monocis_idx) == tuple:
                import ipdb; ipdb.set_trace()
            t_expressed.append(sum(mRNA_counts[sim2Name][:, monocis_idx] > 0) / mRNA_counts[sim2Name].shape[0])
            monocis_expression = np.vstack((monocis_expression, mRNA_counts[sim2Name][:, monocis_idx]))

        avg_expression_diff = np.mean(mRNA_counts[sim1Name][:, polycis_idx]) - np.mean(np.mean(monocis_expression, axis=0))


        if np.mean(t_expressed) > 0.3:
            gene_list.append(polycis)
            expression_disparity.append(avg_expression_diff)

    sorted_gene_list = [x for _, x in zip(expression_disparity, gene_list)]

    # -------------------------------------------------------------
    # DO PLOT
    num_plots = min(25, len(sorted_gene_list))

    r = np.ceil(np.sqrt(num_plots))
    c = np.ceil(num_plots/r)



    sim1_rnas = sorted_gene_list
    plt.figure(figsize=(8.5, 11))
    # import ipdb; ipdb.set_trace()
    for p in range(0,num_plots):

        plt.subplot(r,c,p+1)

        for rna in pc_dict[sim1_rnas[p]]:
            rna_idx = mRNA_names[sim2Name].index(rna)
            plt.plot(mRNA_counts[sim2Name][:,rna_idx], color='silver')


        plt.plot(mRNA_counts[sim1Name][:, mRNA_names[sim1Name].index(sim1_rnas[p])], color='crimson', linewidth=2)
        plt.title(sim1_rnas[p], fontsize=3)


    # save plot
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    exportFigure(plt, plotOutDir, sim1Name + '_' + sim2Name + '_mRNA_Trajectories')
    plt.close("all")

# ------------------------------------------------------------
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('sim 1 name', type=str, help='name of first sim')
    parser.add_argument('sim1Dir', type=str, help='directory containing sim (str)')
    parser.add_argument('sim 2 name', type=str, help='name of second sim')
    parser.add_argument('sim2Dir', type=str, help='directory containing master sim (str)')
    parser.add_argument('polycistron file', type=str, help='polycistornic_mrnas_in_model.tsv for sim 1')
    args = vars(parser.parse_args())
    OutDir = args['sim1Dir'].split('simOut')[0] + 'plotOut/'
    mrnaTrajectories(args['sim 1 name'], args['sim1Dir'], args['sim 2 name'],  args['sim2Dir'], args['polycistron file'], OutDir)