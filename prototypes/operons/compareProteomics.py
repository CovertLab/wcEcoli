"""
compare proteomics validation plots between two sims
color proteins which are present in the operons denoted in polycistrons_in_the_model.tsv
"""
import argparse
import os
from six.moves import cPickle

import numpy as np
from matplotlib import pyplot as plt
from scipy.stats import pearsonr

from wholecell.io.tablereader import TableReader
from wholecell.utils.protein_counts import get_simulated_validation_counts


def compareProteomics(sim1Name, sim1Dir, sim2Name, sim2Dir, pc_rnas_path):

    # get path to validation file from sim1
    path = sim2Dir
    if 'out' not in path:
        print('path structure not recognized, validation data file not found')
        print('expected path structure: .../out/...')
        return

    while 1:
        path, folder = os.path.split(path)
        if folder == 'out':
            validationDataFile = os.path.join(path, folder, folderp1, 'kb', 'validationData.cPickle')
            print(validationDataFile)
            break
        else:
            folderp1 = folder


    validation_data = cPickle.load(open(validationDataFile, "rb"))
    wisniewski_ids = validation_data.protein.wisniewski2014Data["monomerId"]
    wisniewski_counts = validation_data.protein.wisniewski2014Data["avgCounts"]

    # get monomer counts from sims
    sim_wisniewski_counts = {}
    monomer_counts = {}
    monomer_names = {}

    for simName, simDir in zip([sim1Name, sim2Name], [sim1Dir, sim2Dir]):

        monomer_counts_reader = TableReader(os.path.join(simDir, "MonomerCounts"))
        monomer_counts[simName] = monomer_counts_reader.readColumn("monomerCounts")
        monomer_names[simName] = monomer_counts_reader.readAttribute("monomerIds")

        sim_wisniewski_counts[simName] = get_simulated_validation_counts(
            monomer_counts[simName], wisniewski_ids, monomer_names[simName])

    # Plot against wisniewski data for both sims
    fig, ax = plt.subplots(2, sharey=True, figsize=(8.5, 11))

    for p, simName in enumerate([sim2Name, sim2Name]):
        ax[p].scatter(
            np.log10(wisniewski_counts + 1),
            np.log10(sim_wisniewski_counts[simName] + 1),
            c='w', edgecolor='k', alpha=.7)
        ax.set_title(simName)

    plt.suptitle('Wisniewski Proteomics Data Comparison')
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    OutDir = sim1Dir.split('simOut')[0] + 'plotOut/'
    print(OutDir)
    plt.savefig(OutDir + sim1Name + '_' + sim2Name + '_wisniewski_comparison.png')

# -----------------------------------------------------------------------------------------------------------------
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('sim 1 name', type=str, help='name of first sim')
    parser.add_argument('sim1Dir', type=str, help='directory containing sim (str)')
    parser.add_argument('sim 2 name', type=str, help='name of second sim')
    parser.add_argument('sim2Dir', type=str, help='directory containing master sim (str)')
    parser.add_argument('polycistron file', type=str, help='polycistornic_mrnas_in_model.tsv for sim 1')
    args = vars(parser.parse_args())
    compareProteomics(args['sim 1 name'], args['sim1Dir'], args['sim 2 name'], args['sim2Dir'], args['polycistron file'])
