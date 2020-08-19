"""
compare proteomics validation plots between two sims
color proteins which are present in the operons denoted in polycistrons_in_the_model.tsv
"""
import argparse
import os
from six.moves import cPickle

import numpy as np
from matplotlib import pyplot as plt
from matplotlib import cm
from scipy.stats import pearsonr

from prototypes.operons.get_operon_rna_list import pc_mc_conversion_dicts

from wholecell.io.tablereader import TableReader
from wholecell.utils.protein_counts import get_simulated_validation_counts

from functools import partial
from reconstruction import spreadsheets



def compareProteomics(sim1Name, sim1Dir, sim2Name, sim2Dir, pc_rnas_path):
    # get polycistron names
    pc_to_mc_dict, mc_to_pc_dict = pc_mc_conversion_dicts(pc_rnas_path)
    cmap = cm.get_cmap('rainbow', len(pc_to_mc_dict.keys()))
    clist = np.linspace(0,1,len(pc_to_mc_dict))

    # rna to protein dictionary
    rna_to_protein_dict = build_rna_to_protein_dict()

    # get path to validation file from sim1
    path = sim1Dir
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

    # load validation data
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

    for p, simName in enumerate([sim1Name, sim2Name]):
        x = np.log10(wisniewski_counts + 1)
        y = np.log10(sim_wisniewski_counts[simName] + 1)
        r, pval = pearsonr(x,y)
        ax[p].scatter(x,y, c='w', edgecolor='k', alpha=.7)
        ax[p].set_title(simName + ' r2 = {r2:.3f}'.format(r2 = r))

        # plot genes in operons with color (probably need to remove when we have more operons!)
        cn=0
        for genes in pc_to_mc_dict.values():
            for gene in genes:
                try:
                    prot = rna_to_protein_dict[gene.strip('[c]')]
                    if prot in wisniewski_ids:
                        idx = np.where(wisniewski_ids == prot)[0][0]
                        x = np.log10(wisniewski_counts[idx] + 1)
                        y = np.log10(sim_wisniewski_counts[simName][idx] + 1)
                        ax[p].scatter(x,y, color=cmap(clist[cn]), edgecolor='k')
                except:
                    import ipdb; ipdb.set_trace()
            cn += 1



    plt.suptitle('Wisniewski Proteomics Data Comparison')
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    OutDir = sim1Dir.split('simOut')[0] + 'plotOut/'
    if not os.path.isdir(OutDir):
        os.mkdir(OutDir)
    plt.savefig(OutDir + sim1Name + '_' + sim2Name + '_wisniewski_comparison.png')


def build_rna_to_protein_dict():
    protein_tsv = 'reconstruction/ecoli/flat/proteins.tsv'
    DIALECT = "excel-tab"
    JsonReader = partial(spreadsheets.JsonReader, dialect=DIALECT)

    tsv_list = []
    with open(protein_tsv) as tsvfile:
        reader = JsonReader(tsvfile)
        fieldnames = reader.fieldnames
        for row in reader:
            tsv_list.append(row)

    rna_to_protein_dict = {}
    for row in tsv_list:
        rna_to_protein_dict[row['rnaId']] = row['id'] + '[' + row['location'][0] + ']'

    return rna_to_protein_dict
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
