import csv
import os
import numpy as np
from collections import defaultdict
from prototypes.operons.get_operon_rna_list import pc_mc_conversion_dicts
import matplotlib.pyplot as plt
import argparse
from wholecell.io.tablereader import TableReader
from functools import partial
from reconstruction import spreadsheets


def protein_rna_expression_changes(sim1Name, sim1Dir, sim2Name, sim2Dir):
    # LOAD RNA AND PROTEIN DNA
    mRNA_names = {}
    mRNA_counts = {}
    monomer_counts = {}
    monomer_names = {}

    for simName, simDir in zip([sim1Name, sim2Name], [sim1Dir, sim2Dir]):
        mRNA_counts_reader = TableReader(os.path.join(simDir, 'mRNACounts'))
        mRNA_names[simName] = mRNA_counts_reader.readAttribute('mRNA_ids')
        mRNA_counts[simName] = mRNA_counts_reader.readColumn('mRNA_counts')

        monomer_counts_reader = TableReader(os.path.join(simDir, "MonomerCounts"))
        monomer_counts[simName] = monomer_counts_reader.readColumn("monomerCounts")
        monomer_names[simName] = monomer_counts_reader.readAttribute("monomerIds")



    m_diff = []
    m_names = []
    m2_counts = []
    m1_counts = []
    for m in monomer_names[sim1Name]:

        try:
            m_names.append(m)
            idx_1 = monomer_names[sim1Name].index(m)
            idx_2 = monomer_names[sim2Name].index(m)
            m2_mean = np.mean(monomer_counts[sim2Name][:, idx_2])
            m1_mean = np.mean(monomer_counts[sim1Name][:, idx_1])
            m2_counts.append(m2_mean)
            m1_counts.append(m1_mean)

            m_diff.append(m1_mean/m2_mean)
        except:
            import ipdb; ipdb.set_trace()
            print(m)


    sorted_idx = np.array(m_diff, dtype=float).argsort()




    OutDir = args['sim2Dir'].split('simOut')[0] + 'plotOut/'

    with open(OutDir + 'protein_diffs.tsv', 'wt') as f:
        tsv_writer = csv.writer(f, delimiter = '\t')
        tsv_writer.writerow(['monomer name', sim1Name, sim2Name, 'fold change'])
        for i in sorted_idx:
            vals = [round(x,2) for x in [m1_counts[i], m2_counts[i], m_diff[i]]]
            row = [m_names[i]] + vals
            tsv_writer.writerow(row)


# -----------------------------------------------------------------------------------------------------------------
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('sim 1 name', type=str, help='name of first sim')
    parser.add_argument('sim1Dir', type=str, help='directory containing sim (str)')
    parser.add_argument('sim 2 name', type=str, help='name of second sim')
    parser.add_argument('sim2Dir', type=str, help='directory containing master sim (str)')
    args = vars(parser.parse_args())
    protein_rna_expression_changes(args['sim 1 name'], args['sim1Dir'], args['sim 2 name'], args['sim2Dir'])
