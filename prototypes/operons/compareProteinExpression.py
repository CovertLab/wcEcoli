from __future__ import absolute_import, division, print_function

import os
from six.moves import cPickle

import numpy as np
from matplotlib import pyplot as plt
from scipy.stats import pearsonr

from wholecell.io.tablereader import TableReader
from wholecell.utils.protein_counts import get_simulated_validation_counts
from wholecell.analysis.analysis_tools import exportFigure
from models.ecoli.analysis import singleAnalysisPlot






def compareProteinExpression(sim1Name, sim1Dir, sim2Name, sim2Dir, plotOutDir):

    monomer_counts = {}
    monomer_names = {}
    for simName, simDir in zip([sim1Name, sim2Name], [sim1Dir, sim2Dir]):
        monomer_counts_reader = TableReader(os.path.join(simDir, "MonomerCounts"))
        monomer_counts[simName] = monomer_counts_reader.readColumn("monomerCounts")
        monomer_names[simName] = monomer_counts_reader.readAttribute("monomerIds")

    if monomer_names[sim1Name] != monomer_names[sim2Name]:
        print("monomer IDs not equivalent")
        # actually do something about this problem

    plt.figure(figsize=(8.5,11))

    plt.scatter(monomer_counts[sim1Name].mean(0), monomer_counts[sim2Name].mean(0), c='w', edgecolor='k', alpha=.7)

    plt.xlabel(sim1Name + ' protein expression')
    plt.ylabel(sim2Name + ' protein expression')
    plt.suptitle('Average Protein Expression Level', fontweight="bold")
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    exportFigure(plt, plotOutDir, sim1name + '_' + sim2name + '_protExpressionScatter')
    plt.close("all")

# ------------------------------------------------------------
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('sim 1 name', type=str, help='name of first sim')
    parser.add_argument('sim1Dir', type=str, help='directory containing sim (str)')
    parser.add_argument('sim 2 name', type=str, help='name of second sim')
    parser.add_argument('sim2Dir', type=str, help='directory containing master sim (str)')
    args = vars(parser.parse_args())
    OutDir = args['sim1Dir'].split('simOut')[0] + 'plotOut/'
    MassPlots(args['sim 1 name'], args['sim1Dir'], args['sim 2 name'],  args['sim2Dir'], OutDir)