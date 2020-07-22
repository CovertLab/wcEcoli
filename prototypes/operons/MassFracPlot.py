from __future__ import absolute_import, division, print_function

import os

import numpy as np
from matplotlib import pyplot as plt

from wholecell.io.tablereader import TableReader
from wholecell.analysis.analysis_tools import exportFigure
from models.ecoli.analysis import singleAnalysisPlot
from six.moves import zip
import argparse



def MassPlots(sim1name, sim1OutDir, sim2name, sim2OutDir, plotOutDir):

    simData = {}
    simData[sim1name] = {}
    simData[sim2name] = {}

    LineColor = {}
    LineColor[sim1name] = 'steelblue'
    LineColor[sim2name] = 'orange'

    LineStyle = {}
    LineStyle[sim1name] = '-'
    LineStyle[sim2name] = '-.'

    for k, simdir in zip([sim1name, sim2name], [sim1OutDir, sim2OutDir]):

        mass = TableReader(os.path.join(simdir, "Mass"))
        main_reader = TableReader(os.path.join(simdir, "Main"))

        simData[k]['cell'] = mass.readColumn("dryMass")
        simData[k]['protein'] = mass.readColumn("proteinMass")
        simData[k]['tRNA'] = mass.readColumn("tRnaMass")
        simData[k]['rRNA'] = mass.readColumn("rRnaMass")
        simData[k]['mRNA'] = mass.readColumn("mRnaMass")
        simData[k]['DNA'] = mass.readColumn("dnaMass")
        simData[k]['Small Molecules'] = mass.readColumn("smallMoleculeMass")

        initialTime = main_reader.readAttribute("initialTime")
        simData[k]['t'] = (main_reader.readColumn("time") - initialTime) / 60.

        simData[k]['masses'] = np.vstack([
            simData[k]['protein'],
            simData[k]['rRNA'],
            simData[k]['tRNA'],
            simData[k]['mRNA'],
            simData[k]['DNA'],
            simData[k]['Small Molecules'],
        ]).T
        simData[k]['fractions'] = (simData[k]['masses']/ simData[k]['cell'][:, None]).mean(axis=0)


    mass_labels = ["cell", "protein", "rRNA", "tRNA", "mRNA", "DNA", "Small Molecules"]

    # PLOT MASS FRACTIONS
    plt.figure(figsize=(11, 8.5))
    for p, mass in enumerate(mass_labels):

        plt.subplot(2, 4, p+1)
        for k in simData:
            plt.plot(simData[k]['t'], simData[k][mass] / simData[k][mass][0], linewidth=2, color=LineColor[k], linestyle=LineStyle[k], label=k)

        plt.title(mass)
        plt.xlabel("Time (min)")
        plt.ylabel("Mass (normalized by t = 0 min)")
        # plt.legend(legend, loc="best")

    plt.legend()
    plt.suptitle("Mass Fractions", fontweight="bold")
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    exportFigure(plt, plotOutDir, sim1name + '_' + sim2name + '_massFrac')

    plt.close("all")

    # PLOT ABSOLUTE MASSES
    plt.figure(figsize=(11, 8.5))
    for p, mass in enumerate(mass_labels):

        plt.subplot(2, 4, p+1)
        for k in simData:
            plt.plot(simData[k]['t'], simData[k][mass], linewidth=2, color=LineColor[k], linestyle=LineStyle[k], label=k)

        plt.title(mass)
        plt.xlabel("Time (min)")
        plt.ylabel("Mass")
        # plt.legend(legend, loc="best")

    plt.legend()
    plt.suptitle('Absolute Masses', fontweight="bold")
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    exportFigure(plt, plotOutDir, sim1name + '_' + sim2name + '_masses')
    plt.close("all")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('sim 1 name', type=str, help='name of first sim')
    parser.add_argument('sim1Dir', type=str, help='directory containing sim (str)')
    parser.add_argument('sim 2 name', type=str, help='name of second sim')
    parser.add_argument('sim2Dir', type=str, help='directory containing master sim (str)')
    args = vars(parser.parse_args())
    OutDir = args['sim1Dir'].split('simOut')[0] + 'plotOut/'
    MassPlots(args['sim 1 name'], args['sim1Dir'], args['sim 2 name'],  args['sim2Dir'], OutDir)
