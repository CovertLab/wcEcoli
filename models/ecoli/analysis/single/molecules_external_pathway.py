"""
Molecule counts for the molecules involved in the external pathway
"""

import pickle
import os
import json
from matplotlib import pyplot as plt
from collections import OrderedDict
import math

# noinspection PyUnresolvedReferences
import numpy as np

from models.ecoli.analysis import singleAnalysisPlot
from wholecell.analysis.analysis_tools import exportFigure, read_bulk_molecule_counts
from wholecell.io.tablereader import TableReader


class Plot(singleAnalysisPlot.SingleAnalysisPlot):
    def do_plot(self, simOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
        sim_data = self.read_pickle_file(simDataFile)
        if not sim_data.process.metabolism_external_pathway.has_external_pathway:
            print("This simulation does not have an external pathway")
            return

        molecule_names = sim_data.process.metabolism_external_pathway.molecule_names

        # Listeners used
        main_reader = TableReader(os.path.join(simOutDir, 'Main'))
        initial_time = main_reader.readAttribute('initialTime')
        bulk_molecules_table = TableReader(os.path.join(simOutDir, 'BulkMolecules'))

        # Load data
        time = main_reader.readColumn('time') - initial_time
        bulk_molecules = bulk_molecules_table.readColumn('counts')

        with open(os.path.join(simOutDir, 'BulkMolecules', 'attributes.json')) as f:
            data = json.load(f, object_pairs_hook=OrderedDict)
            data_b = data['objectNames']
        # Molecule counts
        plt.figure()
        plt.xlabel("Time (min)")
        plt.title("Molecule counts")

        for m in range(len(molecule_names)):
            plt.subplot(int(math.sqrt(len(molecule_names)))+1, int(math.sqrt(len(molecule_names)))+1, m+1)
            molecule_couts = bulk_molecules[:, data_b.index(molecule_names[m])]
            plt.plot(time / 60., molecule_couts)
            plt.ylabel('\n'.join(molecule_names[m].split('-')))

        ### Create Plot ###
        plt.tight_layout()
        exportFigure(plt, plotOutDir, plotOutFileName, metadata)
        plt.close('all')

if __name__ == '__main__':
    Plot().cli()
