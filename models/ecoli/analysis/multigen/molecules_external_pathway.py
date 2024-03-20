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

from models.ecoli.analysis import multigenAnalysisPlot
from wholecell.analysis.analysis_tools import (exportFigure,
	read_stacked_bulk_molecules, read_stacked_columns)
from wholecell.io.tablereader import TableReader

class Plot(multigenAnalysisPlot.MultigenAnalysisPlot):
    def do_plot(self, seedOutDir, plotOutDir, plotOutFileName, simDataFile,
				validationDataFile, metadata):
        sim_data = self.read_pickle_file(simDataFile)

        if not sim_data.process.metabolism_external_pathway.has_external_pathway:
            print("This simulation does not have an external pathway")
            return

        molecule_names = sim_data.process.metabolism_external_pathway.molecule_names

        cell_paths = self.ap.get_cells()
        sim_dir = cell_paths[0]
        plt.figure()
        simOutDir = os.path.join(sim_dir, 'simOut')

        # Load data
        time = read_stacked_columns(cell_paths, 'Main', 'time')
        bulk_molecules = read_stacked_columns(cell_paths, 'BulkMolecules', 'counts')

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
