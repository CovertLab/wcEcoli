"""
Template for multigen analysis plots
"""

import pickle
import os

from matplotlib import pyplot as plt
# noinspection PyUnresolvedReferences
import numpy as np

from models.ecoli.analysis import multigenAnalysisPlot
from wholecell.analysis.analysis_tools import (exportFigure,
    read_bulk_molecule_counts, read_stacked_bulk_molecules, read_stacked_columns)
from wholecell.io.tablereader import TableReader

PLOT_PROTEINS = ["G6890-MONOMER[c]",
                 "PD03938[c]",
                 "G6737-MONOMER[c]",
                 "RPOD-MONOMER[c]",
                 "PD02936[c]",
                 "RED-THIOREDOXIN2-MONOMER[c]",
                 "EG10542-MONOMER[c]"]

class Plot(multigenAnalysisPlot.MultigenAnalysisPlot):
    def do_plot(self, seedOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
        with open(simDataFile, 'rb') as f:
            sim_data = pickle.load(f)

        cell_paths = self.ap.get_cells()
        sim_dir = cell_paths[0]
        simOutDir = os.path.join(sim_dir, 'simOut')

        def check_validity_and_get_compartment(protein_list):
            revised_protein_list = []
            for protein in protein_list:
                if "[" in protein:
                    protein = protein[:-3] # remove compartment
                if sim_data.getter.is_valid_molecule(protein):
                    revised_name = protein + sim_data.getter.get_compartment_tag(protein)
                    revised_protein_list.append(revised_name)

            return revised_protein_list

        # Load data - degradation information stored in the listener
        rawDegradationRate = read_stacked_columns(cell_paths, 'MonomerCounts', "rawDegradationRate")
        rawDegnProtein = read_stacked_columns(cell_paths, 'MonomerCounts', "rawDegradationProtein")
        activeDegradationRate = read_stacked_columns(cell_paths, 'MonomerCounts', "activeDegradationRate")
        activeDegnProtein = read_stacked_columns(cell_paths, 'MonomerCounts', "activeDegnProtein")

        monomer_counts_reader = TableReader( os.path.join(simOutDir, "MonomerCounts") )
        monomer_idx_dict = {monomer: i for i, monomer in
                            enumerate(monomer_counts_reader.readAttribute(
                                'monomerIds'))}
        # Make sure the proteins inputted are valid and have a compartment tag:
        PLOT_PROTEINS_revised = check_validity_and_get_compartment(PLOT_PROTEINS)
        protein_idx = [monomer_idx_dict[protein] for protein in PLOT_PROTEINS_revised]

        POI_rawDegradationRate = rawDegradationRate[protein_idx] #POI==protein of interest
        POI_rawDegnProtein = rawDegnProtein[protein_idx]
        POI_activeDegradationRate = activeDegradationRate[protein_idx]
        POI_activeDegnProtein = activeDegnProtein[protein_idx]

        plt.figure()

        ### Create Plot ###

        plt.tight_layout()
        exportFigure(plt, plotOutDir, plotOutFileName, metadata)
        plt.close('all')


if __name__ == '__main__':
    Plot().cli()
