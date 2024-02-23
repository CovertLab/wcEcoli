"""
Plot fluxes for the external metabolic pathway
"""

import pickle
import os

from matplotlib import pyplot as plt
# noinspection PyUnresolvedReferences
import numpy as np

from models.ecoli.analysis import singleAnalysisPlot
from models.ecoli.processes.metabolism import COUNTS_UNITS, VOLUME_UNITS, TIME_UNITS, MASS_UNITS
from wholecell.analysis.analysis_tools import exportFigure, read_bulk_molecule_counts
from wholecell.io.tablereader import TableReader
from wholecell.utils import units

FLUX_UNITS = COUNTS_UNITS / VOLUME_UNITS / TIME_UNITS
class Plot(singleAnalysisPlot.SingleAnalysisPlot):
    def do_plot(self, simOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
        sim_data = self.read_pickle_file(simDataFile)
        if not sim_data.process.metabolism_external_pathway.has_external_pathway:
            print("This simulation does not have an external pathway")
            return
        cell_density = sim_data.constants.cell_density

        # Listeners used
        main_reader = TableReader(os.path.join(simOutDir, 'Main'))
        initial_time = main_reader.readAttribute('initialTime')
        fba_results = TableReader(os.path.join(simOutDir, 'FBAResults'))
        mass_table = TableReader(os.path.join(simOutDir, 'Mass'))

        # Load data
        time = main_reader.readColumn('time') - initial_time
        fluxes_external_pathway = fba_results.readColumn('externalPathwayFluxes')

        # Calculate coefficients to be used to convert flux units from
        # mM/s to mmol/gCDW/h
        #cell_mass = mass_table.readColumn('cellMass')
        #dry_mass =  mass_table.readColumn('dryMass')
        #conversion_coeffs = (
        #        dry_mass / cell_mass
        #        * cell_density.asNumber(MASS_UNITS / VOLUME_UNITS)
        #)

        # Fluxes
        plt.figure()

        for m in range(fluxes_external_pathway.shape[1]):
            plt.plot(time / 60., fluxes_external_pathway[:, m], label='Flux '+ str(m))
        plt.xlabel("Time (min)")
        plt.ylabel("Fluxes (mmol/L/sec)")
        plt.title("External Metabolic Pathway Fluxes")
        plt.legend()

        ### Create Plot ###
        plt.tight_layout()
        exportFigure(plt, plotOutDir, plotOutFileName, metadata)
        plt.close('all')

if __name__ == '__main__':
    Plot().cli()
