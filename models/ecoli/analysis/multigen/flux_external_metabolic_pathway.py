"""
Plot fluxes for the reactions in the external metabolic pathway across multiple generations
"""

import pickle
import os

from matplotlib import pyplot as plt
# noinspection PyUnresolvedReferences
import numpy as np

from models.ecoli.analysis import multigenAnalysisPlot
from wholecell.analysis.analysis_tools import (exportFigure,
	read_stacked_bulk_molecules, read_stacked_columns)
from wholecell.io.tablereader import TableReader

from models.ecoli.processes.metabolism import (
	COUNTS_UNITS, VOLUME_UNITS, TIME_UNITS, MASS_UNITS)
from wholecell.utils import units


class Plot(multigenAnalysisPlot.MultigenAnalysisPlot):
    def do_plot(self, seedOutDir, plotOutDir, plotOutFileName, simDataFile,
				validationDataFile, metadata):
        with open(simDataFile, 'rb') as f:
            sim_data = pickle.load(f)

        if not sim_data.process.metabolism_external_pathway.has_external_pathway:
            print("This simulation does not have an external pathway")
            return

        cell_paths = self.ap.get_cells()

        sim_dir = cell_paths[0]
        plt.figure()
        #for sim_dir in cell_paths:
        simOutDir = os.path.join(sim_dir, 'simOut')
        # Load data
        time = read_stacked_columns(cell_paths, 'Main', 'time')
        fluxes_external_pathway = read_stacked_columns(cell_paths, 'FBAResults', 'externalPathwayFluxes')
        # Fluxes
        for m in range(fluxes_external_pathway.shape[1]):
            plt.plot(time / 60., fluxes_external_pathway[:, m],
                     label='Flux ' + str(m))
        plt.xlabel("Time (min)")
        plt.ylabel("Fluxes (mmol/L/sec)")

        cell_density = sim_data.constants.cell_density
        plt.figure()

        # Calculate coefficients to be used to convert flux units from
        # mM/s to mmol/gCDW/h
        cell_mass = read_stacked_columns(
            cell_paths, 'Mass', 'cellMass', ignore_exception=True)
        dry_mass = read_stacked_columns(
            cell_paths, 'Mass', 'dryMass', ignore_exception=True)
        conversion_coeffs = (
                dry_mass / cell_mass
                * cell_density.asNumber(MASS_UNITS / VOLUME_UNITS)
        )

        # Load data
        time = read_stacked_columns(cell_paths, 'Main', 'time')
        fluxes_external_pathway = ((COUNTS_UNITS / MASS_UNITS / TIME_UNITS)
				* (read_stacked_columns(cell_paths, 'FBAResults', 'externalPathwayFluxes')/ conversion_coeffs)
				).asNumber(units.mmol / units.g / units.h)
        reaction_id = sim_data.process.metabolism_external_pathway.rxn_ids

        # Fluxes plots
        for m in range(fluxes_external_pathway.shape[1]):
            plt.plot(time / 60., fluxes_external_pathway[:, m],
                     label='Flux ' + str(reaction_id[m]))
        plt.xlabel("Time (min)")
        plt.ylabel("Fluxes (mmol/gCDW/h)")
        plt.title("External Metabolic Pathway Fluxes")
        plt.legend()
        ### Create Plot ###
        plt.tight_layout()
        exportFigure(plt, plotOutDir, plotOutFileName, metadata)
        plt.close('all')


if __name__ == '__main__':
    Plot().cli()
