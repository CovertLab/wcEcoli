"""
This plot allows one to visualize the protein synthesis and loss rates over
then duration of a simulation, along with a plot of the protein's free monomer
counts with a first order degradation rate equation fitted to it.
"""
import pickle
import os
from scipy.optimize import curve_fit
from matplotlib import pyplot as plt
# noinspection PyUnresolvedReferences
import numpy as np

from models.ecoli.analysis import multigenAnalysisPlot
from wholecell.analysis.analysis_tools import (exportFigure,
    read_bulk_molecule_counts, read_stacked_bulk_molecules, read_stacked_columns)
from wholecell.io.tablereader import TableReader
import wholecell.utils.units as units

PLOT_COMPLEXES = ["CPLX0-2881", "CPLX0-3101", "MONOMER0-160", "MONOMER0-155"]


class Plot(multigenAnalysisPlot.MultigenAnalysisPlot):
    def do_plot(self, seedOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
        with open(simDataFile, 'rb') as f:
            sim_data = pickle.load(f)
        cell_paths = self.ap.get_cells()
        sim_dir = cell_paths[0]
        simOutDir = os.path.join(sim_dir, 'simOut')

        # Extract protein indexes for each new gene
        monomer_counts_reader = TableReader(
            os.path.join(simOutDir, "MonomerCounts"))
        monomer_idx_dict = {monomer: i for i, monomer in
                            enumerate(monomer_counts_reader.readAttribute(
                                'monomerIds'))}
        # Extract the protein IDs
        monomerIDs = monomer_counts_reader.readAttribute("monomerIds")

        # Extract the free monomer counts using the monomer counts listener:
        free_monomer_counts = read_stacked_columns(cell_paths, 'MonomerCounts', "freeMonomerCounts")

        # Extract the complexes:
        complex_counts_listener = TableReader(
            os.path.join(simOutDir, "ComplexationListener"))
        complexIDs = complex_counts_listener.readAttribute("complexIDs")
        complex_idx_dict = {complexID: i for i, complexID in
                            enumerate(complexIDs)}
        eq_complex_counts_listener = TableReader(
            os.path.join(simOutDir, "EquilibriumListener"))
        eq_complexIDs = eq_complex_counts_listener.readAttribute("complexIDs")
        eq_complex_idx_dict = {complexID: i for i, complexID in
                            enumerate(eq_complexIDs)}


        # doubling time function from nora (note the normal doubling time extraction is not working):
        def extract_doubling_times(cell_paths):
            # Load simulation time span data
            time = read_stacked_columns(cell_paths, 'Main', 'time').squeeze()
            # Determine doubling time
            doubling_times = read_stacked_columns(cell_paths, 'Main', 'time',
                                                  fun=lambda x: (x[-1] - x[0])).squeeze().astype(
                int)
            end_generation_times = np.cumsum(doubling_times) + time[0]  #
            start_generation_indices = np.searchsorted(time, end_generation_times[:-1],
                                                       side='left').astype(int)
            start_generation_indices = np.insert(start_generation_indices, 0, 0) + np.arange(
                len(doubling_times))
            end_generation_indices = start_generation_indices + doubling_times
            return time, doubling_times, end_generation_times, start_generation_indices, end_generation_indices

        # extract the doubling times:
        time, doubling_times, end_generation_times, start_generation_indices, end_generation_indices = extract_doubling_times(
            cell_paths)

        def check_validity_and_get_compartment(protein_list):
            revised_protein_list = []
            for protein in protein_list:
                if "[" in protein:
                    protein = protein[:-3] # remove compartment
                if sim_data.getter.is_valid_molecule(protein):
                    revised_name = protein + sim_data.getter.get_compartment_tag(protein)
                    revised_protein_list.append(revised_name)

            return revised_protein_list

        # Make sure the proteins inputted are valid and have a compartment tag:
        PLOT_COMPLEXES_revised = check_validity_and_get_compartment(PLOT_COMPLEXES)

        # plot the loss rate and the production rate:
        for complex in PLOT_COMPLEXES_revised:
            # Determine the complex type:
            monomers = {}
            if complex in complex_idx_dict:
                complex_idx = complex_idx_dict[complex]
                # find the complex and its consituent monomers:
                molecules_to_all_downstream_complexes = sim_data.process.complexation.molecules_to_all_downstream_complexes_dict
                # find where the complex shows up as a value in the dict:
                monomers = {k: v for k, v in molecules_to_all_downstream_complexes.items() if complex in v}
                complex_counts = read_stacked_columns(cell_paths, 'ComplexationListener', "complexCounts")[:, complex_idx]
                complexation_events = read_stacked_columns(cell_paths, 'ComplexationListener', "complexationEvents")[:, complex_idx]
                complex_type = "complexation"
                hi = 5
                # TODO: find the correct way to index to these compelxation events from the reaction ids
                complex_reactions = complex_counts_listener.readAttribute("reactionIDs")
                # todo: fix alignment!
            elif complex in eq_complex_idx_dict:
                complex_idx = eq_complex_idx_dict[complex]
                # find the complex and its consituent monomers:
                molecules_to_all_downstream_complexes = sim_data.process.equilibrium.molecules_to_all_downstream_complexes_dict
                # find where the complex shows up as a value in the dict:
                monomers_temp = {k: v for k, v in molecules_to_all_downstream_complexes.items() if
                            complex in v}
                # remove any monomers that are not actually monomers (i.e. ATP):
                for key in monomers_temp.keys():
                    if key in monomerIDs:
                        monomers[key] = monomers_temp[key]

                complex_counts = read_stacked_columns(cell_paths, 'EquilibriumListener',
                                                      "complexCounts")[:, complex_idx]
                hi = 5
                complexation_events = read_stacked_columns(cell_paths, 'EquilibriumListener', "complexationEvents")[:, complex_idx]
                complex_type = "equilibrium"

            else:
                raise ValueError(f"Complex {complex} not found in complexation complex or equilibrium complex  listeners.")

            # Generate the plots:
            fig, (ax1, ax2, ax3) = plt.subplots(3, 1, sharex=True, figsize=(10, 6))
            # Complex Counts plot
            ax1.plot(time, complex_counts, color='lightseagreen', label='Complex Counts')

            # plot the counts of each monomer that makes up the complex:
            for monomer in monomers.keys():
                # get the counts of the monomer in the complex:
                monomer_info = monomers[monomer]
                monomer_complex_info = monomer_info[complex][1]
                monomer_stoich = monomer_complex_info['stoichiometry']
                monomer_complex_counts = complex_counts * monomer_stoich
                ax1.plot(time, monomer_complex_counts, alpha=0.5, label=f'{monomer} counts within {complex} ({monomer_stoich} per)')


            for i in range(len(end_generation_times)):
                dt = end_generation_times[i]
                ax1.axvline(x=dt, linestyle='--', color="yellowgreen")

            ax1.legend(fontsize=5)
            ax1.set_ylabel("Complex Counts")
            ax1.set_title(f"Complex counts over time for {complex_type} {complex}\n Sim ID: {metadata['description']}")

            # Now plot the free monomer counts for each protein that makes up the complex:
            for monomer in monomers.keys():
                protein_idx = monomer_idx_dict[monomer]
                protein_FMC = free_monomer_counts[:, protein_idx]
                ax2.plot(time, protein_FMC, alpha=0.5, label=f'{monomer}')

            # add vertical lines for the end of each generation:
            for i in range(len(end_generation_times)):
                dt = end_generation_times[i]
                ax2.axvline(x=dt, linestyle='--', color="yellowgreen")


            ax2.set_ylabel("Free Monomer Counts")
            ax2.legend(fontsize=5)

            # add a plot of the complexation events over time:
            ax3.plot(time, complexation_events, color='lightseagreen', label=complex)
            ax3.set_xlabel("Time (s)")
            ax3.set_ylabel("Complexation Events")
            ax3.legend(fontsize=5)

            # add vertical lines for the end of each generation:
            for i in range(len(end_generation_times)):
                dt = end_generation_times[i]
                ax3.axvline(x=dt, linestyle='--', color="yellowgreen")




            plt.tight_layout()





            #save the plot:
            file_name = plotOutFileName + "_" +complex
            exportFigure(plt, plotOutDir, file_name, metadata)


if __name__ == '__main__':
    Plot().cli()
