"""
This plot allows one to visualize a specified protein's counts over time,
including both its free monomeric form and its presence within complexes.
"""

import pickle
import os
from matplotlib import pyplot as plt
import numpy as np
from wholecell.io import tsv
import io
from wholecell.utils.filepath import ROOT_PATH
from models.ecoli.analysis import multigenAnalysisPlot
from wholecell.analysis.analysis_tools import (exportFigure,
     read_stacked_columns)
from wholecell.io.tablereader import TableReader

# Specify proteins to plot (each protein will have its own plot):
PLOT_PROTEINS = ["EG11088-MONOMER",]
# TODO: replace common name functions with actual data in sim_data (added in translation.py) when available!
def get_gene_symbols_for_monomer_ids():
    """
    Extracts the gene symbols for each monomer id in the model.

    Returns: a dictionary mapping monomer ids to gene symbols.
    Code adapted from convert_to_flat.py.
    """
    RNAS_FILE = os.path.join(ROOT_PATH, 'reconstruction', 'ecoli',
                                 'flat', 'rnas.tsv')
    with (io.open(RNAS_FILE, 'rb') as f):
        reader = tsv.reader(f, delimiter='\t')
        headers = next(reader)
        while headers[0].startswith('#'):
            headers = next(reader)

        # extract relevant information
        gene_symbol_index = headers.index('common_name')
        protein_id_index = headers.index('monomer_ids')
        monomer_ids_to_gene_symbols = {}
        for line in reader:
            gene_symbol = line[gene_symbol_index]
            protein_id = list(
                line[protein_id_index][2:-2].split('", "'))[0]
            monomer_ids_to_gene_symbols[protein_id] = gene_symbol

    return monomer_ids_to_gene_symbols

def get_common_name(protein_id):
    """
    Obtains the common names for each protein of interest
    Args:
        protein_id: the name of the protein(s) of interest

    Returns: the common name for the protein(s) of interest
    """
    if protein_id == 'NG-GFP-MONOMER[c]':
        return 'GFP'

    else:
        protein = protein_id[:-3]  # subtract the compartment
        common_name = get_gene_symbols_for_monomer_ids()[protein]

    return common_name

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

        # Extract the total and free monomer counts using the monomer counts listener:
        total_monomer_counts = read_stacked_columns(cell_paths, 'MonomerCounts', "monomerCounts")
        free_monomer_counts = read_stacked_columns(cell_paths, 'MonomerCounts', "freeMonomerCounts")

        # Make the complex listeners (but do not use these to read the counts since they only take in generation 0):
        complex_counts_listener = TableReader(
            os.path.join(simOutDir, "ComplexationListener"))
        eq_complex_counts_listener = TableReader(
            os.path.join(simOutDir, "EquilibriumListener"))


        # Extract the doubling times and generation start/end times (using Nora's method):
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

        def check_validity_and_get_compartment(molecule_list):
            revised_molecule_list = []
            for molecule in molecule_list:
                if "[" in molecule:
                    molecule = molecule[:-3]  # remove compartment
                if sim_data.getter.is_valid_molecule(molecule):
                    revised_name = molecule + sim_data.getter.get_compartment_tag(molecule)
                    revised_molecule_list.append(revised_name)

            return revised_molecule_list

        # Make sure the proteins inputted are valid and have a compartment tag:
        PLOT_PROTEINS_revised = check_validity_and_get_compartment(PLOT_PROTEINS)

        # Generate a plot for each specified protein:
        for protein in PLOT_PROTEINS_revised:

            # Generate the plots:
            fig, (ax1, ax2, ax3, ax4, ax5) = plt.subplots(5, 1, sharex=True, figsize=(20, 12))

            # First plot the total monomer counts:
            ax1.plot(time, total_monomer_counts[:,monomer_idx_dict[protein]], color='orange',label=f"{protein}", linewidth=.75)
            ax1.set_ylabel("Total Monomer Counts")
            ax1.legend(fontsize=5, loc="center left", bbox_to_anchor=(1, 0.5))

            # Second, plot the free monomer counts:
            ax2.plot(time, free_monomer_counts[:,monomer_idx_dict[protein]], color='salmon', label=f"{protein}", linewidth=.75)
            ax2.set_ylabel("Free Monomer Counts")
            ax2.legend(fontsize=5, loc="center left", bbox_to_anchor=(1, 0.5))

            # Third, plot the complexed monomer counts:
            complexation_complexes = {}
            eq_complexes = {}

            # find all complexes that the protein is involved in:
            molecules_to_all_downstream_complexes = sim_data.process.complexation.molecules_to_all_downstream_complexes_dict
            molecules_to_all_downstream_eq_complexes = sim_data.process.equilibrium.molecules_to_all_downstream_complexes_dict

            # Find all complexation complexes involving this protein:
            if protein in molecules_to_all_downstream_complexes.keys():
                molecules_to_parent_complexes = sim_data.process.complexation.molecules_to_parent_complexes_dict
                parent_complexes = molecules_to_parent_complexes[protein]
                complexes = molecules_to_all_downstream_complexes[protein]
                complexation_complexes = complexes
                for complex in complexes.keys():
                    # find the complex counts:
                    complex_idx = complex_counts_listener.readAttribute("complexIDs").index(complex)
                    complex_counts = read_stacked_columns(cell_paths, 'ComplexationListener', "complexCounts")[:, complex_idx]
                    line = ax3.plot(time, complex_counts, alpha=0.5, label=f'complexation complex {complex}', linewidth=0.75)
                    # find the number of monomers per this complex:
                    monomer_stoich = complexes[complex][1]['stoichiometry']
                    complexed_monomer_counts = complex_counts * monomer_stoich
                    line_color = line[0].get_color() # use the same color as the complex line
                    # determine if the complex is the direct parent complex of the protein:
                    if complex in parent_complexes.keys():
                        ax3.plot(time, complexed_monomer_counts, alpha=1, label=f'{protein} counts within\n{complex} ({monomer_stoich} per)', linewidth=1, linestyle=':', color=line_color)
                    else:
                        # the complex is not the direct parent complex of the protein, so check to see which is the subunit root for this complex:
                        parent_complex_to_plot = 0
                        for pc in parent_complexes.keys():
                            grandparent_complexes = molecules_to_parent_complexes[pc]
                            if complex in grandparent_complexes.keys():
                                print(f"Note: {protein} is in complex {complex} via {pc}")
                                parent_complex_to_plot = pc
                            else:
                                print(f"Note: {protein} is in complex {complex} via unknown route")
                        if parent_complex_to_plot !=0:
                            ax3.plot(time, complexed_monomer_counts, alpha=1, label=f'{protein} counts within\n{complex} ({monomer_stoich} per, via {parent_complex_to_plot})', linewidth=1, linestyle=':', color=line_color)
                        else:
                            ax3.plot(time, complexed_monomer_counts, alpha=1, label=f'{protein} counts within\n{complex} ({monomer_stoich} per, via unknown route)', linewidth=1, linestyle=':', color=line_color)

            # Find all equilibrium complexes involving this protein:
            if protein in molecules_to_all_downstream_eq_complexes.keys():
                complexes = molecules_to_all_downstream_eq_complexes[protein]
                eq_complexes = complexes
                for complex in complexes.keys():
                    # find the complex counts:
                    complex_idx = eq_complex_counts_listener.readAttribute("complexIDs").index(complex)
                    complex_counts = read_stacked_columns(cell_paths, 'EquilibriumListener', "complexCounts")[:, complex_idx]
                    line = ax3.plot(time, complex_counts, alpha=0.5, label=f'equilibrium complex {complex}', linewidth=0.75)
                    # find the number of monomers per this complex:
                    monomer_stoich = complexes[complex][1]['stoichiometry']
                    complexed_monomer_counts = complex_counts * monomer_stoich
                    line_color = line[0].get_color()
                    ax3.plot(time, complexed_monomer_counts, alpha=1, label=f'{protein} counts within\n{complex} ({monomer_stoich} per)', linewidth=1, linestyle=':', color=line_color)

            ax3.set_ylabel("Complexed Monomer Counts")
            ax3.legend(fontsize=5, loc="center left", bbox_to_anchor=(1, 0.5))

            # Fourth, plot the total, free, and complexed monomer counts together:
            ax4.plot(time, total_monomer_counts[:,[monomer_idx_dict[protein]]], color='orange',label=f"total monomer counts", linewidth=0.9)
            ax4.plot(time, free_monomer_counts[:,[monomer_idx_dict[protein]]], color='salmon', label=f"free monomer counts", linewidth=0.9)

            # calculate complexed monomer counts:
            if len(complexation_complexes) > 0:
                complexation_complexed_monomer_counts = read_stacked_columns(cell_paths, 'ComplexationListener', "complexedMonomerCounts")[:, monomer_idx_dict[protein]]
                ax4.plot(time, complexation_complexed_monomer_counts, color='deepskyblue', label=f"monomers in complexation complexes", linewidth=.5)
            if len(eq_complexes) > 0:
                eq_complexed_monomer_counts = read_stacked_columns(cell_paths, 'EquilibriumListener', "complexedMonomerCounts")[:, monomer_idx_dict[protein]]
                ax4.plot(time, eq_complexed_monomer_counts, color='mediumseagreen', label=f"monomers in equilibrium complexes", linewidth=.5)

            ax4.set_ylabel("All Count Types")
            ax4.legend(fontsize=5, loc="center left", bbox_to_anchor=(1, 0.5))

            # Fifth, plot the complexation events for each reaction involving this protein:
            if len(complexation_complexes) > 0:
                for complex in complexation_complexes.keys():
                    reaction_ID = complexation_complexes[complex][0]['reaction_id']
                    rxn_idx = complex_counts_listener.readAttribute("reactionIDs").index(reaction_ID)
                    complexation_events = read_stacked_columns(cell_paths, 'ComplexationListener', "complexationEvents")[:, rxn_idx]
                    ax5.plot(time, complexation_events, alpha=0.75, label=f'{reaction_ID} ({complex})', linewidth=0.5)
            if len(eq_complexes) > 0:
                for complex in eq_complexes.keys():
                    reaction_ID = eq_complexes[complex][0]['reaction_id']
                    rxn_idx = eq_complex_counts_listener.readAttribute("reactionIDs").index(reaction_ID)
                    eq_events = read_stacked_columns(cell_paths, 'EquilibriumListener', "complexationEvents")[:, rxn_idx]
                    ax5.plot(time, eq_events, alpha=0.75, label=f'{reaction_ID} ({complex})', linewidth=0.5)

            ax5.set_xlabel("Time (s)")
            ax5.set_ylabel("Complexation Events")
            ax5.legend(fontsize=5, loc="center left", bbox_to_anchor=(1, 0.5))

            # Plot specifications:
            plt.tight_layout()
            # TODO: edit how the common name is extracted when sim_data has this info:
            ax1.set_title(f"Monomer count types for {protein} ({get_common_name(protein)})"
                          f"\n Sim ID: {metadata['description']}")
            plt.tight_layout(pad=3.0)

            # add vertical lines for the end of each generation:
            for i in range(len(end_generation_times)):
                dt = end_generation_times[i]
                a = 0.5
                ax1.axvline(x=dt, linestyle='--', color="yellowgreen", alpha=a)
                ax2.axvline(x=dt, linestyle='--', color="yellowgreen", alpha=a)
                ax3.axvline(x=dt, linestyle='--', color="yellowgreen", alpha=a)
                ax4.axvline(x=dt, linestyle='--', color="yellowgreen", alpha=a)
                ax5.axvline(x=dt, linestyle='--', color="yellowgreen", alpha=a)


            # Save the plot:
            file_name = plotOutFileName + "_" + protein
            exportFigure(plt, plotOutDir, file_name, metadata)


if __name__ == '__main__':
    Plot().cli()
