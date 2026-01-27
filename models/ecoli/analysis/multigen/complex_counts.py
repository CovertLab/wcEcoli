"""
This plot allows one to visualize the counts of a specified complex over time,
along with the counts of the free monomers that make up that complex,
and the complexation events that produce that complex.
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

PLOT_COMPLEXES = ["MONOMER0-160",] #
                  #"MONOMER0-160", "MONOMER0-155", ]


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

        def check_validity_and_get_compartment(molecule_list):
            revised_molecule_list = []
            for molecule in molecule_list:
                if "[" in molecule:
                    molecule = molecule[:-3] # remove compartment
                if sim_data.getter.is_valid_molecule(molecule):
                    revised_name = molecule + sim_data.getter.get_compartment_tag(molecule)
                    revised_molecule_list.append(revised_name)

            return revised_molecule_list

        # Make sure the proteins inputted are valid and have a compartment tag:
        PLOT_COMPLEXES_revised = check_validity_and_get_compartment(PLOT_COMPLEXES)

        # Generate a plot for each complex:
        for complex in PLOT_COMPLEXES_revised:
            # Determine the complex type:
            monomers = {}
            parent_complexes = {}
            if complex in complex_idx_dict:
                complex_type = "complexation"
                complex_idx = complex_idx_dict[complex]
                # find the complex and its consituent monomers:
                molecules_to_all_downstream_complexes = sim_data.process.complexation.molecules_to_all_downstream_complexes_dict
                # determine if the complex itself is a subunit of a larger complex:
                # TODO: figure out how to search for this in any complex type! Maybe make one large dictionary that combines both complexation and equilibrium?
                molecules_to_parent_complexes = sim_data.process.complexation.molecules_to_parent_complexes_dict
                if complex in molecules_to_parent_complexes.keys():
                    parent_complexes = molecules_to_parent_complexes[complex]
                # find where the complex shows up as a value in the dict:
                monomers = {k: v for k, v in molecules_to_all_downstream_complexes.items() if complex in v}
                complex_counts = read_stacked_columns(cell_paths, 'ComplexationListener', "complexCounts")[:, complex_idx]

                # TODO: determine if any complexes are generated by multiple reactions (if so, probalby need to generate a dictionary)
                # Find which reactions this complex is involved in:
                complex_reactions = complex_counts_listener.readAttribute("reactionIDs")
                reaction_to_idx_dict = {}
                complex_makeup = ''
                for monomer in monomers:
                    monomer_info = monomers[monomer]
                    monomer_complex_info = monomer_info[complex][0]
                    reaction_id = monomer_complex_info['reaction_id']
                    complex_makeup = monomer_info[complex][3]['complex_type']
                    if reaction_id in complex_reactions:
                        reaction_idx = complex_reactions.index(reaction_id)
                        reaction_to_idx_dict[reaction_id] = reaction_idx

            elif complex in eq_complex_idx_dict:
                complex_type = "equilibrium"
                complex_idx = eq_complex_idx_dict[complex]
                # find the complex and its consituent monomers:
                molecules_to_all_downstream_complexes = sim_data.process.equilibrium.molecules_to_all_downstream_complexes_dict
                # determine if the complex itself is a subunit of a larger complex: # TODO: figure out how to search for this in any complex type! Maybe make one large dictionary that combines both complexation and equilibrium?
                molecules_to_parent_complexes = sim_data.process.equilibrium.molecules_to_parent_complexes_dict
                if complex in molecules_to_parent_complexes.keys():
                    parent_complexes = molecules_to_parent_complexes[complex]
                # find where the complex shows up as a value in the dict:
                monomers_temp = {k: v for k, v in molecules_to_all_downstream_complexes.items() if
                            complex in v}
                # remove any monomers that are not actually monomers (i.e. ATP):
                for key in monomers_temp.keys():
                    if key in monomerIDs:
                        monomers[key] = monomers_temp[key]

                complex_counts = read_stacked_columns(cell_paths, 'EquilibriumListener',
                                                      "complexCounts")[:, complex_idx]

                # find which reactions this complex is involved in:
                complex_reactions = eq_complex_counts_listener.readAttribute("reactionIDs")
                reaction_to_idx_dict = {}
                complex_makeup = ''
                for monomer in monomers:
                    monomer_info = monomers[monomer]
                    monomer_complex_info = monomer_info[complex][0]
                    reaction_id = monomer_complex_info['reaction_id']
                    complex_makeup = monomer_info[complex][3]['complex_type']
                    if reaction_id in complex_reactions:
                        reaction_idx = complex_reactions.index(reaction_id)
                        reaction_to_idx_dict[reaction_id] = reaction_idx

            else:
                # TODO: determine if any complexes are in both complexation and equilibrium (if so, need to handle that here)?
                raise ValueError(f"Complex {complex} not found in complexation complex or equilibrium complex  listeners.")

            # Generate the plots:
            fig, (ax1, ax2, ax3) = plt.subplots(3, 1, sharex=True, figsize=(10, 6))

            # First, make a plot of the complex itself and its subunits (as well as the complexes it forms)
            ax1.plot(time, complex_counts, color='lightseagreen', label=f'{complex} free complex counts', linewidth=.75, alpha=0.75)

            # plot the counts of each monomer that makes up the complex:
            for monomer in monomers.keys():
                # get the counts of the monomer in the complex:
                monomer_info = monomers[monomer]
                monomer_complex_info = monomer_info[complex][1] # use 1 to get the stoichiometry dict
                monomer_stoich = monomer_complex_info['stoichiometry']
                monomer_complex_counts = complex_counts * monomer_stoich
                # check if the monomer is the result of the complex or a different parent complex:
                monomer_parent_complexes = molecules_to_parent_complexes[monomer]
                if complex in monomer_parent_complexes.keys():
                    # if the complex is a parent complex of the monomer, plot with a different style:
                    ax1.plot(time, monomer_complex_counts, alpha=1, label=f'{monomer} counts within\n{complex} ({monomer_stoich} per)', linewidth=1, linestyle=":")
                else:
                    for pc in monomer_parent_complexes.keys():
                        gparent_complexes = molecules_to_parent_complexes[pc]
                        if complex in gparent_complexes.keys():
                            # if the complex is a parent complex of the monomer's parent complex, note that:
                            ax1.plot(time, monomer_complex_counts, alpha=1, label=f'{monomer} counts within\n{complex} ({monomer_stoich} per, via {pc})', linewidth=.8, linestyle="-.")
                            break
                        else:
                            # if the complex is not a parent complex of the monomer or its parent complexes:
                            ax1.plot(time, monomer_complex_counts, alpha=1,
                                     label=f'{monomer} counts within\n{complex} ({monomer_stoich} per, via an unknown intermediate)',
                                     linewidth=.8, linestyle="-.")

            # Next, plot the free monomer counts for each protein that makes up the complex:
            for monomer in monomers.keys():
                protein_idx = monomer_idx_dict[monomer]
                protein_FMC = free_monomer_counts[:, protein_idx]
                ax2.plot(time, protein_FMC, alpha=0.75, label=f'{monomer}', linewidth=0.75)

            ax2.set_ylabel(f"Free Monomer \nSubunit Counts")
            ax2.legend(fontsize=5, loc="center left", bbox_to_anchor=(1, 0.5))

            # Third, add a plot of the complexation events over time:
            for rxn in reaction_to_idx_dict.keys():
                rxn_idx = reaction_to_idx_dict[rxn]
                if complex_type == "equilibrium":
                    # NOTE: equilibrium reactions can go negative (reverse reactions), so the events here may be negative
                    rxn_events = read_stacked_columns(cell_paths, 'EquilibriumListener', "complexationEvents")[:, rxn_idx]
                else:
                    rxn_events = read_stacked_columns(cell_paths, 'ComplexationListener', "complexationEvents")[:, rxn_idx]
                ax3.plot(time, rxn_events, alpha=0.75, label=f'{rxn} \n(generates 1 {complex})', linewidth=0.75)


            # If the complex is itself a subunit of a larger complex, plot the
            # counts of the parent complex and the events that generate it as well:
            if parent_complexes != {}:
                for parent_complex in parent_complexes:
                    if complex_type == "complexation":
                        # find the counts of complex within the parent complex:
                        parent_complex_idx = complex_idx_dict[parent_complex]
                        parent_complex_counts = read_stacked_columns(cell_paths,
                                                                     'ComplexationListener',
                                                                     "complexCounts")[:,
                                                parent_complex_idx]
                        complex_stoich = parent_complexes[parent_complex][1]['stoichiometry']
                        complex_counts_within_parent_complex = parent_complex_counts * complex_stoich

                        # find the events that generate the parent complex:
                        parent_complex_reaction_ID = parent_complexes[parent_complex][0]['reaction_id']
                        reaction_idx = complex_reactions.index(parent_complex_reaction_ID)
                        parent_complex_events = read_stacked_columns(cell_paths,
                                                                     'ComplexationListener',
                                                                     "complexationEvents")[:,
                                                reaction_idx]
                    else:
                        parent_complex_idx = eq_complex_idx_dict[parent_complex]
                        parent_complex_counts = read_stacked_columns(cell_paths,
                                                                     'EquilibriumListener',
                                                                     "complexCounts")[:,
                                                parent_complex_idx]
                        complex_stoich = parent_complexes[parent_complex][1]['stoichiometry']
                        complex_counts_within_parent_complex = parent_complex_counts * complex_stoich

                        # find the events that generate the parent complex:
                        parent_complex_reaction_ID = parent_complexes[parent_complex][0][
                            'reaction_id']
                        reaction_idx = complex_reactions.index(parent_complex_reaction_ID)
                        parent_complex_events = read_stacked_columns(cell_paths,
                                                                     'ComplexationListener',
                                                                     "complexationEvents")[:,
                                                reaction_idx]


                    ax1.plot(time, complex_counts_within_parent_complex, alpha=0.5,
                             label=f'{complex} counts within\n{parent_complex} ({complex_stoich} per)', linewidth=0.75,
                             linestyle="-.")
                    ax3.plot(time, np.negative(parent_complex_events), alpha=0.5,
                             label=f'{parent_complex_reaction_ID} \n(consumes {complex_stoich} {complex} to generate {parent_complex})',
                             linewidth=0.75)


            # Set the axes and legends (since some might have been added with the parent complexes):
            ax1.legend(fontsize=5, loc="center left", bbox_to_anchor=(1, 0.5))
            ax1.set_ylabel("Complex Counts")
            ax1.set_title(
                f"Complex counts for the {complex_makeup} {complex_type} complex {complex}\n Sim ID: {metadata['description']}")
            ax3.set_xlabel("Time (s)")
            ax3.set_ylabel("Complexation Events")
            ax3.legend(fontsize=5, loc="center left", bbox_to_anchor=(1, 0.5))

            # Add vertical lines for the end of each generation:
            for i in range(len(end_generation_times)):
                dt = end_generation_times[i]
                a = 0.5
                ax1.axvline(x=dt, linestyle='--', color="yellowgreen", alpha=a)
                ax2.axvline(x=dt, linestyle='--', color="yellowgreen", alpha=a)
                ax3.axvline(x=dt, linestyle='--', color="yellowgreen", alpha=a)


            plt.tight_layout()

            # Save the plot:
            file_name = plotOutFileName + "_" +complex
            exportFigure(plt, plotOutDir, file_name, metadata)


if __name__ == '__main__':
    Plot().cli()
