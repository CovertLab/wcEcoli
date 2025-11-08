"""
This plot allows one to visualize the spread of proteins over a simulation.
"""
import pickle
import os
from scipy.optimize import curve_fit
from matplotlib import pyplot as plt
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize

import seaborn as sns
# noinspection PyUnresolvedReferences
import numpy as np
import pandas as pd
from models.ecoli.analysis import multigenAnalysisPlot
from wholecell.analysis.analysis_tools import (exportFigure,
    read_bulk_molecule_counts, read_stacked_bulk_molecules, read_stacked_columns)
from wholecell.io.tablereader import TableReader


PLOT_PROTEINS = ["EG10542-MONOMER[c]"]

PLOT_COMPLEXES = ["CPLX0-2881[c]"]


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
        fmcs = read_stacked_columns(cell_paths, 'MonomerCounts', "freeMonomerCounts")
        (free_monomer_counts,) = read_stacked_bulk_molecules(
            cell_paths, monomerIDs)  # todo: check if this is equal to the above! so that way only one method for calling in the data is needed.

        # Extract the complex IDs:
        complex_reader = TableReader(os.path.join(simOutDir, "ComplexationListener"))
        complex_idx_dict = {complex: i for i, complex in
                            enumerate(complex_reader.readAttribute(
                                'complexIDs'))}
        complexIDs = complex_reader.readAttribute("complexIDs")
        (complex_counts,) = read_stacked_bulk_molecules(
            cell_paths, complexIDs)
        # TODO: in the future, add the functionaility to unpack the complexes and plot the monomer counts of the proteins within them as well.


        # read in the counts to molar conversion (at each time point):
        ectm = read_stacked_columns(cell_paths, 'EnzymeKinetics', "countsToMolar")

        # calculate the free monomer concentrations:
        free_monomer_concentrations = fmcs * ectm # mmol/L

        # calculate the complex concentrations:
        complex_concentrations = complex_counts * ectm # mmol/L

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
        PLOT_PROTEINS_revised = check_validity_and_get_compartment(PLOT_PROTEINS)
        PLOT_COMPLEXES_revised = check_validity_and_get_compartment(PLOT_COMPLEXES)
        n_gens = metadata["generations"]

        def make_generation_mask(cell_paths, n_gens):
            time, doubling_times, end_generation_times, start_generation_indices, end_generation_indices = extract_doubling_times(
                cell_paths)

            # initialize a mask:
            gen_mask = np.full(len(time), '', dtype=object)

            # find the time spanned by each generation:
            for i in range(n_gens):
                start_idx = start_generation_indices[i]
                end_idx = end_generation_indices[i]
                generation_name = f"Generation {i+1}"

                # within the generation, assign gen_mask to the generation_name:
                gen_mask[(time >= time[start_idx]) & (time <= time[end_idx])] = generation_name

            return gen_mask

        # define a function that creates a generation mask that has the time start at 0 for each generation:
        def generation_time_mask(cell_paths, n_gens):
            time, doubling_times, end_generation_times, start_generation_indices, end_generation_indices = extract_doubling_times(
                cell_paths)

            # initialize a mask:
            gen_time_mask = np.full(len(time), np.nan)

            # find the time spanned by each generation:
            for i in range(n_gens):
                start_idx = start_generation_indices[i]
                end_idx = end_generation_indices[i]

                # within the generation, assign gen_time_mask to the time since start of generation:
                gen_time_mask[(time >= time[start_idx]) & (time <= time[end_idx])] = time[
                                                                                          (time >= time[start_idx]) & (
                                                                                                  time <= time[
                                                                                                      end_idx])] - time[start_idx]

            return gen_time_mask

        # create the generation name mask and generation time mask:
        gen_mask = make_generation_mask(cell_paths, n_gens)
        gen_time_mask = generation_time_mask(cell_paths, n_gens)

        # plot the loss rate and the production rate:
        for protein in PLOT_PROTEINS_revised:
            # Gather the relevant data:
            protein_idx = monomer_idx_dict[protein]
            fmconc = free_monomer_concentrations[:, protein_idx]
            data = {'generation': gen_mask,
                    'fmconc': fmconc,
                    'all_time': time,
                    'Time (s)': gen_time_mask}
            p_df = pd.DataFrame(data)

            # Generate the plots:
            sns.violinplot(x='generation', y='fmconc', data=p_df, inner="points",
                           color='lightgray', alpha=0.7)

            sns.stripplot(x='generation',
                          y='fmconc',
                          data=p_df,
                          size=3,
                          hue='Time (s)',
                          palette='crest',
                          alpha=0.3,
                          jitter=0.3,
                          #dodge=True # this makes the code take forever
                          )

            # Color mapping and normalization
            norm = Normalize(vmin=p_df['Time (s)'].min(), vmax=p_df['Time (s)'].max())
            sm = ScalarMappable(cmap='crest', norm=norm)  # Create a ScalarMappable

            # Create a color bar
            cbar = plt.colorbar(sm)
            cbar.set_label('Time (s)')

            # Plot specifics:
            sim_name = metadata["description"]
            degradation_rate_combo = sim_data.protein_degradation_combo_selection

            plt.xlabel('Generation')
            plt.ylabel('Free Monomer Concentration (mmol/L)')
            plt.title(
                f"Free monomer concentration for {protein}\n Sim ID: {sim_name}; \nDegradation rate combo: {degradation_rate_combo}",
                size=12)
            plt.tight_layout()

            # save the plot:
            file_name = plotOutFileName + "_" + protein
            exportFigure(plt, plotOutDir, file_name, metadata)

        # also make a plot of the free monomer counts under it:
        for protein in PLOT_PROTEINS_revised:
            # Gather the relevant data:
            protein_idx = monomer_idx_dict[protein]
            protein_FMC = free_monomer_counts[:, protein_idx]
            fmconc = free_monomer_concentrations[:, protein_idx]
            data = {'generation': gen_mask,
                    'fmconc': fmconc,
                    'all_time': time,
                    'Time (s)': gen_time_mask}
            p_df = pd.DataFrame(data)

            # Generate the plots:
            fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(6*n_gens,10))

            # Monomer Counts plot
            ax1.plot(time, protein_FMC, color='lightseagreen', label='Free Monomer Counts')
            for i in range(len(end_generation_times)):
                dt = end_generation_times[i]
                ax1.axvline(x=dt, linestyle='--', color="yellowgreen")

            # Plot specifics:
            sim_name = metadata["description"]
            degradation_rate_combo = sim_data.protein_degradation_combo_selection
            ax1.legend(fontsize=5)
            ax1.set_ylabel("Free Monomer Counts")
            ax1.set_xlabel("Time (s)")
            ax1.set_title(
                f"Sim ID: {sim_name}; Degradation rate combo: {degradation_rate_combo}")

            # violin plot:
            sns.violinplot(x='generation', y='fmconc', data=p_df, inner="points",
                           color='lightgray', alpha=0.7, ax=ax2)

            sns.stripplot(x='generation',
                          y='fmconc',
                          data=p_df,
                          size=2,
                          hue='Time (s)',
                          palette='viridis',
                          alpha=0.2,
                          jitter=0.3,
                          #dodge=True # this makes the code take forever
                          ax=ax2)

            # Color mapping and normalization
            norm = Normalize(vmin=p_df['Time (s)'].min(), vmax=p_df['Time (s)'].max())
            sm = ScalarMappable(cmap='viridis', norm=norm)  # Create a ScalarMappable

            # Create a color bar
            cbar = fig.colorbar(sm, ax=ax2)
            cbar.set_label('Time (s)')

            # Plot specifics:
            ax2.set_xlabel('Generation')
            ax2.set_ylabel('Free Monomer Concentration (mmol/L)', size=10)
            ax2.legend().remove()

            plt.suptitle(
                f"Free monomer counts over time & concentration violin plot for {protein}",
                size=12)
            plt.tight_layout()

            # save the plot:
            file_name = plotOutFileName + "_and_free_monomer_counts_" + protein
            exportFigure(plt, plotOutDir, file_name, metadata)

        # Do the same for complexes:
        for complex in PLOT_COMPLEXES_revised:
            # Gather the relevant data:
            complex_idx = complex_idx_dict[complex]
            cconc = complex_concentrations[:, complex_idx]
            data = {'generation': gen_mask,
                    'cconc': cconc,
                    'all_time': time,
                    'Time (s)': gen_time_mask}
            c_df = pd.DataFrame(data)

            # Generate the plots:
            sns.violinplot(x='generation', y='cconc', data=c_df, inner="points",
                           color='lightgray', alpha=0.7)

            sns.stripplot(x='generation',
                          y='cconc',
                          data=c_df,
                          size=3,
                          hue='Time (s)',
                          palette='crest',
                          alpha=0.3,
                          jitter=0.3,
                          #dodge=True # this makes the code take forever
                          )

            # Color mapping and normalization
            norm = Normalize(vmin=c_df['Time (s)'].min(), vmax=c_df['Time (s)'].max())
            sm = ScalarMappable(cmap='crest', norm=norm)
            # Create a color bar
            cbar = plt.colorbar(sm)
            cbar.set_label('Time (s)')
            # Plot specifics:
            sim_name = metadata["description"]
            degradation_rate_combo = sim_data.protein_degradation_combo_selection
            plt.xlabel('Generation')
            plt.ylabel('Complex Concentration (mmol/L)')
            plt.title(
                f"Complex concentration for {complex}\n Sim ID: {sim_name}; \nDegradation rate combo: {degradation_rate_combo}",
                size=12)
            plt.tight_layout()

            # save the plot:
            file_name = plotOutFileName + "_" + complex
            exportFigure(plt, plotOutDir, file_name, metadata)


        # create a plot of the complex counts under it:
        for complex in PLOT_COMPLEXES_revised:
            # Gather the relevant data:
            complex_idx = complex_reader.readAttribute("complexIDs").index(complex)
            complex_counts_data = complex_counts[:, complex_idx]
            cconc = complex_concentrations[:, complex_idx]
            data = {'generation': gen_mask,
                    'cconc': cconc,
                    'all_time': time,
                    'Time (s)': gen_time_mask}
            c_df = pd.DataFrame(data)

            # Generate the plots:
            fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(6*n_gens,10))

            # Complex Counts plot
            ax1.plot(time, complex_counts_data, color='lightseagreen', label='Complex Counts')
            for i in range(len(end_generation_times)):
                dt = end_generation_times[i]
                ax1.axvline(x=dt, linestyle='--', color="orange")

            # Plot specifics:
            sim_name = metadata["description"]
            degradation_rate_combo = sim_data.protein_degradation_combo_selection
            ax1.legend(fontsize=5)
            ax1.set_ylabel("Complex Counts")
            ax1.set_xlabel("Time (s)")
            ax1.set_title(
                f"Sim ID: {sim_name}; Degradation rate combo: {degradation_rate_combo}")

            # violin plot:
            sns.violinplot(x='generation', y='cconc', data=c_df, inner="points",
                           color='lightgray', alpha=0.7, ax=ax2)

            sns.stripplot(x='generation',
                          y='cconc',
                          data=c_df,
                          size=2,
                          hue='Time (s)',
                          palette='viridis',
                          alpha=0.2,
                          jitter=0.3,
                          #dodge=True # this makes the code take forever
                          ax=ax2)

            # Color mapping and normalization
            norm = Normalize(vmin=c_df['Time (s)'].min(), vmax=c_df['Time (s)'].max())
            sm = ScalarMappable(cmap='viridis', norm=norm)
            # Create a color bar
            cbar = fig.colorbar(sm, ax=ax2)
            cbar.set_label('Time (s)')
            # Plot specifics:
            ax2.set_xlabel('Generation')
            ax2.set_ylabel('Complex Concentration (mmol/L)', size=10)
            ax2.legend().remove()

            avg_conc = np.mean(cconc)
            std_conc = np.std(cconc)

            ax2.axhline(y=avg_conc, linestyle='--', color="orange", label=f'Average Conc: {avg_conc:.2e} mmol/L')

            plt.suptitle(
                f"Complex counts over time & concentration violin plot for {complex} \n"
                   f"Avg Conc: {avg_conc:.2e} mmol/L; Std Dev: {std_conc:.2e} mmol/L",
                size=12)
            plt.tight_layout()
            # save the plot:
            file_name = plotOutFileName + "_and_complex_counts_" + complex
            exportFigure(plt, plotOutDir, file_name, metadata)



if __name__ == '__main__':
    Plot().cli()
