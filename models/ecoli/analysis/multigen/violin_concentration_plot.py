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

        hi = 5
        ezk_reader = TableReader(
            os.path.join(simOutDir, "EnzymeKinetics"))

        ectm = read_stacked_columns(cell_paths, 'EnzymeKinetics', "countsToMolar")
        ctm = ezk_reader.readColumn('countsToMolar')

        # calculate the free monomer concentrations:
        free_monomer_concentrations = fmcs * ectm # mmol/L
        time = read_stacked_columns(cell_paths, 'Main', 'time').squeeze()












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
        n_gens = metadata["generations"]
        # find the maximum generation length for plotting purposes:
        max_gen_length = np.max(doubling_times)

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
                                                                                                      end_idx])] - time[                                                                                  start_idx]

            return gen_time_mask

        # create the generation name mask and generation time mask:
        gen_mask = make_generation_mask(cell_paths, n_gens)
        gen_time_mask = generation_time_mask(cell_paths, n_gens)



        # plot the loss rate and the production rate:
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
            fig = plt.figure(figsize=(5*n_gens,6))

            sns.violinplot(x='generation', y='fmconc', data=p_df, inner="points",
                           color='lightgray', alpha=0.7)
            # Overlay the individual points
            palette = sns.color_palette("husl", len(p_df['Time (s)'].unique()))

            sns.stripplot(x='generation',
                          y='fmconc',
                          data=p_df,
                          size=3,
                          hue='Time (s)',  # Color of the points
                          alpha=0.3,  # Transparency of the points
                          jitter=0.3,
                          )  # Add jittering to the points for better visibility

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


if __name__ == '__main__':
    Plot().cli()
