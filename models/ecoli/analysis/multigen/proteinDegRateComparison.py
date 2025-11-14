"""
Template for multigen analysis plots
"""

import pickle
import os

# noinspection PyUnresolvedReferences
import math
import numpy as np
import pandas as pd
import plotly.graph_objs as go
from plotly.subplots import make_subplots

from models.ecoli.analysis import multigenAnalysisPlot
from models.ecoli.analysis.multigen.violin_concentration_plots_for_complexes import PLOT_COMPLEXES
from wholecell.analysis.analysis_tools import (exportFigure,
                                               read_bulk_molecule_counts, read_stacked_bulk_molecules, read_stacked_columns)
from wholecell.io.tablereader import TableReader

PLOT_PROTEINS = ["G6890-MONOMER[c]",
                 "PD03938[c]",
                 "G6737-MONOMER[c]",
                 "RPOD-MONOMER[c]",
                 "PD02936[c]",
                 "RED-THIOREDOXIN2-MONOMER[c]"]


ONE_PDF = True

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

        # Extract generation and doubling time information:
        (time, doubling_times, end_generation_times,
         start_generation_indices, end_generation_indices) = extract_doubling_times(cell_paths)
        n_gens = metadata["generations"]
        gen_mask = make_generation_mask(cell_paths, n_gens)

        # Load data - degradation information stored in the listener
        rawDegradationRate = read_stacked_columns(cell_paths, 'MonomerCounts', "rawDegradationRate")
        rawDegnProtein = read_stacked_columns(cell_paths, 'MonomerCounts', "rawDegradationProtein")
        activeDegradationRate = read_stacked_columns(cell_paths, 'MonomerCounts', "activeDegradationRate")
        activeDegnProtein = read_stacked_columns(cell_paths, 'MonomerCounts', "activeDegnProtein")

        # Extract the complex IDs:
        complex_reader = TableReader(os.path.join(simOutDir, "ComplexationListener"))
        complex_idx_dict = {complex: i for i, complex in
                            enumerate(complex_reader.readAttribute(
                                'complexIDs'))}
        complexIDs = complex_reader.readAttribute("complexIDs")
        (complex_counts,) = read_stacked_bulk_molecules(
            cell_paths, complexIDs)

        # read in the counts to molar conversion (at each time point):
        ectm = read_stacked_columns(cell_paths, 'EnzymeKinetics', "countsToMolar")

        # calculate the complex concentrations:
        complex_concentrations = complex_counts * ectm # mmol/L

        monomer_counts_reader = TableReader( os.path.join(simOutDir, "MonomerCounts") )
        monomer_idx_dict = {monomer: i for i, monomer in
                            enumerate(monomer_counts_reader.readAttribute(
                                'monomerIds'))}

        # Make sure the proteins inputted are valid and have a compartment tag:
        PLOT_PROTEINS_revised = check_validity_and_get_compartment(PLOT_PROTEINS)
        PLOT_COMPLEXES_revised = check_validity_and_get_compartment(PLOT_COMPLEXES)
        protein_idx = [monomer_idx_dict[protein] for protein in PLOT_PROTEINS_revised]

        # PLOTTING:
        # --- subplot layout: max 3 columns ---
        n_proteins = len(PLOT_PROTEINS_revised)
        n_cols = min(3, n_proteins)
        n_rows = math.ceil(n_proteins / n_cols)

        fig = make_subplots(
            rows=n_rows,
            cols=n_cols,
            shared_yaxes=False,
            y_title = "Dynamic Degradation rate (1/sec)",
            subplot_titles=PLOT_PROTEINS_revised,
            horizontal_spacing=0.1,
            vertical_spacing=0.12,
        )

        # get the two generation values and assign colors
        gen_values = np.unique(gen_mask)
        gen_colors = {
            gen_values[0]: "rgba(190, 225, 201, 0.7)",  #mint
            gen_values[1]: "rgba(85, 133, 125, 0.7)",  # dark mint
        }

        for i, protein_i in enumerate(protein_idx):
            row = i // n_cols + 1
            col = i % n_cols + 1

            POI_rawDegradationRate = rawDegradationRate[:, protein_i]  # POI==protein of interest
            POI_activeDegradationRate = activeDegradationRate[:, protein_i]

            POI_rawDegradationRate = np.trim_zeros(np.unique(POI_rawDegradationRate))

            data = {
                'Dynamic Degradation Rate': POI_activeDegradationRate,
                'generation': gen_mask,
                'all_time': time,
            }
            df_data = pd.DataFrame(data)
            df_filtered = df_data[df_data['Dynamic Degradation Rate'] != 0]

            for gen in gen_values:
                fig.add_trace(
                    go.Violin(
                        x = df_filtered['generation'][df_filtered['generation'] == gen],
                        y = df_filtered['Dynamic Degradation Rate'][df_filtered['generation'] == gen],
                        name=gen,
                        x0 = PLOT_PROTEINS_revised[i],
                        box_visible=True,  # show box (median & IQR)
                        meanline_visible=True,  # show mean line
                        width = 0.8,
                        jitter=0.35,  # spread points
                        scalemode="count",  # violin width ~ number of points
                        marker=dict(size=3),
                        line=dict(color=gen_colors[gen]),
                        fillcolor=gen_colors[gen],
                        showlegend=False,
                    ),
                    row = row, col=col
                )

            # Horizontal reference line across that subplot
            fig.add_hline(
                y = POI_rawDegradationRate[0],
                line_dash="dash",
                line_width=2.5,
                line_color="#00798c",
                row=row, col=col
            )

        # One entry for legend of dotted hline
        fig.add_trace(
            go.Scatter(
                x=[None], y=[None],  # nothing actually plotted
                mode="lines",
                line=dict(dash="dash", width=2.5, color="#00798c"),
                name="Raw degradation rate"  # <-- legend text
            )
        )

        fig.update_layout(
            height=n_rows*300,
            width=n_cols*400,
            margin=dict(l=60, r=20, t=60, b=60),
            template="simple_white",
            title="Dynamic Protein Degradation Rate Compared Static Raw Rate",
        )
        ### Create Plot ###
        print(f'Saved output to {plotOutDir}')
        fig.write_html(f'{plotOutDir}/{plotOutFileName}.html')
        fig.update_layout(
            plot_bgcolor='rgba(0, 0, 0, 0)',  # Transparent plot area background
            paper_bgcolor='rgba(0, 0, 0, 0)'
        )
        fig.write_image(f'{plotOutDir}/{plotOutFileName}.png')


if __name__ == '__main__':
    Plot().cli()
