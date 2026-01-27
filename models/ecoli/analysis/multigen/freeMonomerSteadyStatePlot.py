"""
Template for multigen analysis plots
"""

import pickle
import os
from wholecell.utils import units

from matplotlib import pyplot as plt
# noinspection PyUnresolvedReferences
import numpy as np
import plotly.graph_objects as go
from models.ecoli.analysis import multigenAnalysisPlot
from wholecell.analysis.analysis_tools import (exportFigure,
    read_bulk_molecule_counts, read_stacked_bulk_molecules, read_stacked_columns)
from wholecell.io.tablereader import TableReader
from wholecell.io.tablereader import TableReader
import io
from wholecell.io import tsv
from wholecell.utils.filepath import ROOT_PATH

HIGHLIGHT_IN_RED = ["PD03867[c]", "EG11734-MONOMER[c]"]#['EG10863-MONOMER[c]','DETHIOBIOTIN-SYN-MONOMER[c]','DCUR-MONOMER[c]']
HIGHLIGHT_IN_BLUE = []#['CARBPSYN-SMALL[c]', 'CDPDIGLYSYN-MONOMER[i]','EG10743-MONOMER[c]','GLUTCYSLIG-MONOMER[c]']
HIGHLIGHT_IN_PURPLE = []#['G6890-MONOMER[c]','PD03938[c]','G6737-MONOMER[c]','RPOD-MONOMER[c]','PD02936[c]','RED-THIOREDOXIN2-MONOMER[c]']

# temporary flag to highlight proteins that make complexes:
PLOT_BY_COMPLEX_TYPE = True

# TODO: add average complexed counts and average total counts!


# function to match gene symbols to monomer ids
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
        with open(validationDataFile, 'rb') as f:
            validation_data = pickle.load(f)
        sim_id = metadata['description']
        cell_paths = self.ap.get_cells()
        sim_dir = cell_paths[0]
        simOutDir = os.path.join(sim_dir, 'simOut')

        # Extract protein indexes for each new gene
        monomer_counts_reader = TableReader(
            os.path.join(simOutDir, "MonomerCounts"))
        monomer_idx_dict = {monomer: i for i, monomer in
                            enumerate(monomer_counts_reader.readAttribute(
                                'monomerIds'))}

        # extract the monomer IDs of interest
        monomerIds = monomer_counts_reader.readAttribute(
            "monomerIds")  # this is the one that matches the indexing  I used earlier to construct the listeners!


        def get_protein_data(sim_data, remove_yibQ):
            protein_ids = sim_data.process.translation.monomer_data['id']
            #deg_rate_source = sim_data.process.translation.monomer_data['deg_rate_source']
            degradation_rates = sim_data.process.translation.monomer_data['deg_rate'].asNumber(
                1 / units.s)
            half_lives = np.log(2) / degradation_rates / 60  # in minutes
            if remove_yibQ:
                indices = [i for i, x in enumerate(protein_ids) if x == "EG12298-MONOMER[c]"]
                protein_ids_wo_yibQ = np.delete(protein_ids, indices[0])
                half_lives_wo_yibQ = np.delete(half_lives, indices[0])
                return protein_ids_wo_yibQ, half_lives_wo_yibQ,
            else:
                return protein_ids, half_lives

        protein_ids, half_lives = (
            get_protein_data(sim_data, remove_yibQ=False))

        # Get the common names for each protein:
        common_names = [get_common_name(protein_id) for protein_id in protein_ids]

        # Extract information on complex relationships:
        complexation_data = sim_data.process.complexation.molecules_to_all_downstream_complexes_dict
        equilibrium_data = sim_data.process.equilibrium.molecules_to_all_downstream_complexes_dict

        complexation_complex_info = [complexation_data.get(protein_id, "NA") for protein_id in
                                     protein_ids]
        equilibrium_complex_info = [equilibrium_data.get(protein_id, "NA") for protein_id in
                                    protein_ids]
        hi = 6



        # extract the data over all generations:
        # Load data
        time = read_stacked_columns(cell_paths, 'Main', 'time')
        free_monomer_counts =  read_stacked_columns(cell_paths, 'MonomerCounts', "freeMonomerCounts")

        # doubling time function from nora:
        def extract_doubling_times(cell_paths):
            # Load data
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

        # find how many proteins were removed via dilution for each doubling time:
        diluted_counts = np.zeros(((len(doubling_times) - 1), len(monomerIds)))
        diluted_counts_over_time = np.zeros(((len(time)), len(monomerIds)))
        for i in range(
                len(doubling_times) - 1):  # -1 is to account for not doing the last generation
            end_gen = end_generation_indices[i]
            start_gen = start_generation_indices[i + 1]  # the first start is zero, so skip that
            print(end_gen, start_gen)

            # find the protein counts at the end of the generation:
            monomer_counts_at_gen_end = free_monomer_counts[end_gen,:]  # get this for each protein

            # find the protein counts at the start of the next generation:
            monomer_counts_at_gen_start = free_monomer_counts[start_gen, :]

            # find the difference between the two:
            protein_counts_removed = monomer_counts_at_gen_end - monomer_counts_at_gen_start
            diluted_counts[i, :] = protein_counts_removed
            diluted_counts_over_time[start_gen,:] = protein_counts_removed  # put it at the start of the next gen in terms of time

        print(diluted_counts)


        # Compute how many proteins were removed via degradation over the entire sim length:
        degraded_counts = read_stacked_columns(cell_paths, 'MonomerCounts', "monomersDegraded")
        avg_degraded_counts = np.mean(degraded_counts, axis=0)

        # Compute how many proteins were removed via complexation over the entire sim length:
        complexed_counts = read_stacked_columns(cell_paths, 'ComplexationListener',
                                                "monomersComplexed")
        avg_complexed_counts = np.mean(complexed_counts, axis=0)

        # Compute how many proteins were removed via equilibrium complexation over the entire sim length:
        eq_complexed_counts = read_stacked_columns(cell_paths, 'EquilibriumListener',
                                                   "freeMonomersComplexed")
        avg_eq_complexed_counts = np.mean(eq_complexed_counts, axis=0)

        # Compute how many counts were added via elongation over the entire sim length:
        elongated_counts = read_stacked_columns(cell_paths, 'MonomerCounts',
                                                "monomersElongated")
        avg_elongated_counts = np.mean(elongated_counts, axis=0)


        # compute how many proteins were removed via dilution over the entire sim length:
        # todo: make sure this is the right way to compute the average (when you only have one dilution timepoint (becuase there was that one graph I did where it was different depending on the size of the array)
        total_diluted_counts = np.sum(diluted_counts, axis = 0)
        avg_diluted_counts = total_diluted_counts / len(time) # divide by the number of timesteps to get the average per timestep

        # compute the average loss rate for each protein:
        avg_loss_rate = avg_degraded_counts + avg_diluted_counts + avg_complexed_counts + avg_eq_complexed_counts
        log_avg_loss_rate = np.log10(avg_loss_rate) # todo: consider adding 1 to this to avoid log(0) and all


        # compute how many counts were added via elongation over the entire sim length:
        log_avg_production_rate = np.log10(avg_elongated_counts) # todo: consider adding 1 to this to avoid log(0) and all and being able to see all proteins

        hi = 5
        average_free_monomer_counts = np.mean(free_monomer_counts, axis=0)

        def hover_text(protein_ids, x, y, half_lives, common_names, complexation_complex_info,
                       equilibrium_complex_info, average_free_monomer_counts):
            texts = []

            def complex_info_formatter(complexation_info_entry):
                if isinstance(complexation_info_entry, dict) and complexation_info_entry:
                    complexation_entries = []
                    for comp_id, comp_data in complexation_info_entry.items():
                        if isinstance(comp_data, list) and comp_data:
                            entry_info = {}
                            for entry in comp_data:
                                entry_info.update(entry)

                            # Extract the necessary values
                            reaction_id = entry_info.get('reaction_id', 'N/A')
                            stoichiometry = entry_info.get('stoichiometry', 'N/A')
                            stoich_unknown = entry_info.get('stoich_unknown', 'N/A')
                            complex_type = entry_info.get('complex_type', 'N/A')

                            # Format the info:
                            formatted_info = f"{comp_id}, {reaction_id}, {stoichiometry}, {stoich_unknown}, {complex_type}"
                            complexation_entries.append(formatted_info)

                    # Join the entries into a single string
                    return ';<br>'.join(
                        complexation_entries) if complexation_entries else 'No Complexation Info'
                else:
                    return 'No Complexation Info'

            for i in range(len(protein_ids)):
                # Extracting complexation and equilibrium info:
                complexation_info_str = complex_info_formatter(complexation_complex_info[i])
                equilibrium_info_str = complex_info_formatter(equilibrium_complex_info[i])

                texts.append((
                    f"{protein_ids[i]}<br>"
                    f"Production Rate: {x[i]:.2f}<br>"
                    f"Loss Rate: {y[i]:.2f}<br>"
                    f"Average Free Monomer Count: {average_free_monomer_counts[i]:.2f}<br>"
                    f"Half Life (mins): {half_lives[i]}<br>"
                    f"Common Name: {common_names[i]}<br>"
                    f"<span style='font-size: 10px;'>Complexation Info (ID, rxn, stoich, stoich unknown, type): {complexation_info_str}<br>"
                    f"<span style='font-size: 10px;'>Equilibrium Info (ID, rxn, stoich, stoich unknown, type): {equilibrium_info_str}</span>"
                ))

            return texts

        # plot the loss rate vs the production rate:
        # Create figure
        fig = go.Figure()
        hover_info = hover_text(protein_ids, log_avg_production_rate, log_avg_loss_rate,
                                half_lives, common_names, complexation_complex_info,
                                equilibrium_complex_info, average_free_monomer_counts)

        # Scatter plot for all proteins (grey)
        fig.add_trace(go.Scatter(
            x=log_avg_production_rate,
            y=log_avg_loss_rate,
            mode='markers',
            hovertext=hover_info,
            marker=dict(size=5, color='lightseagreen', opacity=0.3),
            name='All Proteins'
        ))

        if len(HIGHLIGHT_IN_RED) > 0:
            # Scatter plot for red proteins
            red_protein_indices = [monomer_idx_dict[protein] for protein in HIGHLIGHT_IN_RED]
            fig.add_trace(go.Scatter(
                x=log_avg_production_rate[red_protein_indices],
                y=log_avg_loss_rate[red_protein_indices],
                mode='markers',
                hovertext=HIGHLIGHT_IN_RED,
                marker=dict(size=5, color='red', opacity=1),
                name='Red Proteins'
        ))
        if len(HIGHLIGHT_IN_BLUE) > 0:
            # Scatter plot for blue proteins
            blue_protein_indices = [monomer_idx_dict[protein] for protein in HIGHLIGHT_IN_BLUE]
            fig.add_trace(go.Scatter(
                x=log_avg_production_rate[blue_protein_indices],
                y=log_avg_loss_rate[blue_protein_indices],
                mode='markers',
                hovertext=HIGHLIGHT_IN_BLUE,
                marker=dict(size=5, color='blue', opacity=1),
                name='Blue Proteins'
            ))

        if len(HIGHLIGHT_IN_PURPLE) > 0:
            # Scatter plot for purple proteins
            purple_protein_indices = [monomer_idx_dict[protein] for protein in HIGHLIGHT_IN_PURPLE]
            fig.add_trace(go.Scatter(
                x=log_avg_production_rate[purple_protein_indices],
                y=log_avg_loss_rate[purple_protein_indices],
                mode='markers',
                hovertext=HIGHLIGHT_IN_PURPLE,
                marker=dict(size=5, color='hotpink', opacity=1),
                name='Purple Proteins'
            ))
        # Scater plot for purple proteins:

        # Generate line data
        x = np.linspace(-4, 3, 100)

        # y = x line (black)
        fig.add_trace(go.Scatter(
            x=x, y=x, mode='lines', line=dict(color='black', width=2,),
            name='y = x'
        ))

        # y = x + 1 line (green, dashed)
        fig.add_trace(go.Scatter(
            x=x, y=x + 1, mode='lines',
            line=dict(color='green', width=2, dash='dash'), name='y = x + 1'
        ))

        # y = x - 1 line (green, dashed)
        fig.add_trace(go.Scatter(
            x=x, y=x - 1, mode='lines',
            line=dict(color='green', width=2, dash='dash'), name='y = x - 1'
        ))

        # Layout settings
        fig.update_layout(
            title=f"Log10 Average Loss Rate vs Log10 Average Production Rate ({sim_id})",
            xaxis_title="Log10 Average Production Rate",
            yaxis_title="Log10 Average Loss Rate",
            width=700, height=700,
            showlegend=True,
        )

        # Save the plot:
        plot_name = plotOutFileName +"_"+ sim_id +"_red_" + str(HIGHLIGHT_IN_RED) + "_blue_" + str(HIGHLIGHT_IN_BLUE) + "_purple_" + str(HIGHLIGHT_IN_PURPLE) + ".html"
        fig.write_html(os.path.join(plotOutDir, plot_name))

        # Function to add lines:
        def add_lines(fig):
            # Generate line data
            x = np.linspace(-4, 3, 100)

            # y = x line (black)
            fig.add_trace(go.Scatter(
                x=x, y=x, mode='lines', line=dict(color='black', width=2, ),
                name='y = x'
            ))

            # y = x + 1 line (green, dashed)
            fig.add_trace(go.Scatter(
                x=x, y=x + 1, mode='lines',
                line=dict(color='green', width=2, dash='dash'), name='y = x + 1'
            ))

            # y = x - 1 line (green, dashed)
            fig.add_trace(go.Scatter(
                x=x, y=x - 1, mode='lines',
                line=dict(color='green', width=2, dash='dash'), name='y = x - 1'
            ))

        # Plot other plots if specified:
        # TODO: add two-component and ribosome and replisome complexation types
        if PLOT_BY_COMPLEX_TYPE:
            # Find all proteins that make complexes:
            molecules_to_all_downstream_complexes = sim_data.process.complexation.molecules_to_all_downstream_complexes_dict
            monomers_in_complexes = list(molecules_to_all_downstream_complexes.keys())
            molecules_to_all_downstream_eq_complexes = sim_data.process.equilibrium.molecules_to_all_downstream_complexes_dict
            monomers_in_eq_complexes = list(molecules_to_all_downstream_eq_complexes.keys())

            # Find all proteins not in complexes:
            monomers_not_in_complexes = [monomer for monomer in monomerIds if
                                        monomer not in monomers_in_complexes and
                                        monomer not in monomers_in_eq_complexes]

            # Create separate figure for complexation monomers
            fig_complexation = go.Figure()
            complexation_monomer_indices = [monomer_idx_dict[monomer] for monomer in
                                            monomers_in_complexes if
                                            monomer in monomerIds]
            # Accessing attributes with proper indexing
            if len(complexation_monomer_indices) > 0:
                hover_info_complexation = hover_text(
                    [monomerIds[i] for i in complexation_monomer_indices],
                    log_avg_production_rate[complexation_monomer_indices],
                    log_avg_loss_rate[complexation_monomer_indices],
                    half_lives[complexation_monomer_indices],
                    [common_names[i] for i in complexation_monomer_indices],
                    [complexation_complex_info[i] for i in complexation_monomer_indices],
                    [equilibrium_complex_info[i] for i in complexation_monomer_indices],
                    average_free_monomer_counts[complexation_monomer_indices]
                )

                fig_complexation.add_trace(go.Scatter(
                    x=log_avg_production_rate[complexation_monomer_indices],
                    y=log_avg_loss_rate[complexation_monomer_indices],
                    mode='markers',
                    hovertext=hover_info_complexation,
                    marker=dict(size=4, color='orange', opacity=.5),
                    name=f'Monomers in Complexation Complexes (n={len(complexation_monomer_indices)})'
                ))

                add_lines(fig_complexation)

            fig_complexation.update_layout(
                title=f"Free Monomer Average Loss Rate vs Average Production Rate<br>Sim ID: {sim_id}",
                xaxis_title="Log10 Average Production Rate",
                yaxis_title="Log10 Average Loss Rate",
                width=700, height=700
            )

            # Save the complexation plot
            complexation_plot_name = "complexation_monomers_" + sim_id + ".html"
            fig_complexation.write_html(os.path.join(plotOutDir, complexation_plot_name))

            # Create separate figure for equilibrium monomers
            fig_equilibrium = go.Figure()
            equilibrium_monomer_indices = [monomer_idx_dict[monomer] for monomer in
                                           monomers_in_eq_complexes if
                                           monomer in monomerIds]

            if len(equilibrium_monomer_indices) > 0:
                hover_info_equilibrium = hover_text(
                    [monomerIds[i] for i in equilibrium_monomer_indices],
                    log_avg_production_rate[equilibrium_monomer_indices],
                    log_avg_loss_rate[equilibrium_monomer_indices],
                    half_lives[equilibrium_monomer_indices],
                    [common_names[i] for i in equilibrium_monomer_indices],
                    [complexation_complex_info[i] for i in equilibrium_monomer_indices],
                    [equilibrium_complex_info[i] for i in equilibrium_monomer_indices],
                    average_free_monomer_counts[equilibrium_monomer_indices]
                )

                fig_equilibrium.add_trace(go.Scatter(
                    x=log_avg_production_rate[equilibrium_monomer_indices],
                    y=log_avg_loss_rate[equilibrium_monomer_indices],
                    mode='markers',
                    hovertext=hover_info_equilibrium,
                    marker=dict(size=5, color='lightseagreen', opacity=.6),
                    name=f'Monomers in Equilibrium Complexes (n={len(equilibrium_monomer_indices)})'
                ))

                add_lines(fig_equilibrium)

            fig_equilibrium.update_layout(
                title=f"Free Monomer Average Loss Rate vs Average Production Rate<br>Sim ID: {sim_id}",
                xaxis_title="Log10 Average Production Rate",
                yaxis_title="Log10 Average Loss Rate",
                width=700, height=700
            )

            # Save the equilibrium plot
            equilibrium_plot_name = "equilibrium_monomers_" + sim_id + ".html"
            fig_equilibrium.write_html(os.path.join(plotOutDir, equilibrium_plot_name))

            # Create separate figure for the monomer that makes no complexes
            fig_no_complex = go.Figure()

            if len(monomers_not_in_complexes) > 0:
                no_complex_monomer_indices = [monomer_idx_dict[monomer] for monomer in
                                              monomers_not_in_complexes if
                                              monomer in monomerIds]

                hover_info_no_complex = hover_text(
                    [monomerIds[i] for i in no_complex_monomer_indices],
                    log_avg_production_rate[no_complex_monomer_indices],
                    log_avg_loss_rate[no_complex_monomer_indices],
                    half_lives[no_complex_monomer_indices],
                    [common_names[i] for i in no_complex_monomer_indices],
                    [complexation_complex_info[i] for i in no_complex_monomer_indices],
                    [equilibrium_complex_info[i] for i in no_complex_monomer_indices],
                    average_free_monomer_counts[no_complex_monomer_indices]
                )

                fig_no_complex.add_trace(go.Scatter(
                    x=log_avg_production_rate[no_complex_monomer_indices],
                    y=log_avg_loss_rate[no_complex_monomer_indices],
                    mode='markers',
                    hovertext=hover_info_no_complex,
                    marker=dict(size=4, color='salmon', opacity=.6),
                    name=f'monomers in zero complexes (n={len(no_complex_monomer_indices)})'
                ))

                add_lines(fig_no_complex)

            fig_no_complex.update_layout(
                title=f"Free Monomer Average Loss Rate vs Average Production Rate<br>Sim ID: {sim_id}",
                xaxis_title="Log10 Average Production Rate",
                yaxis_title="Log10 Average Loss Rate",
                width=700, height=700
            )

            # Save the no complex plot
            no_complex_plot_name = "no_complex_monomer_" + sim_id + ".html"
            fig_no_complex.write_html(os.path.join(plotOutDir, no_complex_plot_name))

            # Create a combined figure for all three categories
            fig_combined = go.Figure()

            # Plot monomers that make no complexes
            fig_combined.add_trace(go.Scatter(
                x=log_avg_production_rate[no_complex_monomer_indices],
                y=log_avg_loss_rate[no_complex_monomer_indices],
                mode='markers',
                hovertext=hover_info_no_complex,
                marker=dict(size=3, color='salmon', opacity=.5),
                name=f'Monomers never in Complexes (n={len(no_complex_monomer_indices)})'
            ))

            # Monomers in complexation complexes
            fig_combined.add_trace(go.Scatter(
                x=log_avg_production_rate[complexation_monomer_indices],
                y=log_avg_loss_rate[complexation_monomer_indices],
                mode='markers',
                hovertext=hover_info_complexation,
                marker=dict(size=3, color='orange', opacity=0.4),
                name=f'Monomers in Complexation Complexes (n={len(complexation_monomer_indices)})'
            ))

            # Monomers in equilibrium complexes
            fig_combined.add_trace(go.Scatter(
                x=log_avg_production_rate[equilibrium_monomer_indices],
                y=log_avg_loss_rate[equilibrium_monomer_indices],
                mode='markers',
                hovertext=hover_info_equilibrium,
                marker=dict(size=5, color='green', opacity=0.6),
                name=f'Monomers in Equilibrium Complexes (n={len(equilibrium_monomer_indices)})'
            ))


            add_lines(fig_combined)

            # TODO: fix where the legend falls on the plot!

            fig_combined.update_layout(
                title=f"Free Monomer Average Loss Rate vs Average Production Rate<br>Sim ID: {sim_id}",
                xaxis_title="Log10 Average Production Rate",
                yaxis_title="Log10 Average Loss Rate",
                width=700, height=700,
                margin=dict(l=40, r=40, t=50, b=40),
                legend=dict(
                    x=.01,
                    y=0.99,
                    traceorder='normal',
                    font=dict(size=12),
                ),
                title_font=dict(size=16),
                xaxis=dict(title_font=dict(size=14)),
                yaxis=dict(title_font=dict(size=14)),
            )

            # Lock the aspect ratio to keep it square
            fig.update_yaxes(scaleanchor="x", scaleratio=1)


            # Save the combined plot
            combined_plot_name = "combined_monomers_" + sim_id + ".html"
            fig_combined.write_html(os.path.join(plotOutDir, combined_plot_name))


if __name__ == '__main__':
    Plot().cli()
