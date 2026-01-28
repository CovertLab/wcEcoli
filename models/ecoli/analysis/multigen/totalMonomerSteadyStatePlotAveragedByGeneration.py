"""
Template for multigen analysis plots
"""

import pickle
import os
from wholecell.utils import units
from sklearn.metrics import r2_score
from scipy.stats import pearsonr
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

HIGHLIGHT_IN_RED = []#['EG10863-MONOMER[c]','DETHIOBIOTIN-SYN-MONOMER[c]','DCUR-MONOMER[c]']
HIGHLIGHT_IN_BLUE = []#['CARBPSYN-SMALL[c]', 'CDPDIGLYSYN-MONOMER[i]','EG10743-MONOMER[c]','GLUTCYSLIG-MONOMER[c]']
HIGHLIGHT_IN_PURPLE = []#['G6890-MONOMER[c]','PD03938[c]','G6737-MONOMER[c]','RPOD-MONOMER[c]','PD02936[c]','RED-THIOREDOXIN2-MONOMER[c]']

# temporary flag to highlight proteins that make complexes:
PLOT_BY_COMPLEX_TYPE = True

# Indicate if generations at the beginning should be skipped:
SKIP_INITIAL_GENERATIONS = 0

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


def get_plot_stats(x, y):
    # filter out nan or inf values:
    mask = np.isfinite(x) & np.isfinite(y)
    x = x[mask]
    y = y[mask]
    n = len(x)
    # Calculate pearson R:
    pearson_r, pearson_p = pearsonr(x, y)
    # Calculate pearson R^2:
    pearson_r2 = pearson_r ** 2
    # Calculate coefficent of determination R^2:
    r2 = r2_score(x, y)

    return pearson_r, pearson_r2, r2, n, mask


def add_stats_to_plot(x, y, fig):
    pearson_r, pearson_r2, r2, n, mask = get_plot_stats(x, y)
    # Put the data on the plot:
    plt.text(0.96, 0.03,
             f'Pearson R: {round(pearson_r, 3)}\n Pearson $R^2$: '
             f'{round(pearson_r2, 3)}\n'
             f'Coefficent of Determination $R^2$: {round(r2, 3)}',
             ha='right', va='bottom',
             transform=plt.gca().transAxes,
             fontsize=6, color='gray')

    text = (
        f'Stats for n={n} proteins (filtered out NaNs & infs):<br>'
        f'Pearson R: {round(pearson_r, 3)}; Pearson R<sup>2</sup>: {round(pearson_r2, 3)}<br>'
        f'Coefficient of determination R<sup>2</sup>: {round(r2, 3)}'
    )

    # Get the maximum x and minimum y to position the text in the bottom-right
    x_max = x[mask].max()
    y_min = y[mask].min()

    # Adjust x_max and y_min slightly outside the actual graph boundaries:
    text_offset_x = 0.05
    text_offset_y = 0.75

    # Adding text annotation to the bottom right
    fig.add_annotation(
        x=x_max + text_offset_x,
        y=y_min - text_offset_y,
        text=text,
        showarrow=False,
        bgcolor='rgba(255, 255, 255, 0.8)',
        bordercolor='rgba(0, 0, 0, 0.5)',
        borderwidth=1,
        borderpad=4,
        align='right',
        font=dict(size=8, color='gray'),
        xref='x',
        yref='y',
    )


class Plot(multigenAnalysisPlot.MultigenAnalysisPlot):


    def do_plot(self, seedOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
        with open(simDataFile, 'rb') as f:
            sim_data = pickle.load(f)

        sim_id = metadata['description']
        generations = metadata.get('generations', None)
        all_cells = self.ap.get_cells(generation=range(SKIP_INITIAL_GENERATIONS, generations))

        # Get all cell paths excluding initial generations if specified:
        cell_paths = self.ap.get_cells(generation=range(SKIP_INITIAL_GENERATIONS,len(all_cells)-1))
        # Note: dilution is calculated over the range of the proceeding generation and data is not
        # recorded for a theoretical "first time step of the n + 1 generation", there is no way to
        # definitively get the dilution for the last generation since there is no next generation
        # to compare to, so we just exclude that last generation in the analysis.
        sim_dir = cell_paths[0]
        simOutDir = os.path.join(sim_dir, 'simOut')

        # Extract protein indexes for each new gene
        monomer_counts_reader = TableReader(
            os.path.join(simOutDir, "MonomerCounts"))
        monomer_idx_dict = {monomer: i for i, monomer in
                            enumerate(monomer_counts_reader.readAttribute(
                                'monomerIds'))}

        # Extract the monomer IDs within the simulation:
        monomerIds = monomer_counts_reader.readAttribute(
            "monomerIds")  # this is the one that matches the indexing  I used earlier to construct the listeners!


        def get_protein_data(sim_data, remove_yibQ):
            protein_ids = sim_data.process.translation.monomer_data['id']
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

        # Load in the time spanned by the cells being analyzed:
        time = read_stacked_columns(cell_paths, 'Main', 'time').squeeze()

        # Load in the count data:
        total_monomer_counts = read_stacked_columns(cell_paths, 'MonomerCounts', "monomerCounts")
        free_monomer_counts =  read_stacked_columns(cell_paths, 'MonomerCounts', "freeMonomerCounts")
        complexed_monomer_counts = read_stacked_columns(cell_paths, "ComplexationListener", "complexedMonomerCounts")
        eq_complexed_monomer_counts = read_stacked_columns(cell_paths, "EquilibriumListener", "complexedMonomerCounts")

        # Compute how many proteins were removed via degradation over the entire sim length:
        degraded_counts = read_stacked_columns(cell_paths, 'MonomerCounts', "monomersDegraded")
        avg_degraded_counts = np.mean(degraded_counts, axis=0)

        # Compute how many counts were added via elongation over the entire sim length:
        elongated_counts = read_stacked_columns(cell_paths, 'MonomerCounts',
                                                "monomersElongated")
        avg_elongated_counts = np.mean(elongated_counts, axis=0)

        # Doubling time extraction function from Nora (handles the fact that the doubling times are misaligned)
        def extract_doubling_times(all_cells):
            # Load time data for all cells (including the last generation (that is not included in the plot)
            time = read_stacked_columns(all_cells, 'Main', 'time').squeeze()

            # Determine doubling time
            doubling_times = read_stacked_columns(all_cells, 'Main', 'time',
                                                  fun=lambda x: (x[-1] - x[0])).squeeze().astype(
                int)
            end_generation_times = np.cumsum(doubling_times) + time[0]
            start_generation_indices = np.searchsorted(time, end_generation_times[:-1],
                                                       side='left').astype(int)
            start_generation_indices = np.insert(start_generation_indices, 0, 0) + np.arange(
                len(doubling_times))
            end_generation_indices = start_generation_indices + doubling_times

            return time, doubling_times, end_generation_times, start_generation_indices, end_generation_indices

        # Extract the corrected doubling times:
        all_time, doubling_times, end_generation_times, start_generation_indices, end_generation_indices = extract_doubling_times(
            all_cells)

        # Find how many proteins were removed via dilution for each doubling time:
        # NOTE: to distribute dilution properly for calculating average loss rate,
        # the number of proteins diluted between generations is divided over the
        # timepoints spanned by the first gen and attributed to the first generation.
        # Since there is no way to calculate the number of proteins diluted at the
        # end of the last generation, that generation is not included in the averages.
        diluted_counts_per_generation = np.zeros(((len(doubling_times) - 1), len(monomerIds))) # make an array that is num generations (subtracting out the last generation) x num proteins
        diluted_counts = np.zeros(((len(time)), len(monomerIds))) # make an array that is the timepoints x num proteins
        diluted_counts_per_time_step = np.zeros(((len(time)), len(monomerIds))) # make an array that is the timepoints x num proteins
        for i in range(
                len(doubling_times) - 1):  # -1 is to account for not doing the last generation
            end_gen = end_generation_indices[i] # End time of the CURRENT generation
            start_gen = start_generation_indices[i + 1]  # Starts time of the NEXT generation
            print(end_gen, start_gen)

            # momentarily redefine total_monomer_counts to be over all cells:
            total_monomer_counts_temp = read_stacked_columns(all_cells, 'MonomerCounts', "monomerCounts")

            # Find the number of protein counts at the end of the generation:
            monomer_counts_at_gen_end = total_monomer_counts_temp[end_gen,:]

            # Find the protein counts at the start of the next generation:
            monomer_counts_at_gen_start = total_monomer_counts_temp[start_gen, :] # Note: since the loop goes to doubling_times - 1, this will never index out of bounds

            # Find the total change in protein counts between the generations:
            protein_count_change = monomer_counts_at_gen_end - monomer_counts_at_gen_start

            # Add the degraded counts and subtract the elongated counts at the doubling time to isolate dilution:
            protein_counts_removed = protein_count_change + degraded_counts[end_gen, :] - elongated_counts[end_gen, :]

            # Store the diluted counts for the generation:
            diluted_counts_per_generation[i, :] = protein_counts_removed

            # Store where the counts are actually diluted:
            diluted_counts[end_gen,:] = protein_counts_removed  # put it at the end of the generation timepoint

            # Distribute the diluted counts over the time steps spanned by the generation:
            diluted_counts_per_time_step[start_generation_indices[i]:end_gen, :] = (protein_counts_removed / doubling_times[i])  # distribute the diluted counts over the time steps in the generation

        # Calculate the average loss rate and production rate for each generation:
        average_loss_rate_per_generation = np.zeros_like(diluted_counts_per_generation)
        average_production_rate_per_generation = np.zeros_like(diluted_counts_per_generation)
        for i in range(len(doubling_times) - 1):  # -1 is to account for not doing the last generation
            # Initialize the start and end time for the generation:
            start_gen = start_generation_indices[i]
            end_gen = end_generation_indices[i]
            generation_time_span = end_gen - start_gen

            # Determine the slices for the generation:
            degraded_counts_gen = degraded_counts[start_gen:end_gen, :]
            elongated_counts_gen = elongated_counts[start_gen:end_gen, :]
            diluted_counts_gen = diluted_counts_per_time_step[start_gen:end_gen, :]

            # Compute the average rates for the generation:
            average_loss_rate_per_generation[i, :] = (np.sum(degraded_counts_gen, axis=0) + np.sum(diluted_counts_gen, axis=0)) / generation_time_span
            average_production_rate_per_generation[i, :] = np.sum(elongated_counts_gen, axis=0) / generation_time_span

        # Compute the overall average rates across all generations:
        avg_loss_rate = np.mean(average_loss_rate_per_generation, axis=0)
        avg_production_rate = np.mean(average_production_rate_per_generation, axis=0)

        # Take the log10 of the average rates for plotting:
        log_avg_loss_rate = np.log10(avg_loss_rate)  # todo: consider adding 1 to this to avoid log(0) and all
        log_avg_production_rate = np.log10(avg_production_rate)  # todo: consider adding 1 to this to avoid log(0) and all

        # Compute the average diluted and degraded counts over time (for the plot hover text):
        avg_diluted_counts = np.mean(diluted_counts, axis=0)
        avg_degraded_counts = np.mean(degraded_counts, axis=0)

        # Compute the average of all the monomer counts over time (for the plot hover text):
        average_total_monomer_counts = np.mean(total_monomer_counts, axis=0)
        average_free_monomer_counts = np.mean(free_monomer_counts, axis=0)
        average_complexed_monomer_counts = np.mean(complexed_monomer_counts, axis=0)
        average_eq_complexed_monomer_counts = np.mean(eq_complexed_monomer_counts, axis=0)

        # Indicate the plot title and name:
        plot_title = 'Average Total Monomer Loss vs Production Rate (averaged by generation)'

        def hover_text(protein_ids, x, y, half_lives, common_names, complexation_complex_info,
                       equilibrium_complex_info, average_total_monomer_counts,
                       average_free_monomer_counts,
                       average_complexed_monomer_counts, average_eq_complexed_monomer_counts,
                       avg_diluted_counts, avg_degraded_counts):
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
                # Extracting complexation and equilibrium complex info:
                complexation_info_str = complex_info_formatter(complexation_complex_info[i])
                equilibrium_info_str = complex_info_formatter(equilibrium_complex_info[i])

                texts.append((
                    f"{protein_ids[i]}<br>"
                    f"Half Life (mins): {half_lives[i]}<br>"
                    f"Common Name: {common_names[i]}<br>"
                    f"Average Production Rate: {x[i]:.2f}<br>"
                    f"Average Loss Rate: {y[i]:.2f} (dilution: {np.log10(avg_diluted_counts[i])}; degradation: {avg_degraded_counts[i]})<br>"
                    f"Average Total Monomer Count: {average_total_monomer_counts[i]:.2f}<br>"
                    f"Average Free Monomer Count: {average_free_monomer_counts[i]:.2f}<br>"
                    f"Average Complexed Monomer Count: {average_complexed_monomer_counts[i]:.2f}<br>"
                    f"Average Equilibrium Complexed Monomer Count: {average_eq_complexed_monomer_counts[i]:.2f}<br>"
                    f"<span style='font-size: 10px;'>Complexation Info (ID, rxn, stoich, stoich unknown, type):<br> {complexation_info_str}<br>"
                    f"<span style='font-size: 10px;'>Equilibrium Info (ID, rxn, stoich, stoich unknown, type):<br> {equilibrium_info_str}</span>"
                ))

            return texts

        # Create figure
        fig = go.Figure()
        hover_info = hover_text(protein_ids, log_avg_production_rate, log_avg_loss_rate,
                                half_lives, common_names, complexation_complex_info,
                                equilibrium_complex_info, average_total_monomer_counts,
                                average_free_monomer_counts,
                                average_complexed_monomer_counts,
                                average_eq_complexed_monomer_counts,
                                avg_diluted_counts, avg_degraded_counts)

        # Scatter plot for all proteins (lightseagreen)
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

        add_stats_to_plot(log_avg_production_rate, log_avg_loss_rate, fig)

        # Layout settings
        fig.update_layout(
            title=f"{plot_title}<br>Sim ID:{sim_id} ({len(cell_paths)} cells analyzed), n={len(monomerIds)} proteins",
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
                    average_total_monomer_counts[complexation_monomer_indices],
                    average_free_monomer_counts[complexation_monomer_indices],
                    average_complexed_monomer_counts[complexation_monomer_indices],
                    average_eq_complexed_monomer_counts[complexation_monomer_indices],
                    avg_diluted_counts[complexation_monomer_indices],
                    avg_degraded_counts[complexation_monomer_indices]
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
                title=f"{plot_title}<br>Sim ID:{sim_id} ({len(cell_paths)} cells analyzed), n={len(complexation_monomer_indices)} proteins plotted",
                title_font=dict(size=12),
                xaxis_title="Log10 Average Production Rate",
                yaxis_title="Log10 Average Loss Rate",
                width=700, height=700,
                showlegend=True,
            )

            # Save the complexation plot
            complexation_plot_name = plotOutFileName+ "_complexation_monomers_" + sim_id + ".html"
            fig_complexation.write_html(os.path.join(plotOutDir, complexation_plot_name))

            # Create separate figure for monomers in equilibrium complexes
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
                    average_total_monomer_counts[equilibrium_monomer_indices],
                    average_free_monomer_counts[equilibrium_monomer_indices],
                    average_complexed_monomer_counts[equilibrium_monomer_indices],
                    average_eq_complexed_monomer_counts[equilibrium_monomer_indices],
                    avg_diluted_counts[equilibrium_monomer_indices],
                    avg_degraded_counts[equilibrium_monomer_indices]
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
                title=f"{plot_title}<br>Sim ID:{sim_id} ({len(cell_paths)} cells analyzed), n={len(equilibrium_monomer_indices)} proteins plotted",
                title_font=dict(size=12),
                xaxis_title="Log10 Average Production Rate",
                yaxis_title="Log10 Average Loss Rate",
                width=700, height=700
            )

            # Save the equilibrium plot
            equilibrium_plot_name = plotOutFileName + "_" + "equilibrium_monomers_" + sim_id + ".html"
            fig_equilibrium.write_html(os.path.join(plotOutDir, equilibrium_plot_name))

            # Create separate figure for the monomer that makes no complexes
            fig_no_complex = go.Figure()


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
                    average_total_monomer_counts[no_complex_monomer_indices],
                    average_free_monomer_counts[no_complex_monomer_indices],
                    average_complexed_monomer_counts[no_complex_monomer_indices],
                    average_eq_complexed_monomer_counts[no_complex_monomer_indices],
                    avg_diluted_counts[no_complex_monomer_indices],
                    avg_degraded_counts[no_complex_monomer_indices]
                )

            fig_no_complex.add_trace(go.Scatter(
                    x=log_avg_production_rate[no_complex_monomer_indices],
                    y=log_avg_loss_rate[no_complex_monomer_indices],
                    mode='markers',
                    hovertext=hover_info_no_complex,
                    marker=dict(size=4, color='salmon', opacity=.6),
                    name=f'monomers that do not form \nany complexes (n={len(no_complex_monomer_indices)})'
                ))

            add_lines(fig_no_complex)

            fig_no_complex.update_layout(
                title=f"{plot_title}<br>Sim ID:{sim_id} ({len(cell_paths)} cells analyzed), n={len(no_complex_monomer_indices)} proteins plotted",
                title_font=dict(size=12),
                xaxis_title="Log10 Average Production Rate",
                yaxis_title="Log10 Average Loss Rate",
                width=700, height=700
            )

            # Save the no complex plot
            no_complex_plot_name = plotOutFileName + "_non_complex_forming_monomers_" + sim_id + ".html"
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


            fig_combined.update_layout(
                title=f"{plot_title}<br>Sim ID:{sim_id} ({len(cell_paths)} cells analyzed), n={len(monomerIds)} proteins plotted",
                xaxis_title="Log10 Average Production Rate",
                yaxis_title="Log10 Average Loss Rate",
                width=1000, height=700,
                margin=dict(l=40, r=200, t=50, b=40),
                legend=dict(
                    x=1.01,
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
            combined_plot_name = plotOutFileName + "_all_monomer_types_" + sim_id + ".html"
            fig_combined.write_html(os.path.join(plotOutDir, combined_plot_name))


if __name__ == '__main__':
    Plot().cli()
