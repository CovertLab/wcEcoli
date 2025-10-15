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


PLOT_PROTEINS = ['G6890-MONOMER[c]',
                       'PD03938[c]',
                        'G6737-MONOMER[c]',
                        'RPOD-MONOMER[c]',
                        'PD02936[c]',
                        'RED-THIOREDOXIN2-MONOMER[c]',
                 'EG10542-MONOMER[c]'] # ["PD00196", "EG11111-MONOMER", "EG11545-MONOMER", "EG10320-MONOMER", "ADHP-MONOMER", "EG10580-MONOMER", "YJCQ-MONOMER", "G6988-MONOMER", "MOTB-FLAGELLAR-MOTOR-STATOR-PROTEIN"]


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


        # find how many proteins were removed via dilution for each doubling time:
        diluted_counts = np.zeros(((len(doubling_times)-1), len(monomerIDs)))
        diluted_counts_over_time = np.zeros(((len(time)), len(monomerIDs)))
        for i in range(len(doubling_times) -1): # -1 is to account for not doing the last generation
            end_gen = end_generation_indices[i]
            start_gen = start_generation_indices[i+1] # the first start is zero, so skip that
            print(end_gen, start_gen)

            # find the protein counts at the end of the generation:
            monomer_counts_at_gen_end = free_monomer_counts[end_gen,:] # get this for each protein

            # find the protein counts at the start of the next generation:
            monomer_counts_at_gen_start = free_monomer_counts[start_gen,:]

            # find the difference between the two:
            protein_counts_removed = monomer_counts_at_gen_end - monomer_counts_at_gen_start
            diluted_counts[i,:] = protein_counts_removed
            diluted_counts_over_time[start_gen,:] = protein_counts_removed # put it at the start of the next gen in terms of time


        # Compute how many proteins were removed via degradation over the entire sim length:
        degraded_counts = read_stacked_columns(cell_paths, 'MonomerCounts', "monomersDegraded")

        # Compute how many counts were added via elongation over the entire sim length:
        elongated_counts = read_stacked_columns(cell_paths, 'MonomerCounts', "monomersElongated") # can I just use the monomer counts reader like above?

        # Reorganize the data:
        degraded_counts = degraded_counts * -1 # for graphic purposes, make degradation negative
        diluted_counts_over_time = diluted_counts_over_time * -1

        # calculate the net change per time:
        net_rate = degraded_counts + diluted_counts_over_time + elongated_counts

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

        # plot the loss rate and the production rate:
        for protein in PLOT_PROTEINS_revised:
            # Gather the relevant data:
            protein_idx = monomer_idx_dict[protein]
            protein_FMC = free_monomer_counts[:, protein_idx]
            proteins_degraded = (degraded_counts[:, protein_idx])*-1
            monomer_data = sim_data.process.translation.monomer_data[protein_idx]
            deg_rate = monomer_data["deg_rate"]
            measured_HL = (np.log(2) / deg_rate) / 60

            # Generate the plots:
            fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, figsize=(10, 6))

            # Monomer Counts plot
            ax1.plot(time, protein_FMC, color='lightseagreen', label='Free Monomer Counts')
            for i in range(len(end_generation_times)):
                dt = end_generation_times[i]
                ax1.axvline(x=dt, linestyle='--', color="yellowgreen")

                # Obtain the effective HL (total # of proteins degraded / total # of proteins present in a given duration, excluding zeros):
                start_time= start_generation_indices[i]
                end_time = end_generation_indices[i]
                gen_time = time[start_time:end_time]
                gen_protein_counts = protein_FMC[start_time:end_time]
                gen_degraded_counts = proteins_degraded[start_time:end_time] # Question: double check that the counts are being sliced correctly

                # Calculate the effective degradation rate from the effective half life:
                non_zero_counts = np.nonzero(gen_protein_counts)[0]
                k_eff = gen_degraded_counts[non_zero_counts]/gen_protein_counts[non_zero_counts]
                avg_k_eff = np.mean(k_eff) # Question: is it fine to have the mean happen here? otherwise I think python over rounds it to zero
                avg_half_life = (np.log(2) / avg_k_eff)/60
                print(avg_half_life)

                # Generate a plot of the effective degradation rate over time:
                def effective_HL(t, C0, k):
                    return C0 * k * t

                # Plot the effective half life fit:
                k_avg = np.log(2) / (avg_half_life*60)
                C0_fit = gen_protein_counts
                # Question: should I fix this? I think gen_time should not be in here? maybe it should be 1 to end_gen_time?
                time_for_graph = gen_time - np.ones(len(gen_time))*gen_time[0]
                y_data = effective_HL(time_for_graph, C0_fit, k_avg)
                # prevent extra legend entries:
                if i==0:
                    ax1.plot(gen_time, y_data,
                             color='gray', linestyle='--', label="Effective HL fit")
                else:
                    ax1.plot(gen_time, y_data,
                         color='gray', linestyle='--')
                # add text to indicate the average effective half life for the generation:
                y_pos = y_data[0] * .8
                ax1.text(dt-time_for_graph[-1]/2, y_pos, f"effective HL ≈ {avg_half_life:.1f} min", color="black",
                         ha="center", va="bottom", fontsize=5, rotation=0)

                # Plot the measured half life degradation fit as a comparison:
                k_measured = np.log(2)/(measured_HL*60)
                y_data_measured = effective_HL(time_for_graph, C0_fit, k_measured)
                if i==0:
                    ax1.plot(gen_time, y_data_measured,
                         color='orange', linestyle=':', alpha=0.5, label=f"Measured HL ({measured_HL:.1f} min) fit")
                else:
                    ax1.plot(gen_time, y_data_measured,
                         color='orange', linestyle=':', alpha=0.5)

            # Plot specifics:
            sim_name = metadata["description"]
            degradation_rate_combo = sim_data.protein_degradation_combo_selection
            ax1.legend(fontsize=5)
            ax1.set_ylabel("Free Monomer Counts")
            ax1.set_title(f"Monomer counts and influx/efflux behavior over time for {protein}\n Sim ID:{sim_name}; Degradation rate combo: {degradation_rate_combo}")

            # rates
            ax2.plot(time, degraded_counts[:, protein_idx], color='red', alpha=0.5,
                     label='Degradation')
            ax2.plot(time, diluted_counts_over_time[:, protein_idx], color='blue', alpha=0.5,
                     label='Dilution')
            ax2.plot(time, elongated_counts[:, protein_idx], color='green', alpha=0.5,
                     label='Elongation')
            ax2.plot(time, net_rate[:, protein_idx], color='black', linestyle='--', alpha=1,
                     linewidth=0.2, label='Net change')

            ax2.set_xlabel("Time (s)")
            ax2.set_ylabel("Rate (counts/s)") # Question: is it actually counts/s?
            ax2.set_ylim([-10, 10])
            ax2.legend(fontsize=5)
            plt.tight_layout()



            #save the plot:
            file_name = plotOutFileName + "_" +protein
            exportFigure(plt, plotOutDir, file_name, metadata)


if __name__ == '__main__':
    Plot().cli()
