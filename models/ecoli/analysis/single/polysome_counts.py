"""
Histogram of the average polysome counts in the system.
The polysome count analysis takes number of ribosomes attached to mRNA
at each time point, calculates the sum of polysomes across all time points,
and takes the average of polysomes count by dividing with the number of
time points in the system.

The polysome mw analysis takes mass of each unique mRNA, mass of ribosome 30s
and 50s bulk molecules, protein mass attached to each ribosome, and
the number of ribosomes attached to mRNA at each time point to calculate
the molecular weight of polysome through the relationship
polysome_mw = mRNA_mw + protein_mass + n_ribosome * ribosome_mw.
"""

import os

from matplotlib import pyplot as plt
import numpy as np
import pickle

from wholecell.io.tablereader import TableReader
from wholecell.utils import units
from wholecell.analysis.analysis_tools import exportFigure
from models.ecoli.analysis import singleAnalysisPlot

BIN_WIDTH = 2E-3 # Bin width to the polysome mw histogram in fg

class Plot(singleAnalysisPlot.SingleAnalysisPlot):
    def do_plot(self, simOutDir, plotOutDir, plotOutFileName, simDataFile,
                validationDataFile, metadata):
        with open(simDataFile, 'rb') as f:
            sim_data = pickle.load(f)

        # Listeners used
        ribosome_reader = TableReader(os.path.join(simOutDir, 'RibosomeData'))

        # Load data
        mRNA_TU_index = ribosome_reader.readColumn('mRNA_TU_index')
        n_ribosome_on_mRNA = ribosome_reader.readColumn('n_ribosomes_on_each_mRNA')
        protein_mass_on_mRNA = ribosome_reader.readColumn('protein_mass_on_polysomes')
        n_avogadro = sim_data.constants.n_avogadro

        # Initialize variables
        rna_mw = (
            sim_data.process.transcription.rna_data['mw']/
            n_avogadro).asNumber(units.fg)
        ribosome_30S_mass = sim_data.getter.get_mass(
            sim_data.molecule_ids.s30_full_complex)
        ribosome_50S_mass = sim_data.getter.get_mass(
            sim_data.molecule_ids.s50_full_complex)
        ribosome_mass = (
            (ribosome_30S_mass + ribosome_50S_mass)/
            n_avogadro).asNumber(units.fg)

        mRNA_TU_index_flattened = mRNA_TU_index.ravel()
        n_ribosome_on_mRNA_flattened = n_ribosome_on_mRNA.ravel()
        protein_mass_on_mRNA_flattened = protein_mass_on_mRNA.ravel()

        mRNA_TU_index_flattened_clean = mRNA_TU_index_flattened[
            ~(np.isnan(mRNA_TU_index_flattened))].astype(int)
        n_ribosome_on_each_mRNA_flattened_clean = n_ribosome_on_mRNA_flattened[
            ~(np.isnan(n_ribosome_on_mRNA_flattened))].astype(int)
        highest_ribosome_count = n_ribosome_on_each_mRNA_flattened_clean.max()
        protein_mass_on_mRNA_all_time_step = protein_mass_on_mRNA_flattened[
            ~(np.isnan(protein_mass_on_mRNA_flattened))]


        # Count polysomes at all time steps and take average
        bincount_ribosomes_all_time_step = np.bincount(
            n_ribosome_on_each_mRNA_flattened_clean)
        average_polysome_count = (bincount_ribosomes_all_time_step/
            n_ribosome_on_mRNA.shape[0])

        # Calculate polysomes mw at all times steps and take average
        rna_polysome_mw_all_time_step = rna_mw[mRNA_TU_index_flattened_clean]
        ribosome_polysome_mw_all_time_step = (
                ribosome_mass * n_ribosome_on_each_mRNA_flattened_clean)
        polysome_mw_all_time_step = (
            rna_polysome_mw_all_time_step +
            ribosome_polysome_mw_all_time_step +
            protein_mass_on_mRNA_all_time_step)

        polysome_count, polysome_mw_bin= np.histogram(polysome_mw_all_time_step,
            bins = np.arange(min(polysome_mw_all_time_step),
            max(polysome_mw_all_time_step) + BIN_WIDTH, BIN_WIDTH))
        average_polysome_mw_count = polysome_count / mRNA_TU_index.shape[0]

        # Calculation average polysome mass composition across time steps
        average_polysome_mw_across_time = np.mean(polysome_mw_all_time_step)
        percent_rna_polysome_mw_across_time = (
            np.mean(rna_polysome_mw_all_time_step)/
            average_polysome_mw_across_time * 100)
        percent_ribosome_polysome_mw_across_time = (
            np.mean(ribosome_polysome_mw_all_time_step)/
            average_polysome_mw_across_time * 100)
        percent_protein_polysome_mw_across_time = (
            np.mean(protein_mass_on_mRNA_all_time_step)/
            average_polysome_mw_across_time * 100)

        fig, ax = plt.subplots(3, 1, figsize = (8, 10),
            gridspec_kw={'height_ratios': [10,10,1]})
        polysome_count_ax = ax[0]
        polysome_mw_ax = ax[1]
        polysome_mw_composition = ax[2]

        # Plot polysome count
        # add 1 to include the highest ribosome count
        polysome_count_ax.bar(list(range(highest_ribosome_count + 1)),
            average_polysome_count, align = 'center')
        polysome_count_ax.set_xlabel("Count of Ribosomes Attached to Individual mRNA")
        polysome_count_ax.set_ylabel("Count of mRNA")

        # Plot polysome mw
        polysome_mw_ax.bar(polysome_mw_bin[:-1], average_polysome_mw_count,
            width = np.diff(polysome_mw_bin), align = 'edge')
        polysome_mw_ax.set_xlabel("Polysome Molecular Weight [fg]")
        polysome_mw_ax.set_ylabel("Count of polysome")

        # Plot mw composition of polysome
        polysome_mw_composition.set_frame_on(False)
        polysome_mw_composition.barh('polysome mw \n composition',
            percent_ribosome_polysome_mw_across_time,
            label = "Ribosome")
        polysome_mw_composition.barh('polysome mw \n composition',
            percent_rna_polysome_mw_across_time,
            left = percent_ribosome_polysome_mw_across_time,
            label = "mRNA")
        polysome_mw_composition.barh('polysome mw \n composition',
            percent_protein_polysome_mw_across_time,
            left = 100-percent_protein_polysome_mw_across_time,
            label = "Protein")
        polysome_mw_composition.text(percent_ribosome_polysome_mw_across_time/2,
            0, "{:.0f}%".format(percent_ribosome_polysome_mw_across_time),
            ha = 'center')
        polysome_mw_composition.text(
            percent_ribosome_polysome_mw_across_time,
            0, "{:.0f}%".format(percent_rna_polysome_mw_across_time))
        polysome_mw_composition.text(100,
            0, "{:.0f}%".format(percent_protein_polysome_mw_across_time))
        polysome_mw_composition.legend(ncol = 3, bbox_to_anchor = (0.75,-0.4))

        fig.tight_layout()

        exportFigure(plt,plotOutDir, plotOutFileName, metadata)
        plt.close("all")

if __name__ == '__main__':
    Plot().cli()

