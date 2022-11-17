"""
Histogram of the average polysome counts in the system.
The polysome count analysis script takes number of ribosomes attached to mRNA
at each time point, calculates the sum of polysomes across all time points,
and takes the average of polysomes count by dividing with the number of
time points in the system.
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

        mRNA_TU_index_flattened_clean = mRNA_TU_index_flattened[
            ~(np.isnan(mRNA_TU_index_flattened))].astype(int)
        n_ribosome_on_each_mRNA_flattened_clean = n_ribosome_on_mRNA_flattened[
            ~(np.isnan(n_ribosome_on_mRNA_flattened))].astype(int)
        highest_ribosome_count = n_ribosome_on_each_mRNA_flattened_clean.max()


        # Count polysomes at all time steps and take average
        bincount_ribosomes_all_time_step = np.bincount(
            n_ribosome_on_each_mRNA_flattened_clean)
        average_polysome_count = (bincount_ribosomes_all_time_step/
            n_ribosome_on_mRNA.shape[0])

        # Calculate polysomes mw at all times steps and take average
        polysome_mw_all_time_step = (
            rna_mw[mRNA_TU_index_flattened_clean] +
            ribosome_mass * n_ribosome_on_each_mRNA_flattened_clean)
        polysome_mw_hist, bin_edges = np.histogram(polysome_mw_all_time_step,
            bins = np.arange(min(polysome_mw_all_time_step),
            max(polysome_mw_all_time_step) + BIN_WIDTH, BIN_WIDTH))
        average_polysome_mw = polysome_mw_hist / mRNA_TU_index.shape[0]

        fig = plt.figure(figsize = (8, 10))
        # Plot polysome count
        polysome_count_ax = plt.subplot(2,1,1)
        # add 1 to include the highest ribosome count
        polysome_count_ax.bar(list(range(highest_ribosome_count + 1)),
            average_polysome_count, align = 'center')
        polysome_count_ax.set_xlabel("Count of Ribosomes Attached to Individual mRNA")
        polysome_count_ax.set_ylabel("Count of mRNA")

        # Plot polysome mw
        polysome_mw_ax = plt.subplot(2, 1, 2)
        polysome_mw_ax.bar(bin_edges[:-1], average_polysome_mw,
            width = np.diff(bin_edges), align = 'edge')
        polysome_mw_ax.set_xlabel("Polysome Molecular Weight [fg]")
        polysome_mw_ax.set_ylabel("Count of polysome")

        fig.tight_layout()

        exportFigure(plt,plotOutDir, plotOutFileName, metadata)
        plt.close("all")

if __name__ == '__main__':
    Plot().cli()

