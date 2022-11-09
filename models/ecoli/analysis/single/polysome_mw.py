"""
Histogram of the polysome molecular weight (mw) in the system.
The polysome mw analysis script takes mass of each unique mRNA, mass of ribosome bulk molecules,
and the number of ribosomes attached to mRNA at each time point to calculate the molecular weight of polysome through
the relationship polysome mw = mRNA mw + n_ribosome * ribosome mw.
"""

import os

from matplotlib import pyplot as plt
import numpy as np
from six.moves import cPickle

from wholecell.io.tablereader import TableReader
from wholecell.utils import units
from wholecell.analysis.analysis_tools import exportFigure
from models.ecoli.analysis import singleAnalysisPlot

BIN_WIDTH = 5 # Bin width to the polysome mw histogram in kg/mmol

class Plot(singleAnalysisPlot.SingleAnalysisPlot):
    def do_plot(self, simOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
        with open(simDataFile, 'rb') as f:
            sim_data = cPickle.load(f)

        # Listeners used
        ribosome_reader = TableReader(os.path.join(simOutDir, 'RibosomeData'))

        # Load data
        mRNA_TU_index = ribosome_reader.readColumn('mRNA_TU_index')
        n_ribosome_on_each_mRNA = ribosome_reader.readColumn('n_ribosomes_on_each_mRNA')

        # Initialize variables for polysome mw
        rna_mw = sim_data.process.transcription.rna_data['mw'].asNumber(units.kg/units.mmol) # convert unit to kg/mmol
        ribosome_30S_mass = sim_data.getter.get_mass(sim_data.molecule_ids.s30_full_complex)
        ribosome_50S_mass = sim_data.getter.get_mass(sim_data.molecule_ids.s50_full_complex)
        ribosome_mass = (ribosome_30S_mass + ribosome_50S_mass).asNumber(units.kg/units.mmol)  # convert unit to kg/mmol

        mRNA_TU_index_flattened = mRNA_TU_index.ravel()
        mRNA_TU_index_flattened_clean = mRNA_TU_index_flattened[
                                        np.logical_not(np.isnan(mRNA_TU_index_flattened))].astype(int)
        n_ribosome_on_each_mRNA_flattened = n_ribosome_on_each_mRNA.ravel()
        n_ribosome_on_each_mRNA_flattened_clean = n_ribosome_on_each_mRNA_flattened[
                                                  np.logical_not(np.isnan(n_ribosome_on_each_mRNA_flattened))].astype(int)

        # Calculate polysomes mw at all times steps and take average
        polysome_mw_all_time_step = rna_mw[mRNA_TU_index_flattened_clean] + ribosome_mass * n_ribosome_on_each_mRNA_flattened_clean
        polysome_mw_hist, bin_edges = np.histogram(polysome_mw_all_time_step,
                                                   bins=np.arange(min(polysome_mw_all_time_step),
                                                                  max(polysome_mw_all_time_step) + BIN_WIDTH,
                                                                  BIN_WIDTH))
        average_polysome_mw = polysome_mw_hist/mRNA_TU_index.shape[0]

        # Plot mRNA mw
        fig = plt.figure(figsize=(8, 10))
        plt.hist(np.linspace(min(polysome_mw_all_time_step), max(polysome_mw_all_time_step), len(average_polysome_mw)),
                 weights = average_polysome_mw,
                 bins = bin_edges)
        plt.xlabel("Polysome Molecular Weight [kg/mmol]")
        plt.ylabel("Count of polysome")
        plt.ticklabel_format(useOffset = False, style = 'plain')

        fig.tight_layout()

        exportFigure(plt,plotOutDir, plotOutFileName, metadata)
        plt.close("all")

if __name__ == '__main__':
    Plot().cli()

