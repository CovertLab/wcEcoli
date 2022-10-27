"""
Histogram of the average polysome counts in the system.
The polysome count analysis script takes number of ribosomes attached to mRNA at each time point,
calculates the sum of polysomes across all time points, and takes the average of polysomes count by dividing with
the number of time points in the system.
"""

import os

from matplotlib import pyplot as plt
import numpy as np

from wholecell.io.tablereader import TableReader
from wholecell.analysis.analysis_tools import exportFigure
from models.ecoli.analysis import singleAnalysisPlot


class Plot(singleAnalysisPlot.SingleAnalysisPlot):
    def do_plot(self, simOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):

        # Listeners used
        ribosome_reader = TableReader(os.path.join(simOutDir, 'RibosomeData'))

        # Load data
        n_ribosome_on_each_mRNA = ribosome_reader.readColumn('n_ribosomes_on_each_mRNA')

        # Initialize variables
        n_ribosome_on_each_mRNA_flattened = n_ribosome_on_each_mRNA.ravel()
        highest_ribosome_count = int(max(n_ribosome_on_each_mRNA_flattened))+1

        # Count polysomes at all time steps and take average
        n_ribosome_on_each_mRNA_flattened_clean = n_ribosome_on_each_mRNA_flattened[
                                                  np.logical_not(np.isnan(n_ribosome_on_each_mRNA_flattened))]
        bincount_ribosomes_all_time_step = np.bincount(n_ribosome_on_each_mRNA_flattened_clean.tolist())
        average_polysome_count = bincount_ribosomes_all_time_step/len(n_ribosome_on_each_mRNA)


        fig = plt.figure(figsize = (8, 10))

        plt.hist(np.arange(highest_ribosome_count),
                 weights = average_polysome_count,
                 bins = highest_ribosome_count-1,
                 align = 'left')
        plt.xlabel("Count of Ribosomes Attached to Individual mRNA")
        plt.ylabel("Count of mRNA")

        fig.tight_layout()

        exportFigure(plt,plotOutDir, plotOutFileName, metadata)
        plt.close("all")

if __name__ == '__main__':
    Plot().cli()

