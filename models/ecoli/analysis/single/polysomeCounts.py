"""
Plot polysome histogram
"""

from __future__ import absolute_import, division, print_function

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
        num_ribosome_on_mRNA = ribosome_reader.readColumn('num_ribosome_on_mRNA')

        # Initialize variables
        highest_ribosome_count = 50
        bincount_ribosome_time_step = np.zeros((len(num_ribosome_on_mRNA), highest_ribosome_count))
        bincount_ribosome_time_step[:] = np.NaN
        average_polysome_count = np.zeros(highest_ribosome_count)

        # Count polysomes at each timestep
        for time_index in range(len(num_ribosome_on_mRNA)):
            num_ribosome_time_step = num_ribosome_on_mRNA[time_index, :]
            num_ribosome_time_step_clean = num_ribosome_time_step[
                                           np.logical_not(np.isnan(num_ribosome_time_step))]

            bincount_ribosome = np.bincount(num_ribosome_time_step_clean.tolist())
            bincount_ribosome_time_step[time_index, :len(bincount_ribosome)] = bincount_ribosome

        # Get average of polysome count across timesteps
        for col_index in range(highest_ribosome_count):
            polysome_per_num_ribosome = bincount_ribosome_time_step[:, col_index]
            polysome_per_num_ribosome_clean = polysome_per_num_ribosome[
                                              np.logical_not(np.isnan(polysome_per_num_ribosome))]

            if polysome_per_num_ribosome_clean.size == 0:
                average_polysome_count[col_index] = 0
            else:
                average_polysome_count[col_index] = np.mean(polysome_per_num_ribosome_clean)


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

