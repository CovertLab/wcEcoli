"""
Plot polysome histogram
"""

from __future__ import absolute_import, division, print_function

import os

from matplotlib import pyplot as plt
from matplotlib.ticker import FormatStrFormatter
import numpy as np

from wholecell.io.tablereader import TableReader
from wholecell.analysis.analysis_tools import exportFigure
from models.ecoli.analysis import singleAnalysisPlot
from wholecell.utils import units


class Plot(singleAnalysisPlot.SingleAnalysisPlot):
    def do_plot(self, simOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
        # with open(simDataFile, 'rb') as f:
        #     sim_data = cPickle.load(f)

        # Listeners used
        # main_reader = TableReader(os.path.join(simOutDir,'Main'))
        ribosome_reader = TableReader(os.path.join(simOutDir,'RibosomeData'))

        # Load data
        # initial_time = main_reader.readAttribute('initialTime')
        # time = main_reader.readColumn('time') - initial_time
        num_ribosome_on_mRNA = ribosome_reader.readColumn('num_ribosome_on_mRNA')
        mRNA_unique_index = ribosome_reader.readColumn('mRNA_unique_index')

        highest_ribosome_count = 30
        bincount_ribosome_time_step = np.zeros((len(num_ribosome_on_mRNA),highest_ribosome_count))
        bincount_ribosome_time_step[:] = np.NaN
        average_polysome_count = np.zeros(highest_ribosome_count+1)

        for index in range(len(num_ribosome_on_mRNA)):
            num_ribosome_time_step = num_ribosome_on_mRNA[index,:]
            num_ribosome_time_step_nan_removed = num_ribosome_time_step[np.logical_not(np.isnan(num_ribosome_time_step))]

            bincount_ribosome = np.bincount(num_ribosome_time_step_nan_removed.tolist()) # CHECK THIS LINE ALLOCATE BINCOUNT
            bincount_ribosome_time_step[index, :len(bincount_ribosome)] = bincount_ribosome
            # num_ribosome_time_step = [num_ribosome_on_mRNA[index,:][np.logical_not(np.isnan(num_ribosome_on_mRNA[index,:]))] for index in range(len(num_ribosome_on_mRNA))]

        for col_index in range(highest_ribosome_count):
            polysome_distribution_per_num_ribosome = bincount_ribosome_time_step[:, col_index]

            polysome_distribution_per_num_ribosome_clean = polysome_distribution_per_num_ribosome[np.logical_not(np.isnan(polysome_distribution_per_num_ribosome))]
            if polysome_distribution_per_num_ribosome_clean.size == 0:
                average_polysome_count[col_index] = 0
            else:
                average_polysome_count[col_index] = np.mean(polysome_distribution_per_num_ribosome_clean)

        # import ipdb;
        # ipdb.set_trace()

        # average_ribosome_on_mRNA = []
        # print(np.shape(mRNA_unique_index))
        # # print(mRNA_unique_index[:,:10])
        # print(len(np.unique(mRNA_unique_index)))
        # for unique_index in np.unique(mRNA_unique_index)[:10]:
        #     mRNA_indexes = np.where(mRNA_unique_index == unique_index)
        #     average_ribosome_on_mRNA.append(sum(num_ribosome_on_mRNA[mRNA_indexes])/len(mRNA_unique_index[mRNA_indexes]))
        # print(average_ribosome_on_mRNA)
        # print(np.unique(mRNA_unique_index)[:100])
        # average_ribosome_on_mRNA = [sum(num_ribosome_on_mRNA[np.where(mRNA_unique_index == unique_index)])/sum(mRNA_unique_index[np.where(mRNA_unique_index == unique_index)]) for unique_index in np.unique(mRNA_unique_index)]
        # FRIDAY: NO NEED TO PRETAIN MRNA INFORMATION, JUST TAKE AVERAGE OF RIBOSOME COUNTS AT EVERY TIME STEP

        fig = plt.figure(figsize = (8.5, 15))

        plt.subplot(2, 1, 1) # NOT DONE WITH PLOTTING YET
        plt.hist(np.arange(highest_ribosome_count+1), weights = average_polysome_count, bins = highest_ribosome_count)
        plt.xlabel("Count of Ribosomes Attached to Individual mRNA")
        plt.ylabel("Count of mRNA")

        # plt.subplot(2, 1, 2)
        # plt.hist(np.arange(len(avg_polysome_counts)), weights = avg_polysome_counts, bins = range(len(avg_polysome_counts)))
        # plt.xlabel("Count of Ribosomes Attached to Individual mRNA")
        # plt.ylabel("Count of mRNA")

        fig.tight_layout()

        exportFigure(plt,plotOutDir, plotOutFileName, metadata)
        plt.close("all")


if __name__ == '__main__':
    Plot().cli()

