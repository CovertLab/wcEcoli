from __future__ import absolute_import, division, print_function

import os

import numpy as np
from matplotlib import pyplot as plt

from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from wholecell.io.tablereader import TableReader, TableReaderError

from wholecell.utils.sparkline import whitePadSparklineAxis
from wholecell.analysis.analysis_tools import exportFigure
from models.ecoli.analysis import variantAnalysisPlot
from wholecell.utils import filepath as fp

FONT_SIZE=9


class Plot(variantAnalysisPlot.VariantAnalysisPlot):
	def do_plot(self, inputDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		ap = AnalysisPaths(inputDir, all_variant_plot = True)

		if ap.n_generation <= 1:
			print(f"Plot {__name__} needs more generations of data to plot.")
			return

		fig = plt.figure()
		fig.set_figwidth(15)
		fig.set_figheight(5)

		doublingTimeVariants = [44, 100, 22]

		for variant_composite_index in ap.variants:
			op, varIdx = fp.split_variant_index(variant_composite_index)

			if varIdx == 0:
				plotIdx = 1
			elif varIdx == 1:
				plotIdx = 0
			elif varIdx == 2:
				plotIdx = 2
			else:
				continue

			all_cells = ap.get_cells(
				generation=[2, 3],
				variant=[op * fp.OPERON_PART + varIdx])
			if len(all_cells) == 0:
				continue

			doublingTimes = np.zeros(len(all_cells))
			for idx, simDir in enumerate(all_cells):
				try:
					simOutDir = os.path.join(simDir, "simOut")
					time = TableReader(os.path.join(simOutDir, "Main")).readColumn("time")
				except (TableReaderError, EnvironmentError) as e:
					print('Error reading data for %s; %s' % (simDir, e))
					raise

				doublingTimes[idx] = (time[-1] - time[0]) / 60.

			bins = 16
			ax = plt.subplot2grid((1, 3), (0, plotIdx))
			ax.hist(doublingTimes, bins)
			ax.axvline(x = doublingTimeVariants[varIdx], color = "r", linestyle = "--")

			ax.set_title("%i min" % (doublingTimeVariants[varIdx]), fontsize = FONT_SIZE)

			ax.set_xlabel("Doubling Time (min)", fontsize = FONT_SIZE)

			plt.subplots_adjust(bottom = 0.2)

			whitePadSparklineAxis(ax)

			for tick in ax.yaxis.get_major_ticks():
				tick.label.set_fontsize(FONT_SIZE)
			for tick in ax.xaxis.get_major_ticks():
				tick.label.set_fontsize(FONT_SIZE)

		exportFigure(plt, plotOutDir, plotOutFileName, metadata)


if __name__ == "__main__":
	Plot().cli()
