"""
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 9/13/2019
"""

from __future__ import absolute_import

import os
import matplotlib.pyplot as plt

from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from wholecell.io.tablereader import TableReader
from wholecell.analysis.analysis_tools import exportFigure
from wholecell.analysis.analysis_tools import read_bulk_molecule_counts
from models.ecoli.analysis import multigenAnalysisPlot

TITLE = "tolC"

class Plot(multigenAnalysisPlot.MultigenAnalysisPlot):
	def do_plot(self, seedOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		if not os.path.isdir(seedOutDir):
			raise Exception, "seedOutDir does not currently exist as a directory"

		if not os.path.exists(plotOutDir):
			os.mkdir(plotOutDir)

		mrna_id = "EG11009_RNA[c]"
		monomer_id = "EG11009-MONOMER[e]"
		complex_id = "CPLX0-7964[e]"

		# Get all cells
		ap = AnalysisPaths(seedOutDir, multi_gen_plot=True)
		all_dir = ap.get_cells()

		generation_ticks = [0.]

		fig, axes_list = plt.subplots(3, 1, figsize=(11, 8.5))
		ax0, ax1, ax2 = axes_list

		for gen, simDir in enumerate(all_dir):
			simOutDir = os.path.join(simDir, "simOut")

			time = TableReader(os.path.join(simOutDir, "Main")).readColumn("time")
			generation_ticks.append(time[-1])

			(counts,) = read_bulk_molecule_counts(simOutDir, ([mrna_id, monomer_id, complex_id],))

			ax0.plot(time, counts[:, 0])
			ax1.plot(time, counts[:, 1])
			ax2.plot(time, counts[:, 2])

		for ax, ax_title in zip(
				axes_list,
				["mRNA: {}".format(mrna_id),
				 "Monomer: {}".format(monomer_id),
				 "Complex: {}".format(complex_id)]):
			ax.set_xticks(generation_ticks)
			ax.set_xticklabels(range(len(generation_ticks)))
			ax.set_title(ax_title)

		ax2.set_xlabel("Generation")
		plt.suptitle(TITLE)
		plt.subplots_adjust(hspace=0.4)
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close("all")

if __name__ == "__main__":
	Plot().cli()
