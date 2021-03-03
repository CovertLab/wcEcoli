"""
Cohort analysis identifying seeds and generations in which
stalled transcript elongation occurred
"""

import pickle
import os

from matplotlib import pyplot as plt
import numpy as np

from models.ecoli.analysis import cohortAnalysisPlot
from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from wholecell.analysis.analysis_tools import (exportFigure,
	read_bulk_molecule_counts, read_stacked_bulk_molecules, read_stacked_columns)
from wholecell.io.tablereader import TableReader


class Plot(cohortAnalysisPlot.CohortAnalysisPlot):
	def do_plot(self, variantDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		ap = AnalysisPaths(variantDir, cohort_plot=True)
		cell_paths = ap.get_cells()

		seed_gen_stalled = []
		for sim_dir in cell_paths:
			simOutDir = os.path.join(sim_dir, 'simOut')

			# Listeners used
			main_reader = TableReader(os.path.join(simOutDir, 'Main'))
			stallRNAP_counts_reader = TableReader(os.path.join(simOutDir, 'RnapData'))

			# Load data
			time = main_reader.readColumn('time')
			stallRNAP_counts = stallRNAP_counts_reader.readColumn('didStall')

			# add to result if stalled elongation occurred
			if np.sum(stallRNAP_counts) > 0:
				paths = os.path.normpath(simOutDir).split(os.sep)
				seed = paths[-4]
				gen = paths[-3].split('_')[1]
				seed_gen_stalled.append([seed, gen, simOutDir])

		seed_gen_stalled.sort()
		np.savetxt(fname=f"{plotOutDir}/rnapStalledSims.tsv",
				   X=seed_gen_stalled, delimiter='\t', fmt='%s',
				   comments='', header='seed\tgeneration\tsimOutDir')


if __name__ == '__main__':
	Plot().cli()
