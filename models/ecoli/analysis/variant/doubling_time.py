"""
Compare cell cycle times across variants.

@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 4/26/20
"""

from __future__ import absolute_import, division, print_function

import cPickle
import os

from matplotlib import pyplot as plt
import numpy as np

from models.ecoli.analysis import variantAnalysisPlot
from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from wholecell.analysis.analysis_tools import exportFigure
from wholecell.io.tablereader import TableReader
from wholecell.utils import filepath


def remove_border(ax, bottom=False):
	ax.spines['top'].set_visible(False)
	ax.spines['right'].set_visible(False)
	if bottom:
		ax.spines['bottom'].set_visible(False)
		ax.set_xticks([])
	ax.tick_params(axis='y', labelsize=6)

class Plot(variantAnalysisPlot.VariantAnalysisPlot):
	def do_plot(self, inputDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		if not os.path.isdir(inputDir):
			raise Exception('inputDir does not currently exist as a directory')

		filepath.makedirs(plotOutDir)

		ap = AnalysisPaths(inputDir, variant_plot=True)
		variants = ap.get_variants()

		with open(simDataFile) as f:
			sim_data = cPickle.load(f)

		variant_lengths = []
		variant_counts = []
		labels = []
		for variant in variants:
			lengths = []
			count = 0
			for sim_dir in ap.get_cells(variant=[variant]):
				try:
					sim_out_dir = os.path.join(sim_dir, "simOut")

					# Listeners used
					main_reader = TableReader(os.path.join(sim_out_dir, 'Main'))

					# Load data
					time = main_reader.readColumn('time')
					cycle_length = time[-1] - time[0]

				except:
					continue

				lengths.append(cycle_length / 60)
				count += 1

			variant_lengths.append(np.mean(lengths))
			variant_counts.append(count)
			labels.append(sim_data.moleculeGroups.aaIDs[variant])

		plt.figure()

		plt.subplot(2, 1, 1)
		plt.bar(variants, variant_lengths)
		plt.ylabel('Average cell cycle length (min)', fontsize=8)
		remove_border(plt.gca(), bottom=True)

		plt.subplot(2, 1, 2)
		plt.bar(variants, variant_counts)
		plt.ylabel('Number of variants', fontsize=8)
		plt.xticks(variants, labels, rotation=45, fontsize=6, ha='right')
		remove_border(plt.gca())

		plt.tight_layout()
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close('all')


if __name__ == "__main__":
	Plot().cli()
