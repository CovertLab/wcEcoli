"""
Plots the ribosome abundances for all seeds and all generations as a violin plot.

Notes
-----
There are hard-coded validation values:
	- one represents the average cell (defined as 44% along the cell cycle in
	  age - Neidhardt et al. Physiology of the Bacterial Cell. 1931). To be
	  comparable with the validation data, the molecule abundances are taken
	  from 44% along the cell cycle.
	- ribosome abundance reported in Bremer & Dennis 2008, Table 3, for the
	  40 minute and 60 minute doubling cells.
"""

from __future__ import absolute_import

import os
import cPickle
from multiprocessing import Pool
from matplotlib import pyplot as plt

from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from wholecell.io.tablereader import TableReader
from wholecell.analysis.analysis_tools import exportFigure
from wholecell.analysis.analysis_tools import read_bulk_molecule_counts
from wholecell.utils import units, parallelization

from models.ecoli.analysis import cohortAnalysisPlot

# First generation (counting from zero) from which to gather doubling time
# values.  If fewer generations were run, this script quits early without
# plotting anything.
FIRST_GENERATION = 2

CELL_CYCLE_FRACTION = 0.44 # Average cell is 44% along its cell cycle length
RIBO_VALIDATION = [26e3, 15e3] # Bremer & Dennis 2008, Table 3
FIGSIZE = (5, 5)
COUNTS_BOUNDS = [0, 3e4]


def mp_worker(sim_dir):
	sim_out_dir = os.path.join(sim_dir, 'simOut')
	ribosome_count_avg_cell = None

	try:
		(ribosome_30s_count, ribosome_50s_count) = read_bulk_molecule_counts(
			sim_out_dir, (
				[ribosome_30s_id],
				[ribosome_50s_id]))

		unique_molecule_reader = TableReader(os.path.join(sim_out_dir, 'UniqueMoleculeCounts'))
		unique_molecule_ids = unique_molecule_reader.readAttribute('uniqueMoleculeIds')
		unique_molecule_counts = unique_molecule_reader.readColumn('uniqueMoleculeCounts')
		unique_molecule_reader.close()

		index_ribosome = unique_molecule_ids.index('activeRibosome')
		ribosome_active_count = unique_molecule_counts[:, index_ribosome]

		index_average_cell = int(len(ribosome_active_count) * CELL_CYCLE_FRACTION)
		ribosome_count_avg_cell = ribosome_active_count[index_average_cell] + min(
			ribosome_30s_count[index_average_cell],
			ribosome_50s_count[index_average_cell])

	except Exception as e:
		print('Excluded from analysis due to broken files: {}'.format(sim_out_dir))

	return ribosome_count_avg_cell


class Plot(cohortAnalysisPlot.CohortAnalysisPlot):
	def do_plot(self, variantDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		if not os.path.isdir(variantDir):
			raise Exception, 'variantDir does not currently exist as a directory'

		if not os.path.exists(plotOutDir):
			os.mkdir(plotOutDir)

		analysis_paths = AnalysisPaths(variantDir, cohort_plot = True)
		n_gens = analysis_paths.n_generation

		# Check for sufficient generations
		if n_gens - 1 < FIRST_GENERATION:
			print 'Not enough generations to plot.'
			return

		sim_dirs = analysis_paths.get_cells(
			generation=range(FIRST_GENERATION, n_gens))

		sim_data = cPickle.load(open(simDataFile, 'rb'))

		global ribosome_30s_id
		global ribosome_50s_id

		ribosome_30s_id = sim_data.moleculeIds.s30_fullComplex
		ribosome_50s_id = sim_data.moleculeIds.s50_fullComplex

		p = Pool(parallelization.cpus())
		output = p.map(mp_worker, sim_dirs)
		p.close()
		p.join()

		# Filter output from broken files
		ribosome_counts = [x for x in output if x]

		if not len(ribosome_counts):
			print('Skipping plot due to no viable sims.')
			return

		# Plot
		fig, ax = plt.subplots(1, 1, figsize=FIGSIZE)
		ax.violinplot(ribosome_counts)
		for val in RIBO_VALIDATION:
			ax.axhline(val, color='r', lw=2)
		ax.set_ylim(*COUNTS_BOUNDS)
		ax.set_xticks([])
		ax.set_yticklabels([])
		exportFigure(plt, plotOutDir, '{}__clean'.format(plotOutFileName), metadata)

		ax.set_title('n = {}'.format(len(ribosome_counts)))
		ax.set_ylabel('Molecule abundance (counts)')
		ax.set_yticklabels(ax.get_yticks())
		plt.subplots_adjust(left=0.2, bottom=0.2, right=0.8, top=0.8)
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close("all")


if __name__ == "__main__":
	Plot().cli()