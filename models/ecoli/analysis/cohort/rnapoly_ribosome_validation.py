"""
Plots the RNA polymerase and ribosome abundances for all seeds and all
generations as a violin plot.

Notes
-----
There are hard-coded validation values:
	- one represents the average cell (defined as 44% along the cell cycle in
	age - Neidhardt et al. Physiology of the Bacterial Cell. 1931). To be
	comparable with the validation data, the molecule abundances are taken from
	44% along the cell cycle.
	- RNA polymerase and ribosome abundances reported in Bremer & Dennis 1996,
	Table 3, for the 40 minute doubling cell.

"""

from __future__ import absolute_import

import os
import cPickle
from multiprocessing import Pool
import numpy as np
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

FIGSIZE = (3.5, 3.5)
X_RNAP = 1
X_RIBO = 4
X_RNAP_VAL = 0
X_RIBO_VAL = 3

RNAP_VALIDATION = 5e3    # Bremer & Dennis 1996, Table 3, 40 min doubling
RIBO_VALIDATION = 26.3e3 # Bremer & Dennis 1996, Table 3, 40 min doubling


def mp_worker(sim_dir):
	sim_out_dir = os.path.join(sim_dir, 'simOut')
	rnap_count_avg_cell = None
	ribosome_count_avg_cell = None

	try:
		(rnap_count, ribosome_30s_count, ribosome_50s_count) = read_bulk_molecule_counts(
			sim_out_dir, (
				[rnap_id],
				[ribosome_30s_id],
				[ribosome_50s_id]))

		unique_molecule_reader = TableReader(os.path.join(sim_out_dir, 'UniqueMoleculeCounts'))
		unique_molecule_ids = unique_molecule_reader.readAttribute('uniqueMoleculeIds')
		unique_molecule_counts = unique_molecule_reader.readColumn('uniqueMoleculeCounts')
		unique_molecule_reader.close()

		index_ribosome = unique_molecule_ids.index('activeRibosome')
		index_rnap = unique_molecule_ids.index('activeRnaPoly')
		ribosome_active_count = unique_molecule_counts[:, index_ribosome]
		rnap_active_count = unique_molecule_counts[:, index_rnap]


		index_average_cell = int(len(rnap_count) * CELL_CYCLE_FRACTION)
		rnap_count_avg_cell = rnap_count[index_average_cell] + rnap_active_count[index_average_cell]
		ribosome_count_avg_cell = ribosome_active_count[index_average_cell] + min(
			ribosome_30s_count[index_average_cell],
			ribosome_50s_count[index_average_cell])

	except Exception as e:
		print('Excluded from analysis due to broken files: {}'.format(sim_out_dir))

	return (rnap_count_avg_cell, ribosome_count_avg_cell)


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

		global rnap_id
		global ribosome_30s_id
		global ribosome_50s_id

		rnap_id = sim_data.moleculeIds.rnapFull
		ribosome_30s_id = sim_data.moleculeIds.s30_fullComplex
		ribosome_50s_id = sim_data.moleculeIds.s50_fullComplex

		p = Pool(parallelization.cpus())
		output = p.map(mp_worker, sim_dirs)
		p.close()
		p.join()

		# Filter output from broken files
		rnap_counts = [x[0] for x in output if x[0]]
		ribosome_counts = [x[1] for x in output if x[1]]

		if not len(rnap_counts) or not len(ribosome_counts):
			print('Skipping plot due to no viable sims.')
			return

		# Plot
		fig, ax = plt.subplots(1, 1, figsize=FIGSIZE)

		ax.violinplot(rnap_counts, [X_RNAP])
		ax.scatter(X_RNAP_VAL, RNAP_VALIDATION)

		ax.violinplot(ribosome_counts, [X_RIBO])
		ax.scatter(X_RIBO_VAL, RIBO_VALIDATION)

		ax.set_xticks([X_RNAP_VAL, X_RNAP, X_RIBO_VAL, X_RIBO])
		ax.set_xticklabels([])
		ax.set_yticklabels([])
		exportFigure(plt, plotOutDir, '{}__clean'.format(plotOutFileName), metadata)
		
		ax.set_xticklabels([
			'Validation\n(RNAP)', 'Simulated\n(RNAP)',
			'Validation\n(RIBO)', 'Simulated\n(RIBO)'], fontsize=6)
		ax.set_title('n = {}'.format(len(rnap_counts)))
		ax.set_ylabel('Molecule abundance (counts)')
		ax.set_yticklabels(ax.get_yticks(), fontsize=6)
		ax.tick_params(which="both", direction="out", right=False, top=False)
		plt.subplots_adjust(left=0.2, bottom=0.2, right=0.8, top=0.8)
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)

		plt.close("all")


if __name__ == "__main__":
	Plot().cli()
