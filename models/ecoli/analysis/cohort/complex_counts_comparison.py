"""
Analysis script to compare counts of protein complexes between different sets
of sims.
"""
import numpy as np
import os
import pickle

from matplotlib import pyplot as plt

from models.ecoli.analysis import cohortAnalysisPlot
from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from wholecell.analysis.analysis_tools import (exportFigure,
	read_stacked_bulk_molecules, read_stacked_columns)

MASTER_VARIANT_DIR = '/home/ggsun/projects/wcEcoli/out/operon_complex_test_master/wildtype_000000/'

class Plot(cohortAnalysisPlot.CohortAnalysisPlot):
	def do_plot(self, variantDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		with open(simDataFile, 'rb') as f:
			sim_data = pickle.load(f)

		ap_master = AnalysisPaths(MASTER_VARIANT_DIR, cohort_plot=True)
		ap_operon = AnalysisPaths(variantDir, cohort_plot=True)

		cell_paths_master = ap_master.get_cells()
		cell_paths_operon = ap_operon.get_cells()

		# Load data
		## Simple stacking functions for data from all cells
		complex_names = sim_data.process.complexation.ids_complexes
		(complex_counts_operon, ) = read_stacked_bulk_molecules(cell_paths_operon, (complex_names, ))
		(complex_counts_master, ) = read_stacked_bulk_molecules(cell_paths_master, (complex_names, ))

		log_diff = np.log10(
			(complex_counts_operon.mean(axis=0) + 1)/(complex_counts_master.mean(axis=0) + 1))
		log_diff_mean = log_diff.mean()
		log_diff_std = log_diff.std()

		# Get mask for complexes where log10 diff is higher than mean + 2*std
		operon_higher_mask = log_diff > log_diff_mean + 2*log_diff_std
		operon_higher_complex_names = [complex_names[i] for i in np.nonzero(operon_higher_mask)[0]]

		monomer_counts_operon = read_stacked_columns(
			cell_paths_operon, 'MonomerCounts', 'monomerCounts')
		monomer_counts_master = read_stacked_columns(
			cell_paths_master, 'MonomerCounts', 'monomerCounts')

		plt.figure(figsize=(8, 8))
		plt.plot([1e-2, 5e4], [1e-2, 5e4], linestyle='--', c='gray')
		plt.scatter(
			complex_counts_master.mean(axis=0) + 1,
			complex_counts_operon.mean(axis=0) + 1,
			color='k', s=5)
		# Highlight counts where operon is +2SD higher
		plt.scatter(
			complex_counts_master.mean(axis=0)[operon_higher_mask] + 1,
			complex_counts_operon.mean(axis=0)[operon_higher_mask] + 1,
			color='r', s=5)
		plt.xscale('log')
		plt.yscale('log')
		plt.xlim([0.5, 5e4])
		plt.ylim([0.5, 5e4])
		plt.xlabel('log10(complex counts + 1) (master)')
		plt.ylabel('log10(complex counts + 1) (operon)')
		plt.tight_layout()
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close('all')

		monomer_id_to_rna_set = {
			monomer['id']: monomer['rna_set']
			for monomer in sim_data.process.translation.monomer_data
			}

		with open(os.path.join(plotOutDir, 'complexes_with_higher_counts_in_operon.tsv'), 'w') as f:
			f.write('id\tmonomer_list\trna_list\tlog10_diff\n')
			for i in np.nonzero(operon_higher_mask)[0]:
				monomers = [x for x in sim_data.process.complexation.get_monomers(complex_names[i])['subunitIds']]
				rna_list = []
				for monomer in monomers:
					rna_list.extend(monomer_id_to_rna_set[monomer])

				f.write('%s\t%s\t%s\t%.2f\n' % (
					complex_names[i],
					monomers,
					list(set(rna_list)),
					log_diff[i])
					)

		plt.figure(figsize=(8, 8))
		plt.plot([1e-2, 1e6], [1e-2, 1e6], linestyle='--', c='gray')
		plt.scatter(monomer_counts_master.mean(axis=0) + 1,
			monomer_counts_operon.mean(axis=0) + 1, color='k', s=5)
		plt.xscale('log')
		plt.yscale('log')
		plt.xlim([0.5, 1e6])
		plt.ylim([0.5, 1e6])
		plt.xlabel('log10(monomer counts + 1) (master)')
		plt.ylabel('log10(monomer counts + 1) (operon)')
		plt.tight_layout()
		exportFigure(plt, plotOutDir, plotOutFileName + '_monomer', metadata)
		plt.close('all')


if __name__ == '__main__':
	Plot().cli()
