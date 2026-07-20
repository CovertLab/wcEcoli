"""
Template for cohort analysis plots
"""

import pickle
import os

from matplotlib import pyplot as plt
# noinspection PyUnresolvedReferences
import numpy as np
from matplotlib import cm

import csv

from wholecell.utils import units
from models.ecoli.analysis import cohortAnalysisPlot
from models.ecoli.analysis.cohort import subgen_common as sc
from wholecell.analysis.analysis_tools import (exportFigure, stacked_cell_identification,
	read_bulk_molecule_counts, read_stacked_bulk_molecules, read_stacked_columns)
from wholecell.io.tablereader import TableReader
from wholecell.containers.bulk_objects_container import BulkObjectsContainer

IGNORE_FIRST_N_GENS = 8
SEED_RANGE = np.arange(0, 60)
BATCH_SIZE = 50


class Plot(cohortAnalysisPlot.CohortAnalysisPlot):
	def do_plot(self, variantDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		with open(simDataFile, 'rb') as f:
			sim_data = pickle.load(f)
			# Ignore data from predefined number of generations per seed
		if self.ap.n_generation <= IGNORE_FIRST_N_GENS:
			print('Skipping analysis - not enough generations run.')
			return
		cell_paths = self.ap.get_cells(
			generation=np.arange(IGNORE_FIRST_N_GENS, self.ap.n_generation), seed = SEED_RANGE,
			only_successful=True)

		# Restrict to strict-successful lineages (completed every generation and
		# no cell at the 180-min doubling cap).
		success = sc.compute_lineage_success(self.ap, self.ap.n_generation)
		cell_paths = sc.filter_cells_to_successful(
			cell_paths, success['successful_seeds'])

		if len(cell_paths) == 0:
			print('No valid cell paths found for this variant. Skipping analysis.')
			return
		
		print('Analyzing %d cells...' % len(cell_paths))

		# There are 4346 mRNA ids with counts
		RNA_reader = TableReader(
				os.path.join(cell_paths[0], 'simOut', 'RNACounts'))
		mRNA_ids = RNA_reader.readAttribute('mRNA_cistron_ids')
		RNA_reader.close()

		mRNA_id_to_index = {
			cistron_id: i for (i, cistron_id)
			in enumerate(mRNA_ids)
		}

		# There are 4539 mRNA ids total w/ gene namesclau
		cistron_id_to_gene_id = {
			cistron['id']: cistron['gene_id']
			for cistron in sim_data.process.transcription.cistron_data
		}

		# There are 4310 mRNA ids with associated protein/monomer ids
		protein_id_to_cistron_id = {
			protein['id']:protein['cistron_id']
			for protein in sim_data.process.translation.monomer_data
		}

		monomer_reader = TableReader(
			os.path.join(cell_paths[0], 'simOut', 'MonomerCounts'))
		monomer_ids = monomer_reader.readAttribute('monomerIds')
		monomer_reader.close()


		# Get indexes of monomers in this subcolumn
		monomer_id_to_index = {
			monomer_id: i for (i, monomer_id)
			in enumerate(monomer_ids)
		}

		monomer_indices = np.array([
			monomer_id_to_index[monomer_id] for monomer_id in monomer_ids
		])

		# order cistrons in the order of monomers ids
		cistron_ids_in_order = np.array([
			protein_id_to_cistron_id[monomer_id] for monomer_id in monomer_ids
		])

		# order gene names in the order of monomers ids
		gene_ids_in_order = np.array([
			cistron_id_to_gene_id[cistron_id] for cistron_id in cistron_ids_in_order
		])

		# Get indices of cistron_ids_in_order
		mRNA_ids_indices = np.array([
			mRNA_id_to_index[cistron_id] for cistron_id
			in cistron_ids_in_order
		])

		mean_monomer_counts = np.zeros(len(monomer_indices))
		mean_mRNA_counts = np.zeros(len(mRNA_ids_indices))
		n_cells = 0

		for i in range(0, len(cell_paths), BATCH_SIZE):
			batch_paths = cell_paths[i:i+BATCH_SIZE]
			print(f'Processing batch {i//BATCH_SIZE + 1}/{(len(cell_paths)-1)//BATCH_SIZE + 1}...')
			
			# Read batch data
			batch_monomer = read_stacked_columns(
				batch_paths, 'MonomerCounts', 'monomerCounts',
				ignore_exception=True)[:, monomer_indices]
			batch_mRNA = read_stacked_columns(
				batch_paths, 'RNACounts', 'mRNA_cistron_counts',
				ignore_exception=True)[:, mRNA_ids_indices]
			
			# Accumulate sums
			mean_monomer_counts += batch_monomer.sum(axis=0)
			mean_mRNA_counts += batch_mRNA.sum(axis=0)
			n_cells += batch_monomer.shape[0]
			
			# Free memory
			del batch_monomer, batch_mRNA

		# Calculate means once, after all batches have been accumulated. (This
		# was previously inside the loop, which divided the running sums by a
		# growing n_cells every batch and corrupted the means.)
		mean_monomer_counts /= n_cells
		mean_mRNA_counts /= n_cells

		# Write data to table
		with open(os.path.join(plotOutDir, plotOutFileName + '_40_seeds_last_11_gens.tsv'), 'w') as f:
			writer = csv.writer(f, delimiter='\t')
			writer.writerow([
				'gene_name', 'cistron_name', 'monomer_name',
				'mean_mRNA_count', 
				'mean_protein_count', 
			])

			for i in monomer_indices:
				writer.writerow([
					gene_ids_in_order[i], cistron_ids_in_order[i], monomer_ids[i],
					mean_mRNA_counts[i], 
					mean_monomer_counts[i],					
				])


if __name__ == '__main__':
	Plot().cli()
