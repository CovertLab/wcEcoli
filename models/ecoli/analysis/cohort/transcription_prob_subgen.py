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
from wholecell.analysis.analysis_tools import (exportFigure, stacked_cell_identification,
	read_bulk_molecule_counts, read_stacked_bulk_molecules, read_stacked_columns)
from wholecell.io.tablereader import TableReader
from wholecell.containers.bulk_objects_container import BulkObjectsContainer


IGNORE_FIRST_N_GENS = 8
SEED_RANGE = np.arange(0, 60)
BATCH_SIZE = 100  # Process cells in batches to avoid memory issues

class Plot(cohortAnalysisPlot.CohortAnalysisPlot):
    def do_plot(self, variantDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
        with open(simDataFile, 'rb') as f:
            sim_data = pickle.load(f)

        # Get mRNA cistron data
        transcription = sim_data.process.transcription
        cistron_ids = transcription.cistron_data["id"]
        is_mRNA = transcription.cistron_data['is_mRNA']
        mRNA_cistron_indexes = np.where(is_mRNA)[0]
        mRNA_cistron_ids = np.array([cistron_ids[x] for x in mRNA_cistron_indexes])

        # Ignore data from predefined number of generations per seed
        if self.ap.n_generation <= IGNORE_FIRST_N_GENS:
            print('Skipping analysis - not enough generations run.')
            return
        cell_paths = self.ap.get_cells(
            generation=np.arange(IGNORE_FIRST_N_GENS, self.ap.n_generation),seed=SEED_RANGE,
            only_successful=True)

        print('Analyzing %d cells...' % len(cell_paths))

        # There are 4346 mRNA ids with counts
        RNA_reader = TableReader(
            os.path.join(cell_paths[0], 'simOut', 'RNACounts'))
        mRNA_cistron_ids_counts_table = RNA_reader.readAttribute('mRNA_cistron_ids')
        RNA_reader.close()

        # Process simulatedSynthProb in batches using online mean calculation
        print('Computing simulated synthesis probabilities...')
        running_sum = None
        total_count = 0
        
        for i in range(0, len(cell_paths), BATCH_SIZE):
            batch_paths = cell_paths[i:i+BATCH_SIZE]
            print(f'  Processing batch {i//BATCH_SIZE + 1}/{(len(cell_paths)-1)//BATCH_SIZE + 1}')
            
            batch_data = read_stacked_columns(
                batch_paths, 'RnaSynthProb', 'actual_rna_synth_prob_per_cistron',
                remove_first=True)
            
            if running_sum is None:
                running_sum = batch_data.sum(axis=0)
            else:
                running_sum += batch_data.sum(axis=0)
            
            total_count += batch_data.shape[0]
            del batch_data  # Free memory
        
        simulatedSynthProb = (running_sum / total_count)[mRNA_cistron_indexes]
        del running_sum

        # Process mRNA existence in batches
        print('Computing mRNA existence probabilities...')
        mRNA_exists_sum = None
        total_cells = 0
        
        for i in range(0, len(cell_paths), BATCH_SIZE):
            batch_paths = cell_paths[i:i+BATCH_SIZE]
            print(f'  Processing batch {i//BATCH_SIZE + 1}/{(len(cell_paths)-1)//BATCH_SIZE + 1}')
            
            mRNA_exists_batch = read_stacked_columns(
                batch_paths, 'RNACounts', 'mRNA_cistron_counts', 
                ignore_exception=True, fun=lambda x: x.sum(axis=0) > 0)
            
            if mRNA_exists_sum is None:
                mRNA_exists_sum = mRNA_exists_batch.sum(axis=0)
            else:
                mRNA_exists_sum += mRNA_exists_batch.sum(axis=0)
            
            total_cells += mRNA_exists_batch.shape[0]
            del mRNA_exists_batch  # Free memory

        # Divide by total number of cells to get probability
        p_mRNA_exists_in_gen = mRNA_exists_sum / total_cells
        del mRNA_exists_sum

        # Write data to table
        print('Writing output...')
        with open(os.path.join(plotOutDir, plotOutFileName + '.tsv'), 'w') as f:
            writer = csv.writer(f, delimiter='\t')
            writer.writerow([
                'mRNA_cistron_id',
                'actual_rna_synth_prob',
                'prob_mRNA_exists_in_gen',
            ])

            for i in range(len(mRNA_cistron_ids_counts_table)):
                writer.writerow([
                    mRNA_cistron_ids_counts_table[i], simulatedSynthProb[i], p_mRNA_exists_in_gen[i]
                ])
        

if __name__ == '__main__':
    Plot().cli()