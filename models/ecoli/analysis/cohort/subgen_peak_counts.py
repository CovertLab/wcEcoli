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
monomers_of_interest = ['GLYCDEH-MONOMER[c]',  # gldA
                        'BETAGALACTOSID-MONOMER[c]',  # lacZ
                        'RIBULOKIN-MONOMER[c]',  # araB
                        'BAES-MONOMER[i]',  # baeS
                        'G6504-MONOMER[o]',  # gfcE
                        'EG11250-MONOMER[c]',  # chpS
                        'EG11222-MONOMER[c]',  # alkA
                        'G7263-MONOMER[c]',  # murQ
                        'EG11249-MONOMER[c]',  # mazF const
                        'EG10466-MONOMER[c]'  # hupA const
                        ]

monomers_of_interest_name_dict = {'GLYCDEH-MONOMER[c]': 'gldA',
                                  'BETAGALACTOSID-MONOMER[c]': 'lacZ',
                                  'RIBULOKIN-MONOMER[c]': 'araB',
                                  'BAES-MONOMER[i]': 'baeS',
                                  'G6504-MONOMER[o]': 'gfcE',
                                  'EG11250-MONOMER[c]': 'chpS',
                                  'EG11222-MONOMER[c]': 'alkA',
                                  'G7263-MONOMER[c]': 'murQ',
                                  'EG11249-MONOMER[c]': 'mazF',
                                  'EG10466-MONOMER[c]': 'hupA'
                                  }

class Plot(cohortAnalysisPlot.CohortAnalysisPlot):
    def do_plot(self, variantDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
        with open(simDataFile, 'rb') as f:
            sim_data = pickle.load(f)
            # Ignore data from predefined number of generations per seed
        if self.ap.n_generation <= IGNORE_FIRST_N_GENS:
            print('Skipping analysis - not enough generations run.')
            return
        cell_paths = self.ap.get_cells(
            generation=np.arange(IGNORE_FIRST_N_GENS, self.ap.n_generation), seed=SEED_RANGE,
            only_successful=True)

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

        # There are 4539 mRNA ids total w/ gene names
        cistron_id_to_gene_id = {
            cistron['id']: cistron['gene_id']
            for cistron in sim_data.process.transcription.cistron_data
        }

        # There are 4310 mRNA ids with associated protein/monomer ids
        protein_id_to_cistron_id = {
            protein['id']: protein['cistron_id']
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
            monomer_id_to_index[monomer_id] for monomer_id in monomers_of_interest
        ])

        # order cistrons in the order of monomers ids
        cistron_ids_in_order = np.array([
            protein_id_to_cistron_id[monomer_id] for monomer_id in monomers_of_interest
        ])

        gene_names_in_order = np.array([
            monomers_of_interest_name_dict[monomer_id] for monomer_id in monomers_of_interest
        ])

        # Get indices of cistron_ids_in_order
        mRNA_ids_indices = np.array([
            mRNA_id_to_index[cistron_id] for cistron_id
            in cistron_ids_in_order
        ])

        cell_paths = self.ap.get_cells(
            generation=np.arange(IGNORE_FIRST_N_GENS, self.ap.n_generation), seed=SEED_RANGE,
            only_successful=True)


        # should only be capturing the deltas from 0 to >0 transcripts, otherwise tracking transcript by second
        def count_peaks(time_series_data):
            # Convert counts to a boolean array: True if count > 0, False otherwise
            #is_present = time_series_data > 0
            # look at all new mRNA appearances

            # Use np.diff to find the difference between adjacent time steps (rows, axis=0).
            transition_deltas = np.diff(time_series_data.astype(int), axis=0)

            # Count only the '1's (the False -> True transitions) along the time axis (axis=0).
            on_event_count = (transition_deltas == 1).sum(axis=0)
            return on_event_count

        
        cistron_peaks_per_gen_full = read_stacked_columns(
            cell_paths, 
            'RNACounts', 
            'mRNA_cistron_counts',
            ignore_exception=True, 
            fun=count_peaks)[:, mRNA_ids_indices]

        tabel_cols = ['cell_id'] + monomers_of_interest
        # Write data to table so that the first col is the cell id and the rest are the counts per monomer
        with open(os.path.join(plotOutDir, plotOutFileName + '_count_transcipt_peaks.tsv'), 'w') as f:
            writer = csv.writer(f, delimiter='\t')
            writer.writerow(tabel_cols)

            for i in np.arange(0, len(cell_paths)):
                cell_id = cell_paths[i]
                counts_row = cistron_peaks_per_gen_full[i]
                full_row = [cell_id] + counts_row.tolist()
                writer.writerow(full_row)


if __name__ == '__main__':
    Plot().cli()