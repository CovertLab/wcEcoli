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
from wholecell.utils.sparkline import whitePadSparklineAxis

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

IGNORE_FIRST_N_GENS = 2
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
            generation=np.arange(IGNORE_FIRST_N_GENS, self.ap.n_generation),
            only_successful=True)

        print('Analyzing %d cells...' % len(cell_paths))

        # There are 4346 mRNA ids with counts
        RNA_reader = TableReader(
            os.path.join(cell_paths[0], 'simOut', 'RNACounts'))
        mRNA_cistron_ids_counts_table = RNA_reader.readAttribute('mRNA_cistron_ids')
        RNA_reader.close()

        # 4346
        simulatedSynthProb = read_stacked_columns(
            cell_paths, 'RnaSynthProb', 'actual_rna_synth_prob_per_cistron',
            remove_first=True).mean(axis=0) [mRNA_cistron_indexes]

        # 4346, order of mrnas is the same as in mRnaIndexes/mRnaIds
        mRNA_exists_in_gen = read_stacked_columns(
            cell_paths, 'RNACounts', 'mRNA_cistron_counts', 
            ignore_exception=True, fun=lambda x: x.sum(axis=0) > 0)

        # Divide by total number of cells to get probability
        p_mRNA_exists_in_gen = (
            mRNA_exists_in_gen.sum(axis=0) / mRNA_exists_in_gen.shape[0]
        )

        # Write data to table
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
