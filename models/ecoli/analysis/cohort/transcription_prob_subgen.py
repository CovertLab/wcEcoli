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
        # Get mRNA data
        transcription = sim_data.process.transcription
        # 3269
        rnaIds = transcription.rna_data["id"]
        isMRna = transcription.rna_data['is_mRNA']
        # 3126
        mRnaIndexes = np.where(isMRna)[0]
        # 3126
        mRnaIds = np.array([rnaIds[x] for x in mRnaIndexes])

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
        #3126
        mRNA_ids_counts_table = RNA_reader.readAttribute('mRNA_ids')
        RNA_reader.close()

        # dictionary to map proteins of interest to cistron ids
        protein_id_to_cistron_id = {
            protein['id']: protein['cistron_id']
            for protein in sim_data.process.translation.monomer_data
        }

        # corresponding cistron ids in order of monomers of interest
        cistron_ids_in_order = np.array([
            protein_id_to_cistron_id[monomer_id] for monomer_id in monomers_of_interest
        ])


        # corresponding TU IDs in order of monomers of interest
        cistron_TU_index_dict = {
            cistron_id: transcription.cistron_id_to_rna_indexes(cistron_id) 
            for cistron_id in cistron_ids_in_order
            }

        all_TU_ids = np.array(mRNA_ids_counts_table)

        cistron_TU_ids_dict = {
            cistron_id: all_TU_ids[cistron_TU_index_dict[cistron_id]]
            for cistron_id in cistron_ids_in_order
        }
        
        # 3269
        simulatedSynthProb = read_stacked_columns(
            cell_paths, 'RnaSynthProb', 'actual_rna_synth_prob',
            remove_first=True).mean(axis=0)[mRnaIndexes]
        # 3126, order of mrnas is the same as in mRnaIndexes/mRnaIds
        mRNACounts_sumOverTime = read_stacked_columns(
            cell_paths, 'RNACounts', 'full_mRNA_counts',
            ignore_exception=True).sum(axis = 0)

        mRnasTranscribed = np.array([x != 0 for x in mRNACounts_sumOverTime])

        # Write data to table
        with open(os.path.join(plotOutDir, plotOutFileName + '.tsv'), 'w') as f:
            writer = csv.writer(f, delimiter='\t')
            writer.writerow([
                'mRNA_name',
                'actual_rna_synth_prob',
                'mRNA_counts_aggregated',
                'transcribed'
            ])

            for i in range(len(mRnaIds)):
                writer.writerow([
                    mRnaIds[i], simulatedSynthProb[i], mRNACounts_sumOverTime[i],
                    mRnasTranscribed[i],
                ])

if __name__ == '__main__':
    Plot().cli()
