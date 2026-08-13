"""
Per-cell net change in monomer counts for  10 curated genes.

Writes <plotOutFileName>_new_monomers_per_gen.tsv: one row per cell generation,
first column the cell path, then one column per curated monomer holding
monomerCounts[last] - monomerCounts[first] for that generation.
"""

import os

# noinspection PyUnresolvedReferences
import numpy as np

import csv


from models.ecoli.analysis import cohortAnalysisPlot
from models.ecoli.analysis.cohort import subgen_common as sc
from wholecell.analysis.analysis_tools import read_stacked_columns
from wholecell.io.tablereader import TableReader

IGNORE_FIRST_N_GENS = sc.IGNORE_FIRST_N_GENS
SEED_RANGE = sc.SEED_RANGE
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
        # Ignore data from predefined number of generations per seed
        if self.ap.n_generation <= IGNORE_FIRST_N_GENS:
            print('Skipping analysis - not enough generations run.')
            return
        cell_paths = self.ap.get_cells(
            generation=np.arange(IGNORE_FIRST_N_GENS, self.ap.n_generation), seed=SEED_RANGE,
            only_successful=True)

        print('Analyzing %d cells...' % len(cell_paths))

        # Restrict to strict-successful lineages (completed every generation and
        # no cell at the 180-min doubling cap).
        success = sc.compute_lineage_success(self.ap, self.ap.n_generation)
        cell_paths = sc.filter_cells_to_successful(
            cell_paths, success['successful_seeds'])
        print('Analyzing %d cells from successful lineages...' % len(cell_paths))
        if len(cell_paths) == 0:
            print('No successful-lineage cells found. Skipping.')
            return

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

        def net_new_monomers(time_series_data):
            """Net change in monomer count over the cell generation: last - first.

            Equal to the old np.diff(...).sum(axis=0), written explicitly. Row 0 is the
            daughter's inherited count (Listener.initialUpdate runs before the logger's
            first write, so it is a real state, not the allocate() zeros) and the last row
            is the pre-division count, so this is a NET quantity: synthesis minus
            degradation, which can be negative. Under balanced growth it is close to the
            inherited count, i.e. about half the final count -- it is not gross synthesis.

            There is no per-monomer count of *completed* proteins in the listeners
            (RibosomeData/didTerminate and actualElongations are scalars over all
            monomers). The nearest gross per-monomer quantity is
            RibosomeData/ribosome_init_event_per_monomer (translation initiation events,
            monomerIds order). This script deliberately reproduces the original net-delta
            quantity that protein_distribution_new_monomers_per_gen.tsv has always held.

            Cast to int64 before subtracting: monomerCounts is int64 today (BulkMolecules
            uses BulkObjectsContainer's default dtype), but on an unsigned dtype the
            negative values this quantity depends on would silently wrap to huge positives.
            """
            return (time_series_data[-1].astype(np.int64)
                - time_series_data[0].astype(np.int64))

        # Read one cell at a time so each delta row stays paired with its own
        # cell path. (A single read_stacked_columns(ignore_exception=True) call
        # silently DROPS unreadable cells, which would shift every later row off
        # its cell_id label and eventually IndexError.)
        delta_rows = []
        kept_cell_ids = []
        for cell_path in cell_paths:
            try:
                cell_delta = read_stacked_columns(
                    [cell_path], 'MonomerCounts', 'monomerCounts',
                    fun=net_new_monomers)
            except Exception as e:
                print('  Warning: could not read %s: %s' % (cell_path, e))
                continue
            delta_rows.append(cell_delta[0][monomer_indices])
            kept_cell_ids.append(cell_path)

        tabel_cols = ['cell_id'] + monomers_of_interest
        # Write data to table so that the first col is the cell id and the rest are the counts per monomer
        with open(os.path.join(plotOutDir, plotOutFileName + '_new_monomers_per_gen.tsv'), 'w') as f:
            writer = csv.writer(f, delimiter='\t')
            writer.writerow(tabel_cols)

            for cell_id, counts_row in zip(kept_cell_ids, delta_rows):
                full_row = [cell_id] + counts_row.tolist()
                writer.writerow(full_row)


if __name__ == '__main__':
    Plot().cli()
