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

IGNORE_FIRST_N_GENS = sc.IGNORE_FIRST_N_GENS
SEED_RANGE = sc.SEED_RANGE
TIMEPOINTS_TO_SAMPLE = 10000
SAMPLE_PER_SEED = TIMEPOINTS_TO_SAMPLE // len(SEED_RANGE)


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

        # Strict-successful lineages (completed every generation and no cell at
        # the 180-min doubling cap).
        success = sc.compute_lineage_success(self.ap, self.ap.n_generation)
        # Seed the RNG once for reproducible timepoint subsampling.
        np.random.seed(0)

        # Load from sim_data
        transcription = sim_data.process.transcription
        cistron_data = transcription.cistron_data
        # 4539 cistrons total
        cistron_ids = cistron_data['id']
        uncharged_trna_names = transcription.uncharged_trna_names
        charged_trna_names = transcription.charged_trna_names
        aa_from_trna = transcription.aa_from_trna  # (n_aa, n_trna) 0/1 map; sc.collapse_trna_to_aa transposes it
        molecule_groups = sim_data.molecule_groups
        aa_ids = molecule_groups.amino_acids

        # Filter list for cistron IDs with associated protein ids
        cistron_id_to_protein_id = {
            protein['cistron_id']: protein['id']
            for protein in sim_data.process.translation.monomer_data
            }
        # 4310 cistrons with associated protein ids
        mRNA_cistron_ids = [
            cistron_id for cistron_id in cistron_ids
            if cistron_id in cistron_id_to_protein_id]
        RNA_reader = TableReader(
            os.path.join(cell_paths[0], 'simOut', 'RNACounts'))
        # 4346 mRNA cistron IDs with counts
        mRNA_ids_counts_table = RNA_reader.readAttribute('mRNA_cistron_ids')
        # 22 cistrons corresponding to rRNA genes
        rRNA_ids_counts_table = RNA_reader.readAttribute('rRNA_cistron_ids')
        RNA_reader.close()
        
        # Dictionary with index corresponding to each mRNA cistron with counts
        mRNA_cistron_id_to_index = {
            cistron_id: i for (i, cistron_id)
            in enumerate(mRNA_ids_counts_table)
            }
        # Get indices of cistrons with associated protein ids in the counts table
        mRNA_cistron_indices = np.array([
            mRNA_cistron_id_to_index[cistron_id] for cistron_id
            in mRNA_cistron_ids
            ]) 
        # Get IDs of associated monomers and genes
        monomer_ids = [
            cistron_id_to_protein_id.get(cistron_id, None)
            for cistron_id in mRNA_cistron_ids]
        cistron_id_to_gene_id = {
            cistron['id']: cistron['gene_id'] for cistron in cistron_data
            }
        gene_ids = [
            cistron_id_to_gene_id[cistron_id]
            for cistron_id in mRNA_cistron_ids]
        
        total_seed_ids = []
        total_random_time_steps = []
        total_gen_start_times = []
        total_mRNA_counts = np.empty((0, len(gene_ids)), dtype=np.float64)
        total_rRNA_counts = np.empty((0, len(rRNA_ids_counts_table)), dtype=np.float64)
        total_trna_counts = np.empty((0, len(aa_ids)), dtype=np.float64)

        for seed in SEED_RANGE:
            cell_paths_per_seed = self.ap.get_cells(
                generation=np.arange(IGNORE_FIRST_N_GENS, self.ap.n_generation), seed=[seed],
                only_successful=True)
            
            if len(cell_paths_per_seed) == 0:
                continue

            if seed not in success['successful_seeds']:
                continue

            # Draw random timepoints and align them to generation starts (shared
            # helper: clamps the draw to the available timesteps).
            sample = sc.subsample_seed_timepoints(cell_paths_per_seed, SAMPLE_PER_SEED)
            if sample is None:
                continue
            random_time_indices = sample['time_indices']
            random_time_steps = sample['time_steps']
            aligned_start_times = sample['gen_start_times']
            time_step_rows = random_time_indices[:, np.newaxis]

            # Get counts of mRNAs for each gene across random timepoints
            mRNA_counts = read_stacked_columns(
                cell_paths_per_seed, 'RNACounts', 'mRNA_cistron_counts',
                remove_first=True)[time_step_rows, mRNA_cistron_indices]
            # Partial rRNAs count, the model doesn't seem to have just "rRNA cistron counts", get counts for each rRNA gene
            partial_rRNA_cistron_counts = read_stacked_columns(
                cell_paths_per_seed,'RNACounts', 'partial_rRNA_cistron_counts', remove_first=True)[random_time_indices]

            # Get counts of tRNAs by amino acid type
            (uncharged_trna_counts, charged_trna_counts, ) = read_stacked_bulk_molecules(
                cell_paths_per_seed, (uncharged_trna_names, charged_trna_names, ),
                remove_first=True)
            full_trna_counts = sc.collapse_trna_to_aa(
                charged_trna_counts[random_time_indices],
                uncharged_trna_counts[random_time_indices], aa_from_trna)
            total_seed_ids = total_seed_ids + [seed] * len(random_time_steps)
            total_random_time_steps = total_random_time_steps + random_time_steps.tolist()
            total_gen_start_times = total_gen_start_times + aligned_start_times.tolist()
            total_mRNA_counts = np.vstack([total_mRNA_counts, mRNA_counts])
            total_rRNA_counts = np.vstack([total_rRNA_counts, partial_rRNA_cistron_counts])
            total_trna_counts =  np.vstack([total_trna_counts, full_trna_counts])
        table_cols = ['seed', 'time_step', 'generation_start_time'] + gene_ids + rRNA_ids_counts_table + aa_ids

        # Write data to table
        with open(os.path.join(plotOutDir, plotOutFileName + '.tsv'), 'w') as f:
            writer = csv.writer(f, delimiter='\t')
            writer.writerow(table_cols)
            for i in np.arange(0, len(total_random_time_steps)):
                seed = total_seed_ids[i]
                time_step = total_random_time_steps[i]
                gen_start = total_gen_start_times[i]
                mRNA_counts_row = total_mRNA_counts[i].tolist()
                rRNA_counts_row = total_rRNA_counts[i].tolist()
                trna_counts_row = total_trna_counts[i].tolist()
                counts_row = [seed]+ [time_step]+ [gen_start] + mRNA_counts_row + rRNA_counts_row + trna_counts_row
                writer.writerow(counts_row)
                #import ipdb; ipdb.set_trace()

if __name__ == '__main__':
	Plot().cli()
