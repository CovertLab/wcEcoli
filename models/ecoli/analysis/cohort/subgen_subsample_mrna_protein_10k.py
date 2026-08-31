"""
Pseudo-single-cell subsampling: ~10k random timepoints drawn
across strict-successful seeds, with mRNA, monomer, rRNA and tRNA counts at each.

The budget is split evenly over the seeds ACTUALLY sampled (sc.sample_per_seed),
so failed seeds no longer silently undersample the cohort.

No subgen definition is applied here -- this is data collection, and every gene
is reported.

Writes, all sharing the same (seed, time_step) rows in the same order:

  <name>.tsv                        seed, time_step, generation_start_time,
                                    then one column per gene holding TOTAL mRNA
                                    counts (full + partial), then the rRNA
                                    cistrons, then tRNA counts collapsed by
                                    amino acid.

  <name>_full_mrna.tsv              complete transcripts only, per gene.

  <name>_partial_mrna.tsv           nascent (incomplete) transcripts only.

  <name>_corresponding_monomers.tsv seed, time_step, generation,
                                    generation_start_time, dry_mass_fg, then one
                                    monomer count column per gene, in the same
                                    gene order as the mRNA tables.
  <name>_monomer_mws.tsv            monomer_id, mw_g_per_mol, in that same order.

On the mRNA split: RNACounts/mRNA_cistron_counts includes both partial and full
transcripts, so it is the TOTAL, and total == full + partial exactly (the
cistron-to-TU mapping is linear). Only two of the three are independent; all
three are written because the total is the historical column.
"""

import pickle
import os

import numpy as np

import csv

from models.ecoli.analysis import cohortAnalysisPlot
from models.ecoli.analysis.cohort import subgen_helper_functions as sc
from wholecell.analysis.analysis_tools import (
    read_stacked_bulk_molecules, read_stacked_columns)
from wholecell.io.tablereader import TableReader

IGNORE_FIRST_N_GENS = sc.IGNORE_FIRST_N_GENS
SEED_RANGE = sc.SEED_RANGE
TIMEPOINTS_TO_SAMPLE = sc.TIMEPOINTS_TO_SAMPLE
# Cells per block when reading the heavy per-timestep tables. Bounds memory:
# six columns are now sampled per seed instead of two.
BATCH_SIZE = 100

# only successful seeds are used for sampling to divide equal number of
# data points coming from each seed. 


def _sampled_rows(cell_paths, table, column, row_indices, offsets, n_cols,
        col_indices=None):
    """Read `column` at just the sampled rows, one BATCH_SIZE block at a time.

    Every table here is stacked with remove_first=True, so `offsets` -- built
    from Main/time row counts under that same flag -- maps a global sampled row
    index onto the right block and local row for any of them. A block that
    contains no sampled timestep is never read at all.
    """
    out = np.empty((len(row_indices), n_cols))
    for b in range(0, len(cell_paths), BATCH_SIZE):
        batch = cell_paths[b:b + BATCH_SIZE]
        row_start = offsets[b]
        row_end = offsets[b + len(batch)]
        mask = (row_indices >= row_start) & (row_indices < row_end)
        if not mask.any():
            continue
        local_idx = row_indices[mask] - row_start
        block = read_stacked_columns(
            batch, table, column, remove_first=True)
        if col_indices is None:
            out[mask] = block[local_idx].reshape(len(local_idx), n_cols)
        else:
            out[mask] = block[local_idx][:, col_indices]
    return out


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

        # Strict-successful seeds 
        success = sc.compute_seed_success(self.ap, self.ap.n_generation)
        seeds_to_sample = sorted(
            s for s in success['successful_seeds'] if s in set(SEED_RANGE.tolist()))
        sample_per_seed = sc.sample_per_seed(len(seeds_to_sample))
        print('Sampling %d timepoints from each of %d successful seeds '
            '(target %d total).'
            % (sample_per_seed, len(seeds_to_sample), TIMEPOINTS_TO_SAMPLE))
        if not seeds_to_sample:
            print('No successful seeds found. Skipping.')
            return

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

        # Get subcolumn for monomer IDs in monomer counts table
        monomer_counts_reader = TableReader(
            os.path.join(cell_paths[0], 'simOut', 'MonomerCounts'))
        monomer_ids_monomer_counts_table = monomer_counts_reader.readAttribute('monomerIds')
        monomer_counts_reader.close()
        
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
        
        # Get indeces of monomers in this subcolumn
        monomer_id_to_index = {
            monomer_id: i for (i, monomer_id)
            in enumerate(monomer_ids_monomer_counts_table)
            }
        
        monomer_indexes = np.array([monomer_id_to_index[monomer_id] for monomer_id in monomer_ids])

        # Molecular weights, in the same gene order as every count column, so the
        # sidecar joins to the monomer table by position as well as by id.
        id_to_mw = {m['id']: m['mw'] for m in sim_data.process.translation.monomer_data}
        monomer_mws = [id_to_mw[monomer_id] for monomer_id in monomer_ids]

        n_genes = len(gene_ids)
        total_seed_ids = []
        total_random_time_steps = []
        total_gen_start_times = []
        total_generation_ids = []
        total_dry_mass = []
        total_mRNA_counts = np.empty((0, n_genes), dtype=np.float64)
        total_full_mRNA_counts = np.empty((0, n_genes), dtype=np.float64)
        total_partial_mRNA_counts = np.empty((0, n_genes), dtype=np.float64)
        total_monomer_counts = np.empty((0, n_genes), dtype=np.float64)
        total_rRNA_counts = np.empty((0, len(rRNA_ids_counts_table)), dtype=np.float64)
        total_trna_counts = np.empty((0, len(aa_ids)), dtype=np.float64)

        # Seed the RNG once (not per seed) so each seed draws independently.
        np.random.seed(0)
        for seed in SEED_RANGE:
            cell_paths_per_seed = self.ap.get_cells(
                generation=np.arange(IGNORE_FIRST_N_GENS, self.ap.n_generation), seed=[seed],
                only_successful=True)

            if seed not in success['successful_seeds']:
                continue
            if len(cell_paths_per_seed) == 0:
                continue

            # Draw random timepoints and align them to generation starts 
            sample = sc.subsample_seed_timepoints(
                cell_paths_per_seed, sample_per_seed)
            if sample is None:
                continue
            random_time_indices = sample['time_indices']
            random_time_steps = sample['time_steps']
            aligned_start_times = sample['gen_start_times']

            # Map each sampled time step to its generation number.
            gen_indices_per_cell = [
                int(os.path.basename(os.path.dirname(cp))[-6:])
                for cp in cell_paths_per_seed
                ]
            generation_ids = np.array(gen_indices_per_cell)[sample['gen_index']]

            # Per-cell row counts (remove_first) give the row offsets used to map
            # global sampled indices back into per-batch reads of the heavy tables.
            per_cell_rows = read_stacked_columns(
                cell_paths_per_seed, 'Main', 'time', remove_first=True,
                fun=lambda x: len(x)).flatten().astype(int)
            cell_row_offsets = np.concatenate([[0], np.cumsum(per_cell_rows)])

            def sampled(table, column, n_cols, col_indices=None):
                return _sampled_rows(
                    cell_paths_per_seed, table, column, random_time_indices,
                    cell_row_offsets, n_cols, col_indices)

            # Get counts of mRNAs for each gene across random timepoints. The
            # three RNA columns share the mRNA_cistron_ids subcolumn key, so one
            # index map serves all of them.
            mRNA_counts = sampled(
                'RNACounts', 'mRNA_cistron_counts', n_genes, mRNA_cistron_indices)
            full_mRNA_counts = sampled(
                'RNACounts', 'full_mRNA_cistron_counts', n_genes, mRNA_cistron_indices)
            partial_mRNA_counts = sampled(
                'RNACounts', 'partial_mRNA_cistron_counts', n_genes, mRNA_cistron_indices)

            # Get monomer counts for each gene at the same timepoints
            # (read from this seed's cells, not the whole cohort).
            monomer_counts = sampled(
                'MonomerCounts', 'monomerCounts', n_genes, monomer_indexes)

            # Dry mass at the same timepoints, in femtograms (listener unit).
            dry_mass = sampled('Mass', 'dryMass', 1).flatten()

            # Partial rRNAs count, the model doesn't seem to have just "rRNA cistron counts", get counts for each rRNA gene
            partial_rRNA_cistron_counts = sampled(
                'RNACounts', 'partial_rRNA_cistron_counts',
                len(rRNA_ids_counts_table))

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
            total_generation_ids = total_generation_ids + generation_ids.tolist()
            total_dry_mass = total_dry_mass + dry_mass.tolist()
            total_mRNA_counts = np.vstack([total_mRNA_counts, mRNA_counts])
            total_full_mRNA_counts = np.vstack([total_full_mRNA_counts, full_mRNA_counts])
            total_partial_mRNA_counts = np.vstack([total_partial_mRNA_counts, partial_mRNA_counts])
            total_monomer_counts = np.vstack([total_monomer_counts, monomer_counts])
            total_rRNA_counts = np.vstack([total_rRNA_counts, partial_rRNA_cistron_counts])
            total_trna_counts =  np.vstack([total_trna_counts, full_trna_counts])

        n_rows = len(total_random_time_steps)
        table_cols = ['seed', 'time_step', 'generation_start_time'] + gene_ids + rRNA_ids_counts_table + aa_ids

        # Write data to table
        with open(os.path.join(plotOutDir, plotOutFileName + '.tsv'), 'w') as f:
            writer = csv.writer(f, delimiter='\t')
            writer.writerow(table_cols)
            for i in np.arange(0, n_rows):
                seed = total_seed_ids[i]
                time_step = total_random_time_steps[i]
                gen_start = total_gen_start_times[i]
                mRNA_counts_row = total_mRNA_counts[i].tolist()
                rRNA_counts_row = total_rRNA_counts[i].tolist()
                trna_counts_row = total_trna_counts[i].tolist()
                counts_row = [seed]+ [time_step]+ [gen_start] + mRNA_counts_row + rRNA_counts_row + trna_counts_row
                writer.writerow(counts_row)

        # The full/partial split, as sidecars rather than extra column blocks:
        # the main table is already ~4300 columns wide, and these rows join to it
        # on (seed, time_step).
        for suffix, matrix in (
                ('_full_mrna.tsv', total_full_mRNA_counts),
                ('_partial_mrna.tsv', total_partial_mRNA_counts)):
            path = os.path.join(plotOutDir, plotOutFileName + suffix)
            print('Writing %s' % path)
            with open(path, 'w') as f:
                writer = csv.writer(f, delimiter='\t')
                writer.writerow(['seed', 'time_step'] + gene_ids)
                for i in np.arange(0, n_rows):
                    writer.writerow(
                        [total_seed_ids[i], total_random_time_steps[i]]
                        + matrix[i].tolist())

        monomer_table_cols = (
            ['seed', 'time_step', 'generation', 'generation_start_time',
                'dry_mass_fg']
            + monomer_ids)
        with open(os.path.join(plotOutDir, plotOutFileName + '_corresponding_monomers.tsv'), 'w') as f:
            writer = csv.writer(f, delimiter='\t')
            writer.writerow(monomer_table_cols)
            for i in np.arange(0, n_rows):
                seed = total_seed_ids[i]
                time_step = total_random_time_steps[i]
                generation = total_generation_ids[i]
                gen_start = total_gen_start_times[i]
                dry_mass_row = total_dry_mass[i]
                monomer_counts_row = total_monomer_counts[i].tolist()
                counts_row2 = ([seed] + [time_step] + [generation] + [gen_start]
                    + [dry_mass_row] + monomer_counts_row)
                writer.writerow(counts_row2)

        mws_path = os.path.join(plotOutDir, plotOutFileName + '_monomer_mws.tsv')
        print('Writing %s' % mws_path)
        with open(mws_path, 'w') as f:
            writer = csv.writer(f, delimiter='\t')
            writer.writerow(['monomer_id', 'mw_g_per_mol'])
            for monomer_id, mw in zip(monomer_ids, monomer_mws):
                writer.writerow([monomer_id, mw])

        sim_metadata_path, sim_metadata = sc.load_sim_metadata(variantDir)
        sc.write_run_metadata(
            os.path.join(plotOutDir, plotOutFileName + '_run_metadata.json'),
            script=os.path.basename(__file__),
            parameters={
                'ignore_first_n_gens': IGNORE_FIRST_N_GENS,
                'max_doubling_min': sc.MAX_DOUBLING_MIN,
                'timepoints_to_sample': TIMEPOINTS_TO_SAMPLE,
                'sample_per_seed': sample_per_seed,
                'batch_size': BATCH_SIZE,
                'numpy_seed': 0,
                },
            sim_metadata_path=sim_metadata_path,
            sim_metadata=sim_metadata,
            extra={
                'seeds': {
                    'n_successful': len(seeds_to_sample),
                    'successful_seeds': seeds_to_sample,
                    },
                'timepoints': {'n_sampled': n_rows},
                })

if __name__ == '__main__':
	Plot().cli()
