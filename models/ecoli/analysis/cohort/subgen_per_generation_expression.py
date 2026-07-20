"""
Per-generation expression audit (Task 4).

For every successful lineage (strict filter: completed every generation and no
cell hit the 180-min doubling cap), and for every generation of that lineage,
determine which protein-coding genes had at least one SUCCESSFUL (completed) mRNA
transcription event that generation. A gene counts as "expressed" for the whole
generation if it produced >= 1 completed transcript at any timestep during it.

Completed transcripts come from TranscriptElongationListener/
countRnaCistronSynthesized (attenuation-excluded), summed over the generation;
expressed := (sum > 0). This reads the pre-computed raw extraction rather than
simulation output, so run subgen_raw_extract.py first.

Outputs (to plotOutDir, prefixed by plotOutFileName):
  <name>.tsv          long: seed, generation, gene_name, cistron_id,
                      expressed (0/1), synth_count.
  <name>_metrics.tsv  per-gene summary across successful lineages/generations:
                      n_generations, n_generations_expressed, frac_expressed,
                      n_seeds, n_seeds_ever_expressed.
  <name>_per_generation_summary.tsv  per (seed, generation): n_genes_expressed,
                      n_genes_total, frac_genes_expressed.
"""

import csv
import os

import numpy as np

from models.ecoli.analysis import cohortAnalysisPlot
from models.ecoli.analysis.cohort import subgen_common as sc


class Plot(cohortAnalysisPlot.CohortAnalysisPlot):
	def do_plot(self, variantDir, plotOutDir, plotOutFileName, simDataFile,
			validationDataFile, metadata):
		# Load the raw extraction (errors clearly if it has not been run).
		seeds, generations, is_successful, gene_ids, synth = sc.load_raw_synth(
			plotOutDir)
		gene_ids_key, cistron_ids, _ = sc.load_raw_genes(plotOutDir)
		if gene_ids_key != gene_ids:
			print('WARNING: gene column order differs between synth matrix and '
				'gene key; falling back to the synth-matrix gene order.')
			# Map cistron ids by gene id where possible.
			gene_to_cistron = dict(zip(gene_ids_key, cistron_ids))
			cistron_ids = [gene_to_cistron.get(g, '') for g in gene_ids]

		# Restrict to strict-successful lineages.
		mask = is_successful
		n_total_cells = len(seeds)
		seeds = seeds[mask]
		generations = generations[mask]
		synth = synth[mask]
		n_cells, n_genes = synth.shape
		print('Successful cells: %d / %d ; genes: %d'
			% (n_cells, n_total_cells, n_genes))
		if n_cells == 0:
			print('No successful-lineage cells found. Skipping.')
			return

		expressed = synth > 0  # (n_cells, n_genes) bool

		prefix = os.path.join(plotOutDir, plotOutFileName)

		# --- Long per-(seed, generation, gene) table ---
		long_path = prefix + '.tsv'
		print('Writing %s (%d rows)' % (long_path, n_cells * n_genes))
		with open(long_path, 'w') as f:
			w = csv.writer(f, delimiter='\t')
			w.writerow(['seed', 'generation', 'gene_name', 'cistron_id',
				'expressed', 'synth_count'])
			for ci in range(n_cells):
				seed = int(seeds[ci])
				gen = int(generations[ci])
				exp_row = expressed[ci]
				synth_row = synth[ci]
				for gi in range(n_genes):
					w.writerow([seed, gen, gene_ids[gi], cistron_ids[gi],
						int(exp_row[gi]), int(synth_row[gi])])

		# --- Per-gene metrics across successful lineages/generations ---
		n_gens_expressed = expressed.sum(axis=0)             # over all cells
		frac_expressed = n_gens_expressed / n_cells
		# Per-seed "ever expressed": did the gene fire in >=1 gen of that seed?
		unique_seeds = sorted(set(int(s) for s in seeds))
		ever_by_seed = np.zeros(n_genes, dtype=int)
		for s in unique_seeds:
			srows = seeds == s
			ever_by_seed += expressed[srows].any(axis=0).astype(int)

		metrics_path = prefix + '_metrics.tsv'
		print('Writing %s' % metrics_path)
		with open(metrics_path, 'w') as f:
			w = csv.writer(f, delimiter='\t')
			w.writerow(['gene_name', 'cistron_id', 'n_generations',
				'n_generations_expressed', 'frac_generations_expressed',
				'n_seeds', 'n_seeds_ever_expressed'])
			for gi in range(n_genes):
				w.writerow([gene_ids[gi], cistron_ids[gi], n_cells,
					int(n_gens_expressed[gi]), '%.6g' % frac_expressed[gi],
					len(unique_seeds), int(ever_by_seed[gi])])

		# --- Per-(seed, generation) summary ---
		summary_path = prefix + '_per_generation_summary.tsv'
		print('Writing %s' % summary_path)
		genes_expressed_per_cell = expressed.sum(axis=1)
		with open(summary_path, 'w') as f:
			w = csv.writer(f, delimiter='\t')
			w.writerow(['seed', 'generation', 'n_genes_expressed',
				'n_genes_total', 'frac_genes_expressed'])
			for ci in range(n_cells):
				w.writerow([int(seeds[ci]), int(generations[ci]),
					int(genes_expressed_per_cell[ci]), n_genes,
					'%.6g' % (genes_expressed_per_cell[ci] / n_genes)])
		print('Done.')


if __name__ == '__main__':
	Plot().cli()
