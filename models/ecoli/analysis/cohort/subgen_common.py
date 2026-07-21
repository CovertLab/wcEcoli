"""
Shared helpers for the subgenerational-expression analyses.

This module is NOT an analysis Plot (it has no `Plot` class and is not listed in
`__init__.py`), so the analysis runner never imports or runs it directly. It
consolidates logic that was previously duplicated across
`subgen_expression_definitions.py` and `subgen_definition5_lineage_ci.py`, and
provides the canonical implementation of:

  * the gene set (protein-coding mRNA cistrons, in cistron order),
  * the STRICT "successful lineage" filter (completed every generation AND no
    cell hit the 180-min doubling cap),
  * Definition 5 classifiers -- Form B (per-lineage confidence interval, the
    canonical subgen label) and Form A (pooled-cell threshold),
  * small TSV IO helpers for the raw-extraction pipeline.

Definition 5 is the mean number of *completed* mRNA transcripts per generation,
per gene, read from TranscriptElongationListener/countRnaCistronSynthesized.
"""

import csv
import json
import os
import subprocess

import numpy as np

from wholecell.io.tablereader import TableReader
from wholecell.utils import constants


# --- Shared constants -------------------------------------------------------

# Skip the first few generations (transient from initial conditions).
IGNORE_FIRST_N_GENS = 8
# Default cohort seed range. The sims run 0..127 seeds; scripts pass this to
# ap.get_cells(seed=...) (absent seeds are tolerated) and then gate on the strict
# successful-lineage set, so no seed is silently excluded by a too-small range.
N_SEEDS = 128
SEED_RANGE = np.arange(0, N_SEEDS)
# Max simulation length is lengthSec = 3*60*60 s = 180 min; a generation whose
# doubling time reaches this never divided within the simulation window.
MAX_DOUBLING_MIN = 180.0
DOUBLING_AT_MAX_TOL = 0.01  # treat >= 179.99 min as "hit the cap"
# 95% confidence interval from the standard error of the mean (normal approx;
# t vs z is negligible at ~100 lineages).
CI_Z = 1.96
CONFIDENCE = 0.95

SYNTH_TABLE = 'TranscriptElongationListener'
SYNTH_COLUMN = 'countRnaCistronSynthesized'

# Definition-5 categories and a validated categorical palette (dataviz
# validator, light + dark). Grey is the neutral null category.
CATEGORIES = ['subgen', 'possibly_subgen', 'not_subgen', 'never_expressed']
PALETTE = {
	'subgen': '#1667B8',
	'possibly_subgen': '#109C9C',
	'not_subgen': '#C7362F',
	'never_expressed': '#8A828C',
	}
CAT_LABEL = {
	'subgen': 'Subgenerational (CI < 1)',
	'possibly_subgen': 'Possibly subgen (CI includes 1)',
	'not_subgen': 'Not subgen (CI > 1)',
	'never_expressed': 'Never expressed (mean = 0)',
	}
# Neutral figure colors (shared by the def-5 figure scripts).
INK = '#1b2530'
MUTED = '#5b6672'
GRID = '#e4e8ec'
SURF = '#fcfcfb'


# --- Path / provenance helpers ---------------------------------------------

def parse_cell_id(cell_path):
	"""Split a cell path into (seed_int, generation_int).

	Cell paths look like .../<variant>/<seed>/generation_<gen>/<daughter>.
	Returns (-1, -1) if the path does not match.
	"""
	parts = cell_path.rstrip(os.sep).split(os.sep)
	for i, part in enumerate(parts):
		if part.startswith('generation_'):
			seed = int(parts[i - 1]) if i > 0 and parts[i - 1].isdigit() else -1
			return seed, int(part.split('_')[1])
	return -1, -1


def git_info(repo_dir):
	"""Return the current git hash, branch, and dirty flag for repo_dir."""
	def run(args):
		return subprocess.check_output(
			['git', '-C', repo_dir] + args,
			stderr=subprocess.DEVNULL).decode().strip()
	try:
		return {
			'git_hash': run(['rev-parse', 'HEAD']),
			'git_branch': run(['rev-parse', '--abbrev-ref', 'HEAD']),
			'git_dirty': bool(run(['status', '--porcelain'])),
			}
	except Exception as e:
		return {'git_hash': None, 'git_branch': None, 'git_dirty': None,
			'error': str(e)}


def load_sim_metadata(variant_dir):
	"""Load the simulation's metadata.json (git hash, run time, options).

	The sim-level metadata directory sits one level above the variant directory;
	fall back to a metadata directory inside the variant directory.
	"""
	candidates = [
		os.path.join(os.path.dirname(variant_dir),
			constants.METADATA_DIR, constants.JSON_METADATA_FILE),
		os.path.join(variant_dir,
			constants.METADATA_DIR, constants.JSON_METADATA_FILE),
		]
	for path in candidates:
		if os.path.isfile(path):
			with open(path) as f:
				return path, json.load(f)
	return None, {}


# --- Gene set ---------------------------------------------------------------

def get_mrna_gene_set(sim_data):
	"""Canonical gene set: protein-coding mRNA cistrons, in cistron order.

	Returns (mRNA_cistron_ids, monomer_ids, gene_ids) as parallel lists, one
	entry per mRNA cistron that has an associated protein/monomer.
	"""
	cistron_data = sim_data.process.transcription.cistron_data
	cistron_id_to_protein_id = {
		protein['cistron_id']: protein['id']
		for protein in sim_data.process.translation.monomer_data
		}
	cistron_id_to_gene_id = {
		cistron['id']: cistron['gene_id'] for cistron in cistron_data
		}
	mRNA_cistron_ids = [
		cistron_id for cistron_id in cistron_data['id']
		if cistron_id in cistron_id_to_protein_id]
	monomer_ids = [
		cistron_id_to_protein_id[cistron_id]
		for cistron_id in mRNA_cistron_ids]
	gene_ids = [
		cistron_id_to_gene_id[cistron_id]
		for cistron_id in mRNA_cistron_ids]
	return mRNA_cistron_ids, monomer_ids, gene_ids


def cistron_index_map(cell_path, table, ids):
	"""Map `ids` into the subcolumn order of `table`'s `cistron_ids` attribute.

	Returns (True, index_array) on success, else (False, None). Used to align
	our gene list with a listener's on-disk column order.
	"""
	try:
		reader = TableReader(os.path.join(cell_path, 'simOut', table))
		full_ids = reader.readAttribute('cistron_ids')
		id_to_index = {cid: i for i, cid in enumerate(full_ids)}
		return True, np.array([id_to_index[c] for c in ids])
	except Exception:
		return False, None


# --- Strict "successful lineage" filter ------------------------------------

def compute_lineage_success(ap, n_generation, total_init_sims=None):
	"""Compute the STRICT successful-lineage classification for a cohort.

	A lineage (seed) is "successful" iff it completed EVERY generation AND no
	cell hit the 180-minute doubling cap. This mirrors the rule in
	subgen_definition5_lineage_ci.py (stricter than
	subgen_expression_definitions_complete.py, which only checks that the last
	successful generation is the final one).

	Returns a dict with:
	  seeds            sorted list of seeds that produced any directory
	  doubling         {seed: np.array(n_generation)} doubling time (min), -1 if
	                   the generation did not run / could not be read
	  n_at_180         {seed: int} number of cells at the 180-min cap
	  gens_at_180      {seed: list[int]} which generations hit the cap
	  successful_gens  {seed: set[int]} generations that completed successfully
	  completed_all    {seed: bool} completed every generation
	  in_successful    {seed: bool} strict successful flag
	  successful_seeds set[int] of seeds with in_successful True
	  all_seed_ids     sorted seeds incl. never-ran seeds (if total_init_sims set)
	"""
	seeds = sorted(int(s) for s in ap.get_seeds())
	all_seed_ids = sorted(set(seeds) | set(range(total_init_sims))) \
		if total_init_sims else list(seeds)

	doubling = {s: -np.ones(n_generation) for s in all_seed_ids}
	successful_gens = {s: set() for s in all_seed_ids}
	n_at_180 = {s: 0 for s in all_seed_ids}
	gens_at_180 = {s: [] for s in all_seed_ids}

	for s in seeds:
		for cell_path in ap.get_cells(seed=[s], only_successful=False):
			_, gen = parse_cell_id(cell_path)
			if 0 <= gen < n_generation:
				try:
					time = TableReader(
						os.path.join(cell_path, 'simOut', 'Main')
						).readColumn('time')
					doubling[s][gen] = (time[-1] - time[0]) / 60.0
				except Exception:
					pass
		at_max = np.where(
			doubling[s] >= MAX_DOUBLING_MIN - DOUBLING_AT_MAX_TOL)[0]
		n_at_180[s] = int(at_max.size)
		gens_at_180[s] = at_max.tolist()
		successful_gens[s] = set(
			parse_cell_id(cp)[1]
			for cp in ap.get_cells(seed=[s], only_successful=True))

	completed_all = {
		s: set(range(n_generation)).issubset(successful_gens[s])
		for s in all_seed_ids}
	in_successful = {
		s: (completed_all[s] and n_at_180[s] == 0) for s in all_seed_ids}
	successful_seeds = {s for s in all_seed_ids if in_successful[s]}

	return {
		'seeds': seeds,
		'all_seed_ids': all_seed_ids,
		'doubling': doubling,
		'n_at_180': n_at_180,
		'gens_at_180': gens_at_180,
		'successful_gens': successful_gens,
		'completed_all': completed_all,
		'in_successful': in_successful,
		'successful_seeds': successful_seeds,
		}


def filter_cells_to_successful(cell_paths, successful_seeds):
	"""Keep only cells whose seed is in the strict-successful set."""
	return [cp for cp in cell_paths
		if parse_cell_id(cp)[0] in successful_seeds]


# --- Definition-5 classifiers ----------------------------------------------

def _classify_ci(mean, ci_low, ci_high):
	"""Assign each gene to a Def-5 category from its mean and CI (Form B)."""
	n = len(mean)
	out = np.empty(n, dtype=object)
	for i in range(n):
		if mean[i] == 0:
			out[i] = 'never_expressed'
		elif ci_high[i] < 1:
			out[i] = 'subgen'
		elif ci_low[i] > 1:
			out[i] = 'not_subgen'
		else:
			out[i] = 'possibly_subgen'
	return out


def classify_def5_ci(lambda_matrix, row_indices, n_genes):
	"""Form B: per-gene mean/std/se/CI/category over selected lineage rows.

	`lambda_matrix` is (n_lineages, n_genes) of per-lineage Def-5 rates.
	`row_indices` selects which lineage rows to include (e.g. successful ones).
	The 95% CI is the normal-approx SE of the mean; a gene is `subgen` iff its
	CI upper bound is below 1 completed transcript per generation.
	"""
	rows = lambda_matrix[row_indices, :]
	n = rows.shape[0]
	mean = rows.mean(axis=0) if n else np.full(n_genes, np.nan)
	std = rows.std(axis=0, ddof=1) if n >= 2 else np.zeros(n_genes)
	se = std / np.sqrt(n) if n else np.full(n_genes, np.nan)
	ci_low = np.maximum(0.0, mean - CI_Z * se)
	ci_high = mean + CI_Z * se
	cat = _classify_ci(mean, ci_low, ci_high)
	return {'n': n, 'mean': mean, 'std': std, 'se': se,
		'ci_low': ci_low, 'ci_high': ci_high, 'cat': cat}


def category_counts(cat):
	"""Count genes in each Def-5 category."""
	return {c: int(np.sum(cat == c)) for c in CATEGORIES}


def build_lineage_lambda(seeds_arr, synth_matrix, restrict_seeds=None):
	"""Collapse a per-cell synth matrix to per-lineage Def-5 rates (Form B input).

	`seeds_arr` is (n_cells,) of seed ints; `synth_matrix` is (n_cells, n_genes)
	of completed-transcript counts per cell (already burned-in). For each seed,
	the lineage rate is the mean over that seed's cells. Returns
	(lambda_matrix (n_lineages, n_genes), lineage_seeds list).
	"""
	seeds_arr = np.asarray(seeds_arr)
	unique_seeds = sorted(set(int(s) for s in seeds_arr))
	rows = []
	lineage_seeds = []
	for s in unique_seeds:
		if restrict_seeds is not None and s not in restrict_seeds:
			continue
		mask = seeds_arr == s
		if not mask.any():
			continue
		rows.append(synth_matrix[mask].mean(axis=0))
		lineage_seeds.append(s)
	lambda_matrix = np.array(rows) if rows \
		else np.zeros((0, synth_matrix.shape[1]))
	return lambda_matrix, lineage_seeds


def classify_def5_threshold(values, threshold=1.0):
	"""Form A: pooled-cell subgen mask, subgen iff 0 < value < threshold.

	`values` is the per-gene pooled mean of completed transcripts per cell.
	"""
	values = np.asarray(values, dtype=float)
	return (values > 0) & (values < threshold)


def pooled_mean(synth_matrix, row_mask=None):
	"""Form A pooled value: cell-weighted mean completed transcripts per gene."""
	if row_mask is not None:
		synth_matrix = synth_matrix[row_mask]
	if synth_matrix.shape[0] == 0:
		return np.full(synth_matrix.shape[1], np.nan)
	return synth_matrix.mean(axis=0)


# --- Timepoint-subsampling helpers -----------------------------------------

def subsample_seed_timepoints(cell_paths_per_seed, sample_per_seed):
	"""Draw up to `sample_per_seed` random timesteps (no replacement) from one
	lineage's cells and align each sampled timestep to its generation.

	Reads only Main/time. The caller MUST seed numpy's RNG once before the loop
	over seeds (``np.random.seed(0)``) so the subsample is reproducible; this
	helper does not reseed. Returns None if the lineage has no timesteps.

	The returned ``time_indices`` are row indices into this seed's
	``read_stacked_columns(..., remove_first=True)`` output for any per-timestep
	table (RNACounts, MonomerCounts, Mass, ...), so callers read their own count
	tables with ``remove_first=True`` and index with ``time_indices``.

	Returns dict:
	  time_indices     (S,) row indices into remove_first-stacked per-timestep tables
	  time_steps       (S,) sampled absolute time values (seconds)
	  gen_start_times  (S,) generation start time (s) for each sampled timestep
	  gen_index        (S,) index of the containing cell/generation within
	                   cell_paths_per_seed
	"""
	from wholecell.analysis.analysis_tools import read_stacked_columns
	time = read_stacked_columns(
		cell_paths_per_seed, 'Main', 'time', remove_first=True).flatten()
	if time.size == 0:
		return None
	gen_starts = read_stacked_columns(
		cell_paths_per_seed, 'Main', 'time', remove_first=True,
		fun=lambda x: x[0]).flatten()
	n = min(sample_per_seed, time.size)
	idx = np.random.choice(time.size, size=n, replace=False)
	steps = time[idx]
	# side='right' minus 1 -> the cell whose start time is <= the sample time,
	# i.e. the generation containing that timestep (time is monotonic across a
	# lineage, and gen_starts is increasing).
	gen_index = np.clip(
		np.searchsorted(gen_starts, steps, side='right') - 1, 0, None)
	return {
		'time_indices': idx,
		'time_steps': steps,
		'gen_start_times': gen_starts[gen_index],
		'gen_index': gen_index,
		}


def collapse_trna_to_aa(charged_counts, uncharged_counts, aa_from_trna):
	"""Sum charged+uncharged tRNA counts per amino acid.

	`aa_from_trna` is sim_data's transcription.aa_from_trna, a (n_aa, n_trna)
	0/1 map (pass it as-is, not transposed). `charged_counts`/`uncharged_counts`
	are (rows, n_trna) in uncharged_trna_names/charged_trna_names order. Returns
	(rows, n_aa).
	"""
	m = aa_from_trna.T
	return charged_counts @ m + uncharged_counts @ m


# --- Raw-extraction TSV IO --------------------------------------------------

def write_per_cell_matrix(path, meta_header, meta_rows, gene_ids, matrix,
		value_fmt='%.6g'):
	"""Write a wide per-cell matrix: `meta_header` index cols + one col per gene.

	`meta_rows` is a list of index tuples (one per cell), `matrix` is
	(n_cells, n_genes). Values are formatted with `value_fmt` (set to None to
	write raw ints).
	"""
	with open(path, 'w') as f:
		w = csv.writer(f, delimiter='\t')
		w.writerow(list(meta_header) + list(gene_ids))
		for meta, vals in zip(meta_rows, matrix):
			if value_fmt is None:
				formatted = [int(v) for v in vals]
			else:
				formatted = [value_fmt % v for v in vals]
			w.writerow(list(meta) + formatted)


def read_per_cell_matrix(path, n_meta):
	"""Read a wide per-cell matrix written by write_per_cell_matrix.

	Returns (meta_header list, meta_rows list-of-lists, gene_ids list,
	matrix np.ndarray of floats).
	"""
	with open(path) as f:
		r = csv.reader(f, delimiter='\t')
		header = next(r)
		meta_header = header[:n_meta]
		gene_ids = header[n_meta:]
		meta_rows = []
		data = []
		for row in r:
			if not row:
				continue
			meta_rows.append(row[:n_meta])
			data.append([float(x) for x in row[n_meta:]])
	return meta_header, meta_rows, gene_ids, np.array(data)


# Fixed basename for the raw-extraction outputs. Downstream analyses look for
# these files under this name (NOT the consumer's own plotOutFileName), so the
# extraction and its consumers agree on one place regardless of how each plot is
# invoked.
RAW_EXTRACT_BASENAME = 'subgen_raw_extract'
# Standard file-name suffixes for the raw-extraction outputs.
SYNTH_PER_CELL_SUFFIX = '_synth_per_cell.tsv'
MAX_MRNA_PER_CELL_SUFFIX = '_max_mrna_per_cell.tsv'
MAX_PROTEIN_PER_CELL_SUFFIX = '_max_protein_per_cell.tsv'
LINEAGE_SUCCESS_SUFFIX = '_lineage_success.tsv'
GENES_SUFFIX = '_genes.tsv'
# Number of leading index columns in the per-cell matrices.
PER_CELL_N_META = 3  # seed, generation, is_successful_lineage


def raw_extract_prefix(plot_out_dir):
	"""Path prefix the raw-extraction writes to (and consumers read from)."""
	return os.path.join(plot_out_dir, RAW_EXTRACT_BASENAME)


def _require_raw_file(plot_out_dir, suffix):
	path = raw_extract_prefix(plot_out_dir) + suffix
	if not os.path.isfile(path):
		raise FileNotFoundError(
			'Raw-extraction file not found: %s\nRun the extraction first:\n'
			'  python runscripts/manual/analysisCohort.py '
			'--plot subgen_raw_extract.py <sim_dir>' % path)
	return path


def load_raw_synth(plot_out_dir):
	"""Load the per-cell completed-transcript matrix from the raw extraction.

	Returns (seeds, generations, is_successful, gene_ids, matrix) where the
	first three are per-cell arrays, gene_ids is the column key, and matrix is
	(n_cells, n_genes) of completed transcripts per generation.
	"""
	path = _require_raw_file(plot_out_dir, SYNTH_PER_CELL_SUFFIX)
	_, meta_rows, gene_ids, matrix = read_per_cell_matrix(path, PER_CELL_N_META)
	seeds = np.array([int(r[0]) for r in meta_rows])
	generations = np.array([int(r[1]) for r in meta_rows])
	is_successful = np.array([int(r[2]) for r in meta_rows], dtype=bool)
	return seeds, generations, is_successful, gene_ids, matrix


def load_raw_max(plot_out_dir, which):
	"""Load a per-cell max-count matrix ('mrna' or 'protein') from the raw
	extraction. Returns (seeds, generations, is_successful, gene_ids, matrix)."""
	suffix = MAX_MRNA_PER_CELL_SUFFIX if which == 'mrna' \
		else MAX_PROTEIN_PER_CELL_SUFFIX
	path = _require_raw_file(plot_out_dir, suffix)
	_, meta_rows, gene_ids, matrix = read_per_cell_matrix(path, PER_CELL_N_META)
	seeds = np.array([int(r[0]) for r in meta_rows])
	generations = np.array([int(r[1]) for r in meta_rows])
	is_successful = np.array([int(r[2]) for r in meta_rows], dtype=bool)
	return seeds, generations, is_successful, gene_ids, matrix


def canonical_def5_classification(plot_out_dir):
	"""Canonical subgen classification (Form B) from the raw extraction.

	Restricts to strict-successful lineages, collapses to per-lineage Def-5
	rates, and classifies each gene by its 95% CI vs 1 transcript/gen. Also
	returns Form A (pooled cell-weighted mean over successful cells) and
	Definition 4 (fraction of successful cells with >= 1 completed transcript).

	Returns a dict:
	  gene_ids        column key (list)
	  stats           Form B dict (mean/std/se/ci_low/ci_high/cat/n) -- 'cat' is
	                  the canonical per-gene subgen label
	  formA           per-gene pooled mean completed transcripts (Form A)
	  p_any_synth     per-gene Definition-4 probability
	  n_lineages      number of successful lineages used
	  lineage_seeds   the successful seeds used, in row order
	  seeds/generations/is_successful/synth  the raw per-cell arrays
	"""
	seeds, generations, is_successful, gene_ids, synth = load_raw_synth(
		plot_out_dir)
	n_genes = len(gene_ids)
	successful_seeds = {int(s) for s in seeds[is_successful]}
	lambda_matrix, lineage_seeds = build_lineage_lambda(
		seeds, synth, restrict_seeds=successful_seeds)
	stats = classify_def5_ci(
		lambda_matrix, np.arange(len(lineage_seeds)), n_genes)
	formA = pooled_mean(synth, row_mask=is_successful)
	succ_synth = synth[is_successful]
	p_any_synth = (succ_synth > 0).mean(axis=0) if succ_synth.shape[0] \
		else np.full(n_genes, np.nan)
	return {
		'gene_ids': gene_ids,
		'stats': stats,
		'formA': formA,
		'p_any_synth': p_any_synth,
		'n_lineages': len(lineage_seeds),
		'lineage_seeds': lineage_seeds,
		'successful_seeds': successful_seeds,
		'seeds': seeds,
		'generations': generations,
		'is_successful': is_successful,
		'synth': synth,
		}


def load_raw_genes(plot_out_dir):
	"""Load the gene key (gene_id, cistron_id, monomer_id) as parallel lists."""
	path = _require_raw_file(plot_out_dir, GENES_SUFFIX)
	gene_ids, cistron_ids, monomer_ids = [], [], []
	with open(path) as f:
		r = csv.reader(f, delimiter='\t')
		next(r)  # header
		for row in r:
			if not row:
				continue
			gene_ids.append(row[0])
			cistron_ids.append(row[1])
			monomer_ids.append(row[2])
	return gene_ids, cistron_ids, monomer_ids
