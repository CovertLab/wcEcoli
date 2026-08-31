"""
Shared helpers for the subgenerational-expression analyses.

This module is NOT an analysis Plot (it has no `Plot` class and is not listed in
`__init__.py`), so the analysis runner never imports or runs it directly. It
consolidates logic providing:

  * the gene set (protein-coding mRNA cistrons, in cistron order),
  * the STRICT "successful seed" filter (completed every generation AND no
    cell hit the 180-min doubling cap),
  * THE Definition-5 classifier: classify_def5_ci
  * readers for a per-gene def-5 table (load_def5_categories),
  * curated monomer panels, run-metadata provenance, subsampling, and small TSV
    IO helpers for the raw-extraction pipeline.

Definition 5 is the mean number of *completed* mRNA transcripts per generation,
per gene, read from TranscriptElongationListener/countRnaCistronSynthesized. The more strict definiton 
includes using CI if the gene falls in the 95% CI interval of its per-seed rate lies entirely below 1 transcript per generation 


The pooled cell-weighted mean (`pooled_mean`) is
reported alongside as a descriptive statistic, but a point estimate compared
against 1 is a diagnostic and must never be used to label genes: on sim set 1 the
per-seed point estimate called 1652-1687 genes subgen where the CI rule calls
1859, and inflated `never_expressed` from 135 to ~330.
"""

import csv
import json
import os
import subprocess
from datetime import datetime

import numpy as np


def _table_reader():
	from wholecell.io.tablereader import TableReader
	return TableReader


def _constants():
	from wholecell.utils import constants
	return constants


# Shared constants 
IGNORE_FIRST_N_GENS = 8
N_SEEDS = 128
SEED_RANGE = np.arange(0, N_SEEDS)
MAX_DOUBLING_MIN = 180.0 # mins
DOUBLING_AT_MAX_TOL = 0.01  
# treat >= 179.99 min as "hit the cap"
CI_Z = 1.96
CONFIDENCE = 0.95
SYNTH_TABLE = 'TranscriptElongationListener'
SYNTH_COLUMN = 'countRnaCistronSynthesized'
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
# Neutral figure colors 
INK = '#1b2530'
MUTED = '#5b6672'
GRID = '#e4e8ec'
SURF = '#fcfcfb'

# pseudo-single-cell timepoints the subsample scripts aim to draw across the whole cohort. 
TIMEPOINTS_TO_SAMPLE = 10000

def sample_per_seed(n_seeds, total=None):
	"""Per-seed timepoint draw for the subsample scripts.

	Divides the cohort-wide budget by the number of seeds ACTUALLY sampled, not
	by len(SEED_RANGE). 
	"""
	total = TIMEPOINTS_TO_SAMPLE if total is None else total
	if n_seeds <= 0:
		return 0
	return max(1, int(total // n_seeds))


# Hand-picked subgen genes that strains were made out of The genes are used by subgen_peak_counts.py and subgen_protein_distribution.py.

PANEL_CURATED10 = [
	'GLYCDEH-MONOMER[c]',        # gldA
	'BETAGALACTOSID-MONOMER[c]', # lacZ
	'RIBULOKIN-MONOMER[c]',      # araB
	'BAES-MONOMER[i]',           # baeS
	'G6504-MONOMER[o]',          # gfcE
	'EG11250-MONOMER[c]',        # chpS
	'EG11222-MONOMER[c]',        # alkA
	'G7263-MONOMER[c]',          # murQ
	'EG11249-MONOMER[c]',        # mazF (constitutive control)
	'EG10466-MONOMER[c]',        # hupA (constitutive control)
	]
PANEL_CURATED10_NAMES = {
	'GLYCDEH-MONOMER[c]': 'gldA',
	'BETAGALACTOSID-MONOMER[c]': 'lacZ',
	'RIBULOKIN-MONOMER[c]': 'araB',
	'BAES-MONOMER[i]': 'baeS',
	'G6504-MONOMER[o]': 'gfcE',
	'EG11250-MONOMER[c]': 'chpS',
	'EG11222-MONOMER[c]': 'alkA',
	'G7263-MONOMER[c]': 'murQ',
	'EG11249-MONOMER[c]': 'mazF',
	'EG10466-MONOMER[c]': 'hupA',
	}

# The anaerobic-respiration panel used by subgen_monomer_dynamics.py. One
# representative subunit per complex; the sibling subunits are in the name map so
# a panel can be widened without re-deriving the symbols.
PANEL_ANAEROBIC = [
	'CYTOCHROMEC-MONOMER[p]',       # nrfB
	'DMSA-MONOMER[i]',              # dmsA
	'EG11800-MONOMER[i]',           # hybB
	'EG11815-MONOMER[p]',           # torC
	'EG12244-MONOMER[i]',           # ccp
	'FDNG-MONOMER[m]',              # fdnG
	'FDOG-MONOMER[c]',              # fdoG
	'FORMATEDEHYDROGH-MONOMER[c]',  # fdhF
	'FUM-FE-S[c]',                  # frdB
	'G6848-MONOMER[i]',             # ynfH
	'G7022-MONOMER[p]',             # torZ
	'HYAA-MONOMER[i]',              # hyaA
	'HYCBSMALL-MONOMER[c]',         # hycB
	'MONOMER0-141[i]',              # hyfD
	'NARG-MONOMER[m]',              # narG
	'NARV-MONOMER[m]',              # narV
	'NRFC-MONOMER[c]',              # nrfC
	'TORA-MONOMER[p]',              # torA
	]
PANEL_ANAEROBIC_NAMES = {
	'CYTOCHROMEC-MONOMER[p]': 'nrfB',
	'CYTOCHROMEC552-MONOMER[p]': 'nrfA',
	'DMSA-MONOMER[i]': 'dmsA',
	'DMSB-MONOMER[i]': 'dmsB',
	'DMSC-MONOMER[i]': 'dmsC',
	'EG11800-MONOMER[i]': 'hybB',
	'EG11815-MONOMER[p]': 'torC',
	'EG12244-MONOMER[i]': 'ccp',
	'FDNG-MONOMER[m]': 'fdnG',
	'FDNH-MONOMER[m]': 'fdnH',
	'FDNI-MONOMER[m]': 'fdnI',
	'FDOG-MONOMER[c]': 'fdoG',
	'FDOH-MONOMER[m]': 'fdoH',
	'FDOI-MONOMER[i]': 'fdoI',
	'FORMATEDEHYDROGH-MONOMER[c]': 'fdhF',
	'FUM-FE-S[c]': 'frdB',
	'FUM-FLAVO[c]': 'frdA',
	'FUM-MEMB1[m]': 'frdC',
	'FUM-MEMB2[m]': 'frdD',
	'G6848-MONOMER[i]': 'ynfH',
	'G7022-MONOMER[p]': 'torZ',
	'G7023-MONOMER[m]': 'torY',
	'HYAA-MONOMER[i]': 'hyaA',
	'HYAB-MONOMER[i]': 'hyaB',
	'HYAC-MONOMER[i]': 'hyaC',
	'HYCBSMALL-MONOMER[c]': 'hycB',
	'HYCC-MONOMER[i]': 'hycC',
	'HYCD-MONOMER[i]': 'hycD',
	'HYCELARGE-MONOMER[c]': 'hycE',
	'HYCF-MONOMER[m]': 'hycF',
	'HYCG-MONOMER[i]': 'hycG',
	'MONOMER0-141[i]': 'hyfD',
	'MONOMER0-142[i]': 'hyfE',
	'MONOMER0-143[i]': 'hyfF',
	'MONOMER0-150[c]': 'hyfG',
	'MONOMER0-153[i]': 'hyfB',
	'MONOMER0-154[m]': 'hyfC',
	'NARG-MONOMER[m]': 'narG',
	'NARH-MONOMER[m]': 'narH',
	'NARI-MONOMER[m]': 'narI',
	'NARV-MONOMER[m]': 'narV',
	'NARY-MONOMER[m]': 'narY',
	'NARZ-MONOMER[m]': 'narZ',
	'NRFC-MONOMER[c]': 'nrfC',
	'NRFD-MONOMER[i]': 'nrfD',
	'TORA-MONOMER[p]': 'torA',
	}

CURATED_PANELS = {
	'curated10': (PANEL_CURATED10, PANEL_CURATED10_NAMES),
	'anaerobic': (PANEL_ANAEROBIC, PANEL_ANAEROBIC_NAMES),
	}


def curated_panel(name):
	"""Return (monomer_ids, {monomer_id: gene_symbol}) for a named panel."""
	if name not in CURATED_PANELS:
		raise KeyError('Unknown curated panel %r; available: %s'
			% (name, ', '.join(sorted(CURATED_PANELS))))
	return CURATED_PANELS[name]


# path helpers 
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
	constants = _constants()
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


def sim_metadata_block(sim_metadata_path, sim_metadata):
	"""The `simulation` metadata block shared by every run_metadata.json."""
	return {
		'metadata_source': sim_metadata_path,
		'git_hash': sim_metadata.get('git_hash'),
		'git_branch': sim_metadata.get('git_branch'),
		'run_time': sim_metadata.get('time'),
		'description': sim_metadata.get('description'),
		'variant': sim_metadata.get('variant'),
		'total_gens': sim_metadata.get('total_gens'),
		'total_init_sims': sim_metadata.get('total_init_sims'),
		}


def write_run_metadata(path, script, parameters, sim_metadata_path=None,
		sim_metadata=None, run_time=None, extra=None, quiet=False):
	"""Write a subgen run_metadata.json with the standard metadata layout.

	`script` is the analysis module's filename, `parameters` its run parameters,
	and `extra` any additional top-level blocks (e.g. `seeds`, `cells`,
	`category_counts`) merged in as-is. The git block is read from this module's
	repository, and `run_time` defaults to now.
	"""
	repo_dir = os.path.dirname(os.path.abspath(__file__))
	meta = {
		'analysis': {
			'script': script,
			'run_time': run_time or datetime.now().isoformat(timespec='seconds'),
			'git': git_info(repo_dir),
			'parameters': parameters,
			},
		'simulation': sim_metadata_block(
			sim_metadata_path, sim_metadata or {}),
		}
	meta.update(extra or {})
	with open(path, 'w') as f:
		json.dump(meta, f, indent=2)
	if not quiet:
		print('Wrote %s' % path)
	return meta


#  Gene set 
def get_mrna_gene_set(sim_data):
	"""protein-coding mRNA cistrons, in cistron order.

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
	TableReader = _table_reader()
	try:
		reader = TableReader(os.path.join(cell_path, 'simOut', table))
		full_ids = reader.readAttribute('cistron_ids')
		id_to_index = {cid: i for i, cid in enumerate(full_ids)}
		return True, np.array([id_to_index[c] for c in ids])
	except Exception:
		return False, None


#  Strict "successful seed" filter
def compute_seed_success(ap, n_generation, total_init_sims=None):


	"""A seed is "successful" if it completed EVERY generation AND no
	cell hit the 180-minute doubling cap. 

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

	TableReader = _table_reader()
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


# Definition-5 classifiers 
def _classify_ci(mean, ci_low, ci_high):
	"""Assign each gene to a Def-5 category from its mean and CI (def5_CI)."""
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
	"""def5_CI: per-gene mean/std/se/CI/category over selected seed rows.

	`lambda_matrix` is (n_seeds, n_genes) of per-seed Def-5 rates.
	`row_indices` selects which seed rows to include (e.g. successful ones).
	The 95% CI is the normal-approx SE of the mean; a gene is `subgen` if its
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


def build_seed_lambda(seeds_arr, synth_matrix, restrict_seeds=None):
	"""Collapse a per-cell synth matrix to per-seed Def-5 rates (def5_CI input).

	`seeds_arr` is (n_cells,) of seed ints; `synth_matrix` is (n_cells, n_genes)
	of completed-transcript counts per cell (already burned-in). For each seed,
	the seed rate is the mean over that seed's cells. Returns
	(lambda_matrix (n_seeds, n_genes), seeds_used list).
	"""
	seeds_arr = np.asarray(seeds_arr)
	unique_seeds = sorted(set(int(s) for s in seeds_arr))
	rows = []
	seeds_used = []
	for s in unique_seeds:
		if restrict_seeds is not None and s not in restrict_seeds:
			continue
		mask = seeds_arr == s
		if not mask.any():
			continue
		rows.append(synth_matrix[mask].mean(axis=0))
		seeds_used.append(s)
	lambda_matrix = np.array(rows) if rows \
		else np.zeros((0, synth_matrix.shape[1]))
	return lambda_matrix, seeds_used


def pooled_mean(synth_matrix, row_mask=None):
	"""def5 pooled value: cell-weighted mean completed transcripts per gene."""
	if row_mask is not None:
		synth_matrix = synth_matrix[row_mask]
	if synth_matrix.shape[0] == 0:
		return np.full(synth_matrix.shape[1], np.nan)
	return synth_matrix.mean(axis=0)


# Timepoint-subsampling helpers
def subsample_seed_timepoints(cell_paths_per_seed, sample_per_seed):
	"""randomly select timepoints in a selected seed to pool randomly selected data.
	currently using np.random.seed(0) so it is reproducible 


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
	# seed, and gen_starts is increasing).
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


# Raw-extraction TSV IO 

def write_per_cell_matrix(path, meta_header, meta_rows, gene_ids, matrix,
		value_fmt='%.6g'):
	"""Write a wide per-cell matrix: index cols + one col per gene.

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


# harcdcoded output names so that they can be found later by downstream analysis scripts
EXTRACT_BASENAME = 'subgen_extract'
SYNTH_PER_CELL_SUFFIX = '_synth_per_cell.tsv'
MAX_MRNA_PER_CELL_SUFFIX = '_max_mrna_per_cell.tsv'
MAX_PROTEIN_PER_CELL_SUFFIX = '_max_protein_per_cell.tsv'
FRAC_PROTEIN_ZERO_PER_CELL_SUFFIX = '_frac_protein_zero_per_cell.tsv'
SEED_SUCCESS_SUFFIX = '_seed_success.tsv'
DOUBLING_TIMES_SUFFIX = '_doubling_times.tsv'
GENES_SUFFIX = '_genes.tsv'
# Number of index columns in the per-cell matrices.
PER_CELL_N_META = 3  # seed, generation, is_successful_seed


def extract_prefix(plot_out_dir):
	"""Path prefix the raw-extraction writes to (and consumers read from)."""
	return os.path.join(plot_out_dir, EXTRACT_BASENAME)


def _require_raw_file(plot_out_dir, suffix):
	path = extract_prefix(plot_out_dir) + suffix
	if not os.path.isfile(path):
		raise FileNotFoundError(
			'Raw-extraction file not found: %s\nRun the extraction first:\n'
			'  python runscripts/manual/analysisCohort.py '
			'--plot subgen_extract.py <sim_dir>' % path)
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
	"""subgen classification (def5_CI) from raw extraction.

	Restricts to strict-successful seeds, collapses to per-seed Def-5
	rates, and classifies each gene by its 95% CI vs 1 transcript/gen. Also
	returns def5 (pooled cell-weighted mean over successful cells) and
	Definition 4 (fraction of successful cells with >= 1 completed transcript).

	Returns a dict:
	  gene_ids        column key (list)
	  stats           def5_CI dict (mean/std/se/ci_low/ci_high/cat/n) -- 'cat' is
	                  the canonical per-gene subgen label
	  def5            per-gene pooled mean completed transcripts (def5)
	  p_any_synth     per-gene Definition-4 probability
	  n_seeds      number of successful seeds used
	  seeds_used   the successful seeds used, in row order
	  seeds/generations/is_successful/synth  the raw per-cell arrays
	"""
	seeds, generations, is_successful, gene_ids, synth = load_raw_synth(
		plot_out_dir)
	n_genes = len(gene_ids)
	successful_seeds = {int(s) for s in seeds[is_successful]}
	lambda_matrix, seeds_used = build_seed_lambda(
		seeds, synth, restrict_seeds=successful_seeds)
	stats = classify_def5_ci(
		lambda_matrix, np.arange(len(seeds_used)), n_genes)
	def5 = pooled_mean(synth, row_mask=is_successful)
	succ_synth = synth[is_successful]
	p_any_synth = (succ_synth > 0).mean(axis=0) if succ_synth.shape[0] \
		else np.full(n_genes, np.nan)
	return {
		'gene_ids': gene_ids,
		'stats': stats,
		'def5': def5,
		'p_any_synth': p_any_synth,
		'n_seeds': len(seeds_used),
		'seeds_used': seeds_used,
		'successful_seeds': successful_seeds,
		'seeds': seeds,
		'generations': generations,
		'is_successful': is_successful,
		'synth': synth,
		}


# Reading a per-gene def-5 table back in 
# The two tables carry the SAME def5_CI classification under different
# column names, so every reader must accept either spelling:
#   subgen_seed_ci_pergene_*.tsv -> gene_id / cistron_id / category
#   subgen_expression_table.tsv   -> gene_name / cistron_name /
#                                                  def5_CI_category
# (subgen_expression_table.py now also emits the first spelling as
# aliases, but older outputs do not, and the seed_ci tables never will.)
DEF5_GENE_COLUMNS = ('gene_id', 'gene_name')
DEF5_CISTRON_COLUMNS = ('cistron_id', 'cistron_name')
DEF5_CATEGORY_COLUMNS = ('category', 'def5_CI_category')


def _first_column(row, names, path, what):
	for name in names:
		if name in row:
			return name
	raise KeyError('%s has no %s column (looked for %s). Is it a def-5 per-gene '
		'table?' % (path, what, ' / '.join(names)))


def load_def5_table(path):
	"""Load a per-gene def-5 table as (rows, column_names).

	`rows` is a list of dicts straight from the file; `column_names` is a dict
	with the resolved `gene`, `cistron` and `category` keys for this file's
	spelling, so callers do not have to care which of the two tables
	they were handed.
	"""
	with open(path) as f:
		rows = list(csv.DictReader(f, delimiter='\t'))
	if not rows:
		raise ValueError('%s has no data rows.' % path)
	cols = {
		'gene': _first_column(rows[0], DEF5_GENE_COLUMNS, path, 'gene-id'),
		'cistron': _first_column(
			rows[0], DEF5_CISTRON_COLUMNS, path, 'cistron-id'),
		'category': _first_column(
			rows[0], DEF5_CATEGORY_COLUMNS, path, 'def-5 category'),
		}
	return rows, cols


def load_def5_categories(path, category='subgen'):
	"""Genes in one def-5 category, from either per-gene table.

	`category` defaults to 'subgen', i.e. the CI-based subgenerational set: a gene
	whose 95% CI on the per-seed completed-transcript rate lies entirely below
	1 transcript/generation. Definition 5 always means the CI form -- a point
	estimate compared against 1 is a diagnostic, never a classification.

	Returns a dict:
	  gene_ids     ordered list of matching gene ids
	  cistron_ids  ordered list of matching cistron ids (parallel to gene_ids)
	  by_gene      {gene_id: full row dict} for the matching rows
	  n_total      number of rows in the file (all categories)
	  columns      the resolved column names (see load_def5_table)
	"""
	rows, cols = load_def5_table(path)
	gene_ids, cistron_ids, by_gene = [], [], {}
	for row in rows:
		if row[cols['category']].strip() != category:
			continue
		gene_id = row[cols['gene']].strip()
		gene_ids.append(gene_id)
		cistron_ids.append(row[cols['cistron']].strip())
		by_gene[gene_id] = row
	return {
		'gene_ids': gene_ids,
		'cistron_ids': cistron_ids,
		'by_gene': by_gene,
		'n_total': len(rows),
		'columns': cols,
		}


def load_seed_success_rows(plot_out_dir):
	"""Load the raw extraction's seed table as {seed: {column: value}}.

	Column keys are the file's own header names (seed, n_gens_ran,
	reached_final_gen, completed_all_gens, n_cells_at_180, gens_at_180,
	max_doubling_min, is_successful); values are the raw strings apart from the
	int seed key. Returns {} when the file is absent.
	"""
	path = extract_prefix(plot_out_dir) + SEED_SUCCESS_SUFFIX
	if not os.path.isfile(path):
		return {}
	out = {}
	with open(path) as f:
		for row in csv.DictReader(f, delimiter='\t'):
			out[int(row['seed'])] = row
	return out


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
