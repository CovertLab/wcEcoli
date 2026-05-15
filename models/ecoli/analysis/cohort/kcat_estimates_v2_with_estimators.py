"""
v2-with-estimators: v2 distribution analysis plus two candidate kcat
estimators overlaid on each kcat panel.

Same per-cell read + accumulation pipeline as kcat_estimates_v2.py, with the
following additions:

  1. Top-K gap estimator.  Keeps a larger top-K of kcat samples per pair
     (cap 5000; effective K = max(100, min(5000, 0.002 * n_samples)) so
     resolution scales with cohort size).  Sorts descending, finds the
     largest log-gap between adjacent samples; if it exceeds the threshold
     (~3x ratio), returns the value just below the gap, else returns max.
     Trims spike-shaped distributions; falls back to max on well-behaved
     ones, so it is never more punishing than the existing max estimator.

  2. Smoothed-max estimator (short window).  Per cell, masks kcat <= 0
     to NaN, runs a scipy median_filter with W=5 timesteps, then takes the
     nanmax.  Aggregates across cells with nanmax.  W=5 kills 1-2-sample
     numerical glitches but preserves any spike sustained for >= 3
     consecutive samples -- much less aggressive than the old W=100.

  3. Recommended estimator: max(gap, smoothed_max).  Each of the two
     above is conservative in a different dimension (gap can't see modes
     larger than K; smoothed can't see scattered modes), so taking the max
     gives the benefit of the doubt -- if either estimator believes a
     value is real, the recommended estimate keeps it.

PDF pages are sorted by max / recommended descending so the reactions
where the estimator most aggressively trims max appear first.  All three
estimator values are written to the kcat subplot's suptitle and drawn as
distinct dashed vertical lines (purple = gap, green = smoothed_max,
black solid = recommended).

Outputs in plotOutDir (two parallel "views" of the data):
  default view -- each quantity filtered to its own positive-finite samples
    (enzyme distributions include every timestep where the enzyme exists).
  active view -- every quantity restricted to timesteps where the kcat
    estimate is finite and > 0 (no zeros, infs, or NaNs).

  Per-view CSVs (always produced):
	* kcat_estimates_v2_with_estimators_summary.csv  -- four rows per pair,
	      one per quantity, plus three estimator columns on the kcat row
	      (kcat_gap_estimate, kcat_smoothed_max_estimate, kcat_recommended).
	      Rows sorted by kcat max/recommended descending.
	* kcat_estimates_v2_with_estimators_summary_active.csv

  Per-view PDFs (controlled by PLOT_STYLES_ENABLED; rug is on by default):
	* kcat_estimates_v2_with_estimators_dist_logy{,_active}.pdf
	* kcat_estimates_v2_with_estimators_dist_ecdf{,_active}.pdf
	* kcat_estimates_v2_with_estimators_dist_logy_rug{,_active}.pdf

Toggle MAKE_DISTRIBUTION_PDFS = False to produce only the CSVs.

Memory: log-spaced histograms (20 bins/decade) per pair plus running
min/max/sum/n and two top-K accumulators (small one for the rug, larger
one for the gap detector).  Resident footprint scales with n_pairs, not
with the number of cells.
"""

import contextlib
import csv
import datetime
import os
import pickle

import numpy as np
from scipy.ndimage import median_filter
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

from models.ecoli.analysis import cohortAnalysisPlot
from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from models.ecoli.processes.metabolism import (
	COUNTS_UNITS, VOLUME_UNITS, TIME_UNITS, MASS_UNITS)
from wholecell.io.tablereader import TableReader
from wholecell.utils import units


# ============================================================================
# Configuration
# ============================================================================

IGNORE_FIRST_N_GENS = 4

# Toggle off to skip plotting and produce only the CSV (much faster).
MAKE_DISTRIBUTION_PDFS = True

# When True (and _select_kcat returns a value), emit a second set of three
# PDFs with the chosen kcat overlaid as a dashed vertical line.  Stays off
# until the v2 estimator below is implemented.
MAKE_KCAT_OVERLAY_PDFS = False

# Number of largest sample values to retain per (pair, quantity) for the rug
# plot in the "_dist_logy_rug" PDF.
TOP_K_RUG = 20

# ---------------------------------------------------------------------------
# Estimator parameters
# ---------------------------------------------------------------------------
# Top-K gap detector: keep up to TOP_K_GAP_CAP largest kcat samples per pair.
# At estimation time the effective K used is
#   k_eff = max(TOP_K_GAP_FLOOR, min(TOP_K_GAP_CAP, int(TOP_K_GAP_FRACTION * n)))
# which gives 100-sample minimum resolution for rare reactions and scales
# proportionally for cohorts with millions of samples.
TOP_K_GAP_CAP = 5000
TOP_K_GAP_FLOOR = 100
TOP_K_GAP_FRACTION = 0.002

# Log10 ratio threshold for "this is a gap, not a continuation of the tail".
# 0.5 corresponds to a ~3x jump between adjacent sorted samples; anything
# below this is treated as part of the same tail.
GAP_THRESHOLD_LOG10 = 0.5

# Minimum positive-finite kcat samples required to trust the gap detector.
# Below this we fall back to max (sparse data => no reliable gap inference).
GAP_MIN_SAMPLES = 100

# Median-filter window (timesteps) for the smoothed-max estimator.  W=5
# kills any 1-2-sample numerical spike but preserves anything sustained
# for >=3 consecutive samples.  The old kcat_estimations.py used W=100,
# which erased real metabolic transients of 10-50 s -- avoid going long.
SMOOTH_WINDOW = 5

# Which plot styles to emit (per view).  The full set of available styles is
# {'logy', 'ecdf', 'logy_rug'}; trim or extend this list to control which
# PDFs are produced.  The rug variant carries the most information for the
# v2 estimator (it shows whether the max is one outlier or a real high mode),
# so it is the only style on by default.  Each style is emitted for every
# view that is enabled, so 'logy_rug' here produces both
# kcat_estimates_v2_dist_logy_rug.pdf and ..._dist_logy_rug_active.pdf.
PLOT_STYLES_AVAILABLE = ('logy', 'ecdf', 'logy_rug')
PLOT_STYLES_ENABLED = ('logy_rug',)

_REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(
	os.path.dirname(os.path.abspath(__file__))))))
below_line_directory = os.path.join(
	_REPO_ROOT, "reconstruction", "ecoli", "scripts",
	"new_gene_below_line_proteome_ids")
below_line_essential_monomer_ids_filepath = os.path.join(
	below_line_directory, "below_line_essential_monomer_ids_variant16.csv")
below_line_essential_complex_ids_filepath = os.path.join(
	below_line_directory, "below_line_essential_complex_ids_variant16.csv")


def _logbins(lo_decade, hi_decade, per_decade=20):
	"""Log-spaced bin edges spanning [10**lo_decade, 10**hi_decade]."""
	n_edges = (hi_decade - lo_decade) * per_decade + 1
	return np.logspace(lo_decade, hi_decade, n_edges)


# Per-quantity bin edges and axis labels.  Ranges are chosen wide enough to
# cover plausible values for E. coli at this scale; values outside the range
# are clipped into the edge bins so the spike is still represented.  20 bins
# per decade leaves single-sample maxima distinguishable from the bulk.
QUANTITY_BINS = {
	'kcat':          (_logbins(-3, 8),   'kcat estimate (1/s)'),
	'flux':          (_logbins(-6, 4),   'flux (mmol/g DCW/h)'),
	'enzyme_conc':   (_logbins(-10, -1), 'enzyme concentration (M)'),
	'enzyme_counts': (_logbins(0, 7),    'enzyme counts'),
}
QUANTITY_ORDER = ['kcat', 'flux', 'enzyme_conc', 'enzyme_counts']
# 2x2 grid position (row, col) for each quantity in the per-pair figure.
SUBPLOT_POSITIONS = {
	'kcat':          (0, 0),
	'flux':          (0, 1),
	'enzyme_conc':   (1, 0),
	'enzyme_counts': (1, 1),
}

PERCENTILES = [
	('p05', 0.05), ('p25', 0.25), ('median', 0.50),
	('p75', 0.75), ('p95', 0.95), ('p99', 0.99),
]

# Two parallel views of the data:
#   * 'default'      -- each quantity is filtered to its own positive-finite
#                       samples (current behavior; conc/counts include all
#                       timesteps where the enzyme exists).
#   * 'active_only'  -- every quantity is restricted to timesteps where flux
#                       is finite and > 0, so the enzyme distributions reflect
#                       only when the reaction is actually active.
VIEW_DEFAULT = 'default'
VIEW_ACTIVE = 'active'
VIEWS = [VIEW_DEFAULT, VIEW_ACTIVE]
VIEW_SUFFIX = {VIEW_DEFAULT: '', VIEW_ACTIVE: '_active'}


# ============================================================================
# Helpers
# ============================================================================

def _read_csv_ids(filepath):
	"""Read the first column of a CSV, skipping the header row."""
	ids = []
	with open(filepath, 'r') as f:
		reader = csv.reader(f)
		next(reader)
		for row in reader:
			ids.append(row[0])
	return ids


def _percentile_from_hist(hist, edges, q):
	"""Recover a percentile from a log-binned histogram by linear-in-log interp.

	Used to derive median/p25/p75/p95/p99 without retaining raw samples.
	Returns NaN if the histogram is empty.
	"""
	total = hist.sum()
	if total == 0:
		return np.nan
	target = q * total
	cum = np.cumsum(hist)
	i = int(np.searchsorted(cum, target, side='left'))
	if i >= len(hist):
		return edges[-1]
	prev_cum = cum[i - 1] if i > 0 else 0
	bin_count = cum[i] - prev_cum
	frac = (target - prev_cum) / bin_count if bin_count > 0 else 0.0
	log_lo = np.log10(edges[i])
	log_hi = np.log10(edges[i + 1])
	return 10.0 ** (log_lo + frac * (log_hi - log_lo))


def _merge_topk(existing, new_vals, k):
	"""Merge two arrays of candidate top-K samples, keep the k largest."""
	if existing.size == 0:
		combined = np.asarray(new_vals, dtype=float)
	elif np.asarray(new_vals).size == 0:
		combined = existing
	else:
		combined = np.concatenate([existing, np.asarray(new_vals, dtype=float)])
	if combined.size <= k:
		combined.sort()
		return combined
	out = np.partition(combined, -k)[-k:]
	out.sort()
	return out


def _gap_estimate(top_k_largest, n_samples, max_val, p99_floor=None,
		threshold_log10=GAP_THRESHOLD_LOG10):
	"""Top-K gap estimator for kcat.

	`top_k_largest` is a 1-D array of up to TOP_K_GAP_CAP largest positive
	kcat samples observed across all cells (ascending or unsorted).  Sort
	descending, compute log10 ratios between adjacent samples, and find the
	largest gap.  If the largest gap exceeds `threshold_log10`, return the
	value just below the gap (i.e. trim the spike).  Otherwise return
	`max_val` unchanged so the estimator is never more punishing than max.

	Safeguards (so the detector cannot recommend values below the bulk):
	  * Requires n_samples >= GAP_MIN_SAMPLES; below that, return max_val.
	  * If a candidate trimmed value sits below `p99_floor` (the bulk's
	    high end), return max_val -- a "trim" that drops below p99 is
	    almost certainly a sparse-data artifact, not a real spike-vs-bulk
	    boundary.

	K_effective = max(TOP_K_GAP_FLOOR, min(top_k_largest.size,
	                                       int(TOP_K_GAP_FRACTION * n_samples)))
	scales detector depth with cohort size: small cohorts get the floor,
	million-sample cohorts use thousands of top values.
	"""
	if not np.isfinite(max_val):
		return np.nan
	if n_samples < GAP_MIN_SAMPLES:
		return float(max_val)
	if top_k_largest is None or top_k_largest.size < 5:
		return float(max_val)

	k_eff = max(
		TOP_K_GAP_FLOOR,
		min(int(top_k_largest.size),
			int(TOP_K_GAP_FRACTION * n_samples)),
	)
	k_eff = min(k_eff, int(top_k_largest.size))

	desc = np.sort(top_k_largest)[::-1][:k_eff]
	desc = desc[np.isfinite(desc) & (desc > 0)]
	if desc.size < 2:
		return float(max_val)

	log_desc = np.log10(desc)
	gaps = log_desc[:-1] - log_desc[1:]  # length k_eff - 1
	max_gap_idx = int(np.argmax(gaps))
	if gaps[max_gap_idx] <= threshold_log10:
		return float(max_val)
	trimmed = float(desc[max_gap_idx + 1])
	# Bulk-floor safeguard: never drop below p99.
	if (p99_floor is not None
			and np.isfinite(p99_floor)
			and trimmed < p99_floor):
		return float(max_val)
	return trimmed


def _select_kcat(pair_entry):
	"""Return the recommended kcat estimate (1/s) for this pair.

	Recommended = max(gap_estimate, smoothed_max_estimate).  Each underlying
	estimator is conservative in a different dimension (gap can't see modes
	larger than K; smoothed can't see scattered modes), so taking the max
	gives the benefit of the doubt: if either estimator believes a value is
	real, keep it.  Returns None if neither estimator produced a finite
	value.
	"""
	kc = pair_entry.get('quantities', {}).get('kcat')
	if kc is None:
		return None
	candidates = [
		kc.get('gap_estimate'),
		kc.get('smoothed_max_estimate'),
	]
	finite = [c for c in candidates
		if c is not None and np.isfinite(c) and c > 0]
	if finite:
		return max(finite)
	# Both estimators declined; fall back to max so downstream code never
	# inherits a zero / NaN "estimate".
	max_val = kc.get('max')
	return float(max_val) if (max_val is not None and np.isfinite(max_val)) else None


# ============================================================================
# Per-cell accumulation
# ============================================================================

def _update_accumulators(arr, edges, hist_acc, top_k_acc, stats_acc,
		external_mask=None):
	"""Update histogram + running stats + top-K from one cell's samples.

	`arr` is shape (T, n_pairs).  Only positive, finite samples are counted.
	If `external_mask` (same shape as arr) is provided, it is AND-ed with the
	per-quantity positivity check -- used by the "active-only" view to
	restrict every quantity to timesteps where the reaction flux was nonzero.
	Vectorized across pairs for speed.
	"""
	valid = np.isfinite(arr) & (arr > 0)
	if external_mask is not None:
		valid = valid & external_mask
	# Per-pair scalar stats
	n_per_pair = valid.sum(axis=0)
	any_valid = n_per_pair > 0

	if not np.any(any_valid):
		return

	masked_for_sum = np.where(valid, arr, 0.0)
	sum_per_pair = masked_for_sum.sum(axis=0)
	masked_for_max = np.where(valid, arr, -np.inf)
	max_per_pair = masked_for_max.max(axis=0)
	masked_for_min = np.where(valid, arr, np.inf)
	min_per_pair = masked_for_min.min(axis=0)

	stats_acc['n'] += n_per_pair
	stats_acc['sum'] += sum_per_pair
	stats_acc['n_cells'] += any_valid.astype(np.int64)
	stats_acc['min'] = np.minimum(stats_acc['min'], min_per_pair)
	stats_acc['max'] = np.maximum(stats_acc['max'], max_per_pair)

	# Histogram update via np.add.at over the flattened (sample, pair) grid.
	n_pairs = arr.shape[1]
	n_bins = len(edges) - 1
	clipped = np.clip(arr, edges[0], edges[-1])
	bin_idx = np.searchsorted(edges, clipped, side='right') - 1
	bin_idx = np.clip(bin_idx, 0, n_bins - 1)
	flat_pair = np.broadcast_to(np.arange(n_pairs), arr.shape)
	mask = valid.ravel()
	np.add.at(
		hist_acc,
		(flat_pair.ravel()[mask], bin_idx.ravel()[mask]),
		1,
	)

	# Top-K per pair from this cell.
	if arr.shape[0] >= TOP_K_RUG:
		# argpartition along time axis, then take_along_axis to pull values.
		part = np.argpartition(masked_for_max, -TOP_K_RUG, axis=0)[-TOP_K_RUG:, :]
		cell_top = np.take_along_axis(masked_for_max, part, axis=0)  # (K, n_pairs)
	else:
		cell_top = masked_for_max  # (T, n_pairs)
	for p in np.where(any_valid)[0]:
		col = cell_top[:, p]
		col = col[np.isfinite(col)]
		if col.size:
			top_k_acc[p] = _merge_topk(top_k_acc[p], col, TOP_K_RUG)


# ============================================================================
# Main analysis
# ============================================================================

class Plot(cohortAnalysisPlot.CohortAnalysisPlot):
	def do_plot(self, variantDir, plotOutDir, plotOutFileName, simDataFile,
				validationDataFile, metadata):
		with open(simDataFile, 'rb') as f:
			sim_data = pickle.load(f)

		# --- Below-line essential catalyst set (same as kcat_estimations.py) ---
		below_line_essential_monomer_ids = _read_csv_ids(
			below_line_essential_monomer_ids_filepath)
		below_line_essential_complex_ids = _read_csv_ids(
			below_line_essential_complex_ids_filepath)
		below_line_essential_ids = set(
			below_line_essential_monomer_ids + below_line_essential_complex_ids)

		ap = AnalysisPaths(variantDir, cohort_plot=True)
		cell_paths = ap.get_cells(
			generation=np.arange(IGNORE_FIRST_N_GENS, ap.n_generation),
			only_successful=True)

		if len(cell_paths) == 0:
			print('Skipping analysis -- not enough simulations run.')
			return

		# --- Reaction-catalyst pair identification ---
		fba_reader = TableReader(
			os.path.join(cell_paths[0], 'simOut', 'FBAResults'))
		listener_fba_reaction_ids = fba_reader.readAttribute('reactionIDs')
		listener_catalyst_ids = fba_reader.readAttribute('catalyst_ids')

		reaction_id_to_catalyst_ids_dict = (
			sim_data.process.metabolism.reaction_catalysts)
		catalyst_id_to_reaction_ids_dict = {}
		for reaction_id, cat_ids in reaction_id_to_catalyst_ids_dict.items():
			for cid in cat_ids:
				catalyst_id_to_reaction_ids_dict.setdefault(cid, []).append(
					reaction_id)

		below_line_essential_catalyst_ids_unique = list(
			set(below_line_essential_ids) & set(listener_catalyst_ids))

		pair_reaction_ids = []
		pair_catalyst_ids = []
		for catalyst_id in below_line_essential_catalyst_ids_unique:
			if catalyst_id not in catalyst_id_to_reaction_ids_dict:
				print(f"Warning: catalyst {catalyst_id} not in "
					f"reaction_catalysts mapping.")
				continue
			for reaction_id in catalyst_id_to_reaction_ids_dict[catalyst_id]:
				pair_reaction_ids.append(reaction_id)
				pair_catalyst_ids.append(catalyst_id)

		# Drop reactions that have additional non-essential catalysts -- the
		# flux through them cannot be cleanly attributed to a single enzyme.
		reaction_id_to_idx = {r: i for i, r in enumerate(listener_fba_reaction_ids)}
		catalyst_id_to_idx = {c: i for i, c in enumerate(listener_catalyst_ids)}

		reactions_to_filter_out = set()
		for reaction_id in set(pair_reaction_ids):
			all_cats = reaction_id_to_catalyst_ids_dict.get(reaction_id, [])
			essential_cat_count = sum(
				1 for c in all_cats if c in below_line_essential_ids)
			if len(all_cats) > 1 and essential_cat_count < len(all_cats):
				reactions_to_filter_out.add(reaction_id)
		print(f"Reactions with mixed catalyst sets (filtered out): "
			f"{len(reactions_to_filter_out)}")

		rxn_indexes = []
		cat_indexes = []
		valid_pair_reaction_ids = []
		valid_pair_catalyst_ids = []
		for rxn_id, cat_id in zip(pair_reaction_ids, pair_catalyst_ids):
			if rxn_id in reactions_to_filter_out:
				continue
			if rxn_id not in reaction_id_to_idx:
				print(f"Warning: reaction {rxn_id} not in FBA reaction IDs.")
				continue
			if cat_id not in catalyst_id_to_idx:
				print(f"Warning: catalyst {cat_id} not in catalyst IDs.")
				continue
			rxn_indexes.append(reaction_id_to_idx[rxn_id])
			cat_indexes.append(catalyst_id_to_idx[cat_id])
			valid_pair_reaction_ids.append(rxn_id)
			valid_pair_catalyst_ids.append(cat_id)

		rxn_indexes = np.array(rxn_indexes)
		cat_indexes = np.array(cat_indexes)
		n_pairs = len(rxn_indexes)
		print(f"{n_pairs} valid reaction-catalyst pairs across {len(cell_paths)} cells.")

		if n_pairs == 0:
			print('No valid reaction-catalyst pairs found.')
			return

		cell_density = sim_data.constants.cell_density
		cell_density_value = cell_density.asNumber(MASS_UNITS / VOLUME_UNITS)

		# --- Allocate accumulators (one independent set per view) ---
		accumulators = {v: {} for v in VIEWS}
		top_k = {v: {} for v in VIEWS}
		running_stats = {v: {} for v in VIEWS}
		for view in VIEWS:
			for q in QUANTITY_ORDER:
				edges = QUANTITY_BINS[q][0]
				accumulators[view][q] = np.zeros(
					(n_pairs, len(edges) - 1), dtype=np.int64)
				top_k[view][q] = [
					np.empty(0, dtype=float) for _ in range(n_pairs)]
				running_stats[view][q] = {
					'min':     np.full(n_pairs, np.inf),
					'max':     np.full(n_pairs, -np.inf),
					'sum':     np.zeros(n_pairs),
					'n':       np.zeros(n_pairs, dtype=np.int64),
					'n_cells': np.zeros(n_pairs, dtype=np.int64),
				}

		# Extra kcat-only accumulators for the estimators.  Both are global
		# across views (kcat distribution itself is the same in default and
		# active views; the estimators are time-series / value-set properties
		# of the raw nonzero kcat samples, not of a view-specific filter).
		top_k_gap = [np.empty(0, dtype=float) for _ in range(n_pairs)]
		per_pair_smoothed_max = np.full(n_pairs, -np.inf, dtype=float)

		# --- Per-cell read + accumulate ---
		cells_processed = 0
		for cell_path in cell_paths:
			sim_out = os.path.join(cell_path, 'simOut')
			try:
				counts_to_molar = TableReader(
					os.path.join(sim_out, 'EnzymeKinetics')
				).readColumn('countsToMolar', squeeze=False)[1:]
				cell_mass = TableReader(
					os.path.join(sim_out, 'Mass')
				).readColumn('cellMass', squeeze=True)[1:]
				dry_mass = TableReader(
					os.path.join(sim_out, 'Mass')
				).readColumn('dryMass', squeeze=True)[1:]
				catalyst_counts = TableReader(
					os.path.join(sim_out, 'FBAResults')
				).readColumn('catalyst_counts',
					indices=cat_indexes, squeeze=False)[1:]
				reaction_fluxes = TableReader(
					os.path.join(sim_out, 'FBAResults')
				).readColumn('reactionFluxes',
					indices=rxn_indexes, squeeze=False)[1:]
			except Exception as e:
				print(f"Ignored exception reading {sim_out}: {e!r}")
				continue

			conversion_coeffs = (
				dry_mass / cell_mass * cell_density_value
			)  # shape (T,)
			fluxes = (
				(COUNTS_UNITS / MASS_UNITS / TIME_UNITS)
				* (reaction_fluxes / conversion_coeffs[:, np.newaxis])
			).asNumber(units.mmol / units.g / units.h)
			fluxes[~np.isfinite(fluxes)] = np.nan

			if counts_to_molar.ndim > 1:
				conc = catalyst_counts * counts_to_molar
			else:
				conc = catalyst_counts * counts_to_molar[:, np.newaxis]

			with np.errstate(divide='ignore', invalid='ignore'):
				kcat = fluxes / conc
			kcat[~np.isfinite(kcat)] = np.nan

			sample_arrays = {
				'kcat':          kcat,
				'flux':          fluxes,
				'enzyme_conc':   conc,
				'enzyme_counts': catalyst_counts.astype(float),
			}
			# Joint mask for the active-only view: timesteps where the kcat
			# estimate itself is finite and > 0 (so flux > 0 AND conc > 0,
			# no infs or NaNs from divide-by-zero).  Same shape as the
			# per-quantity arrays so it AND-s in cleanly.
			active_mask = np.isfinite(kcat) & (kcat > 0)
			for view in VIEWS:
				mask = None if view == VIEW_DEFAULT else active_mask
				for q in QUANTITY_ORDER:
					_update_accumulators(
						sample_arrays[q],
						QUANTITY_BINS[q][0],
						accumulators[view][q],
						top_k[view][q],
						running_stats[view][q],
						external_mask=mask,
					)

			# --- Estimator accumulators (kcat only, global across views) ---
			kcat_valid = np.isfinite(kcat) & (kcat > 0)
			n_t = kcat.shape[0]
			# Top-K for the gap detector: per pair, merge this cell's top
			# samples into the global top-K (cap TOP_K_GAP_CAP).
			masked_kcat = np.where(kcat_valid, kcat, -np.inf)
			if n_t >= TOP_K_GAP_CAP:
				part = np.argpartition(
					masked_kcat, -TOP_K_GAP_CAP, axis=0)[-TOP_K_GAP_CAP:, :]
				cell_top = np.take_along_axis(masked_kcat, part, axis=0)
			else:
				cell_top = masked_kcat
			for p in range(n_pairs):
				col = cell_top[:, p]
				col = col[np.isfinite(col)]
				if col.size:
					top_k_gap[p] = _merge_topk(
						top_k_gap[p], col, TOP_K_GAP_CAP)

			# Smoothed max: per-pair median filter (W=5) on the kcat time
			# series with non-positive / non-finite samples zeroed out so
			# they pull the median down for isolated spikes.  Result is
			# masked back to NaN at originally-invalid positions and the
			# nanmax taken across the cell, then merged into the global
			# per-pair smoothed max via np.maximum.
			kcat_for_filter = np.where(kcat_valid, kcat, 0.0)
			for p in range(n_pairs):
				if not np.any(kcat_valid[:, p]):
					continue
				smoothed = median_filter(
					kcat_for_filter[:, p], size=SMOOTH_WINDOW,
					mode='nearest')
				smoothed[~kcat_valid[:, p]] = np.nan
				if np.any(~np.isnan(smoothed)):
					cell_max = float(np.nanmax(smoothed))
				else:
					cell_max = -np.inf
				# smoothed = 0 means "no kcat value persisted for >= 3
				# consecutive timesteps" (the filter erased every spike).
				# Treat that as "no signal" so it won't pollute the
				# recommended estimate downstream.
				if cell_max <= 0:
					cell_max = -np.inf
				if cell_max > per_pair_smoothed_max[p]:
					per_pair_smoothed_max[p] = cell_max

			cells_processed += 1

		if cells_processed == 0:
			print('No cells successfully processed.')
			return
		print(f"Processed {cells_processed} cells.")

		# --- Build per-pair summary stats (per view) ---
		def _ratio(numer, denom):
			if denom is None or not np.isfinite(denom) or denom <= 0:
				return np.nan
			return numer / denom

		pair_stats_by_view = {}
		sorted_indices_by_view = {}
		for view in VIEWS:
			view_stats = []
			for p in range(n_pairs):
				entry = {
					'reaction_id': valid_pair_reaction_ids[p],
					'catalyst_id': valid_pair_catalyst_ids[p],
					'quantities':  {},
				}
				for q in QUANTITY_ORDER:
					stats = running_stats[view][q]
					n = int(stats['n'][p])
					if n == 0:
						entry['quantities'][q] = None
						continue
					edges = QUANTITY_BINS[q][0]
					hist = accumulators[view][q][p]
					mean_val = stats['sum'][p] / n
					percentiles = {
						label: _percentile_from_hist(hist, edges, qval)
						for label, qval in PERCENTILES
					}
					max_val = float(stats['max'][p])
					min_val = float(stats['min'][p])
					p99 = percentiles['p99']
					median = percentiles['median']
					entry['quantities'][q] = {
						'n':       n,
						'n_cells': int(stats['n_cells'][p]),
						'min':     min_val,
						'max':     max_val,
						'mean':    float(mean_val),
						**{k: float(v) for k, v in percentiles.items()},
						'max_over_p99':    _ratio(max_val, p99),
						'max_over_median': _ratio(max_val, median),
						'max_over_mean':   _ratio(max_val, mean_val),
					}
					if q == 'kcat':
						p99_for_floor = (percentiles['p99']
							if np.isfinite(percentiles['p99']) else None)
						gap_est = _gap_estimate(
							top_k_gap[p], n, max_val,
							p99_floor=p99_for_floor)
						sm_raw = per_pair_smoothed_max[p]
						sm_est = (float(sm_raw)
							if np.isfinite(sm_raw) and sm_raw > 0
							else np.nan)
						entry['quantities'][q]['gap_estimate'] = gap_est
						entry['quantities'][q][
							'smoothed_max_estimate'] = sm_est
						# Recommended = max of valid positive estimators;
						# fall back to max if neither produced a usable
						# value, so sparse pairs are never punished.
						finite = [v for v in (gap_est, sm_est)
							if v is not None and np.isfinite(v) and v > 0]
						if finite:
							rec = max(finite)
						elif np.isfinite(max_val):
							rec = float(max_val)
						else:
							rec = np.nan
						entry['quantities'][q]['recommended_estimate'] = rec
						entry['quantities'][q]['max_over_recommended'] = (
							_ratio(max_val, rec))
						# Sort key: max/gap_estimate surfaces every pair
						# where the gap detector flagged a bimodality,
						# even when smoothed_max disagreed and kept the
						# spike in `recommended`.  Falls back to 1.0
						# when gap == max (no trim).
						entry['quantities'][q]['max_over_gap'] = (
							_ratio(max_val, gap_est))
				view_stats.append(entry)
			pair_stats_by_view[view] = view_stats

			def _sort_key(p, _stats=view_stats):
				kc = _stats[p]['quantities'].get('kcat')
				if kc is None:
					return -np.inf
				r = kc.get('max_over_gap')
				return r if (r is not None and np.isfinite(r)) else -np.inf

			sorted_indices_by_view[view] = sorted(
				range(n_pairs), key=_sort_key, reverse=True)

		# --- Write summary CSVs (one per view) ---
		timestamp = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
		n_gens_used = ap.n_generation - IGNORE_FIRST_N_GENS
		view_note = {
			VIEW_DEFAULT: 'per-quantity positive-finite filter',
			VIEW_ACTIVE:  'active-only: all quantities restricted to timesteps where the kcat estimate is finite and > 0',
		}

		def _fmt(x):
			return f"{x:.6g}" if (x is not None and np.isfinite(x)) else ''

		for view in VIEWS:
			provenance = (
				f'# Generated by models/ecoli/analysis/cohort/kcat_estimates_v2.py'
				f' on {timestamp}.'
				f' View: {view} ({view_note[view]}).'
				f' {ap.n_seed} seed(s), gens {IGNORE_FIRST_N_GENS}-{ap.n_generation - 1}'
				f' ({n_gens_used} gens used),'
				f' {cells_processed}/{len(cell_paths)} cells processed.\n'
			)
			csv_path = os.path.join(
				plotOutDir,
				f'kcat_estimates_v2_with_estimators_summary'
				f'{VIEW_SUFFIX[view]}.csv')
			pair_stats = pair_stats_by_view[view]
			sorted_indices = sorted_indices_by_view[view]
			with open(csv_path, 'w', newline='', encoding='utf-8') as fh:
				fh.write(provenance)
				fh.write(
					f'# kcat_gap_estimate: top-K gap detector '
					f'(K_eff up to {TOP_K_GAP_CAP}, scaled to sample count;'
					f' gap threshold {GAP_THRESHOLD_LOG10:.2f} log10;'
					f' falls back to max if n<{GAP_MIN_SAMPLES} or'
					f' trim < p99).\n'
					f'# kcat_smoothed_max_estimate: median filter W='
					f'{SMOOTH_WINDOW} timesteps on nonzero kcat,'
					f' nanmax across cells (empty when no kcat persists'
					f' for >=3 consecutive samples).\n'
					f'# kcat_recommended: max(gap, smoothed_max), falling'
					f' back to max if neither estimator produced a'
					f' valid positive value.\n'
					f'# Sort: rows sorted by kcat max/gap_estimate'
					f' descending, so pairs the gap detector trimmed'
					f' appear first regardless of what smoothed_max did.\n')
				writer = csv.writer(fh)
				writer.writerow([
					'reaction_id', 'catalyst_id', 'quantity',
					'n_samples', 'n_cells_with_data',
					'min', 'p05', 'p25', 'median', 'mean',
					'p75', 'p95', 'p99', 'max',
					'max_over_p99', 'max_over_median', 'max_over_mean',
					'kcat_gap_estimate', 'kcat_smoothed_max_estimate',
					'kcat_recommended',
					'max_over_gap', 'max_over_recommended',
				])
				for p in sorted_indices:
					entry = pair_stats[p]
					for q in QUANTITY_ORDER:
						s = entry['quantities'].get(q)
						if s is None:
							writer.writerow([
								entry['reaction_id'], entry['catalyst_id'], q,
								0, 0,
								'', '', '', '', '', '', '', '', '',
								'', '', '',
								'', '', '', '', '',
							])
							continue
						extra_cols = ['', '', '', '', '']
						if q == 'kcat':
							extra_cols = [
								_fmt(s.get('gap_estimate')),
								_fmt(s.get('smoothed_max_estimate')),
								_fmt(s.get('recommended_estimate')),
								_fmt(s.get('max_over_gap')),
								_fmt(s.get('max_over_recommended')),
							]
						writer.writerow([
							entry['reaction_id'], entry['catalyst_id'], q,
							s['n'], s['n_cells'],
							_fmt(s['min']), _fmt(s['p05']), _fmt(s['p25']),
							_fmt(s['median']), _fmt(s['mean']),
							_fmt(s['p75']), _fmt(s['p95']), _fmt(s['p99']),
							_fmt(s['max']),
							_fmt(s['max_over_p99']), _fmt(s['max_over_median']),
							_fmt(s['max_over_mean']),
							*extra_cols,
						])
			print(f"Wrote {csv_path}")

		# --- PDF generation (per view) ---
		# v2_with_estimators always draws the three estimator lines on the
		# kcat panel directly, so the legacy MAKE_KCAT_OVERLAY_PDFS path
		# from v2 is unused here.
		if MAKE_DISTRIBUTION_PDFS:
			for view in VIEWS:
				self._make_pdfs(
					plotOutDir,
					sorted_indices_by_view[view],
					pair_stats_by_view[view],
					accumulators[view], top_k[view],
					chosen_kcats=None, suffix=VIEW_SUFFIX[view])

	# ------------------------------------------------------------------------
	# Plotting
	# ------------------------------------------------------------------------

	def _make_pdfs(self, plotOutDir, sorted_indices, pair_stats,
				   accumulators, top_k, chosen_kcats=None, suffix=''):
		# Honor module-level enable list, but preserve the canonical style
		# order so output ordering is deterministic regardless of how the
		# user wrote PLOT_STYLES_ENABLED.
		styles = [s for s in PLOT_STYLES_AVAILABLE if s in PLOT_STYLES_ENABLED]
		if not styles:
			print("No plot styles enabled; skipping PDF generation.")
			return
		paths = {
			s: os.path.join(
				plotOutDir,
				f'kcat_estimates_v2_with_estimators_dist_{s}{suffix}.pdf')
			for s in styles
		}
		with contextlib.ExitStack() as stack:
			pdfs = {s: stack.enter_context(PdfPages(paths[s]))
				for s in styles}
			for p in sorted_indices:
				entry = pair_stats[p]
				kc = entry['quantities'].get('kcat')
				if kc is None:
					continue
				chosen = chosen_kcats[p] if chosen_kcats is not None else None
				for s in styles:
					fig = self._render_figure(
						entry, accumulators, top_k, p, s, chosen)
					pdfs[s].savefig(fig, bbox_inches='tight')
					plt.close(fig)
		for path in paths.values():
			print(f"Wrote {path}")

	def _render_figure(self, entry, accumulators, top_k, p, style, chosen_kcat):
		fig, axes = plt.subplots(2, 2, figsize=(13, 8))

		kc = entry['quantities']['kcat']
		def _fmt_ratio(x):
			return f"{x:.3g}" if (x is not None and np.isfinite(x)) else "n/a"
		gap_est = kc.get('gap_estimate')
		sm_est = kc.get('smoothed_max_estimate')
		rec_est = kc.get('recommended_estimate')
		suptitle = (
			f"{entry['reaction_id']}\n"
			f"catalyst: {entry['catalyst_id']}   |   "
			f"kcat max/p99 = {_fmt_ratio(kc['max_over_p99'])},  "
			f"max/median = {_fmt_ratio(kc['max_over_median'])},  "
			f"max/mean = {_fmt_ratio(kc['max_over_mean'])}\n"
			f"gap = {_fmt_ratio(gap_est)},  "
			f"smoothed_max = {_fmt_ratio(sm_est)},  "
			f"recommended = {_fmt_ratio(rec_est)}  "
			f"(max/gap = {_fmt_ratio(kc.get('max_over_gap'))},  "
			f"max/recommended = {_fmt_ratio(kc.get('max_over_recommended'))})"
		)
		fig.suptitle(suptitle, fontsize=10)

		for q in QUANTITY_ORDER:
			r, c = SUBPLOT_POSITIONS[q]
			ax = axes[r, c]
			edges, xlabel = QUANTITY_BINS[q]
			hist = accumulators[q][p]
			s = entry['quantities'].get(q)

			ax.set_xscale('log')
			ax.set_xlabel(xlabel)

			if s is None or hist.sum() == 0:
				ax.text(0.5, 0.5, 'no data', ha='center', va='center',
					transform=ax.transAxes, fontsize=10, color='grey')
				ax.set_title(q, fontsize=10)
				continue

			centers = np.sqrt(edges[:-1] * edges[1:])
			widths = edges[1:] - edges[:-1]

			if style == 'ecdf':
				ax.bar(centers, hist, width=widths, align='center',
					color='#1f77b4', alpha=0.7, edgecolor='none')
				ax.set_ylabel('count')
				ax2 = ax.twinx()
				cum = np.cumsum(hist) / hist.sum()
				ax2.step(edges[1:], cum, where='post',
					color='#d62728', lw=1.2)
				ax2.set_ylim(0, 1.02)
				ax2.set_ylabel('ECDF', color='#d62728')
				ax2.tick_params(axis='y', colors='#d62728')
			else:
				# logy and logy_rug share the histogram body.
				ax.bar(centers, hist, width=widths, align='center',
					color='#1f77b4', alpha=0.7, edgecolor='none')
				ax.set_yscale('log')
				ax.set_ylabel('count (log)')
				if style == 'logy_rug':
					tk = top_k[q][p]
					if tk.size:
						ymin, ymax = ax.get_ylim()
						y_rug = ymax / 1.4
						ax.scatter(tk, np.full(tk.size, y_rug),
							marker='|', s=120, c='k', linewidths=1.0)
						ax.text(0.99, 0.97, f'top {tk.size} values',
							transform=ax.transAxes, ha='right', va='top',
							fontsize=7, color='k')

			ax.set_title(
				f"{q} -- n={s['n']:,}, max/p99={_fmt_ratio(s['max_over_p99'])}",
				fontsize=9)

			# Reference / estimator lines sit at zorder 0.5 so the
			# histogram bars (default zorder 1) are drawn on top -- the
			# distribution stays readable even where the recommended line
			# crosses a populated bin.  The bars use alpha=0.7 so the line
			# colors still show through in the populated regions.
			LINE_Z = 0.5
			if np.isfinite(s['median']):
				ax.axvline(s['median'], color='grey',
					ls='--', lw=0.8, label='median', zorder=LINE_Z)
			if np.isfinite(s['p99']):
				ax.axvline(s['p99'], color='orange',
					ls='--', lw=0.8, label='p99', zorder=LINE_Z)
			if np.isfinite(s['max']):
				ax.axvline(s['max'], color='red',
					ls='--', lw=0.8, label='max', zorder=LINE_Z)

			if q == 'kcat':
				# Estimator overlays: purple = top-K gap, green = smoothed
				# max (W=5), black solid = recommended (max of the two).
				if gap_est is not None and np.isfinite(gap_est):
					ax.axvline(gap_est, color='#7e1ec6', ls='--', lw=1.0,
						label='gap est.', zorder=LINE_Z)
				if sm_est is not None and np.isfinite(sm_est):
					ax.axvline(sm_est, color='#2ca02c', ls='--', lw=1.0,
						label='smoothed_max', zorder=LINE_Z)
				if rec_est is not None and np.isfinite(rec_est):
					ax.axvline(rec_est, color='black', ls='-', lw=1.6,
						label='recommended', zorder=LINE_Z)

			if q == 'kcat':
				# On the rug variant, the rug + caption sit in the upper-right
				# corner; push the legend to upper-left to avoid overlap.
				legend_loc = 'upper left' if style == 'logy_rug' else 'upper right'
				ax.legend(loc=legend_loc, fontsize=6, framealpha=0.85,
					ncol=2)

		# Slightly more top margin than v2 (extra suptitle line for estimators).
		fig.tight_layout(rect=[0, 0, 1, 0.90])
		return fig


if __name__ == '__main__':
	Plot().cli()
