"""
Demonstration of the smoothed_max (v2_smoothed_max) kcat estimator.

The smoothed_max estimator reduces a noisy per-timestep kcat = flux/[E] series
to a single bound.  Aggregated across cells it is opaque, so this analysis makes
every step visible for a few representative reaction *types* (auto-selected):

  Step 1  kcat = flux / [E]            (per timestep, per cell)
  Step 2  drop kcat <= 0 / non-finite  (validity mask)
  Step 3  zero-fill invalid samples     (so isolated spikes get pulled down)
  Step 4  median_filter(W=SMOOTH_WINDOW, mode='nearest')
  Step 5  re-mask invalid positions back to NaN
  Step 6  per-cell nanmax of the smoothed trace
  Step 7  nanmax across all cells (gens >= IGNORE_FIRST_N_GENS) -> smoothed_max

The per-cell smoothing and aggregation are imported / reproduced exactly from
kcat_estimates_v2_with_estimators.py, so the smoothed_max value printed and drawn
here equals the production estimate.

Reaction selection (auto): a first lightweight pass computes max + smoothed_max
+ valid-sample count per (reaction, catalyst) pair, then picks
  * spike-dominated  -- smoothed_max / max smallest (smoothing trimmed the most)
  * moderate-trim    -- smoothed_max / max between the spike and clean cutoffs
  * clean/sustained  -- smoothed_max / max ~ 1 with many valid samples
  * sparse           -- fewest valid samples
  * filtered-out     -- positive raw max but smoothed_max collapses to 0 (no
                        kcat persisted >= ~3 consecutive timesteps), so the
                        pair contributes no bound
  * intermittent     -- flux drops out for stretches (active a fraction of the
                        time) yet smoothed_max ~ max: the zero stretches do not
                        drag the (max-based) estimate down
Set REACTION_CATALYST_PAIRS below to a non-empty list to plot specific pairs
instead (skips the selection pass).

Outputs in plotOutDir:
  * kcat_smoothed_max_demo.pdf   -- one reaction per page (3 annotated panels)
  * kcat_smoothed_max_demo.html  -- self-contained report: a "how to read this"
        guide + the 7-step pipeline, then one section per reaction (embedded
        figure + step captions + table)

Run:
  python runscripts/manual/analysisCohort.py <run_dir> -p kcat_smoothed_max_demo
"""

import base64
import html
import io
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
# Import the production constants + below-line id readers so this demo can never
# drift from the estimator it documents.
from models.ecoli.analysis.cohort.kcat_estimates_v2_with_estimators import (
	SMOOTH_WINDOW,
	IGNORE_FIRST_N_GENS,
	_read_csv_ids,
	below_line_essential_monomer_ids_filepath,
	below_line_essential_complex_ids_filepath,
)
from wholecell.io.tablereader import TableReader
from wholecell.utils import units


# ============================================================================
# Configuration
# ============================================================================

# Non-empty -> plot exactly these (reaction_id, catalyst_id) pairs and skip the
# auto-selection pass.
REACTION_CATALYST_PAIRS = []

# Auto-selection: how many of each representative type to show.
N_SPIKE = 3
N_MODERATE = 1
N_CLEAN = 2
N_SPARSE = 1
# "filtered-out": positive raw max but smoothed_max collapses to 0 (no kcat
# persisted >= ~3 consecutive timesteps), so the pair contributes no bound --
# the spike-rejection working at its extreme.
N_FILTERED = 1
# "intermittent": flux drops out for stretches (enzyme/intermediate gone) but
# the sustained on-periods still survive, so smoothed_max ~ max -- shows the
# zero stretches do NOT drag the estimate down (it's a max, not an average).
N_INTERMITTENT = 1
# A pair is "intermittent" if its active fraction falls in this band: off
# enough to clearly show dropouts (<= MAX) but active enough to have real
# sustained on-periods rather than scattered noise (>= MIN), while still being
# clean (ratio ~ 1).
INTERMITTENT_FRAC_MIN = 0.10
INTERMITTENT_FRAC_MAX = 0.85

# smoothed_max / max <= this -> "spike-dominated" (smoothing trimmed a lot).
SPIKE_RATIO_MAX = 0.5
# smoothed_max / max >= this (with enough samples) -> "clean / sustained".
# Pairs between SPIKE_RATIO_MAX and CLEAN_RATIO_MIN are "moderate-trim".
CLEAN_RATIO_MIN = 0.9
# Minimum valid samples for a pair to qualify as a well-sampled "clean" case.
MIN_VALID_FOR_CLEAN = 500
# Prefer spike/clean demo reactions whose raw max is above this floor (matches
# MIN_KCAT_ESTIMATE in metabolism.py -- below it the kcat is treated as noise
# and never used as a bound).  Falls back to any magnitude if none qualify, so
# the demo still works on tiny test cohorts.
MIN_MAX_FOR_DEMO = 1.0

# Seeds drawn in the overview panels (the estimate itself always uses ALL
# cells, matching production).
MAX_SEEDS_TO_PLOT = 4

# Timesteps on each side of the raw-max sample shown in the single-cell zoom.
ZOOM_HALF_WIDTH = 40


# ============================================================================
# Read + smoothing helpers (mirrors kcat_estimates_v2_with_estimators.py)
# ============================================================================

def _read_cell_kcat(sim_out, rxn_indexes, cat_indexes, cell_density_value):
	"""Return (time, kcat) for one cell, kcat shape (T-1, n_sel).

	Reproduces the production read: flux/[E] with the dry/cell-mass + density
	unit conversion, first timestep dropped.
	"""
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
	).readColumn('catalyst_counts', indices=cat_indexes, squeeze=False)[1:]
	reaction_fluxes = TableReader(
		os.path.join(sim_out, 'FBAResults')
	).readColumn('reactionFluxes', indices=rxn_indexes, squeeze=False)[1:]
	t = TableReader(
		os.path.join(sim_out, 'Main')
	).readColumn('time', squeeze=True)[1:]

	conversion_coeffs = dry_mass / cell_mass * cell_density_value
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
	return t, kcat


def _smooth_cell_column(kcat_col, valid_col):
	"""Return the W=SMOOTH_WINDOW median-filtered, re-masked series (NaN where
	invalid).  Exactly the production transform (Steps 3-5)."""
	filled = np.where(valid_col, kcat_col, 0.0)
	smoothed = median_filter(filled, size=SMOOTH_WINDOW, mode='nearest')
	smoothed = smoothed.astype(float)
	smoothed[~valid_col] = np.nan
	return smoothed


def _fmt_smax(smax_val):
	"""Display string for a smoothed_max value (or its filtered-out state)."""
	if np.isfinite(smax_val) and smax_val > 0:
		return f'{smax_val:.6g}'
	return 'filtered out (0)'


def _cell_smoothed_max(kcat_col, valid_col):
	"""Per-cell smoothed nanmax (Step 6), with the production <=0 -> -inf rule."""
	if not np.any(valid_col):
		return -np.inf
	smoothed = _smooth_cell_column(kcat_col, valid_col)
	if np.any(~np.isnan(smoothed)):
		cell_max = float(np.nanmax(smoothed))
	else:
		cell_max = -np.inf
	if cell_max <= 0:
		cell_max = -np.inf
	return cell_max


# ============================================================================
# Main analysis
# ============================================================================

class Plot(cohortAnalysisPlot.CohortAnalysisPlot):
	def do_plot(self, variantDir, plotOutDir, plotOutFileName, simDataFile,
				validationDataFile, metadata):
		with open(simDataFile, 'rb') as f:
			sim_data = pickle.load(f)

		# --- Below-line essential catalyst set + (reaction, catalyst) pairs ---
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
			print('Skipping kcat_smoothed_max_demo -- not enough sims run.')
			return

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
				continue
			for reaction_id in catalyst_id_to_reaction_ids_dict[catalyst_id]:
				pair_reaction_ids.append(reaction_id)
				pair_catalyst_ids.append(catalyst_id)

		reaction_id_to_idx = {
			r: i for i, r in enumerate(listener_fba_reaction_ids)}
		catalyst_id_to_idx = {
			c: i for i, c in enumerate(listener_catalyst_ids)}

		reactions_to_filter_out = set()
		for reaction_id in set(pair_reaction_ids):
			all_cats = reaction_id_to_catalyst_ids_dict.get(reaction_id, [])
			essential_cat_count = sum(
				1 for c in all_cats if c in below_line_essential_ids)
			if len(all_cats) > 1 and essential_cat_count < len(all_cats):
				reactions_to_filter_out.add(reaction_id)

		rxn_indexes = []
		cat_indexes = []
		valid_pair_reaction_ids = []
		valid_pair_catalyst_ids = []
		for rxn_id, cat_id in zip(pair_reaction_ids, pair_catalyst_ids):
			if rxn_id in reactions_to_filter_out:
				continue
			if rxn_id not in reaction_id_to_idx:
				continue
			if cat_id not in catalyst_id_to_idx:
				continue
			rxn_indexes.append(reaction_id_to_idx[rxn_id])
			cat_indexes.append(catalyst_id_to_idx[cat_id])
			valid_pair_reaction_ids.append(rxn_id)
			valid_pair_catalyst_ids.append(cat_id)

		rxn_indexes = np.array(rxn_indexes)
		cat_indexes = np.array(cat_indexes)
		n_pairs = len(rxn_indexes)
		if n_pairs == 0:
			print('No valid reaction-catalyst pairs found.')
			return

		cell_density_value = sim_data.constants.cell_density.asNumber(
			MASS_UNITS / VOLUME_UNITS)

		pair_index_by_id = {
			(valid_pair_reaction_ids[p], valid_pair_catalyst_ids[p]): p
			for p in range(n_pairs)}

		# --- Pass 1: per-pair max / smoothed_max / n_valid (+ the cell that
		# holds the global raw max, for the zoom panel) ---
		pair_max = np.full(n_pairs, -np.inf)
		pair_smax = np.full(n_pairs, -np.inf)
		pair_nvalid = np.zeros(n_pairs, dtype=np.int64)
		pair_ncells = np.zeros(n_pairs, dtype=np.int64)
		pair_max_cell = [None] * n_pairs
		# Total timesteps seen (same for every pair -- all share the cells), so
		# valid_fraction = pair_nvalid / total_timesteps measures how often each
		# reaction actually carries flux (low => frequent dropouts).
		total_timesteps = 0

		for cell_path in cell_paths:
			sim_out = os.path.join(cell_path, 'simOut')
			try:
				_, kcat = _read_cell_kcat(
					sim_out, rxn_indexes, cat_indexes, cell_density_value)
			except Exception as e:
				print(f'Ignored exception reading {sim_out}: {e!r}')
				continue
			total_timesteps += kcat.shape[0]
			valid = np.isfinite(kcat) & (kcat > 0)
			for p in range(n_pairs):
				vc = valid[:, p]
				if not np.any(vc):
					continue
				col = kcat[:, p]
				cell_raw_max = float(np.max(col[vc]))
				if cell_raw_max > pair_max[p]:
					pair_max[p] = cell_raw_max
					pair_max_cell[p] = cell_path
				sm = _cell_smoothed_max(col, vc)
				if sm > pair_smax[p]:
					pair_smax[p] = sm
				pair_nvalid[p] += int(np.count_nonzero(vc))
				pair_ncells[p] += 1

		valid_frac = np.array([
			pair_nvalid[p] / total_timesteps if total_timesteps else 0.0
			for p in range(n_pairs)])

		# --- Select which pairs to demonstrate ---
		selected = self._select_pairs(
			n_pairs, pair_max, pair_smax, pair_nvalid, pair_ncells,
			valid_frac, valid_pair_reaction_ids, valid_pair_catalyst_ids,
			pair_index_by_id)
		if not selected:
			print('kcat_smoothed_max_demo: no demonstrable pairs found.')
			return

		# --- Pass 2: full traces for the selected pairs (capped seeds) ---
		seeds_all = list(ap.get_seeds())
		seeds_plot = seeds_all[:MAX_SEEDS_TO_PLOT]
		sel_idx = [p for (p, _label) in selected]
		sel_rxn = rxn_indexes[sel_idx]
		sel_cat = cat_indexes[sel_idx]

		# traces[seed_pos][k] -> dict(time, kcat, smoothed, gen_bounds)
		traces = []
		for seed in seeds_plot:
			seed_cells = sorted(ap.get_cells(
				seed=[seed],
				generation=np.arange(IGNORE_FIRST_N_GENS, ap.n_generation),
				only_successful=True))
			per_pair = [dict(time=[], kcat=[], smoothed=[], gen_bounds=[])
				for _ in sel_idx]
			t_offset = 0.0
			for cell_path in seed_cells:
				sim_out = os.path.join(cell_path, 'simOut')
				try:
					t, kcat = _read_cell_kcat(
						sim_out, sel_rxn, sel_cat, cell_density_value)
				except Exception:
					continue
				t_rel = t - t[0] + t_offset
				valid = np.isfinite(kcat) & (kcat > 0)
				for k in range(len(sel_idx)):
					col = kcat[:, k]
					vc = valid[:, k]
					sm = _smooth_cell_column(col, vc)
					raw = np.where(vc, col, np.nan)
					per_pair[k]['time'].append(t_rel)
					per_pair[k]['kcat'].append(raw)
					per_pair[k]['smoothed'].append(sm)
					per_pair[k]['gen_bounds'].append(t_offset)
				t_offset += (t[-1] - t[0])
			traces.append((seed, per_pair))

		# --- Render: one figure per selected pair -> PDF pages + HTML PNGs ---
		pdf_path = os.path.join(plotOutDir, plotOutFileName + '.pdf')
		html_sections = []
		with PdfPages(pdf_path) as pdf:
			for k, (p, label) in enumerate(selected):
				fig = self._render_reaction(
					k, p, label,
					valid_pair_reaction_ids[p], valid_pair_catalyst_ids[p],
					pair_max[p], pair_smax[p], pair_nvalid[p], pair_ncells[p],
					valid_frac[p], pair_max_cell[p], sel_rxn[k], sel_cat[k],
					cell_density_value, traces)
				pdf.savefig(fig, dpi=200)
				html_sections.append(self._figure_html(
					fig, p, label,
					valid_pair_reaction_ids[p], valid_pair_catalyst_ids[p],
					pair_max[p], pair_smax[p], pair_nvalid[p], pair_ncells[p],
					valid_frac[p]))
				plt.close(fig)
		print(f'Wrote {pdf_path}')

		html_path = os.path.join(plotOutDir, plotOutFileName + '.html')
		with open(html_path, 'w', encoding='utf-8') as fh:
			fh.write(self._html_report(html_sections, len(cell_paths),
				len(seeds_all), len(seeds_plot)))
		print(f'Wrote {html_path}')

		# Numeric cross-check line per reaction.
		for p, label in selected:
			print(f'  [{label}] {valid_pair_reaction_ids[p]} / '
				f'{valid_pair_catalyst_ids[p]}: max={pair_max[p]:.6g} '
				f'smoothed_max={_fmt_smax(pair_smax[p])} '
				f'active={valid_frac[p]:.0%}')

	# ------------------------------------------------------------------------
	# Selection
	# ------------------------------------------------------------------------

	def _select_pairs(self, n_pairs, pair_max, pair_smax, pair_nvalid,
			pair_ncells, valid_frac, rxn_ids, cat_ids, pair_index_by_id):
		"""Return an ordered list of (pair_index, type_label)."""
		# Explicit override.
		if REACTION_CATALYST_PAIRS:
			out = []
			for rxn_id, cat_id in REACTION_CATALYST_PAIRS:
				p = pair_index_by_id.get((rxn_id, cat_id))
				if p is None:
					print(f'Override pair not found / has no flux: '
						f'{rxn_id} / {cat_id}')
					continue
				out.append((p, 'specified'))
			return out

		has_signal = np.array([
			pair_max[p] > 0 and np.isfinite(pair_smax[p]) and pair_smax[p] > 0
			for p in range(n_pairs)])
		ratio = np.full(n_pairs, np.nan)
		for p in range(n_pairs):
			if has_signal[p]:
				ratio[p] = pair_smax[p] / pair_max[p]

		chosen = []
		used = set()

		def _take(cands, label, count):
			taken = 0
			for p in cands:
				if taken >= count:
					break
				if p in used:
					continue
				used.add(p)
				chosen.append((p, label))
				taken += 1

		def _prefer_above_floor(cands):
			"""Prefer candidates whose raw max clears the noise floor; fall
			back to all of them if none do (tiny test cohorts)."""
			above = [p for p in cands if pair_max[p] >= MIN_MAX_FOR_DEMO]
			return above if above else cands

		# spike-dominated: smallest ratio below the spike threshold
		spike_cands = [p for p in range(n_pairs)
			if has_signal[p] and ratio[p] <= SPIKE_RATIO_MAX]
		spike_cands.sort(key=lambda p: ratio[p])
		_take(_prefer_above_floor(spike_cands), 'spike-dominated', N_SPIKE)

		# filtered-out: positive raw max but no smoothed signal at all (every
		# valid run was shorter than the filter could keep).  Prefer the ones
		# with the most scattered samples so the erasure is visible.
		filtered_cands = [p for p in range(n_pairs)
			if np.isfinite(pair_max[p]) and pair_max[p] > 0
			and not has_signal[p]]
		filtered_cands.sort(key=lambda p: -pair_nvalid[p])
		_take(_prefer_above_floor(filtered_cands), 'filtered-out', N_FILTERED)

		# moderate-trim: partial smoothing (between the spike and clean cutoffs)
		mid = 0.5 * (SPIKE_RATIO_MAX + CLEAN_RATIO_MIN)
		moderate_cands = [p for p in range(n_pairs)
			if has_signal[p]
			and SPIKE_RATIO_MAX < ratio[p] < CLEAN_RATIO_MIN]
		moderate_cands.sort(key=lambda p: abs(ratio[p] - mid))
		_take(_prefer_above_floor(moderate_cands), 'moderate-trim', N_MODERATE)

		# intermittent: flux drops out for stretches (low valid_frac) yet the
		# sustained on-periods survive (ratio ~ 1) -- proof the zero stretches
		# do not drag smoothed_max down.  Prefer the most dropout-heavy.
		intermittent_cands = [p for p in range(n_pairs)
			if has_signal[p] and ratio[p] >= CLEAN_RATIO_MIN
			and pair_nvalid[p] >= MIN_VALID_FOR_CLEAN
			and INTERMITTENT_FRAC_MIN <= valid_frac[p] <= INTERMITTENT_FRAC_MAX]
		# Highest-magnitude one makes the most convincing demo.
		intermittent_cands.sort(key=lambda p: -pair_max[p])
		_take(_prefer_above_floor(intermittent_cands), 'intermittent',
			N_INTERMITTENT)

		# clean / sustained: ratio ~ 1 with many valid samples
		clean_cands = [p for p in range(n_pairs)
			if has_signal[p] and ratio[p] >= CLEAN_RATIO_MIN
			and pair_nvalid[p] >= MIN_VALID_FOR_CLEAN]
		clean_cands.sort(key=lambda p: -pair_nvalid[p])
		_take(_prefer_above_floor(clean_cands), 'clean/sustained', N_CLEAN)

		# sparse: fewest valid samples (but with some signal)
		sparse_cands = [p for p in range(n_pairs)
			if has_signal[p] and pair_nvalid[p] > 0]
		sparse_cands.sort(key=lambda p: pair_nvalid[p])
		_take(sparse_cands, 'sparse', N_SPARSE)

		return chosen

	# ------------------------------------------------------------------------
	# Rendering
	# ------------------------------------------------------------------------

	def _render_reaction(self, k, p, label, rxn_id, cat_id,
			max_val, smax_val, n_valid, n_cells, valid_frac, max_cell,
			rxn_idx, cat_idx, cell_density_value, traces):
		fig, axes = plt.subplots(3, 1, figsize=(13, 12))
		colors = plt.cm.tab10(np.linspace(0, 1, 10))

		# --- Panel 1: raw kcat per timestep, all plotted seeds ---
		ax = axes[0]
		for s_pos, (seed, per_pair) in enumerate(traces):
			d = per_pair[k]
			if not d['time']:
				continue
			t = np.concatenate(d['time']) / 60.
			raw = np.concatenate(d['kcat'])
			# Rasterize the dense trace (tens of thousands of points) so the
			# vector PDF stays light and renders in all viewers; axes/text
			# remain vector.  Markers make isolated valid samples visible for
			# sparse / filtered-out pairs (where a pure line draws nothing).
			ax.plot(t, raw, lw=0.6, marker='.', ms=1.5,
				color=colors[s_pos % 10], label=f'seed {seed}',
				rasterized=True)
			for gb in d['gen_bounds'][1:]:
				ax.axvline(gb / 60., color='grey', lw=0.4, ls=':', alpha=0.4)
		ax.axhline(max_val, color='red', lw=1.0, ls='--',
			label=f'raw max = {max_val:.3g}')
		ax.set_yscale('log')
		ax.set_xlabel('Time (min, concatenated over generations)')
		ax.set_ylabel('kcat (1/s, log)')
		ax.set_title(
			f'Step 1-2: raw kcat = flux/[E] per timestep, then drop '
			f'kcat<=0 / non-finite (gaps). Faint dotted = generation '
			f'boundaries.')
		ax.legend(fontsize=7, ncol=2, loc='lower right')

		# --- Panel 2: single-cell zoom on the raw-max event ---
		ax = axes[1]
		self._plot_zoom(ax, max_cell, rxn_idx, cat_idx, cell_density_value,
			max_val)

		# --- Panel 3: estimate construction ---
		ax = axes[2]
		for s_pos, (seed, per_pair) in enumerate(traces):
			d = per_pair[k]
			if not d['time']:
				continue
			t = np.concatenate(d['time']) / 60.
			sm = np.concatenate(d['smoothed'])
			ax.plot(t, sm, lw=0.7, color=colors[s_pos % 10],
				label=f'seed {seed} (smoothed)', rasterized=True)
			# per-cell maxima dots
			for cell_t, cell_sm in zip(d['time'], d['smoothed']):
				if np.any(~np.isnan(cell_sm)):
					j = int(np.nanargmax(cell_sm))
					ax.plot(cell_t[j] / 60., cell_sm[j], 'o',
						color=colors[s_pos % 10], ms=4,
						markeredgecolor='k', markeredgewidth=0.4)
		ax.axhline(max_val, color='red', lw=1.0, ls='--',
			label=f'raw max = {max_val:.3g}')
		is_filtered = not (np.isfinite(smax_val) and smax_val > 0)
		if not is_filtered:
			ax.axhline(smax_val, color='green', lw=1.5,
				label=f'smoothed_max = {smax_val:.3g}')
		ax.set_yscale('log')
		ax.set_xlabel('Time (min, concatenated over generations)')
		ax.set_ylabel('smoothed kcat (1/s, log)')
		if is_filtered:
			ax.text(0.5, 0.5,
				'smoothed_max FILTERED OUT (= 0)\nno kcat persisted for >= '
				'~%d consecutive timesteps,\nso the median filter erased every '
				'spike -> this pair contributes no bound' % SMOOTH_WINDOW,
				transform=ax.transAxes, ha='center', va='center',
				fontsize=11, color='#b30000',
				bbox=dict(boxstyle='round', fc='#fff0f0', ec='#b30000'))
			ax.set_title(
				'Step 6-7: every cell smoothed-max <= 0 -> smoothed_max = 0 '
				'(filtered out)')
		else:
			ratio = smax_val / max_val if max_val > 0 else np.nan
			ax.set_title(
				f'Step 6-7: per-cell nanmax (dots) -> max across all cells = '
				f'smoothed_max (green). smoothed_max/max = {ratio:.3g}')
		if label == 'intermittent' and not is_filtered:
			ax.text(0.5, 0.06,
				'Flux is OFF %.0f%% of the time (gaps above), yet '
				'smoothed_max/max = %.3g:\nthe zero stretches only touch the '
				'~%d steps at each dropout edge; smoothed_max is the MAX over '
				'sustained\non-periods, so dropouts do NOT drag it down.'
				% (100 * (1 - valid_frac), smax_val / max_val
					if max_val > 0 else float('nan'), SMOOTH_WINDOW // 2),
				transform=ax.transAxes, ha='center', va='bottom',
				fontsize=9, color='#0a5', bbox=dict(
					boxstyle='round', fc='#f0fff4', ec='#0a5'))
		ax.legend(fontsize=7, ncol=2, loc='lower right')

		fig.suptitle(
			f'[{label}]  {rxn_id}\ncatalyst: {cat_id}   |   '
			f'n_valid={n_valid:,} over {n_cells} cells   |   '
			f'active {valid_frac:.0%} of timesteps   |   '
			f'W={SMOOTH_WINDOW}',
			fontsize=11)
		fig.tight_layout(rect=[0, 0, 1, 0.95])
		return fig

	def _plot_zoom(self, ax, max_cell, rxn_idx, cat_idx, cell_density_value,
			max_val):
		"""Single-cell zoom showing zero-fill -> median filter -> re-mask."""
		if max_cell is None:
			ax.text(0.5, 0.5, 'no max cell', ha='center', va='center',
				transform=ax.transAxes)
			return
		try:
			_, kcat = _read_cell_kcat(
				os.path.join(max_cell, 'simOut'),
				np.array([rxn_idx]), np.array([cat_idx]), cell_density_value)
		except Exception as e:
			ax.text(0.5, 0.5, f'read error: {e!r}', ha='center', va='center',
				transform=ax.transAxes, fontsize=7)
			return
		col = kcat[:, 0]
		valid = np.isfinite(col) & (col > 0)
		filled = np.where(valid, col, 0.0)
		smoothed = median_filter(filled, size=SMOOTH_WINDOW, mode='nearest')
		smoothed_masked = smoothed.astype(float)
		smoothed_masked[~valid] = np.nan

		if np.any(valid):
			center = int(np.nanargmax(np.where(valid, col, np.nan)))
		else:
			center = len(col) // 2
		lo = max(0, center - ZOOM_HALF_WIDTH)
		hi = min(len(col), center + ZOOM_HALF_WIDTH + 1)
		x = np.arange(lo, hi)

		raw_plot = np.where(valid[lo:hi], col[lo:hi], np.nan)
		ax.plot(x, filled[lo:hi], color='#cccccc', lw=1.0,
			label='Step 3: zero-filled')
		ax.plot(x, raw_plot, 'o', color='#1f77b4', ms=3,
			label='Step 1-2: raw valid kcat')
		ax.plot(x, smoothed_masked[lo:hi], color='#2ca02c', lw=1.5,
			label=f'Step 4-5: median filter W={SMOOTH_WINDOW}, re-masked')
		ax.axvline(center, color='red', lw=0.8, ls=':')
		ax.annotate('raw max sample\n(spike)',
			xy=(center, col[center]),
			xytext=(0.98, 0.95), textcoords='axes fraction',
			ha='right', va='top', fontsize=8, color='red',
			arrowprops=dict(arrowstyle='->', color='red', lw=0.8))
		ax.set_yscale('log')
		ax.set_xlabel('Timestep index (within the cell holding the raw max)')
		ax.set_ylabel('kcat (1/s, log)')
		ax.set_title(
			'Steps 3-5: zero-fill invalid -> median filter (W=%d) -> re-mask. '
			'The filter trims spikes shorter than ~W/2.' % SMOOTH_WINDOW)
		ax.legend(fontsize=7, loc='lower right')

	# ------------------------------------------------------------------------
	# HTML report
	# ------------------------------------------------------------------------

	def _figure_html(self, fig, p, label, rxn_id, cat_id,
			max_val, smax_val, n_valid, n_cells, valid_frac):
		buf = io.BytesIO()
		fig.savefig(buf, format='png', dpi=110, bbox_inches='tight')
		buf.seek(0)
		b64 = base64.b64encode(buf.read()).decode('ascii')
		is_filtered = not (np.isfinite(smax_val) and smax_val > 0)
		ratio_str = ('&mdash;' if is_filtered
			else f'{smax_val / max_val:.3g}' if max_val > 0 else 'n/a')
		return (
			f'<section>'
			f'<h2>[{html.escape(label)}] {html.escape(rxn_id)}</h2>'
			f'<p class="cat">catalyst: {html.escape(cat_id)}</p>'
			f'<table>'
			f'<tr><th>raw max</th><th>smoothed_max</th>'
			f'<th>smoothed_max / max</th><th>active fraction</th>'
			f'<th>n_valid</th><th>n_cells</th></tr>'
			f'<tr><td>{max_val:.6g}</td><td>{html.escape(_fmt_smax(smax_val))}</td>'
			f'<td>{ratio_str}</td><td>{valid_frac:.0%}</td>'
			f'<td>{n_valid:,}</td><td>{n_cells}</td></tr>'
			f'</table>'
			f'<img src="data:image/png;base64,{b64}" />'
			f'</section>')

	def _html_report(self, sections, n_cells, n_seeds, n_seeds_plotted):
		steps = [
			('1', 'kcat = flux / [E] computed per timestep, per cell.'),
			('2', 'Drop kcat &le; 0 and non-finite samples (validity mask).'),
			('3', 'Zero-fill invalid samples so isolated spikes get pulled '
				'down by the filter.'),
			('4', f'Apply a median filter of width W={SMOOTH_WINDOW} timesteps '
				'(mode="nearest").'),
			('5', 'Re-mask the originally-invalid positions back to NaN.'),
			('6', 'Take the per-cell nanmax of the smoothed trace.'),
			('7', f'Take the nanmax across all cells (gens &ge; '
				f'{IGNORE_FIRST_N_GENS}, successful only) &rarr; '
				'<b>smoothed_max</b>.'),
		]
		steps_html = ''.join(
			f'<li><b>Step {n}:</b> {text}</li>' for n, text in steps)
		body = ''.join(sections)
		guide = self._reading_guide_html()
		return (
			'<!DOCTYPE html><html><head><meta charset="utf-8">'
			'<title>smoothed_max estimator demonstration</title>'
			'<style>'
			'body{font-family:sans-serif;margin:2em;max-width:1100px;'
			'line-height:1.5;}'
			'section{border-top:2px solid #ccc;margin-top:2em;padding-top:1em;}'
			'img{max-width:100%;height:auto;border:1px solid #eee;}'
			'table{border-collapse:collapse;margin:0.5em 0;}'
			'th,td{border:1px solid #ccc;padding:4px 10px;text-align:right;}'
			'.cat{color:#555;font-family:monospace;}'
			'.guide{background:#f6f8fa;border:1px solid #ddd;border-radius:6px;'
			'padding:0.5em 1.5em;margin:1em 0;}'
			'ol,ul{line-height:1.5;}'
			'code{background:#eef;padding:0 3px;}'
			'</style></head><body>'
			'<h1>smoothed_max (v2_smoothed_max) estimator demonstration</h1>'
			f'<p>Reproduced over {n_cells} cells from {n_seeds} seeds '
			f'(overview panels show up to {n_seeds_plotted}). The '
			'<b>smoothed_max</b> values below are identical to the production '
			'estimator in '
			'<code>kcat_estimates_v2_with_estimators.py</code>.</p>'
			f'<h3>What this is</h3>{guide}'
			f'<h3>The 7-step pipeline</h3><ol>{steps_html}</ol>'
			f'<h3>Reactions</h3>'
			f'{body}'
			'</body></html>')

	def _reading_guide_html(self):
		"""Prose explanation of how to read each per-reaction figure."""
		return (
			'<div class="guide">'
			'<p>Each reaction kcat (= reaction flux / catalyst concentration) '
			'is a noisy per-timestep signal. The plain <b>max</b> estimator '
			'would take its single highest sample as the rate bound &mdash; '
			'which a single weird timestep can inflate by orders of magnitude. '
			'The <b>smoothed_max</b> estimator instead median-filters each '
			'cell\'s trace first, so a value only survives if it persisted for '
			'&ge; ~3 consecutive timesteps. This report shows that happening '
			'for a few representative reaction <i>types</i>.</p>'
			'<p><b>Each reaction has three stacked panels:</b></p>'
			'<ul>'
			'<li><b>Top &mdash; raw kcat, all seeds.</b> Every valid '
			'per-timestep kcat sample, one color per seed, generations '
			'concatenated end to end (faint dotted lines = generation '
			'boundaries). The red dashed line is the single highest sample, '
			'i.e. what the plain <code>max</code> estimator would return. '
			'Downward excursions are timesteps where flux fell.</li>'
			'<li><b>Middle &mdash; single-cell zoom (the key view).</b> Zooms '
			'into the one cell that holds that highest sample. Blue dots = raw '
			'valid kcat; grey line = after invalid samples are zero-filled; '
			'green line = after the W=5 median filter (re-masked). If the peak '
			'was a 1&ndash;2 timestep spike, the green line stays well below '
			'it &mdash; that is the smoothing rejecting a transient the plain '
			'max would have trusted. If the peak was sustained, the green line '
			'rides up to meet it.</li>'
			'<li><b>Bottom &mdash; estimate construction.</b> The smoothed '
			'trace per seed; each dot is one cell\'s max of its smoothed trace '
			'(Step 6); the solid green line is the max of those across all '
			'cells (Step 7) = <b>smoothed_max</b>. Compare it to the red '
			'dashed <code>max</code>: a large gap means smoothing rejected '
			'spikes; green &asymp; red means the high value was real and was '
			'kept.</li>'
			'</ul>'
			'<p><b>Reaction types shown</b> (by <code>smoothed_max / max</code>):</p>'
			'<ul>'
			'<li><b>spike-dominated</b> &mdash; ratio &lt;&lt; 1: the max was '
			'a transient spike; smoothing trims it hard.</li>'
			'<li><b>moderate-trim</b> &mdash; ratio in between: the filter '
			'shaves a partially-sustained peak.</li>'
			'<li><b>clean/sustained</b> &mdash; ratio &asymp; 1: the high '
			'value persisted, so smoothed_max &asymp; max (the estimator '
			'falls back to max).</li>'
			'<li><b>sparse</b> &mdash; few valid timesteps: edge behavior '
			'where there is little to smooth.</li>'
			'<li><b>filtered-out</b> &mdash; positive raw max but '
			'<code>smoothed_max = 0</code>: the flux never persisted for '
			'&ge; ~3 consecutive timesteps, so the median filter erased every '
			'spike and the pair contributes <i>no</i> kcat bound. This is the '
			'spike-rejection at its extreme &mdash; exactly the single-weird-'
			'timestep case the estimator exists to drop.</li>'
			'<li><b>intermittent</b> &mdash; flux drops out for stretches '
			'(enzyme degraded / missing intermediate) so the reaction is '
			'active only part of the time, yet <code>smoothed_max &asymp; '
			'max</code>: the sustained on-periods survive the filter.</li>'
			'</ul>'
			'<p><b>Do the zero stretches drag the estimate down?</b> No. '
			'<code>smoothed_max</code> is a <i>max</i>, not an average. A flux '
			'dropout only pulls down the ~W/2 timesteps right at each dropout '
			'<i>edge</i> (their median window is mostly zeros); the interior '
			'of every sustained on-period is untouched, and the estimate is '
			'the max over the whole trace. So a reaction can be off most of '
			'the time and still get the correct bound from its on-periods '
			'&mdash; see the <b>intermittent</b> example, where the bottom '
			'panel\'s smoothed trace drops out during the gaps but the green '
			'<code>smoothed_max</code> line still sits at the sustained '
			'level.</p>'
			'<p><b>Sanity check:</b> smoothed_max should sit at the top of the '
			'<i>sustained</i> band of the data &mdash; not pinned to an '
			'isolated spike, and not dropped below the bulk of the trace.</p>'
			'</div>')


if __name__ == '__main__':
	Plot().cli()
