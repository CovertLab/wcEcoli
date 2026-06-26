"""
Estimate implied kcats for EVERY catalyzed reaction (not just the curated
below-line list), then compare them to the experimentally measured kcats in the
model.

Intended for a no-kinetics run (the metabolism_kinetic_objective_weight index-0
sim, kinetic_objective_weight = 0): with the kinetic objective off, FBA flux is
set by homeostatic mass balance rather than by any kcat bound, so
implied_kcat = flux / [enzyme] is an unbiased "what turnover would this enzyme
need to support the flux it actually carried" -- not a value the bound itself
shaped.

This reuses the established v2 estimator logic
(kcat_estimates_v2_with_estimators.py) but swaps the curated-CSV reaction
selection for the full set of catalyzed reactions from
sim_data.process.metabolism.reaction_catalysts.

Per (reaction, catalyst) pair, over the settled window (last gens, all seeds), it
computes the implied kcat in BOTH unit bases:
  * per-enzyme turnover, 1/s          (= raw reactionFluxes / [enzyme]);
                                       directly comparable to literature kcats.
  * effective capacity, L/g DCW/h     (= GDCW-converted flux / [enzyme]);
                                       matches the existing v2 estimator TSVs and
                                       is pluggable as a model kcat bound.
and reduces the per-timestep samples with the v2 estimators (no others):
  max, gap (v2 gap detector), smoothed_max (v2 median-filter W=5), and
  selected = max(gap, smoothed_max).

Multiple catalysts (isozymes): one (reaction, catalyst) row per enzyme, each =
flux / THAT enzyme's concentration -- i.e. the kcat the enzyme would need if it
alone carried the reaction's flux (an UPPER bound when there are isozymes).  Rows
whose reaction has >1 catalyst are flagged (multi_catalyst) so the upper-bound
caveat is visible.

Filtering: estimates that are non-finite, <= 0, below a floor (MIN_KCAT_*), or
backed by too few active samples (MIN_SAMPLES) are dropped from the TSVs and the
comparison (the summary CSV keeps every pair with its n_samples).

Outputs (in plotOutDir):
  kcat_estimates_all_catalysts_summary.csv          -- one row per pair, both
      unit bases x {max, gap, smoothed_max, selected}, n_samples, multi-catalyst
      flag, and the matched measured kcat (1/s) where one exists.
  kcat_estimates_all_catalysts_{max,gap,smoothed_max,selected}.tsv  -- L/g DCW/h,
      in the existing v2 TSV schema; drop into
      reconstruction/ecoli/flat/kcat_estimates/ to use as model bounds.
  kcat_vs_measured.csv / kcat_vs_measured.pdf       -- implied (selected, 1/s) vs
      measured (1/s) for the pairs present in both, with the implied/measured
      ratio distribution.

Reads only recorded listeners (FBAResults, EnzymeKinetics, Mass) + the measured
kcat flat file -- no re-simulation.
"""

import csv
import datetime
import json
import os

import numpy as np
from scipy.ndimage import median_filter
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from models.ecoli.analysis import cohortAnalysisPlot
from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from models.ecoli.analysis.cohort.kcat_estimates_v2_with_estimators import (
	IGNORE_FIRST_N_GENS,
	SMOOTH_WINDOW,
	TOP_K_GAP_CAP,
	QUANTITY_BINS,
	_gap_estimate,
	_merge_topk,
	_percentile_from_hist,
)
from models.ecoli.processes.metabolism import (
	COUNTS_UNITS, VOLUME_UNITS, TIME_UNITS, MASS_UNITS)
from reconstruction.ecoli.dataclasses.process.metabolism import (
	MIN_KCAT_ESTIMATE)
from reconstruction.spreadsheets import CSV_DIALECT, JsonWriter
from wholecell.io.tablereader import TableReader
from wholecell.utils import units


# Per-enzyme-turnover floor (1/s) below which an implied kcat is treated as
# near-zero-flux noise and dropped from the comparison.
MIN_KCAT_1S = 0.01
# Effective-capacity floor (L/g DCW/h); mirrors the existing pipeline so the
# TSVs match how MIN_KCAT_ESTIMATE filters kcat_estimates at parca load time.
MIN_KCAT_LGH = MIN_KCAT_ESTIMATE
# Minimum active samples (flux > 0 and [enzyme] > 0) for a pair to be trusted.
MIN_SAMPLES = 20

# Cap on retained top-K kcat samples per pair per basis (bounds memory for the
# full catalyst set; the v2 default 5000 x both bases x thousands of pairs is
# unnecessarily large).  k_eff in _gap_estimate still scales with sample count.
TOP_K_CAP = min(TOP_K_GAP_CAP, 2000)

# Unit bases: label -> human description for provenance.
BASES = ('per_s', 'Lgh')

# Measured kcat flat file (literature turnover numbers, 1/s).
_REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(
	os.path.dirname(os.path.abspath(__file__))))))
METABOLISM_KINETICS_TSV = os.path.join(
	_REPO_ROOT, 'reconstruction', 'ecoli', 'flat', 'metabolism_kinetics.tsv')

REVERSE_TAG = ' (reverse)'


def _strip_location(catalyst_id):
	"""Drop a trailing location tag, e.g. 'CPLX0-235[c]' -> 'CPLX0-235'."""
	return catalyst_id.split('[')[0]


def _base_reaction(reaction_id):
	"""Normalize a model reaction id to the measured-file namespace.

	Removes the ' (reverse)' direction tag and any '__ENZYME' suffix added when a
	multi-enzyme reaction is split, so model reaction ids can be matched against
	the base reactionID column in metabolism_kinetics.tsv.
	"""
	base = reaction_id
	if base.endswith(REVERSE_TAG):
		base = base[:-len(REVERSE_TAG)]
	if '__' in base:
		base = base.split('__')[0]
	return base


def _temperature_adjusted(kcat_value, temp):
	"""Adjust a measured kcat (1/s) to 37 C with Q10 = 2 (blank temp -> 25).

	Mirrors Metabolism.temperature_adjusted_kcat without needing a units.Unum.
	"""
	if temp is None or temp == '':
		temp = 25.0
	return 2 ** ((37.0 - float(temp)) / 10.0) * float(kcat_value)


def _load_measured_kcats():
	"""Return {(base_reaction, base_enzyme): kcat_1s} from the measured file.

	kcat (1/s) lists are temperature-adjusted to 37 C; when a row lists several
	kcats (one per substrate) the max is kept as the enzyme's peak turnover, and
	when several rows map to the same key (e.g. forward/reverse directions) the
	max across them is kept.
	"""
	measured = {}
	if not os.path.exists(METABOLISM_KINETICS_TSV):
		print(f'Measured kcat file not found: {METABOLISM_KINETICS_TSV}')
		return measured
	with open(METABOLISM_KINETICS_TSV, newline='', encoding='utf-8') as fh:
		lines = [ln for ln in fh if not ln.startswith('#')]
	reader = csv.DictReader(lines, delimiter='\t')
	kcat_col = next((c for c in reader.fieldnames if c.startswith('kcat')), None)
	for row in reader:
		raw = (row.get(kcat_col) or '').strip()
		if not raw:
			continue
		try:
			kcats = json.loads(raw)
		except (ValueError, TypeError):
			continue
		if not kcats:
			continue
		temp = (row.get('Temp') or '').strip()
		adj = max(_temperature_adjusted(k, temp) for k in kcats)
		rxn = _base_reaction(str(row['reactionID']).strip('"'))
		enz = _strip_location(str(row['enzymeID']).strip('"'))
		key = (rxn, enz)
		if key not in measured or adj > measured[key]:
			measured[key] = adj
	return measured


def _select(gap_est, smoothed_est, max_val):
	"""Recommended estimate = max(gap, smoothed_max), else max (v2 rule)."""
	finite = [v for v in (gap_est, smoothed_est)
		if v is not None and np.isfinite(v) and v > 0]
	if finite:
		return max(finite)
	return float(max_val) if np.isfinite(max_val) else np.nan


class Plot(cohortAnalysisPlot.CohortAnalysisPlot):
	def do_plot(self, variantDir, plotOutDir, plotOutFileName, simDataFile,
				validationDataFile, metadata):
		import pickle
		with open(simDataFile, 'rb') as f:
			sim_data = pickle.load(f)

		ap = AnalysisPaths(variantDir, cohort_plot=True)
		cell_paths = ap.get_cells(
			generation=np.arange(IGNORE_FIRST_N_GENS, ap.n_generation),
			only_successful=True)
		if len(cell_paths) == 0:
			print('Skipping -- not enough simulations run.')
			return

		# --- Enumerate ALL (reaction, catalyst) pairs present in FBAResults ---
		fba_reader = TableReader(
			os.path.join(cell_paths[0], 'simOut', 'FBAResults'))
		listener_rxn_ids = fba_reader.readAttribute('reactionIDs')
		listener_cat_ids = fba_reader.readAttribute('catalyst_ids')
		rxn_id_to_idx = {r: i for i, r in enumerate(listener_rxn_ids)}
		cat_id_to_idx = {c: i for i, c in enumerate(listener_cat_ids)}

		reaction_catalysts = sim_data.process.metabolism.reaction_catalysts

		rxn_indexes, cat_indexes = [], []
		pair_reaction_ids, pair_catalyst_ids = [], []
		pair_n_catalysts = []
		for reaction_id, cat_ids in reaction_catalysts.items():
			if reaction_id not in rxn_id_to_idx:
				continue
			for catalyst_id in cat_ids:
				if catalyst_id not in cat_id_to_idx:
					continue
				rxn_indexes.append(rxn_id_to_idx[reaction_id])
				cat_indexes.append(cat_id_to_idx[catalyst_id])
				pair_reaction_ids.append(reaction_id)
				pair_catalyst_ids.append(catalyst_id)
				pair_n_catalysts.append(len(cat_ids))
		rxn_indexes = np.array(rxn_indexes)
		cat_indexes = np.array(cat_indexes)
		n_pairs = len(rxn_indexes)
		print(f'{n_pairs} reaction-catalyst pairs (all catalysts) across '
			f'{len(cell_paths)} cells.')
		if n_pairs == 0:
			print('No catalyzed reactions found in FBAResults.')
			return

		cell_density_value = sim_data.constants.cell_density.asNumber(
			MASS_UNITS / VOLUME_UNITS)
		edges = QUANTITY_BINS['kcat'][0]
		n_bins = len(edges) - 1

		# --- Accumulators, one independent set per unit basis ---
		hist = {b: np.zeros((n_pairs, n_bins), dtype=np.int64) for b in BASES}
		stat_max = {b: np.full(n_pairs, -np.inf) for b in BASES}
		top_k = {b: [np.empty(0, dtype=float) for _ in range(n_pairs)]
			for b in BASES}
		smoothed_max = {b: np.full(n_pairs, -np.inf) for b in BASES}
		n_active = np.zeros(n_pairs, dtype=np.int64)  # shared across bases
		conv_sum, conv_n = 0.0, 0

		pair_arange = np.arange(n_pairs)

		def _accumulate(basis, kcat_arr, valid):
			"""Fold one cell's (T, n_pairs) kcat samples into basis accumulators."""
			# Histogram (for the gap p99 floor) + running max.
			vals = kcat_arr[valid]
			pairs = np.broadcast_to(pair_arange, kcat_arr.shape)[valid]
			if vals.size:
				bins = np.searchsorted(edges, vals, side='right') - 1
				np.clip(bins, 0, n_bins - 1, out=bins)
				np.add.at(hist[basis], (pairs, bins), 1)
			masked = np.where(valid, kcat_arr, -np.inf)
			cell_max = masked.max(axis=0)
			np.maximum(stat_max[basis], cell_max, out=stat_max[basis])
			# Only pairs active in this cell need the per-pair loops below;
			# inactive pairs (no flux / no enzyme) contribute nothing.  With the
			# full catalyst set this is ~900 of ~18000 pairs per cell.
			active_pairs = np.where(valid.any(axis=0))[0]
			# Top-K for the gap detector.
			if kcat_arr.shape[0] >= TOP_K_CAP:
				part = np.argpartition(
					masked, -TOP_K_CAP, axis=0)[-TOP_K_CAP:, :]
				cell_top = np.take_along_axis(masked, part, axis=0)
			else:
				cell_top = masked
			for p in active_pairs:
				col = cell_top[:, p]
				col = col[np.isfinite(col)]
				if col.size:
					top_k[basis][p] = _merge_topk(top_k[basis][p], col, TOP_K_CAP)
			# Smoothed-max: per-pair median filter, mask invalid, nanmax.
			kcat_for_filter = np.where(valid, kcat_arr, 0.0)
			for p in active_pairs:
				vp = valid[:, p]
				sm = median_filter(
					kcat_for_filter[:, p], size=SMOOTH_WINDOW, mode='nearest')
				sm[~vp] = np.nan
				cm = float(np.nanmax(sm)) if np.any(~np.isnan(sm)) else -np.inf
				if cm <= 0:
					cm = -np.inf
				if cm > smoothed_max[basis][p]:
					smoothed_max[basis][p] = cm

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
				print(f'Ignored exception reading {sim_out}: {e!r}')
				continue

			conv = dry_mass / cell_mass * cell_density_value  # (T,) g/L
			if counts_to_molar.ndim > 1:
				conc = catalyst_counts * counts_to_molar
			else:
				conc = catalyst_counts * counts_to_molar[:, np.newaxis]

			# per-enzyme turnover (1/s): raw flux (mmol/L/s) / conc (mmol/L).
			with np.errstate(divide='ignore', invalid='ignore'):
				kcat_s = reaction_fluxes / conc
			kcat_s[~np.isfinite(kcat_s)] = np.nan

			# effective capacity (L/g DCW/h): GDCW-converted flux / conc.
			flux_gdcw = (
				(COUNTS_UNITS / MASS_UNITS / TIME_UNITS)
				* (reaction_fluxes / conv[:, np.newaxis])
			).asNumber(units.mmol / units.g / units.h)
			with np.errstate(divide='ignore', invalid='ignore'):
				kcat_l = flux_gdcw / conc
			kcat_l[~np.isfinite(kcat_l)] = np.nan

			valid = np.isfinite(kcat_s) & (kcat_s > 0)
			n_active += valid.sum(axis=0)
			conv_sum += float(np.nansum(conv))
			conv_n += int(np.isfinite(conv).sum())

			_accumulate('per_s', kcat_s, valid)
			_accumulate('Lgh', kcat_l, valid)
			cells_processed += 1

		if cells_processed == 0:
			print('No cells successfully processed.')
			return
		print(f'Processed {cells_processed} cells.')
		conv_mean = conv_sum / conv_n if conv_n else np.nan

		# --- Reduce to per-pair estimates (both bases) ---
		def _estimates(basis, p, n):
			max_val = float(stat_max[basis][p])
			if n == 0 or not np.isfinite(max_val):
				return dict(max=np.nan, gap=np.nan, smoothed_max=np.nan,
					selected=np.nan)
			p99 = _percentile_from_hist(hist[basis][p], edges, 0.99)
			gap = _gap_estimate(top_k[basis][p], n, max_val,
				p99_floor=p99 if np.isfinite(p99) else None)
			sm = smoothed_max[basis][p]
			sm = float(sm) if np.isfinite(sm) and sm > 0 else np.nan
			return dict(max=max_val, gap=gap, smoothed_max=sm,
				selected=_select(gap, sm, max_val))

		measured = _load_measured_kcats()
		print(f'{len(measured)} measured (reaction, enzyme) kcats loaded.')

		rows = []
		for p in range(n_pairs):
			n = int(n_active[p])
			est = {b: _estimates(b, p, n) for b in BASES}
			rxn_id = pair_reaction_ids[p]
			cat_id = pair_catalyst_ids[p]
			meas = measured.get(
				(_base_reaction(rxn_id), _strip_location(cat_id)), np.nan)
			rows.append({
				'reaction_id': rxn_id,
				'catalyst_id': cat_id,
				'n_catalysts': pair_n_catalysts[p],
				'multi_catalyst': int(pair_n_catalysts[p] > 1),
				'n_samples': n,
				'per_s': est['per_s'],
				'Lgh': est['Lgh'],
				'measured_kcat_1s': meas,
			})

		# Sort by selected per-enzyme turnover, descending.
		rows.sort(
			key=lambda r: (r['per_s']['selected']
				if np.isfinite(r['per_s']['selected']) else -np.inf),
			reverse=True)

		timestamp = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
		n_gens_used = ap.n_generation - IGNORE_FIRST_N_GENS
		self._write_summary(plotOutDir, rows, timestamp, ap, n_gens_used,
			cells_processed, len(cell_paths), conv_mean)
		self._write_tsvs(plotOutDir, rows, timestamp, variantDir)
		self._write_comparison(plotOutDir, rows, timestamp)

	# ------------------------------------------------------------------ outputs
	@staticmethod
	def _write_summary(plotOutDir, rows, timestamp, ap, n_gens_used,
			cells_processed, n_cells, conv_mean):
		ests = ('max', 'gap', 'smoothed_max', 'selected')
		header = (['reaction_id', 'catalyst_id', 'n_catalysts', 'multi_catalyst',
			'n_samples']
			+ [f'kcat_per_s_{e}' for e in ests]
			+ [f'kcat_Lgh_{e}' for e in ests]
			+ ['measured_kcat_1s', 'selected_per_s_over_measured'])

		def _fmt(x):
			return f'{x:.6g}' if (x is not None and np.isfinite(x)) else ''

		path = os.path.join(
			plotOutDir, 'kcat_estimates_all_catalysts_summary.csv')
		with open(path, 'w', newline='', encoding='utf-8') as fh:
			fh.write(
				f'# Generated by kcat_estimates_all_catalysts.py on {timestamp}.'
				f' {ap.n_seed} seed(s), gens {IGNORE_FIRST_N_GENS}-'
				f'{ap.n_generation - 1} ({n_gens_used} gens),'
				f' {cells_processed}/{n_cells} cells.'
				f' mean GDCW conversion (dryMass/cellMass*density) = '
				f'{conv_mean:.6g} g/L.'
				f' per-enzyme kcat in 1/s; effective kcat in L/g DCW/h.\n')
			writer = csv.writer(fh)
			writer.writerow(header)
			# Skip pairs that never carried flux (no implied kcat to report).
			for r in (r for r in rows if r['n_samples'] > 0):
				meas = r['measured_kcat_1s']
				sel = r['per_s']['selected']
				ratio = (sel / meas if np.isfinite(meas) and meas > 0
					and np.isfinite(sel) else np.nan)
				writer.writerow(
					[r['reaction_id'], r['catalyst_id'], r['n_catalysts'],
						r['multi_catalyst'], r['n_samples']]
					+ [_fmt(r['per_s'][e]) for e in ests]
					+ [_fmt(r['Lgh'][e]) for e in ests]
					+ [_fmt(meas), _fmt(ratio)])
		print(f'Wrote {path} ({len(rows)} pairs)')

	@staticmethod
	def _write_tsvs(plotOutDir, rows, timestamp, variantDir):
		"""Emit one L/g DCW/h TSV per estimator in the v2 flat-file schema."""
		for est in ('max', 'gap', 'smoothed_max', 'selected'):
			valid = [r for r in rows
				if np.isfinite(r['Lgh'][est]) and r['Lgh'][est] >= MIN_KCAT_LGH
				and r['n_samples'] >= MIN_SAMPLES]
			out_path = os.path.join(
				plotOutDir, f'kcat_estimates_all_catalysts_{est}.tsv')
			with open(out_path, 'w', newline='', encoding='utf-8') as fh:
				fh.write(
					f'# Generated by kcat_estimates_all_catalysts.py from '
					f'{os.path.basename(os.path.normpath(variantDir))} on '
					f'{timestamp}. Estimator: {est} (L/g DCW/h). '
					f'{len(valid)} pairs with kcat >= {MIN_KCAT_LGH} and '
					f'>= {MIN_SAMPLES} samples.\n')
				writer = JsonWriter(
					fh, ['reaction_id', 'catalyst_id', 'kcat_estimate'],
					dialect=CSV_DIALECT)
				writer.writeheader()
				for r in valid:
					writer.writerow({
						'reaction_id': r['reaction_id'],
						'catalyst_id': r['catalyst_id'],
						'kcat_estimate': r['Lgh'][est]})
			print(f'Wrote {out_path} ({len(valid)} rows)')

	@staticmethod
	def _write_comparison(plotOutDir, rows, timestamp):
		"""Compare implied (selected, 1/s) vs measured (1/s): CSV + scatter."""
		pairs = [r for r in rows
			if np.isfinite(r['measured_kcat_1s']) and r['measured_kcat_1s'] > 0
			and np.isfinite(r['per_s']['selected'])
			and r['per_s']['selected'] >= MIN_KCAT_1S
			and r['n_samples'] >= MIN_SAMPLES]

		csv_path = os.path.join(plotOutDir, 'kcat_vs_measured.csv')
		with open(csv_path, 'w', newline='', encoding='utf-8') as fh:
			fh.write(f'# implied (selected, 1/s) vs measured (1/s, 37C) on '
				f'{timestamp}. {len(pairs)} matched pairs.\n')
			writer = csv.writer(fh)
			writer.writerow(['reaction_id', 'catalyst_id', 'measured_kcat_1s',
				'implied_kcat_1s', 'implied_over_measured', 'n_samples',
				'multi_catalyst'])
			for r in pairs:
				meas = r['measured_kcat_1s']
				imp = r['per_s']['selected']
				writer.writerow([r['reaction_id'], r['catalyst_id'],
					f'{meas:.6g}', f'{imp:.6g}', f'{imp / meas:.6g}',
					r['n_samples'], r['multi_catalyst']])
		print(f'Wrote {csv_path} ({len(pairs)} matched pairs)')

		if not pairs:
			print('No measured/implied matches -- skipping scatter.')
			return
		meas = np.array([r['measured_kcat_1s'] for r in pairs])
		imp = np.array([r['per_s']['selected'] for r in pairs])
		multi = np.array([bool(r['multi_catalyst']) for r in pairs])
		ratio = imp / meas

		fig, ax = plt.subplots(figsize=(7, 7))
		lo = min(meas.min(), imp.min()) * 0.5
		hi = max(meas.max(), imp.max()) * 2.0
		ax.plot([lo, hi], [lo, hi], ls='--', color='#888', lw=1.0, label='y = x')
		ax.scatter(meas[~multi], imp[~multi], s=18, alpha=0.6, color='#1f77b4',
			label='single catalyst')
		ax.scatter(meas[multi], imp[multi], s=18, alpha=0.6, color='#d62728',
			label='multi-catalyst (upper bound)')
		ax.set_xscale('log')
		ax.set_yscale('log')
		ax.set_xlim(lo, hi)
		ax.set_ylim(lo, hi)
		ax.set_xlabel('measured kcat (1/s, 37C)')
		ax.set_ylabel('implied kcat (selected, 1/s)')
		med = float(np.median(ratio))
		q1, q3 = np.percentile(ratio, [25, 75])
		ax.set_title(
			'Implied vs measured kcat (no-kinetics run)\n'
			f'n = {len(pairs)} pairs; implied/measured median {med:.2g} '
			f'(IQR {q1:.2g}-{q3:.2g})')
		ax.legend(fontsize=8, loc='best')
		fig.tight_layout()
		pdf_path = os.path.join(plotOutDir, 'kcat_vs_measured.pdf')
		fig.savefig(pdf_path, dpi=200)
		plt.close('all')
		print(f'Wrote {pdf_path}')


if __name__ == '__main__':
	Plot().cli()
