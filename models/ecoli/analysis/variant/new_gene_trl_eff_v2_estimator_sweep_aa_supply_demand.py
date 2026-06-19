"""
Amino-acid supply / demand / concentration / per-AA charging diagnosis for the
new_gene_trl_eff_v2_estimator_sweep (GFP burden across kcat regimes).

This is the listener-level companion to ..._growth_limitation.py.  It opens up the
GrowthLimits listener to ask, AA by AA, whether GFP expression (or the kcat bound)
ever leaves translation short of amino acids -- the only route by which metabolism
could feed back on growth through the ppGpp / stringent controller.

The supply-limited signature is: AA supply falling below demand (aasUsed), free AA
pools depleting, charged-tRNA fraction dropping, and ppGpp rising.  In the GFP
sweep we expect this NOT to happen (supply tracks the reduced demand, charging and
ppGpp stay ~flat, the three kcat blocks overlap) -- the positive contrast is the
wildtype kcat-tightening sweep (kcat_v2_estimator_sweep_aa_supply_demand.py), which
reuses this same machinery.

All quantities come from GrowthLimits (per-AA columns labelled by aaIds):
  aa_supply / aa_synthesis / aa_import / aa_export  -- supply rates
  aasUsed                                           -- translation demand
  aa_conc                                           -- free AA concentration
  charged_trna_conc, uncharged_trna_conc            -- per-AA charged fraction
  ppgpp_conc                                        -- stringent-response signal

Figures (all averaged over the settled window, last 8 gens, all seeds):
  <name>.pdf                     -- 2x3 overview: total supply, total demand,
                                    supply/demand balance (each rel. wildtype),
                                    mean charged-tRNA fraction, total free AA
                                    concentration (rel. wildtype), ppGpp.
  <name>_supply_composition.pdf  -- AA synthesis / import / export (absolute).
  <name>_charging_by_aa.pdf      -- per-AA charged-tRNA fraction heatmap.
  <name>_conc_by_aa.pdf          -- per-AA free AA concentration, fold vs WT.
  <name>_supplydemand_by_aa.pdf  -- per-AA supply/demand, fold vs WT.

The collect/render helpers are sweep-agnostic (driven by an AALayout) so the
wildtype sweep analysis imports and reuses them.  Reads only recorded listeners --
no re-simulation.
"""

import collections
import os

import numpy as np
from matplotlib import pyplot as plt
from matplotlib.colors import TwoSlopeNorm

from models.ecoli.analysis import variantAnalysisPlot
from models.ecoli.sim.variants.new_gene_trl_eff_v2_estimator_sweep import (
	BLOCK_ESTIMATORS,
	BLOCK_SHORT_NAMES,
	TRL_EFF_VALUES,
	is_wildtype,
	variant_to_block,
)
from wholecell.analysis.analysis_tools import exportFigure, read_stacked_columns
from wholecell.io.tablereader import TableReader


BLOCK_COLORS = ('#7f7f7f', '#1f77b4', '#2ca02c')  # gfp_only, gap, smoothed_max
WT_COLOR = '#222222'

# Sweep-agnostic plot layout shared by both the GFP and wildtype AA analyses.
#   groups      -- list of (label, color), one per line/heatmap series.
#   assign      -- {variant_index: (group_idx, x_idx)} for non-wildtype variants.
#   x_values    -- swept-parameter values (trl-eff or kcat multiplier).
#   xlabel      -- x-axis label.
#   invert_x    -- reverse the x-axis on line plots (tighter bound to the right).
#   wt          -- wildtype variant index (or None).
#   title_prefix-- sweep name for figure suptitles.
AALayout = collections.namedtuple('AALayout', [
	'groups', 'assign', 'x_values', 'xlabel', 'invert_x', 'wt', 'title_prefix'])


# ---------------------------------------------------------------------------
# Listener reading
# ---------------------------------------------------------------------------

def _read_aa_labels(cell):
	"""Return the amino-acid labels (location tag stripped) for a cell, or None."""
	try:
		ids = TableReader(
			os.path.join(cell, 'simOut', 'GrowthLimits')
			).readAttribute('aaIds')
		return [str(a).split('[')[0] for a in ids]
	except Exception:
		return None


def _per_aa_mean(cells, table, column):
	"""Return the time-mean of a per-AA column (1D, length n_aa), or None."""
	if len(cells) == 0:
		return None
	try:
		data = read_stacked_columns(
			cells, table, column, remove_first=True, ignore_exception=True)
		if data is None or data.size == 0:
			return None
		return np.nanmean(data, axis=0)
	except Exception:
		return None


def _scalar_mean(cells, table, column):
	"""Return the time-mean of a scalar column as a float, or NaN."""
	if len(cells) == 0:
		return np.nan
	try:
		data = read_stacked_columns(
			cells, table, column, remove_first=True, ignore_exception=True)
		return float(np.nanmean(data))
	except Exception:
		return np.nan


def _charged_fraction(charged, uncharged):
	"""Return per-AA charged fraction charged / (charged + uncharged), or None."""
	if charged is None or uncharged is None:
		return None
	total = charged + uncharged
	with np.errstate(divide='ignore', invalid='ignore'):
		return np.where(total > 0, charged / total, np.nan)


def _safe_ratio(numer, denom):
	"""Return elementwise numer / denom (per AA), NaN where denom is 0/invalid."""
	if numer is None or denom is None:
		return None
	with np.errstate(divide='ignore', invalid='ignore'):
		return np.where(
			(denom != 0) & np.isfinite(denom), numer / denom, np.nan)


def collect_aa_metrics(ap, variant_indexes, window):
	"""Collect per-variant AA supply/demand/concentration/charging from
	GrowthLimits, averaged over the settled window across all seeds.

	Returns a dict whose per-AA entries (supply, demand, conc, charged,
	synthesis, importr, exportr) map variant_index -> 1D array (length n_aa),
	plus scalar ppgpp (variant_index -> float) and aa_labels (list or None).
	"""
	metrics = dict(
		supply={}, demand={}, conc={}, charged={},
		synthesis={}, importr={}, exportr={}, ppgpp={})
	aa_labels = None
	for vi in variant_indexes:
		cells = ap.get_cells(
			variant=[vi], generation=window, only_successful=True)
		if aa_labels is None and len(cells) > 0:
			aa_labels = _read_aa_labels(cells[0])
		metrics['supply'][vi] = _per_aa_mean(cells, 'GrowthLimits', 'aa_supply')
		metrics['demand'][vi] = _per_aa_mean(cells, 'GrowthLimits', 'aasUsed')
		metrics['conc'][vi] = _per_aa_mean(cells, 'GrowthLimits', 'aa_conc')
		metrics['synthesis'][vi] = _per_aa_mean(
			cells, 'GrowthLimits', 'aa_synthesis')
		metrics['importr'][vi] = _per_aa_mean(
			cells, 'GrowthLimits', 'aa_import')
		metrics['exportr'][vi] = _per_aa_mean(
			cells, 'GrowthLimits', 'aa_export')
		charged = _per_aa_mean(cells, 'GrowthLimits', 'charged_trna_conc')
		uncharged = _per_aa_mean(cells, 'GrowthLimits', 'uncharged_trna_conc')
		metrics['charged'][vi] = _charged_fraction(charged, uncharged)
		metrics['ppgpp'][vi] = _scalar_mean(cells, 'GrowthLimits', 'ppgpp_conc')
	metrics['aa_labels'] = aa_labels
	return metrics


# ---------------------------------------------------------------------------
# Aggregation helpers
# ---------------------------------------------------------------------------

def _variant_list(metrics):
	"""Return the variant indexes collected (keys of a per-AA metric dict)."""
	return list(metrics['supply'].keys())


def _nan_sum(arr):
	return float(np.nansum(arr)) if arr is not None else np.nan


def _nan_mean(arr):
	return (float(np.nanmean(arr))
		if arr is not None and np.size(arr) > 0 else np.nan)


def _norm_to_wt(metric, wt):
	"""Return metric normalized so the wildtype value is 1.0 (NaN if no base)."""
	if wt is None or wt not in metric:
		return {vi: np.nan for vi in metric}
	base = metric.get(wt)
	if base is None or not np.isfinite(base) or base == 0:
		return {vi: np.nan for vi in metric}
	return {vi: metric[vi] / base for vi in metric}


def _grouped_series(metric, layout):
	"""Return {group_idx: [value-by-x_idx]} from a scalar per-variant metric."""
	n_x = len(layout.x_values)
	out = {g: [np.nan] * n_x for g in range(len(layout.groups))}
	for vi, (g, xi) in layout.assign.items():
		out[g][xi] = metric.get(vi, np.nan)
	return out


# ---------------------------------------------------------------------------
# Rendering
# ---------------------------------------------------------------------------

def _draw_lines(ax, metric, layout, ylabel, title, wt_hline=True):
	"""Plot one line per group for a scalar per-variant metric on ax."""
	series = _grouped_series(metric, layout)
	for g, (label, color) in enumerate(layout.groups):
		ax.plot(layout.x_values, series[g], marker='o', lw=1.5,
			color=color, label=label)
	if (wt_hline and layout.wt is not None
			and np.isfinite(metric.get(layout.wt, np.nan))):
		ax.axhline(metric[layout.wt], color=WT_COLOR, ls='--', lw=1.0,
			label='wildtype')
	ax.set_xlabel(layout.xlabel)
	ax.set_ylabel(ylabel)
	ax.set_title(title)
	if layout.invert_x:
		ax.invert_xaxis()


def render_overview(plotOutDir, plotOutFileName, metadata, metrics, layout):
	"""Render the 2x3 supply / demand / charging / ppGpp overview figure."""
	vlist = _variant_list(metrics)
	total_supply = {vi: _nan_sum(metrics['supply'][vi]) for vi in vlist}
	total_demand = {vi: _nan_sum(metrics['demand'][vi]) for vi in vlist}
	ratio = {}
	for vi in vlist:
		s, d = total_supply[vi], total_demand[vi]
		ratio[vi] = (s / d
			if np.isfinite(s) and np.isfinite(d) and d != 0 else np.nan)
	total_conc = {vi: _nan_sum(metrics['conc'][vi]) for vi in vlist}
	mean_charged = {vi: _nan_mean(metrics['charged'][vi]) for vi in vlist}
	ppgpp = {vi: metrics['ppgpp'][vi] for vi in vlist}

	supply_n = _norm_to_wt(total_supply, layout.wt)
	demand_n = _norm_to_wt(total_demand, layout.wt)
	ratio_n = _norm_to_wt(ratio, layout.wt)
	conc_n = _norm_to_wt(total_conc, layout.wt)

	fig, axes = plt.subplots(2, 3, figsize=(17, 9))
	_draw_lines(axes[0, 0], supply_n, layout,
		'AA supply (rel. WT)', 'Total AA supply')
	axes[0, 0].legend(fontsize=7, loc='best')
	_draw_lines(axes[0, 1], demand_n, layout,
		'AA demand / used (rel. WT)', 'Total AA demand (aasUsed)')
	_draw_lines(axes[0, 2], ratio_n, layout,
		'supply / demand (rel. WT)', 'AA supply/demand balance')
	_draw_lines(axes[1, 0], mean_charged, layout,
		'charged tRNA fraction', 'Mean charged tRNA fraction')
	_draw_lines(axes[1, 1], conc_n, layout,
		'free AA conc (rel. WT)', 'Total free AA concentration')
	_draw_lines(axes[1, 2], ppgpp, layout, 'ppGpp (uM)', 'ppGpp')

	fig.suptitle(
		f'{layout.title_prefix}: AA supply / demand / charging overview\n'
		'Supply-limited signature = supply/demand below WT, charging down, '
		'ppGpp up (per-variant means over the settled window, all seeds).',
		fontsize=12)
	fig.tight_layout(rect=[0, 0, 1, 0.94])
	exportFigure(plt, plotOutDir, plotOutFileName, metadata)
	plt.close('all')


def render_supply_composition(plotOutDir, plotOutFileName, metadata, metrics,
		layout):
	"""Render the AA synthesis / import / export composition figure (absolute).

	In minimal media synthesis dominates and import/export are ~0; the panels
	make that explicit rather than hiding it behind a normalization.
	"""
	vlist = _variant_list(metrics)
	synthesis = {vi: _nan_sum(metrics['synthesis'][vi]) for vi in vlist}
	importr = {vi: _nan_sum(metrics['importr'][vi]) for vi in vlist}
	exportr = {vi: _nan_sum(metrics['exportr'][vi]) for vi in vlist}

	fig, axes = plt.subplots(1, 3, figsize=(16, 5))
	_draw_lines(axes[0], synthesis, layout,
		'synthesis (counts/step)', 'AA synthesis')
	axes[0].legend(fontsize=7, loc='best')
	_draw_lines(axes[1], importr, layout,
		'import (counts/step)', 'AA import')
	_draw_lines(axes[2], exportr, layout,
		'export (counts/step)', 'AA export')
	fig.suptitle(
		f'{layout.title_prefix}: AA supply composition '
		'(synthesis / import / export, absolute)',
		fontsize=12)
	fig.tight_layout(rect=[0, 0, 1, 0.92])
	exportFigure(plt, plotOutDir, plotOutFileName + '_supply_composition',
		metadata)
	plt.close('all')


def _heatmap_figure(plotOutDir, plotOutFileName, metadata, per_aa_metric,
		layout, aa_labels, title, suffix, mode,
		vmin=0.5, vcenter=1.0, vmax=2.0):
	"""Render a per-group per-AA heatmap (rows = AA, cols = swept param).

	mode 'fraction': values shown directly on [0, 1] (viridis).
	mode 'fold':     values divided by the wildtype per-AA value and shown on a
	                 diverging scale centered at 1.0 (RdBu_r).
	"""
	n_groups = len(layout.groups)
	n_aa = len(aa_labels)
	n_x = len(layout.x_values)
	wt_arr = (per_aa_metric.get(layout.wt)
		if layout.wt is not None else None)

	fig, axes = plt.subplots(
		1, n_groups, figsize=(5.0 * n_groups, max(6.0, 0.32 * n_aa)),
		squeeze=False)
	if mode == 'fold':
		norm = TwoSlopeNorm(vmin=vmin, vcenter=vcenter, vmax=vmax)
	for g, (label, _color) in enumerate(layout.groups):
		grid = np.full((n_aa, n_x), np.nan)
		for vi, (gg, xi) in layout.assign.items():
			if gg != g:
				continue
			arr = per_aa_metric.get(vi)
			if arr is None:
				continue
			col = arr
			if mode == 'fold':
				if wt_arr is None:
					continue
				with np.errstate(divide='ignore', invalid='ignore'):
					col = np.where(
						(wt_arr != 0) & np.isfinite(wt_arr), arr / wt_arr,
						np.nan)
			grid[:, xi] = col
		ax = axes[0][g]
		if mode == 'fraction':
			im = ax.imshow(grid, aspect='auto', cmap='viridis',
				vmin=0.0, vmax=1.0)
		else:
			im = ax.imshow(grid, aspect='auto', cmap='RdBu_r', norm=norm)
		ax.set_xticks(range(n_x))
		ax.set_xticklabels([f'{v:g}' for v in layout.x_values],
			rotation=90, fontsize=6)
		ax.set_yticks(range(n_aa))
		ax.set_yticklabels(aa_labels, fontsize=6)
		ax.set_xlabel(layout.xlabel)
		ax.set_title(label, fontsize=10)
		fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)

	fig.suptitle(f'{layout.title_prefix}: {title}', fontsize=12)
	fig.tight_layout(rect=[0, 0, 1, 0.95])
	exportFigure(plt, plotOutDir, plotOutFileName + suffix, metadata)
	plt.close('all')


def render_per_aa_heatmaps(plotOutDir, plotOutFileName, metadata, metrics,
		layout):
	"""Render the per-AA charging, concentration, and supply/demand heatmaps."""
	aa_labels = metrics.get('aa_labels')
	if not aa_labels:
		print('aa_supply_demand: no AA labels available; skipping heatmaps.')
		return

	_heatmap_figure(plotOutDir, plotOutFileName, metadata, metrics['charged'],
		layout, aa_labels,
		'per-AA charged tRNA fraction (low = AA-limited charging)',
		'_charging_by_aa', 'fraction')

	_heatmap_figure(plotOutDir, plotOutFileName, metadata, metrics['conc'],
		layout, aa_labels,
		'per-AA free AA concentration, fold vs wildtype '
		'(blue = depletion, red = accumulation)',
		'_conc_by_aa', 'fold', vmin=0.5, vcenter=1.0, vmax=2.0)

	vlist = _variant_list(metrics)
	supply_demand = {vi: _safe_ratio(metrics['supply'][vi], metrics['demand'][vi])
		for vi in vlist}
	_heatmap_figure(plotOutDir, plotOutFileName, metadata, supply_demand,
		layout, aa_labels,
		'per-AA supply/demand, fold vs wildtype (blue = supply falling behind)',
		'_supplydemand_by_aa', 'fold', vmin=0.5, vcenter=1.0, vmax=2.0)


def render_all(plotOutDir, plotOutFileName, metadata, metrics, layout):
	"""Render every AA figure for the given metrics and layout."""
	render_overview(plotOutDir, plotOutFileName, metadata, metrics, layout)
	render_supply_composition(
		plotOutDir, plotOutFileName, metadata, metrics, layout)
	render_per_aa_heatmaps(plotOutDir, plotOutFileName, metadata, metrics, layout)


class Plot(variantAnalysisPlot.VariantAnalysisPlot):
	def do_plot(self, inputDir, plotOutDir, plotOutFileName, simDataFile,
				validationDataFile, metadata):
		variant_indexes = self.ap.get_variants()
		n_total_gens = self.ap.n_generation
		ignore_first_n_gens = max(n_total_gens - 8, 0)
		window = np.arange(ignore_first_n_gens, n_total_gens)

		metrics = collect_aa_metrics(self.ap, variant_indexes, window)

		wt = next((vi for vi in variant_indexes if is_wildtype(vi)), None)
		assign = {}
		for vi in variant_indexes:
			if is_wildtype(vi):
				continue
			block_idx, te_idx = variant_to_block(vi)
			assign[vi] = (block_idx, te_idx)
		groups = [(BLOCK_SHORT_NAMES[b], BLOCK_COLORS[b])
			for b in range(len(BLOCK_ESTIMATORS))]
		layout = AALayout(
			groups=groups, assign=assign, x_values=np.array(TRL_EFF_VALUES),
			xlabel='GFP translation efficiency', invert_x=False, wt=wt,
			title_prefix='new_gene_trl_eff_v2_estimator_sweep')

		render_all(plotOutDir, plotOutFileName, metadata, metrics, layout)


if __name__ == '__main__':
	Plot().cli()
