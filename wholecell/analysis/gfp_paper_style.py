"""
Shared drawing and saving for the GFP burden paper's figure panels.

One place for the things that must agree across all forty-odd panels: the
colour of each entity, the fonts, the panel sizes, how a design whose lineages
died is drawn, and how a panel is written to disk.

Every panel is saved as PNG, PDF and SVG with editable text, plus a JSON
sidecar recording the code and data it came from, so any figure can be traced
back and rebuilt.

Typical use inside an analysis script:

	from wholecell.analysis import gfp_paper_style as style

	style.apply()
	fig, ax = plt.subplots(figsize=style.panel_size('single'))
	...
	style.save_panel(fig, plotOutDir, 'Figure_3F_ribosomes_per_mrna',
		script=__file__, batch=inputDir, sim_data_path=simDataFile,
		variants=variants, seeds=seeds)

Colours are an Okabe-Ito derived set checked for colour-vision deficiency
across every pair, not just neighbours: the worst pair separates by dE 8.6
under simulated deutan/tritan vision and 15.6 for normal vision. Amber and
sky blue sit below 3:1 contrast on white, so they are for fills and thick
marks rather than thin lines or text.
"""

from __future__ import annotations

import json
import os
import platform
import subprocess
import sys
from datetime import datetime, timezone
from typing import Dict, Optional, Sequence, Union

import numpy as np

from wholecell.utils import filepath

# --------------------------------------------------------------------------
# Colour: one meaning per colour across the whole paper
COLORS: Dict[str, str] = {
	'gfp': '#009E73',          # the inserted gene and its product
	'rnap': '#0072B2',         # RNA polymerase
	'ribosome': '#D55E00',     # ribosomes
	'rrna': '#882255',         # rRNA operons
	'enzyme': '#E69F00',       # metabolic enzymes
	'replication': '#56B4E9',  # replication, gene dosage
	'control': '#8C8C8C',      # the no-GFP control
}

# Increasing GFP expression: one hue, light to dark. Never a rainbow.
GFP_RAMP: Sequence[str] = (
	'#B8E6D4', '#6CCBA8', '#009E73', '#00704F', '#00402D')

# Ink, for all text and axes. Series colour never carries text.
INK = '#18201D'
INK_MUTED = '#5A6661'
GRID = '#E3E6E4'

# Status colours are reserved and never reused as a series colour.
STATUS = {'good': '#1F7A4D', 'warning': '#A3480A', 'critical': '#8A1C1C'}


def gfp_ramp(n: int) -> list:
	"""n colours from the sequential GFP ramp, light to dark."""
	if n <= 1:
		return [GFP_RAMP[len(GFP_RAMP) // 2]]
	idx = np.linspace(0, len(GFP_RAMP) - 1, n)
	return [GFP_RAMP[int(round(i))] for i in idx]


# RBS strength is a magnitude, so it gets the sequential GFP ramp, and each
# RBS value keeps one colour in every panel (colour follows the entity, not
# its rank within a panel). The lightest ramp step is too pale for a line on
# white, so lines use the rest of the ramp; RBS 0 (transcribed, never
# translated) takes the pale step and is only ever drawn as a fill.
RBS_ZERO_COLOR = GFP_RAMP[0]


def rbs_color(rbs: float, all_rbs: Optional[Sequence[float]] = None) -> str:
	"""
	The colour for one RBS strength, fixed across panels.

	Args:
		rbs: translation efficiency
		all_rbs: the full set of RBS values in the paper's grid, which fixes
			each value's place on the ramp; defaults to the grid in
			gfp_paper_constants
	"""
	from matplotlib.colors import LinearSegmentedColormap, to_hex
	if rbs <= 0:
		return RBS_ZERO_COLOR
	if all_rbs is None:
		from wholecell.analysis import gfp_paper_constants
		all_rbs = gfp_paper_constants.NEW_GENE_TRANSLATION_EFFICIENCY_VALUES
	positive = np.log10(sorted(v for v in all_rbs if v > 0))
	lo, hi = positive[0], positive[-1]
	position = 0.5 if hi == lo else (np.log10(rbs) - lo) / (hi - lo)
	ramp = LinearSegmentedColormap.from_list('gfp_lines', GFP_RAMP[1:])
	return to_hex(ramp(float(np.clip(position, 0, 1))))


def promoter_axis(ax, factors: Sequence[float], label: str = 'promoter strength'):
	"""
	Promoter strength on the x axis at its true spacing.

	Expression factors are already log10 multipliers, so plotting at the
	factor itself is a log axis: 8.33 and 8.66 sit a third and two thirds of
	the way from 8 to 9, as they should.
	"""
	from wholecell.analysis import gfp_paper_constants
	ticks = sorted(set(float(f) for f in factors))
	ax.set_xticks(ticks)
	ax.set_xticklabels([gfp_paper_constants.promoter_label(f) for f in ticks])
	pad = 0.25
	ax.set_xlim(ticks[0] - pad, ticks[-1] + pad)
	ax.set_xlabel(label)


def log_axis(ax, which: str, lo: float, hi: float, data=None):
	"""
	A log axis with plain-number ticks on a 1-3-10 series.

	The range is fixed by the caller so a panel keeps the same scale from one
	batch to the next, and only widens, never narrows, if the data fall
	outside it. Ticks read 30, 100, 300 rather than 3x10^1.

	Args:
		ax: axes
		which: 'x' or 'y'
		lo, hi: default limits
		data: values plotted on this axis, used only to widen the range
	"""
	from matplotlib.ticker import FixedLocator, FuncFormatter, NullLocator
	if data is not None:
		values = np.asarray(data, dtype=float)
		values = values[np.isfinite(values) & (values > 0)]
		if values.size:
			lo = min(lo, values.min() / 1.3)
			hi = max(hi, values.max() * 1.3)
	ticks = [m * 10.0 ** e for e in range(int(np.floor(np.log10(lo))) - 1,
		int(np.ceil(np.log10(hi))) + 1) for m in (1, 3)]
	ticks = [t for t in ticks if lo <= t <= hi]
	axis = ax.xaxis if which == 'x' else ax.yaxis
	(ax.set_xscale if which == 'x' else ax.set_yscale)('log')
	(ax.set_xlim if which == 'x' else ax.set_ylim)(lo, hi)
	axis.set_major_locator(FixedLocator(ticks))
	axis.set_minor_locator(NullLocator())
	axis.set_major_formatter(FuncFormatter(lambda v, _: '{:,g}'.format(v)))


def promoter_lines(ax, designs, column: str, rbs_values: Sequence[float]):
	"""
	One line per RBS strength through the promoter strengths, for `column`.

	Args:
		ax: axes, whose x axis should then be set with promoter_axis
		designs: the design table (gfp_paper_data.design_table), induced
			designs only
		column: the quantity on the y axis
		rbs_values: RBS strengths to draw, in legend order

	Returns:
		(the rows drawn, whether any of them is a failed design)
	"""
	drawn = []
	for rbs in rbs_values:
		line = designs[np.isclose(designs['rbs'], rbs)].sort_values('promoter')
		line = line[np.isfinite(line[column])]
		if line.empty:
			continue
		color = rbs_color(rbs)
		from wholecell.analysis import gfp_paper_constants
		ax.plot(line['promoter'], line[column], color=color, zorder=2)
		design_scatter(ax, line['promoter'], line[column], color,
			failed=line['failed'], zorder=3,
			label=gfp_paper_constants.rbs_label(rbs))
		drawn.append(line)
	return drawn, any(bool(line['failed'].any()) for line in drawn)


def legend_with_failed(ax, any_failed: bool, **kwargs):
	"""The axes legend, plus the hollow-marker entry when a design failed."""
	handles, labels = ax.get_legend_handles_labels()
	if any_failed:
		handles.append(failed_legend_handle())
		labels.append(FAILED_LABEL)
	return ax.legend(handles, labels, **kwargs)


def reference_line(ax, y: float, label: str, horizontal: bool = True):
	"""A dashed physical limit, labelled at its right-hand end in muted ink."""
	if horizontal:
		ax.axhline(y, color=INK_MUTED, linestyle=(0, (4, 3)), linewidth=0.8,
			zorder=1)
		ax.annotate(label, xy=(1, y), xycoords=('axes fraction', 'data'),
			xytext=(4, 0), textcoords='offset points', va='center',
			ha='left', color=INK_MUTED, fontsize=FONT_SIZES['annotation'])
	else:
		ax.axvline(y, color=INK_MUTED, linestyle=(0, (4, 3)), linewidth=0.8,
			zorder=1)


# --------------------------------------------------------------------------
# Size: panels are drawn at final print size, so a 7 pt label really is 7 pt
MM_PER_INCH = 25.4
COLUMN_WIDTH_MM = {'single': 89.0, 'onehalf': 120.0, 'double': 183.0}
DEFAULT_ASPECT = 0.75                      # height / width


def panel_size(width: str = 'single', aspect: float = DEFAULT_ASPECT,
		fraction: float = 1.0) -> tuple:
	"""
	(width, height) in inches for a panel occupying `fraction` of a column.

	Args:
		width: 'single', 'onehalf' or 'double' journal column
		aspect: height / width
		fraction: share of the column this panel takes, for panels that sit
			side by side
	"""
	if width not in COLUMN_WIDTH_MM:
		raise ValueError('unknown column width %r; expected one of %s'
			% (width, sorted(COLUMN_WIDTH_MM)))
	inches = COLUMN_WIDTH_MM[width] / MM_PER_INCH * fraction
	return (inches, inches * aspect)


# Type sizes in points, at final size.
FONT_SIZES = {
	'tick': 7.0, 'label': 8.0, 'legend': 7.0, 'annotation': 7.0,
	'panel_letter': 9.0,
}


def apply() -> None:
	"""
	Set the rcParams every paper panel shares.

	Text stays as text in both SVG and PDF (svg.fonttype 'none' and
	pdf.fonttype 42), so panels remain editable in Illustrator instead of
	arriving as outlines.
	"""
	import matplotlib
	matplotlib.rcParams.update({
		'font.family': 'sans-serif',
		'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
		'font.size': FONT_SIZES['label'],
		'axes.labelsize': FONT_SIZES['label'],
		'axes.titlesize': FONT_SIZES['label'],
		'xtick.labelsize': FONT_SIZES['tick'],
		'ytick.labelsize': FONT_SIZES['tick'],
		'legend.fontsize': FONT_SIZES['legend'],
		'axes.labelcolor': INK,
		'axes.edgecolor': INK_MUTED,
		'text.color': INK,
		'xtick.color': INK_MUTED,
		'ytick.color': INK_MUTED,
		'axes.linewidth': 0.6,
		'xtick.major.width': 0.6,
		'ytick.major.width': 0.6,
		'axes.spines.top': False,
		'axes.spines.right': False,
		'axes.grid': False,
		'grid.color': GRID,
		'grid.linewidth': 0.5,
		'lines.linewidth': 1.2,
		'lines.markersize': 3.5,
		'legend.frameon': False,
		# Constrained layout keeps axis labels inside the figure; a tight
		# bounding box alone clips a y label taller than its axes.
		'figure.constrained_layout.use': True,
		'figure.dpi': 150,
		'savefig.dpi': 600,
		'savefig.bbox': 'tight',
		'savefig.pad_inches': 0.02,
		'mathtext.default': 'regular',  # exponents in the body font
		'svg.fonttype': 'none',        # keep SVG text editable
		'pdf.fonttype': 42,            # embed TrueType, not Type 3
		'ps.fonttype': 42,
	})


# --------------------------------------------------------------------------
# Failed lineages: one rule, applied everywhere
#
# Designs where some lineages died stay in the figures, drawn hollow, with
# "some lineages failed" in the legend. A hollow point's mean comes only from
# the lineages that survived, so it is biased toward the healthier ones; the
# survivor count goes into the panel's metadata for that reason.
FAILED_LABEL = 'some lineages failed'


def design_scatter(ax, x, y, color, failed=None, size=None, label=None,
		**kwargs):
	"""
	Scatter one design series, drawing incomplete designs as hollow markers.

	Args:
		ax: axes to draw on
		x, y: coordinates, one point per design
		color: the entity colour for this series
		failed: boolean per point, True where some lineages did not finish
		size: marker area; defaults to the rcParam marker size squared
		label: legend label for the filled points
	"""
	x = np.asarray(x, dtype=float)
	y = np.asarray(y, dtype=float)
	if failed is None:
		failed = np.zeros(x.shape, dtype=bool)
	failed = np.asarray(failed, dtype=bool)
	if failed.shape != x.shape:
		raise ValueError('failed mask has shape %s but there are %s points'
			% (failed.shape, x.shape))
	import matplotlib
	if size is None:
		size = matplotlib.rcParams['lines.markersize'] ** 2
	whole = ~failed
	ax.scatter(x[whole], y[whole], s=size, c=color, linewidths=0,
		label=label, **kwargs)
	ax.scatter(x[failed], y[failed], s=size, facecolors='none',
		edgecolors=color, linewidths=0.9, **kwargs)
	return ax


def failed_legend_handle(color=INK_MUTED):
	"""A hollow marker handle, for adding FAILED_LABEL to a legend once."""
	from matplotlib.lines import Line2D
	return Line2D([], [], linestyle='none', marker='o', markerfacecolor='none',
		markeredgecolor=color, markeredgewidth=0.9, label=FAILED_LABEL)


# --------------------------------------------------------------------------
# Statistics
#
# Fixed before the data exist, so the test is not chosen after seeing a curve.
def holm_adjust(p_values: Union[Sequence[float], np.ndarray]) -> np.ndarray:
	"""
	Holm-Bonferroni adjusted p-values, in the input order.

	Holm controls the family-wise error rate without assuming independence,
	which matters here because consecutive generations are correlated.
	"""
	p = np.asarray(p_values, dtype=float)
	n = p.size
	order = np.argsort(p)
	adjusted = np.empty(n, dtype=float)
	running = 0.0
	for rank, idx in enumerate(order):
		running = max(running, (n - rank) * p[idx])
		adjusted[idx] = min(1.0, running)
	return adjusted


def paired_divergence(treated: np.ndarray, control: np.ndarray,
		alpha: float = 0.05) -> dict:
	"""
	Per-generation paired test of treated against its seed-matched control.

	The design is seed-matched, so each seed contributes a difference and the
	test asks whether those differences are centred on zero. Used for the
	Figure 1B divergence bar.

	Args:
		treated: (n_seeds, n_generations) per-cell-cycle values with GFP
		control: (n_seeds, n_generations) the same seeds without GFP
		alpha: family-wise significance level after Holm adjustment

	Returns:
		dict with per-generation mean difference, raw and adjusted p-values,
		a significance mask, and `first_sustained`, the first generation that
		is significant and stays significant to the end (None if there is no
		such generation). The bar in 1B starts there.
	"""
	from scipy import stats
	treated = np.asarray(treated, dtype=float)
	control = np.asarray(control, dtype=float)
	if treated.shape != control.shape:
		raise ValueError('treated %s and control %s must be seed-matched'
			% (treated.shape, control.shape))
	difference = treated - control
	n_gens = difference.shape[1]
	p_raw = np.ones(n_gens)
	for gen in range(n_gens):
		column = difference[:, gen]
		column = column[np.isfinite(column)]
		if column.size < 2:
			continue
		# Real per-seed differences always vary, so this branch exists only to
		# keep a degenerate column from producing a meaningless t-statistic:
		# with no spread the t-test divides by ~0 and returns p ~ 0 or p ~ 1
		# depending on float noise alone.
		if np.allclose(column, column[0], rtol=0, atol=1e-12):
			p_raw[gen] = 1.0 if abs(column[0]) <= 1e-12 else 0.0
			continue
		p_raw[gen] = stats.ttest_1samp(column, 0.0).pvalue
	p_adj = holm_adjust(p_raw)
	significant = p_adj < alpha
	first_sustained = None
	for gen in range(n_gens):
		if significant[gen] and significant[gen:].all():
			first_sustained = gen
			break
	return {
		'mean_difference': np.nanmean(difference, axis=0),
		'p_raw': p_raw,
		'p_adjusted': p_adj,
		'significant': significant,
		'first_sustained': first_sustained,
		'n_seeds': difference.shape[0],
		'alpha': alpha,
		'test': 'paired t-test per generation, Holm-adjusted across generations',
	}


# --------------------------------------------------------------------------
# Saving, with provenance
PANEL_SUBDIR = 'paper_figures'
FORMATS = ('png', 'pdf', 'svg')


def _git_state() -> Dict[str, object]:
	try:
		commit = filepath.git_hash()
	except Exception:
		commit = 'unknown'
	try:
		dirty = bool(subprocess.run(
			['git', 'status', '--porcelain'], capture_output=True, text=True,
			cwd=os.path.dirname(os.path.abspath(__file__)),
			timeout=30).stdout.strip())
	except Exception:
		dirty = None
	return {'commit': commit, 'working_tree_dirty': dirty}


def _package_versions() -> Dict[str, str]:
	versions = {'python': sys.version.split()[0]}
	for name in ('numpy', 'scipy', 'matplotlib'):
		try:
			versions[name] = __import__(name).__version__
		except Exception:
			versions[name] = 'unavailable'
	return versions


PANELS_TSV = os.path.join(filepath.ROOT_PATH, 'runscripts', 'gfp_paper',
	'panels.tsv')
SCRIPT_PREFIX = 'gfp_paper_'


def panel_name(script: str, part: Optional[str] = None,
		manifest: str = PANELS_TSV) -> str:
	"""
	Output stem for a panel, e.g. 'Figure_3F_ribosomes_per_gfp_mrna'.

	Scripts are named for their content and the figure label is looked up in
	runscripts/gfp_paper/panels.tsv, so reordering the paper is one manifest
	edit rather than a rename of every script.

	Args:
		script: __file__ of the calling analysis
		part: suffix for scripts that write several files for one panel,
			such as the two halves of a split panel
		manifest: path to panels.tsv

	Raises:
		ValueError: if the script is not in the manifest, or is listed for
			more than one panel, since then the label would be a guess.
	"""
	base = os.path.basename(script)
	labels = []
	with open(manifest) as f:
		header = f.readline().rstrip('\n').split('\t')
		col_panel, col_script = header.index('panel'), header.index('script')
		for line in f:
			fields = line.rstrip('\n').split('\t')
			if len(fields) > col_script and os.path.basename(
					fields[col_script]) == base:
				labels.append(fields[col_panel])
	if len(labels) != 1:
		raise ValueError('%s is listed for %s in %s; expected exactly one'
			' panel' % (base, labels or 'no panel', manifest))
	stem = os.path.splitext(base)[0]
	if stem.startswith(SCRIPT_PREFIX):
		stem = stem[len(SCRIPT_PREFIX):]
	name = 'Figure_%s_%s' % (labels[0], stem)
	return name if part is None else '%s_%s' % (name, part)


def panel_dir(plot_out_dir: str) -> str:
	"""
	The one folder a batch's paper panels are written to.

	Multigen analyses normally write beside their own generation, which would
	scatter panels across the batch. Everything lands in
	<batch>/plotOut/paper_figures/ instead, so one download gets the set.
	"""
	return filepath.makedirs(plot_out_dir, PANEL_SUBDIR)


def save_panel(fig, plot_out_dir: str, name: str, script: Optional[str] = None,
		batch: Optional[str] = None, sim_data_path: Optional[str] = None,
		variants: Optional[Sequence] = None, seeds: Optional[Sequence] = None,
		generation_window: Optional[Sequence[int]] = None,
		n: Optional[int] = None, n_lineages_total: Optional[int] = None,
		n_lineages_survived: Optional[int] = None,
		notes: Optional[str] = None, subdir: Optional[str] = None,
		close: bool = True) -> str:
	"""
	Write one panel as PNG, PDF and SVG plus a JSON provenance sidecar.

	Args:
		fig: the matplotlib figure
		plot_out_dir: the analysis script's plotOutDir
		name: output stem, e.g. 'Figure_3F_ribosomes_per_mrna'
		script: __file__ of the calling analysis
		batch: the batch directory the data came from
		sim_data_path: path to the sim_data used
		variants, seeds: what went into this panel
		generation_window: (first, last) inclusive generation indices
		n: sample size behind the panel, in whatever unit the panel reports
		n_lineages_total, n_lineages_survived: recorded whenever any design
			lost lineages, since a hollow point averages survivors only
		notes: anything a reader of the sidecar would need
		subdir: optional subfolder, for panels with one file per variant/seed
		close: close the figure afterwards

	Returns:
		The directory the files were written to.
	"""
	_warn_long_labels(fig, name)
	out_dir = panel_dir(plot_out_dir)
	if subdir:
		out_dir = filepath.makedirs(out_dir, subdir)
	for extension in FORMATS:
		fig.savefig(os.path.join(out_dir, '%s.%s' % (name, extension)))

	metadata = {
		'panel': name,
		'written_utc': datetime.now(timezone.utc).isoformat(timespec='seconds'),
		'script': os.path.basename(script) if script else None,
		'batch': batch,
		'sim_data': sim_data_path,
		'sim_data_sha256': _sha256(sim_data_path) if sim_data_path else None,
		'variants': list(variants) if variants is not None else None,
		'seeds': list(seeds) if seeds is not None else None,
		'generation_window': (list(generation_window)
			if generation_window is not None else None),
		'n': n,
		'n_lineages_total': n_lineages_total,
		'n_lineages_survived': n_lineages_survived,
		'git': _git_state(),
		'packages': _package_versions(),
		'host': {'node': platform.node(), 'processor': platform.processor(),
			'platform': platform.platform()},
		'notes': notes,
	}
	with open(os.path.join(out_dir, '%s.json' % name), 'w') as f:
		json.dump(metadata, f, indent=1, sort_keys=True, default=str)

	if close:
		import matplotlib.pyplot as plt
		plt.close(fig)
	return out_dir


def _warn_long_labels(fig, name: str) -> None:
	"""Warn when a y label is longer than its axes, which reads badly."""
	import warnings
	fig.canvas.draw()
	renderer = fig.canvas.get_renderer()
	for ax in fig.axes:
		label = ax.yaxis.label
		if not label.get_text():
			continue
		if label.get_window_extent(renderer).height > 1.05 * ax.bbox.height:
			warnings.warn('%s: y label %r is taller than its axes; shorten it'
				' or give the panel more height' % (name, label.get_text()))


def _sha256(path: str) -> Optional[str]:
	"""Hash of the sim_data pickle, so two batches can be told apart."""
	import hashlib
	try:
		digest = hashlib.sha256()
		with open(path, 'rb') as f:
			for block in iter(lambda: f.read(1 << 20), b''):
				digest.update(block)
		return digest.hexdigest()
	except Exception:
		return None
