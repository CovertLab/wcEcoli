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
from typing import Dict, Optional, Sequence

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
		'figure.dpi': 150,
		'savefig.dpi': 600,
		'savefig.bbox': 'tight',
		'savefig.pad_inches': 0.02,
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
def holm_adjust(p_values: Sequence[float]) -> np.ndarray:
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
