"""
Plot translation efficiency vs every numeric y-metric from the scatterplot
gallery CSVs. One PNG per y-column is written to
plotOutDir/trl_eff_vs_all/.

x encoding:
  - expression variants: x = TRL_EFF_VALUES[local_idx - 1]
  - controls:            x = -0.5 (a sentinel placement)

Visual encoding (shared with the other trl_eff_sweep scripts):
  - Color -> kcat category (CATEGORY_COLORS[cat_idx])
  - Marker -> EXPRESSION_MARKER for expression variants, CONTROL_MARKER for
    controls.

Also emits the 5 multi-y overlays (rRNA, RNAP, ribosome subunits, limiting
30s, limiting 50s) and the 4 ratio plots (active/total ribosome, etc.) the
original standalone script produced.

Reads CSVs produced by new_gene_counts_save_dt / save_omes /
save_transcriptome (located in plotOutDir).
"""

import os

import matplotlib.lines as mlines
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from models.ecoli.analysis import variantAnalysisPlot
from models.ecoli.sim.variants.new_gene_trl_eff_sweep import (
	TRL_EFF_VALUES,
	KCAT_MULTIPLIERS,
	category_label,
	is_control,
	variant_to_category,
)

# x sentinel for controls
CONTROL_X = -0.5

# Canonical category visual encoding
CATEGORY_COLORS = {
	0: (27/255, 132/255, 198/255),   # no_kcat -- blue
	1: (202/255,   0/255,  32/255),  # 1.0x    -- red
	2: (230/255, 120/255,  25/255),  # 0.8x    -- orange
	3: (120/255,  60/255, 155/255),  # 0.6x    -- purple
}
CONTROL_MARKER = '^'
EXPRESSION_MARKER = 'o'

MULTI_Y_COLORS = [
	(136/255, 205/255, 240/255),
	(188/255, 140/255, 191/255),
	(66/255, 170/255, 154/255),
	(221/255, 203/255, 119/255),
	(27/255, 132/255, 198/255),
]

MULTI_Y_GROUPS = [
	('rRNA Counts', [
		'Avg Total Ribosome Counts',
		'Avg rRNA Counts: RRSA-RRNA[c]',
		'Avg rRNA Counts: RRLA-RRNA[c]',
		'Avg rRNA Counts: RRFA-RRNA[c]',
	]),
	('RNA Polymerase Counts', [
		'Avg Total RNA Polymerase Counts',
		'Avg Monomer Counts (RNAP): EG10893-MONOMER[c]',
		'Avg Monomer Counts (RNAP): RPOC-MONOMER[c]',
		'Avg Monomer Counts (RNAP): RPOB-MONOMER[c]',
	]),
	('Ribosomal Subunit Counts', [
		'Avg Total Ribosome Counts',
		'Avg Total Counts (Ribosomal 30s)',
		'Avg Total Counts (Ribosomal 50s)',
	]),
	('Limiting Counts Ribosomal 30s', [
		'Avg Total Counts (Ribosomal 30s)',
		'Avg Limiting Protein Counts (Ribosomal 30s)',
		'Avg 16s rRNA Counts (Ribosomal 30s)',
	]),
	('Limiting Counts Ribosomal 50s', [
		'Avg Total Counts (Ribosomal 50s)',
		'Avg Limiting Protein Counts (Ribosomal 50s)',
		'Avg 23s rRNA Counts (Ribosomal 50s)',
		'Avg 5s rRNA Counts (Ribosomal 50s)',
	]),
]

RATIO_PAIRS = [
	('Avg Active Ribosome Counts', 'Avg Total Ribosome Counts'),
	('Avg Active RNA Polymerase Counts', 'Avg Total RNA Polymerase Counts'),
	('Avg Total Ribosome Counts', 'Avg Total RNA Polymerase Counts'),
	('Avg Total Ribosome Concentration (uM)',
	 'Avg Total RNA Polymerase Concentration (uM)'),
]

SKIP_COLS = {'Variant Index', 'trl_eff', 'is_ctrl', 'cat_idx'}


def _safe_filename(label):
	"""Strip characters that are unsafe in filenames."""
	return (label.replace(' ', '_')
		.replace('/', '_')
		.replace(':', '')
		.replace('\\', '_'))


def _variant_to_x(index):
	"""Return x placement for a variant index."""
	cat_idx, local_idx = variant_to_category(int(index))
	if local_idx == 0:
		return CONTROL_X
	return TRL_EFF_VALUES[local_idx - 1]


def _group_key(cat_idx, is_ctrl):
	suffix = 'control' if is_ctrl else 'expression'
	return f'{category_label(cat_idx)} {suffix}'


def _group_style(cat_idx, is_ctrl):
	marker = CONTROL_MARKER if is_ctrl else EXPRESSION_MARKER
	return CATEGORY_COLORS[cat_idx], marker


def _group_masks(df):
	"""Return a dict mapping group label -> boolean mask."""
	cats = df['cat_idx'].astype(int).values
	ctrls = df['is_ctrl'].astype(bool).values
	masks = {}
	for cat_idx in range(len(KCAT_MULTIPLIERS)):
		masks[_group_key(cat_idx, False)] = (cats == cat_idx) & ~ctrls
		masks[_group_key(cat_idx, True)] = (cats == cat_idx) & ctrls
	return masks


def _add_legend(ax):
	"""Add a legend with one entry per category x (expression, control)."""
	handles = []
	for cat_idx in range(len(KCAT_MULTIPLIERS)):
		for is_ctrl in (False, True):
			color, marker = _group_style(cat_idx, is_ctrl)
			handles.append(mlines.Line2D(
				[], [], color=color, marker=marker, linestyle='None',
				markersize=6, alpha=0.7,
				label=_group_key(cat_idx, is_ctrl)))
	ax.legend(handles=handles, fontsize=7, loc='best')


def _plot_single(df, y_col, out_dir):
	fig, ax = plt.subplots(figsize=(10, 8))
	x = df['trl_eff'].values
	y = df[y_col].values
	for label, mask in _group_masks(df).items():
		if not np.any(mask):
			continue
		cat_idx = int(df['cat_idx'].values[np.where(mask)[0][0]])
		is_ctrl = bool(df['is_ctrl'].values[np.where(mask)[0][0]])
		color, marker = _group_style(cat_idx, is_ctrl)
		ax.scatter(x[mask], y[mask], color=color, marker=marker,
				   alpha=0.7, s=40)
	ax.set_xlabel(f'Translation Efficiency (controls at x={CONTROL_X:.1f})')
	ax.set_ylabel(y_col)
	ax.tick_params(axis='x', rotation=0)
	_add_legend(ax)
	plt.tight_layout()
	fname = 'trl_eff_vs_' + _safe_filename(y_col) + '.png'
	plt.savefig(os.path.join(out_dir, fname), dpi=300)
	plt.close('all')


def _plot_multi_y(df, y_cols, y_axis_label, out_dir):
	available = [c for c in y_cols if c in df.columns]
	if not available:
		return
	fig, ax = plt.subplots(figsize=(10, 8))
	x = df['trl_eff'].values
	masks = _group_masks(df)
	first_expression_label = _group_key(0, False)

	for i, y_col in enumerate(available):
		y = df[y_col].values
		color = MULTI_Y_COLORS[i % len(MULTI_Y_COLORS)]
		for label, mask in masks.items():
			if not np.any(mask):
				continue
			cat_idx = int(df['cat_idx'].values[np.where(mask)[0][0]])
			is_ctrl = bool(df['is_ctrl'].values[np.where(mask)[0][0]])
			_, marker = _group_style(cat_idx, is_ctrl)
			ax.scatter(
				x[mask], y[mask], color=color, marker=marker, alpha=0.5, s=40,
				label=(y_col if label == first_expression_label else None))

	ax.set_xlabel(f'Translation Efficiency (controls at x={CONTROL_X:.1f})')
	ax.set_ylabel(y_axis_label)
	handles, _ = ax.get_legend_handles_labels()
	handles.append(mlines.Line2D(
		[], [], color='gray', marker=EXPRESSION_MARKER, linestyle='None',
		markersize=6, label='Expression'))
	handles.append(mlines.Line2D(
		[], [], color='gray', marker=CONTROL_MARKER, linestyle='None',
		markersize=6, label='Control'))
	ax.legend(handles=handles, fontsize=7, loc='best')
	plt.tight_layout()
	fname = 'trl_eff_vs_' + _safe_filename(y_axis_label) + '.png'
	plt.savefig(os.path.join(out_dir, fname), dpi=300)
	plt.close('all')


def _plot_ratio(df, y_col, z_col, out_dir):
	if y_col not in df.columns or z_col not in df.columns:
		return
	fig, ax = plt.subplots(figsize=(10, 8))
	x = df['trl_eff'].values
	with np.errstate(divide='ignore', invalid='ignore'):
		y = df[y_col].values / df[z_col].values
	for label, mask in _group_masks(df).items():
		if not np.any(mask):
			continue
		cat_idx = int(df['cat_idx'].values[np.where(mask)[0][0]])
		is_ctrl = bool(df['is_ctrl'].values[np.where(mask)[0][0]])
		color, marker = _group_style(cat_idx, is_ctrl)
		ax.scatter(x[mask], y[mask], color=color, marker=marker,
				   alpha=0.7, s=40)
	y_label = y_col + ' divided by ' + z_col
	ax.set_xlabel(f'Translation Efficiency (controls at x={CONTROL_X:.1f})')
	ax.set_ylabel(y_label)
	_add_legend(ax)
	plt.tight_layout()
	fname = 'trl_eff_vs_' + _safe_filename(y_label) + '.png'
	plt.savefig(os.path.join(out_dir, fname), dpi=300)
	plt.close('all')


def _load_data(plot_out_dir):
	"""Merge the three gallery CSVs into a single dataframe keyed by variant.

	Adds trl_eff / is_ctrl / cat_idx columns derived from Variant Index.
	"""
	dt_csv = os.path.join(plot_out_dir, 'new_gene_counts_save_dt.csv')
	rnap_csv = os.path.join(
		plot_out_dir,
		'new_gene_counts_save_dt_ribosome_and_rnap_components.csv')
	ribo_csv = os.path.join(
		plot_out_dir,
		'new_gene_counts_save_dt_ribosome_components_by_subunit.csv')

	if not os.path.exists(dt_csv):
		return None

	df = pd.read_csv(dt_csv)

	if os.path.exists(rnap_csv):
		df_rnap = pd.read_csv(rnap_csv).drop(
			columns=['Variants'], errors='ignore')
		alpha = 'Avg Monomer Counts (RNAP): EG10893-MONOMER[c]'
		if alpha in df_rnap.columns:
			df_rnap[alpha] = df_rnap[alpha] / 2
		df = pd.concat([df, df_rnap], axis=1)

	if os.path.exists(ribo_csv):
		df_ribo = pd.read_csv(ribo_csv).drop(
			columns=['Variants'], errors='ignore')
		df = pd.concat([df, df_ribo], axis=1)

	df['trl_eff'] = df['Variant Index'].apply(_variant_to_x)
	df['is_ctrl'] = df['Variant Index'].apply(
		lambda i: is_control(int(i)))
	df['cat_idx'] = df['Variant Index'].apply(
		lambda i: variant_to_category(int(i))[0])
	return df


class Plot(variantAnalysisPlot.VariantAnalysisPlot):
	def do_plot(self, inputDir, plotOutDir, plotOutFileName, simDataFile,
				validationDataFile, metadata):

		df = _load_data(plotOutDir)
		if df is None:
			print(f'CSV not found in {plotOutDir} -- '
				  'run new_gene_counts_save_dt first.')
			return

		out_dir = os.path.join(plotOutDir, 'trl_eff_vs_all')
		os.makedirs(out_dir, exist_ok=True)

		n_single = 0
		for col in df.columns:
			if col in SKIP_COLS:
				continue
			if not np.issubdtype(df[col].dtype, np.number):
				continue
			_plot_single(df, col, out_dir)
			n_single += 1

		for label, y_cols in MULTI_Y_GROUPS:
			_plot_multi_y(df, y_cols, label, out_dir)

		for y_col, z_col in RATIO_PAIRS:
			_plot_ratio(df, y_col, z_col, out_dir)

		print(f'Wrote {n_single} single-y plots + '
			  f'{len(MULTI_Y_GROUPS)} multi-y + {len(RATIO_PAIRS)} ratio '
			  f'plots to {out_dir}')


if __name__ == '__main__':
	Plot().cli()
