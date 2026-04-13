"""
Scatterplot gallery for the new_gene_trl_eff_sweep variant comparing cell
properties with and without kcat constraints at matching translation
efficiencies.

Reads CSVs produced by new_gene_counts_save_dt, save_omes, and
save_transcriptome.  Points are colored by kcat status (blue = no kcat,
red = kcat) and shaped by control status (triangle = control/knockout,
circle = expression variant).

Variant index layout (44 indices total):
  No-kcat half (0-21): 0 = control, 1-21 = expression variants
  Kcat half (22-43): 22 = control, 23-43 = expression variants
"""

import os

import matplotlib.lines as mlines
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from models.ecoli.analysis import variantAnalysisPlot
from models.ecoli.sim.variants.new_gene_trl_eff_sweep import KCAT_HALF_START
from wholecell.analysis.analysis_tools import exportFigure

# -- Visual encoding constants ------------------------------------------------

POSTER_BLUE = (27/255, 132/255, 198/255)
POSTER_RED = (202/255, 0/255, 32/255)

CONTROL_INDICES = {0, KCAT_HALF_START}

OUTPUT_FOLDER_NAME = 'scatterplot_gallery'

# Colors used to distinguish multiple y-series in multi_y plots
MULTI_Y_COLORS = [
	(136/255, 205/255, 240/255),  # light blue
	(188/255, 140/255, 191/255),  # light purple
	(66/255, 170/255, 154/255),   # light green
	(221/255, 203/255, 119/255),  # gold
	(27/255, 132/255, 198/255),   # dark blue
]

EXTENDED_COLORS = [
	'#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0',
	'#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324',
	'#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075',
	'#808080', '#000000',
]


def _group_masks(variant_indices):
	"""Return boolean masks for the four point groups.

	Returns:
		dict mapping group label -> boolean array of same length as
		variant_indices.
	"""
	vi = np.asarray(variant_indices)
	is_kcat = vi >= KCAT_HALF_START
	is_control = np.isin(vi, list(CONTROL_INDICES))

	return {
		'No-kcat expression': ~is_kcat & ~is_control,
		'No-kcat control':    ~is_kcat &  is_control,
		'Kcat expression':     is_kcat & ~is_control,
		'Kcat control':        is_kcat &  is_control,
	}


def _group_style(label):
	"""Return (color, marker) for a group label."""
	styles = {
		'No-kcat expression': (POSTER_BLUE, 'o'),
		'No-kcat control':    (POSTER_BLUE, '^'),
		'Kcat expression':    (POSTER_RED,  'o'),
		'Kcat control':       (POSTER_RED,  '^'),
	}
	return styles[label]


def _add_legend(ax):
	"""Add a 4-entry legend to the axes."""
	handles = []
	for label in ['No-kcat expression', 'No-kcat control',
				  'Kcat expression', 'Kcat control']:
		color, marker = _group_style(label)
		handles.append(mlines.Line2D(
			[], [], color=color, marker=marker, linestyle='None',
			markersize=6, alpha=0.7, label=label))
	ax.legend(handles=handles, fontsize=7, loc='best')


# -- Plotting functions -------------------------------------------------------

def make_scatterplot(dt_data, variant_indices, x_col, y_col, output_folder,
					 output_to_label_dict):
	"""Single x vs y scatterplot with kcat/control visual encoding."""
	fig, ax = plt.subplots(figsize=(8, 8))
	x_data = dt_data[x_col].values
	y_data = dt_data[y_col].values
	masks = _group_masks(variant_indices)

	for label, mask in masks.items():
		if not np.any(mask):
			continue
		color, marker = _group_style(label)
		ax.scatter(x_data[mask], y_data[mask], color=color, marker=marker,
				   alpha=0.7, s=40)

	ax.set_xlabel(output_to_label_dict.get(x_col, x_col))
	ax.set_ylabel(output_to_label_dict.get(y_col, y_col))
	_add_legend(ax)

	filename = x_col.replace(' ', '_') + '_vs_' + y_col.replace(' ', '_')
	plt.savefig(os.path.join(output_folder, filename + '.png'), dpi=300)
	plt.close('all')


def make_scatterplot_special_y(dt_data, variant_indices, x_data, y_data,
							   x_label, y_label, output_folder):
	"""Scatterplot with pre-computed y data (e.g. ratio plots)."""
	fig, ax = plt.subplots(figsize=(8, 8))
	masks = _group_masks(variant_indices)

	for label, mask in masks.items():
		if not np.any(mask):
			continue
		color, marker = _group_style(label)
		ax.scatter(x_data[mask], y_data[mask], color=color, marker=marker,
				   alpha=0.7, s=40)

	ax.set_xlabel(x_label)
	ax.set_ylabel(y_label)
	_add_legend(ax)

	filename = (x_label.replace(' ', '_') + '_vs_'
				+ y_label.replace(' ', '_'))
	plt.savefig(os.path.join(output_folder, filename + '.png'), dpi=300)
	plt.close('all')


def make_scatterplot_multi_y(dt_data, variant_indices, x_col, y_cols,
							 y_axis_label, output_folder,
							 output_to_label_dict):
	"""Overlay multiple y-series, each with its own color, but all points
	still use kcat/control marker shapes."""
	fig, ax = plt.subplots(figsize=(8, 8))
	x_data = dt_data[x_col].values
	masks = _group_masks(variant_indices)

	for i, y_col in enumerate(y_cols):
		y_data = dt_data[y_col].values
		base_color = MULTI_Y_COLORS[i % len(MULTI_Y_COLORS)]

		for label, mask in masks.items():
			if not np.any(mask):
				continue
			_, marker = _group_style(label)
			ax.scatter(x_data[mask], y_data[mask], color=base_color,
					   marker=marker, alpha=0.5, s=40,
					   label=(output_to_label_dict.get(y_col, y_col)
							  if label == 'No-kcat expression' else None))

	ax.set_xlabel(output_to_label_dict.get(x_col, x_col))
	ax.set_ylabel(y_axis_label)
	# Only show series legend (not per-group), plus the shape legend
	handles_series, labels_series = ax.get_legend_handles_labels()
	# Add marker-shape legend entries
	handles_series.append(mlines.Line2D(
		[], [], color='gray', marker='o', linestyle='None',
		markersize=6, label='Expression'))
	handles_series.append(mlines.Line2D(
		[], [], color='gray', marker='^', linestyle='None',
		markersize=6, label='Control'))
	ax.legend(handles=handles_series, fontsize=7, loc='best')

	filename = (x_col.replace(' ', '_') + '_vs_'
				+ y_axis_label.replace(' ', '_'))
	plt.savefig(os.path.join(output_folder, filename + '.png'), dpi=300)
	plt.close('all')


# -- Comparisons dict ----------------------------------------------------------

COMPARISONS = {
	0: {'x': 'Avg New Gene Protein Counts', 'y': 'Doubling Time (min)'},
	1: {'x': 'Avg New Gene Protein Counts', 'y': 'Avg Active Ribosome Counts'},
	2: {'x': 'Avg New Gene Protein Counts', 'y': 'Avg Active RNA Polymerase Counts'},
	3: {'x': 'Avg New Gene Protein Counts', 'y': 'Avg Total Ribosome Counts'},
	4: {'x': 'Avg New Gene Protein Counts', 'y': 'Avg Total RNA Polymerase Counts'},
	5: {'x': 'Avg New Gene Protein Counts', 'y': 'Avg Total Ribosome Concentration (uM)'},
	6: {'x': 'Avg New Gene Protein Counts', 'y': 'Avg Total RNA Polymerase Concentration (uM)'},
	7: {'x': 'Avg New Gene Protein Counts', 'y': 'Avg Dry Mass (fg)'},
	8: {'x': 'Avg New Gene Protein Counts', 'y': 'Avg Cell Mass (fg)'},
	9: {'x': 'Avg New Gene Protein Counts', 'y': 'Avg mRNA Mass (fg)'},
	10: {'x': 'Avg New Gene Protein Counts', 'y': 'Avg Protein Mass (fg)'},
	11: {'x': 'Avg New Gene Protein Counts', 'y': 'Avg DNA Mass (fg)'},
	12: {'x': 'Avg New Gene Protein Counts', 'y': 'Avg rRNA Mass (fg)'},
	13: {'x': 'Avg New Gene Protein Counts', 'y': 'Avg tRNA Mass (fg)'},
	14: {'x': 'Avg New Gene Protein Counts', 'y': 'Avg Water Mass (fg)'},
	15: {'x': 'Avg New Gene Protein Counts', 'y': 'Avg Small Molecule Mass (fg)'},
	16: {'x': 'Avg New Gene Protein Counts', 'y': 'Avg Membrane Mass (fg)'},
	17: {'x': 'Avg New Gene Protein Counts', 'y': 'Percent Completion (fraction of seeds that reached all gens)'},
	18: {'x': 'Avg Active Ribosome Counts', 'y': 'Avg Active RNA Polymerase Counts'},
	19: {'x': 'Avg Active Ribosome Counts', 'y': 'Avg Total Ribosome Counts'},
	20: {'x': 'Avg Active RNA Polymerase Counts', 'y': 'Avg Total RNA Polymerase Counts'},
	21: {'x': 'Avg Total Ribosome Counts', 'y': 'Avg Total RNA Polymerase Counts'},
	22: {'x': 'Avg Total Ribosome Counts', 'y': 'Avg Total Ribosome Concentration (uM)'},
	23: {'x': 'Avg Total RNA Polymerase Counts', 'y': 'Avg Total RNA Polymerase Concentration (uM)'},
	24: {'x': 'Avg Total Ribosome Concentration (uM)', 'y': 'Avg Total RNA Polymerase Concentration (uM)'},
	25: {'x': 'Avg Active Ribosome Counts', 'y': 'Doubling Time (min)'},
	26: {'x': 'Avg Active RNA Polymerase Counts', 'y': 'Doubling Time (min)'},
	27: {'x': 'Avg Total Ribosome Counts', 'y': 'Doubling Time (min)'},
	28: {'x': 'Avg Total RNA Polymerase Counts', 'y': 'Doubling Time (min)'},
	29: {'x': 'Avg Total Ribosome Concentration (uM)', 'y': 'Doubling Time (min)'},
	30: {'x': 'Avg Total RNA Polymerase Concentration (uM)', 'y': 'Doubling Time (min)'},
	31: {'x': 'Avg New Gene Protein Counts', 'y': 'Avg Active Ribosome Counts', 'z': 'Avg Total Ribosome Counts'},
	32: {'x': 'Avg New Gene Protein Counts', 'y': 'Avg Active RNA Polymerase Counts', 'z': 'Avg Total RNA Polymerase Counts'},
	33: {'x': 'Avg New Gene Protein Counts', 'y': 'Avg Total Ribosome Counts', 'z': 'Avg Total RNA Polymerase Counts'},
	34: {'x': 'Avg New Gene Protein Counts', 'y': 'Avg Total Ribosome Concentration (uM)', 'z': 'Avg Total RNA Polymerase Concentration (uM)'},
	35: {'x': 'Avg Cell Mass (fg)', 'y': 'Avg mRNA Mass (fg)'},
	36: {'x': 'Avg Cell Mass (fg)', 'y': 'Avg Protein Mass (fg)'},
	37: {'x': 'Avg Cell Mass (fg)', 'y': 'Avg DNA Mass (fg)'},
	38: {'x': 'Avg Cell Mass (fg)', 'y': 'Avg rRNA Mass (fg)'},
	39: {'x': 'Avg Cell Mass (fg)', 'y': 'Avg tRNA Mass (fg)'},
	40: {'x': 'Avg Cell Mass (fg)', 'y': 'Avg Water Mass (fg)'},
	41: {'x': 'Avg Cell Mass (fg)', 'y': 'Avg Small Molecule Mass (fg)'},
	42: {'x': 'Avg Cell Mass (fg)', 'y': 'Avg Membrane Mass (fg)'},
	43: {'x': 'Avg Total Ribosome Counts', 'y': 'Percent Completion (fraction of seeds that reached all gens)'},
	44: {'x': 'Avg Total RNA Polymerase Counts', 'y': 'Percent Completion (fraction of seeds that reached all gens)'},
	45: {'x': 'Percent Completion (fraction of seeds that reached all gens)', 'y': 'Doubling Time (min)'},
	46: {
		'x': 'Avg New Gene Protein Counts',
		'y_list': ['Avg Total Ribosome Counts',
			'Avg rRNA Counts: RRSA-RRNA[c]',
			'Avg rRNA Counts: RRLA-RRNA[c]',
			'Avg rRNA Counts: RRFA-RRNA[c]'],
		'y_axis_label': 'rRNA Counts',
	},
	47: {
		'x': 'Avg New Gene Protein Counts',
		'y_list': [
			'Avg Total RNA Polymerase Counts',
			'Avg Monomer Counts (RNAP): EG10893-MONOMER[c]',
			'Avg Monomer Counts (RNAP): RPOC-MONOMER[c]',
			'Avg Monomer Counts (RNAP): RPOB-MONOMER[c]'],
		'y_axis_label': 'RNA Polymerase Counts',
	},
	48: {
		'x': 'Avg New Gene Protein Counts',
		'y_list': [
			'Avg Total Ribosome Counts',
			'Avg Total Counts (Ribosomal 30s)',
			'Avg Total Counts (Ribosomal 50s)'],
		'y_axis_label': 'Ribosomal Subunit Counts',
	},
	49: {
		'x': 'Avg New Gene Protein Counts',
		'y_list': [
			'Avg Total Counts (Ribosomal 30s)',
			'Avg Limiting Protein Counts (Ribosomal 30s)',
			'Avg 16s rRNA Counts (Ribosomal 30s)'],
		'y_axis_label': 'Limiting Counts Ribosomal 30s',
	},
	50: {
		'x': 'Avg New Gene Protein Counts',
		'y_list': [
			'Avg Total Counts (Ribosomal 50s)',
			'Avg Limiting Protein Counts (Ribosomal 50s)',
			'Avg 23s rRNA Counts (Ribosomal 50s)',
			'Avg 5s rRNA Counts (Ribosomal 50s)'],
		'y_axis_label': 'Limiting Counts Ribosomal 50s',
	},
	51: {'x': 'Variant Index', 'y': 'Percent Completion (fraction of seeds that reached all gens)'},
	52: {
		'x': 'Doubling Time (min)',
		'y_list': ['Avg Total Ribosome Counts',
			'Avg rRNA Counts: RRSA-RRNA[c]',
			'Avg rRNA Counts: RRLA-RRNA[c]',
			'Avg rRNA Counts: RRFA-RRNA[c]'],
		'y_axis_label': 'rRNA Counts',
	},
	53: {
		'x': 'Doubling Time (min)',
		'y_list': [
			'Avg Total RNA Polymerase Counts',
			'Avg Monomer Counts (RNAP): EG10893-MONOMER[c]',
			'Avg Monomer Counts (RNAP): RPOC-MONOMER[c]',
			'Avg Monomer Counts (RNAP): RPOB-MONOMER[c]'],
		'y_axis_label': 'RNA Polymerase Counts',
	},
	54: {
		'x': 'Doubling Time (min)',
		'y_list': [
			'Avg Total Ribosome Counts',
			'Avg Total Counts (Ribosomal 30s)',
			'Avg Total Counts (Ribosomal 50s)'],
		'y_axis_label': 'Ribosomal Subunit Counts',
	},
	55: {
		'x': 'Doubling Time (min)',
		'y_list': [
			'Avg Total Counts (Ribosomal 30s)',
			'Avg Limiting Protein Counts (Ribosomal 30s)',
			'Avg 16s rRNA Counts (Ribosomal 30s)'],
		'y_axis_label': 'Limiting Counts Ribosomal 30s',
	},
	56: {
		'x': 'Doubling Time (min)',
		'y_list': [
			'Avg Total Counts (Ribosomal 50s)',
			'Avg Limiting Protein Counts (Ribosomal 50s)',
			'Avg 23s rRNA Counts (Ribosomal 50s)',
			'Avg 5s rRNA Counts (Ribosomal 50s)'],
		'y_axis_label': 'Limiting Counts Ribosomal 50s',
	},
	57: {
		'x': 'Doubling Time (min)',
		'y': 'Avg Total Counts (Ribosomal 30s)',
		'z': 'Avg Total Counts (Ribosomal 50s)',
	},
	58: {'x': 'Average New Gene Proteome Mass Fraction', 'y': 'Percent Completion (fraction of seeds that reached all gens)'},
	59: {'x': 'Average New Gene Proteome Mass Fraction', 'y': 'Doubling Time (min)'},
	60: {'x': 'Average New Gene Proteome Mass Fraction', 'y': 'Avg Total Ribosome Counts'},
	61: {'x': 'Average New Gene Proteome Mass Fraction', 'y': 'Avg Total RNA Polymerase Counts'},
	62: {'x': 'Average New Gene Proteome Mass Fraction', 'y': 'Avg Cell Mass (fg)'},
}


class Plot(variantAnalysisPlot.VariantAnalysisPlot):
	def do_plot(self, inputDir, plotOutDir, plotOutFileName, simDataFile,
				validationDataFile, metadata):

		# -- Load CSVs ---------------------------------------------------------
		dt_filename = os.path.join(
			plotOutDir, 'new_gene_counts_save_dt.csv')
		dt_rnap_filename = os.path.join(
			plotOutDir,
			'new_gene_counts_save_dt_ribosome_and_rnap_components.csv')
		dt_ribosome_filename = os.path.join(
			plotOutDir,
			'new_gene_counts_save_dt_ribosome_components_by_subunit.csv')

		if not os.path.exists(dt_filename):
			print(f'CSV not found: {dt_filename}  -- '
				  'run new_gene_counts_save_dt first.')
			return

		dt_data = pd.read_csv(dt_filename)

		if os.path.exists(dt_rnap_filename):
			dt_data_rnap = pd.read_csv(dt_rnap_filename)
			dt_data_rnap = dt_data_rnap.drop(columns=['Variants'],
											 errors='ignore')
			# Divide RNAP alpha subunit by 2 to account for stoichiometry
			alpha_col = 'Avg Monomer Counts (RNAP): EG10893-MONOMER[c]'
			if alpha_col in dt_data_rnap.columns:
				dt_data_rnap[alpha_col] = dt_data_rnap[alpha_col] / 2
			dt_data = pd.concat([dt_data, dt_data_rnap], axis=1)

		dt_data_ribosome = None
		if os.path.exists(dt_ribosome_filename):
			dt_data_ribosome = pd.read_csv(dt_ribosome_filename)
			dt_data_ribosome = dt_data_ribosome.drop(
				columns=['Variants'], errors='ignore')
			dt_data = pd.concat([dt_data, dt_data_ribosome], axis=1)

		variant_indices = dt_data['Variant Index'].values

		# -- Build label dict --------------------------------------------------
		output_to_label_dict = {
			'Avg New Gene Protein Counts': 'GFP Counts',
		}
		for comp in COMPARISONS.values():
			for key in ('x', 'y', 'z'):
				val = comp.get(key)
				if val and val not in output_to_label_dict:
					output_to_label_dict[val] = val
			for y_col in comp.get('y_list', []):
				if y_col not in output_to_label_dict:
					output_to_label_dict[y_col] = y_col

		# -- Create output folder ----------------------------------------------
		output_folder = os.path.join(plotOutDir, OUTPUT_FOLDER_NAME)
		os.makedirs(output_folder, exist_ok=True)

		# -- Generate scatterplots ---------------------------------------------
		for comp in COMPARISONS.values():
			if 'y_list' in comp:
				make_scatterplot_multi_y(
					dt_data, variant_indices, comp['x'],
					comp['y_list'], comp['y_axis_label'],
					output_folder, output_to_label_dict)
			elif 'z' not in comp:
				make_scatterplot(
					dt_data, variant_indices, comp['x'], comp['y'],
					output_folder, output_to_label_dict)
			else:
				x_data = dt_data[comp['x']].values
				y_data = (dt_data[comp['y']].values
						  / dt_data[comp['z']].values)
				x_label = output_to_label_dict[comp['x']]
				y_label = (output_to_label_dict[comp['y']]
						   + ' Divided by '
						   + output_to_label_dict[comp['z']])
				make_scatterplot_special_y(
					dt_data, variant_indices, x_data, y_data,
					x_label, y_label, output_folder)

		# -- Ribosomal component heatmaps --------------------------------------
		if dt_data_ribosome is not None:
			self._plot_ribosome_heatmaps(
				dt_data, dt_data_ribosome, variant_indices,
				output_folder, output_to_label_dict)

		print(f'Scatterplot gallery saved to {output_folder}')

	def _plot_ribosome_heatmaps(self, dt_data, dt_data_ribosome,
								variant_indices, output_folder,
								output_to_label_dict):
		"""Generate ribosomal component scatter and heatmap plots."""
		ribosome_cols = list(dt_data_ribosome.columns)
		for col in ribosome_cols:
			if col not in output_to_label_dict:
				output_to_label_dict[col] = col

		s30_cols = [c for c in ribosome_cols if '30s' in c]
		s50_cols = [c for c in ribosome_cols if '50s' in c]

		sorted_variant_indices = np.argsort(
			dt_data['Doubling Time (min)'].values)
		x_labels = [str(variant_indices[i]) for i in sorted_variant_indices]

		# -- 30s / 50s component scatterplots ----------------------------------
		special_30s = [
			'Avg Total Counts (Ribosomal 30s)',
			'Avg Limiting Protein Counts (Ribosomal 30s)',
			'Avg 16s rRNA Counts (Ribosomal 30s)',
		]
		s30_protein_cols = [c for c in s30_cols if c not in special_30s]

		special_50s = [
			'Avg Total Counts (Ribosomal 50s)',
			'Avg Limiting Protein Counts (Ribosomal 50s)',
			'Avg 23s rRNA Counts (Ribosomal 50s)',
			'Avg 5s rRNA Counts (Ribosomal 50s)',
		]
		s50_protein_cols = [c for c in s50_cols if c not in special_50s]

		# 30s scatter
		plt.figure(figsize=(10, 10))
		x_data = x_labels
		for i, col in enumerate(s30_protein_cols):
			y_data = dt_data_ribosome[col].values[sorted_variant_indices]
			plt.scatter(x_data, y_data,
						color=EXTENDED_COLORS[i % len(EXTENDED_COLORS)],
						alpha=0.5, label=output_to_label_dict[col])
		for col, marker in [
				('Avg Limiting Protein Counts (Ribosomal 30s)', 'v'),
				('Avg 16s rRNA Counts (Ribosomal 30s)', '*'),
				('Avg Total Counts (Ribosomal 30s)', 'x')]:
			if col in dt_data_ribosome.columns:
				plt.scatter(
					x_data,
					dt_data_ribosome[col].values[sorted_variant_indices],
					color='#e6194b', alpha=0.5, label=col, marker=marker)
		plt.xlabel('Variant Index')
		plt.ylabel('Counts')
		plt.savefig(os.path.join(
			output_folder, 'Variant_Index_vs_30s_Components.png'), dpi=300)
		plt.close('all')

		# 50s scatter
		plt.figure(figsize=(10, 10))
		for i, col in enumerate(s50_protein_cols):
			y_data = dt_data_ribosome[col].values[sorted_variant_indices]
			plt.scatter(x_data, y_data,
						color=EXTENDED_COLORS[i % len(EXTENDED_COLORS)],
						alpha=0.5, label=output_to_label_dict[col])
		for col, marker in [
				('Avg Limiting Protein Counts (Ribosomal 50s)', 'v'),
				('Avg 23s rRNA Counts (Ribosomal 50s)', '*'),
				('Avg 5s rRNA Counts (Ribosomal 50s)', 'd'),
				('Avg Total Counts (Ribosomal 50s)', 'x')]:
			if col in dt_data_ribosome.columns:
				plt.scatter(
					x_data,
					dt_data_ribosome[col].values[sorted_variant_indices],
					color='#e6194b', alpha=0.5, label=col, marker=marker)
		plt.xlabel('Variant Index')
		plt.ylabel('Counts')
		plt.savefig(os.path.join(
			output_folder, 'Variant_Index_vs_50s_Components.png'), dpi=300)
		plt.close('all')

		# -- Rank heatmaps -----------------------------------------------------
		# 30s rank heatmap
		s30_cols_plus_rRNA = (['Avg 16s rRNA Counts (Ribosomal 30s)']
							  + s30_protein_cols)
		labels_30 = ['16s rRNA'] + [
			c.replace('Avg Monomer Counts (Ribosomal 30s): ', '')
			for c in s30_protein_cols]
		if all(c in dt_data_ribosome.columns for c in s30_cols_plus_rRNA):
			ranks = dt_data_ribosome[s30_cols_plus_rRNA].rank(
				axis=1, ascending=False)
			plt.figure(figsize=(16, 8))
			sns.heatmap(
				ranks.iloc[sorted_variant_indices].T,
				cmap='viridis', xticklabels=x_labels,
				yticklabels=labels_30)
			plt.xlabel('Variant Index, Ordered by Doubling Time '
					   '(Rank 1 = Most Abundant)')
			plt.ylabel('Protein or rRNA')
			for i, vi in enumerate(sorted_variant_indices):
				dt = dt_data.loc[vi, 'Doubling Time (min)']
				plt.text(i + 0.5, -0.5, f'{dt:.2f}',
						 ha='center', va='center', fontsize=8)
			plt.tight_layout()
			plt.savefig(os.path.join(
				output_folder,
				'Variant_Index_vs_30s_Components_heatmap_ranked.png'))
			plt.close('all')

		# 50s rank heatmap
		s50_cols_plus_rRNA = ([
			'Avg 5s rRNA Counts (Ribosomal 50s)',
			'Avg 23s rRNA Counts (Ribosomal 50s)'] + s50_protein_cols)
		labels_50 = ['5s rRNA', '23s rRNA'] + [
			c.replace('Avg Monomer Counts (Ribosomal 50s): ', '')
			for c in s50_protein_cols]
		if all(c in dt_data_ribosome.columns for c in s50_cols_plus_rRNA):
			ranks = dt_data_ribosome[s50_cols_plus_rRNA].rank(
				axis=1, ascending=False)
			plt.figure(figsize=(16, 8))
			sns.heatmap(
				ranks.iloc[sorted_variant_indices].T,
				cmap='viridis', xticklabels=x_labels,
				yticklabels=labels_50)
			plt.xlabel('Variant Index, Ordered by Doubling Time '
					   '(Rank 1 = Most Abundant)')
			plt.ylabel('Protein or rRNA')
			for i, vi in enumerate(sorted_variant_indices):
				dt = dt_data.loc[vi, 'Doubling Time (min)']
				plt.text(i + 0.5, -0.5, f'{dt:.2f}',
						 ha='center', va='center', fontsize=8)
			plt.tight_layout()
			plt.savefig(os.path.join(
				output_folder,
				'Variant_Index_vs_50s_Components_heatmap_ranked.png'))
			plt.close('all')

		# -- Limiting-protein rank heatmaps ------------------------------------
		# 30s limiting
		lim_30_cols = [
			'Avg Limiting Protein Counts (Ribosomal 30s)',
			'Avg 16s rRNA Counts (Ribosomal 30s)']
		if all(c in dt_data_ribosome.columns for c in lim_30_cols):
			ranks = dt_data_ribosome[lim_30_cols].rank(
				axis=1, ascending=False)
			plt.figure(figsize=(16, 4))
			sns.heatmap(
				ranks.iloc[sorted_variant_indices].T,
				cmap='viridis', xticklabels=x_labels,
				yticklabels=['Lim.Prot.', '16s rRNA'])
			plt.xlabel('Variant Index, Ordered by Doubling Time '
					   '(Rank 1 = Most Abundant)')
			for i, vi in enumerate(sorted_variant_indices):
				dt = dt_data.loc[vi, 'Doubling Time (min)']
				plt.text(i + 0.5, -0.5, f'{dt:.2f}',
						 ha='center', va='center', fontsize=8)
			plt.tight_layout()
			plt.savefig(os.path.join(
				output_folder,
				'Variant_Index_vs_30s_Components_heatmap_ranked_limiting.png'))
			plt.close('all')

		# 50s limiting
		lim_50_cols = [
			'Avg Limiting Protein Counts (Ribosomal 50s)',
			'Avg 5s rRNA Counts (Ribosomal 50s)',
			'Avg 23s rRNA Counts (Ribosomal 50s)']
		if all(c in dt_data_ribosome.columns for c in lim_50_cols):
			ranks = dt_data_ribosome[lim_50_cols].rank(
				axis=1, ascending=False)
			plt.figure(figsize=(16, 4))
			sns.heatmap(
				ranks.iloc[sorted_variant_indices].T,
				cmap='viridis', xticklabels=x_labels,
				yticklabels=['Lim.Prot.', '5s rRNA', '23s rRNA'])
			plt.xlabel('Variant Index, Ordered by Doubling Time '
					   '(Rank 1 = Most Abundant)')
			for i, vi in enumerate(sorted_variant_indices):
				dt = dt_data.loc[vi, 'Doubling Time (min)']
				plt.text(i + 0.5, -0.5, f'{dt:.2f}',
						 ha='center', va='center', fontsize=8)
			plt.tight_layout()
			plt.savefig(os.path.join(
				output_folder,
				'Variant_Index_vs_50s_Components_heatmap_ranked_limiting.png'))
			plt.close('all')

		# 30s vs 50s total rank heatmap
		total_cols = [
			'Avg Total Counts (Ribosomal 30s)',
			'Avg Total Counts (Ribosomal 50s)']
		if all(c in dt_data_ribosome.columns for c in total_cols):
			ranks = dt_data_ribosome[total_cols].rank(
				axis=1, ascending=False)
			plt.figure(figsize=(16, 3))
			sns.heatmap(
				ranks.iloc[sorted_variant_indices].T,
				cmap='viridis', xticklabels=x_labels,
				yticklabels=['30s', '50s'])
			plt.xlabel('Variant Index, Ordered by Doubling Time '
					   '(Rank 1 = Most Abundant)')
			plt.ylabel('Ribosomal Subunit')
			for i, vi in enumerate(sorted_variant_indices):
				dt = dt_data.loc[vi, 'Doubling Time (min)']
				plt.text(i + 0.5, -0.5, f'{dt:.2f}',
						 ha='center', va='center', fontsize=8)
			plt.tight_layout()
			plt.savefig(os.path.join(
				output_folder,
				'Variant_Index_vs_30s_50s_Components_heatmap_ranked_total.png'))
			plt.close('all')

		# -- Ratio heatmaps ----------------------------------------------------
		# 30s ratio
		if (s30_protein_cols
				and 'Avg Limiting Protein Counts (Ribosomal 30s)'
				in dt_data_ribosome.columns):
			protein_counts = dt_data_ribosome[s30_protein_cols].iloc[
				sorted_variant_indices]
			limiting = dt_data_ribosome[
				'Avg Limiting Protein Counts (Ribosomal 30s)'].iloc[
				sorted_variant_indices]
			ratios = protein_counts.div(limiting, axis=0)
			protein_labels = [
				c.replace('Avg Monomer Counts (Ribosomal 30s): ', '')
				for c in s30_protein_cols]
			plt.figure(figsize=(14, 6))
			sns.heatmap(
				ratios.T, cmap='magma', annot=False,
				cbar_kws={'label': 'Ratio to Limiting Protein'},
				xticklabels=x_labels, yticklabels=protein_labels, vmin=1.0)
			plt.xlabel('Variant Index, 30s Ratio of Protein Counts to '
					   'Limiting Protein Count (per Variant)')
			plt.ylabel('Protein')
			for i, vi in enumerate(sorted_variant_indices):
				dt = dt_data.loc[vi, 'Doubling Time (min)']
				plt.text(i + 0.5, -0.5, f'{dt:.2f}',
						 ha='center', va='center', fontsize=8)
			plt.tight_layout()
			plt.savefig(os.path.join(
				output_folder,
				'Variant_Index_vs_30s_Components_ratio_heatmap.png'))
			plt.close('all')

		# 50s ratio
		if (s50_protein_cols
				and 'Avg Limiting Protein Counts (Ribosomal 50s)'
				in dt_data_ribosome.columns):
			protein_counts = dt_data_ribosome[s50_protein_cols].iloc[
				sorted_variant_indices]
			limiting = dt_data_ribosome[
				'Avg Limiting Protein Counts (Ribosomal 50s)'].iloc[
				sorted_variant_indices]
			ratios = protein_counts.div(limiting, axis=0)
			protein_labels = [
				c.replace('Avg Monomer Counts (Ribosomal 50s): ', '')
				for c in s50_protein_cols]
			plt.figure(figsize=(14, 6))
			sns.heatmap(
				ratios.T, cmap='magma', annot=False,
				cbar_kws={'label': 'Ratio to Limiting Protein'},
				xticklabels=x_labels, yticklabels=protein_labels, vmin=1.0)
			plt.xlabel('Variant Index, 50s Ratio of Protein Counts to '
					   'Limiting Protein Count (per Variant)')
			plt.ylabel('Protein')
			for i, vi in enumerate(sorted_variant_indices):
				dt = dt_data.loc[vi, 'Doubling Time (min)']
				plt.text(i + 0.5, -0.5, f'{dt:.2f}',
						 ha='center', va='center', fontsize=8)
			plt.tight_layout()
			plt.savefig(os.path.join(
				output_folder,
				'Variant_Index_vs_50s_Components_ratio_heatmap.png'))
			plt.close('all')

		# -- Excess heatmaps ---------------------------------------------------
		# 30s excess
		if (s30_protein_cols
				and 'Avg Limiting Protein Counts (Ribosomal 30s)'
				in dt_data_ribosome.columns):
			protein_counts = dt_data_ribosome[s30_protein_cols].iloc[
				sorted_variant_indices]
			limiting = dt_data_ribosome[
				'Avg Limiting Protein Counts (Ribosomal 30s)'].iloc[
				sorted_variant_indices]
			differences = protein_counts.subtract(limiting, axis=0)
			protein_labels = [
				c.replace('Avg Monomer Counts (Ribosomal 30s): ', '')
				for c in s30_protein_cols]
			plt.figure(figsize=(14, 6))
			sns.heatmap(
				differences.T, cmap='coolwarm', center=0, annot=False,
				xticklabels=x_labels, yticklabels=protein_labels,
				cbar_kws={
					'label': 'Protein Count \u2013 Limiting Protein Count'})
			plt.xlabel('Variant Index, 30s Excess Protein Counts Relative '
					   'to Limiting Protein (per Variant)')
			plt.ylabel('Protein')
			for i, vi in enumerate(sorted_variant_indices):
				dt = dt_data.loc[vi, 'Doubling Time (min)']
				plt.text(i + 0.5, -0.5, f'{dt:.2f}',
						 ha='center', va='center', fontsize=8)
			plt.tight_layout()
			plt.savefig(os.path.join(
				output_folder,
				'Variant_Index_vs_30s_Components_excess_heatmap.png'))
			plt.close('all')

		# 50s excess
		if (s50_protein_cols
				and 'Avg Limiting Protein Counts (Ribosomal 50s)'
				in dt_data_ribosome.columns):
			protein_counts = dt_data_ribosome[s50_protein_cols].iloc[
				sorted_variant_indices]
			limiting = dt_data_ribosome[
				'Avg Limiting Protein Counts (Ribosomal 50s)'].iloc[
				sorted_variant_indices]
			differences = protein_counts.subtract(limiting, axis=0)
			protein_labels = [
				c.replace('Avg Monomer Counts (Ribosomal 50s): ', '')
				for c in s50_protein_cols]
			plt.figure(figsize=(14, 6))
			sns.heatmap(
				differences.T, cmap='coolwarm', center=0, annot=False,
				xticklabels=x_labels, yticklabels=protein_labels,
				cbar_kws={
					'label': 'Protein Count \u2013 Limiting Protein Count'})
			plt.xlabel('Variant Index, 50s Excess Protein Counts Relative '
					   'to Limiting Protein (per Variant)')
			plt.ylabel('Protein')
			for i, vi in enumerate(sorted_variant_indices):
				dt = dt_data.loc[vi, 'Doubling Time (min)']
				plt.text(i + 0.5, -0.5, f'{dt:.2f}',
						 ha='center', va='center', fontsize=8)
			plt.tight_layout()
			plt.savefig(os.path.join(
				output_folder,
				'Variant_Index_vs_50s_Components_excess_heatmap.png'))
			plt.close('all')


if __name__ == '__main__':
	Plot().cli()
