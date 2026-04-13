"""
Enzyme count and concentration reduction under kcat constraints for the
new_gene_trl_eff_sweep variant.

For each non-control variant, compares the mean enzyme counts and
concentrations (for enzymes with kcat > 1 in the max kcat estimates) against
the matched control variant (variant 0 for no-kcat, variant 22 for kcat).

Outputs:
  - Detail CSV: per-(variant, enzyme) counts, concentrations, and ratios
  - Scatter plot: doubling time vs mean enzyme count ratio
  - Scatter plot: translation efficiency vs mean enzyme count ratio
  - Heatmap: per-enzyme count ratios across kcat variants
"""

import csv
import os
import pickle

import numpy as np
from matplotlib import pyplot as plt
from matplotlib.colors import TwoSlopeNorm

from models.ecoli.analysis import variantAnalysisPlot
from models.ecoli.sim.variants.new_gene_trl_eff_sweep import (
	TRL_EFF_VALUES, N_TRL_EFF, KCAT_HALF_START)
from wholecell.analysis.analysis_tools import exportFigure, read_stacked_columns
from wholecell.io.tablereader import TableReader

# Use last 8 generations
N_GENS_TO_USE = 8


class Plot(variantAnalysisPlot.VariantAnalysisPlot):
	def do_plot(self, inputDir, plotOutDir, plotOutFileName, simDataFile,
				validationDataFile, metadata):

		variant_indexes = self.ap.get_variants()
		available = set(variant_indexes)
		n_total_gens = self.ap.n_generation
		ignore_first_n_gens = max(n_total_gens - N_GENS_TO_USE, 0)
		gen_range = np.arange(ignore_first_n_gens, n_total_gens)

		# -------------------------------------------------------------- #
		# Load kcat estimates and build enzyme index maps                  #
		# -------------------------------------------------------------- #
		# Find a kcat variant to load sim_data from
		kcat_variant = None
		for vi in variant_indexes:
			if vi >= KCAT_HALF_START and (vi - KCAT_HALF_START) != 0:
				kcat_variant = vi
				break
		if kcat_variant is None:
			print('Skipping: no non-control kcat variants found.')
			return

		kb_path = self.ap.get_variant_kb(kcat_variant)
		with open(kb_path, 'rb') as f:
			sim_data = pickle.load(f)

		metabolism = sim_data.process.metabolism
		if not hasattr(metabolism, 'kcat_estimates'):
			print('Skipping: sim_data does not have kcat_estimates.')
			return

		# Filter for kcat > 1
		all_kcat_estimates = metabolism.kcat_estimates['max']
		filtered_cat_ids = sorted(set(
			cat_id for (rxn_id, cat_id), kcat in all_kcat_estimates.items()
			if kcat > 1))

		if not filtered_cat_ids:
			print('Skipping: no enzymes with kcat > 1.')
			return

		# Get listener catalyst_ids from a representative cell
		ref_cells = self.ap.get_cells(
			variant=[kcat_variant], generation=[0], only_successful=True)
		if not ref_cells:
			print('Skipping: no cells for reference variant.')
			return

		fba_reader = TableReader(
			os.path.join(ref_cells[0], 'simOut', 'FBAResults'))
		listener_cat_ids = fba_reader.readAttribute('catalyst_ids')
		cat_id_to_idx = {c: i for i, c in enumerate(listener_cat_ids)}

		# Keep only enzymes present in the listener
		target_cat_ids = [c for c in filtered_cat_ids if c in cat_id_to_idx]
		target_cat_indices = np.array([cat_id_to_idx[c] for c in target_cat_ids])
		n_enzymes = len(target_cat_ids)
		print(f'{n_enzymes} enzymes with kcat > 1 found in listener.')

		# -------------------------------------------------------------- #
		# Helper: compute mean dt, counts, and concentrations for a variant#
		# -------------------------------------------------------------- #
		def get_variant_data(vi):
			"""Return (mean_dt, mean_counts, mean_conc) or Nones."""
			cells = self.ap.get_cells(
				variant=[vi], generation=gen_range, only_successful=True)
			if len(cells) == 0:
				return np.nan, None, None

			# Doubling times
			dt_arr = read_stacked_columns(
				cells, 'Main', 'time',
				fun=lambda x: (x[-1] - x[0]) / 60.)
			mean_dt = np.nanmean(dt_arr)

			# Enzyme counts and concentrations per cell
			cell_counts = []
			cell_concs = []
			for cell_path in cells:
				sim_out = os.path.join(cell_path, 'simOut')
				try:
					cat_counts = TableReader(
						os.path.join(sim_out, 'FBAResults')
					).readColumn('catalyst_counts',
								 indices=target_cat_indices,
								 squeeze=False)[1:]

					counts_to_molar = TableReader(
						os.path.join(sim_out, 'EnzymeKinetics')
					).readColumn('countsToMolar', squeeze=True)[1:]

					# Mean over timesteps
					cell_counts.append(np.mean(cat_counts, axis=0))
					# Concentration: counts * countsToMolar per timestep, then mean
					conc = cat_counts * counts_to_molar[:, np.newaxis]
					cell_concs.append(np.mean(conc, axis=0))
				except Exception as e:
					print(f'  Ignored exception reading {sim_out}: {e!r}')
					continue

			if not cell_counts:
				return mean_dt, None, None

			return (
				mean_dt,
				np.mean(cell_counts, axis=0),
				np.mean(cell_concs, axis=0),
			)

		# -------------------------------------------------------------- #
		# Compute control data                                             #
		# -------------------------------------------------------------- #
		print('Computing control data...')
		ctrl_no_kcat_vi = 0
		ctrl_kcat_vi = KCAT_HALF_START  # 22

		ctrl_dt = {}
		ctrl_counts = {}
		ctrl_conc = {}
		for vi in [ctrl_no_kcat_vi, ctrl_kcat_vi]:
			if vi not in available:
				print(f'  Control variant {vi} not available.')
				continue
			dt, counts, conc = get_variant_data(vi)
			ctrl_dt[vi] = dt
			ctrl_counts[vi] = counts
			ctrl_conc[vi] = conc
			print(f'  Variant {vi}: dt={dt:.1f} min')

		# -------------------------------------------------------------- #
		# Compute data for all non-control variants                        #
		# -------------------------------------------------------------- #
		print('Computing variant data...')
		results = []
		for vi in variant_indexes:
			is_kcat = (vi >= KCAT_HALF_START)
			local_idx = vi - KCAT_HALF_START if is_kcat else vi
			if local_idx == 0:
				continue  # skip controls

			ctrl_vi = ctrl_kcat_vi if is_kcat else ctrl_no_kcat_vi
			if ctrl_counts.get(ctrl_vi) is None:
				continue

			trl_eff = TRL_EFF_VALUES[local_idx - 1]
			print(f'  Variant {vi} (trl_eff={trl_eff}, kcat={is_kcat})...')

			dt, counts, conc = get_variant_data(vi)
			if counts is None:
				continue

			with np.errstate(divide='ignore', invalid='ignore'):
				count_ratio = counts / ctrl_counts[ctrl_vi]
				conc_ratio = conc / ctrl_conc[ctrl_vi]
				count_ratio[~np.isfinite(count_ratio)] = np.nan
				conc_ratio[~np.isfinite(conc_ratio)] = np.nan

			results.append({
				'vi': vi,
				'is_kcat': is_kcat,
				'local_idx': local_idx,
				'trl_eff': trl_eff,
				'dt': dt,
				'ctrl_dt': ctrl_dt[ctrl_vi],
				'counts': counts,
				'conc': conc,
				'ctrl_counts': ctrl_counts[ctrl_vi],
				'ctrl_conc': ctrl_conc[ctrl_vi],
				'count_ratio': count_ratio,
				'conc_ratio': conc_ratio,
				'mean_count_ratio': np.nanmean(count_ratio),
				'mean_conc_ratio': np.nanmean(conc_ratio),
			})

		if not results:
			print('Skipping: no variant data collected.')
			return

		# -------------------------------------------------------------- #
		# Export detail CSV                                                #
		# -------------------------------------------------------------- #
		csv_path = os.path.join(plotOutDir, 'kcat_enzyme_reduction_detail.csv')
		with open(csv_path, 'w', newline='', encoding='utf-8') as fh:
			writer = csv.writer(fh)
			writer.writerow([
				'variant_index', 'trl_eff', 'is_kcat', 'enzyme_id',
				'count_variant', 'count_control', 'count_ratio',
				'conc_variant', 'conc_control', 'conc_ratio',
				'dt_variant', 'dt_control',
			])
			for r in results:
				for ei, eid in enumerate(target_cat_ids):
					writer.writerow([
						r['vi'], r['trl_eff'], r['is_kcat'], eid,
						f'{r["counts"][ei]:.2f}',
						f'{r["ctrl_counts"][ei]:.2f}',
						f'{r["count_ratio"][ei]:.4f}'
						if np.isfinite(r['count_ratio'][ei]) else 'NaN',
						f'{r["conc"][ei]:.6e}',
						f'{r["ctrl_conc"][ei]:.6e}',
						f'{r["conc_ratio"][ei]:.4f}'
						if np.isfinite(r['conc_ratio'][ei]) else 'NaN',
						f'{r["dt"]:.2f}',
						f'{r["ctrl_dt"]:.2f}',
					])
		print(f'Wrote {csv_path}')

		# -------------------------------------------------------------- #
		# Plot 1: Scatter — doubling time vs mean enzyme count ratio       #
		# -------------------------------------------------------------- #
		fig, ax = plt.subplots(figsize=(8, 6))

		for is_kcat, color, label in [(False, 'tab:blue', 'No kcat'),
									  (True, 'tab:red', 'With kcat')]:
			subset = [r for r in results if r['is_kcat'] == is_kcat]
			if not subset:
				continue
			dts = [r['dt'] for r in subset]
			ratios = [r['mean_count_ratio'] for r in subset]
			ax.scatter(dts, ratios, c=color, label=label, s=40, alpha=0.8)
			for r in subset:
				ax.annotate(
					f'{r["trl_eff"]:.2f}',
					(r['dt'], r['mean_count_ratio']),
					fontsize=6, textcoords='offset points', xytext=(4, 4))

		ax.axhline(1.0, color='gray', ls='--', lw=0.8)
		ax.set_xlabel('Doubling Time (min)')
		ax.set_ylabel('Mean Enzyme Count Ratio (variant / control)')
		ax.set_title(
			'Enzyme count reduction vs doubling time\n'
			f'({n_enzymes} enzymes with kcat > 1, last {N_GENS_TO_USE} gens)')
		ax.legend()
		plt.tight_layout()
		exportFigure(plt, plotOutDir, plotOutFileName + '_scatter_dt', metadata)
		plt.close('all')

		# -------------------------------------------------------------- #
		# Plot 2: Scatter — trl_eff vs mean enzyme count ratio             #
		# -------------------------------------------------------------- #
		fig, ax = plt.subplots(figsize=(8, 6))

		for is_kcat, color, label in [(False, 'tab:blue', 'No kcat'),
									  (True, 'tab:red', 'With kcat')]:
			subset = [r for r in results if r['is_kcat'] == is_kcat]
			if not subset:
				continue
			trl_effs = [r['trl_eff'] for r in subset]
			ratios = [r['mean_count_ratio'] for r in subset]
			ax.scatter(trl_effs, ratios, c=color, label=label, s=40, alpha=0.8)

		ax.axhline(1.0, color='gray', ls='--', lw=0.8)
		ax.set_xlabel('Translation Efficiency')
		ax.set_ylabel('Mean Enzyme Count Ratio (variant / control)')
		ax.set_title(
			'Enzyme count reduction vs translation efficiency\n'
			f'({n_enzymes} enzymes with kcat > 1, last {N_GENS_TO_USE} gens)')
		ax.legend()
		plt.tight_layout()
		exportFigure(plt, plotOutDir, plotOutFileName + '_scatter_trl_eff',
					 metadata)
		plt.close('all')

		# -------------------------------------------------------------- #
		# Plot 3: Heatmap — per-enzyme ratios for kcat variants            #
		# -------------------------------------------------------------- #
		kcat_results = [r for r in results if r['is_kcat']]
		if not kcat_results:
			print('Skipping heatmap: no kcat variant data.')
			return

		# Sort columns by doubling time
		kcat_results.sort(key=lambda r: r['dt'])
		ratio_matrix = np.array(
			[r['count_ratio'] for r in kcat_results]).T  # (n_enzymes, n_variants)
		col_labels = [
			f'trl={r["trl_eff"]:.2f}\ndt={r["dt"]:.0f}'
			for r in kcat_results]

		# Sort rows by mean ratio (most reduced first)
		row_means = np.nanmean(ratio_matrix, axis=1)
		row_order = np.argsort(row_means)
		ratio_matrix = ratio_matrix[row_order]
		sorted_enzyme_ids = [target_cat_ids[i] for i in row_order]

		fig_h = max(8, n_enzymes * 0.22 + 2)
		fig_w = max(10, len(kcat_results) * 0.6 + 4)
		fig, ax = plt.subplots(figsize=(fig_w, fig_h))

		vmin = np.nanmin(ratio_matrix)
		vmax = np.nanmax(ratio_matrix)
		# Ensure vcenter=1.0 is within range
		if vmin >= 1.0:
			vmin = 0.9
		if vmax <= 1.0:
			vmax = 1.1
		norm = TwoSlopeNorm(vcenter=1.0, vmin=vmin, vmax=vmax)
		cmap = plt.cm.RdBu_r.copy()
		cmap.set_bad('#cccccc')

		im = ax.imshow(
			ratio_matrix, aspect='auto', cmap=cmap, norm=norm,
			interpolation='nearest')

		cbar = fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02)
		cbar.set_label('Count Ratio (variant / control)', fontsize=9)

		ax.set_xticks(np.arange(len(kcat_results)))
		ax.set_xticklabels(col_labels, fontsize=6, rotation=45, ha='right')
		ax.set_yticks(np.arange(n_enzymes))
		ax.set_yticklabels(sorted_enzyme_ids, fontsize=5)

		ax.set_xlabel('Kcat Variant (sorted by doubling time)', fontsize=10)
		ax.set_ylabel('Enzyme (sorted by mean ratio)', fontsize=10)
		ax.set_title(
			'Per-enzyme count ratio (variant / kcat control)\n'
			f'({n_enzymes} enzymes with kcat > 1, last {N_GENS_TO_USE} gens)',
			fontsize=11)

		plt.tight_layout()
		exportFigure(plt, plotOutDir, plotOutFileName + '_heatmap', metadata)
		plt.close('all')


if __name__ == '__main__':
	Plot().cli()
