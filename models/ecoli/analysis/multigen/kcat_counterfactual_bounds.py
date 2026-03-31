"""
Counterfactual kcat bound visualization for wildtype simulations.

Runs on wildtype sims (no kcat bounds active) to show what the smoothed_max
bounds *would* cut off.  For each reaction, overlays the actual flux with
the hypothetical bound (kcat * [E]) to reveal where the bound would bite.

Reactions are automatically classified into behavioral profiles:
  - Spiking:          max(flux) / median(flux) > SPIKE_RATIO
  - Always saturated: mean(flux/bound) > SATURATED_THRESHOLD
  - Sparse:           >50% of timesteps have flux = 0
  - Normal:           everything else

Produces a multi-page PDF grouped by category, with N_PER_CATEGORY
representative reactions from each profile.  Each page has two panels:
  Top:    Flux and hypothetical bound over time, shaded where flux > bound.
  Bottom: flux/bound ratio over time, horizontal line at 1.0.
"""

import os
import pickle

import numpy as np
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

from models.ecoli.analysis import multigenAnalysisPlot
from models.ecoli.processes.metabolism import (
	COUNTS_UNITS, VOLUME_UNITS, TIME_UNITS, MASS_UNITS)
from wholecell.io.tablereader import TableReader
from wholecell.utils import units

# Number of representative reactions per behavioral category.
N_PER_CATEGORY = 15

# Classification thresholds.
SPIKE_RATIO = 5.0
SATURATED_THRESHOLD = 0.8
SPARSE_ZERO_FRACTION = 0.5

# Ignore first N generations (initialization transient).
IGNORE_FIRST_N_GENS = 4

# Path to the smoothed_max kcat estimates TSV, relative to project root.
KCAT_TSV_RELPATH = os.path.join(
	'reconstruction', 'ecoli', 'flat', 'kcat_estimates',
	'kcat_estimates_smoothed_max.tsv')


def _load_kcat_estimates_from_tsv(tsv_path):
	"""Parse the smoothed_max kcat TSV into {(rxn_id, cat_id): kcat}.

	The TSV has a comment line (starting with '#') followed by a
	tab-separated header with quoted column names, then data rows with
	quoted strings and unquoted floats.
	"""
	estimates = {}
	with open(tsv_path) as f:
		# Skip comment lines
		for line in f:
			if not line.startswith('#'):
				header_line = line
				break
		else:
			return estimates

		# Parse header: strip quotes
		headers = [h.strip().strip('"') for h in header_line.strip().split('\t')]
		ri_col = headers.index('reaction_id')
		ci_col = headers.index('catalyst_id')
		kcat_col = headers.index('kcat_estimate')

		for line in f:
			parts = line.strip().split('\t')
			if len(parts) < 3:
				continue
			rxn_id = parts[ri_col].strip('"')
			cat_id = parts[ci_col].strip('"')
			kcat = float(parts[kcat_col].strip('"'))
			estimates[(rxn_id, cat_id)] = kcat
	return estimates


class Plot(multigenAnalysisPlot.MultigenAnalysisPlot):
	def do_plot(self, seedOutDir, plotOutDir, plotOutFileName, simDataFile,
				validationDataFile, metadata):

		with open(simDataFile, 'rb') as f:
			sim_data = pickle.load(f)

		# Find the TSV file relative to the wcEcoli project root.
		# Walk up from simDataFile to find the project root containing
		# 'reconstruction/ecoli/flat/'.
		tsv_path = None
		candidate = os.path.dirname(os.path.abspath(__file__))
		for _ in range(10):
			candidate = os.path.dirname(candidate)
			test = os.path.join(candidate, KCAT_TSV_RELPATH)
			if os.path.exists(test):
				tsv_path = test
				break

		if tsv_path is None:
			print(f'Skipping: cannot find {KCAT_TSV_RELPATH}')
			return

		kcat_estimates = _load_kcat_estimates_from_tsv(tsv_path)
		if not kcat_estimates:
			print('Skipping: no kcat estimates found in TSV.')
			return

		cell_paths = self.ap.get_cells(
			generation=np.arange(IGNORE_FIRST_N_GENS, self.ap.n_generation),
			only_successful=True)

		if len(cell_paths) == 0:
			print('Skipping analysis -- not enough simulations run.')
			return

		# --- Build per-reaction lookup ------------------------------------
		fba_reader = TableReader(
			os.path.join(cell_paths[0], 'simOut', 'FBAResults'))
		listener_rxn_ids = fba_reader.readAttribute('reactionIDs')
		listener_cat_ids = fba_reader.readAttribute('catalyst_ids')

		rxn_id_to_idx = {r: i for i, r in enumerate(listener_rxn_ids)}
		cat_id_to_idx = {c: i for i, c in enumerate(listener_cat_ids)}

		# rxn_to_pairs: rxn_id -> [(cat_listener_idx, kcat), ...]
		rxn_to_pairs = {}
		for (rxn_id, cat_id), kcat in kcat_estimates.items():
			if rxn_id not in rxn_id_to_idx or cat_id not in cat_id_to_idx:
				continue
			rxn_to_pairs.setdefault(rxn_id, []).append(
				(cat_id_to_idx[cat_id], kcat))

		if not rxn_to_pairs:
			print('Skipping: no kcat-estimated reactions found in listener IDs.')
			return

		constrained_rxn_ids = list(rxn_to_pairs.keys())
		rxn_indexes = np.array([rxn_id_to_idx[r] for r in constrained_rxn_ids])

		all_cat_indexes = sorted({ci for pairs in rxn_to_pairs.values()
								  for ci, _ in pairs})
		cat_idx_local = {ci: li for li, ci in enumerate(all_cat_indexes)}

		cell_density = sim_data.constants.cell_density
		n_rxns = len(constrained_rxn_ids)

		# --- Accumulate time traces ---------------------------------------
		time_all = []
		flux_all = []
		bound_all = []
		gen_boundary_times = []

		time_offset = 0.0

		for cell_path in cell_paths:
			sim_out = os.path.join(cell_path, 'simOut')
			try:
				t_raw = TableReader(
					os.path.join(sim_out, 'Main')
				).readColumn('time', squeeze=True)[1:]

				counts_to_molar = TableReader(
					os.path.join(sim_out, 'EnzymeKinetics')
				).readColumn('countsToMolar', squeeze=True)[1:]

				cell_mass = TableReader(
					os.path.join(sim_out, 'Mass')
				).readColumn('cellMass', squeeze=True)[1:]

				dry_mass = TableReader(
					os.path.join(sim_out, 'Mass')
				).readColumn('dryMass', squeeze=True)[1:]

				cat_counts = TableReader(
					os.path.join(sim_out, 'FBAResults')
				).readColumn('catalyst_counts',
							 indices=np.array(all_cat_indexes),
							 squeeze=False)[1:]

				rxn_fluxes_raw = TableReader(
					os.path.join(sim_out, 'FBAResults')
				).readColumn('reactionFluxes',
							 indices=rxn_indexes,
							 squeeze=False)[1:]

			except Exception as e:
				print(f'Ignored exception reading {sim_out}: {e!r}')
				continue

			t_sec = t_raw - t_raw[0]
			t_h = (t_sec + time_offset) / 3600.
			gen_boundary_times.append(time_offset / 3600.)
			time_offset += (t_raw[-1] - t_raw[0])

			# Reaction fluxes -> mmol/g DCW/h
			conversion_coeffs = (dry_mass / cell_mass
								 * cell_density.asNumber(MASS_UNITS / VOLUME_UNITS))
			fluxes = (
				(COUNTS_UNITS / MASS_UNITS / TIME_UNITS)
				* (rxn_fluxes_raw / conversion_coeffs[:, np.newaxis])
			).asNumber(units.mmol / units.g / units.h)
			fluxes[np.isinf(fluxes)] = np.nan

			# Hypothetical bounds -> mmol/g DCW/h: max over catalysts of kcat * [E]
			bounds = np.zeros_like(fluxes)
			for ri, rxn_id in enumerate(constrained_rxn_ids):
				for cat_idx, kcat in rxn_to_pairs[rxn_id]:
					li = cat_idx_local[cat_idx]
					conc = cat_counts[:, li] * counts_to_molar
					bounds[:, ri] = np.maximum(bounds[:, ri], kcat * conc)

			time_all.append(t_h)
			flux_all.append(fluxes)
			bound_all.append(bounds)

		if not time_all:
			print('No cells processed successfully.')
			return

		time_cat = np.concatenate(time_all)
		flux_cat = np.vstack(flux_all)
		bound_cat = np.vstack(bound_all)

		# --- Classify reactions into behavioral profiles ------------------
		categories = {}  # category -> list of (ri, sort_key)
		for ri in range(n_rxns):
			f = flux_cat[:, ri]
			b = bound_cat[:, ri]
			valid = (b > 0) & np.isfinite(f)

			if not valid.any():
				continue

			f_valid = f[valid]
			b_valid = b[valid]

			median_f = np.median(np.abs(f_valid))
			max_f = np.max(np.abs(f_valid))
			frac_zero = np.mean(np.abs(f_valid) < 1e-12)

			with np.errstate(divide='ignore', invalid='ignore'):
				ratio = f_valid / b_valid
			ratio = ratio[np.isfinite(ratio)]
			mean_ratio = np.mean(ratio) if len(ratio) > 0 else 0.0

			if median_f > 0 and max_f / median_f > SPIKE_RATIO:
				cat = 'spiking'
				sort_key = max_f / median_f
			elif mean_ratio > SATURATED_THRESHOLD:
				cat = 'always_saturated'
				sort_key = mean_ratio
			elif frac_zero > SPARSE_ZERO_FRACTION:
				cat = 'sparse'
				sort_key = frac_zero
			else:
				cat = 'normal'
				sort_key = mean_ratio

			categories.setdefault(cat, []).append((ri, sort_key))

		# Sort each category by sort_key descending, pick top N
		selected = {}  # category -> [ri, ...]
		for cat, entries in categories.items():
			entries.sort(key=lambda x: -x[1])
			selected[cat] = [ri for ri, _ in entries[:N_PER_CATEGORY]]

		# Print summary
		print('\n' + '=' * 80)
		print('Counterfactual bound classification summary')
		print('-' * 80)
		for cat in ['always_saturated', 'spiking', 'sparse', 'normal']:
			if cat not in categories:
				continue
			total = len(categories[cat])
			shown = len(selected.get(cat, []))
			print(f'  {cat:<20s}: {total:4d} reactions total, showing top {shown}')
			for ri in selected.get(cat, []):
				print(f'    - {constrained_rxn_ids[ri]}')
		print('=' * 80 + '\n')

		# --- Multi-page PDF -----------------------------------------------
		category_order = ['always_saturated', 'spiking', 'sparse', 'normal']
		category_labels = {
			'always_saturated': 'Always Saturated (mean flux/bound > {:.1f})'.format(SATURATED_THRESHOLD),
			'spiking': 'Spiking (max/median > {:.0f}x)'.format(SPIKE_RATIO),
			'sparse': 'Sparse (>{:.0%} zero-flux timesteps)'.format(SPARSE_ZERO_FRACTION),
			'normal': 'Normal',
		}

		pdf_path = os.path.join(plotOutDir, plotOutFileName + '.pdf')
		with PdfPages(pdf_path) as pdf:
			for cat in category_order:
				if cat not in selected or not selected[cat]:
					continue

				# Category header page
				fig = plt.figure(figsize=(12, 2))
				fig.text(0.5, 0.5, category_labels.get(cat, cat),
						 ha='center', va='center', fontsize=18, fontweight='bold')
				fig.text(0.5, 0.25,
						 f'{len(categories[cat])} reactions in this category, '
						 f'showing top {len(selected[cat])}',
						 ha='center', va='center', fontsize=12, color='gray')
				pdf.savefig(fig)
				plt.close(fig)

				for ri in selected[cat]:
					rxn_id = constrained_rxn_ids[ri]
					fig, (ax_top, ax_bot) = plt.subplots(
						2, 1, figsize=(12, 6), sharex=True)

					fig.suptitle(
						f'{rxn_id}  [{cat}]', fontsize=10)

					f = flux_cat[:, ri]
					b = bound_cat[:, ri]

					# Top panel: flux and hypothetical bound
					ax_top.plot(time_cat, f,
								lw=0.6, color='#1f77b4', alpha=0.85, label='flux')
					ax_top.plot(time_cat, b,
								lw=0.6, color='#d62728', alpha=0.85,
								ls='--', label='hypothetical bound')

					# Shade where flux > bound (would be clipped)
					exceeds = f > b
					if exceeds.any():
						ax_top.fill_between(
							time_cat, b, f,
							where=exceeds, alpha=0.2, color='#d62728',
							label='would be clipped')

					ax_top.set_ylabel('mmol/g DCW/h', fontsize=9)
					ax_top.legend(fontsize=8, loc='upper right')

					# Bottom panel: flux/bound ratio
					with np.errstate(divide='ignore', invalid='ignore'):
						ratio = f / b
					ratio[~np.isfinite(ratio)] = np.nan
					ax_bot.plot(time_cat, ratio,
								lw=0.6, color='#2ca02c', alpha=0.85)
					ax_bot.axhline(1.0, color='k', lw=0.8, ls='--', alpha=0.5)
					ax_bot.set_ylabel('flux / bound', fontsize=9)
					ax_bot.set_xlabel('Time (h)', fontsize=9)

					finite_ratio = ratio[np.isfinite(ratio)]
					if len(finite_ratio) > 0:
						ax_bot.set_ylim(-0.05, max(1.2, np.nanmax(ratio) * 1.1))
					else:
						ax_bot.set_ylim(-0.05, 1.2)

					# Generation boundaries
					for gt in gen_boundary_times[1:]:
						ax_top.axvline(gt, color='gray', lw=0.5, ls='--', alpha=0.5)
						ax_bot.axvline(gt, color='gray', lw=0.5, ls='--', alpha=0.5)

					for ax in (ax_top, ax_bot):
						ax.set_xlim(time_cat[0], time_cat[-1])
						ax.tick_params(labelsize=7)

					plt.tight_layout()
					pdf.savefig(fig)
					plt.close(fig)

		print(f'PDF written to {pdf_path}')


if __name__ == '__main__':
	Plot().cli()
