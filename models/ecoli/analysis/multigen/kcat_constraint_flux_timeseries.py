"""
Per-reaction flux/bound/enzyme time-series for all kcat-constrained reactions.

For sims with active kcat constraints (new_gene_internal_shift_with_kcat
variant), shows for each constrained reaction its flux, kcat upper bound,
enzyme concentration, and enzyme counts over time across multiple generations.

Page 1: GFP protein counts and concentration overview.
Pages 2+: One page per (rxn, catalyst) pair with three panels (flux vs bound,
enzyme concentration, enzyme counts), sorted by saturation rate descending.

Skips gracefully if selected_kcat_estimates is absent.
"""

import os
import pickle

import numpy as np
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

from models.ecoli.analysis import multigenAnalysisPlot
from models.ecoli.processes.metabolism import (
	COUNTS_UNITS, VOLUME_UNITS, TIME_UNITS, MASS_UNITS)
from models.ecoli.sim.variants.new_gene_internal_shift import (
	determine_new_gene_ids_and_indices)
from wholecell.analysis.analysis_tools import read_stacked_bulk_molecules
from wholecell.io.tablereader import TableReader
from wholecell.utils import units

# Only include reactions with kcat >= this threshold.
MIN_KCAT_ESTIMATE = 1.0

# Bound is considered "active" when flux exceeds this fraction of the bound.
SATURATION_THRESHOLD = 0.9

IGNORE_FIRST_N_GENS = 0


class Plot(multigenAnalysisPlot.MultigenAnalysisPlot):
	def do_plot(self, seedOutDir, plotOutDir, plotOutFileName, simDataFile,
				validationDataFile, metadata):

		with open(simDataFile, 'rb') as f:
			sim_data = pickle.load(f)

		metabolism = sim_data.process.metabolism
		if not hasattr(metabolism, 'selected_kcat_estimates'):
			print('Skipping: sim_data does not have selected_kcat_estimates.'
				  ' Run with a kcat constraint variant.')
			return

		selected_estimates = metabolism.selected_kcat_estimates

		# --- Build target reaction set ------------------------------------
		targets = {
			pair for pair, kcat in selected_estimates.items()
			if kcat >= MIN_KCAT_ESTIMATE
		}

		if not targets:
			print('No target reactions found.')
			return

		# --- Get new gene monomer IDs for GFP page ------------------------
		try:
			(new_gene_mRNA_ids, _, new_gene_monomer_ids,
			 _) = determine_new_gene_ids_and_indices(sim_data)
		except Exception as e:
			print(f'Could not determine new gene IDs: {e!r}')
			new_gene_monomer_ids = []

		# --- Get cell paths -----------------------------------------------
		cell_paths = self.ap.get_cells(
			generation=np.arange(IGNORE_FIRST_N_GENS, self.ap.n_generation),
			only_successful=True)
		if len(cell_paths) == 0:
			print('Skipping analysis -- not enough simulations run.')
			return

		# --- Build index lookups from first cell --------------------------
		fba_reader = TableReader(
			os.path.join(cell_paths[0], 'simOut', 'FBAResults'))
		listener_rxn_ids = fba_reader.readAttribute('reactionIDs')
		listener_cat_ids = fba_reader.readAttribute('catalyst_ids')

		rxn_id_to_idx = {r: i for i, r in enumerate(listener_rxn_ids)}
		cat_id_to_idx = {c: i for i, c in enumerate(listener_cat_ids)}

		valid_targets = []
		target_kcats = {}
		for rxn_id, cat_id in sorted(targets):
			if rxn_id in rxn_id_to_idx and cat_id in cat_id_to_idx:
				valid_targets.append((rxn_id, cat_id))
				target_kcats[(rxn_id, cat_id)] = selected_estimates[
					(rxn_id, cat_id)]

		if not valid_targets:
			print('No valid target reactions found in listener.')
			return

		n_targets = len(valid_targets)
		rxn_indexes = np.array([rxn_id_to_idx[r] for r, _ in valid_targets])
		cat_indexes = np.array([cat_id_to_idx[c] for _, c in valid_targets])

		cell_density = sim_data.constants.cell_density

		# --- Read GFP data ------------------------------------------------
		gfp_counts_stacked = None
		gfp_conc_all = []
		gfp_time_all = []
		gfp_gen_boundaries = []
		gfp_time_offset = 0.0

		if new_gene_monomer_ids:
			for cell_path in cell_paths:
				sim_out = os.path.join(cell_path, 'simOut')
				try:
					t_raw = TableReader(
						os.path.join(sim_out, 'Main')
					).readColumn('time', squeeze=True)[1:]
					counts_to_molar = TableReader(
						os.path.join(sim_out, 'EnzymeKinetics')
					).readColumn('countsToMolar', squeeze=True)[1:]
				except Exception:
					continue

				t_sec = t_raw - t_raw[0]
				t_h = (t_sec + gfp_time_offset) / 3600.
				gfp_gen_boundaries.append(gfp_time_offset / 3600.)
				gfp_time_offset += (t_raw[-1] - t_raw[0])
				gfp_time_all.append(t_h)

			# Read stacked bulk molecules for GFP counts
			(gfp_counts_stacked,) = read_stacked_bulk_molecules(
				cell_paths, new_gene_monomer_ids)
			# Remove first timestep per cell to match [1:] slicing
			# read_stacked_bulk_molecules doesn't skip first timestep,
			# so we need to handle the time array accordingly
			gfp_time_for_counts = np.concatenate([
				TableReader(os.path.join(cp, 'simOut', 'Main')
				).readColumn('time', squeeze=True)
				for cp in cell_paths])
			# Recompute GFP time to match stacked counts (no [1:] skip)
			gfp_time_offset2 = 0.0
			gfp_time_all2 = []
			gfp_gen_boundaries2 = []
			for cell_path in cell_paths:
				sim_out = os.path.join(cell_path, 'simOut')
				try:
					t_raw = TableReader(
						os.path.join(sim_out, 'Main')
					).readColumn('time', squeeze=True)
				except Exception:
					continue
				t_sec = t_raw - t_raw[0]
				t_h = (t_sec + gfp_time_offset2) / 3600.
				gfp_gen_boundaries2.append(gfp_time_offset2 / 3600.)
				gfp_time_offset2 += (t_raw[-1] - t_raw[0])
				gfp_time_all2.append(t_h)

			gfp_time = np.concatenate(gfp_time_all2)
			gfp_gen_boundaries = gfp_gen_boundaries2

			# Compute GFP concentration: counts * countsToMolar
			# Need countsToMolar without [1:] skip to match
			ctm_all = []
			for cell_path in cell_paths:
				sim_out = os.path.join(cell_path, 'simOut')
				try:
					ctm = TableReader(
						os.path.join(sim_out, 'EnzymeKinetics')
					).readColumn('countsToMolar', squeeze=True)
					ctm_all.append(ctm)
				except Exception:
					continue
			counts_to_molar_stacked = np.concatenate(ctm_all)
			if gfp_counts_stacked.ndim == 1:
				gfp_conc = gfp_counts_stacked * counts_to_molar_stacked
			else:
				gfp_conc = (gfp_counts_stacked
							* counts_to_molar_stacked[:, np.newaxis])

		# --- Accumulate flux/enzyme time traces per cell ------------------
		time_all = []
		flux_all = []
		cat_counts_all = []
		enzyme_conc_all = []
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
							 indices=cat_indexes,
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

			# Convert fluxes to mmol/g DCW/h
			conversion_coeffs = (dry_mass / cell_mass
								 * cell_density.asNumber(MASS_UNITS / VOLUME_UNITS))
			fluxes = (
				(COUNTS_UNITS / MASS_UNITS / TIME_UNITS)
				* (rxn_fluxes_raw / conversion_coeffs[:, np.newaxis])
			).asNumber(units.mmol / units.g / units.h)
			fluxes[np.isinf(fluxes)] = np.nan

			# Enzyme concentration (molar)
			enzyme_conc = cat_counts * counts_to_molar[:, np.newaxis]

			time_all.append(t_h)
			flux_all.append(fluxes)
			cat_counts_all.append(cat_counts)
			enzyme_conc_all.append(enzyme_conc)

		if not time_all:
			print('No cells processed successfully.')
			return

		time_cat = np.concatenate(time_all)
		flux_cat = np.vstack(flux_all)
		cat_counts_cat = np.vstack(cat_counts_all)
		enzyme_conc_cat = np.vstack(enzyme_conc_all)

		# --- Compute bounds and saturation rates --------------------------
		# Bound = kcat * [E] for each (rxn, catalyst) pair
		bound_cat = np.full_like(flux_cat, np.nan)
		saturation_rates = np.full(n_targets, np.nan)

		for ti in range(n_targets):
			pair = valid_targets[ti]
			kcat = target_kcats[pair]
			bound_cat[:, ti] = kcat * enzyme_conc_cat[:, ti]

			# Saturation rate: fraction of timesteps where flux > 0.9 * bound
			valid = (bound_cat[:, ti] > 0) & np.isfinite(flux_cat[:, ti])
			if valid.any():
				saturation_rates[ti] = np.mean(
					flux_cat[valid, ti] > SATURATION_THRESHOLD * bound_cat[valid, ti])

		# Sort by saturation rate descending
		has_data = ~np.isnan(saturation_rates)
		sort_order = np.argsort(saturation_rates[has_data])[::-1]
		sorted_idxs = np.where(has_data)[0][sort_order]
		# Include targets without data at the end
		no_data_idxs = np.where(~has_data)[0]
		sorted_idxs = np.concatenate([sorted_idxs, no_data_idxs])

		# --- Write multi-page PDF -----------------------------------------
		pdf_path = os.path.join(plotOutDir, plotOutFileName + '.pdf')
		n_pages = 0

		with PdfPages(pdf_path) as pdf:
			# Page 1: GFP overview
			if new_gene_monomer_ids and gfp_counts_stacked is not None:
				fig, (ax_counts, ax_conc) = plt.subplots(
					2, 1, figsize=(14, 8), sharex=True)

				fig.suptitle('GFP Protein Overview', fontsize=12,
							 fontweight='bold')

				# GFP counts
				if gfp_counts_stacked.ndim == 1:
					ax_counts.plot(gfp_time, gfp_counts_stacked,
								   lw=0.8, color='green', alpha=0.85,
								   label=new_gene_monomer_ids[0])
				else:
					for m, mid in enumerate(new_gene_monomer_ids):
						ax_counts.plot(gfp_time, gfp_counts_stacked[:, m],
									   lw=0.8, alpha=0.85, label=mid)
				ax_counts.set_ylabel('Protein counts (molecules)', fontsize=9)
				ax_counts.legend(fontsize=7, loc='upper left')

				# GFP concentration
				if gfp_conc.ndim == 1:
					ax_conc.plot(gfp_time, gfp_conc,
								 lw=0.8, color='green', alpha=0.85,
								 label=new_gene_monomer_ids[0])
				else:
					for m, mid in enumerate(new_gene_monomer_ids):
						ax_conc.plot(gfp_time, gfp_conc[:, m],
									 lw=0.8, alpha=0.85, label=mid)
				ax_conc.set_ylabel('Protein concentration (M)', fontsize=9)
				ax_conc.set_xlabel('Time (h)', fontsize=9)
				ax_conc.legend(fontsize=7, loc='upper left')

				# Generation boundaries
				for gt in gfp_gen_boundaries[1:]:
					ax_counts.axvline(gt, color='gray', lw=0.5, ls='--',
									  alpha=0.5)
					ax_conc.axvline(gt, color='gray', lw=0.5, ls='--',
									alpha=0.5)

				for ax in (ax_counts, ax_conc):
					ax.set_xlim(gfp_time[0], gfp_time[-1])
					ax.tick_params(labelsize=7)

				plt.tight_layout()
				pdf.savefig(fig)
				plt.close(fig)
				n_pages += 1

			# Per-reaction pages
			for ti in sorted_idxs:
				rxn_id, cat_id = valid_targets[ti]
				kcat = target_kcats[(rxn_id, cat_id)]
				sat_rate = saturation_rates[ti]
				sat_str = f'{sat_rate:.1%}' if np.isfinite(sat_rate) else 'N/A'

				fig, (ax_flux, ax_conc, ax_counts) = plt.subplots(
					3, 1, figsize=(14, 10), sharex=True,
					gridspec_kw={'height_ratios': [2, 1, 1]})

				fig.suptitle(
					f'{rxn_id}  /  {cat_id}\n'
					f'kcat = {kcat:.2f}  |  saturation rate = {sat_str}',
					fontsize=10)

				# Panel 1: Flux vs bound
				ax_flux.plot(time_cat, bound_cat[:, ti],
							 lw=1.2, color='red', alpha=0.8, ls='--',
							 label=f'kcat bound (kcat={kcat:.2f})')
				ax_flux.plot(time_cat, flux_cat[:, ti],
							 lw=0.8, color='blue', alpha=0.85,
							 label='flux', zorder=10)
				ax_flux.set_ylabel('mmol/g DCW/h', fontsize=9)
				ax_flux.legend(fontsize=7, loc='upper right')

				# Panel 2: Enzyme concentration
				ax_conc.plot(time_cat, enzyme_conc_cat[:, ti],
							 lw=0.8, color='teal', alpha=0.85)
				ax_conc.set_ylabel('Enzyme conc (M)', fontsize=9)

				# Panel 3: Enzyme counts
				ax_counts.plot(time_cat, cat_counts_cat[:, ti],
							   lw=0.8, color='teal', alpha=0.85)
				ax_counts.set_ylabel('Enzyme counts (molecules)', fontsize=9)
				ax_counts.set_xlabel('Time (h)', fontsize=9)

				# Generation boundaries
				for gt in gen_boundary_times[1:]:
					for ax in (ax_flux, ax_conc, ax_counts):
						ax.axvline(gt, color='gray', lw=0.5, ls='--',
								   alpha=0.5)

				for ax in (ax_flux, ax_conc, ax_counts):
					ax.set_xlim(time_cat[0], time_cat[-1])
					ax.tick_params(labelsize=7)

				plt.tight_layout()
				pdf.savefig(fig)
				plt.close(fig)
				n_pages += 1

		print(f'PDF written to {pdf_path}  ({n_pages} pages, '
			  f'{n_targets} reactions)')


if __name__ == '__main__':
	Plot().cli()
