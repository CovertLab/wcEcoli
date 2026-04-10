"""
Diagnose whether KETOGLUTREDUCT-RXN forms a futile cycle by showing forward,
reverse, and net flux decomposition for the reaction family catalyzed by
PGLYCDEHYDROG-CPLX (SerA), along with enzyme levels and metabolite mass
balances.

8 panels (4 rows x 2 cols):
  Row 1: KETOGLUTREDUCT-RXN fwd/rev/net | RXN-16701 fwd/rev/net
  Row 2: PGLYCDEHYDROG-RXN fwd/rev      | RXN-14932 flux
  Row 3: Enzyme counts over time         | Enzyme concentrations (mmol/L)
  Row 4: R-2-HYDROXYGLUTARATE mass bal   | 2-KETOGLUTARATE top-10 reactions

x-axis is simulated time in hours across generations; vertical dashed lines
mark generation boundaries.
"""

import os
import pickle

import numpy as np
from matplotlib import pyplot as plt

from models.ecoli.analysis import multigenAnalysisPlot
from models.ecoli.processes.metabolism import (
	COUNTS_UNITS, VOLUME_UNITS, TIME_UNITS, MASS_UNITS)
from wholecell.analysis.analysis_tools import exportFigure
from wholecell.io.tablereader import TableReader
from wholecell.utils import units

# Reactions to trace individually (panels 1-4).
REACTION_PAIRS = [
	('KETOGLUTREDUCT-RXN', 'KETOGLUTREDUCT-RXN (reverse)'),
	('RXN-16701', 'RXN-16701 (reverse)'),
	('PGLYCDEHYDROG-RXN', 'PGLYCDEHYDROG-RXN (reverse)'),
	('RXN-14932', None),  # no reverse entry expected
]

# Enzymes to track (panels 5-6).
ENZYME_IDS = [
	'PGLYCDEHYDROG-CPLX[c]',
	'CPLX0-9749[c]',
	'2OXOGLUTARATEDEH-CPLX[c]',
]
ENZYME_LABELS = ['SerA (PGLYCDEHYDROG-CPLX)', 'YdiJ (CPLX0-9749)',
				 'SucAB/Lpd (2OXOGLUTARATEDEH-CPLX)']

# Metabolites for mass-balance panels.
MASS_BAL_METABOLITE = 'R-2-HYDROXYGLUTARATE[c]'
TOP_FLUX_METABOLITE = '2-KETOGLUTARATE[c]'
TOP_N = 10

IGNORE_FIRST_N_GENS = 0


class Plot(multigenAnalysisPlot.MultigenAnalysisPlot):
	def do_plot(self, seedOutDir, plotOutDir, plotOutFileName, simDataFile,
				validationDataFile, metadata):

		with open(simDataFile, 'rb') as f:
			sim_data = pickle.load(f)

		cell_density = sim_data.constants.cell_density
		reaction_stoich = sim_data.process.metabolism.reaction_stoich

		# --- Find reactions touching target metabolites ---
		def reactions_for_metabolite(met_id):
			"""Return {base_rxn_id: stoich_coeff} for every reaction_stoich
			entry that contains met_id.  Forward and reverse entries for the
			same base reaction are collapsed: the stoich_coeff is taken from
			the forward entry (reverse has negated stoich already)."""
			result = {}
			for rxn_id, stoich in reaction_stoich.items():
				if met_id in stoich:
					# Use the base reaction id (strip ' (reverse)' suffix)
					if rxn_id.endswith(' (reverse)'):
						continue  # skip reverse; net flux handles sign
					result[rxn_id] = stoich[met_id]
			return result

		r2hg_rxns = reactions_for_metabolite(MASS_BAL_METABOLITE)
		akg_rxns = reactions_for_metabolite(TOP_FLUX_METABOLITE)

		# --- Get cell paths ---
		cell_paths = self.ap.get_cells(
			generation=np.arange(IGNORE_FIRST_N_GENS, self.ap.n_generation),
			only_successful=True)
		if len(cell_paths) == 0:
			print('Skipping analysis -- no cells found.')
			return

		# --- Build index lookups from the first cell ---
		fba0 = TableReader(os.path.join(cell_paths[0], 'simOut', 'FBAResults'))
		listener_rxn_ids = list(fba0.readAttribute('reactionIDs'))
		listener_cat_ids = list(fba0.readAttribute('catalyst_ids'))
		base_rxn_ids = list(fba0.readAttribute('base_reaction_ids'))

		rxn_id_to_idx = {r: i for i, r in enumerate(listener_rxn_ids)}
		cat_id_to_idx = {c: i for i, c in enumerate(listener_cat_ids)}
		base_id_to_idx = {r: i for i, r in enumerate(base_rxn_ids)}

		# Resolve reaction indexes for panels 1-4.
		fwd_rev_indexes = []  # list of (fwd_idx_or_None, rev_idx_or_None)
		for fwd_id, rev_id in REACTION_PAIRS:
			fi = rxn_id_to_idx.get(fwd_id)
			ri = rxn_id_to_idx.get(rev_id) if rev_id else None
			fwd_rev_indexes.append((fi, ri))

		# Resolve enzyme indexes for panels 5-6.
		enz_indexes = []
		for eid in ENZYME_IDS:
			idx = cat_id_to_idx.get(eid)
			if idx is None:
				print(f'Warning: enzyme {eid} not found in catalyst_ids.')
			enz_indexes.append(idx)

		# Resolve base-reaction indexes for panels 7-8.
		r2hg_base_idxs = {}
		for rxn_id, coeff in r2hg_rxns.items():
			bi = base_id_to_idx.get(rxn_id)
			if bi is not None:
				r2hg_base_idxs[rxn_id] = (bi, coeff)

		akg_base_idxs = {}
		for rxn_id, coeff in akg_rxns.items():
			bi = base_id_to_idx.get(rxn_id)
			if bi is not None:
				akg_base_idxs[rxn_id] = (bi, coeff)

		# Collect all needed base-reaction column indexes.
		all_base_cols = sorted(set(
			[bi for bi, _ in r2hg_base_idxs.values()] +
			[bi for bi, _ in akg_base_idxs.values()]
		))
		base_col_remap = {bi: i for i, bi in enumerate(all_base_cols)}

		# Collect all needed reactionFlux column indexes.
		all_rxn_cols = sorted(set(
			idx for pair in fwd_rev_indexes for idx in pair if idx is not None
		))
		rxn_col_remap = {ci: i for i, ci in enumerate(all_rxn_cols)}

		# Collect all needed enzyme column indexes.
		all_enz_cols = sorted(set(i for i in enz_indexes if i is not None))
		enz_col_remap = {ci: i for i, ci in enumerate(all_enz_cols)}

		# --- Accumulate traces across generations ---
		time_all = []
		flux_all = []          # (T, n_rxn_cols) — reactionFluxes in mmol/gDCW/h
		base_flux_all = []     # (T, n_base_cols) — base_reaction_fluxes
		enz_counts_all = []    # (T, n_enz_cols)
		enz_conc_all = []      # (T, n_enz_cols) mmol/L
		gen_boundary_times = []
		time_offset = 0.0

		for cell_path in cell_paths:
			sim_out = os.path.join(cell_path, 'simOut')
			try:
				t_raw = TableReader(os.path.join(sim_out, 'Main')
					).readColumn('time', squeeze=True)[1:]
				cell_mass = TableReader(os.path.join(sim_out, 'Mass')
					).readColumn('cellMass', squeeze=True)[1:]
				dry_mass = TableReader(os.path.join(sim_out, 'Mass')
					).readColumn('dryMass', squeeze=True)[1:]
				counts_to_molar = TableReader(os.path.join(sim_out, 'EnzymeKinetics')
					).readColumn('countsToMolar', squeeze=True)[1:]

				fba_reader = TableReader(os.path.join(sim_out, 'FBAResults'))
				rxn_fluxes_raw = fba_reader.readColumn(
					'reactionFluxes', indices=all_rxn_cols, squeeze=False)[1:]
				base_fluxes_raw = fba_reader.readColumn(
					'base_reaction_fluxes', indices=all_base_cols, squeeze=False)[1:]
				cat_counts = fba_reader.readColumn(
					'catalyst_counts', indices=all_enz_cols, squeeze=False)[1:]
			except Exception as e:
				print(f'Skipped cell {sim_out}: {e!r}')
				continue

			# Time
			t_sec = t_raw - t_raw[0]
			t_h = (t_sec + time_offset) / 3600.0
			gen_boundary_times.append(time_offset / 3600.0)
			time_offset += (t_raw[-1] - t_raw[0])

			# Flux conversion: counts/volume/time -> mmol/gDCW/h
			conversion = (dry_mass / cell_mass
						  * cell_density.asNumber(MASS_UNITS / VOLUME_UNITS))
			fluxes = (
				(COUNTS_UNITS / MASS_UNITS / TIME_UNITS)
				* (rxn_fluxes_raw / conversion[:, np.newaxis])
			).asNumber(units.mmol / units.g / units.h)
			fluxes[np.isinf(fluxes)] = np.nan

			base_fluxes = (
				(COUNTS_UNITS / MASS_UNITS / TIME_UNITS)
				* (base_fluxes_raw / conversion[:, np.newaxis])
			).asNumber(units.mmol / units.g / units.h)
			base_fluxes[np.isinf(base_fluxes)] = np.nan

			# Enzyme concentrations
			enz_conc = cat_counts * counts_to_molar[:, np.newaxis]

			time_all.append(t_h)
			flux_all.append(fluxes)
			base_flux_all.append(base_fluxes)
			enz_counts_all.append(cat_counts)
			enz_conc_all.append(enz_conc)

		if not time_all:
			print('No cells processed.')
			return

		time_cat = np.concatenate(time_all)
		flux_cat = np.vstack(flux_all)
		base_flux_cat = np.vstack(base_flux_all)
		enz_counts_cat = np.vstack(enz_counts_all)
		enz_conc_cat = np.vstack(enz_conc_all)

		# --- Plot ---
		fig, axes = plt.subplots(4, 2, figsize=(14, 16), sharex=True)

		def add_gen_lines(ax):
			for gt in gen_boundary_times[1:]:
				ax.axvline(gt, color='gray', lw=0.5, ls='--', alpha=0.4)

		# --- Panels 1-4: reaction fluxes ---
		titles = [
			'KETOGLUTREDUCT-RXN', 'RXN-16701',
			'PGLYCDEHYDROG-RXN', 'RXN-14932',
		]
		for pi, ((fi, ri), title) in enumerate(zip(fwd_rev_indexes, titles)):
			row, col = divmod(pi, 2)
			ax = axes[row, col]

			if fi is not None:
				fwd = flux_cat[:, rxn_col_remap[fi]]
				ax.plot(time_cat, fwd, lw=0.7, color='#1f77b4', alpha=0.8,
						label='forward')
			else:
				fwd = np.zeros(len(time_cat))

			if ri is not None:
				rev = flux_cat[:, rxn_col_remap[ri]]
				ax.plot(time_cat, rev, lw=0.7, color='#d62728', alpha=0.8,
						label='reverse')
			else:
				rev = np.zeros(len(time_cat))

			net = fwd - rev
			ax.plot(time_cat, net, lw=0.7, color='black', ls='--', alpha=0.7,
					label='net (fwd-rev)')

			ax.set_title(title, fontsize=9)
			ax.set_ylabel('Flux (mmol/gDCW/h)', fontsize=7)
			ax.legend(fontsize=6, loc='best')
			ax.tick_params(labelsize=7)
			add_gen_lines(ax)

		# --- Panel 5: enzyme counts ---
		ax5 = axes[2, 0]
		colors = ['#1f77b4', '#ff7f0e', '#2ca02c']
		for ei, (eidx, label, c) in enumerate(
				zip(enz_indexes, ENZYME_LABELS, colors)):
			if eidx is not None:
				ax5.plot(time_cat, enz_counts_cat[:, enz_col_remap[eidx]],
						 lw=0.7, color=c, alpha=0.8, label=label)
		ax5.set_title('Enzyme counts', fontsize=9)
		ax5.set_ylabel('Counts', fontsize=7)
		ax5.legend(fontsize=5, loc='best')
		ax5.tick_params(labelsize=7)
		add_gen_lines(ax5)

		# --- Panel 6: enzyme concentrations ---
		ax6 = axes[2, 1]
		for ei, (eidx, label, c) in enumerate(
				zip(enz_indexes, ENZYME_LABELS, colors)):
			if eidx is not None:
				ax6.plot(time_cat, enz_conc_cat[:, enz_col_remap[eidx]],
						 lw=0.7, color=c, alpha=0.8, label=label)
		ax6.set_title('Enzyme concentrations', fontsize=9)
		ax6.set_ylabel('Concentration (mmol/L)', fontsize=7)
		ax6.legend(fontsize=5, loc='best')
		ax6.tick_params(labelsize=7)
		add_gen_lines(ax6)

		# --- Panel 7: R-2-HYDROXYGLUTARATE mass balance ---
		ax7 = axes[3, 0]
		if r2hg_base_idxs:
			for rxn_id, (bi, coeff) in sorted(r2hg_base_idxs.items()):
				contribution = base_flux_cat[:, base_col_remap[bi]] * coeff
				mean_val = np.nanmean(contribution)
				ax7.plot(time_cat, contribution, lw=0.7, alpha=0.7,
						 label=f'{rxn_id} (mean={mean_val:.3g})')
			ax7.axhline(0, color='black', lw=0.5, alpha=0.3)
		ax7.set_title(f'{MASS_BAL_METABOLITE} mass balance', fontsize=9)
		ax7.set_ylabel('Net flux contribution\n(mmol/gDCW/h)', fontsize=7)
		ax7.legend(fontsize=5, loc='best')
		ax7.tick_params(labelsize=7)
		ax7.set_xlabel('Time (h)', fontsize=8)
		add_gen_lines(ax7)

		# --- Panel 8: 2-KETOGLUTARATE top-10 reactions ---
		ax8 = axes[3, 1]
		if akg_base_idxs:
			# Compute mean absolute contribution for ranking.
			rxn_means = {}
			for rxn_id, (bi, coeff) in akg_base_idxs.items():
				contribution = base_flux_cat[:, base_col_remap[bi]] * coeff
				rxn_means[rxn_id] = (np.nanmean(np.abs(contribution)),
									 np.nanmean(contribution), coeff)

			top_rxns = sorted(rxn_means.keys(),
							  key=lambda r: rxn_means[r][0], reverse=True)[:TOP_N]
			labels = []
			means = []
			for rxn_id in top_rxns:
				abs_mean, signed_mean, _ = rxn_means[rxn_id]
				labels.append(rxn_id)
				means.append(signed_mean)

			y_pos = np.arange(len(labels))
			bar_colors = ['#2ca02c' if m > 0 else '#d62728' for m in means]
			ax8.barh(y_pos, means, color=bar_colors, alpha=0.8)
			ax8.set_yticks(y_pos)
			ax8.set_yticklabels(labels, fontsize=5)
			ax8.axvline(0, color='black', lw=0.5, alpha=0.3)
		ax8.set_title(f'{TOP_FLUX_METABOLITE} top-{TOP_N} by |flux|', fontsize=9)
		ax8.set_xlabel('Mean net flux contribution\n(mmol/gDCW/h)', fontsize=7)
		ax8.tick_params(labelsize=7)
		add_gen_lines(ax8)

		fig.suptitle('KETOGLUTREDUCT-RXN futile cycle diagnostic', fontsize=12,
					 y=1.002)
		plt.tight_layout()
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close('all')


if __name__ == '__main__':
	Plot().cli()
