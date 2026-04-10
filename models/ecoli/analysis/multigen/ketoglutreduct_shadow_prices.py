"""
Shadow prices and reduced costs for metabolites and reactions involved in the
KETOGLUTREDUCT-RXN futile cycle.

Shadow prices reveal metabolite scarcity/surplus in the FBA dual; reduced costs
show whether reactions are at their optimality bounds.

4 panels (2 rows x 2 cols):
  Row 1: Shadow prices (all target metabolites) | Shadow prices (excl. PROTON)
  Row 2: Reduced costs over time                | Reduced cost vs flux scatter

Also exports a TSV summary to plotOutDir.
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

# Metabolites for shadow price panels.
SHADOW_METABOLITES = [
	'PROTON[c]',
	'NAD[c]',
	'NADH[c]',
	'2-KETOGLUTARATE[c]',
	'R-2-HYDROXYGLUTARATE[c]',
]
# Classification: (H)omeostatic, (K)inetic-only, (N)either.
MET_CLASS = {
	'PROTON[c]': 'N',
	'NAD[c]': 'H',
	'NADH[c]': 'H',
	'2-KETOGLUTARATE[c]': 'H',
	'R-2-HYDROXYGLUTARATE[c]': 'K',
}

# Reactions for reduced cost panels.
REDUCED_COST_RXNS = [
	'KETOGLUTREDUCT-RXN',
	'KETOGLUTREDUCT-RXN (reverse)',
	'RXN-16701',
	'RXN-16701 (reverse)',
	'PGLYCDEHYDROG-RXN',
	'PGLYCDEHYDROG-RXN (reverse)',
	'RXN-14932',
]

IGNORE_FIRST_N_GENS = 0


class Plot(multigenAnalysisPlot.MultigenAnalysisPlot):
	def do_plot(self, seedOutDir, plotOutDir, plotOutFileName, simDataFile,
				validationDataFile, metadata):

		with open(simDataFile, 'rb') as f:
			sim_data = pickle.load(f)

		cell_density = sim_data.constants.cell_density

		cell_paths = self.ap.get_cells(
			generation=np.arange(IGNORE_FIRST_N_GENS, self.ap.n_generation),
			only_successful=True)
		if len(cell_paths) == 0:
			print('Skipping analysis -- no cells found.')
			return

		# --- Build index lookups from first cell ---
		fba0 = TableReader(os.path.join(cell_paths[0], 'simOut', 'FBAResults'))
		output_mol_ids = list(fba0.readAttribute('outputMoleculeIDs'))
		listener_rxn_ids = list(fba0.readAttribute('reactionIDs'))

		# Flexibly find shadow-price column indexes (substring match).
		met_indexes = {}
		for met in SHADOW_METABOLITES:
			for i, oid in enumerate(output_mol_ids):
				if met in oid or oid == met:
					met_indexes[met] = i
					break
			else:
				print(f'Warning: {met} not found in outputMoleculeIDs.')

		# Reaction indexes for reduced costs and fluxes.
		rxn_id_to_idx = {r: i for i, r in enumerate(listener_rxn_ids)}
		rc_indexes = {}
		for rxn in REDUCED_COST_RXNS:
			if rxn in rxn_id_to_idx:
				rc_indexes[rxn] = rxn_id_to_idx[rxn]
			else:
				print(f'Warning: {rxn} not found in reactionIDs.')

		# Column indexes to read.
		sp_cols = sorted(met_indexes.values())
		sp_col_remap = {c: i for i, c in enumerate(sp_cols)}
		rc_cols = sorted(rc_indexes.values())
		rc_col_remap = {c: i for i, c in enumerate(rc_cols)}

		# --- Accumulate traces ---
		time_all = []
		sp_all = []       # shadow prices
		rc_all = []       # reduced costs
		flux_all = []     # reaction fluxes (for scatter)
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

				fba = TableReader(os.path.join(sim_out, 'FBAResults'))
				shadow = fba.readColumn('shadowPrices',
										indices=sp_cols, squeeze=False)[1:]
				rcosts = fba.readColumn('reducedCosts',
										indices=rc_cols, squeeze=False)[1:]
				rxn_fluxes = fba.readColumn('reactionFluxes',
											indices=rc_cols, squeeze=False)[1:]
			except Exception as e:
				print(f'Skipped cell {sim_out}: {e!r}')
				continue

			t_sec = t_raw - t_raw[0]
			t_h = (t_sec + time_offset) / 3600.0
			gen_boundary_times.append(time_offset / 3600.0)
			time_offset += (t_raw[-1] - t_raw[0])

			# Convert fluxes to mmol/gDCW/h.
			conversion = (dry_mass / cell_mass
						  * cell_density.asNumber(MASS_UNITS / VOLUME_UNITS))
			fluxes = (
				(COUNTS_UNITS / MASS_UNITS / TIME_UNITS)
				* (rxn_fluxes / conversion[:, np.newaxis])
			).asNumber(units.mmol / units.g / units.h)
			fluxes[np.isinf(fluxes)] = np.nan

			time_all.append(t_h)
			sp_all.append(shadow)
			rc_all.append(rcosts)
			flux_all.append(fluxes)

		if not time_all:
			print('No cells processed.')
			return

		time_cat = np.concatenate(time_all)
		sp_cat = np.vstack(sp_all)
		rc_cat = np.vstack(rc_all)
		flux_cat = np.vstack(flux_all)

		def add_gen_lines(ax):
			for gt in gen_boundary_times[1:]:
				ax.axvline(gt, color='gray', lw=0.5, ls='--', alpha=0.4)

		# --- Plot ---
		fig, axes = plt.subplots(2, 2, figsize=(14, 10))

		# Panel 1: all shadow prices
		ax1 = axes[0, 0]
		colors = ['#d62728', '#1f77b4', '#ff7f0e', '#2ca02c', '#9467bd']
		for mi, met in enumerate(SHADOW_METABOLITES):
			if met not in met_indexes:
				continue
			ci = sp_col_remap[met_indexes[met]]
			cls = MET_CLASS.get(met, '?')
			ax1.plot(time_cat, sp_cat[:, ci], lw=0.7, alpha=0.8,
					 color=colors[mi % len(colors)],
					 label=f'{met} ({cls})')
		ax1.set_title('Shadow prices (all)', fontsize=9)
		ax1.set_ylabel('Shadow price', fontsize=7)
		ax1.legend(fontsize=5, loc='best')
		ax1.tick_params(labelsize=7)
		add_gen_lines(ax1)

		# Panel 2: shadow prices excluding PROTON[c]
		ax2 = axes[0, 1]
		for mi, met in enumerate(SHADOW_METABOLITES):
			if met == 'PROTON[c]' or met not in met_indexes:
				continue
			ci = sp_col_remap[met_indexes[met]]
			cls = MET_CLASS.get(met, '?')
			ax2.plot(time_cat, sp_cat[:, ci], lw=0.7, alpha=0.8,
					 color=colors[mi % len(colors)],
					 label=f'{met} ({cls})')
		ax2.set_title('Shadow prices (excl. PROTON)', fontsize=9)
		ax2.set_ylabel('Shadow price', fontsize=7)
		ax2.legend(fontsize=5, loc='best')
		ax2.tick_params(labelsize=7)
		add_gen_lines(ax2)

		# Panel 3: reduced costs over time
		ax3 = axes[1, 0]
		rc_colors = plt.cm.tab10(np.linspace(0, 1, len(REDUCED_COST_RXNS)))
		for ri, rxn in enumerate(REDUCED_COST_RXNS):
			if rxn not in rc_indexes:
				continue
			ci = rc_col_remap[rc_indexes[rxn]]
			ax3.plot(time_cat, rc_cat[:, ci], lw=0.7, alpha=0.8,
					 color=rc_colors[ri], label=rxn[:30])
		ax3.set_title('Reduced costs', fontsize=9)
		ax3.set_ylabel('Reduced cost', fontsize=7)
		ax3.set_xlabel('Time (h)', fontsize=8)
		ax3.legend(fontsize=4.5, loc='best')
		ax3.tick_params(labelsize=7)
		add_gen_lines(ax3)

		# Panel 4: reduced cost vs flux magnitude scatter
		ax4 = axes[1, 1]
		for ri, rxn in enumerate(REDUCED_COST_RXNS):
			if rxn not in rc_indexes:
				continue
			ci = rc_col_remap[rc_indexes[rxn]]
			ax4.scatter(np.abs(flux_cat[:, ci]), rc_cat[:, ci],
						s=1, alpha=0.3, color=rc_colors[ri],
						label=rxn[:30])
		ax4.set_title('Reduced cost vs |flux|', fontsize=9)
		ax4.set_xlabel('|Flux| (mmol/gDCW/h)', fontsize=7)
		ax4.set_ylabel('Reduced cost', fontsize=7)
		ax4.legend(fontsize=4.5, loc='best')
		ax4.tick_params(labelsize=7)

		fig.suptitle('KETOGLUTREDUCT-RXN shadow prices & reduced costs',
					 fontsize=12, y=1.002)
		plt.tight_layout()
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close('all')

		# --- TSV export ---
		tsv_path = os.path.join(plotOutDir,
								'ketoglutreduct_shadow_prices_data.tsv')
		header_cols = ['time_h']
		data_cols = [time_cat]

		for met in SHADOW_METABOLITES:
			if met in met_indexes:
				header_cols.append(f'SP_{met}')
				data_cols.append(sp_cat[:, sp_col_remap[met_indexes[met]]])

		for rxn in REDUCED_COST_RXNS:
			if rxn in rc_indexes:
				header_cols.append(f'RC_{rxn}')
				data_cols.append(rc_cat[:, rc_col_remap[rc_indexes[rxn]]])
				header_cols.append(f'Flux_{rxn}')
				data_cols.append(flux_cat[:, rc_col_remap[rc_indexes[rxn]]])

		data_matrix = np.column_stack(data_cols)
		np.savetxt(tsv_path, data_matrix, delimiter='\t',
				   header='\t'.join(header_cols), comments='')
		print(f'Wrote {tsv_path}')


if __name__ == '__main__':
	Plot().cli()
