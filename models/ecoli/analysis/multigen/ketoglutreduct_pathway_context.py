"""
Broader pathway context for the KETOGLUTREDUCT-RXN futile cycle: what fraction
of each metabolite's total turnover comes from the cycle reactions?

4 panels (2 rows x 2 cols):
  Row 1: Bar chart of all R-2-HYDROXYGLUTARATE reactions |
         Bar chart of top-20 2-KETOGLUTARATE reactions
  Row 2: Bar chart of top-20 NAD reactions               |
         Summary text table (production, consumption, cycle share)

Green bars = net production of target metabolite, red = net consumption.
Star markers highlight KETOGLUTREDUCT-RXN and RXN-16701.
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

TARGET_METABOLITES = [
	'R-2-HYDROXYGLUTARATE[c]',
	'2-KETOGLUTARATE[c]',
	'NAD[c]',
	'NADH[c]',
	'PROTON[c]',
]

# Classification for summary table.
MET_CLASS = {
	'R-2-HYDROXYGLUTARATE[c]': 'K',
	'2-KETOGLUTARATE[c]': 'H',
	'NAD[c]': 'H',
	'NADH[c]': 'H',
	'PROTON[c]': 'N',
}

# Reactions that form the futile cycle (for highlighting).
CYCLE_RXNS = {'KETOGLUTREDUCT-RXN', 'RXN-16701'}

TOP_N = 20
IGNORE_FIRST_N_GENS = 0


class Plot(multigenAnalysisPlot.MultigenAnalysisPlot):
	def do_plot(self, seedOutDir, plotOutDir, plotOutFileName, simDataFile,
				validationDataFile, metadata):

		with open(simDataFile, 'rb') as f:
			sim_data = pickle.load(f)

		cell_density = sim_data.constants.cell_density
		reaction_stoich = sim_data.process.metabolism.reaction_stoich

		# --- Find base reactions per metabolite ---
		# Only consider forward entries (skip ' (reverse)' to avoid double-
		# counting since base_reaction_fluxes already nets fwd-rev).
		met_rxns = {}  # {met_id: {base_rxn_id: stoich_coeff}}
		for met_id in TARGET_METABOLITES:
			rxns = {}
			for rxn_id, stoich in reaction_stoich.items():
				if rxn_id.endswith(' (reverse)'):
					continue
				if met_id in stoich:
					rxns[rxn_id] = stoich[met_id]
			met_rxns[met_id] = rxns

		# --- Get cell paths ---
		cell_paths = self.ap.get_cells(
			generation=np.arange(IGNORE_FIRST_N_GENS, self.ap.n_generation),
			only_successful=True)
		if len(cell_paths) == 0:
			print('Skipping analysis -- no cells found.')
			return

		# --- Build base-reaction index lookup ---
		fba0 = TableReader(os.path.join(cell_paths[0], 'simOut', 'FBAResults'))
		base_rxn_ids = list(fba0.readAttribute('base_reaction_ids'))
		base_id_to_idx = {r: i for i, r in enumerate(base_rxn_ids)}

		# Collect all base-reaction column indexes needed.
		all_needed_idxs = set()
		met_base_info = {}  # {met_id: {rxn_id: (col_idx, stoich_coeff)}}
		for met_id, rxns in met_rxns.items():
			info = {}
			for rxn_id, coeff in rxns.items():
				bi = base_id_to_idx.get(rxn_id)
				if bi is not None:
					info[rxn_id] = (bi, coeff)
					all_needed_idxs.add(bi)
			met_base_info[met_id] = info

		all_cols = sorted(all_needed_idxs)
		col_remap = {c: i for i, c in enumerate(all_cols)}

		# --- Accumulate base_reaction_fluxes ---
		flux_accum = []
		n_timesteps = 0

		for cell_path in cell_paths:
			sim_out = os.path.join(cell_path, 'simOut')
			try:
				cell_mass = TableReader(os.path.join(sim_out, 'Mass')
					).readColumn('cellMass', squeeze=True)[1:]
				dry_mass = TableReader(os.path.join(sim_out, 'Mass')
					).readColumn('dryMass', squeeze=True)[1:]

				base_fluxes_raw = TableReader(os.path.join(sim_out, 'FBAResults')
					).readColumn('base_reaction_fluxes',
								 indices=all_cols, squeeze=False)[1:]
			except Exception as e:
				print(f'Skipped cell {sim_out}: {e!r}')
				continue

			conversion = (dry_mass / cell_mass
						  * cell_density.asNumber(MASS_UNITS / VOLUME_UNITS))
			fluxes = (
				(COUNTS_UNITS / MASS_UNITS / TIME_UNITS)
				* (base_fluxes_raw / conversion[:, np.newaxis])
			).asNumber(units.mmol / units.g / units.h)
			fluxes[np.isinf(fluxes)] = np.nan

			flux_accum.append(fluxes)
			n_timesteps += fluxes.shape[0]

		if not flux_accum:
			print('No cells processed.')
			return

		flux_all = np.vstack(flux_accum)  # (total_T, n_cols)

		# --- Compute mean contribution per reaction per metabolite ---
		# contribution = mean(net_flux) * stoich_coeff
		met_contributions = {}  # {met_id: {rxn_id: mean_contribution}}
		for met_id, info in met_base_info.items():
			contribs = {}
			for rxn_id, (bi, coeff) in info.items():
				ci = col_remap[bi]
				mean_flux = np.nanmean(flux_all[:, ci])
				contribs[rxn_id] = mean_flux * coeff
			met_contributions[met_id] = contribs

		# --- Plot ---
		fig, axes = plt.subplots(2, 2, figsize=(16, 14))

		def plot_bar(ax, met_id, top_n=None, title_extra=''):
			"""Horizontal bar chart of reactions sorted by |contribution|."""
			contribs = met_contributions.get(met_id, {})
			if not contribs:
				ax.text(0.5, 0.5, 'No data', ha='center', va='center')
				return

			sorted_rxns = sorted(contribs.keys(),
								 key=lambda r: abs(contribs[r]),
								 reverse=True)
			if top_n:
				sorted_rxns = sorted_rxns[:top_n]

			labels = []
			vals = []
			for rxn in sorted_rxns:
				labels.append(rxn)
				vals.append(contribs[rxn])

			y_pos = np.arange(len(labels))
			bar_colors = ['#2ca02c' if v > 0 else '#d62728' for v in vals]
			ax.barh(y_pos, vals, color=bar_colors, alpha=0.8)

			# Star markers on cycle reactions.
			for yi, rxn in enumerate(labels):
				if rxn in CYCLE_RXNS:
					ax.plot(vals[yi], yi, marker='*', color='gold',
							markersize=10, markeredgecolor='black',
							markeredgewidth=0.5, zorder=5)

			ax.set_yticks(y_pos)
			ax.set_yticklabels(labels, fontsize=5)
			ax.axvline(0, color='black', lw=0.5, alpha=0.3)
			n_label = f' (top {top_n})' if top_n else ''
			ax.set_title(f'{met_id}{n_label}{title_extra}', fontsize=8)
			ax.set_xlabel('Mean net flux contribution (mmol/gDCW/h)', fontsize=7)
			ax.tick_params(labelsize=6)
			ax.invert_yaxis()

		# Panel 1: R-2-HYDROXYGLUTARATE — all reactions
		plot_bar(axes[0, 0], 'R-2-HYDROXYGLUTARATE[c]')

		# Panel 2: 2-KETOGLUTARATE — top 20
		plot_bar(axes[0, 1], '2-KETOGLUTARATE[c]', top_n=TOP_N)

		# Panel 3: NAD — top 20
		plot_bar(axes[1, 0], 'NAD[c]', top_n=TOP_N)

		# Panel 4: summary text table
		ax4 = axes[1, 1]
		ax4.axis('off')

		rows = []
		header = ['Metabolite', 'Class', 'Total Prod', 'Total Cons',
				  'Cycle Contrib', 'Cycle %']
		for met_id in TARGET_METABOLITES:
			contribs = met_contributions.get(met_id, {})
			cls = MET_CLASS.get(met_id, '?')
			total_prod = sum(v for v in contribs.values() if v > 0)
			total_cons = sum(v for v in contribs.values() if v < 0)
			cycle_contrib = sum(contribs.get(rxn, 0) for rxn in CYCLE_RXNS)
			turnover = total_prod + abs(total_cons)
			pct = (abs(cycle_contrib) / (turnover / 2) * 100
				   if turnover > 0 else 0)
			rows.append([
				met_id[:25], cls,
				f'{total_prod:.3g}', f'{total_cons:.3g}',
				f'{cycle_contrib:.3g}', f'{pct:.1f}%',
			])

		table = ax4.table(cellText=rows, colLabels=header, loc='center',
						  cellLoc='center')
		table.auto_set_font_size(False)
		table.set_fontsize(6)
		table.scale(1.0, 1.4)
		# Style header.
		for j in range(len(header)):
			table[0, j].set_facecolor('#d4e6f1')
			table[0, j].set_text_props(weight='bold')
		ax4.set_title('Futile cycle turnover summary', fontsize=9, pad=15)

		fig.suptitle('KETOGLUTREDUCT-RXN pathway context — metabolite turnover',
					 fontsize=12, y=1.002)
		plt.tight_layout()
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close('all')


if __name__ == '__main__':
	Plot().cli()
