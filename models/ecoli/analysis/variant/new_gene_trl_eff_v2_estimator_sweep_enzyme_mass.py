"""
Enzyme-reduction and mass-category evidence for the
new_gene_trl_eff_v2_estimator_sweep.

This complements the growth_limitation analysis (which covers ribosome/RNAP
machinery and supply signals).  Here we confirm the premise of the kcat
experiment -- that GFP ribosome competition reduces metabolic catalyst levels --
and show to what extent that should tighten the kcat bounds (bound = kcat*[E],
so the ceiling scales directly with catalyst concentration).

Across the trl-eff gradient, averaged over the settled window (last 8 gens, all
seeds), it plots (2x3):

  1 dry-mass category fractions (protein/RNA/DNA/small mol), gfp_only block
  2 RNA sub-fractions (rRNA/tRNA/mRNA), gfp_only block -- rRNA tracks ribosomes
  3 mean metabolic-catalyst concentration, normalized to wildtype, per block
      -> the direct "enzymes are reduced" / "bound is reduced" curve
  4 mean metabolic-catalyst concentration (absolute), per block
  5 mean reaction flux of the kcat-constrained pairs, per block -> demand
  6 normalized catalyst conc vs normalized flux overlaid (gfp_only) -> if they
      fall in parallel, the bound (proportional to [E]) tracks demand and
      saturation stays ~constant (bounds do not trigger harder under stress)

Catalyst concentration uses the same (reaction, catalyst) pairs and
conc = catalyst_counts * countsToMolar that define the v2 kcat bounds, so it is
exactly the [E] that sets each bound.  Reads only recorded listeners -- no
re-simulation.
"""

import os
import pickle

import numpy as np
from matplotlib import pyplot as plt

from models.ecoli.analysis import variantAnalysisPlot
from models.ecoli.sim.variants.new_gene_trl_eff_v2_estimator_sweep import (
	BLOCK_ESTIMATORS,
	BLOCK_SHORT_NAMES,
	TRL_EFF_VALUES,
	N_PER_BLOCK,
	is_wildtype,
	variant_to_block,
)
from models.ecoli.analysis.cohort.kcat_estimates_v2_with_estimators import (
	_read_csv_ids,
	below_line_essential_monomer_ids_filepath,
	below_line_essential_complex_ids_filepath,
)
from models.ecoli.processes.metabolism import (
	COUNTS_UNITS, VOLUME_UNITS, TIME_UNITS, MASS_UNITS)
from wholecell.analysis.analysis_tools import exportFigure, read_stacked_columns
from wholecell.io.tablereader import TableReader
from wholecell.utils import units


BLOCK_COLORS = ('#7f7f7f', '#1f77b4', '#2ca02c')  # gfp_only, gap, smoothed_max
WILDTYPE_COLOR = '#222222'
CATEGORY_COLORS = {
	'protein': '#1f77b4', 'RNA': '#2ca02c', 'DNA': '#d62728',
	'small mol': '#ff7f0e', 'rRNA': '#2ca02c', 'tRNA': '#9467bd',
	'mRNA': '#8c564b'}


def _series_by_block(metric, variant_indexes):
	"""Return {block_idx: list-of-value-by-trl_eff-index}."""
	by_block = {b: [np.nan] * N_PER_BLOCK
		for b in range(len(BLOCK_ESTIMATORS))}
	for vi in variant_indexes:
		if is_wildtype(vi):
			continue
		block_idx, te_idx = variant_to_block(vi)
		by_block[block_idx][te_idx] = metric.get(vi, np.nan)
	return by_block


def _col(cells, table, column):
	"""Stacked column over the window (remove_first), or None if unreadable."""
	if len(cells) == 0:
		return None
	try:
		return read_stacked_columns(
			cells, table, column, remove_first=True, ignore_exception=True)
	except Exception:
		return None


def _build_pair_indices(sim_data, example_cell):
	"""Map the below-line-essential (reaction, catalyst) pairs to FBAResults
	column indices -- the same pairs the v2 kcat bounds constrain."""
	essential = set(
		_read_csv_ids(below_line_essential_monomer_ids_filepath)
		+ _read_csv_ids(below_line_essential_complex_ids_filepath))
	fba = TableReader(os.path.join(example_cell, 'simOut', 'FBAResults'))
	listener_rxn_ids = fba.readAttribute('reactionIDs')
	listener_cat_ids = fba.readAttribute('catalyst_ids')
	r2c = sim_data.process.metabolism.reaction_catalysts
	c2r = {}
	for rxn, cats in r2c.items():
		for cat in cats:
			c2r.setdefault(cat, []).append(rxn)
	ridx = {r: i for i, r in enumerate(listener_rxn_ids)}
	cidx = {c: i for i, c in enumerate(listener_cat_ids)}

	# Drop reactions with extra non-essential catalysts (flux not cleanly
	# attributable), matching the v2 estimator pipeline.
	filter_out = set()
	for rxn in {r for cat in (essential & set(listener_cat_ids))
			for r in c2r.get(cat, [])}:
		all_cats = r2c.get(rxn, [])
		ess_n = sum(1 for c in all_cats if c in essential)
		if len(all_cats) > 1 and ess_n < len(all_cats):
			filter_out.add(rxn)

	rxn_indexes, cat_indexes = [], []
	for cat in (essential & set(listener_cat_ids)):
		for rxn in c2r.get(cat, []):
			if rxn in filter_out or rxn not in ridx or cat not in cidx:
				continue
			rxn_indexes.append(ridx[rxn])
			cat_indexes.append(cidx[cat])
	return np.array(rxn_indexes), np.array(cat_indexes)


class Plot(variantAnalysisPlot.VariantAnalysisPlot):
	def do_plot(self, inputDir, plotOutDir, plotOutFileName, simDataFile,
				validationDataFile, metadata):
		variant_indexes = self.ap.get_variants()
		n_total_gens = self.ap.n_generation
		ignore_first_n_gens = max(n_total_gens - 8, 0)
		window = np.arange(ignore_first_n_gens, n_total_gens)

		with open(simDataFile, 'rb') as f:
			sim_data = pickle.load(f)

		all_cells = self.ap.get_cells(only_successful=True)
		if len(all_cells) == 0:
			print('enzyme_mass: no successful cells; skipping.')
			return
		rxn_indexes, cat_indexes = _build_pair_indices(sim_data, all_cells[0])
		cell_density_value = sim_data.constants.cell_density.asNumber(
			MASS_UNITS / VOLUME_UNITS)

		# Per-variant means.
		frac = {c: {} for c in ('protein', 'RNA', 'DNA', 'small mol',
			'rRNA', 'tRNA', 'mRNA')}
		conc, flux = {}, {}
		for vi in variant_indexes:
			cells = self.ap.get_cells(
				variant=[vi], generation=window, only_successful=True)
			dry = _col(cells, 'Mass', 'dryMass')
			if dry is not None and np.nanmean(dry) > 0:
				dm = float(np.nanmean(dry))
				for key, col in (
						('protein', 'proteinMass'), ('RNA', 'rnaMass'),
						('DNA', 'dnaMass'), ('small mol', 'smallMoleculeMass'),
						('rRNA', 'rRnaMass'), ('tRNA', 'tRnaMass'),
						('mRNA', 'mRnaMass')):
					arr = _col(cells, 'Mass', col)
					frac[key][vi] = (float(np.nanmean(arr)) / dm
						if arr is not None else np.nan)
			else:
				for key in frac:
					frac[key][vi] = np.nan

			# Catalyst concentration = catalyst_counts * countsToMolar -- the
			# [E] that sets each kcat bound.
			conc[vi] = np.nan
			flux[vi] = np.nan
			if len(cat_indexes) > 0:
				cc = _col(cells, 'FBAResults', 'catalyst_counts')
				c2m = _col(cells, 'EnzymeKinetics', 'countsToMolar')
				if cc is not None and c2m is not None:
					c2m1 = c2m.reshape(c2m.shape[0], -1)[:, :1]
					conc_ts = cc[:, cat_indexes] * c2m1
					conc[vi] = float(np.nanmean(conc_ts))

				rf = _col(cells, 'FBAResults', 'reactionFluxes')
				cm = _col(cells, 'Mass', 'cellMass')
				if rf is not None and dry is not None and cm is not None:
					conv = (dry.reshape(-1) / cm.reshape(-1)
						* cell_density_value)
					with np.errstate(divide='ignore', invalid='ignore'):
						f = ((COUNTS_UNITS / MASS_UNITS / TIME_UNITS)
							* (rf[:, rxn_indexes] / conv[:, np.newaxis])
							).asNumber(units.mmol / units.g / units.h)
					f[~np.isfinite(f)] = np.nan
					flux[vi] = float(np.nanmean(f))

		wt = next((vi for vi in variant_indexes if is_wildtype(vi)), None)

		def _norm(metric):
			base = metric.get(wt) if wt is not None else None
			if base is None or not np.isfinite(base) or base == 0:
				return {vi: np.nan for vi in metric}
			return {vi: metric[vi] / base for vi in metric}

		conc_norm = _norm(conc)
		flux_norm = _norm(flux)
		x = np.array(TRL_EFF_VALUES)

		fig, axes = plt.subplots(2, 3, figsize=(17, 9))

		def _categories(ax, keys, title, ylabel):
			for key in keys:
				series = _series_by_block(frac[key], variant_indexes)[0]
				ax.plot(x, series, marker='o', lw=1.5,
					color=CATEGORY_COLORS[key], label=key)
			ax.set_xlabel('GFP translation efficiency')
			ax.set_ylabel(ylabel)
			ax.set_title(title)
			ax.legend(fontsize=7)

		def _by_block(ax, metric, title, ylabel, wt_fmt=None):
			for b in range(len(BLOCK_ESTIMATORS)):
				ax.plot(x, _series_by_block(metric, variant_indexes)[b],
					marker='o', lw=1.5, color=BLOCK_COLORS[b],
					label=BLOCK_SHORT_NAMES[b])
			if wt is not None and np.isfinite(metric.get(wt, np.nan)):
				lbl = ('wildtype' if wt_fmt is None
					else f'wildtype ({wt_fmt.format(metric[wt])})')
				ax.axhline(metric[wt], color=WILDTYPE_COLOR, ls='--', lw=1.0,
					label=lbl)
			ax.set_xlabel('GFP translation efficiency')
			ax.set_ylabel(ylabel)
			ax.set_title(title)
			ax.legend(fontsize=7)

		_categories(axes[0, 0], ('protein', 'RNA', 'DNA', 'small mol'),
			'Dry-mass category fractions (gfp_only)', 'fraction of dry mass')
		_categories(axes[0, 1], ('rRNA', 'tRNA', 'mRNA'),
			'RNA sub-fractions (gfp_only) -- rRNA tracks ribosomes',
			'fraction of dry mass')
		_by_block(axes[0, 2], conc_norm,
			'Metabolic catalyst conc / wildtype\n'
			'(= kcat-bound ceiling / wildtype)', 'fraction of wildtype')

		_by_block(axes[1, 0], conc,
			'Mean metabolic catalyst conc (absolute)',
			'conc (counts x countsToMolar)')
		_by_block(axes[1, 1], flux,
			'Mean flux of kcat-constrained reactions',
			'flux (mmol/g/h)')

		# Panel 6: normalized conc vs flux (gfp_only) -- parallel => bound
		# tracks demand, saturation ~constant.
		ax = axes[1, 2]
		conc_n_series = _series_by_block(conc_norm, variant_indexes)[0]
		flux_n_series = _series_by_block(flux_norm, variant_indexes)[0]
		ax.plot(x, conc_n_series, marker='o', lw=1.6, color='#1f77b4',
			label='catalyst conc (= bound)')
		ax.plot(x, flux_n_series, marker='s', lw=1.6, color='#d62728',
			label='reaction flux (= demand)')
		ax.set_xlabel('GFP translation efficiency')
		ax.set_ylabel('fraction of wildtype (gfp_only)')
		ax.set_title('Bound vs demand (gfp_only)\n'
			'parallel -> saturation ~constant, bound not triggered harder')
		ax.legend(fontsize=7)

		fig.suptitle(
			'new_gene_trl_eff_v2_estimator_sweep: GFP reduces metabolic '
			'enzyme levels (and the kcat ceiling) -- but demand falls with it\n'
			f'(averaged over gens {ignore_first_n_gens}-{n_total_gens - 1}, '
			'all seeds)',
			fontsize=13)
		fig.tight_layout(rect=[0, 0, 1, 0.92])
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close('all')


if __name__ == '__main__':
	Plot().cli()
