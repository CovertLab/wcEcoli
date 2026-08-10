"""
One row per variant: the basic statistics every other figure and conversation
needs, so nobody has to re-derive them from a specialised analysis.

Why this exists
---------------
Doubling time, construct mRNA and protein counts, the construct's share of the
proteome, cell mass and the machinery pools are quoted constantly and have so
far been pulled out of whichever analysis happened to compute them as a side
effect. That is how the same quantity ends up with two values in two reports.
This script is the single source for them, it runs in a couple of minutes, and
it is meant to be run on every batch.

It is deliberately descriptive. Nothing here is a decomposition or a test; it
is the table you put beside a figure.

What this measures
------------------
Per variant, averaged over the analysis window and reduced per cell so that
standard errors are across cells rather than across autocorrelated timesteps:

  tau                    doubling time, birth to division, minutes
  growth_rate            ln(2) / tau, per hour
  dry_mass               Mass/dryMass, fg
  protein_mass           Mass/proteinMass, fg
  rna_mass               Mass/rnaMass, fg
  gfp_mrna               construct mRNA count, RNACounts/mRNA_counts
  gfp_protein            construct monomer count, MonomerCounts
  gfp_mrna_frac          construct share of all mRNA counts
  gfp_protein_frac       construct share of all monomer counts
  gfp_proteome_mass_pct  construct monomer mass / total protein mass, percent
  active_rnap            UniqueMoleculeCounts
  active_ribosome        UniqueMoleculeCounts
  n_cells                cells contributing to the row

The proteome percentage is mass-based, not count-based, because that is the
quantity comparable to a measured proteomics fraction: monomer count times
monomer molecular weight, over Mass/proteinMass. Both are reported since the
count fraction is what heatmap analyses in this repo use.

Emits a CSV alongside the PDF. The CSV is the point; the PDF is a glance.
"""

import csv
import os
import pickle

from matplotlib import pyplot as plt
import numpy as np

from models.ecoli.analysis import variantAnalysisPlot
from wholecell.analysis.analysis_tools import (exportFigure,
	first_cell_with_table, read_stacked_columns)
from wholecell.io.tablereader import TableReader
from wholecell.utils import units

# Matches every other analysis in this study.
IGNORE_FIRST_N_GENS = 16

AVOGADRO = 6.02214076e23
FG_PER_G = 1e15


def _window(n_generation):
	"""Generations to analyse: everything past the burn-in."""
	if n_generation <= IGNORE_FIRST_N_GENS:
		return None
	return np.arange(IGNORE_FIRST_N_GENS, n_generation)


def _sem(values):
	"""Standard error of the mean across cells. NaN if fewer than two."""
	v = np.asarray(values, dtype=float).ravel()
	v = v[np.isfinite(v)]
	if v.size < 2:
		return float('nan')
	return float(np.std(v, ddof=1) / np.sqrt(v.size))


def _tau_minutes(x):
	"""Reduce one cell's time column to its doubling time in minutes."""
	return np.array([[(x[-1, 0] - x[0, 0]) / 60.0]])


def _cell_mean(x):
	"""Per cell: mean over that cell's timesteps, one row per cell."""
	return np.array([[float(np.mean(x))]])


def _sum_over(idx):
	"""Per cell: sum the given subcolumns, then average over time."""
	def fn(x):
		if idx is None or len(idx) == 0:
			return np.array([[0.0]])
		return np.array([[float(np.mean(x[:, idx].sum(axis=1)))]])
	return fn


def _sum_all(x):
	"""Per cell: sum across all subcolumns, then average over time."""
	return np.array([[float(np.mean(x.sum(axis=1)))]])


class Plot(variantAnalysisPlot.VariantAnalysisPlot):
	def do_plot(self, inputDir, plotOutDir, plotOutFileName, simDataFile,
			validationDataFile, metadata):
		variants = sorted(self.ap.get_variants())
		if not variants:
			print('No variants found.')
			return

		generations = _window(self.ap.n_generation)
		if generations is None:
			print('Run has only %d generations, fewer than the %d-generation '
				'burn-in; using all of them.'
				% (self.ap.n_generation, IGNORE_FIRST_N_GENS))
			generations = np.arange(self.ap.n_generation)

		with open(simDataFile, 'rb') as handle:
			sim_data = pickle.load(handle)
		construct = self._construct_ids(sim_data)

		rows = []
		for variant in variants:
			row = self._measure(variant, generations, construct)
			if row is None:
				print('No usable cells for variant %d; skipping.' % variant)
				continue
			row['variant'] = variant
			rows.append(row)

		if not rows:
			print('Nothing measurable.')
			return

		self._write_csv(plotOutDir, plotOutFileName, rows)
		self._report(rows)
		self._plot(plotOutDir, plotOutFileName, rows, metadata)

	def _construct_ids(self, sim_data):
		"""
		The construct's mRNA (cistron) ids, monomer ids and monomer molecular
		weights, or empty lists if there is no new gene in this sim_data.
		"""
		transcription = sim_data.process.transcription
		cistron_data = transcription.cistron_data.struct_array
		monomer_data = sim_data.process.translation.monomer_data.struct_array

		cistrons = sorted(set(
			cistron_data[cistron_data['is_new_gene']]['id'].tolist()))
		if not cistrons:
			print('No new gene in this sim_data; construct columns will be 0.')
			return dict(cistrons=[], monomers=[], mw={})

		mapping = dict(zip(monomer_data['cistron_id'], monomer_data['id']))
		monomers = [str(mapping[c]) for c in cistrons if c in mapping]

		mw_units = sim_data.process.translation.monomer_data['mw']
		mws = np.asarray(mw_units.asNumber(units.g / units.mol), dtype=float)
		ids = [str(i) for i in monomer_data['id']]
		mw = {m: float(mws[ids.index(m)]) for m in monomers if m in ids}

		print('Construct cistrons: %s' % cistrons)
		print('Construct monomers: %s' % monomers)
		return dict(cistrons=cistrons, monomers=monomers, mw=mw)

	def _measure(self, variant, generations, construct):
		"""Everything for one variant."""
		cell_paths = self.ap.get_cells(
			variant=[variant], generation=generations)
		if len(cell_paths) == 0:
			return None

		main_cell = first_cell_with_table(cell_paths, 'Main')
		if main_cell is None:
			return None

		taus = read_stacked_columns(cell_paths, 'Main', 'time',
			ignore_exception=True, fun=_tau_minutes)
		if taus.size == 0:
			return None
		tau = float(np.mean(taus))
		row = dict(tau=tau, tau_sem=_sem(taus), n_cells=int(taus.size),
			growth_rate=float(np.log(2) / (tau / 60.0)) if tau else float('nan'))

		for label, column in (('dry_mass', 'dryMass'),
				('protein_mass', 'proteinMass'), ('rna_mass', 'rnaMass')):
			values = read_stacked_columns(cell_paths, 'Mass', column,
				ignore_exception=True, fun=_cell_mean)
			row[label] = float(np.mean(values)) if values.size else float('nan')
			row['%s_sem' % label] = _sem(values)

		# ---- construct mRNA ----------------------------------------------
		row['gfp_mrna'] = 0.0
		row['gfp_mrna_frac'] = 0.0
		rna_cell = first_cell_with_table(cell_paths, 'RNACounts')
		if rna_cell and construct['cistrons']:
			reader = TableReader(os.path.join(rna_cell, 'simOut', 'RNACounts'))
			names = list(reader.readAttribute('mRNA_cistron_ids'))
			idx = np.array([names.index(c) for c in construct['cistrons']
				if c in names], dtype=int)
			if idx.size:
				counts = read_stacked_columns(cell_paths, 'RNACounts',
					'mRNA_cistron_counts', ignore_exception=True,
					fun=_sum_over(idx))
				total = read_stacked_columns(cell_paths, 'RNACounts',
					'mRNA_cistron_counts', ignore_exception=True,
					fun=_sum_all)
				if counts.size:
					row['gfp_mrna'] = float(np.mean(counts))
					row['gfp_mrna_sem'] = _sem(counts)
				if total.size and float(np.mean(total)):
					row['gfp_mrna_frac'] = (
						float(np.mean(counts)) / float(np.mean(total)))

		# ---- construct protein -------------------------------------------
		row['gfp_protein'] = 0.0
		row['gfp_protein_frac'] = 0.0
		row['gfp_proteome_mass_pct'] = 0.0
		mono_cell = first_cell_with_table(cell_paths, 'MonomerCounts')
		if mono_cell and construct['monomers']:
			reader = TableReader(os.path.join(
				mono_cell, 'simOut', 'MonomerCounts'))
			names = list(reader.readAttribute('monomerIds'))
			idx = np.array([names.index(m) for m in construct['monomers']
				if m in names], dtype=int)
			if idx.size:
				counts = read_stacked_columns(cell_paths, 'MonomerCounts',
					'monomerCounts', ignore_exception=True, fun=_sum_over(idx))
				total = read_stacked_columns(cell_paths, 'MonomerCounts',
					'monomerCounts', ignore_exception=True, fun=_sum_all)
				if counts.size:
					row['gfp_protein'] = float(np.mean(counts))
					row['gfp_protein_sem'] = _sem(counts)
				if total.size and float(np.mean(total)):
					row['gfp_protein_frac'] = (
						float(np.mean(counts)) / float(np.mean(total)))

				# Mass-based proteome share, the quantity comparable to a
				# measured proteomics fraction. Counts times molecular weight
				# over the cell's total protein mass.
				mw = np.mean([construct['mw'][m] for m in construct['monomers']
					if m in construct['mw']]) if construct['mw'] else 0.0
				grams = float(np.mean(counts)) * mw / AVOGADRO
				fg = grams * FG_PER_G
				if row.get('protein_mass'):
					row['gfp_proteome_mass_pct'] = 100.0 * fg / row['protein_mass']

		# ---- machinery pools ---------------------------------------------
		umc_cell = first_cell_with_table(cell_paths, 'UniqueMoleculeCounts')
		row['active_rnap'] = float('nan')
		row['active_ribosome'] = float('nan')
		if umc_cell:
			reader = TableReader(os.path.join(
				umc_cell, 'simOut', 'UniqueMoleculeCounts'))
			unique_ids = reader.readAttribute('uniqueMoleculeIds')
			counts = read_stacked_columns(cell_paths, 'UniqueMoleculeCounts',
				'uniqueMoleculeCounts', ignore_exception=True)
			for label, key in (('active_rnap', 'active_RNAP'),
					('active_ribosome', 'active_ribosome')):
				if key in unique_ids and counts.size:
					row[label] = float(
						np.mean(counts[:, unique_ids.index(key)]))

		return row

	FIELDS = ['variant', 'n_cells', 'tau', 'tau_sem', 'growth_rate',
		'dry_mass', 'dry_mass_sem', 'protein_mass', 'protein_mass_sem',
		'rna_mass', 'rna_mass_sem',
		'gfp_mrna', 'gfp_mrna_sem', 'gfp_mrna_frac',
		'gfp_protein', 'gfp_protein_sem', 'gfp_protein_frac',
		'gfp_proteome_mass_pct', 'active_rnap', 'active_ribosome']

	def _write_csv(self, plot_out_dir, plot_out_file_name, rows):
		path = os.path.join(plot_out_dir, plot_out_file_name + '.csv')
		with open(path, 'w', newline='') as handle:
			writer = csv.DictWriter(handle, fieldnames=self.FIELDS)
			writer.writeheader()
			for row in rows:
				writer.writerow({k: row.get(k, '') for k in self.FIELDS})

	def _report(self, rows):
		print('\nPer-variant summary')
		print('  %-4s %6s %7s %9s %11s %11s %9s %9s'
			% ('var', 'cells', 'tau', 'dry mass', 'GFP mRNA', 'GFP protein',
				'proteome', 'act RNAP'))
		for r in rows:
			print('  %-4d %6d %7.1f %9.1f %11.1f %11.4g %8.2f%% %9.0f'
				% (r['variant'], r['n_cells'], r['tau'], r['dry_mass'],
					r['gfp_mrna'], r['gfp_protein'],
					r['gfp_proteome_mass_pct'], r['active_rnap']))
		print()

	def _plot(self, plot_out_dir, plot_out_file_name, rows, metadata):
		variants = [r['variant'] for r in rows]
		fig, axes = plt.subplots(1, 4, figsize=(19, 4))

		axes[0].errorbar(variants, [r['tau'] for r in rows],
			yerr=[r['tau_sem'] for r in rows], marker='o')
		axes[0].set_ylabel('doubling time (min)')

		axes[1].plot(variants, [r['gfp_mrna'] for r in rows], marker='o')
		axes[1].set_ylabel('construct mRNA count')

		axes[2].plot(variants, [r['gfp_protein'] for r in rows], marker='o')
		axes[2].set_ylabel('construct monomer count')

		axes[3].plot(variants, [r['gfp_proteome_mass_pct'] for r in rows],
			marker='o')
		axes[3].set_ylabel('construct % of proteome mass')

		for ax in axes:
			ax.set_xlabel('variant')
		plt.tight_layout()
		exportFigure(plt, plot_out_dir, plot_out_file_name, metadata)
		plt.close('all')
