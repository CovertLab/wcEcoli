"""
Figure B (PLAN.md §6.2) — the honesty and validation panel.

Mean copy number against growth rate, one series each for the new gene, the
rRNA operon mean, and a terminus-proximal control gene, with
`get_average_copy_number` overlaid at each variant's own measured doubling time.
A second panel reports the residual between simulated and analytical.

The point of this figure is to be checkable. If the simulation departs from
Cooper-Helmstetter, every downstream dosage claim inherits that error, so the
residual belongs on the page rather than in a footnote. On Batch 1 the median
residual is 4.6-6.0% across variants, with no systematic trend in growth rate.

Single-batch analysis. Reductions happen per cell (see `fun=` below) to keep
memory flat.
"""

import csv
import os
import pickle

from matplotlib import pyplot as plt
import numpy as np

from models.ecoli.analysis import variantAnalysisPlot
from models.ecoli.analysis.variant.dosage_channel_decomposition import _window
from wholecell.analysis.analysis_tools import (exportFigure,
	first_cell_with_table, read_stacked_columns)
from wholecell.io.tablereader import TableReader

# A terminus-proximal reference gene. Chosen for position, not function: it is
# the control for "does this gene move because of burden, or because of where
# it sits". Resolved by replichore fraction at run time rather than hard-coded
# by id, so it stays valid if the genome changes.
TERMINUS_TARGET_FRACTION = 0.95


def _cell_mean(x):
	return x.mean(axis=0)[None, :]


def _tau_minutes(x):
	return np.array([[(x[-1, 0] - x[0, 0]) / 60.0]])


class Plot(variantAnalysisPlot.VariantAnalysisPlot):
	def do_plot(self, inputDir, plotOutDir, plotOutFileName, simDataFile,
			validationDataFile, metadata):
		variants = sorted(self.ap.get_variants())
		if len(variants) < 2:
			print('Need at least two variants. Found %d.' % len(variants))
			return

		generations = _window(self.ap.n_generation)
		if generations is None:
			generations = np.arange(self.ap.n_generation)

		with open(simDataFile, 'rb') as handle:
			sim_data = pickle.load(handle)
		cistron_data = sim_data.process.transcription.cistron_data.struct_array
		replication = sim_data.process.replication

		coord_of = dict(zip(
			cistron_data['id'], cistron_data['replication_coordinate']))
		rrna_of = dict(zip(cistron_data['id'], cistron_data['is_rRNA']))
		new_of = dict(zip(cistron_data['id'], cistron_data['is_new_gene']))

		rows = []
		series = {}
		cistron_ids = None
		selections = None

		for variant in variants:
			cell_paths = self.ap.get_cells(
				variant=[variant], generation=generations)
			if len(cell_paths) == 0:
				continue

			if cistron_ids is None:
				readable = first_cell_with_table(cell_paths, 'RnaSynthProb')
				if readable is None:
					print('No readable RnaSynthProb for variant %d.' % variant)
					continue
				rsp = TableReader(os.path.join(
					readable, 'simOut', 'RnaSynthProb'))
				cistron_ids = rsp.readAttribute('cistron_ids')
				coords = np.array(
					[coord_of.get(c, 0) for c in cistron_ids], dtype=float)
				arm = np.where(
					coords >= 0,
					replication.replichore_lengths[0],
					replication.replichore_lengths[1])
				frac = np.abs(coords) / arm
				is_rrna = np.array(
					[bool(rrna_of.get(c, False)) for c in cistron_ids])
				is_new = np.array(
					[bool(new_of.get(c, False)) for c in cistron_ids])
				# Terminus reference: the gene closest to the target fraction
				# that is neither rRNA nor the construct.
				eligible = ~(is_rrna | is_new)
				dist = np.where(
					eligible, np.abs(frac - TERMINUS_TARGET_FRACTION), np.inf)
				term_idx = int(np.argmin(dist))
				is_term = np.zeros_like(is_rrna)
				is_term[term_idx] = True
				selections = [
					('new_gene', is_new),
					('rrna_operons', is_rrna),
					('terminus_control', is_term),
					]
				print('terminus reference gene: %s (f=%.3f)'
					% (cistron_ids[term_idx], frac[term_idx]))

			copy_number = read_stacked_columns(
				cell_paths, 'RnaSynthProb', 'gene_copy_number',
				ignore_exception=True, fun=_cell_mean)
			taus = read_stacked_columns(
				cell_paths, 'Main', 'time',
				ignore_exception=True, fun=_tau_minutes)
			if copy_number.size == 0 or taus.size == 0:
				continue

			mean_cn = copy_number.mean(axis=0)
			tau = float(np.mean(taus))
			growth_rate = 60.0 * np.log(2) / tau  # per hour

			for name, sel in selections:
				if not sel.any():
					continue
				sim = float(mean_cn[sel].mean())
				analytical = float(replication.get_average_copy_number(
					tau, coords[sel]).mean())
				rows.append(dict(
					variant=variant, tau_min=tau, growth_rate_per_h=growth_rate,
					gene=name, coord=float(coords[sel].mean()),
					replichore_fraction=float(frac[sel].mean()),
					copy_number_sim=sim, copy_number_ch=analytical,
					residual_pct=100.0 * (sim - analytical) / analytical,
					n_cells=copy_number.shape[0],
					))
				series.setdefault(name, []).append(
					(growth_rate, sim, analytical))

		if not rows:
			print('No usable data.')
			return

		self._write_csv(plotOutDir, plotOutFileName, rows)
		self._report(rows)
		self._plot(plotOutDir, plotOutFileName, metadata, series)

	def _write_csv(self, plotOutDir, plotOutFileName, rows):
		fields = ['variant', 'tau_min', 'growth_rate_per_h', 'gene', 'coord',
			'replichore_fraction', 'copy_number_sim', 'copy_number_ch',
			'residual_pct', 'n_cells']
		path = os.path.join(plotOutDir, plotOutFileName + '.csv')
		with open(path, 'w', newline='') as handle:
			writer = csv.DictWriter(handle, fieldnames=fields)
			writer.writeheader()
			for row in rows:
				writer.writerow(row)
		print('Wrote %s' % path)

	def _report(self, rows):
		print('')
		print('=' * 62)
		print('COPY NUMBER vs GROWTH RATE (PLAN.md §6.2)')
		print('=' * 62)
		resid = np.array([abs(r['residual_pct']) for r in rows])
		print('  |simulated - analytical| across all points: '
			'median %.2f%%, max %.2f%%' % (np.median(resid), resid.max()))
		for name in sorted({r['gene'] for r in rows}):
			sub = [r for r in rows if r['gene'] == name]
			sub.sort(key=lambda r: r['tau_min'])
			print('  %-17s f=%.3f  cn %.3f -> %.3f across tau %.0f -> %.0f min'
				% (name, sub[0]['replichore_fraction'],
					sub[0]['copy_number_sim'], sub[-1]['copy_number_sim'],
					sub[0]['tau_min'], sub[-1]['tau_min']))
		print('=' * 62)

	def _plot(self, plotOutDir, plotOutFileName, metadata, series):
		fig, axes = plt.subplots(2, 1, figsize=(6.5, 6.6), sharex=True,
			gridspec_kw={'height_ratios': [2.2, 1]})
		ax, ax_r = axes

		for i, (name, points) in enumerate(sorted(series.items())):
			points.sort()
			g = np.array([p[0] for p in points])
			sim = np.array([p[1] for p in points])
			ch = np.array([p[2] for p in points])
			colour = 'C%d' % i
			ax.plot(g, sim, 'o-', color=colour, lw=1.8, ms=5, label=name)
			ax.plot(g, ch, '--', color=colour, lw=1.1, alpha=0.65)
			ax_r.plot(g, 100 * (sim - ch) / ch, 'o-', color=colour,
				lw=1.4, ms=4)

		ax.set_ylabel('mean copy number')
		ax.set_title('Copy number vs growth rate (dashed = Cooper-Helmstetter)',
			fontsize='medium')
		ax.legend(fontsize='x-small', frameon=False)
		ax_r.axhline(0, color='k', lw=0.8)
		ax_r.set_ylabel('residual (%)')
		ax_r.set_xlabel('growth rate (doublings / h)')
		for a in axes:
			a.spines['top'].set_visible(False)
			a.spines['right'].set_visible(False)

		plt.tight_layout()
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close('all')


if __name__ == '__main__':
	Plot().cli()
