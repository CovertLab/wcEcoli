"""
Figure F (PLAN.md §6.5) — the mechanism in one lineage.

Stacked rows across induction for a single seed:

	new gene copy number
	mean rRNA operon copy number
	a terminus-proximal control gene's copy number
	number of origins (replication rounds in flight)
	cell mass
	doubling time

The point is visual and specific: as burden accumulates, the origin-proximal
rows fall further than the terminus-proximal row, and they fall in step with
the loss of overlapping replication rounds. It extends the existing cascade
figure rather than replacing it.

This is a multigen analysis, so it runs against one seed of one variant. Run it
on the highest-burden variant to see the effect; on the control it should be
flat, which is itself a useful negative.

Time is plotted as a continuous axis across generations so the division events
are visible as sawteeth rather than being averaged away.
"""

import os
import pickle

from matplotlib import pyplot as plt
import numpy as np

from models.ecoli.analysis import multigenAnalysisPlot
from wholecell.analysis.analysis_tools import (exportFigure,
	first_cell_with_table, read_stacked_columns)
from wholecell.io.tablereader import TableReader

# Terminus-proximal reference, resolved by position at run time.
TERMINUS_TARGET_FRACTION = 0.95


class Plot(multigenAnalysisPlot.MultigenAnalysisPlot):
	def do_plot(self, seedOutDir, plotOutDir, plotOutFileName, simDataFile,
			validationDataFile, metadata):
		cell_paths = self.ap.get_cells()
		if len(cell_paths) == 0:
			print('No cells found for this seed.')
			return

		with open(simDataFile, 'rb') as handle:
			sim_data = pickle.load(handle)
		cistron_data = sim_data.process.transcription.cistron_data.struct_array
		replication = sim_data.process.replication

		readable = first_cell_with_table(cell_paths, 'RnaSynthProb')
		if readable is None:
			print('No cell in this lineage has a readable RnaSynthProb.')
			return
		rsp = TableReader(os.path.join(readable, 'simOut', 'RnaSynthProb'))
		cistron_ids = rsp.readAttribute('cistron_ids')

		coord_of = dict(zip(
			cistron_data['id'], cistron_data['replication_coordinate']))
		rrna_of = dict(zip(cistron_data['id'], cistron_data['is_rRNA']))
		new_of = dict(zip(cistron_data['id'], cistron_data['is_new_gene']))

		coords = np.array(
			[coord_of.get(c, 0) for c in cistron_ids], dtype=float)
		arm = np.where(
			coords >= 0,
			replication.replichore_lengths[0],
			replication.replichore_lengths[1])
		frac = np.abs(coords) / arm

		is_rrna = np.array([bool(rrna_of.get(c, False)) for c in cistron_ids])
		is_new = np.array([bool(new_of.get(c, False)) for c in cistron_ids])
		if not is_new.any():
			print('No new gene in this simulation; this plot expects the new '
				'gene option to be enabled.')
			return

		eligible = ~(is_rrna | is_new)
		dist = np.where(
			eligible, np.abs(frac - TERMINUS_TARGET_FRACTION), np.inf)
		term_idx = int(np.argmin(dist))
		print('terminus reference gene: %s (f=%.3f)'
			% (cistron_ids[term_idx], frac[term_idx]))

		new_idx = np.where(is_new)[0]
		rrna_idx = np.where(is_rrna)[0]

		# Per-time-point series, stacked across the lineage.
		time = read_stacked_columns(
			cell_paths, 'Main', 'time', ignore_exception=True).squeeze()
		cn_new = read_stacked_columns(
			cell_paths, 'RnaSynthProb', 'gene_copy_number',
			ignore_exception=True,
			fun=lambda x: x[:, new_idx].mean(axis=1)[:, None]).squeeze()
		cn_rrna = read_stacked_columns(
			cell_paths, 'RnaSynthProb', 'gene_copy_number',
			ignore_exception=True,
			fun=lambda x: x[:, rrna_idx].mean(axis=1)[:, None]).squeeze()
		cn_term = read_stacked_columns(
			cell_paths, 'RnaSynthProb', 'gene_copy_number',
			ignore_exception=True,
			fun=lambda x: x[:, [term_idx]]).squeeze()
		n_oric = read_stacked_columns(
			cell_paths, 'ReplicationData', 'numberOfOric',
			ignore_exception=True).squeeze()
		cell_mass = read_stacked_columns(
			cell_paths, 'Mass', 'cellMass', ignore_exception=True).squeeze()

		if time.size == 0:
			print('No readable time series.')
			return

		hours = (time - time[0]) / 3600.0

		# Per-generation doubling times, drawn as a step series.
		gen_edges, gen_taus = [], []
		for path in cell_paths:
			try:
				t = TableReader(
					os.path.join(path, 'simOut', 'Main')).readColumn('time')
			except Exception:
				continue
			if len(t) < 2:
				continue
			gen_edges.append((t[0] - time[0]) / 3600.0)
			gen_taus.append((t[-1] - t[0]) / 60.0)

		self._plot(plotOutDir, plotOutFileName, metadata, hours, cn_new,
			cn_rrna, cn_term, n_oric, cell_mass, gen_edges, gen_taus,
			cistron_ids[term_idx], frac)
		self._report(cn_new, cn_rrna, cn_term, gen_taus)

	def _report(self, cn_new, cn_rrna, cn_term, gen_taus):
		print('')
		print('=' * 62)
		print('LINEAGE TRACE (PLAN.md §6.5)')
		print('=' * 62)
		if len(gen_taus) >= 2:
			print('  doubling time: %.1f -> %.1f min over %d generations'
				% (gen_taus[0], gen_taus[-1], len(gen_taus)))
		# Compare the first and last tenth of the trace.
		def ends(v):
			k = max(1, len(v) // 10)
			return float(np.mean(v[:k])), float(np.mean(v[-k:]))
		for name, series in (('new gene', cn_new), ('rRNA operons', cn_rrna),
				('terminus ref', cn_term)):
			if series.size:
				a, b = ends(series)
				print('  %-13s copy number %.3f -> %.3f  (lost %.3f)'
					% (name, a, b, a - b))
		print('=' * 62)

	def _plot(self, plotOutDir, plotOutFileName, metadata, hours, cn_new,
			cn_rrna, cn_term, n_oric, cell_mass, gen_edges, gen_taus,
			term_id, frac):
		panels = [
			('new gene\ncopy number', cn_new, 'C0'),
			('rRNA operons\nmean copy number', cn_rrna, 'C1'),
			('%s (terC)\ncopy number' % term_id, cn_term, 'C2'),
			('origins\nper cell', n_oric, 'C3'),
			('cell mass\n(fg)', cell_mass, 'C4'),
			]
		fig, axes = plt.subplots(len(panels) + 1, 1, figsize=(9, 11),
			sharex=True)

		for ax, (label, series, colour) in zip(axes, panels):
			if series is None or np.size(series) == 0:
				ax.text(0.5, 0.5, 'no data', ha='center', va='center',
					transform=ax.transAxes)
			else:
				ax.plot(hours[:len(series)], series, color=colour, lw=0.8)
			ax.set_ylabel(label, fontsize='x-small')
			ax.spines['top'].set_visible(False)
			ax.spines['right'].set_visible(False)

		ax = axes[-1]
		if gen_edges:
			ax.step(gen_edges, gen_taus, where='post', color='C5', lw=1.4)
		ax.set_ylabel('doubling\ntime (min)', fontsize='x-small')
		ax.set_xlabel('time (h)')
		ax.spines['top'].set_visible(False)
		ax.spines['right'].set_visible(False)

		axes[0].set_title(
			'Dosage cascade through induction — origin-proximal rows fall '
			'furthest', fontsize='medium')

		plt.tight_layout()
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close('all')


if __name__ == '__main__':
	Plot().cli()
