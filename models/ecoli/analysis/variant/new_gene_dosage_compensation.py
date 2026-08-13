"""
Measure how completely the dosage normaliser compensates gene copy-number loss.

Why this exists
---------------
`transcription.synth_prob_from_ppgpp` divides every gene's synthesis probability
by the copy number that gene is expected to have:

	growth = fit(ppGpp)          <- calibrated on five wildtype conditions
	tau    = ln(2) / growth
	n_avg  = get_average_copy_number(tau, wt_replication_coordinate)
	factor = (growth + deg_rate) / n_avg

The share a transcription unit ends up with is proportional to
`a * expression * loss` normalised over all promoters, where

	a = n_actual / n_avg

so if the normaliser's expectation matched reality, `a` would be 1 for every
gene, copy-number loss would cancel, and there would be no dosage feedback at
all. Everything therefore turns on one question, which was previously answered
by assumption in both directions:

	**Is that tau the cell's realised doubling time, or not?**

It is not. It is inferred from the current ppGpp concentration, and under
heterologous burden the ppGpp sensor under-reports the slowdown. The
uncompensated remainder is a real copy-number feedback, and because `a` depends
on the gene's replichore fraction the remainder does NOT cancel through the
normalisation -- origin-proximal genes lose more of it, which is exactly the
asymmetry the rRNA operons and the RNAP subunit genes sit in.

Nothing here needs a model change or a re-run. `GrowthLimits/ppgpp_conc` is
already logged in the same units `synth_prob_from_ppgpp` consumes, so the
inferred growth rate, the inferred tau, `n_avg`, `factor` and `a` are all
reconstructible from data already on disk.

What this measures
------------------
Per variant, reduced per cell so that standard errors are across cells:

  tau_real        measured doubling time, birth to division
  ppgpp_conc      mean `GrowthLimits/ppgpp_conc`
  tau_inferred    ln(2) / fit(ppgpp_conc), the tau the normaliser actually used
  sees_frac       (tau_inferred span) / (tau_real span), 1.0 = a faithful
                  sensor and therefore complete compensation

and per gene class:

  n_actual        measured `RnaSynthProb/promoter_copy_number`, summed
  n_avg           the normaliser's expectation at tau_inferred, summed
  a               n_actual / n_avg -- the availability term. This is the number
                  the whole question reduces to
  n_ch_real       Cooper-Helmstetter at tau_REAL, as a cross-check that actual
                  copy number tracks the realised doubling time
  factor_mean     (growth_inferred + deg_rate) / n_avg

Reading the result
------------------
  a flat at 1 across the ladder      -> compensation complete, no dosage
                                        feedback, and the reallocation seen
                                        elsewhere is entirely the dilution term
  a falls, and falls MORE for the    -> live copy-number feedback. The spread
  origin-proximal classes               between classes is its size, since only
                                        the differential survives normalisation
  sees_frac near 0                   -> the sensor is blind to the burden, so
                                        the dilution term is also switched off;
                                        expect the stable-RNA signature in
                                        `new_gene_vs_unregulated_genes` to
                                        vanish in the same run

Note the two arms are gated by the same number: `loss` uses the inferred growth
too, so a sensor that misses the slowdown suppresses the dilution arm at the
same time as it releases the dosage arm.

Aggregation: `a` is reported as sum(n_actual) / sum(n_avg) over the class, not
the mean of per-gene ratios, because the share allocation is copy-weighted.

Class overlap warning: the RNAP-subunit and ribosomal-protein classes share
transcription units (E. coli co-transcribes rpoA and rpoBC inside r-protein
operons), so those two rows must not be summed.

Emits a CSV alongside the PDF.
"""

import csv
import os
import pickle

from matplotlib import pyplot as plt
import numpy as np

from models.ecoli.analysis import variantAnalysisPlot
from wholecell.analysis.analysis_tools import (exportFigure,
	read_stacked_columns)
from wholecell.utils import fitting, units

IGNORE_FIRST_N_GENS = 16

# Classes reported, in the order the CSV lists them. rnap_subunits and
# ribosomal_proteins overlap; see the module docstring.
CLASS_ORDER = ['rrna', 'rnap_subunits', 'ribosomal_proteins', 'trna',
	'mrna_not_construct', 'construct', 'genome']


def _window(n_generation):
	if n_generation <= IGNORE_FIRST_N_GENS:
		return None
	return np.arange(IGNORE_FIRST_N_GENS, n_generation)


def _sem(values):
	"""Return the standard error of the mean. NaN if fewer than two values."""
	v = np.asarray(values, dtype=float).ravel()
	v = v[np.isfinite(v)]
	if v.size < 2:
		return float('nan')
	return float(np.std(v, ddof=1) / np.sqrt(v.size))


def _tau_minutes(x):
	"""Reduce one cell's time column to its doubling time in minutes."""
	return np.array([[(x[-1, 0] - x[0, 0]) / 60.0]])


def _time_mean(x):
	"""Per cell: mean over that cell's timesteps, one row per cell."""
	return np.asarray(x).mean(axis=0, keepdims=True)


class Plot(variantAnalysisPlot.VariantAnalysisPlot):
	def do_plot(self, inputDir, plotOutDir, plotOutFileName, simDataFile,
			validationDataFile, metadata):
		variants = sorted(self.ap.get_variants())
		if not variants:
			print('No variants found.')
			return

		generations = _window(self.ap.n_generation)
		if generations is None:
			print('Run has only %d generations; using all of them.'
				% self.ap.n_generation)
			generations = np.arange(self.ap.n_generation)

		with open(simDataFile, 'rb') as handle:
			sim_data = pickle.load(handle)

		masks, frac = self._classes(sim_data)
		transcription = sim_data.process.transcription
		self.wt_coords = np.asarray(
			transcription.rna_data['wt_replication_coordinate'])
		self.deg_rate = transcription.rna_data['deg_rate'].asNumber(
			1 / units.s)
		self.growth_params = transcription._ppgpp_growth_parameters
		self.get_n_avg = sim_data.process.replication.get_average_copy_number
		self.masks, self.frac = masks, frac

		rows = []
		for variant in variants:
			measured = self._measure(variant, generations)
			if measured is None:
				print('No usable cells for variant %d; skipping.' % variant)
				continue
			rows.extend(measured)

		if not rows:
			print('Nothing measurable.')
			return

		self._write_csv(plotOutDir, plotOutFileName, rows)
		self._report(rows)
		self._plot(plotOutDir, plotOutFileName, rows, metadata)

	# ------------------------------------------------------------ classes ----

	def _classes(self, sim_data):
		"""TU-level masks and each unit's replichore fraction."""
		transcription = sim_data.process.transcription
		rna_data = transcription.rna_data
		cistron_data = transcription.cistron_data.struct_array
		ids = list(rna_data['id'])

		coords = np.asarray(rna_data['replication_coordinate'])
		lengths = sim_data.process.replication.replichore_lengths
		frac = np.where(coords > 0, coords / lengths[0], -coords / lengths[1])

		new_tu = np.zeros(len(ids), dtype=bool)
		try:
			new_cistrons = sorted(set(
				cistron_data[cistron_data['is_new_gene']]['id'].tolist()))
			if new_cistrons:
				matrix = transcription.cistron_tu_mapping_matrix.toarray()
				cistron_ids = list(cistron_data['id'])
				sel = [cistron_ids.index(c) for c in new_cistrons
					if c in cistron_ids]
				new_tu = np.asarray(
					matrix[np.array(sel), :].sum(axis=0)).ravel() > 0
		except Exception as exc:  # noqa: BLE001 - construct row is optional
			print('Could not locate the construct (%s); continuing.' % exc)

		is_rnap = np.asarray(rna_data['includes_RNAP'], dtype=bool)
		is_rprot = np.asarray(
			rna_data['includes_ribosomal_protein'], dtype=bool)
		is_mrna = np.asarray(rna_data['is_mRNA'], dtype=bool)
		masks = {
			'rrna': np.asarray(rna_data['is_rRNA'], dtype=bool),
			'rnap_subunits': is_rnap,
			'ribosomal_proteins': is_rprot,
			'trna': np.asarray(rna_data['is_tRNA'], dtype=bool),
			'mrna_not_construct': is_mrna & ~new_tu,
			'construct': new_tu,
			'genome': np.ones(len(ids), dtype=bool),
		}
		overlap = int((is_rnap & is_rprot).sum())
		if overlap:
			print('NOTE: %d transcription units are in both the RNAP and '
				'ribosomal-protein classes; do not sum those rows.' % overlap)
		return masks, frac

	# ------------------------------------------------------------ measure ----

	def _inferred(self, ppgpp):
		"""Growth rate and doubling time the normaliser derives from ppGpp."""
		y = fitting.interpolate_linearized_fit(ppgpp, *self.growth_params)
		growth = max(float(y), 0.0)
		if growth <= 0:
			return 0.0, float('inf')
		return growth, float(np.log(2) / growth / 60)

	def _measure(self, variant, generations):
		cell_paths = self.ap.get_cells(
			variant=[variant], generation=generations)
		if len(cell_paths) == 0:
			return None

		# All three reads must cover the same cells, so restrict to cells that
		# have every listener before reading any of them.
		cell_paths = np.array([p for p in cell_paths
			if all(os.path.exists(os.path.join(p, 'simOut', t,
				'attributes.json'))
				for t in ('Main', 'GrowthLimits', 'RnaSynthProb'))])
		if len(cell_paths) == 0:
			print('  variant %d: no cell has all three listeners.' % variant)
			return None

		taus = read_stacked_columns(cell_paths, 'Main', 'time',
			ignore_exception=True, fun=_tau_minutes).ravel()
		ppgpp = read_stacked_columns(cell_paths, 'GrowthLimits', 'ppgpp_conc',
			ignore_exception=True, fun=_time_mean).ravel()
		copies = read_stacked_columns(cell_paths, 'RnaSynthProb',
			'promoter_copy_number', ignore_exception=True, fun=_time_mean)

		if taus.size == 0 or ppgpp.size == 0 or copies.size == 0:
			return None
		if not (taus.size == ppgpp.size == copies.shape[0]):
			print('  variant %d: read lengths disagree (%d/%d/%d); skipping.'
				% (variant, taus.size, ppgpp.size, copies.shape[0]))
			return None
		if np.allclose(ppgpp, 0):
			print('  variant %d: ppgpp_conc is all zero -- ppGpp regulation is '
				'probably off, and this analysis does not apply.' % variant)
			return None

		# Per cell: the tau the normaliser used, and its expectation for every
		# transcription unit. One call per cell rather than per timestep; the
		# within-cycle variation in ppGpp is second order for this purpose.
		n_avg = np.zeros_like(copies)
		growths = np.zeros(taus.size)
		taus_inf = np.zeros(taus.size)
		for i, p in enumerate(ppgpp):
			growths[i], taus_inf[i] = self._inferred(float(p))
			n_avg[i, :] = np.asarray(
				self.get_n_avg(taus_inf[i], self.wt_coords))
		n_ch_real = np.array([np.asarray(self.get_n_avg(t, self.wt_coords))
			for t in taus])

		rows = []
		for name in CLASS_ORDER:
			mask = self.masks.get(name)
			if mask is None or not mask.any():
				continue
			# Copy-weighted, because the share allocation is copy-weighted.
			act = copies[:, mask].sum(axis=1)
			exp = n_avg[:, mask].sum(axis=1)
			ref = n_ch_real[:, mask].sum(axis=1)
			with np.errstate(divide='ignore', invalid='ignore'):
				a = np.where(exp > 0, act / exp, np.nan)
				a_ref = np.where(ref > 0, act / ref, np.nan)
			loss = growths[:, None] + self.deg_rate[None, mask]
			with np.errstate(divide='ignore', invalid='ignore'):
				factor = np.where(n_avg[:, mask] > 0,
					loss / n_avg[:, mask], np.nan)
			rows.append(dict(
				variant=variant, gene_class=name,
				n_cells=int(taus.size), n_tus=int(mask.sum()),
				mean_replichore_fraction=float(self.frac[mask].mean()),
				tau_real=float(np.mean(taus)), tau_real_sem=_sem(taus),
				ppgpp_conc=float(np.mean(ppgpp)),
				tau_inferred=float(np.mean(taus_inf)),
				tau_inferred_sem=_sem(taus_inf),
				n_actual=float(np.mean(act)),
				n_avg=float(np.mean(exp)),
				n_ch_real=float(np.mean(ref)),
				a=float(np.nanmean(a)), a_sem=_sem(a),
				a_vs_real_tau=float(np.nanmean(a_ref)),
				factor_mean=float(np.nanmean(factor)),
				))
		return rows

	# ------------------------------------------------------------- output ----

	FIELDS = ['variant', 'gene_class', 'n_cells', 'n_tus',
		'mean_replichore_fraction', 'tau_real', 'tau_real_sem', 'ppgpp_conc',
		'tau_inferred', 'tau_inferred_sem', 'n_actual', 'n_avg', 'n_ch_real',
		'a', 'a_sem', 'a_vs_real_tau', 'factor_mean']

	def _write_csv(self, plot_out_dir, plot_out_file_name, rows):
		path = os.path.join(plot_out_dir, plot_out_file_name + '.csv')
		with open(path, 'w') as handle:
			w = csv.DictWriter(handle, fieldnames=self.FIELDS)
			w.writeheader()
			for r in rows:
				w.writerow({k: r.get(k, '') for k in self.FIELDS})
		print('\nWrote %s' % path)

	def _report(self, rows):
		variants = sorted({r['variant'] for r in rows})
		gen = [r for r in rows if r['gene_class'] == 'genome']
		print('\nThe doubling time the normaliser used, against the real one')
		print('  %-8s %-10s %-12s %-12s %s'
			% ('variant', 'ppGpp uM', 'tau_real', 'tau_inferred', 'ratio'))
		for r in gen:
			print('  %-8d %-10.2f %-12.2f %-12.2f %.3f'
				% (r['variant'], r['ppgpp_conc'], r['tau_real'],
					r['tau_inferred'], r['tau_inferred'] / r['tau_real']))
		if len(gen) >= 2:
			lo, hi = gen[1] if len(gen) > 2 else gen[0], gen[-1]
			d_real = hi['tau_real'] - lo['tau_real']
			d_inf = hi['tau_inferred'] - lo['tau_inferred']
			sees = d_inf / d_real if d_real else float('nan')
			print('\n  The sensor sees %.0f%% of the realised slowdown '
				'(%.1f of %.1f min).' % (100 * sees, d_inf, d_real))
			if sees > 0.9:
				print('  VERDICT: compensation is essentially complete. Gene '
					'dosage loss cancels, so it cannot drive a feedback loop, '
					'and any reallocation must come from the dilution term.')
			elif sees > 0.4:
				print('  VERDICT: compensation is PARTIAL. A real copy-number '
					'feedback survives, sized by the spread in `a` below.')
			else:
				print('  VERDICT: the sensor is effectively BLIND to this '
					'burden. Gene dosage loss is barely compensated and the '
					'dosage arm is live -- but the dilution term, which uses '
					'the same inferred growth, is suppressed in the same '
					'measure. Expect stable RNA to show no extra loss in '
					'new_gene_vs_unregulated_genes on this run.')

		print('\nAvailability term a = n_actual / n_avg, by class')
		hdr = '  %-22s %-6s' % ('class', 'f')
		for v in variants:
			hdr += ' %8s' % ('v%d' % v)
		print(hdr + '   change')
		for name in CLASS_ORDER:
			sub = [r for r in rows if r['gene_class'] == name]
			if not sub:
				continue
			line = '  %-22s %-6.3f' % (name, sub[0]['mean_replichore_fraction'])
			for r in sub:
				line += ' %8.4f' % r['a']
			first = sub[1]['a'] if len(sub) > 2 else sub[0]['a']
			line += '   %+.1f%%' % (100 * (sub[-1]['a'] / first - 1))
			print(line)
		print('\n  Only the SPREAD between classes survives the L1 '
			'normalisation; a common shift cancels.')

		# The control was previously described but never printed, which made it
		# look like a missing block. Print the numbers.
		print('\nControl: a_vs_real_tau = n_actual / Cooper-Helmstetter at the '
			'REALISED tau')
		hdr = '  %-22s %-6s' % ('class', 'f')
		for v in variants:
			hdr += ' %8s' % ('v%d' % v)
		print(hdr)
		for name in CLASS_ORDER:
			sub = [r for r in rows if r['gene_class'] == name]
			if not sub:
				continue
			line = '  %-22s %-6.3f' % (name, sub[0]['mean_replichore_fraction'])
			for r in sub:
				line += ' %8.4f' % r['a_vs_real_tau']
			print(line)
		ctrl = [r['a_vs_real_tau'] for r in rows if r['variant'] >= 1
			and np.isfinite(r['a_vs_real_tau'])]
		if ctrl:
			print('  range %.4f - %.4f over every class and variant'
				% (min(ctrl), max(ctrl)))
			print('  This should stay near 1 and be FLAT across classes. Flat '
				'means actual copy\n  number tracks the realised doubling '
				'time uniformly, so any f-dependence in\n  `a` above is the '
				'ppGpp sensor rather than the cell. If this column develops '
				'its\n  own ordering in f, that interpretation fails.')

	def _plot(self, plot_out_dir, plot_out_file_name, rows, metadata):
		gen = [r for r in rows if r['gene_class'] == 'genome']
		fig, axes = plt.subplots(1, 3, figsize=(13, 3.6))

		v = [r['variant'] for r in gen]
		axes[0].plot(v, [r['tau_real'] for r in gen], 'o-', color='#c2410c',
			label='realised')
		axes[0].plot(v, [r['tau_inferred'] for r in gen], 's--',
			color='#4a3aa7', label='inferred from ppGpp')
		axes[0].set_ylabel('doubling time (min)')
		axes[0].set_title('What the normaliser thinks')

		for name in ('rrna', 'rnap_subunits', 'mrna_not_construct', 'genome'):
			sub = [r for r in rows if r['gene_class'] == name]
			if not sub:
				continue
			axes[1].plot([r['variant'] for r in sub], [r['a'] for r in sub],
				'o-', label='%s (f=%.2f)'
				% (name, sub[0]['mean_replichore_fraction']))
		axes[1].axhline(1.0, color='#767469', lw=0.8, ls=':')
		axes[1].set_ylabel('a = n_actual / n_avg')
		axes[1].set_title('Uncompensated dosage')

		for name in ('rrna', 'mrna_not_construct', 'genome'):
			sub = [r for r in rows if r['gene_class'] == name]
			if not sub:
				continue
			axes[2].plot([r['variant'] for r in sub],
				[r['a_vs_real_tau'] for r in sub], 'o-', label=name)
		axes[2].axhline(1.0, color='#767469', lw=0.8, ls=':')
		axes[2].set_ylabel('n_actual / C-H at realised tau')
		axes[2].set_title('Control: copies track real tau')

		for ax in axes:
			ax.set_xlabel('variant')
			ax.legend(fontsize=7, frameon=False)
		plt.tight_layout()
		exportFigure(plt, plot_out_dir, plot_out_file_name, metadata)
		plt.close('all')


if __name__ == '__main__':
	Plot().cli()
