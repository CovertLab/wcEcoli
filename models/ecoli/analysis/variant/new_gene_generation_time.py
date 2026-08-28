"""
Doubling time, generation by generation, per lineage.

The simplest possible test of whether the copy-number feedback does anything,
and deliberately the least processed analysis in this study.

Two paths run from GFP induction to a slower cell. GFP translation competes for
ribosomes immediately, within a timestep. Gene dosage loss cannot act until
growth has already slowed, an initiation has been delayed, and the copies have
gone -- at least one full cycle, realistically one to two. So:

	one drop after induction, then recovery  ->  only the direct path
	a drop, then a SECOND drop a generation
	or two later                             ->  the dosage path carries flux

The frozen batch is the contrast. Its dosage channel is roughly twice as strong,
so a second drop should be deeper there. Matched arms, matched medium.

Why this rather than something richer
-------------------------------------
Growth rate is set by ribosome content, so ribosome content and doubling time
are close to the same variable, and mass accumulation is not an independent
check on either. Earlier designs in this study measured per-event windows,
concentrations and derivatives, and each one added a confound: cell-cycle phase
coherence, division halving, window-versus-cycle geometry, control-trace
overcorrection. Division times have none of those. They need no alignment,
because induction happens at a fixed generation index; no smoothing, because
there is one value per cell; and no normalisation, because minutes are minutes.

Ribosome content per generation is emitted as a secondary column, since the
growth law need not hold during a transient and a divergence between the two
would itself be informative. The doubling time is the headline.

Output
------
	<name>.csv              one row per (variant, seed, generation)
	<name>_by_gen.csv       per (variant, generation): mean, SEM, n lineages
	<name>.pdf / .png       doubling time against generation, one line per rung
"""

import csv
import os
import pickle

from matplotlib import pyplot as plt
import numpy as np

from models.ecoli.analysis import variantAnalysisPlot
from wholecell.analysis.analysis_tools import (exportFigure,
	first_cell_with_table, read_stacked_bulk_molecules)
from wholecell.io.tablereader import TableReader

# Generation at which new gene expression is switched on. Must match
# NEW_GENE_INDUCTION_GEN in models/ecoli/sim/variants/new_gene_internal_shift.py.
INDUCTION_GEN = 8

# Timesteps averaged at each end of a cycle for the ribosome-content columns. A
# single timestep at the division boundary can be read mid-update.
EDGE_STEPS = 5


def _sem(values):
	v = np.asarray(values, dtype=float)
	v = v[np.isfinite(v)]
	if v.size < 2:
		return float('nan')
	return float(np.std(v, ddof=1) / np.sqrt(v.size))


class Plot(variantAnalysisPlot.VariantAnalysisPlot):
	def do_plot(self, inputDir, plotOutDir, plotOutFileName, simDataFile,
			validationDataFile, metadata):
		with open(simDataFile, 'rb') as handle:
			sim_data = pickle.load(handle)

		variants = sorted(self.ap.get_variants())
		if not variants:
			print('No variants found.')
			return

		rows = []
		for variant in variants:
			for seed in sorted(self.ap.get_seeds(variant=variant)):
				rows.extend(self._one_lineage(variant, seed, sim_data))
		if not rows:
			print('No complete generations found.')
			return

		by_gen = self._aggregate(rows)
		self._write(plotOutDir, plotOutFileName, rows, by_gen)
		self._print_summary(by_gen)
		self._plot(plotOutDir, plotOutFileName, by_gen)

	# ---- extraction ------------------------------------------------------

	def _one_lineage(self, variant, seed, sim_data):
		"""
		One row per COMPLETED generation.

		A generation counts only if the next one exists, which is how we know
		the cell divided rather than being truncated by the end of the run or by
		stalling. That exclusion is survivorship-biased -- the cells that stall
		are the burdened ones -- so the count of contributing lineages is
		reported per generation and must be read alongside the means.
		"""
		present = {}
		for gen in range(self.ap.n_generation):
			paths = self.ap.get_cells(variant=[variant], seed=[seed],
				generation=[gen])
			if len(paths):
				present[gen] = paths[0]

		out = []
		for gen, path in sorted(present.items()):
			if gen + 1 not in present:
				continue          # did not divide, or end of the run
			try:
				t = TableReader(os.path.join(path, 'simOut',
					'Main')).readColumn('time')
			except Exception as exc:  # noqa: BLE001 - skip the generation
				print('Variant %d seed %d gen %d: no time column (%s).'
					% (variant, seed, gen, exc))
				continue
			if t.size < 2 * EDGE_STEPS:
				continue
			row = dict(variant=variant, seed=seed, generation=gen,
				tau_min=float((t[-1] - t[0]) / 60.0),
				n_timesteps=int(t.size),
				post_induction=int(gen >= INDUCTION_GEN))
			row.update(self._content(path, sim_data))
			out.append(row)
		return out

	@staticmethod
	def _content(path, sim_data):
		"""
		Ribosome content at the start and end of the cycle.

		Secondary to the doubling time. Reported because the growth law -- that
		growth rate tracks ribosome content -- holds at balanced growth and need
		not hold during a transient, so a divergence between the two is
		informative on its own.
		"""
		blank = dict(ribosome_birth=float('nan'), ribosome_div=float('nan'),
			mass_birth=float('nan'), mass_div=float('nan'),
			content_ratio=float('nan'))
		try:
			umc = TableReader(os.path.join(path, 'simOut',
				'UniqueMoleculeCounts'))
			ids = umc.readAttribute('uniqueMoleculeIds')
			counts = umc.readColumn('uniqueMoleculeCounts')
			mass = TableReader(os.path.join(path, 'simOut',
				'Mass')).readColumn('cellMass')
			if 'active_ribosome' not in ids:
				return blank
			active = counts[:, ids.index('active_ribosome')]
			(s30, s50) = read_stacked_bulk_molecules([path],
				([sim_data.molecule_ids.s30_full_complex],
					[sim_data.molecule_ids.s50_full_complex]),
				ignore_exception=True)
			ribo = active + np.minimum(np.asarray(s30).squeeze(),
				np.asarray(s50).squeeze())
		except Exception as exc:  # noqa: BLE001 - the headline does not need it
			print('  ribosome content unavailable (%s)' % exc)
			return blank

		e = EDGE_STEPS
		rb, rd = float(np.mean(ribo[:e])), float(np.mean(ribo[-e:]))
		mb, md = float(np.mean(mass[:e])), float(np.mean(mass[-e:]))
		ratio = ((rd / rb) / (md / mb)
			if rb > 0 and mb > 0 and md > 0 else float('nan'))
		return dict(ribosome_birth=rb, ribosome_div=rd, mass_birth=mb,
			mass_div=md, content_ratio=ratio)

	# ---- aggregation and output ------------------------------------------

	@staticmethod
	def _aggregate(rows):
		out = []
		keys = sorted({(r['variant'], r['generation']) for r in rows})
		for variant, gen in keys:
			sub = [r for r in rows
				if r['variant'] == variant and r['generation'] == gen]
			taus = [r['tau_min'] for r in sub]
			cont = [r['content_ratio'] for r in sub
				if np.isfinite(r['content_ratio'])]
			out.append(dict(variant=variant, generation=gen,
				n_lineages=len(sub),
				tau_mean=float(np.mean(taus)), tau_sem=_sem(taus),
				tau_median=float(np.median(taus)),
				content_ratio_mean=float(np.mean(cont)) if cont
					else float('nan'),
				content_ratio_sem=_sem(cont)))
		return out

	@staticmethod
	def _write(plot_out_dir, name, rows, by_gen):
		p1 = os.path.join(plot_out_dir, name + '.csv')
		fields = ['variant', 'seed', 'generation', 'post_induction', 'tau_min',
			'n_timesteps', 'ribosome_birth', 'ribosome_div', 'mass_birth',
			'mass_div', 'content_ratio']
		with open(p1, 'w', newline='') as h:
			w = csv.DictWriter(h, fieldnames=fields)
			w.writeheader()
			for r in rows:
				w.writerow(r)
		p2 = os.path.join(plot_out_dir, name + '_by_gen.csv')
		with open(p2, 'w', newline='') as h:
			w = csv.DictWriter(h, fieldnames=list(by_gen[0].keys()))
			w.writeheader()
			for r in by_gen:
				w.writerow(r)
		print('Wrote %s\nWrote %s' % (p1, p2))

	@staticmethod
	def _print_summary(by_gen):
		variants = sorted({r['variant'] for r in by_gen})
		gens = sorted({r['generation'] for r in by_gen})
		look = {(r['variant'], r['generation']): r for r in by_gen}

		print('\nDoubling time by generation, minutes. Induction at generation '
			'%d.' % INDUCTION_GEN)
		print('A SECOND fall one to two generations after induction is the '
			'dosage path.')
		print('An immediate fall and then recovery is the direct path alone.')
		print()
		header = '%-4s' % 'gen'
		for v in variants:
			header += ' %13s' % ('v%d' % v)
		print(header)
		for gen in gens:
			line = '%-4d' % gen
			for v in variants:
				r = look.get((v, gen))
				line += ' %13s' % ('%.1f+/-%.1f' % (r['tau_mean'],
					r['tau_sem']) if r else '-')
			print(line + ('   <- induction' if gen == INDUCTION_GEN else ''))

		print('\nContributing lineages per generation. Read this WITH the '
			'means above:')
		print('cells that stall never divide, and the ones that stall are the '
			'burdened ones,')
		print('so a falling count means the surviving means are optimistic.')
		print()
		print(header)
		for gen in gens:
			line = '%-4d' % gen
			for v in variants:
				r = look.get((v, gen))
				line += ' %13s' % (r['n_lineages'] if r else '-')
			print(line)

	@staticmethod
	def _plot(plot_out_dir, name, by_gen):
		variants = sorted({r['variant'] for r in by_gen})
		look = {(r['variant'], r['generation']): r for r in by_gen}
		gens = sorted({r['generation'] for r in by_gen})
		fig, axes = plt.subplots(1, 2, figsize=(11, 4))
		cmap = plt.get_cmap('viridis')
		shade = {v: cmap(0.1 + 0.75 * i / max(1, len(variants) - 1))
			for i, v in enumerate(variants)}

		for key, ax, ylabel in (('tau_mean', axes[0], 'doubling time (min)'),
				('content_ratio_mean', axes[1],
					'ribosome content ratio\n(1 = kept pace with mass)')):
			for v in variants:
				xs = [g for g in gens if (v, g) in look]
				ys = [look[(v, g)][key] for g in xs]
				es = [look[(v, g)][key.replace('_mean', '_sem')] for g in xs]
				ax.errorbar(xs, ys, yerr=es, marker='o', ms=3.5, lw=1.3,
					color=shade[v], label='v%d' % v, capsize=2)
			ax.axvline(INDUCTION_GEN - 0.5, color='C3', lw=1.4)
			ax.set_xlabel('generation')
			ax.set_ylabel(ylabel)
		axes[1].axhline(1.0, color='0.5', lw=0.8, ls='--')
		axes[0].legend(fontsize='xx-small', ncol=3, frameon=False)
		axes[0].set_title('One fall or two?', fontsize='small')
		axes[1].set_title('Secondary: did ribosomes keep pace with mass?',
			fontsize='small')
		plt.tight_layout()
		exportFigure(plt, plot_out_dir, name, metadata=None)
		plt.close('all')


if __name__ == '__main__':
	Plot().cli()
