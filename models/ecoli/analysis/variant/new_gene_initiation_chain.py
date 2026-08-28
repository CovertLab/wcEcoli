"""
The transcription initiation chain, per generation, per seed, per gene class.

Every earlier analysis in this study averages over generations 16-23, which
hides the one thing that would settle whether gene dosage loss matters: GFP's
direct proteome cost lands immediately, whereas the copy-number route cannot act
until growth has slowed, a replication round has been lost, and the copies have
gone. Those two effects are separated in TIME, not in magnitude.

This traces the whole chain generation by generation:

	expression x loss / EXPECTED copies  -> p           (per transcription unit)
	ACTUAL copies x p, renormalised      -> share       (per promoter)
	RNAP budget x share                  -> initiations (per timestep)
	                                     -> pools, growth

Most of that is logged rather than reconstructed. Only five quantities are not:
f_ppgpp, the blended expression vector, inferred growth and tau, the loss term,
and the expected copy number. All five are validated in one shot -- see below.

Do these numbers match what the process actually saw?
-----------------------------------------------------
Not automatically, and the script tests rather than assumes.

  ppGpp. GrowthLimits/ppgpp_conc is written by the elongation model; transcript
  initiation computes its own from the same molecule count. Two computations, one
  timestep. Units are identical (both umol/L) but the values need not be.

  Promoter copy number. The process reads the promoter container twice, and the
  read at evolveState is what builds target_rna_synth_prob. The listener reads it
  again in its own update(), AFTER every process has run, so replication or
  division in between makes the logged value differ from the one used.

Both are covered by the identity checks. Identity 1 compares our reconstruction
against basal_prob_ppgpp_synth_prob, which is the process's own output, so
passing it means we recovered the process's inputs -- not merely that our
arithmetic is self-consistent. Identity 4 does the same for copy number, and if
it fails, `copies_implied` (= target / p, proportional to the copy-number vector
the process actually multiplied by) is emitted as the authoritative fallback.

Verified before this script was written: the reconstruction reproduces
synth_prob_from_ppgpp to a maximum relative error of 5e-16, with and without the
frozen-expectation override. The gate was then tested against deliberately
broken inputs -- a frozen batch reconstructed as if unfrozen, which is the
silent failure that defeated the same fix in new_gene_dosage_compensation, trips
identity 1 at 5.2e-2, seven orders above tolerance.

One limitation to know. Identity 4 is blind to a UNIFORM error in the copy-number
column, because a scaling common to every gene cancels in the promoter-level
normalisation -- halving every copy number leaves it at 0.0. It catches per-gene
errors: halving only the rRNA operons trips it at 4.7e-1. So a globally stale
copy-number column would pass identity 4 while making every `a` wrong by a
constant factor. The spread between classes, which is what the argument turns
on, would still be right. `copies_implied` is emitted so that constant can be
recovered if it matters.

Output
------
	<name>_cell.csv      one row per (variant, seed, generation, aggregation)
	<name>_class.csv     ... x gene_class
	<name>_by_gen.csv    aggregated across seeds, with n_lineages
	<name>_identity.csv  the six checks
"""

import csv
import os
import pickle

from matplotlib import pyplot as plt
import numpy as np

from models.ecoli.analysis import variantAnalysisPlot
from wholecell.analysis.analysis_tools import (exportFigure,
	first_cell_with_table, read_bulk_molecule_counts)
from wholecell.io.tablereader import TableReader
from wholecell.utils import fitting, units
from wholecell.utils.fitting import normalize

# Generation at which new gene expression is switched on. Must match
# NEW_GENE_INDUCTION_GEN in models/ecoli/sim/variants/new_gene_internal_shift.py.
INDUCTION_GEN = 8

# Timesteps averaged at each end for the birth and mid-cycle aggregations.
EDGE_STEPS = 5

# Replichore fraction bounds for the two positional reference classes.
#
# terminus_ref follows the existing mask idiom in new_gene_machinery_allocation
# rather than the single-gene idiom elsewhere, because one gene is noisy.
#
# origin_ref does not exist anywhere in the repo and has to exclude rRNA, RNAP
# subunits AND ribosomal proteins: all three cluster near the origin, so an
# unfiltered origin reference would just be the machinery again, and the whole
# point of the pair is to be a positional control that carries no machinery.
TERMINUS_FRACTION = 0.90
ORIGIN_FRACTION = 0.10

# Timesteps sampled per cell for the per-transcription-unit identity checks.
# The full arrays are (n_timesteps, 3266); reducing to class sums first keeps
# memory bounded, so the identities are checked on a subsample instead.
IDENTITY_SAMPLES = 5

AGGREGATIONS = ('time_mean', 'birth', 'mid_cycle')

CLASS_ORDER = ('rnap_subunits', 'rrna', 'ribosomal_proteins', 'construct',
	'mrna_not_construct', 'origin_ref', 'terminus_ref', 'genome')


def _sem(values):
	v = np.asarray(values, dtype=float)
	v = v[np.isfinite(v)]
	if v.size < 2:
		return float('nan')
	return float(np.std(v, ddof=1) / np.sqrt(v.size))


def _slices(n):
	"""Timestep slices for the three aggregations."""
	e = min(EDGE_STEPS, max(1, n // 4))
	mid = max(0, n // 2 - e // 2)
	return dict(time_mean=slice(None), birth=slice(0, e),
		mid_cycle=slice(mid, mid + e))


class Plot(variantAnalysisPlot.VariantAnalysisPlot):
	def do_plot(self, inputDir, plotOutDir, plotOutFileName, simDataFile,
			validationDataFile, metadata):
		with open(simDataFile, 'rb') as handle:
			sim_data = pickle.load(handle)

		variants = sorted(self.ap.get_variants())
		if not variants:
			print('No variants found.')
			return

		ctx = self._context(sim_data)
		self._report_classes(ctx)

		cell_rows, class_rows, ident_rows = [], [], []
		self._dropped = dict(stale=0, unreadable=0, incomplete=0)
		for variant in variants:
			ctx['frozen_tau'] = self._frozen_tau(variant)
			for path, seed, gen in self._cells(variant):
				out = self._one_cell(path, variant, seed, gen, ctx)
				if out is None:
					continue
				cell_rows.extend(out[0])
				class_rows.extend(out[1])
				ident_rows.extend(out[2])

		if not class_rows:
			print('No readable generations.')
			return

		print('\nTimesteps dropped: %d stale (budget zero or no chromosome).'
			% self._dropped['stale'])
		print('Cells dropped: %d unreadable, %d incomplete (did not divide).'
			% (self._dropped['unreadable'], self._dropped['incomplete']))

		gate = self._print_identities(ident_rows)
		by_gen = self._aggregate(class_rows, cell_rows)
		self._write(plotOutDir, plotOutFileName, cell_rows, class_rows,
			ident_rows, by_gen)
		self._print_summary(by_gen)
		if not gate:
			print('\n*** IDENTITY 1 FAILED. The reconstruction does not match '
				'what the process used.')
			print('*** Chain numbers above are NOT trustworthy. Fix before '
				'interpreting.')
		self._plot(plotOutDir, plotOutFileName, by_gen)

	# ---- context from sim_data -------------------------------------------

	@staticmethod
	def _context(sim_data):
		"""Masks, rate constants and fit parameters, all TU-indexed."""
		tr = sim_data.process.transcription
		rna = tr.rna_data.struct_array
		cistron = tr.cistron_data.struct_array
		ids = [str(i) for i in rna['id']]

		coords = np.asarray(rna['replication_coordinate'], dtype=float)
		lengths = sim_data.process.replication.replichore_lengths
		frac = np.where(coords > 0, coords / lengths[0], -coords / lengths[1])

		is_rnap = np.asarray(rna['includes_RNAP'], dtype=bool)
		is_rprot = np.asarray(rna['includes_ribosomal_protein'], dtype=bool)
		is_rrna = np.asarray(rna['is_rRNA'], dtype=bool)
		is_mrna = np.asarray(rna['is_mRNA'], dtype=bool)

		# The construct has no TU-level flag; reach it through the cistron map.
		new_tu = np.zeros(len(ids), dtype=bool)
		try:
			new_cistrons = sorted(set(
				cistron[cistron['is_new_gene']]['id'].tolist()))
			if new_cistrons:
				matrix = tr.cistron_tu_mapping_matrix.toarray()
				cids = list(cistron['id'])
				sel = [cids.index(c) for c in new_cistrons if c in cids]
				new_tu = np.asarray(
					matrix[np.array(sel), :].sum(axis=0)).ravel() > 0
		except Exception as exc:  # noqa: BLE001 - construct is optional
			print('Could not locate the construct (%s).' % exc)

		machinery = is_rrna | is_rnap | is_rprot
		masks = {
			'rnap_subunits': is_rnap,
			'rrna': is_rrna,
			'ribosomal_proteins': is_rprot,
			'construct': new_tu,
			'mrna_not_construct': is_mrna & ~new_tu,
			'origin_ref': (frac < ORIGIN_FRACTION) & ~machinery & ~new_tu,
			'terminus_ref': (frac > TERMINUS_FRACTION) & ~machinery & ~new_tu,
			'genome': np.ones(len(ids), dtype=bool),
			}

		# New gene monomer, for the GFP protein trace.
		monomer = sim_data.process.translation.monomer_data.struct_array
		c2m = dict(zip(monomer['cistron_id'], monomer['id']))
		new_monomers = [c2m[c] for c in
			cistron[cistron['is_new_gene']]['id'].tolist() if c in c2m]

		return dict(ids=ids, masks=masks, frac=frac,
			deg=tr.rna_data['deg_rate'].asNumber(1 / units.s),
			wt_coords=np.asarray(tr.rna_data['wt_replication_coordinate']),
			exp_free=tr.exp_free, exp_ppgpp=tr.exp_ppgpp,
			km_sq=tr._ppgpp_km_squared,
			growth_params=tr._ppgpp_growth_parameters,
			is_rrna=is_rrna,
			get_n_avg=sim_data.process.replication.get_average_copy_number,
			new_monomers=new_monomers,
			mol=sim_data.molecule_ids, frozen_tau=None)

	@staticmethod
	def _report_classes(ctx):
		print('Gene classes, TU-indexed:')
		print('%-22s %7s %9s' % ('class', 'n_tus', 'mean f'))
		for name in CLASS_ORDER:
			m = ctx['masks'][name]
			print('%-22s %7d %9.3f'
				% (name, int(m.sum()),
					float(ctx['frac'][m].mean()) if m.any() else float('nan')))
		if int(ctx['masks']['origin_ref'].sum()) < 20:
			print('WARNING: origin_ref has few members. It excludes rRNA, RNAP '
				'subunits and ribosomal proteins,')
			print('which is deliberate -- but a thin reference should be '
				'reported as weak rather than quoted.')
		overlap = int((ctx['masks']['rnap_subunits']
			& ctx['masks']['ribosomal_proteins']).sum())
		if overlap:
			print('NOTE: %d TUs are in BOTH the RNAP and ribosomal-protein '
				'classes (co-transcribed operons).' % overlap)
			print('Do not sum those two rows.')

	def _frozen_tau(self, variant):
		"""
		Frozen expectation tau, from the PER-VARIANT pickle.

		simDataFile is the base pickle and variant functions mutate a copy, so
		reading the attribute from the base one silently returns None and takes
		the unfrozen path -- the exact failure that defeated the same fix in
		new_gene_dosage_compensation.
		"""
		try:
			with open(self.ap.get_variant_kb(variant), 'rb') as handle:
				sd = pickle.load(handle)
		except Exception as exc:  # noqa: BLE001 - assume unfrozen
			print('Variant %d: could not read per-variant sim_data (%s); '
				'assuming the expectation is not frozen.' % (variant, exc))
			return None
		return getattr(sd.process.transcription, 'frozen_expectation_tau', None)

	def _cells(self, variant):
		"""(path, seed, generation) for every COMPLETED generation."""
		out = []
		by_seed = {}
		for path in self.ap.get_cells(variant=[variant]):
			seed = int(self.ap.get_cell_seed(path))
			gen = int(self.ap.get_cell_generation(path))
			by_seed.setdefault(seed, {})[gen] = path
		for seed, gens in sorted(by_seed.items()):
			for gen, path in sorted(gens.items()):
				# A generation counts only if the next exists, i.e. the cell
				# divided. Survivorship-biased -- the cells that stall are the
				# burdened ones -- so n_lineages is reported per generation.
				if gen + 1 not in gens:
					self._dropped['incomplete'] += 1
					continue
				out.append((path, seed, gen))
		return out

	# ---- per cell --------------------------------------------------------

	def _one_cell(self, path, variant, seed, gen, ctx):
		sim_out = os.path.join(path, 'simOut')
		try:
			t = TableReader(os.path.join(sim_out, 'Main')).readColumn('time')
			rsp = TableReader(os.path.join(sim_out, 'RnaSynthProb'))
			rnap = TableReader(os.path.join(sim_out, 'RnapData'))
			gl = TableReader(os.path.join(sim_out, 'GrowthLimits'))
			ppgpp = np.asarray(gl.readColumn('ppgpp_conc'), dtype=float)
			budget = np.asarray(rnap.readColumn('didInitialize'), dtype=float)
			max_p = np.asarray(rsp.readColumn('max_p'), dtype=float)
		except Exception as exc:  # noqa: BLE001 - a bad cell is skippable
			print('v%d s%d g%d: unreadable (%s).' % (variant, seed, gen, exc))
			self._dropped['unreadable'] += 1
			return None
		n = t.size
		if n < 2 * EDGE_STEPS or ppgpp.size != n:
			self._dropped['unreadable'] += 1
			return None

		# Stale timesteps. When the budget is zero transcript_initiation returns
		# early, after writing target_rna_synth_prob but before max_p,
		# actual_rna_synth_prob, tu_is_overcrowded, rnaInitEvent and
		# didInitialize -- those five carry the previous timestep's values.
		good = (budget > 0) & np.isfinite(ppgpp) & (ppgpp > 0)
		self._dropped['stale'] += int((~good).sum())
		if good.sum() < 2 * EDGE_STEPS:
			return None

		# Reduce over transcription units FIRST, then over time: the raw
		# columns are (n_timesteps, 3266) and reducing the other way round
		# holds several hundred MB per cell.
		series = {}
		for key, table, col in (
				('copies', rsp, 'promoter_copy_number'),
				('p_logged', rsp, 'basal_prob_ppgpp_synth_prob'),
				('p_attn', rsp, 'basal_prob_updated'),
				('share_target', rsp, 'target_rna_synth_prob'),
				('share_actual', rsp, 'actual_rna_synth_prob'),
				('overcrowded', rsp, 'tu_is_overcrowded'),
				('init', rnap, 'rnaInitEvent')):
			try:
				arr = np.asarray(table.readColumn(col), dtype=float)
			except Exception as exc:  # noqa: BLE001
				print('v%d s%d g%d: %s missing (%s).'
					% (variant, seed, gen, col, exc))
				self._dropped['unreadable'] += 1
				return None
			series[key] = {c: arr[:, ctx['masks'][c]].sum(axis=1)
				for c in CLASS_ORDER}
			if key == 'p_logged':
				ident = self._identities(arr, rsp, rnap, ppgpp, budget, good,
					ctx, variant, seed, gen)
			del arr

		cell = self._cell_scalars(sim_out, t, ppgpp, budget, max_p, good, ctx)
		rows_cell, rows_class = [], []
		sl = _slices(int(good.sum()))
		for agg, s in sl.items():
			idx = np.where(good)[0][s]
			base = dict(variant=variant, seed=seed, generation=gen,
				aggregation=agg, post_induction=int(gen >= INDUCTION_GEN))
			rows_cell.append(dict(base, n_timesteps=len(idx),
				**{k: float(np.mean(v[idx])) for k, v in cell.items()}))
			rows_class.extend(
				self._class_rows(base, idx, series, ppgpp, ctx))
		return rows_cell, rows_class, ident

	def _cell_scalars(self, sim_out, t, ppgpp, budget, max_p, good, ctx):
		"""Per-timestep cell-level quantities, all same length as t."""
		f_ppgpp = ppgpp ** 2 / (ctx['km_sq'] + ppgpp ** 2)
		growth = np.array([max(float(fitting.interpolate_linearized_fit(
			p, *ctx['growth_params'])), 0.0) for p in ppgpp])
		with np.errstate(divide='ignore', invalid='ignore'):
			tau_inf = np.where(growth > 0, np.log(2) / growth / 60, np.nan)
		tau_used = (np.full_like(tau_inf, ctx['frozen_tau'])
			if ctx['frozen_tau'] is not None else tau_inf)
		out = dict(ppgpp_conc=ppgpp, f_ppgpp=f_ppgpp, growth_inferred=growth,
			tau_inferred=tau_inf, tau_used=tau_used, budget=budget,
			max_p=max_p)
		out['tau_realised'] = np.full_like(ppgpp, (t[-1] - t[0]) / 60.0)
		out['frozen_tau'] = np.full_like(
			ppgpp, ctx['frozen_tau'] if ctx['frozen_tau'] is not None
			else np.nan)

		try:
			umc = TableReader(os.path.join(sim_out, 'UniqueMoleculeCounts'))
			uids = umc.readAttribute('uniqueMoleculeIds')
			counts = umc.readColumn('uniqueMoleculeCounts')
			out['active_rnap'] = counts[:, uids.index('active_RNAP')].astype(float)
			out['active_ribosome'] = counts[
				:, uids.index('active_ribosome')].astype(float)
		except Exception:  # noqa: BLE001 - pools are informative only
			out['active_rnap'] = np.full_like(ppgpp, np.nan)
			out['active_ribosome'] = np.full_like(ppgpp, np.nan)
		try:
			mol = ctx['mol']
			names = [mol.full_RNAP, mol.s30_full_complex, mol.s50_full_complex]
			free = read_bulk_molecule_counts(sim_out, (names,))[0]
			out['inactive_rnap'] = np.asarray(free[:, 0], dtype=float)
			out['inactive_ribosome'] = np.minimum(
				np.asarray(free[:, 1], dtype=float),
				np.asarray(free[:, 2], dtype=float))
		except Exception:  # noqa: BLE001
			out['inactive_rnap'] = np.full_like(ppgpp, np.nan)
			out['inactive_ribosome'] = np.full_like(ppgpp, np.nan)
		out['total_rnap'] = out['active_rnap'] + out['inactive_rnap']
		out['total_ribosome'] = out['active_ribosome'] + out['inactive_ribosome']
		try:
			gfp = read_bulk_molecule_counts(sim_out, (ctx['new_monomers'],))[0]
			out['gfp_protein'] = np.asarray(gfp, dtype=float).sum(axis=1)
		except Exception:  # noqa: BLE001
			out['gfp_protein'] = np.full_like(ppgpp, np.nan)
		try:
			mass = TableReader(os.path.join(sim_out, 'Mass'))
			out['cell_mass'] = np.asarray(mass.readColumn('cellMass'), float)
			out['dry_mass'] = np.asarray(mass.readColumn('dryMass'), float)
		except Exception:  # noqa: BLE001
			out['cell_mass'] = np.full_like(ppgpp, np.nan)
			out['dry_mass'] = np.full_like(ppgpp, np.nan)
		return {k: v for k, v in out.items() if v.size == ppgpp.size}

	def _class_rows(self, base, idx, series, ppgpp, ctx):
		"""One row per gene class, for one aggregation."""
		mean_p = float(np.mean(ppgpp[idx]))
		f_ppgpp = mean_p ** 2 / (ctx['km_sq'] + mean_p ** 2)
		growth = max(float(fitting.interpolate_linearized_fit(
			mean_p, *ctx['growth_params'])), 0.0)
		tau_inf = np.log(2) / growth / 60 if growth > 0 else np.nan
		tau_used = ctx['frozen_tau'] if ctx['frozen_tau'] is not None else tau_inf
		loss = growth + ctx['deg']
		n_avg = ctx['get_n_avg'](tau_used, ctx['wt_coords'])
		expr = ctx['exp_free'] * (1 - f_ppgpp) + ctx['exp_ppgpp'] * f_ppgpp
		p_recon = normalize(expr * loss / n_avg)
		p_recon[ctx['is_rrna']] = p_recon[ctx['is_rrna']].mean()

		rows = []
		for name in CLASS_ORDER:
			m = ctx['masks'][name]
			if not m.any():
				continue
			take = lambda k: float(np.mean(series[k][name][idx]))
			copies, p_log = take('copies'), take('p_logged')
			share_t, init = take('share_target'), take('init')
			exp_copies = float(n_avg[m].sum())
			rows.append(dict(base, gene_class=name, n_tus=int(m.sum()),
				mean_replichore_fraction=float(ctx['frac'][m].mean()),
				expression=float(expr[m].sum()),
				loss_mean=float(loss[m].mean()),
				loss_growth_frac=float(np.mean(growth / loss[m])),
				loss_deg_frac=float(np.mean(ctx['deg'][m] / loss[m])),
				expected_copies=exp_copies,
				actual_copies=copies,
				copies_implied=share_t / p_log if p_log else float('nan'),
				a=copies / exp_copies if exp_copies else float('nan'),
				p_logged=p_log, p_recon=float(p_recon[m].sum()),
				p_after_attenuation=take('p_attn'),
				share_target=share_t, share_actual=take('share_actual'),
				initiations=init,
				init_per_copy=init / copies if copies else float('nan'),
				overcrowded=take('overcrowded')))
		return rows

	@staticmethod
	def _identities(p_arr, rsp, rnap, ppgpp, budget, good, ctx, variant, seed,
			gen):
		"""The six checks, on a subsample of timesteps."""
		idx = np.where(good)[0]
		if idx.size == 0:
			return []
		pick = idx[np.linspace(0, idx.size - 1, min(IDENTITY_SAMPLES,
			idx.size)).astype(int)]
		target = np.asarray(rsp.readColumn('target_rna_synth_prob'), float)
		copies = np.asarray(rsp.readColumn('promoter_copy_number'), float)
		init = np.asarray(rnap.readColumn('rnaInitEvent'), float)
		rows = []
		for i in pick:
			pp = float(ppgpp[i])
			f = pp ** 2 / (ctx['km_sq'] + pp ** 2)
			g = max(float(fitting.interpolate_linearized_fit(
				pp, *ctx['growth_params'])), 0.0)
			tau = np.log(2) / g / 60 if g > 0 else np.nan
			tu = ctx['frozen_tau'] if ctx['frozen_tau'] is not None else tau
			n_avg = ctx['get_n_avg'](tu, ctx['wt_coords'])
			expr = ctx['exp_free'] * (1 - f) + ctx['exp_ppgpp'] * f
			mine = normalize(expr * (g + ctx['deg']) / n_avg)
			mine[ctx['is_rrna']] = mine[ctx['is_rrna']].mean()
			ref = p_arr[i]
			nz = ref > 0
			id1 = float(np.abs(mine[nz] / ref[nz] - 1).max()) if nz.any() else np.nan
			w = copies[i] * ref
			pred = w / w.sum() if w.sum() > 0 else w
			tgt = target[i]
			nzt = tgt > 0
			id4 = float(np.abs(pred[nzt] / tgt[nzt] - 1).max()) if nzt.any() else np.nan
			rows.append(dict(variant=variant, seed=seed, generation=gen,
				timestep=int(i),
				id1_p_recon_vs_logged=id1,
				id2_sum_p=float(ref.sum()),
				id3_sum_target=float(tgt.sum()),
				id4_target_vs_copies_x_p=id4,
				id5_init_sum_vs_budget=float(init[i].sum() - budget[i]),
				id6_init_vs_budget_x_share=float(
					np.abs(init[i] - budget[i] * tgt).sum() / max(budget[i], 1))))
		return rows

	# ---- aggregation and output ------------------------------------------

	@staticmethod
	def _aggregate(class_rows, cell_rows):
		keys = ('expected_copies', 'actual_copies', 'copies_implied', 'a',
			'p_logged', 'p_recon', 'share_target', 'share_actual',
			'initiations', 'init_per_copy', 'expression', 'loss_mean',
			'loss_growth_frac', 'overcrowded')
		cell_keys = ('tau_realised', 'ppgpp_conc', 'tau_inferred', 'tau_used',
			'budget', 'active_rnap', 'total_rnap', 'active_ribosome',
			'total_ribosome', 'gfp_protein', 'dry_mass')
		cells = {}
		for r in cell_rows:
			cells.setdefault((r['variant'], r['generation'],
				r['aggregation']), []).append(r)
		out = []
		groups = {}
		for r in class_rows:
			groups.setdefault((r['variant'], r['generation'],
				r['aggregation'], r['gene_class']), []).append(r)
		for (v, g, agg, cls), rows in sorted(groups.items()):
			row = dict(variant=v, generation=g, aggregation=agg,
				gene_class=cls, n_lineages=len(rows),
				post_induction=rows[0]['post_induction'],
				n_tus=rows[0]['n_tus'])
			for k in keys:
				vals = [r[k] for r in rows if np.isfinite(r[k])]
				row[k] = float(np.mean(vals)) if vals else float('nan')
				row[k + '_sem'] = _sem(vals)
			for k in cell_keys:
				vals = [r[k] for r in cells.get((v, g, agg), [])
					if k in r and np.isfinite(r[k])]
				row[k] = float(np.mean(vals)) if vals else float('nan')
			out.append(row)
		return out

	# The origin-minus-terminus share transfer is the quantity that survives
	# normalisation, so it is computed here rather than left to the reader.
	@staticmethod
	def _transfer(by_gen):
		look = {(r['variant'], r['generation'], r['aggregation'],
			r['gene_class']): r for r in by_gen}
		out = []
		keys = sorted({(r['variant'], r['generation'], r['aggregation'])
			for r in by_gen})
		for v, g, agg in keys:
			o = look.get((v, g, agg, 'origin_ref'))
			t = look.get((v, g, agg, 'terminus_ref'))
			if not o or not t or not t['share_target']:
				continue
			out.append(dict(variant=v, generation=g, aggregation=agg,
				origin_share=o['share_target'], terminus_share=t['share_target'],
				origin_over_terminus=o['share_target'] / t['share_target'],
				origin_a=o['a'], terminus_a=t['a'],
				a_spread=1 - o['a'] / t['a'] if t['a'] else float('nan')))
		return out

	def _write(self, out_dir, name, cell_rows, class_rows, ident_rows, by_gen):
		for tag, rows in (('_cell', cell_rows), ('_class', class_rows),
				('_identity', ident_rows), ('_by_gen', by_gen),
				('_transfer', self._transfer(by_gen))):
			if not rows:
				continue
			path = os.path.join(out_dir, name + tag + '.csv')
			fields = list(rows[0].keys())
			with open(path, 'w', newline='') as h:
				w = csv.DictWriter(h, fieldnames=fields, extrasaction='ignore')
				w.writeheader()
				for r in rows:
					w.writerow(r)
			print('Wrote %s (%d rows)' % (path, len(rows)))

	@staticmethod
	def _print_identities(rows):
		if not rows:
			print('\nNo identity checks ran.')
			return False
		print('\nIDENTITY CHECKS, worst case over %d sampled timesteps'
			% len(rows))
		spec = [('id1_p_recon_vs_logged', 'recon p vs logged p', 1e-9, 'max rel'),
			('id2_sum_p', 'sum of p over TUs', None, 'want 1'),
			('id3_sum_target', 'sum of share over TUs', None, 'want 1'),
			('id4_target_vs_copies_x_p', 'share vs copies x p', 1e-6, 'max rel'),
			('id5_init_sum_vs_budget', 'init sum minus budget', None, 'want 0'),
			('id6_init_vs_budget_x_share', 'init vs budget x share', None,
				'multinomial noise')]
		gate = True
		for key, label, tol, note in spec:
			vals = np.array([r[key] for r in rows if np.isfinite(r[key])])
			if vals.size == 0:
				print('  %-24s no finite values' % label)
				continue
			if key in ('id2_sum_p', 'id3_sum_target'):
				worst = float(np.abs(vals - 1).max())
				print('  %-24s worst |x - 1| = %.3e   (%s)'
					% (label, worst, note))
			else:
				worst = float(np.abs(vals).max())
				verdict = ''
				if tol is not None:
					ok = worst < tol
					verdict = '   PASS' if ok else '   *** FAIL (tol %.0e)' % tol
					if key == 'id1_p_recon_vs_logged':
						gate = ok
				print('  %-24s worst %.3e   (%s)%s'
					% (label, worst, note, verdict))
		return gate

	@staticmethod
	def _print_summary(by_gen):
		rows = [r for r in by_gen if r['aggregation'] == 'time_mean']
		if not rows:
			return
		variants = sorted({r['variant'] for r in rows})
		gens = sorted({r['generation'] for r in rows})
		look = {(r['variant'], r['generation'], r['gene_class']): r
			for r in rows}
		print('\nrRNA operons, time-mean, by generation. Induction at gen %d.'
			% INDUCTION_GEN)
		print('%-5s %-4s %9s %9s %8s %9s %9s %8s'
			% ('v', 'gen', 'actual', 'expected', 'a', 'share', 'init', 'lin'))
		for v in variants:
			for g in gens:
				r = look.get((v, g, 'rrna'))
				if not r:
					continue
				print('%-5d %-4d %9.3f %9.3f %8.4f %9.5f %9.4f %8d'
					% (v, g, r['actual_copies'], r['expected_copies'], r['a'],
						r['share_target'], r['initiations'], r['n_lineages']))
			print()
		print('n_lineages falls where cells stalled; the survivors are the')
		print('healthier ones, so read the means against that column.')

	@staticmethod
	def _plot(out_dir, name, by_gen):
		rows = [r for r in by_gen if r['aggregation'] == 'time_mean']
		if not rows:
			return
		variants = sorted({r['variant'] for r in rows})
		gens = sorted({r['generation'] for r in rows})
		look = {(r['variant'], r['generation'], r['gene_class']): r
			for r in rows}
		panels = [('actual_copies', 'actual copies'),
			('expected_copies', 'expected copies'),
			('a', 'a = actual / expected'),
			('share_target', 'share of initiation'),
			('initiations', 'initiation events'),
			('init_per_copy', 'initiations per copy')]
		fig, axes = plt.subplots(len(panels), 1,
			figsize=(9, 1.7 * len(panels)), sharex=True)
		cmap = plt.get_cmap('viridis')
		shade = {v: cmap(0.1 + 0.75 * i / max(1, len(variants) - 1))
			for i, v in enumerate(variants)}
		for ax, (key, label) in zip(np.atleast_1d(axes), panels):
			for v in variants:
				xs = [g for g in gens if (v, g, 'rrna') in look]
				ys = [look[(v, g, 'rrna')][key] for g in xs]
				pre = [y for x, y in zip(xs, ys)
					if x < INDUCTION_GEN and np.isfinite(y)]
				base = np.mean(pre) if pre else np.nan
				ax.plot(xs, np.array(ys) / base if base else ys, marker='o',
					ms=3, lw=1.3, color=shade[v], label='v%d' % v)
			ax.axvline(INDUCTION_GEN - 0.5, color='C3', lw=1.3)
			ax.axhline(1.0, color='0.7', lw=0.7, ls='--')
			ax.set_ylabel(label, fontsize='x-small')
			ax.tick_params(labelsize='x-small')
		np.atleast_1d(axes)[0].legend(fontsize='xx-small', ncol=4,
			frameon=False)
		np.atleast_1d(axes)[0].set_title('rRNA operons, relative to '
			'pre-induction. Red line: induction.', fontsize='small')
		np.atleast_1d(axes)[-1].set_xlabel('generation', fontsize='small')
		plt.tight_layout()
		exportFigure(plt, out_dir, name, metadata=None)
		plt.close('all')


if __name__ == '__main__':
	Plot().cli()
