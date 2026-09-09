"""
Transcription-factor promoter occupancy under new-gene burden, and what the
freed sites actually buy.

Motivation
----------
A lab member running antibiotic simulations saw LexA protein counts fall and
inferred: fewer LexA -> fewer occupied binding sites -> de-repression -> more
expression of LexA targets. The first two steps are exactly how the model
works. The third does not follow, and this script measures the gap.

LexA is a "0CS" transcription factor, so TfBinding hard-codes
pPromoterBound = 1.0 and its binding reduces to

	bound sites = min(available promoter sites, free LexA dimers)

with no Kd, no ligand and no residence time -- TfUnbinding frees every TF each
timestep and TfBinding re-draws from scratch. Bound state then enters
transcription additively per promoter copy,

	promoter_init_probs = basal_prob + ppgpp_scale * (delta_prob . bound_TF)

after which the whole vector is CLIPPED AT ZERO AND NORMALISED TO SUM TO 1
(transcript_initiation.py:177-178). That normalisation is the crux: losing a
repressor buys a larger *share* of a fixed initiation budget, not additional
absolute transcription. If burden shrinks the budget faster than de-repression
grows the share, realised transcription of the regulon falls even though every
target's synthesis probability rose.

So the quantity of interest is a product of two separately measurable factors:

	regulon initiations  =  regulon share of budget  x  total initiation budget
	     (net)                  (de-repression)            (RNAP scarcity)

That is an identity, not a fit, because initiations are a multinomial draw over
the normalised probability vector. Both factors are reported, at pooled grain
where the identity is exact and at per-lineage grain where it carries a CI.

Why the estimator is not an event count
---------------------------------------
Realised initiations for the LexA regulon are about 89 events per generation
across all 50 target transcription units, and 1-10 per gene. Most LexA targets
are deeply subgenerational. A single-seed comparison of raw counts returns 89
against 89 -- pure Poisson noise -- while the probability-weighted estimator

	expected initiations = sum_t actual_rna_synth_prob[t, TU] * didInitialize[t]

resolves the same contrast at 96.6 against 80.5. It is the conditional
expectation of the multinomial draw given the probability vector, so it is
Rao-Blackwellised: same mean, far lower variance. Both are emitted; the
expected one is the estimator, the realised one is the honesty check.

All 23 transcription factors, not just LexA
-------------------------------------------
The occupancy read is identical for every TF, so it costs nothing to cover all
of them -- and LexA is not the largest exposure. Shares of the basal
transcription budget, measured from the shipped sim_data:

	IHF 5.45%   Fis 4.66%   Lrp 3.94%   H-NS 3.07%
	ArcA 2.12%  CRP 2.09%   FNR 1.93%   ...   LexA 0.39%

IHF, Fis, H-NS and LeuO are 0CS, so they sit under exactly the same
min(sites, TFs) cap as LexA and together hold 13.2% of the budget -- 34x
LexA's share. CRP, FNR and ArcA are mostly activators, so losing THEIR
occupancy pushes expression down rather than up. Whichever way the LexA result
lands, that table is where a larger effect would hide.

What is NOT here
----------------
There is no SOS response in this model: no RecA-mediated LexA autocleavage, no
ssDNA or DNA-damage species, no damage sensing, and no stress condition in
condition_defs.tsv. LexA's only "inactive" state in the reconstruction is a
genetic perturbation zeroing EG10533_RNA. Any LexA movement measured here is a
dilution/synthesis consequence of altered growth, and must not be reported as
SOS induction.

A readout trap, recorded here because it is easy to hit
------------------------------------------------------
MonomerCounts and BulkMolecules both report the FREE pool only. Bound dimers
are decremented out of the bulk count (tf_binding.py:160) and live solely as a
bound_TF flag on the promoter unique-molecule. In minimal media the free LexA
pool is ~0.4 while ~82 dimers are bound and invisible, so

	total LexA monomer equivalents = MonomerCounts[PD00205] + 2 * nActualBound

A fall in the free pool alone is equally consistent with MORE binding, since
anything that raises chromosome copy number raises site count and sequesters
more LexA. Both the free and the reconstructed total are emitted, per TF.

Output
------
	<name>_tf_occupancy.csv    (variant, seed, generation, tf)
	<name>_decomposition.csv   (variant, seed, generation) -- the identity
	<name>_per_gene.csv        (variant, seed, generation, LexA target TU)
	<name>_gates.csv           the mechanism checks
	<name>_manifest.txt        coverage, gates, dropped-cell census
"""

import csv
import json
import os
import pickle
from collections import Counter

from matplotlib import pyplot as plt
import numpy as np

from models.ecoli.analysis import variantAnalysisPlot
from wholecell.analysis.analysis_tools import (exportFigure,
	first_cell_with_table, read_bulk_molecule_counts)
from wholecell.io.tablereader import TableReader

# LexA. The TF the model binds is the DIMER (2 PD00205 -> 1 PC00010,
# complexation_reactions.tsv:1023); tf_binding only ever views the dimer.
LEXA_TF_ID = 'PC00010'
LEXA_MONOMER_ID = 'PD00205[c]'

# Fallback when metadata carries no induction generation. The burden-ladder and
# internal-shift variants both import NEW_GENE_INDUCTION_GEN from
# new_gene_internal_shift, but the value has differed between batches, so it is
# read from the run's own metadata and only falls back to this.
DEFAULT_INDUCTION_GEN = 8

# Floating tolerance for the decomposition identity. It is exact arithmetic on
# pooled sums, so this is generous.
IDENTITY_TOL = 1e-9


class Plot(variantAnalysisPlot.VariantAnalysisPlot):

	def do_plot(self, inputDir, plotOutDir, plotOutFileName, simDataFile,
			validationDataFile, metadata):

		self._dropped = Counter()

		with open(simDataFile, 'rb') as handle:
			sim_data = pickle.load(handle)

		induction_gen = self._induction_gen(metadata)
		variants = self.ap.get_variants()
		print('Variants found: %s' % (variants,))
		print('Induction generation: %d' % induction_gen)

		probe = first_cell_with_table(self.ap.get_cells(), 'RnaSynthProb')
		if probe is None:
			print('No cell has a readable RnaSynthProb table. Nothing to do.')
			return
		ctx = self._context(sim_data, probe)
		self._report_tf_table(ctx)

		tf_rows = []
		decomp_rows = []
		gene_rows = []
		gate_rows = []

		for variant in variants:
			for path, seed, gen in self._cells(variant):
				got = self._one_cell(path, variant, seed, gen, ctx)
				if got is None:
					continue
				tf_rows.extend(got['tf'])
				decomp_rows.append(got['decomp'])
				gene_rows.extend(got['gene'])
				gate_rows.append(got['gate'])

		if not decomp_rows:
			print('No readable cells. Nothing written.')
			return

		name = plotOutFileName
		self._write(os.path.join(plotOutDir, name + '_tf_occupancy.csv'),
			tf_rows)
		self._write(os.path.join(plotOutDir, name + '_decomposition.csv'),
			decomp_rows)
		self._write(os.path.join(plotOutDir, name + '_per_gene.csv'),
			gene_rows)
		self._write(os.path.join(plotOutDir, name + '_gates.csv'), gate_rows)

		gates = self._check_gates(gate_rows, decomp_rows, induction_gen)
		self._manifest(os.path.join(plotOutDir, name + '_manifest.txt'),
			decomp_rows, tf_rows, gene_rows, gates, ctx, induction_gen)
		self._figure(plotOutDir, name, decomp_rows, tf_rows, induction_gen)

	# ---- setup -----------------------------------------------------------

	@staticmethod
	def _induction_gen(metadata):
		"""
		Induction generation from the run's own metadata.

		Read rather than hardcoded: the value has differed between batches
		(1 in the local gfp_shift_minimal regression run, 8 in the burden
		ladder), and getting it wrong silently mislabels the pre-induction
		baseline as burdened.
		"""
		for key in ('new_gene_induction_gen', 'induction_gen'):
			value = (metadata or {}).get(key)
			if value is not None:
				try:
					return int(value)
				except (TypeError, ValueError):
					pass
		print('WARNING: metadata carries no new_gene_induction_gen; assuming '
			'%d.' % DEFAULT_INDUCTION_GEN)
		return DEFAULT_INDUCTION_GEN

	@staticmethod
	def _context(sim_data, probe):
		"""Static per-run indexing, built once."""
		reg = sim_data.process.transcription_regulation
		transcription = sim_data.process.transcription
		common = sim_data.common_names

		rsp = TableReader(os.path.join(probe, 'simOut', 'RnaSynthProb'))
		rna_ids = list(rsp.readAttribute('rnaIds'))
		tf_ids = list(rsp.readAttribute('tf_ids'))
		n_tu = len(rna_ids)
		n_tf = len(tf_ids)

		delta = reg.delta_prob
		basal = np.asarray(reg.basal_prob, dtype=float)

		# Derive each TF's targets exactly as TfBinding does, so the TUs whose
		# fitted deltaV is zero are still counted as available promoter sites
		# (they are, in the process) while contributing no effect.
		targets = {}
		deltas = {}
		for j, tf in enumerate(tf_ids):
			sel = delta['deltaJ'] == j
			tu = np.asarray(delta['deltaI'][sel], dtype=int)
			dv = np.asarray(delta['deltaV'][sel], dtype=float)
			order = np.argsort(tu)
			targets[tf] = tu[order]
			deltas[tf] = dv[order]

		lexa_tu = targets.get(LEXA_TF_ID, np.array([], dtype=int))
		lexa_dv = deltas.get(LEXA_TF_ID, np.array([], dtype=float))

		def tu_label(index):
			rna = rna_ids[index]
			bare = rna.split('[')[0]
			label = str(common.get_common_name(bare.replace('_RNA', '')))
			if label == bare.replace('_RNA', ''):
				label = str(common.get_common_name(bare))
			return label

		# LexA target TU -> its monomers, so protein can be reported per gene.
		monomer_data = sim_data.process.translation.monomer_data
		cistron_to_monomers = {}
		for row in monomer_data:
			cistron_to_monomers.setdefault(row['cistron_id'], []).append(
				row['id'])
		cistron_ids = list(transcription.cistron_data['id'])

		lexa_monomers = {}
		for index in lexa_tu:
			mons = []
			try:
				for ci in transcription.rna_id_to_cistron_indexes(
						rna_ids[index]):
					mons.extend(cistron_to_monomers.get(cistron_ids[ci], []))
			except Exception:  # noqa: BLE001 - a TU with no cistron mapping
				pass
			lexa_monomers[int(index)] = mons

		# Flat column indices for n_bound_TF_per_TU, which is stored as
		# (t, n_TU * n_TF) in C order. readSubcolumn is broken for it (its
		# declared label list is n_TU long against a width of n_TU * n_TF), and
		# reshaping the full column holds ~430 MB per cell, so the LexA slice is
		# pulled by index instead. Verified equal to the reshape, and its
		# per-timestep sum equals nActualBound.
		lexa_tf_index = tf_ids.index(LEXA_TF_ID) if LEXA_TF_ID in tf_ids \
			else None
		lexa_bound_flat = (lexa_tu * n_tf + lexa_tf_index
			if lexa_tf_index is not None else np.array([], dtype=int))

		return dict(
			rna_ids=rna_ids, tf_ids=tf_ids, n_tu=n_tu, n_tf=n_tf,
			targets=targets, deltas=deltas, basal=basal,
			tf_to_gene=reg.tf_to_gene_id, tf_type=reg.tf_to_tf_type,
			active_to_bound=reg.active_to_bound,
			lexa_tu=lexa_tu, lexa_dv=lexa_dv, lexa_tf_index=lexa_tf_index,
			lexa_bound_flat=np.asarray(lexa_bound_flat, dtype=int),
			lexa_monomers=lexa_monomers,
			tu_label={int(i): tu_label(int(i)) for i in lexa_tu},
			basal_total=float(basal.sum()))

	@staticmethod
	def _report_tf_table(ctx):
		print('')
		print('%-20s %-8s %-6s %7s %12s' % ('TF', 'gene', 'type', 'n_TU',
			'basal_share'))
		rows = []
		for tf in ctx['tf_ids']:
			tu = ctx['targets'][tf]
			share = (ctx['basal'][tu].sum() / ctx['basal_total'] * 100
				if ctx['basal_total'] > 0 else float('nan'))
			rows.append((share, tf, tu))
		for share, tf, tu in sorted(rows, reverse=True):
			print('%-20s %-8s %-6s %7d %11.2f%%'
				% (tf, ctx['tf_to_gene'].get(tf, '?'),
					ctx['tf_type'].get(tf, '?'), len(tu), share))
		print('')
		if LEXA_TF_ID not in ctx['tf_ids']:
			print('WARNING: LexA (%s) is not in tf_ids for this run. The '
				'decomposition and per-gene tables will be empty.'
				% LEXA_TF_ID)

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
				# divided. This is survivorship-biased -- the lineages that
				# stall are the burdened ones -- so n_lineages is reported per
				# generation in the manifest rather than assumed constant.
				if gen + 1 not in gens:
					self._dropped['incomplete'] += 1
					continue
				out.append((path, seed, gen))
		return out

	# ---- per cell --------------------------------------------------------

	def _one_cell(self, path, variant, seed, gen, ctx):
		sim_out = os.path.join(path, 'simOut')
		try:
			time = TableReader(os.path.join(sim_out, 'Main')
				).readColumn('time')
			rsp = TableReader(os.path.join(sim_out, 'RnaSynthProb'))
			rnap = TableReader(os.path.join(sim_out, 'RnapData'))
			mass = TableReader(os.path.join(sim_out, 'Mass'))
			umc = TableReader(os.path.join(sim_out, 'UniqueMoleculeCounts'))

			bound = np.asarray(rsp.readColumn('nActualBound', squeeze=False),
				dtype=float)
			sites = np.asarray(
				rsp.readColumn('n_available_promoters', squeeze=False),
				dtype=float)
			p_bound = np.asarray(rsp.readColumn('pPromoterBound',
				squeeze=False), dtype=float)
			budget = np.asarray(rnap.readColumn('didInitialize'), dtype=float)
			dry_mass = np.asarray(mass.readColumn('dryMass'), dtype=float)
		except Exception as exc:  # noqa: BLE001 - preempted/empty simOut
			self._dropped['unreadable'] += 1
			print('skip v%d s%d g%d: %s' % (variant, seed, gen, exc))
			return None

		# When didInitialize == 0, transcript_initiation returns early after
		# writing target_rna_synth_prob but BEFORE max_p,
		# actual_rna_synth_prob, tu_is_overcrowded, rnaInitEvent and
		# didInitialize, so those five carry the PREVIOUS timestep's values.
		# Everything derived from them is masked to the timesteps that ran.
		good = budget > 0
		if not good.any():
			self._dropped['no_initiation'] += 1
			return None

		free_active = self._free_active(sim_out, ctx)
		total_tf = free_active + bound

		# --- per TF ---
		tf_out = []
		share_actual_by_tf = {}
		share_target_by_tf = {}
		actual = np.asarray(rsp.readColumn('actual_rna_synth_prob',
			squeeze=False), dtype=float)
		target = np.asarray(rsp.readColumn('target_rna_synth_prob',
			squeeze=False), dtype=float)
		for j, tf in enumerate(ctx['tf_ids']):
			tu = ctx['targets'][tf]
			sa = actual[:, tu].sum(axis=1)
			st = target[:, tu].sum(axis=1)
			share_actual_by_tf[tf] = sa
			share_target_by_tf[tf] = st
			has_sites = sites[:, j] > 0
			tf_out.append(dict(
				variant=variant, seed=seed, generation=gen,
				tf=tf, gene=ctx['tf_to_gene'].get(tf, ''),
				tf_type=ctx['tf_type'].get(tf, ''),
				n_target_tu=int(len(tu)),
				sites=float(sites[:, j].mean()),
				bound=float(bound[:, j].mean()),
				occupancy=float(np.divide(
					bound[has_sites, j], sites[has_sites, j]).mean())
					if has_sites.any() else float('nan'),
				free_active=float(free_active[:, j].mean()),
				total_active=float(total_tf[:, j].mean()),
				p_promoter_bound=float(p_bound[:, j].mean()),
				# TF-limited means the free pool is exhausted: every copy is on
				# DNA. Site-limited means every available site is occupied.
				# Both can hold at once when the two counts coincide.
				frac_tf_limited=float(
					(bound[:, j] >= total_tf[:, j] - 0.5).mean()),
				frac_site_limited=float(
					(bound[:, j] >= sites[:, j] - 0.5).mean()),
				share_actual=float(sa[good].mean()),
				share_target=float(st[good].mean()),
				basal_share=float(ctx['basal'][tu].sum() / ctx['basal_total'])
					if ctx['basal_total'] > 0 else float('nan')))

		# --- the decomposition, LexA ---
		decomp = self._decompose(rsp, rnap, umc, mass, ctx, good, budget,
			bound, sites, free_active, total_tf, share_actual_by_tf,
			share_target_by_tf, time, dry_mass, variant, seed, gen, sim_out)

		# --- per gene ---
		gene_out = self._per_gene(rsp, rnap, sim_out, ctx, good, budget,
			actual, dry_mass, variant, seed, gen)

		# --- mechanism gates ---
		j = ctx['lexa_tf_index']
		if j is None:
			gate = dict(variant=variant, seed=seed, generation=gen,
				p_bound_is_binary=True, cap_respected=True,
				frac_stale_timesteps=float((~good).mean()))
		else:
			uniq = np.unique(p_bound[:, j])
			gate = dict(
				variant=variant, seed=seed, generation=gen,
				# 0CS means pPromoterBound is hard-coded to 1.0, and 0.0 only
				# when the TF count is zero. Anything else means the TF was
				# reclassified and the whole premise changes.
				p_bound_is_binary=bool(
					np.all(np.isin(np.round(uniq, 12), [0.0, 1.0]))),
				# bound = min(sites, free dimers) -- the cap that makes the
				# mechanism work at all.
				cap_respected=bool(np.all(
					bound[:, j] <= np.minimum(sites[:, j], total_tf[:, j])
					+ 0.5)),
				frac_stale_timesteps=float((~good).mean()))

		return dict(tf=tf_out, decomp=decomp, gene=gene_out, gate=gate)

	@staticmethod
	def _free_active(sim_out, ctx):
		"""
		Free (unbound) active-TF counts, one column per TF.

		This is the FREE pool: TfBinding decrements bound copies out of the
		bulk count, so DNA-bound TF is invisible here. The reconstructed total
		is free + nActualBound, which is what total_active reports.
		"""
		ids = ['%s[c]' % tf for tf in ctx['tf_ids']]
		try:
			counts, = read_bulk_molecule_counts(sim_out, (ids,))
			return np.asarray(counts, dtype=float).reshape(
				-1, len(ctx['tf_ids']))
		except Exception:  # noqa: BLE001 - missing bulk ids
			reader = TableReader(os.path.join(sim_out, 'BulkMolecules'))
			names = list(reader.readAttribute('objectNames'))
			out = []
			for mol in ids:
				if mol in names:
					out.append(reader.readColumn(
						'counts', [names.index(mol)], squeeze=False)[:, 0])
				else:
					out.append(np.full(len(names) and 0 or 0, np.nan))
			return np.asarray(out, dtype=float).T

	def _decompose(self, rsp, rnap, umc, mass, ctx, good, budget, bound,
			sites, free_active, total_tf, share_actual_by_tf,
			share_target_by_tf, time, dry_mass, variant, seed, gen, sim_out):
		"""
		The identity, per cell.

		expected_inits = sum_t share[t] * budget[t], and pooled_share is
		defined as expected_inits / total_budget so that

			expected_inits == pooled_share * total_budget

		holds exactly. Ratios between variants are formed downstream from the
		pooled sums, where the product identity therefore also holds exactly.
		"""
		j = ctx['lexa_tf_index']
		share = (share_actual_by_tf[LEXA_TF_ID] if j is not None
			else np.zeros_like(budget))
		share_t = (share_target_by_tf[LEXA_TF_ID] if j is not None
			else np.zeros_like(budget))

		expected = float((share[good] * budget[good]).sum())
		total_budget = float(budget[good].sum())

		init = np.asarray(rnap.readColumn('rnaInitEvent', squeeze=False),
			dtype=float)
		realised = float(init[np.ix_(good, ctx['lexa_tu'])].sum()) \
			if len(ctx['lexa_tu']) else float('nan')
		total_realised = float(init[good].sum())

		unique_ids = list(umc.readAttribute('uniqueMoleculeIds'))
		unique = np.asarray(umc.readColumn('uniqueMoleculeCounts'),
			dtype=float)

		def unique_mean(key):
			return (float(unique[:, unique_ids.index(key)].mean())
				if key in unique_ids else float('nan'))

		try:
			ppgpp = float(np.nanmean(TableReader(
				os.path.join(sim_out, 'GrowthLimits')
				).readColumn('ppgpp_conc')))
		except Exception:  # noqa: BLE001 - column absent in older vintages
			ppgpp = float('nan')

		lexa_free_monomer = self._monomer_mean(sim_out, [LEXA_MONOMER_ID])

		return dict(
			variant=variant, seed=seed, generation=gen,
			doubling_time_min=float((time[-1] - time[0]) / 60.),
			dry_mass_fg=float(dry_mass.mean()),
			protein_mass_fg=float(np.asarray(
				mass.readColumn('proteinMass'), dtype=float).mean()),
			ppgpp_conc=ppgpp,
			active_rnap=unique_mean('active_RNAP'),
			active_ribosome=unique_mean('active_ribosome'),
			# de-repression factor
			lexa_sites=float(sites[:, j].mean()) if j is not None
				else float('nan'),
			lexa_bound=float(bound[:, j].mean()) if j is not None
				else float('nan'),
			lexa_free_sites=float((sites[:, j] - bound[:, j]).mean())
				if j is not None else float('nan'),
			lexa_occupancy=float(np.divide(
				bound[:, j], np.maximum(sites[:, j], 1)).mean())
				if j is not None else float('nan'),
			lexa_free_dimers=float(free_active[:, j].mean())
				if j is not None else float('nan'),
			lexa_total_dimers=float(total_tf[:, j].mean())
				if j is not None else float('nan'),
			lexa_free_monomer=lexa_free_monomer,
			# The trap: MonomerCounts omits DNA-bound LexA entirely.
			lexa_total_monomer_equivalents=(
				lexa_free_monomer + 2. * float(bound[:, j].mean())
				if j is not None else float('nan')),
			pooled_share=expected / total_budget if total_budget > 0
				else float('nan'),
			pooled_share_target=float(
				(share_t[good] * budget[good]).sum() / total_budget)
				if total_budget > 0 else float('nan'),
			# throughput factor
			total_budget=total_budget,
			total_realised_inits=total_realised,
			# net
			expected_regulon_inits=expected,
			realised_regulon_inits=realised)

	@staticmethod
	def _monomer_mean(sim_out, ids):
		try:
			reader = TableReader(os.path.join(sim_out, 'MonomerCounts'))
			names = list(reader.readAttribute('monomerIds'))
			idx = [names.index(m) for m in ids if m in names]
			if not idx:
				return float('nan')
			return float(reader.readColumn(
				'monomerCounts', idx, squeeze=False).sum(axis=1).mean())
		except Exception:  # noqa: BLE001
			return float('nan')

	def _per_gene(self, rsp, rnap, sim_out, ctx, good, budget, actual,
			dry_mass, variant, seed, gen):
		"""One row per LexA target transcription unit."""
		tu = ctx['lexa_tu']
		if not len(tu):
			return []

		bound_per_tu = np.asarray(rsp.readColumn(
			'n_bound_TF_per_TU', ctx['lexa_bound_flat'], squeeze=False),
			dtype=float)
		copies = np.asarray(rsp.readColumn('promoter_copy_number',
			squeeze=False), dtype=float)[:, tu]
		init = np.asarray(rnap.readColumn('rnaInitEvent', squeeze=False),
			dtype=float)[:, tu]

		try:
			mono = TableReader(os.path.join(sim_out, 'MonomerCounts'))
			mono_names = list(mono.readAttribute('monomerIds'))
			mono_counts = np.asarray(mono.readColumn('monomerCounts'),
				dtype=float)
		except Exception:  # noqa: BLE001
			mono_names, mono_counts = [], None

		mass_mean = float(dry_mass.mean())
		out = []
		for k, index in enumerate(tu):
			index = int(index)
			share = actual[:, index]
			mons = ctx['lexa_monomers'].get(index, [])
			idx = [mono_names.index(m) for m in mons if m in mono_names]
			protein = (float(mono_counts[:, idx].sum(axis=1).mean())
				if idx and mono_counts is not None else float('nan'))
			out.append(dict(
				variant=variant, seed=seed, generation=gen,
				tu_id=ctx['rna_ids'][index], gene=ctx['tu_label'][index],
				basal_prob=float(ctx['basal'][index]),
				delta_v=float(ctx['lexa_dv'][k]),
				# For 13 of these TUs deltaV == -basal_prob exactly, so a bound
				# promoter copy is driven to ZERO initiation and unbinding is
				# an on/off switch rather than a graded change.
				derepression_ceiling=float(
					abs(ctx['lexa_dv'][k]) / ctx['basal'][index])
					if ctx['basal'][index] > 0 else float('nan'),
				bound_copies=float(bound_per_tu[:, k].mean()),
				promoter_copies=float(copies[:, k].mean()),
				frac_copies_bound=float(np.divide(
					bound_per_tu[:, k],
					np.maximum(copies[:, k], 1)).mean()),
				share_actual=float(share[good].mean()),
				expected_inits=float((share[good] * budget[good]).sum()),
				realised_inits=float(init[good, k].sum()),
				protein_counts=protein,
				protein_per_fg=protein / mass_mean if mass_mean > 0
					else float('nan')))
		return out

	# ---- gates, manifest, figure ----------------------------------------

	@staticmethod
	def _check_gates(gate_rows, decomp_rows, induction_gen):
		gates = {}
		gates['p_promoter_bound_binary'] = all(
			r['p_bound_is_binary'] for r in gate_rows)
		gates['binding_cap_respected'] = all(
			r['cap_respected'] for r in gate_rows)

		# The identity, on pooled sums, per variant.
		worst = 0.0
		by_variant = {}
		for row in decomp_rows:
			acc = by_variant.setdefault(row['variant'], [0.0, 0.0, 0.0])
			acc[0] += row['expected_regulon_inits']
			acc[1] += row['total_budget']
			acc[2] += row['pooled_share'] * row['total_budget']
		for variant, (expected, total, reconstructed) in by_variant.items():
			denom = max(abs(expected), 1e-30)
			worst = max(worst, abs(expected - reconstructed) / denom)
		gates['identity_max_rel_error'] = worst
		gates['identity_holds'] = worst < IDENTITY_TOL

		# Pre-induction generations must be indistinguishable across variants:
		# every variant runs the same baseline until induction, so a spread
		# there means the variants differ for a reason this analysis has not
		# accounted for.
		pre = [r for r in decomp_rows if r['generation'] < induction_gen]
		spread = float('nan')
		if pre:
			means = {}
			for row in pre:
				means.setdefault(row['variant'], []).append(
					row['lexa_occupancy'])
			vals = [np.nanmean(v) for v in means.values()]
			if len(vals) > 1 and np.isfinite(np.nanmean(vals)):
				spread = float((np.nanmax(vals) - np.nanmin(vals))
					/ max(abs(np.nanmean(vals)), 1e-30))
		gates['pre_induction_occupancy_spread'] = spread
		gates['pre_induction_flat'] = (
			bool(spread < 0.05) if np.isfinite(spread) else False)

		gates['max_frac_stale_timesteps'] = float(max(
			(r['frac_stale_timesteps'] for r in gate_rows), default=0.0))
		return gates

	def _manifest(self, path, decomp_rows, tf_rows, gene_rows, gates, ctx,
			induction_gen):
		lines = []
		add = lines.append
		add('TF DE-REPRESSION UNDER NEW-GENE BURDEN')
		add('')
		add('INDUCTION GENERATION: %d  (generations < this are the '
			'pre-induction baseline)' % induction_gen)
		add('')
		add('SEED COVERAGE')
		by_variant = {}
		for row in decomp_rows:
			by_variant.setdefault(row['variant'], set()).add(row['seed'])
		gens = {}
		for row in decomp_rows:
			gens.setdefault(row['variant'], set()).add(row['generation'])
		for variant in sorted(by_variant):
			add('  variant %-3d  n_seeds=%-4d generations %d-%d  '
				'n_lineage_gens=%d'
				% (variant, len(by_variant[variant]),
					min(gens[variant]), max(gens[variant]),
					sum(1 for r in decomp_rows
						if r['variant'] == variant)))
		add('')
		add('DROPPED CELLS')
		if self._dropped:
			for key, count in sorted(self._dropped.items()):
				add('  %-16s %d' % (key, count))
		else:
			add('  none')
		add('')
		add('GATES')
		add('  pPromoterBound[lexA] in {0,1} (0CS path)     %s'
			% ('PASS' if gates['p_promoter_bound_binary'] else 'FAIL'))
		add('  bound <= min(sites, dimers) (the cap)        %s'
			% ('PASS' if gates['binding_cap_respected'] else 'FAIL'))
		add('  share x budget == expected inits            %s  '
			'(max rel err %.2e)'
			% ('PASS' if gates['identity_holds'] else 'FAIL',
				gates['identity_max_rel_error']))
		add('  pre-induction occupancy flat across variants %s  '
			'(spread %.3f)'
			% ('PASS' if gates['pre_induction_flat'] else 'FAIL',
				gates['pre_induction_occupancy_spread']))
		add('  max fraction of stale timesteps              %.5f'
			% gates['max_frac_stale_timesteps'])
		add('')
		add('If either of the first two gates fails, the mechanism described in')
		add('this script no longer holds and nothing downstream is valid.')
		add('')
		add('TF BUDGET SHARES (static, from sim_data)')
		add('  %-20s %-8s %-6s %7s %12s'
			% ('TF', 'gene', 'type', 'n_TU', 'basal_share'))
		rows = []
		for tf in ctx['tf_ids']:
			tu = ctx['targets'][tf]
			rows.append((ctx['basal'][tu].sum() / ctx['basal_total'], tf, tu))
		for share, tf, tu in sorted(rows, reverse=True):
			add('  %-20s %-8s %-6s %7d %11.2f%%'
				% (tf, ctx['tf_to_gene'].get(tf, '?'),
					ctx['tf_type'].get(tf, '?'), len(tu), share * 100))
		add('')
		add('ROW COUNTS')
		add('  tf_occupancy   %d' % len(tf_rows))
		add('  decomposition  %d' % len(decomp_rows))
		add('  per_gene       %d' % len(gene_rows))
		with open(path, 'w') as handle:
			handle.write('\n'.join(lines) + '\n')
		print('')
		print('\n'.join(lines))

	@staticmethod
	def _write(path, rows):
		if not rows:
			print('nothing to write to %s' % path)
			return
		keys = list(rows[0].keys())
		with open(path, 'w') as handle:
			writer = csv.DictWriter(handle, fieldnames=keys)
			writer.writeheader()
			writer.writerows(rows)
		print('wrote %s (%d rows)' % (path, len(rows)))

	@staticmethod
	def _figure(out_dir, name, decomp_rows, tf_rows, induction_gen):
		"""
		Diagnostic only -- the publishable figures are built downstream from
		the CSVs. Panels are per generation so the transient around induction
		is visible, which is the point: LexA dimers are effectively never
		degraded (ProteinDegradation views only the monomer, whose free pool is
		~0.4), so LexA falls by dilution alone and responds with a lag.
		"""
		variants = sorted({r['variant'] for r in decomp_rows})
		gens = sorted({r['generation'] for r in decomp_rows})
		panels = [
			('lexa_occupancy', 'LexA occupancy'),
			('lexa_free_sites', 'unoccupied LexA sites'),
			('lexa_total_dimers', 'total LexA dimers'),
			('pooled_share', 'regulon share of budget'),
			('total_budget', 'total initiation budget'),
			('expected_regulon_inits', 'regulon initiations (expected)'),
			('active_rnap', 'active RNAP'),
			('doubling_time_min', 'doubling time (min)')]

		fig, axes = plt.subplots(len(panels), 1,
			figsize=(9, 1.6 * len(panels)), sharex=True)
		cmap = plt.get_cmap('viridis')
		shade = {v: cmap(0.1 + 0.75 * i / max(1, len(variants) - 1))
			for i, v in enumerate(variants)}
		for ax, (key, label) in zip(np.atleast_1d(axes), panels):
			for variant in variants:
				xs, ys = [], []
				for gen in gens:
					vals = [r[key] for r in decomp_rows
						if r['variant'] == variant
						and r['generation'] == gen
						and np.isfinite(r[key])]
					if vals:
						xs.append(gen)
						ys.append(np.mean(vals))
				if xs:
					ax.plot(xs, ys, marker='o', ms=3, lw=1.3,
						color=shade[variant], label='v%d' % variant)
			ax.axvline(induction_gen - 0.5, color='C3', lw=1.3)
			ax.set_ylabel(label, fontsize='xx-small')
			ax.tick_params(labelsize='xx-small')
		np.atleast_1d(axes)[0].legend(fontsize='xx-small', ncol=8,
			frameon=False)
		np.atleast_1d(axes)[0].set_title(
			'LexA de-repression and what it buys. Red line: induction.',
			fontsize='small')
		np.atleast_1d(axes)[-1].set_xlabel('generation', fontsize='small')
		plt.tight_layout()
		exportFigure(plt, out_dir, name, metadata=None)
		plt.close('all')


if __name__ == '__main__':
	Plot().cli()
