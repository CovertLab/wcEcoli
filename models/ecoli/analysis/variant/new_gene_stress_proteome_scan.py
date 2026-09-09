"""
Do any stress-response proteins move under new-gene burden, and does the
movement mean anything?

What is actually modelled
-------------------------
Only the stringent response. That is not a caveat, it is the headline, and it
was established by an exhaustive pass over models/ecoli/processes/, the
reconstruction dataclasses and the flat files:

  ppGpp / stringent     YES, fully mechanistic. RelA activated by ribosomes
                        with uncharged tRNA in the A site, SpoT synthesis plus
                        uncharged-tRNA-inhibited hydrolysis
                        (polypeptide_elongation.py:786,865-880). ppGpp then
                        rescales BOTH the entire basal synthesis-probability
                        vector and the active RNAP fraction
                        (transcript_initiation.py:129-152), over 393 curated
                        genes. On by default.
  tRNA attenuation      YES, 48 operons, responsive to charged-tRNA levels.
  SOS                   NO. LexA is present as a 0CS repressor but there is no
                        RecA-mediated autocleavage, no ssDNA or DNA-damage
                        species and no damage sensing anywhere.
  heat shock            NO. Zero hits for rpoH/dnaK/groL/ibpA in any Python
                        code, and none of the 16 genes is even ppGpp-linked.
                        There is no misfolded-protein species, so there is no
                        substrate for a heat-shock loop even in principle.
  oxidative             NO. oxyR/soxS/soxR have fold-change data but are absent
                        from condition/tf_condition.tsv, so they never enter
                        tf_ids and contribute exactly zero. No ROS species.
  envelope              NO. cpxR likewise has data but is unused;
                        complexation_reactions_removed.tsv:19 removes the CpxR
                        complex for a missing subunit.
  RpoS / general        NO sigma mechanism -- every holoenzyme complexation
                        reaction is deleted (complexation_reactions_removed.tsv
                        :4-11, "Sigma factors are not used in the current
                        version of the model"), so the only RNAP species is
                        core APORNAP-CPLX and rpoD/rpoS/rpoH protein is inert
                        mass. A large sigma-S-dependent gene block IS
                        ppGpp-regulated, so it moves for stringent reasons.
  chaperones/proteases  NO. protein_degradation.py is 95 lines of first-order
                        Poisson decay; ClpP/Lon/HslV/FtsH/DnaK/GroEL are
                        catalytically inert.

So what this script can honestly measure is whether the model's inert stress
proteins DRIFT in the direction a real cell would move them. That is a
statement about growth-rate coupling and ppGpp, not about stress sensing.

The dilution confound, and why the percentile column exists
-----------------------------------------------------------
This is the methodological crux, and ignoring it would make the whole analysis
worthless. Steady-state protein per cell is roughly

	synthesis / (growth rate + degradation)

and burden slows growth substantially (tau ~54 -> ~87 min across the ladder).
So EVERY poorly-degraded protein rises per cell with no regulation whatsoever.
A heat-shock protein at +40% means nothing if the median protein is +40%.

Each gene is therefore reported three ways -- counts per cell, concentration
(per dry mass), and proteome mass fraction -- and then against the only honest
reference: the fold-change distribution of the WHOLE proteome. The percentile
column says where a gene's fold change falls in that distribution. A stress
gene at the 50th percentile has not responded; it has been diluted like
everything else. Effect size is reported first and the percentile alongside; a
z-score would be wrong here because the fold-change distribution is not normal.

The unbiased half
-----------------
The curated panel can only confirm or deny what we already suspected, so all
~4300 monomers are also scanned and ranked, with pathway annotation applied
afterwards. If something is moving that is not on the panel, that is where it
shows up.

Grain and statistics
--------------------
Protein is a slow integrator, so successive generations of one lineage are
strongly autocorrelated: the effective n is the number of SEEDS, not the number
of lineage-generations. Every summary statistic is therefore a per-lineage mean
first, with the CI taken across lineages.

Output
------
	<name>_panel_by_lineage.csv   (variant, seed, pathway, gene) x 3 norms
	<name>_panel_summary.csv      (variant, pathway, gene) mean, CI, FC, pct
	<name>_proteome_scan.csv      every monomer, all three norms, FC, pct
	<name>_global.csv             (variant, seed) growth, ppGpp, machinery
	<name>_manifest.txt           coverage, the modelled-pathway verdict
"""

import csv
import os
import pickle
from collections import Counter, OrderedDict

from matplotlib import pyplot as plt
import numpy as np

from models.ecoli.analysis import variantAnalysisPlot
from wholecell.analysis.analysis_tools import (exportFigure,
	first_cell_with_table)
from wholecell.io.tablereader import TableReader
from wholecell.utils import units

DEFAULT_INDUCTION_GEN = 8

# How many generations at the end of the run to average over. Matches the
# convention in the sibling new-gene analyses, which use the last 8.
LAST_N_GENS = 8

# The curated panel. Grouped by pathway, with the modelled status recorded in
# PATHWAY_STATUS below so the report cannot quietly imply a mechanism exists.
PANEL = OrderedDict([
	('SOS', [
		'recA', 'recN', 'recX', 'uvrA', 'uvrB', 'uvrC', 'uvrD', 'umuC',
		'umuD', 'dinB', 'dinD', 'dinF', 'dinG', 'dinI', 'polB', 'ruvA',
		'ruvB', 'ssb', 'sulA', 'lexA', 'symE', 'ydjM', 'yebG', 'ftsK',
		'recQ', 'phr', 'cho', 'sbmC']),
	('heat_shock', [
		'dnaK', 'dnaJ', 'grpE', 'groL', 'groS', 'htpG', 'ibpA', 'ibpB',
		'clpB', 'clpP', 'clpX', 'lon', 'hslU', 'hslV', 'ftsH', 'rpoH']),
	('oxidative', [
		'katG', 'katE', 'sodA', 'sodB', 'sodC', 'ahpC', 'ahpF', 'dps',
		'trxA', 'grxA', 'oxyR', 'soxS']),
	('envelope', [
		'degP', 'rpoE', 'rseA', 'rseB', 'rseC', 'spy', 'ompC', 'ompF',
		'cpxR']),
	('cold_general', [
		'cspA', 'cspD', 'rpoS', 'osmY', 'osmB', 'osmC', 'otsA', 'otsB',
		'bolA', 'uspA', 'rmf', 'sra', 'mazE', 'mazF']),
	('stringent', ['relA', 'spoT', 'dksA']),
	])

PATHWAY_STATUS = {
	'SOS': 'repressor present (LexA, 0CS); NO induction mechanism',
	'heat_shock': 'NOT modelled; genes are inert mass; 0/16 ppGpp-linked',
	'oxidative': 'NOT modelled; OxyR/SoxS never enter tf_ids',
	'envelope': 'NOT modelled; CpxR complex removed for a missing subunit',
	'cold_general': 'no sigma mechanism; most genes ARE ppGpp-linked',
	'stringent': 'MODELLED -- the only real stress pathway here',
	}


class Plot(variantAnalysisPlot.VariantAnalysisPlot):

	def do_plot(self, inputDir, plotOutDir, plotOutFileName, simDataFile,
			validationDataFile, metadata):

		self._dropped = Counter()

		with open(simDataFile, 'rb') as handle:
			sim_data = pickle.load(handle)

		induction_gen = self._induction_gen(metadata)
		variants = self.ap.get_variants()
		n_gens = self.ap.n_generation
		first_gen = max(induction_gen, n_gens - LAST_N_GENS)
		print('Variants: %s' % (variants,))
		print('Generations %d-%d averaged (induction at %d, %d total).'
			% (first_gen, n_gens - 1, induction_gen, n_gens))

		probe = first_cell_with_table(self.ap.get_cells(), 'MonomerCounts')
		if probe is None:
			print('No cell has a readable MonomerCounts table. Nothing to do.')
			return
		ctx = self._context(sim_data, probe)
		self._report_panel_coverage(ctx)

		# Per (variant, seed): mean over the analysed generations, so that the
		# lineage is the unit of replication.
		per_lineage = []
		global_rows = []
		for variant in variants:
			for seed in self.ap.get_seeds(variant=variant):
				got = self._one_lineage(variant, seed, first_gen, n_gens, ctx)
				if got is None:
					continue
				per_lineage.append(got)
				global_rows.append(got['global'])

		if not per_lineage:
			print('No readable lineages. Nothing written.')
			return

		name = plotOutFileName
		panel_rows, summary_rows, scan_rows = self._summarise(
			per_lineage, ctx, variants)

		self._write(os.path.join(plotOutDir, name + '_panel_by_lineage.csv'),
			panel_rows)
		self._write(os.path.join(plotOutDir, name + '_panel_summary.csv'),
			summary_rows)
		self._write(os.path.join(plotOutDir, name + '_proteome_scan.csv'),
			scan_rows)
		self._write(os.path.join(plotOutDir, name + '_global.csv'),
			global_rows)
		self._manifest(os.path.join(plotOutDir, name + '_manifest.txt'),
			per_lineage, summary_rows, scan_rows, ctx, first_gen, n_gens,
			induction_gen)
		self._figure(plotOutDir, name, summary_rows, variants)

	# ---- setup -----------------------------------------------------------

	@staticmethod
	def _induction_gen(metadata):
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
		reader = TableReader(os.path.join(probe, 'simOut', 'MonomerCounts'))
		monomer_ids = list(reader.readAttribute('monomerIds'))

		transcription = sim_data.process.transcription
		translation = sim_data.process.translation

		monomer_to_cistron = dict(zip(translation.monomer_data['id'],
			translation.monomer_data['cistron_id']))
		cistron_to_gene = dict(zip(transcription.cistron_data['id'],
			transcription.cistron_data['gene_id']))
		gene_to_symbol = dict(zip(
			sim_data.process.replication.gene_data['name'],
			sim_data.process.replication.gene_data['symbol']))

		symbols = []
		for monomer in monomer_ids:
			cistron = monomer_to_cistron.get(monomer, '')
			gene = cistron_to_gene.get(cistron, '')
			symbols.append(gene_to_symbol.get(gene, monomer))

		# symbol -> every monomer index carrying it. A symbol can map to more
		# than one monomer, so counts are summed rather than assumed unique.
		by_symbol = {}
		for index, symbol in enumerate(symbols):
			by_symbol.setdefault(symbol, []).append(index)

		masses = sim_data.getter.get_masses(monomer_ids).asNumber(
			units.fg / units.count)

		ppgpp_genes = set()
		try:
			ppgpp_genes = {
				str(sim_data.common_names.get_common_name(
					gene.split('[')[0]))
				for gene in transcription.ppgpp_regulated_genes}
		except Exception:  # noqa: BLE001 - attribute absent in old vintages
			print('WARNING: could not read ppgpp_regulated_genes; the '
				'ppgpp_linked column will be blank.')

		panel = OrderedDict()
		for pathway, genes in PANEL.items():
			panel[pathway] = [(g, by_symbol.get(g, [])) for g in genes]

		return dict(monomer_ids=monomer_ids, symbols=symbols,
			by_symbol=by_symbol, masses=np.asarray(masses, dtype=float),
			panel=panel, ppgpp_genes=ppgpp_genes)

	@staticmethod
	def _report_panel_coverage(ctx):
		print('')
		print('CURATED PANEL COVERAGE')
		for pathway, entries in ctx['panel'].items():
			present = [g for g, idx in entries if idx]
			missing = [g for g, idx in entries if not idx]
			linked = [g for g, _ in entries if g in ctx['ppgpp_genes']]
			print('  %-14s %2d/%2d in model, %2d ppGpp-linked -- %s'
				% (pathway, len(present), len(entries), len(linked),
					PATHWAY_STATUS[pathway]))
			if missing:
				print('      absent: %s' % ', '.join(missing))
		print('')

	# ---- per lineage -----------------------------------------------------

	def _one_lineage(self, variant, seed, first_gen, n_gens, ctx):
		"""
		Mean monomer counts over the analysed generations of one lineage.

		The lineage, not the generation, is the unit of replication here:
		protein is a slow integrator, so consecutive generations are heavily
		autocorrelated and treating them as independent would understate every
		interval.
		"""
		cells = self.ap.get_cells(variant=[variant], seed=[seed],
			generation=np.arange(first_gen, n_gens), only_successful=True)
		if not len(cells):
			self._dropped['no_cells'] += 1
			return None

		counts = []
		dry = []
		protein = []
		ppgpp = []
		rnap = []
		ribosome = []
		taus = []
		for path in cells:
			sim_out = os.path.join(path, 'simOut')
			try:
				mono = TableReader(os.path.join(sim_out, 'MonomerCounts')
					).readColumn('monomerCounts')
				mass = TableReader(os.path.join(sim_out, 'Mass'))
				time = TableReader(os.path.join(sim_out, 'Main')
					).readColumn('time')
			except Exception as exc:  # noqa: BLE001 - preempted simOut
				self._dropped['unreadable'] += 1
				print('skip v%d s%d %s: %s' % (variant, seed,
					os.path.basename(path), exc))
				continue
			counts.append(np.asarray(mono, dtype=float).mean(axis=0))
			dry.append(float(np.asarray(
				mass.readColumn('dryMass'), dtype=float).mean()))
			protein.append(float(np.asarray(
				mass.readColumn('proteinMass'), dtype=float).mean()))
			taus.append(float((time[-1] - time[0]) / 60.))
			ppgpp.append(self._mean_col(sim_out, 'GrowthLimits',
				'ppgpp_conc'))
			rnap.append(self._unique_mean(sim_out, 'active_RNAP'))
			ribosome.append(self._unique_mean(sim_out, 'active_ribosome'))

		if not counts:
			self._dropped['no_readable_cells'] += 1
			return None

		mean_counts = np.mean(np.vstack(counts), axis=0)
		dry_mass = float(np.mean(dry))
		protein_mass = float(np.mean(protein))

		return dict(
			variant=variant, seed=seed,
			n_gens=len(counts),
			counts=mean_counts,
			# Three normalisations, because the first one alone cannot separate
			# a real response from slower dilution.
			per_fg=mean_counts / dry_mass if dry_mass > 0
				else np.full_like(mean_counts, np.nan),
			mass_frac=(mean_counts * ctx['masses'] / protein_mass
				if protein_mass > 0 else np.full_like(mean_counts, np.nan)),
			**{'global': dict(
				variant=variant, seed=seed, n_gens=len(counts),
				doubling_time_min=float(np.mean(taus)),
				dry_mass_fg=dry_mass, protein_mass_fg=protein_mass,
				ppgpp_conc=float(np.nanmean(ppgpp)),
				active_rnap=float(np.nanmean(rnap)),
				active_ribosome=float(np.nanmean(ribosome)))})

	@staticmethod
	def _mean_col(sim_out, table, column):
		try:
			return float(np.nanmean(TableReader(
				os.path.join(sim_out, table)).readColumn(column)))
		except Exception:  # noqa: BLE001
			return float('nan')

	@staticmethod
	def _unique_mean(sim_out, key):
		try:
			reader = TableReader(os.path.join(sim_out,
				'UniqueMoleculeCounts'))
			ids = list(reader.readAttribute('uniqueMoleculeIds'))
			if key not in ids:
				return float('nan')
			return float(np.asarray(reader.readColumn('uniqueMoleculeCounts'),
				dtype=float)[:, ids.index(key)].mean())
		except Exception:  # noqa: BLE001
			return float('nan')

	# ---- summarise -------------------------------------------------------

	def _summarise(self, per_lineage, ctx, variants):
		"""
		Per-lineage rows, per-gene summaries, and the whole-proteome scan.

		The control is the lowest variant index present, which is the knockout
		in every new-gene variant layout used here.
		"""
		norms = ('counts', 'per_fg', 'mass_frac')
		control = min(variants)

		by_variant = {}
		for row in per_lineage:
			by_variant.setdefault(row['variant'], []).append(row)

		# Proteome-wide mean per variant and normalisation.
		means = {}
		for variant, rows in by_variant.items():
			for norm in norms:
				means[(variant, norm)] = np.mean(
					np.vstack([r[norm] for r in rows]), axis=0)

		# Whole-proteome fold-change distributions, which are the reference the
		# percentile column is measured against. Monomers absent in the control
		# are excluded rather than given an infinite fold change.
		fc = {}
		for variant in by_variant:
			for norm in norms:
				base = means.get((control, norm))
				here = means[(variant, norm)]
				if base is None:
					fc[(variant, norm)] = np.full_like(here, np.nan)
					continue
				with np.errstate(divide='ignore', invalid='ignore'):
					ratio = np.where(base > 0, here / base, np.nan)
				fc[(variant, norm)] = ratio

		panel_rows = []
		for row in per_lineage:
			for pathway, entries in ctx['panel'].items():
				for gene, idx in entries:
					if not idx:
						continue
					panel_rows.append(dict(
						variant=row['variant'], seed=row['seed'],
						pathway=pathway, gene=gene,
						n_monomers=len(idx),
						counts=float(row['counts'][idx].sum()),
						per_fg=float(row['per_fg'][idx].sum()),
						mass_frac=float(row['mass_frac'][idx].sum())))

		summary_rows = []
		for variant in sorted(by_variant):
			rows = by_variant[variant]
			for pathway, entries in ctx['panel'].items():
				for gene, idx in entries:
					if not idx:
						continue
					entry = dict(
						variant=variant, pathway=pathway, gene=gene,
						pathway_status=PATHWAY_STATUS[pathway],
						ppgpp_linked=int(gene in ctx['ppgpp_genes']),
						n_lineages=len(rows))
					for norm in norms:
						vals = np.array([r[norm][idx].sum() for r in rows],
							dtype=float)
						mean, low, high = self._ci(vals)
						entry['%s_mean' % norm] = mean
						entry['%s_ci_low' % norm] = low
						entry['%s_ci_high' % norm] = high
						base_vals = np.array(
							[r[norm][idx].sum()
								for r in by_variant.get(control, [])],
							dtype=float)
						base = np.nanmean(base_vals) if len(base_vals) \
							else float('nan')
						ratio = (mean / base if base and np.isfinite(base)
							and base > 0 else float('nan'))
						entry['%s_fc' % norm] = ratio
						# Where this gene's fold change sits in the
						# whole-proteome fold-change distribution. 50 means it
						# moved exactly like the median protein, i.e. it did
						# not respond.
						entry['%s_proteome_pct' % norm] = self._percentile(
							fc[(variant, norm)], ratio)
					summary_rows.append(entry)

		scan_rows = []
		for variant in sorted(by_variant):
			rows = by_variant[variant]
			for index, monomer in enumerate(ctx['monomer_ids']):
				entry = dict(
					variant=variant, monomer=monomer,
					gene=ctx['symbols'][index],
					ppgpp_linked=int(
						ctx['symbols'][index] in ctx['ppgpp_genes']),
					pathway=self._pathway_of(ctx['symbols'][index]),
					n_lineages=len(rows))
				for norm in norms:
					vals = np.array([r[norm][index] for r in rows],
						dtype=float)
					mean, low, high = self._ci(vals)
					entry['%s_mean' % norm] = mean
					entry['%s_ci_low' % norm] = low
					entry['%s_ci_high' % norm] = high
					entry['%s_fc' % norm] = float(
						fc[(variant, norm)][index])
					entry['%s_proteome_pct' % norm] = self._percentile(
						fc[(variant, norm)], fc[(variant, norm)][index])
				scan_rows.append(entry)

		return panel_rows, summary_rows, scan_rows

	@staticmethod
	def _pathway_of(gene):
		for pathway, genes in PANEL.items():
			if gene in genes:
				return pathway
		return ''

	@staticmethod
	def _ci(values):
		"""Mean and a 95% normal-approximation CI across lineages."""
		vals = np.asarray(values, dtype=float)
		vals = vals[np.isfinite(vals)]
		if not len(vals):
			return float('nan'), float('nan'), float('nan')
		mean = float(vals.mean())
		if len(vals) < 2:
			return mean, float('nan'), float('nan')
		half = 1.96 * float(vals.std(ddof=1)) / np.sqrt(len(vals))
		return mean, mean - half, mean + half

	@staticmethod
	def _percentile(distribution, value):
		"""Percentile of `value` within the proteome fold-change spread."""
		if value is None or not np.isfinite(value):
			return float('nan')
		dist = np.asarray(distribution, dtype=float)
		dist = dist[np.isfinite(dist)]
		if not len(dist):
			return float('nan')
		return float((dist < value).mean() * 100.)

	# ---- manifest and figure --------------------------------------------

	def _manifest(self, path, per_lineage, summary_rows, scan_rows, ctx,
			first_gen, n_gens, induction_gen):
		lines = []
		add = lines.append
		add('STRESS PROTEOME SCAN UNDER NEW-GENE BURDEN')
		add('')
		add('WHAT IS MODELLED')
		for pathway, status in PATHWAY_STATUS.items():
			add('  %-14s %s' % (pathway, status))
		add('')
		add('Only the stringent response is mechanistically modelled. Movement')
		add('in any other panel is passive growth-rate coupling, and the')
		add('proteome percentile columns are what separate a real change from')
		add('dilution. A gene near the 50th percentile moved like the median')
		add('protein, i.e. it did not respond.')
		add('')
		add('GENERATIONS AVERAGED: %d-%d  (induction at %d, %d total)'
			% (first_gen, n_gens - 1, induction_gen, n_gens))
		add('')
		add('LINEAGE COVERAGE')
		by_variant = {}
		for row in per_lineage:
			by_variant.setdefault(row['variant'], []).append(row)
		for variant in sorted(by_variant):
			rows = by_variant[variant]
			add('  variant %-3d  n_lineages=%-4d gens_per_lineage %d-%d'
				% (variant, len(rows), min(r['n_gens'] for r in rows),
					max(r['n_gens'] for r in rows)))
		add('')
		add('DROPPED')
		if self._dropped:
			for key, count in sorted(self._dropped.items()):
				add('  %-20s %d' % (key, count))
		else:
			add('  none')
		add('')
		add('PANEL COVERAGE')
		for pathway, entries in ctx['panel'].items():
			present = [g for g, idx in entries if idx]
			missing = [g for g, idx in entries if not idx]
			linked = [g for g, _ in entries if g in ctx['ppgpp_genes']]
			add('  %-14s %2d/%2d in model, %2d ppGpp-linked'
				% (pathway, len(present), len(entries), len(linked)))
			if missing:
				add('      absent from model: %s' % ', '.join(missing))
		add('')
		add('ROW COUNTS')
		add('  panel_summary   %d' % len(summary_rows))
		add('  proteome_scan   %d' % len(scan_rows))
		add('')
		add('PROTEOME FOLD-CHANGE REFERENCE (concentration, per dry mass)')
		add('  %-8s %9s %9s %9s' % ('variant', 'median', 'p10', 'p90'))
		for variant in sorted(by_variant):
			vals = np.array([r['per_fg_fc'] for r in scan_rows
				if r['variant'] == variant], dtype=float)
			vals = vals[np.isfinite(vals)]
			if len(vals):
				add('  %-8d %9.3f %9.3f %9.3f'
					% (variant, np.median(vals),
						np.percentile(vals, 10), np.percentile(vals, 90)))
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
			writer = csv.DictWriter(handle, fieldnames=keys,
				extrasaction='ignore')
			writer.writeheader()
			writer.writerows(rows)
		print('wrote %s (%d rows)' % (path, len(rows)))

	@staticmethod
	def _figure(out_dir, name, summary_rows, variants):
		"""
		Diagnostic: each pathway's concentration fold change against the
		whole-proteome median, per variant. If a pathway's band straddles the
		proteome median it has not responded.
		"""
		pathways = list(PANEL.keys())
		shown = [v for v in sorted(variants) if v != min(variants)]
		if not shown:
			shown = sorted(variants)

		fig, axes = plt.subplots(len(pathways), 1,
			figsize=(9, 1.5 * len(pathways)), sharex=True)
		for ax, pathway in zip(np.atleast_1d(axes), pathways):
			for offset, variant in enumerate(shown):
				vals = np.array([r['per_fg_fc'] for r in summary_rows
					if r['variant'] == variant and r['pathway'] == pathway],
					dtype=float)
				vals = vals[np.isfinite(vals)]
				if not len(vals):
					continue
				ax.scatter(np.full(len(vals), offset)
					+ np.random.uniform(-0.12, 0.12, len(vals)),
					vals, s=9, alpha=0.7)
			ax.axhline(1.0, color='0.7', lw=0.8, ls='--')
			ax.set_ylabel('%s\nFC' % pathway, fontsize='xx-small')
			ax.set_yscale('log')
			ax.tick_params(labelsize='xx-small')
		np.atleast_1d(axes)[-1].set_xticks(range(len(shown)))
		np.atleast_1d(axes)[-1].set_xticklabels(
			['v%d' % v for v in shown], fontsize='x-small')
		np.atleast_1d(axes)[0].set_title(
			'Concentration fold change vs control, by pathway. '
			'Dashed line: no change.', fontsize='small')
		plt.tight_layout()
		exportFigure(plt, out_dir, name, metadata=None)
		plt.close('all')


if __name__ == '__main__':
	Plot().cli()
