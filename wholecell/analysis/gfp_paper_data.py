"""
Shared data access for the GFP burden paper's figure panels.

Every data panel needs the same handful of per-cell-cycle numbers: growth,
gfp mRNA and protein, the machinery on gfp, the pools it draws from. Reading
them separately in each panel would re-read thousands of simulations per panel
and, worse, let panels quietly disagree on how a quantity is counted. This
module reads each cell cycle once, applies one set of counting rules, and
caches the result beside the batch, so every panel starts from the same table.

Three levels, each built from the one below:

	cycle_table     one row per (variant, seed, generation): raw time means
	lineage_table   one row per (variant, seed): means over the analysis window
	design_table    one row per variant: means over surviving lineages

Counting rules, fixed here so no panel has to decide them
-----------------------------------------------------------
- **Physical mRNA.** Some batches split the mRNA pool across several
	RNACounts fields (copies being degraded are reported separately from
	intact ones). The count used everywhere is the sum of every field the
	batch provides, so two batches never differ for bookkeeping reasons.
- **RNAP on gfp** is the number of gfp transcripts still being synthesised,
	one polymerase each: `rnap_engaged_mRNA_counts` where the batch has it,
	otherwise `partial_mRNA_counts`.
- **Ribosomes on gfp** is the sum of `n_ribosomes_on_each_mRNA` over the
	mRNAs whose transcription unit is gfp. RibosomeData indexes transcription
	units over all RNAs while RNACounts indexes mRNAs only, so each listener's
	index is resolved from its own id list, never shared.
- **Ratios are ratios of time means within a cycle** (for example ribosomes
	on gfp / gfp mRNA, both averaged over the cycle). That is the per-mRNA
	average weighted by time, and it does not blow up at time steps where the
	denominator happens to be zero.
- **Total ribosomes** are active ribosomes plus free 70S-forming pairs,
	min(30S, 50S), as in the existing heatmap analyses.
- **Growth** is the cycle mean of `Mass.instantaneous_growth_rate`, reported
	per hour (decision D4: growth rate, not doubling time). Doubling time is
	kept as well because some panels and checks still want it.

What counts as a lineage surviving
----------------------------------
A lineage survives if its last analysis generation (ANALYSIS_LAST_GEN) both
exists and divided. Design means use surviving lineages only, so a design
that lost lineages is biased toward its healthier ones; those designs are
flagged `failed` and drawn hollow (gfp_paper_style.design_scatter).

Only cell cycles that divided enter the window means; a cycle that stopped
part-way has a truncated time series and no meaningful growth rate.

Paper panels always read the whole batch. They deliberately ignore the
analysis runner's variant, seed and generation filters, and in particular
--only-successful, which would drop exactly the dying lineages that the
survival counts and hollow markers exist to report. A panel that wants a
subset of designs selects it from the tables.
"""

from __future__ import annotations

import glob
import json
import os
import pickle
import re
import warnings
from dataclasses import asdict, dataclass
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from wholecell.analysis import gfp_paper_constants as C
from wholecell.io.tablereader import TableReader
from wholecell.utils import parallelization

# Bump whenever read_cycle changes what a column means, so stale caches are
# rebuilt rather than silently reused.
CACHE_VERSION = 1
CACHE_SUBDIR = os.path.join('paper_figures', '_cache')
REBUILD_ENV = 'GFP_PAPER_REBUILD_CACHE'

AVOGADRO = 6.02214076e23
SECONDS_PER_HOUR = 3600.0
CONTROL_VARIANT = 0


# --------------------------------------------------------------------------
# What the new gene is, read from sim_data rather than hard-coded
@dataclass(frozen=True)
class GeneInfo:
	"""The inserted gene's ids, sizes and the footprints that cap it."""
	rna_id: str                     # transcription unit id, e.g. NG001_RNA[c]
	monomer_id: str
	monomer_mass_fg: float          # one GFP molecule
	length_nt: float                # transcription unit length
	rnap_footprint_nt: float
	ribosome_footprint_nt: float
	s30_id: str
	s50_id: str

	@property
	def ribosomes_when_full(self) -> float:
		"""Ribosomes on one gfp mRNA packed footprint to footprint."""
		return self.length_nt / self.ribosome_footprint_nt

	@property
	def rnaps_when_full(self) -> float:
		"""RNAPs on one gfp copy packed footprint to footprint."""
		return self.length_nt / self.rnap_footprint_nt


def gene_info(sim_data) -> GeneInfo:
	"""Read the new gene's ids and sizes from sim_data."""
	from wholecell.utils import units
	transcription = sim_data.process.transcription
	translation = sim_data.process.translation

	new_cistrons = transcription.cistron_data['id'][
		transcription.cistron_data['is_new_gene']]
	if len(new_cistrons) != 1:
		raise ValueError('expected exactly one new gene cistron, found %d;'
			' the paper panels assume a single inserted gene'
			% len(new_cistrons))
	cistron = new_cistrons[0]

	rna_ids = list(transcription.rna_data['id'])
	tu_rows = list(transcription.cistron_id_to_rna_indexes(cistron))
	if len(tu_rows) != 1:
		# The panels count gfp mRNA as one transcription unit; an operon
		# containing the new gene would need its own counting rules.
		raise ValueError('expected one transcription unit for %s, found %d'
			% (cistron, len(tu_rows)))
	tu = tu_rows[0]

	monomers = translation.monomer_data
	monomer_rows = np.where(monomers['cistron_id'] == cistron)[0]
	if len(monomer_rows) != 1:
		raise ValueError('expected one monomer for %s, found %d'
			% (cistron, len(monomer_rows)))
	m = monomer_rows[0]
	mw = monomers['mw'][m].asNumber(units.g / units.mol)

	return GeneInfo(
		rna_id=rna_ids[tu],
		monomer_id=monomers['id'][m],
		monomer_mass_fg=float(mw / AVOGADRO * 1e15),
		length_nt=float(transcription.rna_data['length'][tu].asNumber(
			units.nt)),
		rnap_footprint_nt=float(
			transcription.active_rnap_footprint_size.asNumber(units.nt)),
		ribosome_footprint_nt=float(
			translation.active_ribosome_footprint_size.asNumber(units.nt)),
		s30_id=sim_data.molecule_ids.s30_full_complex,
		s50_id=sim_data.molecule_ids.s50_full_complex,
	)


# --------------------------------------------------------------------------
# Which design each variant index is
_DESIGN_LIST_DESC = re.compile(
	r'expression factor ([-\d.eE]+), translation\s+efficiency ([-\d.eE]+)')


def design_of(batch_dir: str, metadata: dict, variant: int
		) -> Tuple[float, float]:
	"""
	(promoter, RBS) for a variant index, as the paper reports them.

	Promoter is the expression factor before exponentiation (8 means
	10 ** 7 times native), RBS is the translation efficiency. The control is
	(0, 0).

	The decoding depends on which variant produced the batch:
	- new_gene_internal_shift and its relatives: a strict factorial, decoded
		from the lists the batch recorded in its metadata, never from the
		module's current lists, which may have changed since the run.
	- new_gene_design_list: decoded from the variant's own recorded
		description, cross-checked against NEW_GENE_DESIGNS.
	"""
	design_index = variant % 1000
	if design_index == CONTROL_VARIANT:
		return 0.0, 0.0
	name = metadata.get('variant', '')

	if name == 'new_gene_design_list':
		from models.ecoli.sim.variants.new_gene_design_list import (
			NEW_GENE_DESIGNS)
		recorded = _recorded_design(batch_dir, variant)
		if design_index >= len(NEW_GENE_DESIGNS):
			if recorded is None:
				raise ValueError('variant %d is outside NEW_GENE_DESIGNS and'
					' has no recorded description' % variant)
			return recorded
		expression, trl_eff = NEW_GENE_DESIGNS[design_index]
		current = (float(expression), float(trl_eff))
		if recorded is not None and not np.allclose(recorded, current):
			raise ValueError(
				'variant %d was run as %s but NEW_GENE_DESIGNS now says %s;'
				' the list was edited after this batch ran'
				% (variant, recorded, current))
		return current

	factors = metadata.get('new_gene_expression_factors')
	efficiencies = metadata.get('new_gene_translation_efficiency_values')
	if factors is None or efficiencies is None:
		raise ValueError('batch metadata has no new gene expression factors'
			' or translation efficiencies, and variant %r is not a design'
			' list; cannot decode variant %d' % (name, variant))
	separator = len(efficiencies)
	trl_eff_index = design_index % separator
	expression_index = design_index // separator + (trl_eff_index != 0)
	return float(factors[expression_index]), float(
		efficiencies[trl_eff_index])


def _recorded_design(batch_dir, variant) -> Optional[Tuple[float, float]]:
	"""The (promoter, RBS) a design-list variant wrote in its description."""
	for path in glob.glob(os.path.join(batch_dir, '*_%06d' % variant,
			'metadata', 'description')):
		with open(path) as f:
			match = _DESIGN_LIST_DESC.search(f.read())
		if match:
			return float(match.group(1)), float(match.group(2))
	return None


# --------------------------------------------------------------------------
# Locating cell cycles
def sim_dirs(batch_dir: str) -> List[Tuple[int, int, int, str]]:
	"""
	[(variant, seed, generation, simOut path)] for every cell in a batch.

	Lineage-only batches (one daughter per generation) are assumed, which is
	how every paper batch is run.
	"""
	out = []
	pattern = os.path.join(batch_dir, '*_[0-9]*', '[0-9]' * 6,
		'generation_[0-9]*', '[0-9]' * 6, 'simOut')
	for path in sorted(glob.glob(pattern)):
		parts = path.split(os.sep)
		variant = re.fullmatch(r'.*_(\d{6})', parts[-5])
		generation = re.fullmatch(r'generation_(\d{6})', parts[-3])
		# Strict matches, so copies such as "generation_000000 11.41.58 AM"
		# left behind by a file browser are never read as a generation.
		if variant is None or generation is None:
			continue
		out.append((int(variant.group(1)), int(parts[-4]),
			int(generation.group(1)), path))
	return out


# --------------------------------------------------------------------------
# Reading one cell cycle
class _Cell:
	"""Lazy, tolerant column access for one simOut directory."""

	def __init__(self, path):
		self.path = path
		self._readers = {}
		self.missing = []

	def reader(self, listener):
		if listener not in self._readers:
			try:
				self._readers[listener] = TableReader(
					os.path.join(self.path, listener))
			except Exception:                               # noqa: BLE001
				self._readers[listener] = None
		return self._readers[listener]

	def column(self, listener, name):
		"""The column, or None (recorded in `missing`) if it is absent."""
		reader = self.reader(listener)
		if reader is None or name not in reader.columnNames():
			self.missing.append('%s/%s' % (listener, name))
			return None
		return reader.readColumn(name)

	def attribute(self, listener, name):
		reader = self.reader(listener)
		if reader is None or name not in reader.attributeNames():
			return None
		return reader.readAttribute(name)


def _mean(values) -> float:
	if values is None or np.size(values) == 0:
		return np.nan
	return float(np.nanmean(values))


def _index(ids, wanted) -> Optional[int]:
	if ids is None:
		return None
	ids = list(ids)
	return ids.index(wanted) if wanted in ids else None


def read_cycle(path: str, gene: GeneInfo) -> Dict[str, object]:
	"""
	Time means of the raw quantities for one cell cycle.

	Only raw means are stored; every ratio is derived later in
	add_derived_columns, so a change of definition does not need a re-read.
	"""
	cell = _Cell(path)
	row: Dict[str, object] = {}

	time = cell.column('Main', 'time')
	row['n_steps'] = 0 if time is None else int(time.size)
	row['doubling_time_min'] = (float(time[-1] - time[0]) / 60.0
		if time is not None and time.size > 1 else np.nan)
	row['divided'] = os.path.exists(
		os.path.join(path, 'Daughter1_inherited_state.cPickle'))

	growth = cell.column('Mass', 'instantaneous_growth_rate')
	# The first step is recorded before any growth has happened.
	row['growth_rate_per_h'] = (_mean(growth[1:]) * SECONDS_PER_HOUR
		if growth is not None and growth.size > 1 else np.nan)
	for name in ('cellMass', 'dryMass', 'proteinMass', 'mRnaMass',
			'rRnaMass'):
		row[name] = _mean(cell.column('Mass', name))

	# GFP protein
	monomers = cell.column('MonomerCounts', 'monomerCounts')
	m = _index(cell.attribute('MonomerCounts', 'monomerIds'), gene.monomer_id)
	if monomers is not None and m is not None:
		gfp_monomers = _mean(monomers[:, m].astype(float))
		row['gfp_monomers'] = gfp_monomers
		row['gfp_protein_mass_fg'] = gfp_monomers * gene.monomer_mass_fg
	else:
		row['gfp_monomers'] = row['gfp_protein_mass_fg'] = np.nan

	# gfp mRNA: the physical count, summing every field the batch provides
	mrna_col = _index(cell.attribute('RNACounts', 'mRNA_ids'), gene.rna_id)
	physical = None
	mrna_fields = []
	for name in ('mRNA_counts', 'decaying_mRNA_counts',
			'decaying_cleared_mRNA_counts'):
		reader = cell.reader('RNACounts')
		if reader is None or name not in reader.columnNames():
			continue
		values = reader.readColumn(name)
		physical = values.astype(float) if physical is None else (
			physical + values)
		mrna_fields.append(name)
	row['mrna_fields'] = ';'.join(mrna_fields)
	if physical is not None and mrna_col is not None:
		row['gfp_mrna'] = _mean(physical[:, mrna_col])
		row['all_mrna'] = _mean(physical.sum(axis=1))
	else:
		row['gfp_mrna'] = row['all_mrna'] = np.nan

	# RNAP on gfp: transcripts still being made, one polymerase each
	engaged = None
	reader = cell.reader('RNACounts')
	for name in ('rnap_engaged_mRNA_counts', 'partial_mRNA_counts'):
		if reader is not None and name in reader.columnNames():
			engaged = reader.readColumn(name)
			row['rnap_field'] = name
			break
	row['rnap_on_gfp'] = (_mean(engaged[:, mrna_col])
		if engaged is not None and mrna_col is not None else np.nan)

	# Promoter copies and crowding flags, in RnaSynthProb's own index
	tu_col = _index(cell.attribute('RnaSynthProb', 'rnaIds'), gene.rna_id)
	copies = cell.column('RnaSynthProb', 'promoter_copy_number')
	row['gfp_promoter_copies'] = (_mean(copies[:, tu_col])
		if copies is not None and tu_col is not None else np.nan)
	crowded = cell.column('RnaSynthProb', 'tu_is_overcrowded')
	row['gfp_rnap_overcrowded_frac'] = (_mean(crowded[:, tu_col])
		if crowded is not None and tu_col is not None else np.nan)
	actual = cell.column('RnaSynthProb', 'actual_rna_synth_prob')
	target = cell.column('RnaSynthProb', 'target_rna_synth_prob')
	row['gfp_actual_synth_prob'] = (_mean(actual[:, tu_col])
		if actual is not None and tu_col is not None else np.nan)
	row['gfp_target_synth_prob'] = (_mean(target[:, tu_col])
		if target is not None and tu_col is not None else np.nan)

	# Ribosomes on gfp mRNAs. mRNA_TU_index runs over every transcription
	# unit, which is RnaSynthProb's rnaIds order.
	tu_index = cell.column('RibosomeData', 'mRNA_TU_index')
	per_mrna = cell.column('RibosomeData', 'n_ribosomes_on_each_mRNA')
	if tu_index is not None and per_mrna is not None and tu_col is not None:
		loaded = np.nan_to_num(per_mrna)
		row['ribosomes_on_gfp'] = _mean(
			np.where(tu_index == tu_col, loaded, 0).sum(axis=1))
		row['ribosomes_on_mrna'] = _mean(loaded.sum(axis=1))
	else:
		row['ribosomes_on_gfp'] = row['ribosomes_on_mrna'] = np.nan
	ribo_crowded = cell.column('RibosomeData', 'mRNA_is_overcrowded')
	ribo_monomer = _index(cell.attribute('RibosomeData', 'monomerIds'),
		gene.monomer_id)
	row['gfp_ribosome_overcrowded_frac'] = (
		_mean(ribo_crowded[:, ribo_monomer])
		if ribo_crowded is not None and ribo_monomer is not None else np.nan)
	row['elongation_rate'] = _mean(
		cell.column('RibosomeData', 'effectiveElongationRate'))

	# Machinery pools
	counts = cell.column('UniqueMoleculeCounts', 'uniqueMoleculeCounts')
	names = cell.attribute('UniqueMoleculeCounts', 'uniqueMoleculeIds')
	for key, molecule in (('active_ribosomes', 'active_ribosome'),
			('active_rnap', 'active_RNAP')):
		i = _index(names, molecule)
		row[key] = (_mean(counts[:, i])
			if counts is not None and i is not None else np.nan)
	bulk = cell.column('BulkMolecules', 'counts')
	bulk_names = cell.attribute('BulkMolecules', 'objectNames')
	i30, i50 = _index(bulk_names, gene.s30_id), _index(bulk_names, gene.s50_id)
	row['free_ribosomes'] = (
		_mean(np.minimum(bulk[:, i30], bulk[:, i50]))
		if bulk is not None and i30 is not None and i50 is not None
		else np.nan)

	row['ppgpp_conc'] = _mean(cell.column('GrowthLimits', 'ppgpp_conc'))
	row['missing_columns'] = ';'.join(sorted(set(cell.missing)))
	return row


def add_derived_columns(df: pd.DataFrame, gene: GeneInfo) -> pd.DataFrame:
	"""
	Ratios the panels plot, each defined once.

	Works on any of the three tables, since each carries the same raw means.
	"""
	df = df.copy()
	# A raw column the batch could not provide becomes NaN, so the ratios
	# built on it are NaN rather than the whole table failing.
	for name in RAW_COLUMNS:
		if name not in df:
			df[name] = np.nan
	with np.errstate(divide='ignore', invalid='ignore'):
		df['total_ribosomes'] = df['active_ribosomes'] + df['free_ribosomes']
		df['gfp_proteome_share'] = df['gfp_protein_mass_fg'] / df[
			'proteinMass']
		df['gfp_mrna_share'] = df['gfp_mrna'] / df['all_mrna']
		df['ribosomes_per_gfp_mrna'] = _safe_ratio(
			df['ribosomes_on_gfp'], df['gfp_mrna'])
		df['ribosome_spacing_nt'] = _safe_ratio(
			gene.length_nt, df['ribosomes_per_gfp_mrna'])
		df['ribosome_ceiling_per_gfp_mrna'] = _safe_ratio(
			df['total_ribosomes'], df['gfp_mrna'])
		df['rnap_per_gfp_copy'] = _safe_ratio(
			df['rnap_on_gfp'], df['gfp_promoter_copies'])
		df['rnap_spacing_nt'] = _safe_ratio(
			gene.length_nt, df['rnap_per_gfp_copy'])
		df['gfp_ribosome_share'] = _safe_ratio(
			df['ribosomes_on_gfp'], df['active_ribosomes'])
		df['gfp_rnap_share'] = _safe_ratio(
			df['rnap_on_gfp'], df['active_rnap'])
		# The RNAP share gfp would hold if every copy were packed footprint to
		# footprint: the ceiling the spacing limit imposes.
		df['gfp_rnap_share_at_spacing_limit'] = _safe_ratio(
			gene.rnaps_when_full * df['gfp_promoter_copies'],
			df['active_rnap'])
	return df


def _safe_ratio(numerator, denominator):
	"""numerator / denominator, NaN wherever the denominator is not positive."""
	numerator = np.asarray(numerator, dtype=float)
	denominator = np.asarray(denominator, dtype=float)
	out = np.full(np.broadcast(numerator, denominator).shape, np.nan)
	ok = np.broadcast_to(denominator > 0, out.shape)
	np.divide(numerator, denominator, out=out, where=ok)
	return out


# --------------------------------------------------------------------------
# The three tables
def _read_one(args):
	path, gene_dict = args
	try:
		return read_cycle(path, GeneInfo(**gene_dict))
	except Exception as e:                                  # noqa: BLE001
		return {'read_error': '%s: %s' % (type(e).__name__, e)}


def cycle_table(batch_dir: str, sim_data, metadata: dict, plot_out_dir: str,
		cpus: int = 1) -> pd.DataFrame:
	"""
	One row per cell cycle in the batch, every generation, with derived
	columns and each variant's (promoter, RBS).

	Cached in <plot_out_dir>/paper_figures/_cache/cycles.csv. Cycles already
	in the cache are not re-read, so adding generations or seeds to a batch
	only reads the new ones. Set GFP_PAPER_REBUILD_CACHE=1 to re-read all.
	"""
	gene = gene_info(sim_data)
	cache_dir = os.path.join(plot_out_dir, CACHE_SUBDIR)
	os.makedirs(cache_dir, exist_ok=True)
	cache_path = os.path.join(cache_dir, 'cycles.csv')
	stamp_path = os.path.join(cache_dir, 'cycles.json')
	stamp = {'cache_version': CACHE_VERSION, 'gene': asdict(gene)}

	cached = pd.DataFrame()
	if os.path.exists(cache_path) and not os.environ.get(REBUILD_ENV):
		try:
			with open(stamp_path) as f:
				if json.load(f) == stamp:
					cached = pd.read_csv(cache_path)
		except Exception:                                   # noqa: BLE001
			cached = pd.DataFrame()

	cells = sim_dirs(batch_dir)
	done = set(cached['path']) if len(cached) else set()
	todo = [c for c in cells if c[3] not in done]
	if todo:
		pool_size = max(1, min(cpus, len(todo)))
		jobs = [(path, asdict(gene)) for _, _, _, path in todo]
		if pool_size > 1:
			pool = parallelization.pool(pool_size)
			rows = pool.map(_read_one, jobs)
			pool.close()
			pool.join()
		else:
			rows = [_read_one(job) for job in jobs]
		fresh = pd.DataFrame([dict(variant=v, seed=s, generation=g,
			path=p, **r) for (v, s, g, p), r in zip(todo, rows)])
		cached = pd.concat([cached, fresh], ignore_index=True)
		# Several paper analyses may build the cache at once. Writing to a
		# private file and renaming it means a reader sees either the old
		# cache or the new one, never half of one; the worst case is
		# duplicated reading, not a corrupt table.
		suffix = '.%d.tmp' % os.getpid()
		cached.to_csv(cache_path + suffix, index=False)
		with open(stamp_path + suffix, 'w') as f:
			json.dump(stamp, f, indent=1, sort_keys=True)
		os.replace(stamp_path + suffix, stamp_path)
		os.replace(cache_path + suffix, cache_path)

	present = {c[3] for c in cells}
	df = cached[cached['path'].isin(present)].copy()
	if 'read_error' in df and df['read_error'].notna().any():
		bad = df[df['read_error'].notna()]
		warnings.warn('%d cell cycles could not be read and are excluded,'
			' e.g. %s: %s' % (len(bad), bad['path'].iloc[0],
			bad['read_error'].iloc[0]))
		df = df[df['read_error'].isna()]
	# After a lineage dies the workflow still creates the later generations'
	# folders, but they hold no data. Those cells never ran: dropping them
	# keeps last_generation honest and keeps them out of the missing-column
	# warning, which is meant for batches that lack a listener.
	df = df[df['n_steps'] > 0]

	designs ={v: design_of(batch_dir, metadata, v)
		for v in sorted(df['variant'].unique())}
	df['promoter'] = df['variant'].map(lambda v: designs[v][0])
	df['rbs'] = df['variant'].map(lambda v: designs[v][1])
	df['divided'] = df['divided'].astype(bool)
	return add_derived_columns(df, gene).sort_values(
		['variant', 'seed', 'generation']).reset_index(drop=True)


# Raw means averaged over the window; derived columns are recomputed from
# these, never averaged, so a lineage's ribosomes per mRNA is the ratio of its
# mean ribosomes to its mean mRNA.
RAW_COLUMNS = (
	'doubling_time_min', 'growth_rate_per_h', 'cellMass', 'dryMass',
	'proteinMass', 'mRnaMass', 'rRnaMass', 'gfp_monomers',
	'gfp_protein_mass_fg', 'gfp_mrna', 'all_mrna', 'rnap_on_gfp',
	'gfp_promoter_copies', 'gfp_rnap_overcrowded_frac',
	'gfp_actual_synth_prob', 'gfp_target_synth_prob', 'ribosomes_on_gfp',
	'ribosomes_on_mrna', 'gfp_ribosome_overcrowded_frac', 'elongation_rate',
	'active_ribosomes', 'active_rnap', 'free_ribosomes', 'ppgpp_conc',
)


def in_window(cycles: pd.DataFrame, first_gen: int = C.ANALYSIS_FIRST_GEN,
		last_gen: int = C.ANALYSIS_LAST_GEN) -> pd.DataFrame:
	"""Cycles inside the analysis window that divided."""
	return cycles[(cycles['generation'] >= first_gen)
		& (cycles['generation'] <= last_gen) & cycles['divided']]


def lineage_table(cycles: pd.DataFrame, gene: GeneInfo,
		first_gen: int = C.ANALYSIS_FIRST_GEN,
		last_gen: int = C.ANALYSIS_LAST_GEN) -> pd.DataFrame:
	"""
	One row per lineage: window means, survival, and growth relative to the
	seed-matched control.
	"""
	last = cycles[cycles['generation'] == last_gen]
	survived = set(map(tuple, last[last['divided']][['variant', 'seed']]
		.values))
	reached = cycles.groupby(['variant', 'seed'])['generation'].max()

	window = in_window(cycles, first_gen, last_gen)
	raw = [c for c in RAW_COLUMNS if c in window]
	lineages = window.groupby(['variant', 'seed', 'promoter', 'rbs'])[raw] \
		.mean().reset_index()
	n_cycles = window.groupby(['variant', 'seed']).size()

	# Lineages with no divided window cycle still count toward n_total.
	everyone = cycles.groupby(['variant', 'seed', 'promoter', 'rbs']).size() \
		.reset_index()[['variant', 'seed', 'promoter', 'rbs']]
	lineages = everyone.merge(lineages, how='left',
		on=['variant', 'seed', 'promoter', 'rbs'])
	keys = list(zip(lineages['variant'], lineages['seed']))
	lineages['n_window_cycles'] = [int(n_cycles.get(k, 0)) for k in keys]
	lineages['last_generation'] = [int(reached.get(k, -1)) for k in keys]
	lineages['survived'] = [k in survived for k in keys]

	lineages = add_derived_columns(lineages, gene)
	control = lineages[lineages['variant'] == CONTROL_VARIANT].set_index(
		'seed')['growth_rate_per_h']
	lineages['control_growth_rate_per_h'] = lineages['seed'].map(control)
	lineages['relative_growth_rate'] = _safe_ratio(
		lineages['growth_rate_per_h'], lineages['control_growth_rate_per_h'])
	return lineages


def relative_cycle_growth(cycles: pd.DataFrame, lineages: pd.DataFrame
		) -> pd.Series:
	"""
	Each cycle's growth rate over its seed-matched control lineage's window
	mean, for per-cycle scatters such as 3H.
	"""
	control = lineages[lineages['variant'] == CONTROL_VARIANT].set_index(
		'seed')['growth_rate_per_h']
	return pd.Series(_safe_ratio(cycles['growth_rate_per_h'],
		cycles['seed'].map(control)), index=cycles.index)


def design_table(lineages: pd.DataFrame, gene: GeneInfo) -> pd.DataFrame:
	"""
	One row per design: means over surviving lineages, with a 95% CI on
	relative growth and the survival counts the hollow markers need.

	Ratios are again recomputed from mean raw quantities.
	"""
	from scipy import stats
	rows = []
	for (variant, promoter, rbs), group in lineages.groupby(
			['variant', 'promoter', 'rbs']):
		alive = group[group['survived']]
		row = dict(variant=variant, promoter=promoter, rbs=rbs,
			n_lineages_total=len(group), n_lineages_survived=len(alive))
		raw = [c for c in RAW_COLUMNS if c in alive]
		row.update(alive[raw].mean().to_dict() if len(alive)
			else {c: np.nan for c in raw})
		rel = alive['relative_growth_rate'].dropna().values
		row['relative_growth_rate'] = float(np.mean(rel)) if rel.size else (
			np.nan)
		if rel.size > 1:
			half = stats.t.ppf(0.975, rel.size - 1) * rel.std(ddof=1) / (
				np.sqrt(rel.size))
			row['relative_growth_rate_ci95'] = float(half)
		else:
			row['relative_growth_rate_ci95'] = np.nan
		rows.append(row)
	designs = add_derived_columns(pd.DataFrame(rows), gene)
	designs['failed'] = designs['n_lineages_survived'] < designs[
		'n_lineages_total']
	return designs.sort_values('variant').reset_index(drop=True)


def lost_designs(designs: pd.DataFrame) -> str:
	"""Designs with no surviving lineage, for a panel's provenance notes."""
	lost = designs[designs['n_lineages_survived'] == 0]
	if lost.empty:
		return 'none'
	return ', '.join('promoter %g RBS %g' % (p, r)
		for p, r in zip(lost['promoter'], lost['rbs']))


# --------------------------------------------------------------------------
# Per-gene profiles, for the transcriptome panel
def _physical_mrna_matrix(path):
	"""(time x mRNA TU) physical counts and their ids, or (None, None)."""
	cell = _Cell(path)
	ids = cell.attribute('RNACounts', 'mRNA_ids')
	reader = cell.reader('RNACounts')
	if ids is None or reader is None:
		return None, None
	total = None
	for name in ('mRNA_counts', 'decaying_mRNA_counts',
			'decaying_cleared_mRNA_counts'):
		if name in reader.columnNames():
			values = reader.readColumn(name).astype(float)
			total = values if total is None else total + values
	return total, list(ids)


def _cycle_tu_means(path):
	try:
		counts, ids = _physical_mrna_matrix(path)
	except Exception:                                       # noqa: BLE001
		return None, None
	if counts is None or counts.size == 0:
		return None, None
	return counts.mean(axis=0), ids


def cistron_profiles(batch_dir: str, sim_data, data: 'PaperData',
		variants: Sequence[int], plot_out_dir: str, cpus: int = 1
		) -> pd.DataFrame:
	"""
	Mean physical mRNA copies per cistron for each requested variant.

	Each divided window cycle gives a time-mean count per mRNA transcription
	unit (summing every RNACounts mRNA field, as everywhere else); these are
	mapped to cistrons through sim_data's cistron-TU matrix, the same mapping
	the RNACounts listener uses, then averaged over a lineage's window and
	over surviving lineages.

	Returns:
		DataFrame indexed by cistron id with columns: gene_id,
		is_ribosomal_protein, is_RNAP, is_new_gene, and one column of mean
		counts per variant. Cached per variant beside the cycle cache.
	"""
	transcription = sim_data.process.transcription
	cistrons = transcription.cistron_data
	mrna_cistron = np.where(cistrons['is_mRNA'])[0]
	rna_index = {rna: i for i, rna in enumerate(transcription.rna_data['id'])}
	mapping = transcription.cistron_tu_mapping_matrix.tocsr()[mrna_cistron]

	out = pd.DataFrame({
		'gene_id': cistrons['gene_id'][mrna_cistron],
		'is_ribosomal_protein': cistrons['is_ribosomal_protein'][mrna_cistron],
		'is_RNAP': cistrons['is_RNAP'][mrna_cistron],
		'is_new_gene': cistrons['is_new_gene'][mrna_cistron],
	}, index=pd.Index(cistrons['id'][mrna_cistron], name='cistron_id'))

	window = in_window(data.cycles)
	alive = data.lineages[data.lineages['survived']]
	cache_dir = os.path.join(plot_out_dir, CACHE_SUBDIR)
	os.makedirs(cache_dir, exist_ok=True)
	for variant in variants:
		seeds = set(alive[alive['variant'] == variant]['seed'])
		rows = window[(window['variant'] == variant)
			& window['seed'].isin(seeds)]
		paths = sorted(rows['path'])
		cache = os.path.join(cache_dir, 'cistrons_%06d.npz' % variant)
		key = np.array([str(CACHE_VERSION)] + paths)
		profile = None
		if os.path.exists(cache) and not os.environ.get(REBUILD_ENV):
			stored = np.load(cache, allow_pickle=False)
			if np.array_equal(stored['key'], key):
				profile = stored['profile']
		if profile is None:
			profile = _variant_profile(rows, mapping, rna_index, cpus)
			tmp = cache + '.%d.tmp.npz' % os.getpid()
			np.savez(tmp, key=key, profile=profile)
			os.replace(tmp, cache)
		out[variant] = profile
	return out


def _variant_profile(rows, mapping, rna_index, cpus):
	"""Mean over surviving lineages of each lineage's window-mean profile."""
	paths = list(rows['path'])
	if not paths:
		return np.full(mapping.shape[0], np.nan)
	if cpus > 1 and len(paths) > 1:
		pool = parallelization.pool(min(cpus, len(paths)))
		results = pool.map(_cycle_tu_means, paths)
		pool.close()
		pool.join()
	else:
		results = [_cycle_tu_means(p) for p in paths]
	per_lineage = {}
	for seed, (means, ids) in zip(rows['seed'], results):
		if means is None:
			continue
		tu = np.zeros(mapping.shape[1])
		tu[[rna_index[i] for i in ids]] = means
		per_lineage.setdefault(seed, []).append(mapping.dot(tu))
	if not per_lineage:
		return np.full(mapping.shape[0], np.nan)
	return np.mean([np.mean(v, axis=0) for v in per_lineage.values()],
		axis=0)


# --------------------------------------------------------------------------
# One call for a panel
@dataclass
class PaperData:
	"""Everything a panel needs, built once per analysis run."""
	gene: GeneInfo
	cycles: pd.DataFrame
	lineages: pd.DataFrame
	designs: pd.DataFrame

	def window_cycles(self) -> pd.DataFrame:
		"""Divided cycles in the analysis window, with relative growth."""
		window = in_window(self.cycles).copy()
		window['relative_growth_rate'] = relative_cycle_growth(
			window, self.lineages)
		return window

	def provenance(self, variants: Optional[Sequence[int]] = None) -> dict:
		"""Keyword arguments for gfp_paper_style.save_panel."""
		lineages = self.lineages
		if variants is not None:
			lineages = lineages[lineages['variant'].isin(variants)]
		return dict(
			variants=sorted(int(v) for v in lineages['variant'].unique()),
			seeds=sorted(int(s) for s in lineages['seed'].unique()),
			generation_window=(C.ANALYSIS_FIRST_GEN, C.ANALYSIS_LAST_GEN),
			n_lineages_total=int(len(lineages)),
			n_lineages_survived=int(lineages['survived'].sum()),
		)


def load(batch_dir: str, sim_data_path: str, metadata: dict,
		plot_out_dir: str, cpus: int = 1) -> PaperData:
	"""Build the three tables for a batch, reading from cache when possible."""
	with open(sim_data_path, 'rb') as f:
		sim_data = pickle.load(f)
	gene = gene_info(sim_data)
	cycles = cycle_table(batch_dir, sim_data, metadata, plot_out_dir, cpus)
	lineages = lineage_table(cycles, gene)
	designs = design_table(lineages, gene)
	missing: set = set()
	for value in cycles.get('missing_columns', pd.Series(dtype=str)).dropna():
		missing.update(x for x in str(value).split(';') if x)
	if missing:
		warnings.warn('this batch lacks listener columns, so the derived'
			' quantities that need them are NaN: %s' % ', '.join(
			sorted(missing)))
	return PaperData(gene, cycles, lineages, designs)
