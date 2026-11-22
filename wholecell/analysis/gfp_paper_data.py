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
CACHE_VERSION = 5
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


@dataclass(frozen=True)
class MachineryInfo:
	"""
	Ids for the host machinery Figure 5 follows, each in the index space of
	the listener that reports it.
	"""
	rpoc_gene_id: str               # RnaSynthProb gene_ids (gene copies)
	rpoc_cistron_id: str            # RnapData cistron_ids (transcription)
	rpoc_monomer_id: str            # MonomerCounts / RibosomeData monomerIds
	full_rnap_id: str               # free RNAP in BulkMolecules
	rnap_formation_reactions: Tuple[str, ...]   # complexation reaction ids
	# (label, complex id, subunit ids, stoichiometry) for limiting-subunit
	# analysis; subunits include rRNAs for the ribosomal subunits
	complexes: Tuple[Tuple[str, str, Tuple[str, ...], Tuple[float, ...]], ...]


RPOC_MONOMER_ID = 'RPOC-MONOMER[c]'


def machinery_info(sim_data) -> MachineryInfo:
	"""Read the machinery ids from sim_data, checking each exists."""
	translation = sim_data.process.translation
	transcription = sim_data.process.transcription
	complexation = sim_data.process.complexation
	full_rnap = sim_data.molecule_ids.full_RNAP

	rnap = complexation.get_monomers(full_rnap)
	if RPOC_MONOMER_ID not in list(rnap['subunitIds']):
		raise ValueError('%s is not a subunit of %s; the RNAP subunit ids'
			' have changed' % (RPOC_MONOMER_ID, full_rnap))
	monomers = list(translation.monomer_data['id'])
	cistron = translation.monomer_data['cistron_id'][
		monomers.index(RPOC_MONOMER_ID)]
	cistrons = list(transcription.cistron_data['id'])
	gene_id = transcription.cistron_data['gene_id'][cistrons.index(cistron)]

	stoich = complexation.stoich_matrix()
	row = list(complexation.molecule_names).index(full_rnap)
	reactions = tuple(complexation.ids_reactions[i]
		for i in np.where(stoich[row, :] > 0)[0])
	if not reactions:
		raise ValueError('no complexation reaction produces %s' % full_rnap)

	complexes = []
	for label, complex_id in (('30S', sim_data.molecule_ids.s30_full_complex),
			('50S', sim_data.molecule_ids.s50_full_complex),
			('RNAP', full_rnap)):
		parts = complexation.get_monomers(complex_id)
		complexes.append((label, complex_id,
			tuple(str(x) for x in parts['subunitIds']),
			tuple(float(x) for x in parts['subunitStoich'])))
	return MachineryInfo(rpoc_gene_id=str(gene_id),
		rpoc_cistron_id=str(cistron), rpoc_monomer_id=RPOC_MONOMER_ID,
		full_rnap_id=full_rnap, rnap_formation_reactions=reactions,
		complexes=tuple(complexes))


def replichore_fraction(sim_data, coords) -> np.ndarray:
	"""
	Position along the replichore, 0 at oriC and 1 at terC, for chromosome
	coordinates relative to oriC (as in replication_coordinate).

	Each arm is normalised by its own length, as Cooper-Helmstetter
	(replication.get_average_copy_number) does.
	"""
	right, left = sim_data.process.replication.replichore_lengths
	coords = np.asarray(coords, dtype=float)
	return np.where(coords >= 0, coords / right, -coords / left)


def nearest_design(designs: pd.DataFrame, share: float) -> pd.Series:
	"""
	The induced design with survivors whose mean GFP proteome share is
	closest to `share`, so "GFP 35%" means the same design in every panel.

	Raises:
		ValueError: if no induced design has a surviving lineage.
	"""
	candidates = designs[(designs['promoter'] > 0)
		& (designs['n_lineages_survived'] > 0)
		& np.isfinite(designs['gfp_proteome_share'])]
	if candidates.empty:
		raise ValueError('no induced design has a surviving lineage')
	gap = (candidates['gfp_proteome_share'] - share).abs()
	return candidates.loc[gap.idxmin()]


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
		description, which is what actually ran. NEW_GENE_DESIGNS is used
		only when no description was recorded; when the two disagree (the
		list was edited for a later batch) the recorded design wins, with a
		warning.
	"""
	design_index = variant % 1000
	if design_index == CONTROL_VARIANT:
		return 0.0, 0.0
	name = metadata.get('variant', '')
	if name == 'new_gene_internal_shift_with_kcat':
		from models.ecoli.sim.variants.new_gene_internal_shift_with_kcat \
			import WILDTYPE_CONTROL_INDEX
		if design_index == WILDTYPE_CONTROL_INDEX:
			return 0.0, 0.0

	if name == 'new_gene_design_list':
		from models.ecoli.sim.variants.new_gene_design_list import (
			NEW_GENE_DESIGNS)
		recorded = _recorded_design(batch_dir, variant)
		current = None
		if design_index < len(NEW_GENE_DESIGNS):
			expression, trl_eff = NEW_GENE_DESIGNS[design_index]
			current = (float(expression), float(trl_eff))
		if recorded is not None:
			if current is not None and not np.allclose(recorded, current):
				warnings.warn('variant %d ran as %s, but NEW_GENE_DESIGNS now'
					' says %s (edited for a later batch); using what ran'
					% (variant, recorded, current))
			return recorded
		if current is None:
			raise ValueError('variant %d is outside NEW_GENE_DESIGNS and'
				' has no recorded description' % variant)
		return current

	factors = metadata.get('new_gene_expression_factors')
	efficiencies = metadata.get('new_gene_translation_efficiency_values')
	if (factors is None or efficiencies is None) and not name.startswith(
			'new_gene'):
		# A batch without new gene designs (wildtype, kcat scans): every
		# variant is a no-GFP condition.
		return 0.0, 0.0
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

	def column(self, listener, name, indices=None):
		"""
		The column, or None (recorded in `missing`) if it is absent.

		`indices` reads only those subcolumns, which matters for wide columns
		such as BulkMolecules counts, where reading every molecule to use two
		dominates the cost of a cycle.
		"""
		reader = self.reader(listener)
		if reader is None or name not in reader.columnNames():
			self.missing.append('%s/%s' % (listener, name))
			return None
		if indices is not None:
			return reader.readColumn(name, indices=np.asarray(indices),
				squeeze=False)
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


def read_cycle(path: str, gene: GeneInfo, machinery: MachineryInfo
		) -> Dict[str, object]:
	"""
	Time means of the raw quantities for one cell cycle.

	Only raw means are stored; every ratio is derived later in
	add_derived_columns, so a change of definition does not need a re-read.
	"""
	cell = _Cell(path)
	row: Dict[str, object] = {}

	time = cell.column('Main', 'time')
	row['n_steps'] = 0 if time is None else int(time.size)
	doubling_time_min = (float(time[-1] - time[0]) / 60.0
		if time is not None and time.size > 1 else np.nan)
	row['doubling_time_min'] = doubling_time_min
	row['divided'] = os.path.exists(
		os.path.join(path, 'Daughter1_inherited_state.cPickle'))

	growth = cell.column('Mass', 'instantaneous_growth_rate')
	# The first step is recorded before any growth has happened.
	row['growth_rate_per_h'] = (_mean(growth[1:]) * SECONDS_PER_HOUR
		if growth is not None and growth.size > 1 else np.nan)
	for name in ('cellMass', 'dryMass', 'proteinMass', 'mRnaMass',
			'rRnaMass', 'cellVolume'):
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
	# Every active ribosome counted by the protein it is making. Unlike the
	# mRNA-attached count above, this includes ribosomes finishing on gfp
	# mRNAs that have left the mRNA list, so it sums to the active ribosome
	# count across monomers and is the basis of the GFP ribosome share.
	row['ribosomes_making_gfp'] = (_mean(cell.column('RibosomeData',
		'n_ribosomes_per_transcript', indices=[ribo_monomer]))
		if ribo_monomer is not None else np.nan)
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
	bulk_names = cell.attribute('BulkMolecules', 'objectNames')
	wanted = [_index(bulk_names, m) for m in (gene.s30_id, gene.s50_id,
		machinery.full_rnap_id)]
	bulk = (cell.column('BulkMolecules', 'counts', indices=wanted)
		if None not in wanted else None)
	row['free_ribosomes'] = (_mean(np.minimum(bulk[:, 0], bulk[:, 1]))
		if bulk is not None else np.nan)
	row['free_rnap'] = _mean(bulk[:, 2]) if bulk is not None else np.nan

	# rpoC and RNAP capacity (Figure 5G, 5H). Event counts become rates per
	# second over the cycle, so cycles of different length compare.
	seconds = doubling_time_min * 60.0

	def rate(listener, column, ids_attribute=None, wanted_id=None,
			ids_list=None):
		if ids_attribute is not None:
			i = _index(cell.attribute(listener, ids_attribute), wanted_id)
			if i is None:
				return np.nan
			values = cell.column(listener, column, indices=[i])
		elif ids_list is not None:
			ids = cell.attribute(listener, 'reactionIDs')
			idx = [_index(ids, r) for r in ids_list]
			if ids is None or None in idx:
				return np.nan
			values = cell.column(listener, column, indices=idx)
		else:
			values = cell.column(listener, column)
		if values is None or not np.isfinite(seconds) or seconds <= 0:
			return np.nan
		return float(np.nansum(values)) / seconds

	g_rpoc = _index(cell.attribute('RnaSynthProb', 'gene_ids'),
		machinery.rpoc_gene_id)
	row['rpoc_gene_copies'] = (_mean(cell.column('RnaSynthProb',
		'gene_copy_number', indices=[g_rpoc])) if g_rpoc is not None
		else np.nan)
	row['rpoc_transcription_per_s'] = rate('RnapData',
		'rna_init_event_per_cistron', 'cistron_ids', machinery.rpoc_cistron_id)
	row['rpoc_translation_per_s'] = rate('RibosomeData',
		'ribosome_init_event_per_monomer', 'monomerIds',
		machinery.rpoc_monomer_id)
	m_rpoc = _index(cell.attribute('MonomerCounts', 'monomerIds'),
		machinery.rpoc_monomer_id)
	row['rpoc_monomers'] = (_mean(monomers[:, m_rpoc])
		if monomers is not None and m_rpoc is not None else np.nan)
	# complexationEvents is indexed by the listener's reactionIDs, which is
	# not the order of its complexIDs attribute; resolve by reaction id.
	row['rnap_formed_per_s'] = rate('ComplexationListener',
		'complexationEvents', ids_list=machinery.rnap_formation_reactions)
	row['rnap_initiations_per_s'] = rate('RnapData', 'didInitialize')
	# gfp's own initiation rates (Figure 6C, 6D), in each listener's index
	row['gfp_transcription_per_s'] = rate('RnapData', 'rnaInitEvent',
		'rnaIds', gene.rna_id)
	row['gfp_translation_per_s'] = rate('RibosomeData',
		'ribosome_init_event_per_monomer', 'monomerIds', gene.monomer_id)

	row['ppgpp_conc'] = _mean(cell.column('GrowthLimits', 'ppgpp_conc'))
	row['oric'] = _mean(cell.column('ReplicationData', 'numberOfOric'))
	# Metabolic enzymes: every FBA catalyst molecule (complexes counted once),
	# summed, for the total enzyme concentration in Figure 1E.
	catalysts = cell.column('FBAResults', 'catalyst_counts')
	row['enzyme_counts'] = (_mean(catalysts.sum(axis=1))
		if catalysts is not None and catalysts.ndim == 2 else np.nan)
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
		df['total_rnap'] = df['active_rnap'] + df['free_rnap']
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
			df['ribosomes_making_gfp'], df['active_ribosomes'])
		# Concentrations in micromolar: count / (Avogadro x volume)
		df['total_ribosomes_uM'] = _safe_ratio(df['total_ribosomes'],
			AVOGADRO * 1e-15 * df['cellVolume'] * 1e-6)
		df['total_rnap_uM'] = _safe_ratio(df['total_rnap'],
			AVOGADRO * 1e-15 * df['cellVolume'] * 1e-6)
		df['enzyme_uM'] = _safe_ratio(df['enzyme_counts'],
			AVOGADRO * 1e-15 * df['cellVolume'] * 1e-6)
		# Amino acids incorporated per second per fL: the translational demand
		# on metabolism (Figure 4G), a concentration rate like enzyme supply.
		df['translation_demand'] = _safe_ratio(
			df['active_ribosomes'] * df['elongation_rate'], df['cellVolume'])
		df['gfp_transcription_share'] = _safe_ratio(
			df['gfp_transcription_per_s'], df['rnap_initiations_per_s'])
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
	path, gene_dict, machinery_dict = args
	try:
		return read_cycle(path, GeneInfo(**gene_dict),
			MachineryInfo(**machinery_dict))
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
	machinery = machinery_info(sim_data)
	cache_dir = os.path.join(plot_out_dir, CACHE_SUBDIR)
	os.makedirs(cache_dir, exist_ok=True)
	cache_path = os.path.join(cache_dir, 'cycles.csv')
	stamp_path = os.path.join(cache_dir, 'cycles.json')
	# Round-tripped through JSON so tuples compare equal to the lists a
	# stored stamp reads back as.
	stamp = json.loads(json.dumps({'cache_version': CACHE_VERSION,
		'gene': asdict(gene), 'machinery': asdict(machinery)}))

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
		jobs = [(path, asdict(gene), asdict(machinery))
			for _, _, _, path in todo]
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
	'free_rnap', 'rpoc_gene_copies', 'rpoc_transcription_per_s',
	'rpoc_translation_per_s', 'rpoc_monomers', 'rnap_formed_per_s',
	'rnap_initiations_per_s', 'cellVolume', 'ribosomes_making_gfp', 'oric',
	'enzyme_counts', 'gfp_transcription_per_s', 'gfp_translation_per_s',
)
# Derived in add_derived_columns from RAW_COLUMNS; listed for readers.


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


def seed_matched_ratio(lineages: pd.DataFrame, variant: int, column: str,
		baseline: int = CONTROL_VARIANT) -> dict:
	"""
	A lineage quantity in one design over the same seed's baseline lineage.

	Only surviving design lineages with a baseline lineage are used. Returns
	the mean ratio, a t-based 95% CI half-width, the per-seed ratios and n.
	"""
	from scipy import stats
	base = lineages[lineages['variant'] == baseline].set_index('seed')[column]
	rows = lineages[(lineages['variant'] == variant) & lineages['survived']]
	ratios = pd.Series(_safe_ratio(rows[column].values,
		rows['seed'].map(base).values), index=rows['seed'].values).dropna()
	n = int(ratios.size)
	half = (float(stats.t.ppf(0.975, n - 1) * ratios.std(ddof=1) / np.sqrt(n))
		if n > 1 else np.nan)
	return dict(mean=float(ratios.mean()) if n else np.nan, ci95=half,
		per_seed=ratios, n=n)


def _ci95(values) -> float:
	"""t-based 95% CI half-width of a mean, NaN for fewer than two values."""
	from scipy import stats
	values = np.asarray(values, dtype=float)
	if values.size < 2:
		return np.nan
	return float(stats.t.ppf(0.975, values.size - 1) * values.std(ddof=1)
		/ np.sqrt(values.size))


def design_table(lineages: pd.DataFrame, gene: GeneInfo) -> pd.DataFrame:
	"""
	One row per design: means over surviving lineages, with a 95% CI on
	relative growth and the survival counts the hollow markers need.

	Ratios are again recomputed from mean raw quantities.
	"""
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
		row['relative_growth_rate_ci95'] = _ci95(rel)
		# The doubling-time versions of growth panels (DOUBLING_TIME_PART)
		# show absolute minutes, so their interval is across lineages too.
		row['doubling_time_min_ci95'] = _ci95(
			alive['doubling_time_min'].dropna().values)
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
# Per-gene profiles: one value per gene (or TU, or monomer) per design
#
# Each kind names a per-gene listener column and the attribute holding its
# ids. Readers are module-level so a process pool can pickle them.
VECTOR_KINDS = {
	# physical mRNA copies per mRNA transcription unit
	'mrna_tu': ('RNACounts', None, 'mRNA_ids'),
	# gene copies (dosage) per gene, the replication-driven quantity
	'gene_copies': ('RnaSynthProb', 'gene_copy_number', 'gene_ids'),
	# protein copies per monomer
	'monomers': ('MonomerCounts', 'monomerCounts', 'monomerIds'),
	# RNAP on each mRNA transcription unit: transcripts still being made
	'rnap_mrna_tu': ('RNACounts', None, 'mRNA_ids'),
	# RNAP on each rRNA operon
	'rnap_rrna_tu': ('RNACounts', 'partial_rRNA_counts', 'rRNA_ids'),
	# active ribosomes by the protein they are making
	'ribosomes_monomer': ('RibosomeData', 'n_ribosomes_per_transcript',
		'monomerIds'),
	# FBA catalysts (enzymes and enzyme complexes), molecules per cell
	'catalysts': ('FBAResults', 'catalyst_counts', 'catalyst_ids'),
	# metabolic flux per base reaction, mmol per g dry weight per hour;
	# `extra` is (cell density in g/L,)
	'base_fluxes': ('FBAResults', 'base_reaction_fluxes',
		'base_reaction_ids'),
	# promoter copies per transcription unit (TU-indexed dosage)
	'promoter_copies': ('RnaSynthProb', 'promoter_copy_number', 'rnaIds'),
	# limiting-subunit credit for one complex; `extra` is its subunits
	'limiting': ('BulkMolecules', 'counts', 'objectNames'),
}


def _read_vector(args):
	"""(time-mean vector, ids) of one per-gene column for one cycle."""
	kind, path, extra = args
	listener, column, id_attribute = VECTOR_KINDS[kind]
	try:
		if kind == 'mrna_tu':
			values, ids = _physical_mrna_matrix(path)
		elif kind == 'rnap_mrna_tu':
			cell = _Cell(path)
			reader = cell.reader('RNACounts')
			name = next((n for n in ('rnap_engaged_mRNA_counts',
				'partial_mRNA_counts') if reader is not None
				and n in reader.columnNames()), None)
			values = None if name is None else reader.readColumn(name)
			ids = cell.attribute('RNACounts', 'mRNA_ids')
			ids = None if ids is None else list(ids)
		elif kind == 'limiting':
			return _limiting_credit(path, *extra)
		elif kind == 'base_fluxes':
			return _base_fluxes(path, *extra)
		else:
			cell = _Cell(path)
			values = cell.column(listener, column)
			ids = cell.attribute(listener, id_attribute)
			ids = None if ids is None else list(ids)
	except Exception:                                       # noqa: BLE001
		return None, None
	if values is None or ids is None or values.size == 0:
		return None, None
	return values.astype(float).mean(axis=0), ids


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


@dataclass
class VectorProfile:
	"""A per-gene quantity for one design, kept per lineage for CIs."""
	ids: List[str]
	per_lineage: Dict[int, np.ndarray]   # seed -> window-mean vector

	@property
	def mean(self) -> np.ndarray:
		"""Mean over surviving lineages (NaN if there are none)."""
		if not self.per_lineage:
			return np.full(len(self.ids), np.nan)
		return np.mean(list(self.per_lineage.values()), axis=0)


def vector_profiles(data: 'PaperData', kind: str, variants: Sequence[int],
		plot_out_dir: str, cpus: int = 1, extra: Optional[tuple] = None,
		label: Optional[str] = None) -> Dict[int, VectorProfile]:
	"""
	Window means of a per-gene column for each variant, per surviving lineage.

	Every divided window cycle gives a time-mean vector; a lineage's value is
	the mean of its cycles, and a design's (VectorProfile.mean) is the mean
	of its surviving lineages, the same rule as the scalar tables. Results
	are cached per (kind, variant) and rebuilt whenever the set of cycles
	behind them changes.

	Args:
		data: the batch's PaperData
		kind: a key of VECTOR_KINDS
		variants: variant indices to profile
		plot_out_dir: the analysis plotOutDir, which holds the cache
		cpus: processes for reading
		extra: kind-specific arguments passed to the reader (for 'limiting',
			the complex's subunit ids and stoichiometry)
		label: distinguishes caches of one kind with different `extra`
	"""
	if kind not in VECTOR_KINDS:
		raise ValueError('unknown vector kind %r; expected one of %s'
			% (kind, sorted(VECTOR_KINDS)))
	window = in_window(data.cycles)
	alive = data.lineages[data.lineages['survived']]
	cache_dir = os.path.join(plot_out_dir, CACHE_SUBDIR)
	os.makedirs(cache_dir, exist_ok=True)
	out = {}
	for variant in variants:
		seeds = set(alive[alive['variant'] == variant]['seed'])
		rows = window[(window['variant'] == variant)
			& window['seed'].isin(seeds)].sort_values(['seed', 'generation'])
		name = kind if label is None else '%s_%s' % (kind, label)
		key = np.array([str(CACHE_VERSION), name, repr(extra)]
			+ list(rows['path']))
		cache = os.path.join(cache_dir, '%s_%06d.npz' % (name, variant))
		profile = None
		if os.path.exists(cache) and not os.environ.get(REBUILD_ENV):
			stored = np.load(cache, allow_pickle=False)
			if np.array_equal(stored['key'], key):
				profile = VectorProfile(list(stored['ids']), dict(zip(
					stored['seeds'].tolist(), stored['values'])))
		if profile is None:
			profile = _build_profile(kind, rows, cpus, extra)
			seeds_arr = np.array(sorted(profile.per_lineage), dtype=int)
			values = (np.array([profile.per_lineage[s] for s in seeds_arr])
				if len(seeds_arr) else np.zeros((0, len(profile.ids))))
			tmp = cache + '.%d.tmp.npz' % os.getpid()
			np.savez(tmp, key=key, ids=np.array(profile.ids), seeds=seeds_arr,
				values=values)
			os.replace(tmp, cache)
		out[variant] = profile
	return out


def _build_profile(kind, rows, cpus, extra=None) -> VectorProfile:
	paths = list(rows['path'])
	jobs = [(kind, p, extra) for p in paths]
	if cpus > 1 and len(jobs) > 1:
		pool = parallelization.pool(min(cpus, len(jobs)))
		results = pool.map(_read_vector, jobs)
		pool.close()
		pool.join()
	else:
		results = [_read_vector(j) for j in jobs]
	ids: Optional[List[str]] = None
	per_seed: Dict[int, list] = {}
	for seed, (means, cycle_ids) in zip(rows['seed'], results):
		if means is None:
			continue
		if ids is None:
			ids = cycle_ids
		elif cycle_ids != ids:
			raise ValueError('%s ids differ between cycles of one batch'
				% kind)
		per_seed.setdefault(int(seed), []).append(means)
	return VectorProfile(ids or [], {seed: np.mean(v, axis=0)
		for seed, v in per_seed.items()})


def cistron_profiles(sim_data, data: 'PaperData', variants: Sequence[int],
		plot_out_dir: str, cpus: int = 1) -> pd.DataFrame:
	"""
	Mean physical mRNA copies per cistron for each requested variant.

	Physical copies per mRNA transcription unit (vector_profiles 'mrna_tu')
	are mapped to cistrons through sim_data's cistron-TU matrix, the same
	mapping the RNACounts listener uses. The mapping is linear, so mapping the
	design mean equals averaging mapped lineages.

	Returns:
		DataFrame indexed by cistron id with columns: gene_id,
		is_ribosomal_protein, is_RNAP, is_new_gene, and one column of mean
		counts per variant.
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
	profiles = vector_profiles(data, 'mrna_tu', variants, plot_out_dir, cpus)
	for variant, profile in profiles.items():
		tu = np.zeros(mapping.shape[1])
		if profile.ids:
			tu[[rna_index[i] for i in profile.ids]] = profile.mean
		else:
			tu[:] = np.nan
		out[variant] = mapping.dot(tu)
	return out


# --------------------------------------------------------------------------
# Limiting subunit (Figure 5F)
TIED = '_tied'          # extra id: share of time steps with a tie for limiting


def _limiting_credit(path, subunit_ids, stoich):
	"""
	Time-mean share of "limiting" credit per subunit of one complex.

	At each time step the limiting subunits are those with the fewest free
	copies per copy needed (free / stoichiometry, from BulkMolecules). When
	several subunits tie, the credit is split equally among them rather than
	given to whichever comes first in the list: the rRNAs come last, so first-
	in-list credit would assign every protein-rRNA tie (common at zero) to
	protein. The share of tied time steps is returned as the extra id TIED.
	"""
	cell = _Cell(path)
	names = cell.attribute('BulkMolecules', 'objectNames')
	idx = [_index(names, m) for m in subunit_ids]
	if names is None or None in idx:
		return None, None
	free = cell.column('BulkMolecules', 'counts', indices=idx)
	if free is None or free.size == 0:
		return None, None
	credit, tied_share = split_limiting_credit(free, stoich)
	return np.append(credit, tied_share), list(subunit_ids) + [TIED]


def split_limiting_credit(free, stoich) -> Tuple[np.ndarray, float]:
	"""
	(time-mean credit per subunit, share of tied time steps).

	Args:
		free: (time x subunit) free copies
		stoich: copies of each subunit one complex needs
	"""
	per_copy = np.asarray(free, dtype=float) / np.asarray(stoich, dtype=float)
	tied = per_copy == per_copy.min(axis=1, keepdims=True)
	n_tied = tied.sum(axis=1, keepdims=True)
	return (tied / n_tied).mean(axis=0), float(np.mean(n_tied[:, 0] > 1))


def is_rrna(molecule_id: str) -> bool:
	"""Whether a complex subunit is an rRNA rather than a protein."""
	return 'RRNA' in molecule_id.upper()


def limiting_profiles(data: 'PaperData', variants: Sequence[int],
		plot_out_dir: str, cpus: int = 1) -> Dict[str, Dict[int, VectorProfile]]:
	"""{complex label: {variant: credit per subunit}} for 30S, 50S and RNAP."""
	if data.machinery is None:
		raise ValueError('PaperData has no machinery ids')
	out = {}
	for label, _, subunits, stoich in data.machinery.complexes:
		out[label] = vector_profiles(data, 'limiting', variants,
			plot_out_dir, cpus, extra=(tuple(subunits), tuple(stoich)),
			label=label)
	return out


# --------------------------------------------------------------------------
# Event-anchored responses around replication initiation (Figure 5D)
EVENT_COLUMNS = ('rpoc_copies', 'rpoc_transcription', 'rpoc_translation',
	'terminus_transcription')


def terminus_tus(sim_data) -> List[str]:
	"""
	Transcription units near the terminus, the control for dosage effects:
	replichore fraction above TERMINUS_MIN_FRACTION, excluding rRNA, RNAP and
	ribosomal-protein units and any unit carrying the inserted gene.
	"""
	rna = sim_data.process.transcription.rna_data
	transcription = sim_data.process.transcription
	new_gene_tus = set()
	cistrons = transcription.cistron_data
	for cistron in cistrons['id'][cistrons['is_new_gene']]:
		new_gene_tus.update(transcription.cistron_id_to_rna_indexes(cistron))
	f = replichore_fraction(sim_data, rna['replication_coordinate'])
	keep = ((f > C.TERMINUS_MIN_FRACTION) & ~rna['is_rRNA']
		& ~rna['includes_RNAP'] & ~rna['includes_ribosomal_protein'])
	return [str(rna['id'][i]) for i in np.where(keep)[0]
		if i not in new_gene_tus]


def _cycle_events(args):
	"""Per-window sums around every replication initiation in one cycle."""
	path, machinery_dict, terminus = args
	machinery = MachineryInfo(**machinery_dict)
	try:
		cell = _Cell(path)
		time = cell.column('Main', 'time')
		oric = cell.column('ReplicationData', 'numberOfOric')
		if time is None or oric is None or time.size < 2:
			return []
		g = _index(cell.attribute('RnaSynthProb', 'gene_ids'),
			machinery.rpoc_gene_id)
		c = _index(cell.attribute('RnapData', 'cistron_ids'),
			machinery.rpoc_cistron_id)
		m = _index(cell.attribute('RibosomeData', 'monomerIds'),
			machinery.rpoc_monomer_id)
		tus = cell.attribute('RnapData', 'rnaIds')
		t_idx = [_index(tus, t) for t in terminus]
		if None in (g, c, m) or tus is None or None in t_idx:
			return []
		series = {
			'rpoc_copies': cell.column('RnaSynthProb', 'gene_copy_number',
				indices=[g])[:, 0],
			'rpoc_transcription': cell.column('RnapData',
				'rna_init_event_per_cistron', indices=[c])[:, 0],
			'rpoc_translation': cell.column('RibosomeData',
				'ribosome_init_event_per_monomer', indices=[m])[:, 0],
			'terminus_transcription': cell.column('RnapData',
				'rnaInitEvent', indices=t_idx).sum(axis=1),
		}
	except Exception:                                       # noqa: BLE001
		return []
	rows = []
	windows = [('pre', C.EVENT_PRE_WINDOW_MIN)] + [
		('post_%g_%g' % w, w) for w in C.EVENT_POST_WINDOWS_MIN]
	for event, step in enumerate(np.where(np.diff(oric) > 0)[0] + 1):
		minutes = (time - time[step]) / 60.0
		event_rows = []
		for name, (lo, hi) in windows:
			inside = (minutes >= lo) & (minutes < hi)
			n = int(inside.sum())
			if n < C.EVENT_MIN_STEPS:
				event_rows = None
				break
			row = dict(event=event, window=name, window_lo=lo, window_hi=hi,
				n_steps=n)
			for key, values in series.items():
				row[key] = float(np.sum(values[inside]))
			event_rows.append(row)
		if event_rows:       # only events with every window complete
			rows.extend(event_rows)
	return rows


def initiation_events(data: 'PaperData', sim_data, variants: Sequence[int],
		plot_out_dir: str, cpus: int = 1) -> pd.DataFrame:
	"""
	Window sums around every replication initiation in the analysis window of
	surviving lineages, one row per (event, window). An event is kept only if
	its pre window and every post window fall inside its own cell cycle, so
	all windows of one event come from the same cell.

	Cached per variant; rebuilt when the set of cycles changes.
	"""
	if data.machinery is None:
		raise ValueError('PaperData has no machinery ids')
	terminus = terminus_tus(sim_data)
	window = in_window(data.cycles)
	alive = data.lineages[data.lineages['survived']]
	cache_dir = os.path.join(plot_out_dir, CACHE_SUBDIR)
	os.makedirs(cache_dir, exist_ok=True)
	frames = []
	for variant in variants:
		seeds = set(alive[alive['variant'] == variant]['seed'])
		rows = window[(window['variant'] == variant)
			& window['seed'].isin(seeds)].sort_values(['seed', 'generation'])
		cache = os.path.join(cache_dir, 'events_%06d.csv' % variant)
		stamp = cache[:-4] + '.json'
		key = json.loads(json.dumps({'version': CACHE_VERSION,
			'paths': list(rows['path']), 'terminus': terminus,
			'pre': C.EVENT_PRE_WINDOW_MIN, 'post': C.EVENT_POST_WINDOWS_MIN,
			'min_steps': C.EVENT_MIN_STEPS}))
		frame = None
		if (os.path.exists(cache) and os.path.exists(stamp)
				and not os.environ.get(REBUILD_ENV)):
			with open(stamp) as f:
				if json.load(f) == key:
					frame = pd.read_csv(cache)
		if frame is None:
			jobs = [(p, asdict(data.machinery), terminus)
				for p in rows['path']]
			if cpus > 1 and len(jobs) > 1:
				pool = parallelization.pool(min(cpus, len(jobs)))
				results = pool.map(_cycle_events, jobs)
				pool.close()
				pool.join()
			else:
				results = [_cycle_events(j) for j in jobs]
			records = []
			for (seed, gen), events in zip(
					zip(rows['seed'], rows['generation']), results):
				for r in events:
					records.append(dict(variant=variant, seed=int(seed),
						generation=int(gen), **r))
			frame = pd.DataFrame(records, columns=['variant', 'seed',
				'generation', 'event', 'window', 'window_lo', 'window_hi',
				'n_steps'] + list(EVENT_COLUMNS))
			suffix = '.%d.tmp' % os.getpid()
			frame.to_csv(cache + suffix, index=False)
			with open(stamp + suffix, 'w') as f:
				json.dump(key, f)
			os.replace(stamp + suffix, stamp)
			os.replace(cache + suffix, cache)
		frames.append(frame)
	return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()


def pooled_response(events: pd.DataFrame, column: str,
		n_boot: int = C.EVENT_BOOTSTRAP, seed: int = 0) -> pd.DataFrame:
	"""
	Each post window's pooled rate over the pre window's, with a bootstrap CI.

	Pooling (summed counts over summed time steps across every event) is used
	for every quantity, events and copies alike, because single events are
	sparse: rpoC sees well under one transcription event per window, so
	per-event or per-lineage ratios are dominated by zeros. The CI resamples
	lineages with replacement, so events from one lineage stay together and
	the interval reflects lineage-to-lineage variation.

	Returns:
		one row per post window: window, mid (minutes), ratio, lo, hi,
		n_events, n_lineages
	"""
	rng = np.random.default_rng(seed)
	pre = events[events['window'] == 'pre']
	posts = [w for w in events['window'].unique() if w != 'pre']
	lineages = np.array(sorted(events['seed'].unique()))
	sums = events.groupby(['window', 'seed'])[[column, 'n_steps']].sum()

	def by_lineage(window):
		"""(counts, steps) per lineage, zero where a lineage has none."""
		table = sums.loc[window] if window in sums.index.get_level_values(
			0) else pd.DataFrame(columns=[column, 'n_steps'])
		counts = np.array([table[column].get(l, 0.0) for l in lineages],
			dtype=float)
		steps = np.array([table['n_steps'].get(l, 0.0) for l in lineages],
			dtype=float)
		return counts, steps

	pre_counts, pre_steps = by_lineage('pre')
	# Each bootstrap draw is a vector of how many times each lineage is
	# resampled, so every window is resampled with the same lineages.
	draws = (rng.multinomial(lineages.size, np.full(lineages.size,
		1 / lineages.size), size=n_boot) if lineages.size > 1
		else np.zeros((0, lineages.size)))

	def ratio(weights, counts, steps):
		with np.errstate(divide='ignore', invalid='ignore'):
			post_rate = (weights @ counts) / (weights @ steps)
			pre_rate = (weights @ pre_counts) / (weights @ pre_steps)
			return post_rate / pre_rate

	out = []
	for window in posts:
		rows = events[events['window'] == window]
		counts, steps = by_lineage(window)
		point = float(ratio(np.ones(lineages.size), counts, steps))
		boots = np.atleast_1d(ratio(draws, counts, steps)) if len(draws) \
			else np.array([])
		boots = boots[np.isfinite(boots)]
		out.append(dict(window=window,
			mid=float((rows['window_lo'].iloc[0] + rows['window_hi'].iloc[0])
				/ 2),
			ratio=point,
			lo=float(np.percentile(boots, 2.5)) if boots.size else np.nan,
			hi=float(np.percentile(boots, 97.5)) if boots.size else np.nan,
			n_events=int(pre[['seed', 'generation', 'event']]
				.drop_duplicates().shape[0]),
			n_lineages=int(lineages.size)))
	return pd.DataFrame(out).sort_values('mid').reset_index(drop=True)



# --------------------------------------------------------------------------
# Machinery traces through the run (Figure 2E)
TRACE_POINTS = 10       # equal fractions of each cell cycle
# Quantities a trace can follow; each is a per-time-step series in one cell.
# Instantaneous growth rate (dry-mass growth / dry mass) has bursts of one to
# two minutes when metabolism refills small-molecule pools: in one 102-step
# burst small-molecule mass rose ~6x its usual amount while protein mass grew
# as normal. Per-cycle growth averages these out; traces show a centred
# rolling median over this many 1 s steps so the bursts do not dominate.
GROWTH_TRACE_SMOOTH_STEPS = 301
TRACE_SERIES = ('rnap', 'ribosomes', 'gfp_mrna', 'gfp_protein',
	'growth_rate', 'cell_mass', 'oric', 'rrna_copies', 'gfp_copies')


def _trace_series(cell, name, gene_dict, full_rnap):
	"""One per-time-step series for a cell, or None if it cannot be read."""
	if name in ('rnap', 'ribosomes'):
		names = cell.attribute('UniqueMoleculeCounts', 'uniqueMoleculeIds')
		bulk_names = cell.attribute('BulkMolecules', 'objectNames')
		molecule = 'active_RNAP' if name == 'rnap' else 'active_ribosome'
		i = _index(names, molecule)
		free = ([full_rnap] if name == 'rnap'
			else [gene_dict['s30_id'], gene_dict['s50_id']])
		idx = [_index(bulk_names, m) for m in free]
		if i is None or None in idx:
			return None
		active = cell.column('UniqueMoleculeCounts', 'uniqueMoleculeCounts',
			indices=[i])[:, 0]
		bulk = cell.column('BulkMolecules', 'counts', indices=idx)
		return active + (bulk[:, 0] if name == 'rnap' else bulk.min(axis=1))
	if name == 'gfp_mrna':
		counts, ids = _physical_mrna_matrix(cell.path)
		if counts is None or gene_dict['rna_id'] not in ids:
			return None
		return counts[:, ids.index(gene_dict['rna_id'])]
	if name == 'gfp_protein':
		m = _index(cell.attribute('MonomerCounts', 'monomerIds'),
			gene_dict['monomer_id'])
		return None if m is None else cell.column('MonomerCounts',
			'monomerCounts', indices=[m])[:, 0]
	if name == 'growth_rate':
		values = cell.column('Mass', 'instantaneous_growth_rate')
		if values is None:
			return None
		values = values.astype(float) * SECONDS_PER_HOUR
		values[0] = np.nan          # recorded before any growth happened
		return pd.Series(values).rolling(GROWTH_TRACE_SMOOTH_STEPS,
			center=True, min_periods=1).median().values
	if name == 'cell_mass':
		return cell.column('Mass', 'cellMass')
	if name == 'oric':
		return cell.column('ReplicationData', 'numberOfOric')
	if name in ('rrna_copies', 'gfp_copies'):
		ids = cell.attribute('RnaSynthProb', 'rnaIds')
		if ids is None:
			return None
		if name == 'gfp_copies':
			idx = [_index(ids, gene_dict['rna_id'])]
		else:
			# rRNA operons by the naming the RNACounts listener reports them
			# under; summing their promoter copies gives operon copies.
			rrna = cell.attribute('RNACounts', 'rRNA_ids')
			idx = [_index(ids, r) for r in (rrna or [])]
		if not idx or None in idx:
			return None
		return cell.column('RnaSynthProb', 'promoter_copy_number',
			indices=idx).sum(axis=1)
	raise ValueError('unknown trace series %r' % name)


def trace_series(sim_out: str, name: str, gene: GeneInfo,
		machinery: MachineryInfo) -> Optional[np.ndarray]:
	"""One per-time-step series (a TRACE_SERIES name) for a cell."""
	return _trace_series(_Cell(sim_out), name, asdict(gene),
		machinery.full_rnap_id)


def _cycle_trace(args):
	"""Mean of each series in each of n_points fractions of one cycle."""
	path, gene_dict, full_rnap, n_points, series = args
	try:
		cell = _Cell(path)
		values = {name: _trace_series(cell, name, gene_dict, full_rnap)
			for name in series}
	except Exception:                                       # noqa: BLE001
		return []
	lengths = {len(v) for v in values.values() if v is not None}
	if len(lengths) != 1 or lengths.pop() < n_points:
		return []
	n = len(next(v for v in values.values() if v is not None))
	edges = np.linspace(0, n, n_points + 1).astype(int)
	rows = []
	for k, (a, b) in enumerate(zip(edges[:-1], edges[1:])):
		row = {'point': k}
		for name, v in values.items():
			# Growth rate is already smoothed; the bin median keeps any
			# remaining burst from setting the bin. Counts use the bin mean.
			summary = np.nanmedian if name == 'growth_rate' else np.nanmean
			row[name] = (float(summary(v[a:b])) if v is not None
				and np.isfinite(v[a:b]).any() else np.nan)
		rows.append(row)
	return rows


def cycle_traces(data: 'PaperData', variants: Sequence[int],
		plot_out_dir: str, cpus: int = 1, n_points: int = TRACE_POINTS,
		series: Sequence[str] = ('rnap', 'ribosomes')) -> pd.DataFrame:
	"""
	Chosen quantities through every generation of each lineage, as the mean
	over each of n_points equal fractions of the cell cycle (the median
	for growth rate, itself smoothed over GROWTH_TRACE_SMOOTH_STEPS).

	Splitting each cycle into equal fractions lines lineages up generation
	by generation even though their cycles differ in length. RNAP and
	ribosomes are totals (active plus free), as in the scalar tables; gfp
	mRNA is the physical count; growth rate is per hour. Every generation is
	included, not just the analysis window.

	Returns:
		columns variant, seed, generation, point, x (generation + fraction of
		the cycle at the bin centre), and one per series. Cached per
		(variant, series, n_points).
	"""
	unknown = set(series) - set(TRACE_SERIES)
	if unknown:
		raise ValueError('unknown trace series %s; expected %s'
			% (sorted(unknown), TRACE_SERIES))
	if data.machinery is None:
		raise ValueError('PaperData has no machinery ids')
	gene_dict = asdict(data.gene)
	cache_dir = os.path.join(plot_out_dir, CACHE_SUBDIR)
	os.makedirs(cache_dir, exist_ok=True)
	tag = '%d_%s' % (n_points, '-'.join(series))
	frames = []
	for variant in variants:
		rows = data.cycles[data.cycles['variant'] == variant].sort_values(
			['seed', 'generation'])
		cache = os.path.join(cache_dir, 'traces_%s_%06d.csv' % (tag, variant))
		stamp = cache[:-4] + '.json'
		key = json.loads(json.dumps({'version': CACHE_VERSION,
			'n_points': n_points, 'series': list(series),
			'growth_smooth': GROWTH_TRACE_SMOOTH_STEPS,
			'paths': list(rows['path'])}))
		frame = None
		if (os.path.exists(cache) and os.path.exists(stamp)
				and not os.environ.get(REBUILD_ENV)):
			with open(stamp) as f:
				if json.load(f) == key:
					frame = pd.read_csv(cache)
		if frame is None:
			jobs = [(p, gene_dict, data.machinery.full_rnap_id, n_points,
				tuple(series)) for p in rows['path']]
			if cpus > 1 and len(jobs) > 1:
				pool = parallelization.pool(min(cpus, len(jobs)))
				results = pool.map(_cycle_trace, jobs)
				pool.close()
				pool.join()
			else:
				results = [_cycle_trace(j) for j in jobs]
			records = [dict(variant=variant, seed=int(seed),
				generation=int(gen), x=gen + (r['point'] + 0.5) / n_points,
				**r) for seed, gen, trace in zip(rows['seed'],
				rows['generation'], results) for r in trace]
			frame = pd.DataFrame(records, columns=['variant', 'seed',
				'generation', 'x', 'point'] + list(series))
			suffix = '.%d.tmp' % os.getpid()
			frame.to_csv(cache + suffix, index=False)
			with open(stamp + suffix, 'w') as f:
				json.dump(key, f)
			os.replace(stamp + suffix, stamp)
			os.replace(cache + suffix, cache)
		frames.append(frame)
	return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()


# --------------------------------------------------------------------------
# Metabolism (Figure 4)
def _base_fluxes(path, cell_density_g_per_L):
	"""
	Time-mean flux per base reaction in mmol per g dry weight per hour.

	The listener records mmol per litre of cell per second; dividing by the
	dry-mass density (cell density x dry mass / cell mass) and converting
	seconds to hours puts conditions with different cell composition on the
	same biomass basis.
	"""
	cell = _Cell(path)
	fluxes = cell.column('FBAResults', 'base_reaction_fluxes')
	ids = cell.attribute('FBAResults', 'base_reaction_ids')
	dry = cell.column('Mass', 'dryMass')
	total = cell.column('Mass', 'cellMass')
	if fluxes is None or ids is None or dry is None or total is None:
		return None, None
	with np.errstate(divide='ignore', invalid='ignore'):
		dry_density = cell_density_g_per_L * dry / total
		per_dry = fluxes * SECONDS_PER_HOUR / dry_density[:, None]
	return np.nanmean(per_dry[1:], axis=0), list(ids)


def kcat_enzymes(sim_data, estimator: str = C.KCAT_ESTIMATOR) -> List[str]:
	"""
	Catalysts with a kcat estimate under `estimator`, the enzymes a kcat
	bound can act on. Read from sim_data, or from the flat file for sim_data
	built before kcat estimates existed.
	"""
	metabolism = sim_data.process.metabolism
	estimates = getattr(metabolism, 'kcat_estimates', {}).get(estimator)
	if estimates is None:
		import csv
		from wholecell.utils import filepath
		path = os.path.join(filepath.ROOT_PATH, 'reconstruction', 'ecoli',
			'flat', 'kcat_estimates', 'kcat_estimates_%s.tsv' % estimator)
		with open(path) as f:
			rows = csv.DictReader((line for line in f
				if not line.startswith('#')), delimiter='\t')
			estimates = {(r['reaction_id'], r['catalyst_id']): float(
				r['kcat_estimate']) for r in rows
				if float(r['kcat_estimate']) >= 1.0}
	return sorted({catalyst for _, catalyst in estimates})


def catalyst_monomers(sim_data, catalysts: Sequence[str]) -> List[str]:
	"""Monomer subunits of the given catalysts (a monomer is its own)."""
	complexation = sim_data.process.complexation
	complexes = set(complexation.ids_complexes)
	out: set = set()
	for catalyst in catalysts:
		if catalyst in complexes:
			out.update(str(m) for m in complexation.get_monomers(catalyst)[
				'subunitIds'])
		else:
			out.add(catalyst)
	return sorted(out)


def catalyst_concentrations(data: 'PaperData', variants: Sequence[int],
		plot_out_dir: str, cpus: int = 1) -> Dict[int, VectorProfile]:
	"""
	Concentration (micromolar) of every FBA catalyst, per surviving lineage.

	Each lineage's window-mean catalyst counts over its window-mean volume.
	"""
	counts = vector_profiles(data, 'catalysts', variants, plot_out_dir, cpus)
	volume = data.lineages.set_index(['variant', 'seed'])['cellVolume']
	out = {}
	for variant, profile in counts.items():
		out[variant] = VectorProfile(profile.ids, {seed: values / (
			AVOGADRO * 1e-15 * volume[(variant, seed)] * 1e-6)
			for seed, values in profile.per_lineage.items()})
	return out


def kcat_conditions(batch_dir: str, variants: Sequence[int],
		plot_out_dir: str) -> pd.DataFrame:
	"""
	Whether each variant ran with a kcat bound, and which, read from the
	variant's own modified sim_data rather than its index, so the panels do
	not depend on how a variant file lays out its indices.

	Returns:
		columns variant, bound (bool), estimator, multiplier (NaN when
		unbound), n_bounds. Cached as kcat_conditions.csv.
	"""
	cache = os.path.join(plot_out_dir, CACHE_SUBDIR, 'kcat_conditions.csv')
	if os.path.exists(cache) and not os.environ.get(REBUILD_ENV):
		cached = pd.read_csv(cache)
		if set(variants) <= set(cached['variant']):
			return cached[cached['variant'].isin(variants)].reset_index(
				drop=True)
	rows = []
	for variant in variants:
		paths = glob.glob(os.path.join(batch_dir, '*_%06d' % variant, 'kb',
			'simData_Modified.cPickle'))
		row = dict(variant=variant, bound=False, estimator=None,
			multiplier=np.nan, n_bounds=0)
		if paths:
			with open(paths[0], 'rb') as f:
				metabolism = pickle.load(f).process.metabolism
			selected = getattr(metabolism, 'selected_kcat_estimates', None)
			if selected:
				row.update(bound=True, n_bounds=len(selected),
					estimator=getattr(metabolism, 'kcat_estimate_quantile',
					None), multiplier=float(getattr(metabolism,
					'kcat_estimate_multiplier', np.nan)))
		rows.append(row)
	frame = pd.DataFrame(rows)
	os.makedirs(os.path.dirname(cache), exist_ok=True)
	tmp = cache + '.%d.tmp' % os.getpid()
	frame.to_csv(tmp, index=False)
	os.replace(tmp, cache)
	return frame


def enzyme_concentration_ratios(data: 'PaperData', sim_data,
		variants: Sequence[int], plot_out_dir: str, cpus: int = 1,
		baseline: int = CONTROL_VARIANT) -> Tuple[pd.DataFrame, List[str]]:
	"""
	Concentration of each kcat-bounded enzyme in each variant over the
	baseline, using design means over surviving lineages.

	Returns:
		(DataFrame indexed by catalyst id with one ratio column per variant,
		the catalyst ids used). Enzymes absent from the baseline are dropped.
	"""
	enzymes = set(kcat_enzymes(sim_data))
	profiles = catalyst_concentrations(data, sorted(set(variants)
		| {baseline}), plot_out_dir, cpus)
	ids = profiles[baseline].ids
	keep = [i for i, c in enumerate(ids) if c in enzymes]
	base = profiles[baseline].mean[keep]
	ok = base > 0
	frame = pd.DataFrame(index=pd.Index(np.array(ids)[keep][ok],
		name='catalyst_id'))
	for variant in variants:
		frame[variant] = profiles[variant].mean[keep][ok] / base[ok]
	return frame, list(frame.index)


def monomer_mass_shares(sim_data, profiles: Dict[int, VectorProfile]
		) -> Tuple[List[str], Dict[int, np.ndarray]]:
	"""
	Each monomer's share of total protein mass (count x molecular weight over
	the sum over monomers), from 'monomers' profiles. Counts include
	monomers inside complexes.
	"""
	from wholecell.utils import units
	monomers = sim_data.process.translation.monomer_data
	mw = dict(zip(monomers['id'], monomers['mw'].asNumber(
		units.g / units.mol)))
	ids = next(iter(profiles.values())).ids
	weights = np.array([mw.get(i, np.nan) for i in ids])
	shares = {}
	for variant, profile in profiles.items():
		mass = profile.mean * weights
		shares[variant] = mass / np.nansum(mass)
	return ids, shares


# --------------------------------------------------------------------------
# One call for a panel
@dataclass
class PaperData:
	"""Everything a panel needs, built once per analysis run."""
	gene: GeneInfo
	cycles: pd.DataFrame
	lineages: pd.DataFrame
	designs: pd.DataFrame
	machinery: Optional[MachineryInfo] = None

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


def load_batch(batch_dir: str, cpus: int = 1) -> 'PaperData':
	"""
	PaperData for a batch directory alone, for scripts that combine batches.

	Uses the batch's own metadata and sim_data and caches in
	<batch>/plotOut, the same place its variant analyses cache.
	"""
	with open(os.path.join(batch_dir, 'metadata', 'metadata.json')) as f:
		metadata = json.load(f)
	return load(batch_dir, os.path.join(batch_dir, 'kb', 'simData.cPickle'),
		metadata, os.path.join(batch_dir, 'plotOut'), cpus)


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
	return PaperData(gene, cycles, lineages, designs,
		machinery_info(sim_data))
