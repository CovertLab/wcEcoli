"""
Sub-generational lineage traces for the copy-number feedback arm.

Every analysis in this study so far has averaged over generations, and that is
what has been hiding the effect. Two things were established the hard way:
the quantities of interest move on a 10-26 minute scale inside 50-90 minute
generations, and once a quantity is averaged over a generation it becomes nearly
collinear with cell mass -- which zeroes every elasticity, including the
unburdened control's. So this analysis does not average anything. It emits the
raw per-timestep trace across the whole lineage, in several normalisations, for
looking at.

The hypothesis being tested has two phases:

  Phase 1  GFP expression -> transcription/translation reallocation -> RNAP and
           ribosome loss -> mass accumulation rate loss -> replication
           initiation delayed (critical mass reached more slowly) -> longer
           doubling time AND gene copy number loss.

  Phase 2  copy number loss -> FURTHER RNAP and ribosome loss -> further mass
           accumulation loss -> further lengthening of the doubling time.

Phase 1 is established. Phase 2 is the open arm.

Why the replication rows are the useful ones
--------------------------------------------
Replication initiation in this model is triggered by exactly the mechanism the
hypothesis names, and the intermediate is logged. From
models/ecoli/processes/chromosome_replication.py:

	massPerOrigin = cellMass / n_oric
	criticalMassPerOriC = massPerOrigin / criticalInitiationMass
	if criticalMassPerOriC >= 1.0:   # initiate

and criticalInitiationMass is derived from nutrientToDoublingTime[media_id], so
within one batch it is a CONSTANT -- the bar does not move as the cell slows.
That makes "mass accumulation loss -> delayed initiation" a plotted line
crossing a fixed threshold of 1.0 rather than an inference, and its slope is the
mass accumulation rate per origin. The constancy is checked, not assumed
(gate 2).

What these plots can and cannot settle
--------------------------------------
They cannot establish Phase 2 causally. Within a lineage, copy number, cell
mass, RNAP and ribosome counts all co-vary because a single growth process sets
all of them, so no amount of looking separates "copy loss caused further
machinery loss" from "both fell because growth slowed."

What they can settle is the one comparison that does not have that problem.
GFP accumulates and then plateaus, so view 6 plots everything against GFP
proteome mass fraction instead of time and asks: after the burden stops growing,
does anything keep getting worse? At constant realised burden the direct route
cannot be strengthening, so continued degradation there is not the direct route.
In the phase plane that is a vertical excursion or an open hysteresis loop at the
plateau, against a single retraced curve if there is no second phase.

The rigorous statistic comes later and will be anchored on replication
initiation events, because an initiation steps origin-proximal copy number up at
approximately constant cell mass -- the only discontinuity available, and the
variation the collinearity destroys everywhere else. Hence loop_lineage_events.csv.

Output
------
Everything lands in ONE folder per batch, <batch>/plotOut/loop_lineage_views/,
so it can be zipped and downloaded in one go. This is a variant analysis rather
than a multigen one purely for that reason: analysisMultigen writes to
<variant>/<seed>/plotOut, which would scatter output over 12 directories per
batch and produce 12 fragments of the tidy CSV.
"""

import csv
import gzip
import os
import pickle

import matplotlib
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np

from models.ecoli.analysis import variantAnalysisPlot
from models.ecoli.sim.variants.new_gene_internal_shift import (
	determine_new_gene_ids_and_indices)
from wholecell.analysis.analysis_tools import read_bulk_molecule_counts
from wholecell.io.tablereader import TableReader
from wholecell.utils import units

# Generation at which new gene expression is switched on. Must match
# NEW_GENE_INDUCTION_GEN in models/ecoli/sim/variants/new_gene_internal_shift.py.
INDUCTION_GEN = 8

VARIANTS = (0, 4, 7)

# Seeds read for the DATA. Statistics on 4 seeds cannot reach significance --
# the smallest attainable two-sided Wilcoxon p at n=4 is 0.125 -- so the
# per-generation summary covers every available seed.
SEEDS = tuple(range(16))

# Seeds rendered as FIGURES, and the only ones written at full per-timestep
# resolution. Figures and the full-resolution CSV are what make the output
# large (~200 MB per batch for 4 seeds), and they are for eyeballing, which
# does not need 16 of them. The per-generation summary that carries the
# statistics is tiny and covers all of SEEDS.
FIGURE_SEEDS = (0, 1, 2, 3)

# Terminus-proximal negative control, resolved by replichore position at run
# time. Reused from multigen/copy_number_lineage_trace.py.
TERMINUS_TARGET_FRACTION = 0.95

OUT_SUBDIR = 'loop_lineage_views'
WINDOWS = (('gens06-14', (6, 14)), ('all-gens', None))
SMOOTH_TIMESTEPS = 7

# Figures only -- the CSV keeps every timestep. The simulation timestep is about
# one second, and nothing here changes on that scale, so drawing every point
# costs ~450 MB of vector PDF per batch and buys nothing visible. Thinning to
# one point per PLOT_STRIDE_SEC leaves ~300 plotted points per inch of axis,
# which is finer than the display and far finer than any mid-generation change.
# Stated in the manifest and on every figure rather than applied silently.
PLOT_STRIDE_SEC = 10.0

# Rows, grouped so that each arrow of the hypothesis sits between adjacent
# blocks. Blocks F and G feed back into block B, which is what closes the loop.
# (key, label, prefer_log)
BLOCKS = (
	('A. driver: GFP', (
		('gfp_gene_copies', 'GFP gene\ncopies', False),
		('gfp_init_events', 'GFP transcript\ninits', False),
		('gfp_mrna', 'GFP mRNA\ncount', False),
		('gfp_protein', 'GFP protein\ncount', False),
		('gfp_proteome_frac', 'GFP proteome\nmass fraction', False),
		)),
	('B. mass accumulation', (
		('growth_rate', 'instantaneous\ngrowth rate', False),
		('cell_mass', 'cell mass\n(fg)', False),
		('dry_mass', 'dry mass\n(fg)', False),
		('protein_frac', 'protein mass\nfraction', False),
		('rrna_mass_frac', 'rRNA mass\nfraction', False),
		)),
	('C. replication timing', (
		('mass_per_origin', 'mass per\norigin (fg)', False),
		('critical_mass_per_oric', 'critical mass\nper oriC', False),
		('critical_init_mass', 'critical init\nmass (fg)', False),
		('n_oric', 'origins\nper cell', False),
		('n_forks', 'replication\nforks', False),
		('doubling_time', 'doubling\ntime (min)', False),
		)),
	('D. gene dosage', (
		('rrna_gene_copies', 'rRNA operon\ncopies', False),
		('rnap_gene_copies', 'RNAP subunit\ngene copies', False),
		('terminus_gene_copies', 'terminus control\ngene copies', False),
		)),
	('E. transcription', (
		('rrna_init_events', 'rRNA transcript\ninits', False),
		('rnap_init_events', 'RNAP subunit\ntranscript inits', False),
		('total_init_events', 'total transcript\ninits', False),
		('rnap_budget', 'RNAP budget\n(activations)', False),
		('free_rnap', 'free RNAP\n(APORNAP)', False),
		)),
	('F. machinery and what limits it', (
		('free_16s', 'free 16S\nrRNA', False),
		('free_23s', 'free 23S\nrRNA', False),
		('free_5s', 'free 5S\nrRNA', False),
		('limiting_rprotein', 'limiting r-protein\n(per stoich)', False),
		('limiting_rprotein_idx', 'WHICH r-protein\nis limiting', False),
		('free_30s', 'free 30S', False),
		('free_50s', 'free 50S', False),
		('active_ribosome', 'active\nribosomes', False),
		('total_ribosome', 'total\nribosomes', False),
		('ribosome_elong_rate', 'ribosome elong.\nrate (aa/s)', False),
		('limiting_rnap_subunit', 'limiting RNAP\nsubunit (per stoich)', False),
		('limiting_rnap_idx', 'WHICH RNAP\nsubunit limits', False),
		('active_rnap', 'active\nRNAP', False),
		('total_rnap', 'total\nRNAP', False),
		('rnap_elong_rate', 'RNAP elong.\nrate (nt/s)', False),
		)),
	('G. context', (
		('ppgpp', 'ppGpp\n(uM)', False),
		)),
	)

KEYS = tuple(k for _, rows in BLOCKS for k, _, _ in rows)
LABEL = {k: lab for _, rows in BLOCKS for k, lab, _ in rows}

# Rows that are counts of a thing and so make sense per unit mass (view 5).
# Fractions, rates, ratios and identities do not.
NOT_PER_MASS = frozenset((
	'gfp_proteome_frac', 'growth_rate', 'protein_frac', 'rrna_mass_frac',
	'critical_mass_per_oric', 'critical_init_mass', 'doubling_time',
	'ribosome_elong_rate', 'rnap_elong_rate', 'ppgpp',
	'limiting_rprotein_idx', 'limiting_rnap_idx', 'mass_per_origin',
	'cell_mass', 'dry_mass',
	))

# Rows that must not be renormalised to 1 at generation start (view 4) because
# they are identities or already normalised.
NOT_RELATIVE = frozenset(('limiting_rprotein_idx', 'limiting_rnap_idx'))

VIEWS = (
	('view1-counts', 'raw counts'),
	('view2-derivative', 'd/dt per minute'),
	('view3-smoothed-derivative', 'd/dt, smoothed over %d timesteps'
		% SMOOTH_TIMESTEPS),
	('view4-relative', 'value / value at generation start'),
	('view5-per-mass', 'value / dry mass'),
	('view6-phaseplane', 'x = GFP proteome mass fraction'),
	)


def _smooth(y, w):
	"""Centred moving average that keeps length and tolerates NaN."""
	if w < 2 or y.size < w:
		return y
	k = np.ones(w) / w
	good = np.isfinite(y).astype(float)
	filled = np.where(np.isfinite(y), y, 0.0)
	num = np.convolve(filled, k, mode='same')
	den = np.convolve(good, k, mode='same')
	with np.errstate(invalid='ignore', divide='ignore'):
		return np.where(den > 0, num / den, np.nan)


class Plot(variantAnalysisPlot.VariantAnalysisPlot):
	def do_plot(self, inputDir, plotOutDir, plotOutFileName, simDataFile,
			validationDataFile, metadata):
		with open(simDataFile, 'rb') as handle:
			sim_data = pickle.load(handle)

		out_dir = os.path.join(plotOutDir, OUT_SUBDIR)
		if not os.path.exists(out_dir):
			os.makedirs(out_dir)

		ctx = self._context(sim_data)
		self._notes = []
		self._gates = {}
		self._bn_rows = []

		available = sorted(set(self.ap.get_variants()) & set(VARIANTS))
		if not available:
			print('None of variants %s are present in this batch.' % (VARIANTS,))
			return
		self._note('variants plotted: %s' % (available,))

		# traces[(variant, seed)] = dict of key -> array, plus 'time',
		# 'gen_bounds', 'init_events'
		traces, events = {}, []
		for variant in available:
			for seed in SEEDS:
				tr = self._lineage(variant, seed, ctx)
				if tr is None:
					continue
				traces[(variant, seed)] = tr
				events.extend(tr.pop('_events'))

		if not traces:
			print('No readable lineages.')
			return

		self._check_gates(traces)
		self._write_by_gen(out_dir, traces)
		self._write_bottleneck(out_dir)
		self._write_timeseries(out_dir, traces)
		self._write_events(out_dir, events)
		self._render(out_dir, traces, available)
		self._write_manifest(out_dir, ctx)
		print('\nWrote %s' % out_dir)

	# ---- context ---------------------------------------------------------

	def _note(self, msg):
		print(msg)
		self._notes.append(msg)

	def _context(self, sim_data):
		"""Resolve every id and index once, against sim_data."""
		transcription = sim_data.process.transcription
		rna_data = transcription.rna_data
		rna_ids = list(rna_data['id'])

		(new_mrna_ids, _, new_monomer_ids, _
			) = determine_new_gene_ids_and_indices(sim_data)

		# TU index of the construct: the TU whose cistron set is the new gene.
		is_new_tu = np.zeros(len(rna_ids), bool)
		cistron_tu = transcription.cistron_tu_mapping_matrix
		cistron_ids = list(transcription.cistron_data['id'])
		for mrna_id in new_mrna_ids:
			ci = cistron_ids.index(mrna_id)
			is_new_tu |= cistron_tu[ci, :].toarray().ravel() > 0

		# Replichore fraction, normalised PER ARM as in
		# new_gene_initiation_chain.py:184 and copy_number_lineage_trace.py:70.
		# Dividing |coordinate| by the sum of both arms instead caps the
		# fraction at ~0.5 and no gene ever reads as terminus-proximal.
		coords = rna_data['replication_coordinate'].astype(float)
		lengths = sim_data.process.replication.replichore_lengths
		frac = np.where(coords > 0, coords / lengths[0], -coords / lengths[1])

		is_rrna = rna_data['is_rRNA'].copy()
		is_rnap = rna_data['includes_RNAP'].copy()

		# Terminus control: the non-rRNA, non-RNAP, non-construct TU closest to
		# TERMINUS_TARGET_FRACTION.
		eligible = ~(is_rrna | is_rnap | is_new_tu)
		cand = np.where(eligible)[0]
		term_idx = cand[np.argmin(
			np.abs(frac[cand] - TERMINUS_TARGET_FRACTION))]

		mg, mi = sim_data.molecule_groups, sim_data.molecule_ids
		complexation = sim_data.process.complexation

		def stoich_of(complex_id, subunit_ids):
			mons = complexation.get_monomers(complex_id)
			lookup = dict(zip(mons['subunitIds'], mons['subunitStoich']))
			return np.array([lookup[s] for s in subunit_ids], float)

		r_prot_ids = list(mg.s30_proteins) + list(mg.s50_proteins)
		r_prot_stoich = np.concatenate((
			stoich_of(mi.s30_full_complex, mg.s30_proteins),
			stoich_of(mi.s50_full_complex, mg.s50_proteins)))

		rnap_sub_ids = list(mg.RNAP_subunits)
		rnap_sub_stoich = stoich_of(mi.full_RNAP, rnap_sub_ids)

		# Per-gene bottleneck set: every ribosomal-protein and RNAP subunit,
		# resolved to the TU carrying its cistron so its promoter copy number
		# can be read alongside its free monomer count. Ribosome and RNAP
		# assembly are bounded by min(free subunit / stoichiometry), so these
		# are the genes through which a dosage effect could actually propagate
		# -- unlike rRNA, which runs in large surplus.
		monomer_data = sim_data.process.translation.monomer_data.struct_array
		mono_to_cistron = dict(zip(monomer_data['id'],
			monomer_data['cistron_id']))
		cistron_index = {c: i for i, c in enumerate(cistron_ids)}
		bottleneck = []
		for group, ids, stoichs in (
				('s30_protein', mg.s30_proteins,
					stoich_of(mi.s30_full_complex, mg.s30_proteins)),
				('s50_protein', mg.s50_proteins,
					stoich_of(mi.s50_full_complex, mg.s50_proteins)),
				('rnap_subunit', rnap_sub_ids, rnap_sub_stoich)):
			for mid, st in zip(list(ids), stoichs):
				cis = mono_to_cistron.get(mid)
				ci = cistron_index.get(cis)
				if ci is None:
					continue
				tus = np.where(cistron_tu[ci, :].toarray().ravel() > 0)[0]
				if not tus.size:
					continue
				bottleneck.append(dict(group=group, monomer=mid, cistron=cis,
					stoich=float(st), tus=tus,
					frac=float(frac[tus].mean())))
		print('Bottleneck set: %d subunit genes resolved to TUs'
			% len(bottleneck))

		# GFP monomer mass, for the proteome fraction.
		gfp_mass = float(sum(
			sim_data.getter.get_mass(m).asNumber(units.fg / units.count)
			for m in new_monomer_ids))

		ctx = dict(
			rna_ids=rna_ids,
			new_tu=np.where(is_new_tu)[0],
			new_mrna_ids=new_mrna_ids,
			new_monomer_ids=new_monomer_ids,
			gfp_monomer_mass_fg=gfp_mass,
			rrna_tu=np.where(is_rrna)[0],
			rnap_tu=np.where(is_rnap)[0],
			term_tu=np.array([term_idx]),
			term_frac=float(frac[term_idx]),
			term_id=rna_ids[term_idx],
			r_prot_ids=r_prot_ids, r_prot_stoich=r_prot_stoich,
			rnap_sub_ids=rnap_sub_ids, rnap_sub_stoich=rnap_sub_stoich,
			s16=list(mg.s30_16s_rRNA), s23=list(mg.s50_23s_rRNA),
			s5=list(mg.s50_5s_rRNA),
			s30=[mi.s30_full_complex], s50=[mi.s50_full_complex],
			free_rnap=[mi.full_RNAP],
			bottleneck=bottleneck,
			)
		print('Construct TUs %s; rRNA TUs %d; RNAP-subunit TUs %d; '
			'terminus control %s at f=%.3f'
			% (ctx['new_tu'], ctx['rrna_tu'].size, ctx['rnap_tu'].size,
				ctx['term_id'], ctx['term_frac']))
		return ctx

	# ---- reading ---------------------------------------------------------

	def _lineage(self, variant, seed, ctx):
		"""Concatenate every generation of one lineage into continuous arrays."""
		by_gen = {}
		for path in self.ap.get_cells(variant=[variant], seed=[seed]):
			by_gen[int(self.ap.get_cell_generation(path))] = path
		if not by_gen:
			return None

		acc = {k: [] for k in KEYS}
		times, bounds, ev = [], [], []
		n_stale = 0
		for gen in sorted(by_gen):
			one = self._one_cell(by_gen[gen], ctx)
			if one is None:
				continue
			t = one.pop('_time')
			n_stale += one.pop('_stale')
			for row in one.pop('_bottleneck', []):
				self._bn_rows.append(dict(variant=variant, seed=seed,
					generation=gen, **row))
			# Doubling time is a per-generation scalar, drawn as a step.
			one['doubling_time'] = np.full(
				t.size, (t[-1] - t[0]) / 60.0 if t.size > 1 else np.nan)
			times.append(t)
			bounds.append((float(t[0]), float(t[-1]), gen))
			for k in KEYS:
				acc[k].append(one.get(k, np.full(t.size, np.nan)))
			# Replication initiations: numberOfOric steps up.
			no = one['n_oric']
			rise = np.where(np.diff(no) > 0)[0] + 1
			for i in rise:
				ev.append(dict(variant=variant, seed=seed, generation=gen,
					kind='replication_initiation', time_min=t[i] / 60.0,
					critical_mass_per_oric=one['critical_mass_per_oric'][i],
					mass_per_origin=one['mass_per_origin'][i],
					n_oric=no[i],
					rrna_gene_copies=one['rrna_gene_copies'][i],
					rnap_gene_copies=one['rnap_gene_copies'][i],
					terminus_gene_copies=one['terminus_gene_copies'][i]))
			ev.append(dict(variant=variant, seed=seed, generation=gen,
				kind='division', time_min=t[-1] / 60.0,
				critical_mass_per_oric=one['critical_mass_per_oric'][-1],
				mass_per_origin=one['mass_per_origin'][-1],
				n_oric=no[-1],
				rrna_gene_copies=one['rrna_gene_copies'][-1],
				rnap_gene_copies=one['rnap_gene_copies'][-1],
				terminus_gene_copies=one['terminus_gene_copies'][-1]))

		if not times:
			return None
		tr = {k: np.concatenate(acc[k]) for k in KEYS}
		tr['time'] = np.concatenate(times) / 60.0
		tr['gen_bounds'] = [(a / 60.0, b / 60.0, g) for a, b, g in bounds]
		tr['_events'] = ev
		if n_stale:
			self._note('v%d seed %d: %d timesteps (about one per generation) '
				'where the replication listener had not yet assigned '
				'criticalInitiationMass; both replication-threshold rows are '
				'blanked there' % (variant, seed, n_stale))
		self._note('v%d seed %d: %d generations, %d timesteps'
			% (variant, seed, len(bounds), tr['time'].size))
		return tr

	def _one_cell(self, path, ctx):
		sim_out = os.path.join(path, 'simOut')
		try:
			t = TableReader(os.path.join(sim_out, 'Main')).readColumn('time')
			mass = TableReader(os.path.join(sim_out, 'Mass'))
			rsp = TableReader(os.path.join(sim_out, 'RnaSynthProb'))
			rnap = TableReader(os.path.join(sim_out, 'RnapData'))
			repl = TableReader(os.path.join(sim_out, 'ReplicationData'))
			ribo = TableReader(os.path.join(sim_out, 'RibosomeData'))
			rnac = TableReader(os.path.join(sim_out, 'RNACounts'))
			mono = TableReader(os.path.join(sim_out, 'MonomerCounts'))
			umc = TableReader(os.path.join(sim_out, 'UniqueMoleculeCounts'))
			gl = TableReader(os.path.join(sim_out, 'GrowthLimits'))
		except Exception as exc:  # noqa: BLE001 - skip unreadable generations
			print('  unreadable: %s (%s)' % (path, exc))
			return None

		t = np.atleast_1d(t).astype(float)
		n = t.size
		out = {'_time': t}

		copies = rsp.readColumn('promoter_copy_number')
		inits = rnap.readColumn('rnaInitEvent')

		take = lambda arr, idx, how: (
			how(arr[:, idx], axis=1) if idx.size else np.full(n, np.nan))

		out['gfp_gene_copies'] = take(copies, ctx['new_tu'], np.sum)
		out['gfp_init_events'] = take(inits, ctx['new_tu'], np.sum)
		out['rrna_gene_copies'] = take(copies, ctx['rrna_tu'], np.mean)
		out['rnap_gene_copies'] = take(copies, ctx['rnap_tu'], np.mean)
		out['terminus_gene_copies'] = take(copies, ctx['term_tu'], np.mean)
		out['rrna_init_events'] = take(inits, ctx['rrna_tu'], np.sum)
		out['rnap_init_events'] = take(inits, ctx['rnap_tu'], np.sum)
		out['total_init_events'] = np.sum(inits, axis=1)
		out['rnap_budget'] = rnap.readColumn('didInitialize').astype(float)

		# determine_new_gene_ids_and_indices returns CISTRON ids, so this must
		# read mRNA_cistron_counts. mRNA_counts is indexed by mRNA_ids, which
		# are TU ids ('TU-8381[c]'), and the lookup silently matches nothing.
		mrna_ids = list(rnac.readAttribute('mRNA_cistron_ids'))
		mrna = rnac.readColumn('mRNA_cistron_counts')
		gi = [mrna_ids.index(i) for i in ctx['new_mrna_ids'] if i in mrna_ids]
		out['gfp_mrna'] = (mrna[:, gi].sum(axis=1) if gi
			else np.full(n, np.nan))

		mono_ids = list(mono.readAttribute('monomerIds'))
		counts = mono.readColumn('monomerCounts')
		mj = [mono_ids.index(i) for i in ctx['new_monomer_ids']
			if i in mono_ids]
		out['gfp_protein'] = (counts[:, mj].sum(axis=1) if mj
			else np.full(n, np.nan))

		protein_mass = mass.readColumn('proteinMass')
		dry = mass.readColumn('dryMass')
		cell = mass.readColumn('cellMass')
		with np.errstate(invalid='ignore', divide='ignore'):
			out['gfp_proteome_frac'] = (
				out['gfp_protein'] * ctx['gfp_monomer_mass_fg'] / protein_mass)
			out['protein_frac'] = protein_mass / dry
			out['rrna_mass_frac'] = mass.readColumn('rRnaMass') / dry
		out['growth_rate'] = mass.readColumn('instantaneous_growth_rate')
		out['cell_mass'] = cell
		out['dry_mass'] = dry

		n_oric = repl.readColumn('numberOfOric').astype(float)
		out['n_oric'] = n_oric
		out['critical_mass_per_oric'] = repl.readColumn('criticalMassPerOriC')
		out['critical_init_mass'] = repl.readColumn('criticalInitiationMass')
		with np.errstate(invalid='ignore', divide='ignore'):
			out['mass_per_origin'] = np.where(n_oric > 0, cell / n_oric, np.nan)
		# criticalMassPerOriC and criticalInitiationMass are initialised to 0. in
		# the listener's allocate() and only assigned inside calculateRequest.
		# On the first timestep of each generation the listener therefore
		# reports 0 for both while numberOfOric is already 4 -- so the stale
		# condition is the VALUE being non-positive, not n_oric == 0 (which
		# never happens here). Blank those out: left in, the one zero per
		# generation makes criticalInitiationMass look variable when the whole
		# point of that row is that it is constant. The count is reported in the
		# manifest, so this is not a silent drop.
		stale = ~(out['critical_init_mass'] > 0)
		out['_stale'] = int(np.count_nonzero(stale))
		out['critical_mass_per_oric'] = np.where(
			stale, np.nan, out['critical_mass_per_oric'])
		out['critical_init_mass'] = np.where(
			stale, np.nan, out['critical_init_mass'])
		forks = repl.readColumn('fork_coordinates')
		out['n_forks'] = (np.isfinite(np.atleast_2d(forks)).sum(axis=1) / 2.0
			if forks.ndim > 1 else np.full(n, np.nan))

		# Bulk molecules: one call for everything, as the helper requires.
		(free_rnap, s30, s50, s16, s23, s5, r_prot, rnap_sub
			) = read_bulk_molecule_counts(sim_out, (
				ctx['free_rnap'], ctx['s30'], ctx['s50'], ctx['s16'],
				ctx['s23'], ctx['s5'], ctx['r_prot_ids'], ctx['rnap_sub_ids']))
		out['free_rnap'] = np.ravel(free_rnap).astype(float)
		out['free_30s'] = np.ravel(s30).astype(float)
		out['free_50s'] = np.ravel(s50).astype(float)
		out['free_16s'] = np.atleast_2d(s16).sum(axis=1).astype(float)
		out['free_23s'] = np.atleast_2d(s23).sum(axis=1).astype(float)
		out['free_5s'] = np.atleast_2d(s5).sum(axis=1).astype(float)

		# Limiting subunit: fewest copies once divided by stoichiometry. The
		# monomer counts include subunits already inside assembled complexes,
		# so this is a supply measure, not a free-pool measure.
		per = np.atleast_2d(r_prot).astype(float) / ctx['r_prot_stoich']
		out['limiting_rprotein'] = per.min(axis=1)
		out['limiting_rprotein_idx'] = per.argmin(axis=1).astype(float)
		per_rnap = (np.atleast_2d(rnap_sub).astype(float)
			/ ctx['rnap_sub_stoich'])
		out['limiting_rnap_subunit'] = per_rnap.min(axis=1)
		out['limiting_rnap_idx'] = per_rnap.argmin(axis=1).astype(float)

		uids = list(umc.readAttribute('uniqueMoleculeIds'))
		ucounts = umc.readColumn('uniqueMoleculeCounts')
		out['active_ribosome'] = ucounts[:, uids.index('active_ribosome')
			].astype(float)
		out['active_rnap'] = ucounts[:, uids.index('active_RNAP')
			].astype(float) if 'active_RNAP' in uids else np.full(n, np.nan)
		out['total_ribosome'] = out['active_ribosome'] + np.minimum(
			out['free_30s'], out['free_50s'])
		out['total_rnap'] = out['active_rnap'] + out['free_rnap']

		out['ribosome_elong_rate'] = ribo.readColumn('effectiveElongationRate')
		dt = np.diff(t, prepend=t[0] - (t[1] - t[0] if n > 1 else 1.0))
		with np.errstate(invalid='ignore', divide='ignore'):
			out['rnap_elong_rate'] = np.where(out['active_rnap'] > 0,
				rnap.readColumn('actualElongations') / dt / out['active_rnap'],
				np.nan)
		out['ppgpp'] = gl.readColumn('ppgpp_conc')

		# Per-gene bottleneck traces. Free monomer counts come from
		# BulkMolecules, so these are FREE pools -- subunits already inside an
		# assembled ribosome or RNAP are excluded, which is what makes
		# free/stoich an assembly bound rather than a total-supply figure.
		bn = ctx['bottleneck']
		if bn:
			(free_sub,) = read_bulk_molecule_counts(sim_out,
				([b['monomer'] for b in bn],))
			free_sub = np.atleast_2d(free_sub).astype(float)
			bn_rows = []
			for j, b in enumerate(bn):
				bn_rows.append(dict(group=b['group'], monomer=b['monomer'],
					stoich=b['stoich'], frac=b['frac'],
					# Summed over every TU variant carrying this cistron, so
					# this is the total PROMOTER copy number driving the gene,
					# not its physical chromosomal dosage. rplB sits in three
					# TU variants, so the two differ by that multiplicity.
					# Promoter copies is the right quantity for transcription;
					# do not read it as gene dosage. n_tus relates the two.
					promoter_copies=float(
						copies[:, b['tus']].sum(axis=1).mean()),
					n_tus=int(b['tus'].size),
					free_count=float(free_sub[:, j].mean()),
					assembly_bound=float(
						(free_sub[:, j] / b['stoich']).mean())))
			out['_bottleneck'] = bn_rows

		for k, v in list(out.items()):
			if k.startswith('_'):
				continue
			v = np.asarray(v, float).ravel()
			if v.size != n:
				out[k] = np.full(n, np.nan)
			else:
				out[k] = v
		return out

	# ---- gates -----------------------------------------------------------

	def _check_gates(self, traces):
		"""The verification gates, recorded so a failure is visible in the
		download rather than only in the console."""
		# Gate: criticalMassPerOriC crosses 1.0 where numberOfOric rises.
		hit = tot = 0
		for tr in traces.values():
			rise = np.where(np.diff(tr['n_oric']) > 0)[0] + 1
			tot += rise.size
			for i in rise:
				lo = max(0, i - 2)
				if np.nanmax(tr['critical_mass_per_oric'][lo:i + 1]) >= 1.0:
					hit += 1
		self._gates['initiation_threshold'] = (
			'%d of %d origin increases had criticalMassPerOriC >= 1.0 in the '
			'preceding 2 timesteps%s' % (hit, tot,
				'' if tot and hit == tot else '   *** REVIEW'))

		# Gate: the threshold is constant within the batch.
		#
		# criticalInitiationMass is initialised to 0. in the listener's
		# allocate() and only assigned inside calculateRequest, so on the first
		# timestep of each generation it reads 0 while numberOfOric is already
		# 4. Those zeros must be excluded or the range reads as the whole mean
		# and this gate false-fails. _one_cell blanks them to NaN; this filter
		# is belt-and-braces so the gate cannot be fooled by either form.
		raw = np.concatenate([tr['critical_init_mass'] for tr in
			traces.values()])
		n_masked = int(np.count_nonzero(~np.isfinite(raw)))
		allv = raw[np.isfinite(raw)]
		n_masked += int(np.count_nonzero(allv <= 0))
		allv = allv[allv > 0]
		spread = float(np.ptp(allv)) if allv.size else float('nan')
		self._gates['critical_mass_constant'] = (
			'criticalInitiationMass range %.6g over %d positive timesteps '
			'(mean %.6g; %d unassigned timesteps excluded)%s'
			% (spread, allv.size,
				float(np.mean(allv)) if allv.size else float('nan'), n_masked,
				'' if spread < 1e-6 else '   *** NOT CONSTANT, section C '
				'readings assume it is'))

		# Gate: rows that came back entirely NaN.
		dead = [k for k in KEYS if all(
			not np.any(np.isfinite(tr[k])) for tr in traces.values())]
		self._gates['all_nan_rows'] = (
			'none' if not dead else '%s   *** these rows carry no data'
			% (dead,))

		for name, msg in self._gates.items():
			self._note('GATE %s: %s' % (name, msg))

	# ---- views -----------------------------------------------------------

	def _stride(self, tr):
		"""Timestep stride giving roughly PLOT_STRIDE_SEC between plotted
		points. Figures only; never applied to the CSVs."""
		t = tr['time']
		if t.size < 3:
			return 1
		dt_min = float(np.median(np.diff(t)))
		if not np.isfinite(dt_min) or dt_min <= 0:
			return 1
		return max(1, int(round(PLOT_STRIDE_SEC / (dt_min * 60.0))))

	def _transform(self, tr, view):
		"""Return (x, {key: y}) for one view of one lineage."""
		t = tr['time']
		if view == 'view1-counts':
			return t, {k: tr[k] for k in KEYS}
		if view in ('view2-derivative', 'view3-smoothed-derivative'):
			dt = np.gradient(t)
			out = {}
			for k in KEYS:
				with np.errstate(invalid='ignore', divide='ignore'):
					d = np.gradient(tr[k]) / dt
				out[k] = (_smooth(d, SMOOTH_TIMESTEPS)
					if view.startswith('view3') else d)
			return t, out
		if view == 'view4-relative':
			out = {}
			for k in KEYS:
				if k in NOT_RELATIVE:
					out[k] = tr[k]
					continue
				y = np.full(t.size, np.nan)
				for a, b, _ in tr['gen_bounds']:
					m = (t >= a) & (t <= b)
					if not np.any(m):
						continue
					seg = tr[k][m]
					ref = next((v for v in seg if np.isfinite(v) and v != 0),
						np.nan)
					with np.errstate(invalid='ignore', divide='ignore'):
						y[m] = seg / ref
				out[k] = y
			return t, out
		if view == 'view5-per-mass':
			out = {}
			for k in KEYS:
				if k in NOT_PER_MASS:
					out[k] = tr[k]
				else:
					with np.errstate(invalid='ignore', divide='ignore'):
						out[k] = tr[k] / tr['dry_mass']
			return t, out
		if view == 'view6-phaseplane':
			return tr['gfp_proteome_frac'], {k: tr[k] for k in KEYS}
		raise ValueError(view)

	def _render(self, out_dir, traces, variants):
		matplotlib.rcParams.update({'font.size': 8,
			'axes.spines.top': False, 'axes.spines.right': False})
		for variant in variants:
			seeds = [s for s in FIGURE_SEEDS if (variant, s) in traces]
			if not seeds:
				continue
			for view, view_desc in VIEWS:
				for wname, wrange in WINDOWS:
					self._one_pdf(out_dir, variant, seeds, traces, view,
						view_desc, wname, wrange)

	def _one_pdf(self, out_dir, variant, seeds, traces, view, view_desc,
			wname, wrange):
		path = os.path.join(out_dir, 'v%02d_%s_%s.pdf' % (variant, view, wname))
		n_rows = len(KEYS)
		with PdfPages(path) as pdf:
			for seed in seeds:
				tr = traces[(variant, seed)]
				x, ys = self._transform(tr, view)
				stride = self._stride(tr)
				keep = np.zeros(x.size, bool)
				keep[::stride] = True
				if wrange is not None:
					gens = np.full(tr['time'].size, -1)
					for a, b, g in tr['gen_bounds']:
						gens[(tr['time'] >= a) & (tr['time'] <= b)] = g
					keep &= (gens >= wrange[0]) & (gens <= wrange[1])
					if not np.any(keep):
						continue
				fig, axes = plt.subplots(n_rows, 1, sharex=True,
					figsize=(18, 1.4 * n_rows))
				phase = view == 'view6-phaseplane'
				for ax, key in zip(axes, KEYS):
					xi, yi = x[keep], ys[key][keep]
					if phase:
						ax.plot(xi, yi, '.', ms=1.6, color='#4a3aa7',
							alpha=.55)
					else:
						ax.plot(xi, yi, lw=.9, color='#4a3aa7')
					ax.set_ylabel(LABEL[key], fontsize=7.5, rotation=0,
						ha='right', va='center', labelpad=8)
					ax.tick_params(labelsize=7)
					ax.margins(x=.01)
					if key == 'critical_mass_per_oric':
						ax.axhline(1.0, color='#b3123c', lw=1.0, ls='-')
					if key in ('gfp_gene_copies', 'gfp_proteome_frac'):
						ax.axhline(0, color='#e2e1dc', lw=.8)
					if view == 'view4-relative' and key not in NOT_RELATIVE:
						ax.axhline(1.0, color='#767469', lw=.7, ls=':')
						ax.axhline(2.0, color='#0d8f63', lw=.9, ls='--')
					if not phase:
						self._mark(ax, tr, wrange)
				axes[-1].set_xlabel('GFP proteome mass fraction' if phase
					else 'time (min, continuous across generations)')
				fig.suptitle('variant %d, seed %d  |  %s  |  %s  |  %s\n'
					'dashed grey = division, dotted orange = replication '
					'initiation%s  |  plotted every %d timesteps '
					'(~%.0f s); the CSV keeps every timestep'
					% (variant, seed, view, view_desc, wname,
						'  (not marked in the phase plane)' if phase else '',
						stride, PLOT_STRIDE_SEC),
					fontsize=10, y=.999)
				fig.tight_layout(rect=(0, 0, 1, .992))
				pdf.savefig(fig)
				if view == 'view1-counts' and wname == 'gens06-14':
					fig.savefig(os.path.join(out_dir,
						'v%02d_seed%02d_counts_gens06-14.png' % (variant, seed)),
						dpi=110, bbox_inches='tight')
				plt.close(fig)

	def _mark(self, ax, tr, wrange):
		"""Division and replication-initiation markers."""
		for a, b, g in tr['gen_bounds']:
			if wrange is not None and not (wrange[0] <= g <= wrange[1]):
				continue
			ax.axvline(b, color='#767469', lw=.7, ls='--', alpha=.8)
		t = tr['time']
		rise = np.where(np.diff(tr['n_oric']) > 0)[0] + 1
		for i in rise:
			g = next((gg for a, b, gg in tr['gen_bounds']
				if a <= t[i] <= b), None)
			if wrange is not None and (g is None
					or not (wrange[0] <= g <= wrange[1])):
				continue
			ax.axvline(t[i], color='#c2410c', lw=.7, ls=':', alpha=.85)

	# ---- writing ---------------------------------------------------------

	def _write_by_gen(self, out_dir, traces):
		"""Per-(variant, seed, generation) summary for EVERY seed.

		This is what the statistics run on. Four seeds cannot support an
		inferential claim -- the smallest attainable two-sided Wilcoxon p at
		n=4 is 0.125, so a real effect and a null are indistinguishable -- and
		this table is small enough (a few thousand rows) to cover all 16
		without the download cost of full per-timestep traces.

		Three aggregations per generation, because copy number oscillates about
		two-fold within a cycle and a mean over the cycle hides where in it a
		quantity sits: the generation mean, the value at birth (first 5% of the
		cycle) and the value at the end (last 5%).
		"""
		path = os.path.join(out_dir, 'loop_lineage_by_gen.csv')
		with open(path, 'w', newline='') as fh:
			w = csv.writer(fh)
			w.writerow(['variant', 'seed', 'generation', 'aggregation',
				'n_timesteps', 'duration_min'] + list(KEYS))
			for (variant, seed), tr in sorted(traces.items()):
				t = tr['time']
				for a, b, g in tr['gen_bounds']:
					m = (t >= a) & (t <= b)
					n = int(np.count_nonzero(m))
					if not n:
						continue
					e = max(1, n // 20)
					idx = np.where(m)[0]
					for label, sel in (('mean', idx),
							('birth', idx[:e]), ('end', idx[-e:])):
						vals = []
						for k in KEYS:
							y = tr[k][sel]
							y = y[np.isfinite(y)]
							vals.append('' if not y.size
								else '%.6g' % float(np.mean(y)))
						w.writerow([variant, seed, g, label, n,
							'%.4f' % (b - a)] + vals)
		print('  %s (%.2f MB)'
			% (os.path.basename(path), os.path.getsize(path) / 1e6))

	def _write_bottleneck(self, out_dir):
		"""Per-gene traces for every ribosomal-protein and RNAP subunit.

		The point of this table: ribosome and RNAP assembly are bounded by
		min(free subunit / stoichiometry), so if a gene-dosage effect
		propagates to the machinery at all, it propagates through these genes.
		It cannot propagate through rRNA, which runs at 100-200x surplus under
		burden. Emitting gene copy number beside the free pool for each subunit
		is what makes the per-gene mediation testable.
		"""
		if not self._bn_rows:
			return
		path = os.path.join(out_dir, 'loop_lineage_bottleneck.csv')
		cols = ['variant', 'seed', 'generation', 'group', 'monomer', 'stoich',
			'frac', 'n_tus', 'promoter_copies', 'free_count',
			'assembly_bound']
		with open(path, 'w', newline='') as fh:
			w = csv.DictWriter(fh, fieldnames=cols)
			w.writeheader()
			for r in self._bn_rows:
				w.writerow({c: ('%.6g' % r[c] if isinstance(r[c], float)
					else r[c]) for c in cols})
		print('  %s (%d rows, %.2f MB)' % (os.path.basename(path),
			len(self._bn_rows), os.path.getsize(path) / 1e6))

	def _write_timeseries(self, out_dir, traces):
		"""One row per timestep, one column per quantity.

		Wide rather than long: long format repeats variant/seed/generation/time
		and the quantity name for all 40 values of every timestep, which came to
		21 MB gzipped for a single variant and two seeds -- about 130 MB for a
		full batch, which defeats the point of a downloadable folder. Wide is
		the same information at roughly a quarter the size and loads directly
		into pandas.
		"""
		path = os.path.join(out_dir, 'loop_lineage_timeseries.csv.gz')
		with gzip.open(path, 'wt', newline='') as fh:
			w = csv.writer(fh)
			w.writerow(['variant', 'seed', 'generation', 'time_min']
				+ list(KEYS))
			for (variant, seed), tr in sorted(traces.items()):
				if seed not in FIGURE_SEEDS:
					continue
				t = tr['time']
				gens = np.full(t.size, -1)
				for a, b, g in tr['gen_bounds']:
					gens[(t >= a) & (t <= b)] = g
				cols = [tr[k] for k in KEYS]
				for i in range(t.size):
					w.writerow([variant, seed, gens[i], '%.4f' % t[i]]
						+ ['' if not np.isfinite(c[i]) else '%.6g' % c[i]
							for c in cols])
		print('  %s (%.1f MB)'
			% (os.path.basename(path), os.path.getsize(path) / 1e6))

	def _write_events(self, out_dir, events):
		path = os.path.join(out_dir, 'loop_lineage_events.csv')
		cols = ['variant', 'seed', 'generation', 'kind', 'time_min',
			'critical_mass_per_oric', 'mass_per_origin', 'n_oric',
			'rrna_gene_copies', 'rnap_gene_copies', 'terminus_gene_copies']
		with open(path, 'w', newline='') as fh:
			w = csv.DictWriter(fh, fieldnames=cols)
			w.writeheader()
			for e in events:
				w.writerow({c: e.get(c, '') for c in cols})
		print('  %s (%d events)' % (os.path.basename(path), len(events)))

	def _write_manifest(self, out_dir, ctx):
		path = os.path.join(out_dir, 'loop_lineage_manifest.txt')
		with open(path, 'w') as fh:
			fh.write('new_gene_loop_lineage_views\n')
			fh.write('=' * 60 + '\n\nVERIFICATION GATES\n')
			for name, msg in self._gates.items():
				fh.write('  %-26s %s\n' % (name, msg))
			fh.write('\nNOTES\n')
			for msg in self._notes:
				fh.write('  %s\n' % msg)
			fh.write('\nTERMINUS CONTROL\n  %s at replichore fraction %.3f\n'
				% (ctx['term_id'], ctx['term_frac']))
			fh.write('\nSEED COVERAGE\n  loop_lineage_by_gen.csv: seeds %s '
				'(statistics)\n  loop_lineage_timeseries.csv.gz and all '
				'figures: seeds %s only\n'
				% (list(SEEDS), list(FIGURE_SEEDS)))
			fh.write('\nFIGURE THINNING\n  figures plot one point per ~%.0f s '
				'(PLOT_STRIDE_SEC); loop_lineage_timeseries.csv.gz keeps every '
				'timestep at full resolution\n' % PLOT_STRIDE_SEC)
			fh.write('\nLIMITING-SUBUNIT INDEX LEGEND\n')
			fh.write('  limiting_rprotein_idx -> index into:\n')
			for i, mid in enumerate(ctx['r_prot_ids']):
				fh.write('    %3d  %s (stoich %g)\n'
					% (i, mid, ctx['r_prot_stoich'][i]))
			fh.write('  limiting_rnap_idx -> index into:\n')
			for i, mid in enumerate(ctx['rnap_sub_ids']):
				fh.write('    %3d  %s (stoich %g)\n'
					% (i, mid, ctx['rnap_sub_stoich'][i]))
			fh.write('\nROWS, IN PLOT ORDER\n')
			for block, rows in BLOCKS:
				fh.write('  %s\n' % block)
				for k, lab, _ in rows:
					fh.write('    %-24s %s\n'
						% (k, lab.replace('\n', ' ')))
		print('  %s' % os.path.basename(path))


if __name__ == '__main__':
	Plot().cli()
