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

# Seeds written at full PER-TIMESTEP resolution. The event-anchored test
# (report 42) needs this resolution, and its events cluster within a lineage --
# so 4 seeds would put it back at the n=4 power ceiling that the 16-seed run
# exists to remove. Kept equal to SEEDS by default for that reason; the cost is
# CSV size, roughly 200 MB per batch, not runtime.
TRACE_SEEDS = SEEDS

# Seeds rendered as FIGURES: 36 multi-page PDFs, one page per seed. This is
# the expensive output -- at 16 seeds it is 576 pages and roughly 850 MB per
# batch, against ~210 MB at 4 seeds. Kept equal to SEEDS anyway: the
# lineage-to-lineage spread is itself part of what these figures are for, and
# judging which lineages are typical from a quarter of them is how the n=4
# reading went wrong before. Cut this back to a subset if the download size
# becomes the binding constraint.
FIGURE_SEEDS = SEEDS

# Named subunits that actually hold the assembly bottleneck, measured in
# loop_lineage_views_minimal_v4. Group sums average these away: RNAP is held by
# rpoC 71-87% of the time, and 2-4 r-proteins carry 80% of the protein-limited
# time once rRNA clears. Each gets its own dosage / transcription / translation
# / free-pool quartet so every step of the chain can be followed on ONE gene.
#
# rpoB and rpoC share an operon and a promoter position (f=0.109), so dosage
# cannot separate them; the difference is production cost (rpoB efficiency
# 0.830 vs rpoC 0.810, length 1342 vs 1407 aa). Both are logged for that
# contrast. Labels are common names; ids are resolved at run time and a
# missing id yields NaN columns rather than an error, since the limiting
# identity is ParCa-dependent and may differ in another batch.
NAMED_LIMITERS = (
	('rpoA', 'EG10893-MONOMER[c]'),
	('rpoB', 'RPOB-MONOMER[c]'),
	('rpoC', 'RPOC-MONOMER[c]'),
	('rpsQ', 'EG10916-MONOMER[c]'),
	('rpsT', 'EG10919-MONOMER[c]'),
	('rplL', 'EG10873-MONOMER[c]'),
	('rplM', 'EG10874-MONOMER[c]'),
	)

# The 36 multi-page PDFs are the expensive output (~850 MB and most of the
# runtime at 16 seeds). Set False to skip them and keep only the CSVs and the
# loop-chain figure, which is what a re-run for new columns needs.
MAKE_BIG_FIGURES = False

# THE LOOP, one row per step, for a single named gene. The last row returns
# to the first: copy number -> ... -> replication timing -> copy number. This
# is the figure to read top to bottom and then wrap around.
#
# Step 4 is deliberately the TOTAL protein count, not the free pool. The free
# pool is a residual of a nearly balanced flux and does not fall under burden
# even for the gene that holds the bottleneck.
# Event-rate rows are sparse -- rpoC initiates transcription about once every
# 300 timesteps -- so at one point per PLOT_STRIDE_SEC they render as isolated
# spikes with no readable trend. They are shown as a centred rolling mean over
# this many timesteps (~2 min), stated on the figure. Level rows are untouched.
LOOP_CHAIN_ROLL = 120

# Generation window for the loop-chain figure. The whole lineage compresses
# 24 generations into one axis, which hides the per-generation structure the
# figure exists to show; this brackets the GFP induction at generation 8 with
# two generations before and six after.
LOOP_CHAIN_WINDOW = (6, 14)

LOOP_CHAIN_GENE = 'rpoC'
LOOP_CHAIN_RATES = frozenset((
	'{g}_trs_init', '{g}_trl_init', 'rnap_formation', 'rrna_init_events',
	's50_formation',
	))

LOOP_CHAIN = (
	('gfp_mrna', '1. GFP mRNA', 'the driver'),
	('gfp_protein', '2. GFP protein', 'the driver'),
	('gfp_proteome_frac', '3. GFP proteome\nfraction', 'realised burden'),
	('{g}_copies', '4. {g} gene\ncopies', 'dosage in'),
	('{g}_trs_init', '5. {g}\ntranscript inits', 'transcription'),
	('{g}_trl_init', '6. {g}\ntranslation inits', 'translation'),
	('{g}_monomer', '7. {g} protein\n(total count)', 'protein supply'),
	('rnap_formation', '8. RNAP formation\nevents', 'machine assembly'),
	('total_rnap', '9. total RNAP', 'machine pool'),
	('rnap_budget', '10. transcription\nbudget', 'global capacity'),
	('rrna_init_events', '11. rRNA\ntranscript inits', 'downstream demand'),
	('s50_formation', '12. 50S formation\nevents', 'ribosome assembly'),
	('total_ribosome', '13. total\nribosomes', 'machine pool'),
	('active_ribosome', '14. active\nribosomes', 'translation capacity'),
	('cell_mass', '15. cell mass', 'the output'),
	('growth_rate', '16. instantaneous\ngrowth rate', 'mass accumulation'),
	('doubling_time', '17. doubling time', 'the phenotype'),
	('critical_mass_per_oric', '18. critical mass\nper oriC', 'replication trigger'),
	('n_oric', '19. origins\nper cell', 'LOOP CLOSES -> row 4'),
	)


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

# Derivative views only. At division every count halves, so d/dt has a huge
# negative spike at each generation junction that dominates the y-axis and
# hides the within-cycle signal. The spike is an artefact of the discontinuity,
# not a rate, so samples within this many seconds of a division boundary are
# blanked. Stated on every derivative figure rather than dropped silently.
DERIV_BLANK_SEC = 20.0

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
		('rprot_gene_copies', 'r-protein ALL\ngene copies', False),
		('rp_origin_gene_copies', 'r-protein ORIGIN\ngene copies (f<0.20)',
			False),
		('rp_term_gene_copies', 'r-protein TERMINUS\ngene copies (f>0.75)',
			False),
		('terminus_gene_copies', 'terminus control\ngene copies', False),
		('term_class_gene_copies', 'terminus CLASS\ngene copies', False),
		)),
	('E. transcription', (
		('rrna_init_events', 'rRNA transcript\ninits', False),
		('rnap_init_events', 'RNAP subunit\ntranscript inits', False),
		('rprot_init_events', 'r-protein ALL\ntranscript inits', False),
		('rp_origin_init_events', 'r-protein ORIGIN\ntranscript inits',
			False),
		('rp_term_init_events', 'r-protein TERMINUS\ntranscript inits',
			False),
		('rp_origin_monomer', 'r-protein ORIGIN\nprotein supply', False),
		('rp_term_monomer', 'r-protein TERMINUS\nprotein supply', False),
		('rnap_monomer', 'RNAP subunit\nprotein supply', False),
		('terminus_init_events', 'terminus control\ntranscript inits', False),
		('term_class_init_events', 'terminus CLASS\ntranscript inits', False),
		('total_init_events', 'total transcript\ninits', False),
		('rnap_budget', 'RNAP budget\n(activations)', False),
		('free_rnap', 'free RNAP\n(APORNAP)', False),
		)),
	('F. machinery and what limits it', (
		('free_16s', 'free 16S\nrRNA', False),
		('free_23s', 'free 23S\nrRNA', False),
		('free_5s', 'free 5S\nrRNA', False),
		('s30_assembly_bound', '30S assembly bound\n(ALL subunits)', False),
		('s30_limiting_idx', 'WHICH 30S subunit\nlimits', False),
		('s50_assembly_bound', '50S assembly bound\n(ALL subunits)', False),
		('s50_limiting_idx', 'WHICH 50S subunit\nlimits', False),
		('limiting_rprotein', 'limiting r-protein\n(PROTEINS only)', False),
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
	('G. rRNA chain: init -> mature -> free -> subunit', (
		('rrna_init_16s', '16S transcription\ninits', False),
		('matured_16s', '16S matured\n(supply flux)', False),
		('s30_formation', '30S formation\nevents (output)', False),
		('rrna_init_23s', '23S transcription\ninits', False),
		('matured_23s', '23S matured\n(supply flux)', False),
		('rrna_init_5s', '5S transcription\ninits', False),
		('matured_5s', '5S matured\n(supply flux)', False),
		('s50_formation', '50S formation\nevents (output)', False),
		('rrna_matured', 'all rRNA matured\n(total supply)', False),
		('rnap_formation', 'RNAP formation\nevents', False),
		)),
	('H. the 70S gate: min(free 30S, free 50S)', (
		('inactive_ribosome', 'inactive ribosome\n= min(30S, 50S)', False),
		('ribo_gate_is_50s', 'is 50S the gate?\n(1 = yes)', False),
		('ribo_gate_margin', 'gate margin\n(frac. above min)', False),
		('translation_inits', '70S formation\n(translation inits)', False),
		('ribosome_terminations', 'ribosome\nterminations', False),
		('ribo_activation_reduced', 'ribosome activation cut\n(mRNA overcrowding)', False),
		('rnap_terminations', 'RNAP\nterminations', False),
		)),
	('I. per-gene chain: transcription -> translation, by position', (
		('rp_origin_trs_init', 'origin r-prot\ntranscript inits', False),
		('rp_origin_trl_init', 'origin r-prot\ntranslation inits', False),
		('rp_term_trs_init', 'terminus r-prot\ntranscript inits', False),
		('rp_term_trl_init', 'terminus r-prot\ntranslation inits', False),
		('rnap_sub_trs_init', 'RNAP subunit\ntranscript inits', False),
		('rnap_sub_trl_init', 'RNAP subunit\ntranslation inits', False),
		)),
	('J. named bottleneck subunits (one row per step, per gene)', (
		('rpoA_copies', 'rpoA\ngene copies', False),
		('rpoA_trs_init', 'rpoA\ntranscript inits', False),
		('rpoA_trl_init', 'rpoA\ntranslation inits', False),
		('rpoA_monomer', 'rpoA\nprotein (total)', False),
		('rpoA_free', 'rpoA\nfree pool', False),
		('rpoB_copies', 'rpoB\ngene copies', False),
		('rpoB_trs_init', 'rpoB\ntranscript inits', False),
		('rpoB_trl_init', 'rpoB\ntranslation inits', False),
		('rpoB_monomer', 'rpoB\nprotein (total)', False),
		('rpoB_free', 'rpoB\nfree pool', False),
		('rpoC_copies', 'rpoC\ngene copies', False),
		('rpoC_trs_init', 'rpoC\ntranscript inits', False),
		('rpoC_trl_init', 'rpoC\ntranslation inits', False),
		('rpoC_monomer', 'rpoC\nprotein (total)', False),
		('rpoC_free', 'rpoC\nfree pool', False),
		('rpsQ_copies', 'rpsQ\ngene copies', False),
		('rpsQ_trs_init', 'rpsQ\ntranscript inits', False),
		('rpsQ_trl_init', 'rpsQ\ntranslation inits', False),
		('rpsQ_monomer', 'rpsQ\nprotein (total)', False),
		('rpsQ_free', 'rpsQ\nfree pool', False),
		('rpsT_copies', 'rpsT\ngene copies', False),
		('rpsT_trs_init', 'rpsT\ntranscript inits', False),
		('rpsT_trl_init', 'rpsT\ntranslation inits', False),
		('rpsT_monomer', 'rpsT\nprotein (total)', False),
		('rpsT_free', 'rpsT\nfree pool', False),
		('rplL_copies', 'rplL\ngene copies', False),
		('rplL_trs_init', 'rplL\ntranscript inits', False),
		('rplL_trl_init', 'rplL\ntranslation inits', False),
		('rplL_monomer', 'rplL\nprotein (total)', False),
		('rplL_free', 'rplL\nfree pool', False),
		('rplM_copies', 'rplM\ngene copies', False),
		('rplM_trs_init', 'rplM\ntranscript inits', False),
		('rplM_trl_init', 'rplM\ntranslation inits', False),
		('rplM_monomer', 'rplM\nprotein (total)', False),
		('rplM_free', 'rplM\nfree pool', False),
		)),
	('K. context', (
		('ppgpp', 'ppGpp\n(uM)', False),
		)),
	)

KEYS = tuple(k for _, rows in BLOCKS for k, _, _ in rows)


def _rolling_mean(y, w):
	"""Centred rolling mean that ignores NaN, without a pandas dependency.

	Used only for the sparse event-rate rows of the loop-chain figure: rpoC
	initiates transcription about once every 300 timesteps, so a raw trace is
	unreadable at any plotting stride.
	"""
	k = np.ones(int(w))
	ok = np.isfinite(y).astype(float)
	num = np.convolve(np.where(np.isfinite(y), y, 0.0), k, mode='same')
	den = np.convolve(ok, k, mode='same')
	return np.where(den > 0, num / np.maximum(den, 1e-9), np.nan)

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
	('view2-derivative', 'd/dt per minute; blanked within %.0f s of division'
		% DERIV_BLANK_SEC),
	('view3-smoothed-derivative',
		'd/dt smoothed over %d timesteps; blanked within %.0f s of division'
		% (SMOOTH_TIMESTEPS, DERIV_BLANK_SEC)),
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
		self._lim_rows = []

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

		# Seeds that do not exist are skipped silently by get_cells, so a batch
		# with fewer seeds than requested yields a quietly smaller result. That
		# is the failure mode this run exists to avoid -- analysing 4 lineages
		# while believing there are 16 puts the statistics back under the n=4
		# two-sided Wilcoxon floor of p=0.125. Report found-vs-requested per
		# variant, and warn if any variant is short.
		self._found = {v: sorted(s for (vv, s) in traces if vv == v)
			for v in available}
		for v in available:
			got = self._found[v]
			missing = [s for s in SEEDS if s not in got]
			msg = 'v%d: %d of %d requested seeds found' % (
				v, len(got), len(SEEDS))
			if missing:
				msg += ('  *** MISSING %s -- statistics will run on fewer '
					'lineages than intended ***' % (missing,))
			self._note(msg)

		if not traces:
			print('No readable lineages.')
			return

		self._check_gates(traces)
		self._write_by_gen(out_dir, traces)
		self._write_bottleneck(out_dir)
		self._write_limiting(out_dir)
		self._write_timeseries(out_dir, traces)
		self._write_events(out_dir, events)
		self._render(out_dir, traces, available)
		self._render_loop_chain(out_dir, traces, available)
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
		# Ribosomal-protein TUs. These OVERLAP the RNAP-subunit set -- four of
		# five RNAP-bearing TUs also carry r-proteins -- so the two classes
		# must never be summed. They are kept separate here because the
		# event-anchored test needs r-protein transcription specifically: the
		# limiting r-protein FREE POOL sits at about one integer copy and is
		# unusable at event resolution, whereas initiation counts are not.
		is_rprot = rna_data['includes_ribosomal_protein'].copy()

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

		# TRUE assembly bounds. The complexes include their rRNAs as subunits
		# -- 30S is 21 proteins PLUS 16S, 50S is 33 proteins PLUS 23S and 5S --
		# so a bound minimised over proteins alone can name the wrong subunit,
		# and does: at the unburdened control free 16S sits at ~1.1 against
		# ~5.3 for the scarcest protein, so the rRNA is the real constraint
		# there. These take the minimum over EVERY subunit of the complex, at
		# the stoichiometry the complexation reaction actually uses.
		# RNAP needs no such fix: its subunit set is only alpha/beta/beta-prime
		# and contains no RNA.
		def full_subunits(cplx):
			m = complexation.get_monomers(cplx)
			return list(m['subunitIds']), np.array(m['subunitStoich'], float)

		s30_all_ids, s30_all_stoich = full_subunits(mi.s30_full_complex)
		s50_all_ids, s50_all_stoich = full_subunits(mi.s50_full_complex)

		# Assembly OUTPUT. Scarcity of a subunit only matters if it actually
		# caps the rate at which the machine is made, and that rate is logged:
		# ComplexationListener/complexationEvents, indexed by the complex
		# produced. Without this column we can say what is scarce but never
		# whether the scarcity binds.
		# ComplexationListener/complexationEvents is REACTION-indexed, but the
		# listener's complexIDs attribute is a DIFFERENT ordering of the same
		# length (1107), so looking a complex up in complexIDs silently returns
		# the wrong column. It happens to be right for 30S and 50S and wrong for
		# RNAP (reaction 99 vs complexID 104), which is the kind of error that
		# passes every smoke test. Resolve the producing reaction from the
		# stoichiometry matrix instead, and sum if more than one produces it.
		cx = sim_data.process.complexation
		cx_mol = list(cx.molecule_names)
		cx_S = cx.stoich_matrix()
		def formation_rxns(cplx):
			if cplx not in cx_mol:
				return np.array([], int)
			row = np.asarray(cx_S[cx_mol.index(cplx), :]).ravel()
			return np.where(row > 0)[0]
		self._cplx_targets = [(formation_rxns(c), lab) for c, lab in (
			(mi.s30_full_complex, 's30_formation'),
			(mi.s50_full_complex, 's50_formation'),
			(mi.full_RNAP, 'rnap_formation'))]
		self._n_cx_rxns = int(cx_S.shape[1])
		for rx, lab in self._cplx_targets:
			print('  %s <- complexation reaction(s) %s' % (lab, list(rx)))
		# rRNA SUPPLY. Free rRNA is set by maturation, not transcription
		# directly, so this separates "scarce because little is made" from
		# "scarce because it is consumed instantly".
		# Split by species, and take EVERY operon's copy of each. 16S feeds
		# 30S while 23S and 5S feed 50S, so a single summed column cannot say
		# which machine a change in supply acts on. An earlier version sliced
		# [:1] from each group and so counted 3 of the 22 rRNA molecules.
		self._mature_rrna = dict(
			matured_16s=list(mg.s30_16s_rRNA),
			matured_23s=list(mg.s50_23s_rRNA),
			matured_5s=list(mg.s50_5s_rRNA))

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

		# Ribosomal proteins SPLIT BY POSITION. The r-protein set spans nearly
		# the whole copy-number channel -- 13 TUs at f<0.20 lose ~9% of
		# promoter share relative to the genome average while 8 TUs at f>0.75
		# GAIN ~9% -- so summing all 42 cancels the signal almost exactly.
		# That is the same L1-normalisation cancellation that has caught this
		# study before.
		#
		# Split by POSITION, never by which gene happens to be limiting.
		# Limiting-ness is an outcome (and an unstable one: 18 different
		# r-proteins hold the title at v7), so selecting on it would bias the
		# analysis. Position is fixed in the genome and cannot respond to
		# burden.
		#
		# The two groups are matched on everything else that matters --
		# co-regulation, expression regime, assembly demand, and the global
		# consumption rate -- which makes origin-vs-terminus r-proteins a
		# tighter control than the generic terminus class.
		# Define the positional groups ONCE, at gene level, and derive both the
		# TU list and the monomer list from the same set. Selecting TUs and
		# monomers independently gave different gene sets (13 vs 10 origin,
		# 8 vs 3 terminus), which would make transcription and protein supply
		# non-comparable -- the whole point is that they describe the same
		# genes.
		rp_genes_o = [b for b in bottleneck
			if b['group'] != 'rnap_subunit' and b['frac'] < 0.20]
		rp_genes_t = [b for b in bottleneck
			if b['group'] != 'rnap_subunit' and b['frac'] > 0.75]
		rp_origin = np.unique(np.concatenate([b['tus'] for b in rp_genes_o])
			) if rp_genes_o else np.array([], int)
		rp_term = np.unique(np.concatenate([b['tus'] for b in rp_genes_t])
			) if rp_genes_t else np.array([], int)

		# Terminus CLASS for the event-anchored control. The single terminus
		# control gene initiates ~once per 1000 timesteps (mean 0.001), so it
		# works as a copy-number control but not as a transcription outcome.
		# A class of terminus-proximal TUs, excluding machinery and the
		# construct, gives a count large enough to compare against.
		term_class = np.where((frac > 0.90)
			& ~(is_rrna | is_rnap | is_rprot | is_new_tu))[0]

		# Monomer-id lists for the positional groups, so PROTEIN SUPPLY can be
		# read the same way transcription is. Free pools are consumed on
		# arrival and sit at ~1 integer copy, which is unusable; TOTAL monomer
		# counts (MonomerCounts, including protein already inside an assembled
		# complex) accumulate, so they integrate the transcription pulse.
		rp_o_mon = [b['monomer'] for b in rp_genes_o]
		rp_t_mon = [b['monomer'] for b in rp_genes_t]
		rnap_mon = [b['monomer'] for b in bottleneck
			if b['group'] == 'rnap_subunit']
		# Index tables for the per-gene FLUX columns. Transcription flux comes
		# from RnapData/rna_init_event_per_cistron, which is cistron-indexed and
		# so needs no TU mapping; translation flux from
		# RibosomeData/ribosome_init_event_per_monomer. Together they give the
		# transcript-initiation -> protein-synthesis step measured as RATES on
		# the same genes the positional contrast is defined on. Monomer COUNTS
		# integrate history and so cannot date a change; these can.
		mono_index = {m: i for i, m in enumerate(monomer_data['id'])}
		rnap_genes = [b for b in bottleneck if b['group'] == 'rnap_subunit']
		# Per-gene index quartet for each named bottleneck subunit. Cistron index
		# drives transcription, monomer index drives translation, the TU set
		# drives promoter copy number and the bulk id drives the free pool.
		named = []
		for lab, mid in NAMED_LIMITERS:
			cis = mono_to_cistron.get(mid)
			ci = cistron_index.get(cis) if cis else None
			mi_ = mono_index.get(mid)
			tus = (np.where(cistron_tu[ci, :].toarray().ravel() > 0)[0]
				if ci is not None else np.array([], int))
			named.append(dict(lab=lab, monomer=mid, cistron=cis, ci=ci, mi=mi_,
				tus=tus, ok=(ci is not None and mi_ is not None and tus.size > 0)))
			print('  named limiter %-5s %-22s %s'
				% (lab, mid, 'resolved (%d TUs)' % tus.size
					if named[-1]['ok'] else '*** NOT FOUND -- columns will be NaN ***'))
		def cis_idx(genes):
			return np.array([cistron_index[b['cistron']] for b in genes
				if b['cistron'] in cistron_index], int)
		def mon_idx(genes):
			return np.array([mono_index[b['monomer']] for b in genes
				if b['monomer'] in mono_index], int)

		print('  positional r-protein groups: origin %d genes / %d TUs, '
			'terminus %d genes / %d TUs (same gene set drives both '
			'transcription and protein supply)'
			% (len(rp_o_mon), rp_origin.size, len(rp_t_mon), rp_term.size))

		# Metadata for the limiting-fraction table: every subunit of every
		# machine, with its stoichiometry, its chromosomal position and the
		# TUs that transcribe it, so "how often is this the bottleneck" sits
		# in the same row as "how many copies of its gene are there".
		mono_to_tus = {b['monomer']: b['tus'] for b in bottleneck}
		mono_to_frac = {b['monomer']: b['frac'] for b in bottleneck}
		def subunit_meta(ids, stoich, cplx):
			out = []
			for k, (mid, st) in enumerate(zip(ids, stoich)):
				is_rna = 'RNA' in mid
				tus = mono_to_tus.get(mid)
				if tus is None and is_rna:
					# rRNA subunit: find the TUs whose cistrons carry it
					cis = mid.split('[')[0]
					ci = cistron_index.get(cis)
					tus = (np.where(
						cistron_tu[ci, :].toarray().ravel() > 0)[0]
						if ci is not None else np.array([], int))
				out.append(dict(cplx=cplx, idx=k, subunit=mid,
					stoich=float(st), is_rna=bool(is_rna),
					frac=mono_to_frac.get(mid, float(frac[tus].mean())
						if tus is not None and len(tus) else float('nan')),
					tus=tus if tus is not None else np.array([], int)))
			return out

		limiting_meta = (subunit_meta(s30_all_ids, s30_all_stoich, '30S')
			+ subunit_meta(s50_all_ids, s50_all_stoich, '50S')
			+ subunit_meta(rnap_sub_ids, rnap_sub_stoich, 'RNAP'))

		ctx = dict(
			rna_ids=rna_ids,
			limiting_meta=limiting_meta,
			term_class_tu=term_class,
			rp_origin_mon=rp_o_mon,
			rp_term_mon=rp_t_mon,
			rnap_mon=rnap_mon,
			rp_origin_cis=cis_idx(rp_genes_o),
			rp_term_cis=cis_idx(rp_genes_t),
			rnap_cis=cis_idx(rnap_genes),
			rp_origin_mon_idx=mon_idx(rp_genes_o),
			rp_term_mon_idx=mon_idx(rp_genes_t),
			rnap_mon_idx=mon_idx(rnap_genes),
			named=named,
			rp_origin_tu=rp_origin,
			rp_term_tu=rp_term,
			new_tu=np.where(is_new_tu)[0],
			new_mrna_ids=new_mrna_ids,
			new_monomer_ids=new_monomer_ids,
			gfp_monomer_mass_fg=gfp_mass,
			rrna_tu=np.where(is_rrna)[0],
			rnap_tu=np.where(is_rnap)[0],
			rprot_tu=np.where(is_rprot)[0],
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
			s30_all_ids=s30_all_ids, s30_all_stoich=s30_all_stoich,
			s50_all_ids=s50_all_ids, s50_all_stoich=s50_all_stoich,
			)
		overlap = np.intersect1d(ctx['rnap_tu'], ctx['rprot_tu']).size
		print('Construct TUs %s; rRNA TUs %d; RNAP-subunit TUs %d; '
			'r-protein TUs %d (%d overlap the RNAP set -- do not sum); '
			'terminus control %s at f=%.3f; terminus class %d TUs at f>0.90; '
			'r-protein split %d origin (f<0.20) / %d terminus (f>0.75)'
			% (ctx['new_tu'], ctx['rrna_tu'].size, ctx['rnap_tu'].size,
				ctx['rprot_tu'].size, overlap, ctx['term_id'],
				ctx['term_frac'], ctx['term_class_tu'].size,
				ctx['rp_origin_tu'].size, ctx['rp_term_tu'].size))
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
			limpack = one.pop('_limiting', None)
			if limpack is not None:
				lim, nts, cop = limpack
				for m in ctx['limiting_meta']:
					key = (m['cplx'], m['idx'])
					if key not in lim:
						continue
					cnt, mean_, med_, frac_lt1 = lim[key]
					pc = (float(cop[:, m['tus']].sum(axis=1).mean())
						if len(m['tus']) else float('nan'))
					self._lim_rows.append(dict(
						variant=variant, seed=seed, generation=gen,
						complex=m['cplx'], subunit=m['subunit'],
						is_rna=int(m['is_rna']), stoich=m['stoich'],
						replichore_frac=m['frac'], n_timesteps=nts,
						limiting_frac=cnt / nts if nts else float('nan'),
						free_per_stoich_mean=mean_,
						free_per_stoich_median=med_,
						frac_below_one=frac_lt1, promoter_copies=pc))
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
		out['rprot_gene_copies'] = take(copies, ctx['rprot_tu'], np.mean)
		out['rrna_init_events'] = take(inits, ctx['rrna_tu'], np.sum)
		out['rnap_init_events'] = take(inits, ctx['rnap_tu'], np.sum)
		out['rprot_init_events'] = take(inits, ctx['rprot_tu'], np.sum)
		# The terminus control's own transcription. At a replication
		# initiation its copy number stays exactly flat while origin-proximal
		# copies step up, so this is the within-cell negative control for the
		# event-anchored test: same cell, same instant, same mass change, no
		# copy step.
		out['terminus_init_events'] = take(inits, ctx['term_tu'], np.sum)
		out['rp_origin_gene_copies'] = take(
			copies, ctx['rp_origin_tu'], np.mean)
		out['rp_term_gene_copies'] = take(copies, ctx['rp_term_tu'], np.mean)
		out['rp_origin_init_events'] = take(
			inits, ctx['rp_origin_tu'], np.sum)
		out['rp_term_init_events'] = take(inits, ctx['rp_term_tu'], np.sum)
		out['term_class_gene_copies'] = take(
			copies, ctx['term_class_tu'], np.mean)
		out['term_class_init_events'] = take(
			inits, ctx['term_class_tu'], np.sum)
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

		# Protein supply per positional group. Total counts, not free pools.
		def mon_sum(ids):
			idx = [mono_ids.index(i) for i in ids if i in mono_ids]
			return (counts[:, idx].sum(axis=1).astype(float) if idx
				else np.full(n, np.nan))

		out['rp_origin_monomer'] = mon_sum(ctx['rp_origin_mon'])
		out['rp_term_monomer'] = mon_sum(ctx['rp_term_mon'])
		out['rnap_monomer'] = mon_sum(ctx['rnap_mon'])

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
		(free_rnap, s30, s50, s16, s23, s5, r_prot, rnap_sub,
				s30_all, s50_all) = read_bulk_molecule_counts(sim_out, (
				ctx['free_rnap'], ctx['s30'], ctx['s50'], ctx['s16'],
				ctx['s23'], ctx['s5'], ctx['r_prot_ids'], ctx['rnap_sub_ids'],
				ctx['s30_all_ids'], ctx['s50_all_ids']))
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

		# The true bounds, over every subunit including the rRNAs.
		p30 = np.atleast_2d(s30_all).astype(float) / ctx['s30_all_stoich']
		p50 = np.atleast_2d(s50_all).astype(float) / ctx['s50_all_stoich']
		out['s30_assembly_bound'] = p30.min(axis=1)
		out['s30_limiting_idx'] = p30.argmin(axis=1).astype(float)
		out['s50_assembly_bound'] = p50.min(axis=1)
		out['s50_limiting_idx'] = p50.argmin(axis=1).astype(float)

		# Per-subunit limiting counts. The argmin columns above are indices,
		# so averaging them over a generation is meaningless; what is wanted
		# is the FRACTION OF TIMESTEPS each subunit holds the minimum. Count
		# it here and aggregate per generation in _write_limiting.
		lim = {}
		for cplx, per, meta_n in (('30S', p30, p30.shape[1]),
				('50S', p50, p50.shape[1]),
				('RNAP', per_rnap, per_rnap.shape[1])):
			am = per.argmin(axis=1)
			cnt = np.bincount(am, minlength=meta_n)
			for k in range(meta_n):
				lim[(cplx, k)] = (int(cnt[k]),
					float(per[:, k].mean()), float(np.median(per[:, k])),
					float((per[:, k] < 1).mean()))
		out['_limiting'] = (lim, n, copies)

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

		# Assembly output rates, and the supply/demand terms around them.
		try:
			cl = TableReader(os.path.join(sim_out, 'ComplexationListener'))
			cev = cl.readColumn('complexationEvents')
			if cev.shape[1] != self._n_cx_rxns:
				raise ValueError('complexationEvents has %d columns but the '
					'complexation network has %d reactions; the '
					'reaction indexing assumption is broken'
					% (cev.shape[1], self._n_cx_rxns))
			for rx, lab in self._cplx_targets:
				out[lab] = (cev[:, rx].sum(axis=1).astype(float) if rx.size
					else np.full(n, np.nan))
		except Exception:  # noqa: BLE001 - listener absent in older batches
			for _, lab in self._cplx_targets:
				out[lab] = np.full(n, np.nan)
		try:
			ml = TableReader(os.path.join(sim_out, 'RnaMaturationListener'))
			mids = list(ml.readAttribute('mature_rna_ids'))
			mgen = ml.readColumn('mature_rnas_generated')
			tot = np.zeros(n)
			for lab, species in self._mature_rrna.items():
				idx = [mids.index(i) for i in species if i in mids]
				v = (mgen[:, idx].sum(axis=1).astype(float) if idx
					else np.full(n, np.nan))
				out[lab] = v
				tot = tot + np.nan_to_num(v)
			out['rrna_matured'] = tot
		except Exception:  # noqa: BLE001
			for lab in self._mature_rrna:
				out[lab] = np.full(n, np.nan)
			out['rrna_matured'] = np.full(n, np.nan)
		out['translation_inits'] = ribo.readColumn('didInitialize').astype(
			float)

		# THE 70S GATE. polypeptide_initiation.py:90 activates
		# activationProb * min(free 30S, free 50S), so the smaller free subunit
		# pool is not merely correlated with 70S formation -- it is the term in
		# the rate law. Log which subunit holds the minimum and how far apart
		# the two are, because the median gap is only a few percent and could
		# flip under burden.
		out['inactive_ribosome'] = np.minimum(out['free_30s'], out['free_50s'])
		out['ribo_gate_is_50s'] = (out['free_50s']
			< out['free_30s']).astype(float)
		with np.errstate(invalid='ignore', divide='ignore'):
			out['ribo_gate_margin'] = np.where(out['inactive_ribosome'] > 0,
				(np.maximum(out['free_30s'], out['free_50s'])
					- out['inactive_ribosome']) / out['inactive_ribosome'],
				np.nan)
		# NOTE: this flag is mRNA OVERCROWDING, not subunit shortage --
		# polypeptide_initiation.py:183 sets it when per-transcript initiation
		# probability exceeds the ribosome-footprint maximum and cannot be
		# rescaled away. It is a second, independent constraint on translation,
		# so it is worth logging, but it does NOT test the 30S/50S gate.
		out['ribo_activation_reduced'] = ribo.readColumn(
			'is_n_ribosomes_to_activate_reduced').astype(float)
		# Termination rates, so formation can be read as NET accumulation
		# rather than gross flux. Ribosomes and RNAPs recycle, so initiation
		# alone does not say whether the pool is growing.
		out['ribosome_terminations'] = ribo.readColumn(
			'didTerminate').astype(float)
		out['rnap_terminations'] = rnap.readColumn('didTerminate').astype(float)
		# rRNA transcription initiation per species, to pair with maturation
		# per species and complete the rRNA chain init -> mature -> free -> 30S/50S.
		for lab, col in (('rrna_init_16s', 'rRNA16S_initiated'),
			('rrna_init_23s', 'rRNA23S_initiated'),
			('rrna_init_5s', 'rRNA5S_initiated')):
			try:
				out[lab] = ribo.readColumn(col).astype(float)
			except Exception:  # noqa: BLE001 - absent in older batches
				out[lab] = np.full(n, np.nan)

		# Per-gene transcription and translation FLUX on the positional groups.
		# This is the transcript -> protein step: if copy loss reaches protein
		# supply, origin-proximal r-protein genes must lose translation
		# initiations relative to terminus-proximal ones, on the same gene set
		# whose transcription is measured directly above.
		try:
			cis_ev = rnap.readColumn('rna_init_event_per_cistron')
			mon_ev = ribo.readColumn('ribosome_init_event_per_monomer')
			for lab in ('rp_origin', 'rp_term', 'rnap_sub'):
				ci = ctx[lab.replace('rnap_sub', 'rnap') + '_cis']
				mo = ctx[lab.replace('rnap_sub', 'rnap') + '_mon_idx']
				out[lab + '_trs_init'] = (
					cis_ev[:, ci].sum(axis=1).astype(float) if ci.size
					else np.full(n, np.nan))
				out[lab + '_trl_init'] = (
					mon_ev[:, mo].sum(axis=1).astype(float) if mo.size
					else np.full(n, np.nan))
		except Exception:  # noqa: BLE001 - absent in older batches
			for lab in ('rp_origin', 'rp_term', 'rnap_sub'):
				out[lab + '_trs_init'] = np.full(n, np.nan)
				out[lab + '_trl_init'] = np.full(n, np.nan)

		# Named bottleneck subunits, one quartet each. These are the columns
		# that let the chain be followed on a single gene instead of a group
		# mean: dosage -> transcription -> translation -> protein.
		named = ctx.get('named') or []
		try:
			cis_ev = rnap.readColumn('rna_init_event_per_cistron')
			mon_ev = ribo.readColumn('ribosome_init_event_per_monomer')
		except Exception:  # noqa: BLE001
			cis_ev = mon_ev = None
		for g in named:
			lab = g['lab']
			if g['ok'] and cis_ev is not None:
				out[lab + '_copies'] = copies[:, g['tus']].sum(axis=1).astype(float)
				out[lab + '_trs_init'] = cis_ev[:, g['ci']].astype(float)
				out[lab + '_trl_init'] = mon_ev[:, g['mi']].astype(float)
			else:
				for suf in ('_copies', '_trs_init', '_trl_init'):
					out[lab + suf] = np.full(n, np.nan)
			# TOTAL monomer count, which accumulates, alongside the FREE pool.
			# The free pool is a residual of a nearly balanced flux and does not
			# track supply -- rpoC is the RNAP bottleneck yet its free pool does
			# not fall under burden. Total count is the supply instrument; the
			# free pool is kept only for the assembly-bound arithmetic.
			out[lab + '_monomer'] = mon_sum([g['monomer']])

		ok_named = [g for g in named if g['ok']]
		if ok_named:
			(free_named,) = read_bulk_molecule_counts(sim_out,
				([g['monomer'] for g in ok_named],))
			free_named = np.atleast_2d(free_named).astype(float)
			for j, g in enumerate(ok_named):
				out[g['lab'] + '_free'] = free_named[:, j]
		for g in named:
			if not g['ok']:
				out[g['lab'] + '_free'] = np.full(n, np.nan)
				out[g['lab'] + '_monomer'] = np.full(n, np.nan)

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
			# Mask the division discontinuities before differentiating, so the
			# halving spike neither appears nor leaks into the smoother.
			near_div = np.zeros(t.size, bool)
			w = DERIV_BLANK_SEC / 60.0
			for a, b, _g in tr['gen_bounds']:
				near_div |= np.abs(t - b) <= w
				near_div |= np.abs(t - a) <= w
			out = {}
			for k in KEYS:
				with np.errstate(invalid='ignore', divide='ignore'):
					d = np.gradient(tr[k]) / dt
				d = np.where(near_div, np.nan, d)
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
		if not MAKE_BIG_FIGURES:
			self._note('MAKE_BIG_FIGURES is False -- skipped the 36 multi-page '
				'PDFs; CSVs and the loop-chain figure are unaffected')
			return
		for variant in variants:
			seeds = [s for s in FIGURE_SEEDS if (variant, s) in traces]
			if not seeds:
				continue
			for view, view_desc in VIEWS:
				for wname, wrange in WINDOWS:
					self._one_pdf(out_dir, variant, seeds, traces, view,
						view_desc, wname, wrange)

	def _render_loop_chain(self, out_dir, traces, variants):
		"""One page per lineage: every step of the loop as its own row.

		Two passes -- raw counts and a smoothed derivative -- because the level
		and the rate answer different questions. The level shows where a
		quantity ends up; the derivative shows when it turned, which is what
		ordering an eight-step chain needs.
		"""
		g = LOOP_CHAIN_GENE
		rows = [(k.format(g=g), lab.format(g=g), note)
			for k, lab, note in LOOP_CHAIN]
		missing = [k for k, _, _ in rows if k not in KEYS]
		if missing:
			self._note('loop-chain figure skipped, missing columns: %s'
				% (missing,))
			return
		matplotlib.rcParams.update({'font.size': 8,
			'axes.spines.top': False, 'axes.spines.right': False})
		for variant in variants:
			seeds = [sd for sd in FIGURE_SEEDS if (variant, sd) in traces]
			if not seeds:
				continue
			for mode in ('counts', 'derivative'):
				path = os.path.join(out_dir, 'v%02d_loopchain_%s_%s_gens%02d-%02d.pdf'
					% (variant, g, mode,
						LOOP_CHAIN_WINDOW[0], LOOP_CHAIN_WINDOW[1]))
				with PdfPages(path) as pdf:
					for seed in seeds:
						fig = self._loop_chain_page(traces[(variant, seed)],
							rows, variant, seed, mode, g)
						if fig is None:
							continue
						pdf.savefig(fig)
						if seed == seeds[0]:
							fig.savefig(os.path.join(out_dir,
								'v%02d_loopchain_%s_%s_gens%02d-%02d.png'
								% (variant, g, mode, LOOP_CHAIN_WINDOW[0],
									LOOP_CHAIN_WINDOW[1])), dpi=110,
								bbox_inches='tight')
						plt.close(fig)
				print('  %s' % os.path.basename(path))

	def _loop_chain_page(self, tr, rows, variant, seed, mode, gene):
		x, ys = self._transform(tr, 'view1-counts' if mode == 'counts'
			else 'view3-smoothed-derivative')
		stride = self._stride(tr)
		keep = np.zeros(x.size, bool)
		keep[::stride] = True
		if LOOP_CHAIN_WINDOW is not None:
			gens = np.full(tr['time'].size, -1)
			for a, b, g in tr['gen_bounds']:
				gens[(tr['time'] >= a) & (tr['time'] <= b)] = g
			keep &= ((gens >= LOOP_CHAIN_WINDOW[0])
				& (gens <= LOOP_CHAIN_WINDOW[1]))
		if not np.any(keep):
			return None
		fig, axes = plt.subplots(len(rows), 1, sharex=True,
			figsize=(16, 1.25 * len(rows)))
		for i, (ax, (key, lab, note)) in enumerate(zip(axes, rows)):
			y = ys[key]
			raw_key = LOOP_CHAIN[i][0]
			if mode == 'counts' and raw_key in LOOP_CHAIN_RATES:
				y = _rolling_mean(y, LOOP_CHAIN_ROLL)
			ax.plot(x[keep], y[keep], lw=.95, color='#4a3aa7')
			ax.set_ylabel(lab, fontsize=7.5, rotation=0, ha='right',
				va='center', labelpad=10)
			ax.tick_params(labelsize=7)
			ax.margins(x=.01)
			# The step name sits inside the axes so the row reads as a stage of
			# the chain rather than as one more unlabelled trace.
			ax.text(.998, .90, note, transform=ax.transAxes, ha='right',
				va='top', fontsize=7, color='#8a7f2f' if i == len(rows) - 1
				else '#96999e',
				fontweight='bold' if i in (0, len(rows) - 1) else 'normal')
			if key == 'critical_mass_per_oric' and mode == 'counts':
				ax.axhline(1.0, color='#b3123c', lw=1.0)
			if mode == 'derivative':
				ax.axhline(0.0, color='#c9ced4', lw=.7)
			self._mark(ax, tr, LOOP_CHAIN_WINDOW, label_induction=(i == 0))
		# Clamp to the plotted span. Masking the data does not move the axis,
		# so without this a 9-generation window renders on a 24-generation
		# axis with the data squeezed into a third of the width.
		xv = x[keep]
		if xv.size:
			pad = 0.01 * (xv.max() - xv.min() or 1.0)
			for ax in axes:
				ax.set_xlim(xv.min() - pad, xv.max() + pad)
		axes[-1].set_xlabel('time (min, continuous across generations)')
		# Build the title in two steps. Adjacent string literals concatenate
		# BEFORE % is applied, so a second % operator inside the literal block
		# tries to fill every placeholder in the whole title.
		roll_note = ('event-rate rows are a %d-timestep rolling mean; '
			'level rows are raw' % LOOP_CHAIN_ROLL) if mode == 'counts' else (
			'derivative smoothed over %d timesteps' % SMOOTH_TIMESTEPS)
		fig.suptitle('THE LOOP, STEP BY STEP -- %s  |  variant %d, seed %d  |  %s'
			'  |  generations %s'
			'\nread top to bottom, then wrap: row 19 feeds row 4.  %s.'
			'\nSOLID RED = GFP induction; dashed grey = division; '
			'dotted orange = replication initiation'
			% (gene, variant, seed, mode,
				'%d-%d' % LOOP_CHAIN_WINDOW if LOOP_CHAIN_WINDOW
				else 'all', roll_note), fontsize=10, y=.997)
		fig.tight_layout(rect=(0, 0, 1, .988))
		return fig

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
						self._mark(ax, tr, wrange,
							label_induction=(key == KEYS[0]))
				axes[-1].set_xlabel('GFP proteome mass fraction' if phase
					else 'time (min, continuous across generations)')
				fig.suptitle('variant %d, seed %d  |  %s  |  %s  |  %s\n'
					'SOLID RED = GFP induction; dashed grey = division; '
					'dotted orange = replication initiation%s  |  '
					'plotted every %d timesteps '
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

	def _induction_time(self, tr):
		"""Time at which the construct switches on, in minutes.

		Induction is scheduled into internal_shift_dict and applied at the
		start of generation INDUCTION_GEN, so it is the left edge of that
		generation's bounds rather than anything logged.
		"""
		for a, b, g in tr['gen_bounds']:
			if g == INDUCTION_GEN:
				return a
		return None

	def _mark(self, ax, tr, wrange, label_induction=False):
		"""Division, replication-initiation and GFP-induction markers."""
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
		# GFP induction, drawn last and heaviest so it reads over the top of
		# the division and replication marks.
		ti = self._induction_time(tr)
		if ti is not None and (wrange is None
				or wrange[0] <= INDUCTION_GEN <= wrange[1]):
			ax.axvline(ti, color='#b3123c', lw=2.2, ls='-', alpha=.9,
				zorder=6)
			if label_induction:
				ax.text(ti, 1.04, ' GFP induced (gen %d)' % INDUCTION_GEN,
					transform=ax.get_xaxis_transform(), fontsize=9,
					color='#b3123c', fontweight='650', ha='left',
					va='bottom', zorder=7)

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

	def _write_limiting(self, out_dir):
		"""How often each subunit is the bottleneck, beside its gene dosage.

		This is the table that answers the three questions together: which
		component limits each machine, what fraction of the time it does so,
		and how many gene copies it has while doing it. One row per
		(variant, seed, generation, complex, subunit) over every subunit of
		30S, 50S and RNAP -- rRNAs included, which the protein-only version
		could not represent.
		"""
		if not self._lim_rows:
			return
		path = os.path.join(out_dir, 'loop_lineage_limiting.csv')
		cols = ['variant', 'seed', 'generation', 'complex', 'subunit',
			'is_rna', 'stoich', 'replichore_frac', 'n_timesteps',
			'limiting_frac', 'free_per_stoich_mean',
			'free_per_stoich_median', 'frac_below_one', 'promoter_copies']
		with open(path, 'w', newline='') as fh:
			w = csv.DictWriter(fh, fieldnames=cols)
			w.writeheader()
			for r in self._lim_rows:
				w.writerow({c: ('%.6g' % r[c] if isinstance(r[c], float)
					else r[c]) for c in cols})
		print('  %s (%d rows, %.2f MB)' % (os.path.basename(path),
			len(self._lim_rows), os.path.getsize(path) / 1e6))

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
				if seed not in TRACE_SEEDS:
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
			# Report the seeds actually READ, per variant, not the ones requested.
			# An earlier version printed list(SEEDS) unconditionally, so a 2-seed
			# batch produced a manifest claiming 16.
			found = getattr(self, '_found', {})
			fh.write('\nSEED COVERAGE\n  requested: %s\n' % (list(SEEDS),))
			for v in sorted(found):
				got = found[v]
				missing = [x for x in SEEDS if x not in got]
				fh.write('  v%d FOUND %d of %d: %s%s\n'
					% (v, len(got), len(SEEDS), got,
						('   *** MISSING %s ***' % (missing,)) if missing else ''))
			fh.write('  per-timestep trace written for the found seeds in %s\n'
				'  figures rendered for the found seeds in %s\n'
				% (list(TRACE_SEEDS), list(FIGURE_SEEDS)))
			fh.write('\nFIGURE THINNING\n  figures plot one point per ~%.0f s '
				'(PLOT_STRIDE_SEC); loop_lineage_timeseries.csv.gz keeps every '
				'timestep at full resolution\n' % PLOT_STRIDE_SEC)
			fh.write('\nLIMITING-SUBUNIT INDEX LEGEND\n')
			fh.write('  limiting_rprotein_idx -> index into:\n')
			for i, mid in enumerate(ctx['r_prot_ids']):
				fh.write('    %3d  %s (stoich %g)\n'
					% (i, mid, ctx['r_prot_stoich'][i]))
			for nm, ids, st in (('s30_limiting_idx', ctx['s30_all_ids'],
						ctx['s30_all_stoich']),
					('s50_limiting_idx', ctx['s50_all_ids'],
						ctx['s50_all_stoich'])):
				fh.write('  %s -> index into (RNA entries marked):\n' % nm)
				for i, mid in enumerate(ids):
					tag = '  <-- rRNA' if 'RNA' in mid else ''
					fh.write('    %3d  %s (stoich %g)%s\n'
						% (i, mid, st[i], tag))
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
