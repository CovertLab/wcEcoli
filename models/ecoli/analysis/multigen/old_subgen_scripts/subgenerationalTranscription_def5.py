"""
Definition-5 rewrite of subgenerationalTranscription.py (Figure 5B/5E/5F/5G).

Three changes from the original, all to align it with Definition 5:

  1. A gene counts as transcribed in a generation when it produced at least one
     COMPLETED transcript (TranscriptElongationListener/
     countRnaCistronSynthesized > 0), not merely when an mRNA molecule was
     present (which also captures transcripts inherited from the mother and
     initiations later lost to tRNA attenuation).
  2. The first IGNORE_FIRST_N_GENS generations are dropped as burn-in, so the
     per-seed frequency is not contaminated by initial-condition transients.
     The original averaged over every generation including startup.
  3. The original ran on the seed-0 seed unconditionally. Seed 0 is not
     necessarily healthy -- a seed that stalls out (cells pinned at the
     180-minute length cap) still passes ap.get_cells(only_successful=True),
     because those cells did write daughter state, and its collapsing
     transcription inflates the subgenerational fraction. This version instead
     plots the first N_SEEDS_TO_PLOT STRICT-successful seeds (completed
     every generation, no cell at the doubling cap), the same gate every other
     def-5 analysis uses.

  4. Genes are labelled by the CANONICAL Definition-5 classification -- the
     cohort-wide def5_CI categories from sc.canonical_def5_classification(), i.e.
     the 95% CI of the per-seed rate against 1 transcript/generation. Earlier
     versions of this script classified by the per-seed POINT ESTIMATE
     (mean == 0 / 0 < mean < 1 / mean >= 1), which silently disagreed with every
     def-5 table: on sim set 1 it called 1652-1687 genes subgen per seed where
     the cohort CI calls 1859, and inflated never_expressed from 135 to ~330
     (a single seed never fires many genes that some seed does). The
     per-seed point estimate is still printed as a diagnostic, but it no longer
     labels anything.

Definition 5 always means the CI form. Categories are sc.CATEGORIES -- subgen /
possibly_subgen / not_subgen / never_expressed -- coloured from sc.PALETTE, so
these panels and the per-gene tables always describe the same gene sets. Because
the categories are cohort-wide, all plotted seeds report the same counts; what
varies per seed is the frequency scatter and the event raster.

Requires the raw extraction (subgen_extract.py) to have been run on the cohort
first: the classification is read from the COHORT plotOut directory, one level
above this multigen plotOut. One figure set per plotted seed, seed-suffixed.
"""

import os

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from wholecell.io.tablereader import TableReader
from wholecell.utils import constants
from wholecell.utils.sparkline import whitePadSparklineAxis
from wholecell.analysis.analysis_tools import exportFigure
from models.ecoli.analysis import multigenAnalysisPlot
from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from models.ecoli.analysis.cohort import subgen_helper_functions as sc

# Number of post-burn-in generations to draw in the transcription-event raster.
RASTER_N_GENS = 5
# How many strict-successful seeds to plot (one figure set each). This figure
# is inherently single-seed, so a few seeds give a sense of the
# seed-to-seed spread. Mirrors N_SEEDS_TO_PLOT in
# subgen_monomer_dynamics_def5.py.
N_SEEDS_TO_PLOT = 3

# Category colors and labels come from subgen_helper_functions (sc.PALETTE / sc.CAT_LABEL)
# so these panels match every other def-5 figure. The 5B raster keeps two rows --
# the confident calls -- because an event raster of `possibly_subgen` genes would
# be read as a claim the CI explicitly declines to make.
RASTER_CATEGORIES = ['not_subgen', 'subgen']


def remove_xaxis(axis):
	axis.spines['bottom'].set_visible(False)
	axis.tick_params(bottom=False, axis='x', labelbottom=False)
	axis.set_xlabel('')


class Plot(multigenAnalysisPlot.MultigenAnalysisPlot):
	def do_plot(self, seedOutDir, plotOutDir, plotOutFileName, simDataFile,
			validationDataFile, metadata):
		# The multigen framework hands us one seed directory, but this figure is
		# only meaningful on a seed that stayed healthy for the whole run, so
		# the seed choice is made here rather than by the caller. self.ap is
		# scoped to seedOutDir; a cohort view over the parent variant directory
		# is needed both to reach the other seeds' simOut and to get the
		# cohort-wide generation count -- self.ap.n_generation counts only the
		# generations THIS seed produced, which would make a truncated seed
		# look like it had completed every generation.
		variant_dir = os.path.dirname(os.path.normpath(seedOutDir))
		cohort_plot_out = os.path.join(variant_dir, constants.PLOTOUT_DIR)
		requested_seed = self._requested_seed(seedOutDir)
		coh_ap = AnalysisPaths(variant_dir, cohort_plot=True)
		n_generation = coh_ap.n_generation
		if n_generation <= 1:
			print('Skipping -- only runs for multigen')
			return

		# The canonical def-5 labels are cohort-wide, so they come from the COHORT
		# plotOut (one level up from this multigen plotOut -- see
		# runscripts/manual/analysisMultigen.py, which appends the seed directory).
		try:
			clf = sc.canonical_def5_classification(cohort_plot_out)
		except FileNotFoundError as e:
			print('Cannot classify genes: %s' % e)
			print('This figure now uses the canonical cohort def5_CI categories, '
				'so the raw extraction must run first:\n'
				'  python runscripts/manual/analysisCohort.py '
				'--plot subgen_extract.py <sim_dir>')
			return
		if clf['n_seeds'] == 0:
			print('No successful seeds in the cohort classification. Skipping.')
			return
		cat_counts = sc.category_counts(clf['stats']['cat'])
		print('Canonical def5_CI classification over %d successful seeds: %s'
			% (clf['n_seeds'],
				', '.join('%s=%d' % (c, cat_counts[c]) for c in sc.CATEGORIES)))

		successful, n_seeds, reasons = self._seed_success(
			variant_dir, coh_ap, n_generation)
		if not successful:
			print('WARNING: no strict-successful seed among %d seeds (needs '
				'all %d generations with no cell at the %g-min doubling cap); '
				'skipping.' % (n_seeds, n_generation, sc.MAX_DOUBLING_MIN))
			return

		ordered = sorted(successful)
		if requested_seed in successful:
			seeds_to_plot = [requested_seed] + [
				s for s in ordered if s != requested_seed]
		else:
			print('WARNING: seed %s is NOT a strict-successful seed%s.'
				% (requested_seed, reasons.get(requested_seed, '')))
			print('  Plotting strict-successful seeds instead; this figure '
				'is not meaningful on a seed that died or stalled.')
			seeds_to_plot = ordered
		seeds_to_plot = seeds_to_plot[:N_SEEDS_TO_PLOT]
		print('Strict-successful seeds: %d of %d seeds. Plotting %s.'
			% (len(successful), n_seeds,
				', '.join(str(s) for s in seeds_to_plot)))

		# Burn-in: drop the first IGNORE_FIRST_N_GENS generations if enough
		# remain. Computed against the cohort generation count, not the number
		# of directories returned, so the boundary is a generation index.
		burn_in = sc.IGNORE_FIRST_N_GENS
		if n_generation <= burn_in + 1:
			print('Only %d generations; using all after a reduced burn-in.'
				% n_generation)
			burn_in = max(0, n_generation - 2)

		sim_data = self.read_pickle_file(simDataFile)
		validation_data = self.read_pickle_file(validationDataFile)

		cistron_ids = sim_data.process.transcription.cistron_data['id']
		is_mRNA = sim_data.process.transcription.cistron_data['is_mRNA']
		mRNA_cistron_indexes = np.where(is_mRNA)[0]
		mRNA_cistron_ids = np.array([cistron_ids[x] for x in mRNA_cistron_indexes])

		# Align the gene set with the classification. `is_mRNA` covers every mRNA
		# cistron, while the canonical def-5 gene set is the protein-coding subset
		# (sc.get_mrna_gene_set: cistrons that have an associated monomer), so a
		# handful of cistrons here have no category. Drop them, so this figure and
		# the per-gene tables describe exactly the same genes.
		cat_by_cistron = self._category_by_cistron(cohort_plot_out, clf)
		keep = np.array([cid in cat_by_cistron for cid in mRNA_cistron_ids])
		n_dropped = int((~keep).size - keep.sum())
		if n_dropped:
			print('Dropping %d of %d mRNA cistrons with no def-5 category '
				'(not protein-coding); %d genes remain.'
				% (n_dropped, keep.size, int(keep.sum())))
		if not keep.any():
			print('No mRNA cistron matched the def-5 gene key. Skipping.')
			return
		mRNA_cistron_indexes = mRNA_cistron_indexes[keep]
		mRNA_cistron_ids = mRNA_cistron_ids[keep]
		categories = np.array(
			[cat_by_cistron[cid] for cid in mRNA_cistron_ids], dtype=object)

		gens = np.arange(burn_in, n_generation)
		for seed in seeds_to_plot:
			freqDir = coh_ap.get_cells(
				seed=[seed], generation=gens, only_successful=True)
			if len(freqDir) <= 1:
				print('Skipping seed %d -- %d post-burn-in generations.'
					% (seed, len(freqDir)))
				continue

			# Detect the completed-transcript column.
			first_out = os.path.join(freqDir[0], 'simOut')
			try:
				TableReader(os.path.join(first_out, sc.SYNTH_TABLE)
					).readColumn(sc.SYNTH_COLUMN)
			except Exception:
				print('WARNING: %s/%s not found; cannot run the Definition-5 '
					'variant on this cohort.' % (sc.SYNTH_TABLE, sc.SYNTH_COLUMN))
				return

			self._plot_seed(plotOutDir, plotOutFileName, metadata,
				validation_data, seed, burn_in, freqDir, mRNA_cistron_indexes,
				mRNA_cistron_ids, categories, clf['n_seeds'])

	def _category_by_cistron(self, cohort_plot_out, clf):
		"""{cistron_id: def5_CI category} from the canonical classification.

		The classification is keyed by gene id, so it is joined to cistron ids
		through the raw extraction's gene key, which was written from the same
		ordered gene list as the synth matrix.
		"""
		gene_ids, cistron_ids, _ = sc.load_raw_genes(cohort_plot_out)
		cat_by_gene = dict(zip(clf['gene_ids'], clf['stats']['cat']))
		out = {}
		for gene_id, cistron_id in zip(gene_ids, cistron_ids):
			if gene_id in cat_by_gene:
				out[cistron_id] = cat_by_gene[gene_id]
		return out

	def _requested_seed(self, seedOutDir):
		"""The seed the framework was pointed at (the seed dir's basename)."""
		name = os.path.basename(os.path.normpath(seedOutDir))
		try:
			return int(name)
		except ValueError:
			return -1

	def _seed_success(self, variant_dir, coh_ap, n_generation):
		"""Strict successful-seed set for this cohort.

		Prefers the flags already persisted by subgen_extract.py (cheap);
		falls back to recomputing them, which reads every cell's Main/time.
		Returns (successful_seeds, n_seeds_seen, {seed: rejection reason}).
		"""
		rows = sc.load_seed_success_rows(
			os.path.join(variant_dir, constants.PLOTOUT_DIR))
		if rows:
			print('Using the strict seed flags from the raw extraction.')
			successful = {s for s, r in rows.items()
				if r['is_successful'].strip() == 'True'}
			reasons = {
				s: ' (completed_all_gens=%s, %s/%d gens ran, n_cells_at_180=%s,'
					' gens_at_180=%s)'
					% (r['completed_all_gens'], r['n_gens_ran'], n_generation,
						r['n_cells_at_180'], r['gens_at_180'] or '-')
				for s, r in rows.items() if s not in successful}
			return successful, len(rows), reasons

		print('No raw-extraction seed table found; recomputing the strict '
			'successful-seed set (reads every cell\'s Main/time).')
		_, sim_metadata = sc.load_sim_metadata(variant_dir)
		success = sc.compute_seed_success(coh_ap, n_generation,
			total_init_sims=sim_metadata.get('total_init_sims'))
		reasons = {}
		for s in success['all_seed_ids']:
			if success['in_successful'][s]:
				continue
			missing = sorted(
				set(range(n_generation)) - success['successful_gens'][s])
			reasons[s] = (' (completed_all_gens=%s, %d of %d gens missing,'
				' n_cells_at_180=%d, gens_at_180=%s)'
				% (success['completed_all'][s], len(missing), n_generation,
					success['n_at_180'][s],
					','.join(str(g) for g in success['gens_at_180'][s]) or '-'))
		return (success['successful_seeds'], len(success['all_seed_ids']),
			reasons)

	def _plot_seed(self, plotOutDir, plotOutFileName, metadata,
			validation_data, seed, burn_in, freqDir, mRNA_cistron_indexes,
			mRNA_cistron_ids, categories, n_seeds):
		"""Produce the 5B/5E/5F/5G panels for one strict-successful seed.

		`categories` is the cohort-wide def5_CI label per gene, parallel to
		`mRNA_cistron_ids` -- it does NOT depend on this seed, so every plotted
		seed reports the same category counts by construction.
		"""
		# Seed-suffixed so one seed never overwrites another's panels.
		name = '%s_seed%06d' % (plotOutFileName, seed)

		transcribedBool = []       # per gen: bool, >=1 completed transcript
		synthPerGen = []           # per gen: completed transcript count per gene
		simulatedSynthProbs = []   # for x-axis ordering only
		time = []
		time_eachGen = []
		transcriptionEvents = None  # raster of completed-transcript timesteps

		for gen, simDir in enumerate(freqDir):
			simOutDir = os.path.join(simDir, 'simOut')

			# Read this generation's tables together; if any is missing or
			# unreadable (e.g. a generation that divided but did not write every
			# listener), skip the whole generation so the per-gen arrays stay
			# aligned rather than crashing the plot.
			try:
				synth = TableReader(os.path.join(simOutDir, sc.SYNTH_TABLE)
					).readColumn(sc.SYNTH_COLUMN)
				# Column order here is the full cistron_ids; subset to mRNAs.
				synth = synth[:, mRNA_cistron_indexes]
				synth_prob = TableReader(os.path.join(simOutDir, 'RnaSynthProb')
					).readColumn('actual_rna_synth_prob_per_cistron'
						)[:, mRNA_cistron_indexes]
			except Exception as e:
				print('  Skipping generation %d (%s): %s' % (gen, simDir, e))
				continue

			synthSum = synth.sum(axis=0)
			transcribedBool.append(synthSum > 0)
			synthPerGen.append(synthSum)
			simulatedSynthProbs.append(np.mean(synth_prob, axis=0))

			if gen < RASTER_N_GENS:
				main_reader = TableReader(os.path.join(simOutDir, 'Main'))
				gen_time = main_reader.readColumn('time')
				time += gen_time.tolist()
				time_eachGen.append(gen_time.tolist()[0])
				events = synth != 0
				transcriptionEvents = events if transcriptionEvents is None \
					else np.vstack((transcriptionEvents, events))

		if not synthPerGen:
			print('Skipping seed %d -- no readable generation.' % seed)
			return

		time = np.array(time)
		if time.size:
			time_eachGen.append(time[-1])
		time_eachGen = np.array(time_eachGen)
		transcribedBool = np.array(transcribedBool)
		synthPerGen = np.array(synthPerGen)
		simulatedSynthProbs = np.array(simulatedSynthProbs)
		n_gens_used = synthPerGen.shape[0]
		print('Seed %d: %d post-burn-in generations (burn-in=%d).'
			% (seed, n_gens_used, burn_in))

		# Order genes by mean simulated synthesis probability (x-axis only).
		indexingOrder = np.argsort(np.mean(simulatedSynthProbs, axis=0))
		freqOrdered = np.mean(transcribedBool, axis=0)[indexingOrder]
		def5MeanOrdered = np.mean(synthPerGen, axis=0)[indexingOrder]
		eventsOrdered = transcriptionEvents[:, indexingOrder] \
			if transcriptionEvents is not None else None
		mRNA_ids_ordered = mRNA_cistron_ids[indexingOrder]

		# Canonical def5_CI categories, reordered onto the x-axis order. These are
		# cohort-wide, so they are identical for every plotted seed.
		catOrdered = categories[indexingOrder]
		idx_by_cat = {c: np.where(catOrdered == c)[0] for c in sc.CATEGORIES}
		colors = np.array(
			[sc.PALETTE[c] for c in catOrdered], dtype=object)
		n_genes = len(freqOrdered)
		print('  def5_CI categories (cohort-wide, %d seeds): %s'
			% (n_seeds, ', '.join(
				'%s=%d (%.1f%%)'
				% (c, len(idx_by_cat[c]), 100. * len(idx_by_cat[c]) / n_genes)
				for c in sc.CATEGORIES)))
		# Diagnostic only: what this single seed's point estimate would have
		# said. Kept visible because the gap is the reason the CI form is canonical,
		# but it labels nothing.
		lin_never = int(np.sum(def5MeanOrdered == 0))
		lin_notsub = int(np.sum(def5MeanOrdered >= 1))
		print('  [diagnostic] this seed\'s point estimate would say: '
			'never=%d, 0<mean<1=%d, mean>=1=%d'
			% (lin_never, n_genes - lin_never - lin_notsub, lin_notsub))

		# --- Figure 5B top: frequency scatter + histogram ---
		fig = plt.figure(figsize=(16, 8))
		scatterAxis = plt.subplot2grid((2, 4), (0, 0), colspan=3, rowspan=2)
		histAxis = plt.subplot2grid((2, 4), (0, 3), colspan=1, rowspan=2,
			sharey=scatterAxis)
		scatterAxis.scatter(np.arange(n_genes), freqOrdered, marker='o',
			facecolors=list(colors), edgecolors='none', s=20)
		scatterAxis.set_xlim([0, n_genes])
		scatterAxis.set_ylim([-.01, 1.01])
		whitePadSparklineAxis(scatterAxis)
		histAxis.hist(freqOrdered, bins=n_gens_used + 1, orientation='horizontal',
			color=sc.PALETTE['subgen'])
		histAxis.set_xscale('log')
		whitePadSparklineAxis(histAxis)
		histAxis.xaxis.tick_bottom()
		plt.suptitle('Frequency of >=1 COMPLETED transcript per generation '
			'(Definition 5, def5_CI categories over %d seeds | seed %d '
			'strict-successful seed, burn-in=%d gens, %d gens)'
			% (n_seeds, seed, burn_in, n_gens_used), fontsize=13)
		scatterAxis.set_xlabel(
			'Genes ordered by simulated synthesis probability', fontsize=12)
		scatterAxis.set_ylabel('Fraction of generations', fontsize=12)
		# Category counts stacked down the right-hand side of the histogram, in the
		# canonical order rather than at fixed y positions (four categories no
		# longer fit the old 0 / 0.5 / 1 anchors).
		x_text = histAxis.get_xlim()[1] * 1.6
		for row, c in enumerate(sc.CATEGORIES):
			histAxis.text(x_text, 1.0 - 0.18 * row, '%d %s\n(%0.1f%%)'
				% (len(idx_by_cat[c]), c,
					100. * len(idx_by_cat[c]) / n_genes),
				fontsize=11, va='center', color=sc.PALETTE[c])
		exportFigure(plt, plotOutDir, name + '_5B_top', metadata)
		plt.close('all')

		# --- Figure 5B bottom: completed-transcript event raster ---
		if eventsOrdered is not None and time.size:
			def event_times(indexes):
				out = []
				for i in indexes:
					mask = eventsOrdered[:, i]
					v = (time[mask] / 3600.).tolist()
					out.append(v if mask.sum() else [-1])
				return out
			fig = plt.figure(figsize=(16, 8))
			topAxis = plt.subplot(2, 1, 1)
			botAxis = plt.subplot(2, 1, 2, sharex=topAxis)
			top_cat, bot_cat = RASTER_CATEGORIES
			for axis, cat in ((topAxis, top_cat), (botAxis, bot_cat)):
				idx = idx_by_cat[cat]
				events = event_times(idx)
				if events:
					axis.eventplot(events, orientation='horizontal',
						linewidths=2., linelengths=4., color=sc.PALETTE[cat])
				axis.set_xlim([0, time[-1] / 3600.])
				axis.set_ylim([-1, max(len(idx), 1)])
				axis.set_yticks([])
				axis.set_ylabel('%s\n(n=%d)' % (sc.CAT_LABEL[cat], len(idx)),
					fontsize=12)
			botAxis.set_xlabel('Time (gens)', fontsize=14)
			botAxis.set_xticks(time_eachGen / 3600.)
			botAxis.set_xticklabels(np.arange(len(time_eachGen)))
			plt.suptitle('Completed-transcript events (Definition 5, def5_CI over '
				'%d seeds | seed %d)' % (n_seeds, seed), fontsize=13)
			exportFigure(plt, plotOutDir, name + '_5B_bottom', metadata)
			plt.close('all')

		# --- Figures 5E/5F/5G: gene-category composition by def5_CI category ---
		self._category_bars(plotOutDir, name, metadata, validation_data,
			mRNA_ids_ordered, catOrdered, n_seeds)

	def _category_bars(self, plotOutDir, name, metadata, validation_data,
			mRNA_ids_ordered, catOrdered, n_seeds):
		"""Composition of a gene group across the four def5_CI categories."""
		xloc = np.arange(len(sc.CATEGORIES))
		width = 0.8
		# Index by id so a gene is looked up once instead of by np.where per hit.
		index_by_id = {str(g): i for i, g in enumerate(mRNA_ids_ordered)}

		def new_counts():
			return {c: 0 for c in sc.CATEGORIES}

		def bar(counts, total, ylabel, xlabel, fname):
			if total == 0:
				return
			fig = plt.figure()
			ax = plt.subplot(1, 1, 1)
			ax.bar(xloc + width,
				[counts[c] / float(total) for c in sc.CATEGORIES], width,
				color=[sc.PALETTE[c] for c in sc.CATEGORIES], edgecolor='none')
			whitePadSparklineAxis(ax)
			ax.set_ylabel(ylabel)
			ax.set_xticks(xloc + 1.5 * width)
			ax.set_xticklabels(
				['never' if c == 'never_expressed' else c.replace('_', '-')
					for c in sc.CATEGORIES], fontsize=8)
			ax.set_xlabel('%s (def5_CI, %d seeds)' % (xlabel, n_seeds))
			plt.subplots_adjust(right=0.9, bottom=0.2, left=0.2, top=0.9)
			exportFigure(plt, plotOutDir, fname, metadata)
			plt.close()

		# 5E: essential genes
		essential = validation_data.essential_genes.essential_cistrons
		counts = new_counts()
		n_ess = 0
		for g in essential:
			i = index_by_id.get(str(g))
			if i is None:
				continue
			counts[catOrdered[i]] += 1
			n_ess += 1
		bar(counts, n_ess, 'Fraction of essential genes',
			'Total essential genes: %s' % n_ess,
			name + '_5E')

		# 5F/5G: gene functions (unknown, resistance)
		geneFunctions = validation_data.geneFunctions.geneFunctions
		unknown = new_counts()
		resistance = new_counts()
		for frameID, function_ in geneFunctions.items():
			# geneFunctions is keyed by gene frame id (e.g. EG10001); the mRNA
			# cistron id is that frame id + '_RNA'. Match exactly (like 5E) rather
			# than by substring, which could false-match a longer variable-length
			# id (e.g. G7263 inside G72631_RNA).
			i = index_by_id.get(frameID + '_RNA')
			if i is None:
				continue
			key = catOrdered[i]
			if function_ in ['Unknown function', 'Unclear/under-characterized']:
				unknown[key] += 1
			elif function_ in ['Antibiotic resistance', 'Toxin/antitoxin']:
				resistance[key] += 1
		bar(unknown, sum(unknown.values()),
			'Fraction of poorly understood genes',
			'Total poorly understood genes: %s' % sum(unknown.values()),
			name + '_5F')
		bar(resistance, sum(resistance.values()),
			'Fraction of antibiotic-related genes',
			'Total antibiotic-related genes: %s' % sum(resistance.values()),
			name + '_5G')


if __name__ == '__main__':
	Plot().cli()
