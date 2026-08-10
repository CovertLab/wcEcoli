"""
Definition-5 rewrite of subgenerationalTranscription.py (Figure 5B/5E/5F/5G).

Three changes from the original, all to align it with Definition 5:

  1. A gene counts as transcribed in a generation when it produced at least one
     COMPLETED transcript (TranscriptElongationListener/
     countRnaCistronSynthesized > 0), not merely when an mRNA molecule was
     present (which also captures transcripts inherited from the mother and
     initiations later lost to tRNA attenuation).
  2. The first IGNORE_FIRST_N_GENS generations are dropped as burn-in, so the
     per-lineage frequency is not contaminated by initial-condition transients.
     The original averaged over every generation including startup.
  3. The original ran on the seed-0 lineage unconditionally. Seed 0 is not
     necessarily healthy -- a lineage that stalls out (cells pinned at the
     180-minute length cap) still passes ap.get_cells(only_successful=True),
     because those cells did write daughter state, and its collapsing
     transcription inflates the subgenerational fraction. This version instead
     plots the first N_LINEAGES_TO_PLOT STRICT-successful lineages (completed
     every generation, no cell at the doubling cap), the same gate every other
     def-5 analysis uses.

Genes are classified by their per-lineage Definition-5 rate (mean completed
transcripts per generation): never (mean == 0), subgen (0 < mean < 1), or
not_subgen (mean >= 1). One figure set per plotted lineage, seed-suffixed.
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
from models.ecoli.analysis.cohort import subgen_common as sc

# Number of post-burn-in generations to draw in the transcription-event raster.
RASTER_N_GENS = 5
# How many strict-successful lineages to plot (one figure set each). This figure
# is inherently single-lineage, so a few lineages give a sense of the
# lineage-to-lineage spread. Mirrors N_SEEDS_TO_PLOT in
# subgen_monomer_dynamics_def5.py.
N_LINEAGES_TO_PLOT = 3

COLOR_NEVER = 'y'      # never expressed (freq/mean == 0)
COLOR_NOTSUB = 'r'     # not subgenerational (mean >= 1)
COLOR_SUB = 'b'        # subgenerational (0 < mean < 1)


def remove_xaxis(axis):
	axis.spines['bottom'].set_visible(False)
	axis.tick_params(bottom=False, axis='x', labelbottom=False)
	axis.set_xlabel('')


class Plot(multigenAnalysisPlot.MultigenAnalysisPlot):
	def do_plot(self, seedOutDir, plotOutDir, plotOutFileName, simDataFile,
			validationDataFile, metadata):
		# The multigen framework hands us one seed directory, but this figure is
		# only meaningful on a lineage that stayed healthy for the whole run, so
		# the lineage choice is made here rather than by the caller. self.ap is
		# scoped to seedOutDir; a cohort view over the parent variant directory
		# is needed both to reach the other seeds' simOut and to get the
		# cohort-wide generation count -- self.ap.n_generation counts only the
		# generations THIS seed produced, which would make a truncated lineage
		# look like it had completed every generation.
		variant_dir = os.path.dirname(os.path.normpath(seedOutDir))
		requested_seed = self._requested_seed(seedOutDir)
		coh_ap = AnalysisPaths(variant_dir, cohort_plot=True)
		n_generation = coh_ap.n_generation
		if n_generation <= 1:
			print('Skipping -- only runs for multigen')
			return

		successful, n_seeds, reasons = self._lineage_success(
			variant_dir, coh_ap, n_generation)
		if not successful:
			print('WARNING: no strict-successful lineage among %d seeds (needs '
				'all %d generations with no cell at the %g-min doubling cap); '
				'skipping.' % (n_seeds, n_generation, sc.MAX_DOUBLING_MIN))
			return

		ordered = sorted(successful)
		if requested_seed in successful:
			seeds_to_plot = [requested_seed] + [
				s for s in ordered if s != requested_seed]
		else:
			print('WARNING: seed %s is NOT a strict-successful lineage%s.'
				% (requested_seed, reasons.get(requested_seed, '')))
			print('  Plotting strict-successful lineages instead; this figure '
				'is not meaningful on a lineage that died or stalled.')
			seeds_to_plot = ordered
		seeds_to_plot = seeds_to_plot[:N_LINEAGES_TO_PLOT]
		print('Strict-successful lineages: %d of %d seeds. Plotting %s.'
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

			self._plot_lineage(plotOutDir, plotOutFileName, metadata,
				validation_data, seed, burn_in, freqDir, mRNA_cistron_indexes,
				mRNA_cistron_ids)

	def _requested_seed(self, seedOutDir):
		"""The seed the framework was pointed at (the seed dir's basename)."""
		name = os.path.basename(os.path.normpath(seedOutDir))
		try:
			return int(name)
		except ValueError:
			return -1

	def _lineage_success(self, variant_dir, coh_ap, n_generation):
		"""Strict successful-lineage set for this cohort.

		Prefers the flags already persisted by subgen_raw_extract.py (cheap);
		falls back to recomputing them, which reads every cell's Main/time.
		Returns (successful_seeds, n_seeds_seen, {seed: rejection reason}).
		"""
		rows = sc.load_lineage_success_rows(
			os.path.join(variant_dir, constants.PLOTOUT_DIR))
		if rows:
			print('Using the strict lineage flags from the raw extraction.')
			successful = {s for s, r in rows.items()
				if r['is_successful'].strip() == 'True'}
			reasons = {
				s: ' (completed_all_gens=%s, %s/%d gens ran, n_cells_at_180=%s,'
					' gens_at_180=%s)'
					% (r['completed_all_gens'], r['n_gens_ran'], n_generation,
						r['n_cells_at_180'], r['gens_at_180'] or '-')
				for s, r in rows.items() if s not in successful}
			return successful, len(rows), reasons

		print('No raw-extraction lineage table found; recomputing the strict '
			'successful-lineage set (reads every cell\'s Main/time).')
		_, sim_metadata = sc.load_sim_metadata(variant_dir)
		success = sc.compute_lineage_success(coh_ap, n_generation,
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

	def _plot_lineage(self, plotOutDir, plotOutFileName, metadata,
			validation_data, seed, burn_in, freqDir, mRNA_cistron_indexes,
			mRNA_cistron_ids):
		"""Produce the 5B/5E/5F/5G panels for one strict-successful lineage."""
		# Seed-suffixed so one lineage never overwrites another's panels.
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

		# Definition-5 categories on this lineage.
		neverIdx = np.where(def5MeanOrdered == 0)[0]
		notSubIdx = np.where(def5MeanOrdered >= 1)[0]
		subIdx = np.array([i for i in np.arange(len(def5MeanOrdered))
			if i not in set(neverIdx.tolist()) | set(notSubIdx.tolist())],
			dtype=int)
		colors = np.repeat(COLOR_SUB, len(freqOrdered))
		colors[neverIdx] = COLOR_NEVER
		colors[notSubIdx] = COLOR_NOTSUB
		print('  never=%d, subgen=%d (%.1f%%), not_subgen=%d'
			% (len(neverIdx), len(subIdx),
				100. * len(subIdx) / len(freqOrdered), len(notSubIdx)))

		# --- Figure 5B top: frequency scatter + histogram ---
		fig = plt.figure(figsize=(16, 8))
		scatterAxis = plt.subplot2grid((2, 4), (0, 0), colspan=3, rowspan=2)
		histAxis = plt.subplot2grid((2, 4), (0, 3), colspan=1, rowspan=2,
			sharey=scatterAxis)
		scatterAxis.scatter(np.arange(len(freqOrdered)), freqOrdered, marker='o',
			facecolors=colors, edgecolors='none', s=20)
		scatterAxis.set_xlim([0, len(freqOrdered)])
		scatterAxis.set_ylim([-.01, 1.01])
		whitePadSparklineAxis(scatterAxis)
		histAxis.hist(freqOrdered, bins=n_gens_used + 1, orientation='horizontal',
			color=COLOR_SUB)
		histAxis.set_xscale('log')
		whitePadSparklineAxis(histAxis)
		histAxis.xaxis.tick_bottom()
		plt.suptitle('Frequency of >=1 COMPLETED transcript per generation '
			'(Definition 5, seed %d strict-successful lineage, burn-in=%d gens, '
			'%d gens)' % (seed, burn_in, n_gens_used), fontsize=14)
		scatterAxis.set_xlabel(
			'Genes ordered by simulated synthesis probability', fontsize=12)
		scatterAxis.set_ylabel('Fraction of generations', fontsize=12)
		histAxis.text(histAxis.get_xlim()[1] * 1.6, 0, '%s never\n(%0.1f%%)'
			% (len(neverIdx), 100. * len(neverIdx) / len(freqOrdered)),
			fontsize=12, va='center', color=COLOR_NEVER)
		histAxis.text(histAxis.get_xlim()[1] * 1.6, 1, '%s not-subgen\n(%0.1f%%)'
			% (len(notSubIdx), 100. * len(notSubIdx) / len(freqOrdered)),
			fontsize=12, va='center', color=COLOR_NOTSUB)
		histAxis.text(histAxis.get_xlim()[1] * 1.6, 0.5, '%s subgen\n(%0.1f%%)'
			% (len(subIdx), 100. * len(subIdx) / len(freqOrdered)),
			fontsize=12, va='center', color=COLOR_SUB)
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
			notSubEvents = event_times(notSubIdx)
			subEvents = event_times(subIdx)
			fig = plt.figure(figsize=(16, 8))
			topAxis = plt.subplot(2, 1, 1)
			botAxis = plt.subplot(2, 1, 2, sharex=topAxis)
			if notSubEvents:
				topAxis.eventplot(notSubEvents, orientation='horizontal',
					linewidths=2., linelengths=4., color=COLOR_NOTSUB)
			topAxis.set_xlim([0, time[-1] / 3600.])
			topAxis.set_ylim([-1, max(len(notSubIdx), 1)])
			topAxis.set_yticks([])
			topAxis.set_ylabel('mean >= 1', fontsize=14)
			if subEvents:
				botAxis.eventplot(subEvents, orientation='horizontal',
					linewidths=2., linelengths=4., color=COLOR_SUB)
			botAxis.set_xlim([0, time[-1] / 3600.])
			botAxis.set_ylim([-1, max(len(subIdx), 1)])
			botAxis.set_yticks([])
			botAxis.set_ylabel('0 < mean < 1', fontsize=14)
			botAxis.set_xlabel('Time (gens)', fontsize=14)
			botAxis.set_xticks(time_eachGen / 3600.)
			botAxis.set_xticklabels(np.arange(len(time_eachGen)))
			plt.suptitle('Completed-transcript events (Definition 5, seed %d)'
				% seed, fontsize=14)
			exportFigure(plt, plotOutDir, name + '_5B_bottom', metadata)
			plt.close('all')

		# --- Figures 5E/5F/5G: gene-category composition by freq group ---
		self._category_bars(plotOutDir, name, metadata, validation_data,
			mRNA_ids_ordered, def5MeanOrdered)

	def _classify_mean(self, m):
		if m == 0:
			return 'r'  # never
		elif m >= 1:
			return 'b'  # not subgen
		return 'g'      # subgen

	def _category_bars(self, plotOutDir, name, metadata, validation_data,
			mRNA_ids_ordered, def5MeanOrdered):
		xloc = np.arange(3)
		width = 0.8
		id_set = set(mRNA_ids_ordered)

		def bar(counts, total, ylabel, xlabel, fname):
			if total == 0:
				return
			fig = plt.figure()
			ax = plt.subplot(1, 1, 1)
			ax.bar(xloc + width,
				[counts['r'] / float(total), counts['g'] / float(total),
					counts['b'] / float(total)], width,
				color=[COLOR_NEVER, COLOR_SUB, COLOR_NOTSUB], edgecolor='none')
			whitePadSparklineAxis(ax)
			ax.set_ylabel(ylabel)
			ax.set_xticks(xloc + 1.5 * width)
			ax.set_xticklabels(['never', 'subgen', 'not-subgen'])
			ax.set_xlabel(xlabel)
			plt.subplots_adjust(right=0.9, bottom=0.15, left=0.2, top=0.9)
			exportFigure(plt, plotOutDir, fname, metadata)
			plt.close()

		# 5E: essential genes
		essential = validation_data.essential_genes.essential_cistrons
		counts = {'r': 0, 'g': 0, 'b': 0}
		n_ess = 0
		for g in essential:
			if str(g) not in id_set:
				continue
			i = np.where(mRNA_ids_ordered == str(g))[0][0]
			counts[self._classify_mean(def5MeanOrdered[i])] += 1
			n_ess += 1
		bar(counts, n_ess, 'Percentage of essential genes',
			'Total essential genes: %s' % n_ess,
			name + '_5E')

		# 5F/5G: gene functions (unknown, resistance)
		geneFunctions = validation_data.geneFunctions.geneFunctions
		unknown = {'r': 0, 'g': 0, 'b': 0}
		resistance = {'r': 0, 'g': 0, 'b': 0}
		for frameID, function_ in geneFunctions.items():
			# geneFunctions is keyed by gene frame id (e.g. EG10001); the mRNA
			# cistron id is that frame id + '_RNA'. Match exactly (like 5E) rather
			# than by substring, which could false-match a longer variable-length
			# id (e.g. G7263 inside G72631_RNA).
			matches = np.where(mRNA_ids_ordered == frameID + '_RNA')[0]
			if matches.size == 0:
				continue
			i = matches[0]
			key = self._classify_mean(def5MeanOrdered[i])
			if function_ in ['Unknown function', 'Unclear/under-characterized']:
				unknown[key] += 1
			elif function_ in ['Antibiotic resistance', 'Toxin/antitoxin']:
				resistance[key] += 1
		bar(unknown, sum(unknown.values()),
			'Percentage of poorly understood genes',
			'Total poorly understood genes: %s' % sum(unknown.values()),
			name + '_5F')
		bar(resistance, sum(resistance.values()),
			'Percentage of antibiotic-related genes',
			'Total antibiotic-related genes: %s' % sum(resistance.values()),
			name + '_5G')


if __name__ == '__main__':
	Plot().cli()
