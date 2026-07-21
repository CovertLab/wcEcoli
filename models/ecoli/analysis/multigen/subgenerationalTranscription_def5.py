"""
Definition-5 rewrite of subgenerationalTranscription.py (Figure 5B/5E/5F/5G).

Two changes from the original, both to align it with Definition 5:

  1. A gene counts as transcribed in a generation when it produced at least one
     COMPLETED transcript (TranscriptElongationListener/
     countRnaCistronSynthesized > 0), not merely when an mRNA molecule was
     present (which also captures transcripts inherited from the mother and
     initiations later lost to tRNA attenuation).
  2. The first IGNORE_FIRST_N_GENS generations are dropped as burn-in, so the
     per-lineage frequency is not contaminated by initial-condition transients.
     The original averaged over every generation including startup.

Genes are classified by their per-lineage Definition-5 rate (mean completed
transcripts per generation): never (mean == 0), subgen (0 < mean < 1), or
not_subgen (mean >= 1). Runs on the seed-0 lineage, like the original.
"""

import os

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from wholecell.io.tablereader import TableReader
from wholecell.utils.sparkline import whitePadSparklineAxis
from wholecell.analysis.analysis_tools import exportFigure
from models.ecoli.analysis import multigenAnalysisPlot
from models.ecoli.analysis.cohort import subgen_common as sc

# Number of post-burn-in generations to draw in the transcription-event raster.
RASTER_N_GENS = 5

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
		if 0 not in self.ap._path_data['seed']:
			print('Skipping -- only runs for seed 0')
			return
		# only_successful=True drops generations that did not finish dividing
		# (and therefore may be missing listener tables), so we never try to read
		# a simOut dir for an incomplete generation.
		allDir = self.ap.get_cells(seed=[0], only_successful=True)
		if len(allDir) <= 1:
			print('Skipping -- only runs for multigen')
			return

		# Burn-in: drop the first IGNORE_FIRST_N_GENS generations if enough remain.
		burn_in = sc.IGNORE_FIRST_N_GENS
		if len(allDir) <= burn_in + 1:
			print('Only %d generations; using all after a reduced burn-in.'
				% len(allDir))
			burn_in = max(0, len(allDir) - 2)
		freqDir = allDir[burn_in:]

		sim_data = self.read_pickle_file(simDataFile)
		validation_data = self.read_pickle_file(validationDataFile)

		cistron_ids = sim_data.process.transcription.cistron_data['id']
		is_mRNA = sim_data.process.transcription.cistron_data['is_mRNA']
		mRNA_cistron_indexes = np.where(is_mRNA)[0]
		mRNA_cistron_ids = np.array([cistron_ids[x] for x in mRNA_cistron_indexes])

		# Detect the completed-transcript column.
		first_out = os.path.join(freqDir[0], 'simOut')
		try:
			TableReader(os.path.join(first_out, sc.SYNTH_TABLE)
				).readColumn(sc.SYNTH_COLUMN)
		except Exception:
			print('WARNING: %s/%s not found; cannot run the Definition-5 '
				'variant on this cohort.' % (sc.SYNTH_TABLE, sc.SYNTH_COLUMN))
			return

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

		time = np.array(time)
		if time.size:
			time_eachGen.append(time[-1])
		time_eachGen = np.array(time_eachGen)
		transcribedBool = np.array(transcribedBool)
		synthPerGen = np.array(synthPerGen)
		simulatedSynthProbs = np.array(simulatedSynthProbs)

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
		histAxis.hist(freqOrdered, bins=len(freqDir) + 1, orientation='horizontal',
			color=COLOR_SUB)
		histAxis.set_xscale('log')
		whitePadSparklineAxis(histAxis)
		histAxis.xaxis.tick_bottom()
		plt.suptitle('Frequency of >=1 COMPLETED transcript per generation '
			'(Definition 5, burn-in=%d gens)' % burn_in, fontsize=14)
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
		exportFigure(plt, plotOutDir, plotOutFileName + '_5B_top', metadata)
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
			plt.suptitle('Completed-transcript events (Definition 5)', fontsize=14)
			exportFigure(plt, plotOutDir, plotOutFileName + '_5B_bottom', metadata)
			plt.close('all')

		# --- Figures 5E/5F/5G: gene-category composition by freq group ---
		self._category_bars(plotOutDir, plotOutFileName, metadata,
			validation_data, mRNA_ids_ordered, def5MeanOrdered)

	def _classify_mean(self, m):
		if m == 0:
			return 'r'  # never
		elif m >= 1:
			return 'b'  # not subgen
		return 'g'      # subgen

	def _category_bars(self, plotOutDir, plotOutFileName, metadata,
			validation_data, mRNA_ids_ordered, def5MeanOrdered):
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
			plotOutFileName + '_5E')

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
			plotOutFileName + '_5F')
		bar(resistance, sum(resistance.values()),
			'Percentage of antibiotic-related genes',
			'Total antibiotic-related genes: %s' % sum(resistance.values()),
			plotOutFileName + '_5G')


if __name__ == '__main__':
	Plot().cli()
