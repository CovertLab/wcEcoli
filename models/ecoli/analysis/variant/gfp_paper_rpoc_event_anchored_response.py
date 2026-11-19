"""
rpoC's response to its own replication, anchored on each initiation event.

# Figure 5D

Around every replication initiation in the analysis window of surviving
high-GFP lineages, each quantity's pooled rate in a post window is divided by
its pooled rate in the pre window (EVENT_PRE_WINDOW_MIN,
EVENT_POST_WINDOWS_MIN), so 1 means no change:
	rpoC gene copies         RnaSynthProb gene_copy_number
	rpoC transcription       RnapData rna_init_event_per_cistron
	rpoC translation         RibosomeData ribosome_init_event_per_monomer
	terminus transcription   RnapData rnaInitEvent summed over the terminus
	                         reference units, the control: its copies do not
	                         change in these minutes

Every line uses the same estimator, summed counts over summed time steps
across all events (gfp_paper_data.pooled_response), because rpoC sees well
under one transcription event per window, so per-event or per-lineage ratios
are mostly zeros. Shaded bands are 95% bootstrap CIs over lineages. The
storyboard version mixed pooled and per-lineage estimators and showed
translation as a share of all initiations; this one uses raw rates for all.

The supporting claim: gene dosage drives machinery expression, since rpoC
transcription follows its copies up while terminus transcription does not.
"""

import pickle

import numpy as np
from matplotlib import pyplot as plt

from models.ecoli.analysis import variantAnalysisPlot
from wholecell.analysis import gfp_paper_constants as C
from wholecell.analysis import gfp_paper_data
from wholecell.analysis import gfp_paper_style as style

LINES = (
	('rpoc_copies', '$\\it{rpoC}$ gene copies', '#7FB8E0'),
	('rpoc_transcription', '$\\it{rpoC}$ transcription', '#0072B2'),
	('rpoc_translation', '$\\it{rpoC}$ translation', '#003D63'),
	('terminus_transcription', 'terminus transcription', None),
)


class Plot(variantAnalysisPlot.VariantAnalysisPlot):
	def do_plot(self, inputDir, plotOutDir, plotOutFileName, simDataFile,
			validationDataFile, metadata):
		data = gfp_paper_data.load(inputDir, simDataFile, metadata,
			plotOutDir, self.cpus)
		with open(simDataFile, 'rb') as f:
			sim_data = pickle.load(f)
		high = gfp_paper_data.nearest_design(data.designs,
			C.FIG5_PROTEOME_SHARES[-1])
		variant = int(high['variant'])
		events = gfp_paper_data.initiation_events(data, sim_data, [variant],
			plotOutDir, self.cpus)
		if events.empty:
			print('No complete initiation events; skipping panel 5D.')
			return
		pre_mid = sum(C.EVENT_PRE_WINDOW_MIN) / 2

		style.apply()
		fig, ax = plt.subplots(figsize=style.panel_size('single', 0.62))
		notes = []
		for column, label, color in LINES:
			color = color or style.COLORS['control']
			r = gfp_paper_data.pooled_response(events, column)
			x = np.concatenate([[pre_mid], r['mid']])
			y = np.concatenate([[1.0], r['ratio']])
			ax.fill_between(r['mid'], r['lo'], r['hi'], color=color,
				alpha=0.18, linewidth=0, zorder=1)
			ax.plot(x, y, color=color, zorder=2, label=label)
			ax.scatter(x, y, s=12, color=color, edgecolors='white',
				linewidths=0.5, zorder=3)
			notes.append('%s %s' % (column, ', '.join('%.3f [%.3f, %.3f]'
				% (a, b, c) for a, b, c in zip(r['ratio'], r['lo'],
				r['hi']))))
		ax.axhline(1.0, color=style.INK_MUTED, linestyle=(0, (4, 3)),
			linewidth=0.8, zorder=0)
		ax.axvline(0.0, color=style.INK_MUTED, linestyle=(0, (4, 3)),
			linewidth=0.8, zorder=0)
		ax.set_xlabel('minutes after replication initiation')
		ax.set_ylabel('relative to before initiation')
		ax.legend(loc='upper left', bbox_to_anchor=(1.0, 0.9))
		n_events = int(events[events['window'] == 'pre'][['seed',
			'generation', 'event']].drop_duplicates().shape[0])

		style.save_panel(fig, plotOutDir, style.panel_name(__file__),
			script=__file__, batch=inputDir, sim_data_path=simDataFile,
			n=n_events,
			notes='Variant %d (GFP %.3f of proteome). %d events from %d'
				' lineages, generations %d-%d. Pre window %s min, post %s.'
				' Pooled ratio [95%% bootstrap CI over lineages, %d'
				' resamples] per post window: %s.' % (variant,
				high['gfp_proteome_share'], n_events,
				events['seed'].nunique(), C.ANALYSIS_FIRST_GEN,
				C.ANALYSIS_LAST_GEN, C.EVENT_PRE_WINDOW_MIN,
				C.EVENT_POST_WINDOWS_MIN, C.EVENT_BOOTSTRAP,
				'; '.join(notes)),
			**data.provenance([variant]))


if __name__ == '__main__':
	Plot().cli()
