"""
Where the cell's active RNA polymerases and ribosomes are working, without
GFP and with high GFP.

# Figure 2G

Left: share of active RNAP on each gene group, from transcripts still being
made (RNACounts rnap_engaged_mRNA_counts, or partial_mRNA_counts in older
batches) and RNAP on rRNA operons (partial_rRNA_counts), over the active RNAP
count. Right: share of active ribosomes making each group's proteins
(RibosomeData n_ribosomes_per_transcript, which sums to the active ribosome
count) over active ribosomes. Bars are means over surviving lineages of each
lineage's window-mean share, with t-based 95% CIs; the value is written above
each bar.

Groups, from sim_data flags: GFP (the inserted gene), rRNA (rRNA operons;
RNAP only), ribosomal proteins and RNAP subunits (transcription units or
monomers flagged as such). "High GFP" is the design nearest the second entry
of FIG5_PROTEOME_SHARES.

The supporting claim: the RNAP and ribosomes the cell has are diverted to GFP
and away from making more RNAP and ribosomes.
"""

import pickle

import numpy as np
from matplotlib import pyplot as plt

from models.ecoli.analysis import variantAnalysisPlot
from wholecell.analysis import gfp_paper_constants as C
from wholecell.analysis import gfp_paper_data
from wholecell.analysis import gfp_paper_style as style

PERCENT = 100.0


def _ci(values):
	from scipy import stats
	values = np.asarray(values, dtype=float)
	n = values.size
	return (float(stats.t.ppf(0.975, n - 1) * values.std(ddof=1)
		/ np.sqrt(n)) if n > 1 else np.nan)


class Plot(variantAnalysisPlot.VariantAnalysisPlot):
	def do_plot(self, inputDir, plotOutDir, plotOutFileName, simDataFile,
			validationDataFile, metadata):
		data = gfp_paper_data.load(inputDir, simDataFile, metadata,
			plotOutDir, self.cpus)
		with open(simDataFile, 'rb') as f:
			sim_data = pickle.load(f)
		high = gfp_paper_data.nearest_design(data.designs,
			C.FIG5_PROTEOME_SHARES[-1])
		control = gfp_paper_data.CONTROL_VARIANT
		variants = [control, int(high['variant'])]
		labels = {control: 'no GFP', int(high['variant']):
			'GFP %.0f%% of proteome' % (PERCENT * high['gfp_proteome_share'])}
		colors = {control: style.COLORS['control'],
			int(high['variant']): style.GFP_RAMP[3]}

		groups = self._groups(sim_data, data.gene)
		mrna = gfp_paper_data.vector_profiles(data, 'rnap_mrna_tu', variants,
			plotOutDir, self.cpus)
		rrna = gfp_paper_data.vector_profiles(data, 'rnap_rrna_tu', variants,
			plotOutDir, self.cpus)
		ribosomes = gfp_paper_data.vector_profiles(data, 'ribosomes_monomer',
			variants, plotOutDir, self.cpus)
		active = data.lineages.set_index(['variant', 'seed'])

		def shares(profile, ids, denominator_column, variant):
			cols = [profile.ids.index(i) for i in ids if i in profile.ids]
			out = []
			for seed, values in profile.per_lineage.items():
				denominator = active.loc[(variant, seed), denominator_column]
				if denominator > 0:
					out.append(values[cols].sum() / denominator)
			return np.array(out)

		style.apply()
		fig, (left, right) = plt.subplots(1, 2, figsize=style.panel_size(
			'double', 0.36), gridspec_kw={'width_ratios': [4, 3]})
		notes = []
		for ax, machine, bars in (
				(left, 'RNAP', [('GFP', 'mrna', groups['gfp_tu']),
					('rRNA', 'rrna', None),
					('ribosomal\nproteins', 'mrna', groups['rprotein_tu']),
					('RNAP\nsubunits', 'mrna', groups['rnap_tu'])]),
				(right, 'ribosome', [('GFP', 'mono', groups['gfp_monomer']),
					('ribosomal\nproteins', 'mono', groups['rprotein_monomer']),
					('RNAP\nsubunits', 'mono', groups['rnap_monomer'])])):
			width = 0.36
			for k, variant in enumerate(variants):
				means, cis = [], []
				for label, source, ids in bars:
					if source == 'rrna':
						profile = rrna[variant]
						ids = profile.ids
					else:
						profile = mrna[variant] if source == 'mrna' else (
							ribosomes[variant])
					values = shares(profile, ids, 'active_rnap'
						if machine == 'RNAP' else 'active_ribosomes', variant)
					means.append(PERCENT * values.mean() if values.size
						else np.nan)
					cis.append(PERCENT * _ci(values))
					notes.append('%s %s %s %.2f%%' % (labels[variant],
						machine, label.replace('\n', ' '), means[-1]))
				x = np.arange(len(bars)) + (k - 0.5) * width
				ax.bar(x, means, width=width * 0.92, color=colors[variant],
					yerr=np.where(np.isfinite(cis), cis, 0),
					error_kw=dict(ecolor=style.INK, elinewidth=0.6,
					capsize=1.2), label=labels[variant], zorder=2)
				for xi, m, c in zip(x, means, np.nan_to_num(cis)):
					if np.isfinite(m):
						ax.annotate('%.2g' % m if m < 1 else '%.0f' % m,
							xy=(xi, m + c), xytext=(0, 2),
							textcoords='offset points', ha='center',
							va='bottom', fontsize=style.FONT_SIZES[
							'annotation'], color=style.INK)
			ax.set_xticks(range(len(bars)))
			ax.set_xticklabels([label for label, _, _ in bars])
			ax.set_ylim(bottom=0)
			ax.set_ylabel('share of active %s (%%)' % ('RNA polymerases'
				if machine == 'RNAP' else 'ribosomes'))
		handles, legend_labels = left.get_legend_handles_labels()
		style.legend_below(fig, handles, legend_labels, ncol=2)

		style.save_panel(fig, plotOutDir, style.panel_name(__file__),
			script=__file__, batch=inputDir, sim_data_path=simDataFile,
			n=int(min(len(p.per_lineage) for p in ribosomes.values())),
			notes='Window-mean shares per surviving lineage, then mean with'
				' t-based 95%% CI. High GFP = variant %d. %s.'
				% (high['variant'], '; '.join(notes)),
			**data.provenance(variants))

	@staticmethod
	def _groups(sim_data, gene):
		"""Ids per gene group, in each listener's id space."""
		transcription = sim_data.process.transcription
		rna = transcription.rna_data
		cistrons = transcription.cistron_data
		monomers = sim_data.process.translation.monomer_data
		cistron_flags = dict(zip(cistrons['id'], zip(
			cistrons['is_ribosomal_protein'], cistrons['is_RNAP'])))
		monomer_ids = list(monomers['id'])
		monomer_cistron = list(monomers['cistron_id'])
		return dict(
			gfp_tu=[gene.rna_id],
			rprotein_tu=list(rna['id'][rna['includes_ribosomal_protein']]),
			rnap_tu=list(rna['id'][rna['includes_RNAP']]),
			gfp_monomer=[gene.monomer_id],
			rprotein_monomer=[m for m, c in zip(monomer_ids, monomer_cistron)
				if cistron_flags.get(c, (False, False))[0]],
			rnap_monomer=[m for m, c in zip(monomer_ids, monomer_cistron)
				if cistron_flags.get(c, (False, False))[1]],
		)


if __name__ == '__main__':
	Plot().cli()
