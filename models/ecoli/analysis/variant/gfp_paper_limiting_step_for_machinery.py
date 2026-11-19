"""
Which component limits building more ribosomal subunits and RNA polymerase.

# Figure 5F

Left: share of time 30S and 50S assembly is limited by a protein subunit
rather than an rRNA. Right: share of time each RNAP subunit is the one in
shortest supply. Bars are for the no-GFP control and the designs nearest
each FIG5_PROTEOME_SHARES, means over surviving lineages with t-based 95%
CIs.

At each time step the limiting subunit is the one with the fewest free copies
(BulkMolecules) per copy the complex needs. Ties split the credit equally:
the storyboard version gave every tie to whichever subunit came first in the
list, and because the rRNAs come last that credited protein-rRNA ties
(common when both are at zero) to protein. The share of tied time steps is in
the sidecar.

The supporting claim: under high burden the bottleneck to more machinery is
protein, and rpoC is the RNAP subunit in shortest supply.
"""

import numpy as np
from matplotlib import pyplot as plt

from models.ecoli.analysis import variantAnalysisPlot
from wholecell.analysis import gfp_paper_constants as C
from wholecell.analysis import gfp_paper_data
from wholecell.analysis import gfp_paper_style as style

PERCENT = 100.0
RNAP_SUBUNITS = (('rpoC', 'RPOC-MONOMER[c]'), ('rpoA', 'EG10893-MONOMER[c]'),
	('rpoB', 'RPOB-MONOMER[c]'))


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
		control = gfp_paper_data.CONTROL_VARIANT
		chosen = [gfp_paper_data.nearest_design(data.designs, share)
			for share in C.FIG5_PROTEOME_SHARES]
		variants = [control] + [int(d['variant']) for d in chosen]
		variants = list(dict.fromkeys(variants))       # drop duplicates
		labels = {control: 'no GFP'}
		for d in chosen:
			labels[int(d['variant'])] = 'GFP %.0f%% of proteome' % (
				PERCENT * d['gfp_proteome_share'])
		# Light to dark with burden, skipping the palest step, so two designs
		# read as clearly different greens.
		steps = style.GFP_RAMP[1:4]
		ramp = [steps[int(round(i))] for i in np.linspace(0, len(steps) - 1,
			max(len(variants) - 1, 1))]
		colors = {control: style.COLORS['control']}
		colors.update(dict(zip(variants[1:], ramp)))
		limiting = gfp_paper_data.limiting_profiles(data, variants,
			plotOutDir, self.cpus)

		def protein_share(profile):
			ids = np.array(profile.ids)
			protein = np.array([not gfp_paper_data.is_rrna(i) for i in ids]) \
				& (ids != gfp_paper_data.TIED)
			return [v[protein].sum() for v in profile.per_lineage.values()]

		def subunit_share(profile, monomer):
			i = profile.ids.index(monomer)
			return [v[i] for v in profile.per_lineage.values()]

		def tied(profile):
			i = profile.ids.index(gfp_paper_data.TIED)
			return float(profile.mean[i])

		style.apply()
		fig, (left, right) = plt.subplots(1, 2, figsize=style.panel_size(
			'double', 0.32), gridspec_kw={'width_ratios': [1, 1.1]})
		width = 0.8 / len(variants)
		notes = []
		for k, variant in enumerate(variants):
			offset = (k - (len(variants) - 1) / 2) * width
			for j, complex_label in enumerate(('30S', '50S')):
				values = protein_share(limiting[complex_label][variant])
				left.barh(j + offset, PERCENT * np.mean(values), height=width
					* 0.9, color=colors[variant], xerr=PERCENT * _ci(values),
					error_kw=dict(ecolor=style.INK, elinewidth=0.7,
					capsize=1.2), label=labels[variant] if j == 0 else None)
				notes.append('%s %s protein-limited %.3f (tied %.3f)'
					% (labels[variant], complex_label, np.mean(values),
					tied(limiting[complex_label][variant])))
			for j, (name, monomer) in enumerate(RNAP_SUBUNITS):
				values = subunit_share(limiting['RNAP'][variant], monomer)
				right.bar(j + offset, PERCENT * np.mean(values), width=width
					* 0.9, color=colors[variant], yerr=PERCENT * _ci(values),
					error_kw=dict(ecolor=style.INK, elinewidth=0.7,
					capsize=1.2))
				notes.append('%s %s holds RNAP %.3f' % (labels[variant],
					name, np.mean(values)))
		left.set_yticks([0, 1])
		left.set_yticklabels(['30S', '50S'])
		left.invert_yaxis()
		left.tick_params(axis='y', length=0)
		left.set_xlim(0, 100)
		left.set_xlabel('time limited by protein (%)')
		right.set_xticks(range(len(RNAP_SUBUNITS)))
		right.set_xticklabels([name for name, _ in RNAP_SUBUNITS],
			style='italic')
		right.set_ylim(0, 100)
		right.set_ylabel('time holding the RNAP bottleneck (%)')
		handles, legend_labels = left.get_legend_handles_labels()
		fig.legend(handles, legend_labels, loc='lower center', ncol=len(
			variants), bbox_to_anchor=(0.5, -0.08))

		style.save_panel(fig, plotOutDir, style.panel_name(__file__),
			script=__file__, batch=inputDir, sim_data_path=simDataFile,
			n=int(data.lineages[data.lineages['variant'].isin(variants)
				& data.lineages['survived']].shape[0]),
			notes='Limiting = fewest free copies per copy needed; ties split.'
				' Designs: %s. %s.' % (', '.join('%d = %s' % (v, labels[v])
				for v in variants), '; '.join(notes)),
			**data.provenance(variants))


if __name__ == '__main__':
	Plot().cli()
