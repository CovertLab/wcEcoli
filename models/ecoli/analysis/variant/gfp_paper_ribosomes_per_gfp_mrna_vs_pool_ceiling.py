"""
Ribosomes per gfp mRNA against gfp mRNA copies, beside the two ceilings that
could stop it rising: the mRNA filling up, and the cell running out of
ribosomes.

# Figure 3F

Each point is one design (promoter x RBS), averaged over the analysis window
of its surviving lineages; one line per RBS strength, running through the
promoter strengths. Hollow points are designs where some lineages failed.

- "mRNA full of ribosomes" is the gfp transcript length over the ribosome
	footprint, both read from sim_data: the most ribosomes one mRNA can hold.
	It is a physical limit only once the footprint is applied per nucleotide
	(the footprint units fix); in batches run before that fix the model's own
	cap sat about three times lower.
- The pool ceiling is every ribosome in the cell (active plus free 30S/50S
	pairs) divided by the gfp mRNA count: the ribosomes per mRNA there would
	be if gfp held all of them. It is drawn for the strongest RBS, the designs
	closest to it.

The claim this panel supports is that the strongest designs meet the pool
ceiling before the mRNA fills, so space on the transcript is not what limits
translation.
"""

import numpy as np
from matplotlib import pyplot as plt

from models.ecoli.analysis import variantAnalysisPlot
from wholecell.analysis import gfp_paper_constants as C
from wholecell.analysis import gfp_paper_data
from wholecell.analysis import gfp_paper_style as style


class Plot(variantAnalysisPlot.VariantAnalysisPlot):
	def do_plot(self, inputDir, plotOutDir, plotOutFileName, simDataFile,
			validationDataFile, metadata):
		data = gfp_paper_data.load(inputDir, simDataFile, metadata,
			plotOutDir, self.cpus)
		designs = data.designs[data.designs['promoter'] > 0]
		gene = data.gene

		style.apply()
		fig, ax = plt.subplots(figsize=style.panel_size('single', 0.62))
		shown = []
		for rbs in C.RBS_LINES_RIBOSOME:
			line = designs[np.isclose(designs['rbs'], rbs)].sort_values(
				'promoter')
			line = line[np.isfinite(line['ribosomes_per_gfp_mrna'])
				& (line['gfp_mrna'] > 0)]
			if line.empty:
				continue
			color = style.rbs_color(rbs)
			ax.plot(line['gfp_mrna'], line['ribosomes_per_gfp_mrna'],
				color=color, zorder=2)
			style.design_scatter(ax, line['gfp_mrna'],
				line['ribosomes_per_gfp_mrna'], color, failed=line['failed'],
				label=C.rbs_label(rbs), zorder=3)
			shown.append(line)

		# The pool ceiling for the strongest RBS shown
		ceiling_rbs = max(C.RBS_LINES_RIBOSOME)
		ceiling = designs[np.isclose(designs['rbs'], ceiling_rbs)].sort_values(
			'promoter')
		ceiling = ceiling[np.isfinite(ceiling['ribosome_ceiling_per_gfp_mrna'])
			& (ceiling['gfp_mrna'] > 0)]
		if not ceiling.empty:
			ax.plot(ceiling['gfp_mrna'], ceiling['ribosome_ceiling_per_gfp_mrna'],
				color=style.COLORS['control'], linestyle=(0, (1.5, 1.5)),
				linewidth=1.0, zorder=1)
			ax.scatter(ceiling['gfp_mrna'],
				ceiling['ribosome_ceiling_per_gfp_mrna'], s=14,
				facecolors='white', edgecolors=style.COLORS['control'],
				linewidths=0.9, zorder=2,
				label='all ribosomes in the cell, %s designs'
					% C.rbs_label(ceiling_rbs))

		style.reference_line(ax, gene.ribosomes_when_full,
			'mRNA full of ribosomes\n(%g nt footprint)'
			% gene.ribosome_footprint_nt)

		# Storyboard ranges, widened only if the data fall outside them
		plotted = shown + [ceiling]
		style.log_axis(ax, 'x', 30, 3000,
			data=np.concatenate([d['gfp_mrna'].values for d in plotted]))
		style.log_axis(ax, 'y', 0.5, 400, data=np.concatenate(
			[d['ribosomes_per_gfp_mrna'].values for d in shown]
			+ [ceiling['ribosome_ceiling_per_gfp_mrna'].values]))
		ax.set_xlabel('gfp mRNA per cell (copies)')
		ax.set_ylabel('ribosomes per gfp mRNA')
		handles, labels = ax.get_legend_handles_labels()
		if any(line['failed'].any() for line in shown):
			handles.append(style.failed_legend_handle())
			labels.append(style.FAILED_LABEL)
		ax.legend(handles, labels, loc='upper center',
			bbox_to_anchor=(0.5, -0.22), ncol=2, handletextpad=0.4,
			columnspacing=1.2)

		variants = sorted(int(v) for v in designs['variant'])
		style.save_panel(fig, plotOutDir, style.panel_name(__file__),
			script=__file__, batch=inputDir, sim_data_path=simDataFile,
			n=int(sum(len(line) for line in shown)),
			notes='One point per design: window means over surviving'
				' lineages, ratio of mean ribosomes on gfp to mean physical'
				' gfp mRNA. mRNA-full line = %g nt / %g nt. Designs with no'
				' surviving lineage are omitted: %s.'
				% (gene.length_nt, gene.ribosome_footprint_nt,
				gfp_paper_data.lost_designs(designs)),
			**data.provenance(variants))


if __name__ == '__main__':
	Plot().cli()
