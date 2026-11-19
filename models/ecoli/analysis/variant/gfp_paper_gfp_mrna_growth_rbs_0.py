"""
gfp mRNA made and growth rate, for designs that transcribe gfp but never
translate it (RBS 0), across promoter strengths.

# Figure 3C

Split into two files, as the storyboard review asked:
	..._mrna     gfp mRNA copies per cell (physical count, window mean)
	..._growth   growth rate relative to the seed-matched no-GFP control,
				 mean over surviving lineages with a 95% CI

The supporting claim: promoter strength raises gfp mRNA up to the RNAP
spacing cap, and transcription alone carries no growth cost however much mRNA
it makes. The dashed line in the growth half is the control (1).

Bars use the pale GFP step, the colour reserved for RBS 0 across the paper.
"""

import numpy as np
from matplotlib import pyplot as plt

from models.ecoli.analysis import variantAnalysisPlot
from wholecell.analysis import gfp_paper_data
from wholecell.analysis import gfp_paper_style as style

RBS = 0
BAR_WIDTH = 0.22               # in expression-factor units on the x axis


class Plot(variantAnalysisPlot.VariantAnalysisPlot):
	def do_plot(self, inputDir, plotOutDir, plotOutFileName, simDataFile,
			validationDataFile, metadata):
		data = gfp_paper_data.load(inputDir, simDataFile, metadata,
			plotOutDir, self.cpus)
		designs = data.designs[(data.designs['promoter'] > 0)
			& np.isclose(data.designs['rbs'], RBS)].sort_values('promoter')
		if designs.empty:
			print('No RBS %g designs in this batch; skipping panel 3C.' % RBS)
			return
		variants = sorted(int(v) for v in designs['variant'])
		provenance = data.provenance(variants)
		lost = gfp_paper_data.lost_designs(designs)
		style.apply()

		# gfp mRNA
		fig, ax = plt.subplots(figsize=style.panel_size('single', 1.1, 0.5))
		self._bars(ax, designs, designs['gfp_mrna'])
		style.promoter_axis(ax, designs['promoter'])
		ax.set_ylim(bottom=0)
		ax.set_ylabel('gfp mRNA (copies)')
		style.save_panel(fig, plotOutDir, style.panel_name(__file__, 'mrna'),
			script=__file__, batch=inputDir, sim_data_path=simDataFile,
			n=len(variants),
			notes='RBS %g designs. Physical gfp mRNA, window means over'
				' surviving lineages. Omitted (no survivors): %s.'
				% (RBS, lost), **provenance)

		# growth relative to control
		fig, ax = plt.subplots(figsize=style.panel_size('single', 1.1, 0.5))
		self._bars(ax, designs, designs['relative_growth_rate'],
			error=designs['relative_growth_rate_ci95'])
		ax.axhline(1.0, color=style.INK_MUTED, linestyle=(0, (4, 3)),
			linewidth=0.8, zorder=1)
		style.promoter_axis(ax, designs['promoter'])
		ax.set_ylim(0, max(1.2, float(np.nanmax(
			designs['relative_growth_rate'] + designs[
			'relative_growth_rate_ci95'].fillna(0))) * 1.05))
		ax.set_ylabel('growth rate relative to no GFP')
		style.save_panel(fig, plotOutDir,
			style.panel_name(__file__, 'growth'), script=__file__,
			batch=inputDir, sim_data_path=simDataFile, n=len(variants),
			notes='RBS %g designs. Per lineage: window-mean growth rate over'
				' the seed-matched control lineage; bars are the mean over'
				' surviving lineages, error bars a t-based 95%% CI. Omitted'
				' (no survivors): %s.' % (RBS, lost), **provenance)

	@staticmethod
	def _bars(ax, designs, heights, error=None):
		"""Bars at each promoter, hollow where some lineages failed."""
		color = style.rbs_color(RBS)
		for (_, row), height in zip(designs.iterrows(), heights):
			if not np.isfinite(height):
				continue
			failed = bool(row['failed'])
			ax.bar(row['promoter'], height, width=BAR_WIDTH,
				color='white' if failed else color, edgecolor=style.COLORS[
				'gfp'], linewidth=0.8, zorder=2)
		if error is not None:
			ok = np.isfinite(heights) & np.isfinite(error)
			ax.errorbar(designs['promoter'][ok], heights[ok],
				yerr=error[ok], fmt='none', ecolor=style.INK, elinewidth=0.7,
				capsize=1.5, zorder=3)


if __name__ == '__main__':
	Plot().cli()
