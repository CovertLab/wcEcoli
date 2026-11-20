"""
GFP made per cell per hour against GFP's share of the proteome.

# Figure 6D

One point per design, coloured by RBS strength, hollow where some lineages
failed. GFP made per cell per hour is the protein copies a cell must make to
hold its GFP level while growing: GFP monomers x (growth rate + GFP
degradation rate), window means over surviving lineages. The gfp translation
initiation rate is recorded in the sidecar as an independent check.

The dashed curve is "share x growth": production proportional to GFP share
times relative growth, with relative growth the straight-line fit of 2C
(control included) and the one proportionality constant fitted by least
squares through the origin. Its peak is where more GFP per cell is outweighed
by slower growth; the peak share is in the sidecar.

There is no growth axis here, so one version only.

The supporting claim: GFP output per cell per hour peaks at moderate
expression, because past that the burden starves GFP's own production.
"""

import pickle

import numpy as np
from matplotlib import pyplot as plt
from matplotlib.lines import Line2D

from models.ecoli.analysis import variantAnalysisPlot
from wholecell.analysis import gfp_paper_constants as C
from wholecell.analysis import gfp_paper_data
from wholecell.analysis import gfp_paper_style as style

PERCENT = 100.0
SECONDS_PER_HOUR = 3600.0


class Plot(variantAnalysisPlot.VariantAnalysisPlot):
	def do_plot(self, inputDir, plotOutDir, plotOutFileName, simDataFile,
			validationDataFile, metadata):
		from wholecell.utils import units
		data = gfp_paper_data.load(inputDir, simDataFile, metadata,
			plotOutDir, self.cpus)
		with open(simDataFile, 'rb') as f:
			sim_data = pickle.load(f)
		monomers = sim_data.process.translation.monomer_data
		k = list(monomers['id']).index(data.gene.monomer_id)
		kdeg = float(monomers['deg_rate'][k].asNumber(1 / units.s)) * (
			SECONDS_PER_HOUR)
		designs = data.designs[(data.designs['promoter'] > 0)
			& (data.designs['n_lineages_survived'] > 0)
			& np.isfinite(data.designs['gfp_proteome_share'])]
		share = designs['gfp_proteome_share'].values
		made = designs['gfp_monomers'].values * (
			designs['growth_rate_per_h'].values + kdeg)
		initiated = designs['gfp_translation_per_s'].values * SECONDS_PER_HOUR

		# share x growth: relative growth = a + b s (as 2C, control at 1)
		complete = ~designs['failed'].values
		fit_s = np.append(0.0, share[complete])
		fit_g = np.append(1.0, designs['relative_growth_rate'].values[
			complete])
		b, a = np.polyfit(fit_s, fit_g, 1) if fit_s.size >= 2 else (np.nan,
			np.nan)
		basis = share * (a + b * share)
		scale = float(np.dot(basis, made) / np.dot(basis, basis)) if np.dot(
			basis, basis) > 0 else np.nan
		grid = np.linspace(0, 0.9, 200)
		curve = scale * grid * (a + b * grid)
		peak = float(-a / (2 * b)) if b < 0 else np.nan

		style.apply()
		fig, ax = plt.subplots(figsize=style.panel_size('single', 0.8))
		keep = curve >= 0
		ax.plot(PERCENT * grid[keep], curve[keep], color=style.INK_MUTED,
			linestyle=(0, (4, 3)), linewidth=1.0, zorder=1)
		handles = []
		for rbs in sorted(designs['rbs'].unique()):
			rows = designs['rbs'].values == rbs
			color = style.rbs_color(rbs)
			style.design_scatter(ax, PERCENT * share[rows], made[rows], color,
				failed=designs['failed'].values[rows], zorder=3)
			handles.append(Line2D([], [], linestyle='none', marker='o',
				color=color, label=C.rbs_label(rbs)))
		if designs['failed'].any():
			handles.append(style.failed_legend_handle())
		handles.append(Line2D([], [], color=style.INK_MUTED,
			linestyle=(0, (4, 3)), label='share × growth fit'))
		style.legend_below(fig, handles, ncol=2)
		ax.set_xlim(0, 90)
		ax.set_ylim(bottom=0)
		ax.set_xlabel('GFP share of proteome (%)')
		ax.set_ylabel('GFP made per cell per hour')
		ax.yaxis.set_major_formatter(plt.FuncFormatter(
			lambda v, _: '%gk' % (v / 1e3) if v else '0'))

		style.save_panel(fig, plotOutDir, style.panel_name(__file__),
			script=__file__, batch=inputDir, sim_data_path=simDataFile,
			n=int(len(designs)),
			notes='GFP made = monomers x (growth rate + degradation %.3g /h).'
				' Fit: relative growth = %.4f + %.4f x share (fraction);'
				' scale %.4g; peak at share %.3f. Per design (variant: share,'
				' made per h, initiations per h): %s.' % (kdeg, a, b, scale,
				peak, ', '.join('%d: %.3f, %.0f, %.0f' % (v, s, m, i) for v,
				s, m, i in zip(designs['variant'], share, made, initiated))),
			**data.provenance(sorted(int(v) for v in designs['variant'])))


if __name__ == '__main__':
	Plot().cli()
