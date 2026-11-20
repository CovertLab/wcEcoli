"""
Growth of cells without GFP when a kcat bound forces a smaller perceived
enzyme concentration.

# Figure 4E

Each variant's kcat condition is read from its own modified sim_data
(gfp_paper_data.kcat_conditions), not from its index, so any variant layout
works. A bound at multiplier m scales every kcat estimate by m, which acts
like a (1 - m) reduction in enzyme concentration; the x axis is
100 x (1 - m). One line per estimator. The reference is the unbound no-GFP
variant with the lowest index.

Points are means over surviving lineages with t-based 95% CIs. Two versions
while decision D4 is open:
	growth rate     each lineage over the same seed's reference lineage;
	                dashed line at 1 (plain name)
	doubling time   minutes, dashed line at the reference mean (suffix
	                DOUBLING_TIME_PART)

The supporting claim: forcing enzymes to look scarcer does slow growth, so
the model can produce a metabolic limit when one is imposed.
"""

import numpy as np
from matplotlib import pyplot as plt

from models.ecoli.analysis import variantAnalysisPlot
from wholecell.analysis import gfp_paper_constants as C
from wholecell.analysis import gfp_paper_data
from wholecell.analysis import gfp_paper_style as style

PERCENT = 100.0


class Plot(variantAnalysisPlot.VariantAnalysisPlot):
	def do_plot(self, inputDir, plotOutDir, plotOutFileName, simDataFile,
			validationDataFile, metadata):
		data = gfp_paper_data.load(inputDir, simDataFile, metadata,
			plotOutDir, self.cpus)
		variants = sorted(int(v) for v in data.designs['variant'])
		conditions = gfp_paper_data.kcat_conditions(inputDir, variants,
			plotOutDir).merge(data.designs[['variant', 'promoter']],
			on='variant')
		unbound = conditions[~conditions['bound'] & (conditions['promoter']
			== 0)]
		bounded = conditions[conditions['bound'] & (conditions['promoter']
			== 0)]
		if unbound.empty or bounded.empty:
			print('Panel 4E needs an unbound and at least one kcat-bounded'
				' no-GFP variant; this batch has %d and %d. Skipping.'
				% (len(unbound), len(bounded)))
			return
		reference = int(unbound['variant'].min())
		alive = data.lineages[data.lineages['survived']]
		reference_dt = float(alive[alive['variant'] == reference][
			'doubling_time_min'].mean())

		for part, column, ylabel, ref_value in (
				(None, 'growth_rate_per_h', 'growth rate relative to wild'
					' type', 1.0),
				(C.DOUBLING_TIME_PART, 'doubling_time_min',
					'doubling time (min)', reference_dt)):
			style.apply()
			fig, ax = plt.subplots(figsize=style.panel_size('single', 0.62))
			notes = []
			estimators = sorted(bounded['estimator'].fillna('unknown').unique())
			for k, estimator in enumerate(estimators):
				rows = bounded[bounded['estimator'].fillna('unknown')
					== estimator].sort_values('multiplier', ascending=False)
				xs, ys, cis = [0.0], [ref_value], [np.nan]
				for _, row in rows.iterrows():
					variant = int(row['variant'])
					if part is None:
						r = gfp_paper_data.seed_matched_ratio(data.lineages,
							variant, column, baseline=reference)
						mean, ci = r['mean'], r['ci95']
					else:
						values = alive[alive['variant'] == variant][column]
						mean = float(values.mean()) if len(values) else np.nan
						ci = gfp_paper_data._ci95(values.dropna().values)
					xs.append(PERCENT * (1 - row['multiplier']))
					ys.append(mean)
					cis.append(ci)
					notes.append('%s x%.2f: %.3f +/- %.3f' % (estimator,
						row['multiplier'], mean, ci))
				color = style.INK if len(estimators) == 1 else (
					style.GROWTH_RAMP[min(4, 2 + k)])
				ax.errorbar(xs, ys, yerr=np.where(np.isfinite(cis), cis, 0),
					color=color, marker='o', markersize=4, linewidth=1.0,
					elinewidth=0.7, capsize=1.5, label=None if len(estimators)
					== 1 else estimator, zorder=3)
			ax.axhline(ref_value, color=style.INK_MUTED,
				linestyle=(0, (4, 3)), linewidth=0.8, zorder=1)
			ax.set_xlabel('forced reduction in perceived enzyme'
				' concentration (%)')
			ax.set_ylabel(ylabel)
			if part is None:
				ax.set_ylim(bottom=0)
			if len(estimators) > 1:
				ax.legend(loc='lower left')
			style.save_panel(fig, plotOutDir,
				style.panel_name(__file__, part), script=__file__,
				batch=inputDir, sim_data_path=simDataFile,
				n=int(len(bounded)),
				notes='Reference = unbound variant %d. Conditions read from'
					' each variant sim_data. %s.' % (reference,
					'; '.join(notes)),
				**data.provenance(variants))


if __name__ == '__main__':
	Plot().cli()
