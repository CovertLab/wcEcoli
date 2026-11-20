#!/usr/bin/env python
"""
ppGpp against growth when growth is slowed by kcat tightening versus by GFP.

# Figure 4H

Two series, each one point per variant (means over surviving lineages of
each lineage's window means):
	kcat tightening, no GFP   the no-GFP variants of the kcat batch (M1): the
	                          unbound reference and every bounded variant
	                          using KCAT_ESTIMATOR (all bounded variants if
	                          none use it); growth over the same seed's
	                          unbound reference
	GFP expression, no bound  the variants of the GFP batch (B1) that run
	                          without a kcat bound, control included; growth
	                          relative to the seed-matched control
Bound or not is read from each variant's sim_data, so a GFP batch that
carries a bound is excluded rather than silently mixed in (the storyboard's
version of this series ran with a 1.0x bound). Lines join points in order of
growth as a guide.

Two versions while decision D4 is open: relative growth rate on x (plain
name) and doubling time in minutes on x (suffix DOUBLING_TIME_PART).

The supporting claim: a real metabolic limit raises ppGpp sharply as growth
slows; GFP slows growth without that signature.

Usage:
	python runscripts/gfp_paper/gfp_paper_ppgpp_vs_growth_kcat_tightening_vs_gfp.py \\
		--kcat-batch M1_DIR --gfp-batch B1_DIR [--out DIR] [--cpus N]
"""

import argparse
import os

import numpy as np
from matplotlib import pyplot as plt

from wholecell.analysis import gfp_paper_constants as C
from wholecell.analysis import gfp_paper_data
from wholecell.analysis import gfp_paper_style as style


def series(batch, cpus, kcat):
	"""(variant, relative growth, doubling time, ppGpp) per kept variant."""
	data = gfp_paper_data.load_batch(batch, cpus)
	variants = sorted(int(v) for v in data.designs['variant'])
	conditions = gfp_paper_data.kcat_conditions(batch, variants,
		os.path.join(batch, 'plotOut')).merge(
		data.designs[['variant', 'promoter']], on='variant')
	if kcat:
		pool = conditions[conditions['promoter'] == 0]
		unbound = pool[~pool['bound']]
		if unbound.empty:
			raise SystemExit('%s has no unbound no-GFP variant' % batch)
		reference = int(unbound['variant'].min())
		bounded = pool[pool['bound']]
		chosen = bounded[bounded['estimator'] == C.KCAT_ESTIMATOR]
		keep = [reference] + sorted(int(v) for v in (chosen if len(chosen)
			else bounded)['variant'])
	else:
		keep = sorted(int(v) for v in conditions[~conditions['bound']][
			'variant'])
		reference = gfp_paper_data.CONTROL_VARIANT
		if reference not in keep:
			raise SystemExit('%s: the control variant runs with a kcat'
				' bound' % batch)
	alive = data.lineages[data.lineages['survived']]
	rows = []
	for variant in keep:
		ratio = gfp_paper_data.seed_matched_ratio(data.lineages, variant,
			'growth_rate_per_h', baseline=reference)
		group = alive[alive['variant'] == variant]
		rows.append((variant, 1.0 if variant == reference else ratio['mean'],
			float(group['doubling_time_min'].mean()),
			float(group['ppgpp_conc'].mean())))
	return np.array(rows, dtype=float), data


def main():
	parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
	parser.add_argument('--kcat-batch', required=True)
	parser.add_argument('--gfp-batch', required=True)
	parser.add_argument('--out', default=None,
		help='plotOut directory; default the kcat batch\'s plotOut')
	parser.add_argument('--cpus', type=int, default=1)
	args = parser.parse_args()
	out = args.out or os.path.join(args.kcat_batch, 'plotOut')
	kcat_rows, _ = series(args.kcat_batch, args.cpus, kcat=True)
	gfp_rows, _ = series(args.gfp_batch, args.cpus, kcat=False)

	for part, col, xlabel in ((None, 1, 'growth rate relative to wild type'),
			(C.DOUBLING_TIME_PART, 2, 'doubling time (min)')):
		style.apply()
		fig, ax = plt.subplots(figsize=style.panel_size('single', 0.7))
		for rows, color, label in ((kcat_rows, style.COLORS['enzyme'],
				'k$_{cat}$ tightening, no GFP'), (gfp_rows,
				style.COLORS['gfp'], 'GFP expression, no tightening')):
			rows = rows[np.isfinite(rows[:, col]) & np.isfinite(rows[:, 3])]
			rows = rows[np.argsort(rows[:, col])]
			ax.plot(rows[:, col], rows[:, 3], color=color, linewidth=1.0,
				zorder=2)
			ax.scatter(rows[:, col], rows[:, 3], s=16, color=color,
				edgecolors='white', linewidths=0.5, label=label, zorder=3)
		ax.set_xlabel(xlabel)
		ax.set_ylabel('ppGpp (µM)')
		ax.set_ylim(bottom=0)
		ax.legend(loc='upper right' if part is None else 'upper left')
		style.save_panel(fig, out, style.panel_name(__file__, part),
			script=__file__, batch='%s;%s' % (args.kcat_batch,
			args.gfp_batch), sim_data_path=os.path.join(args.kcat_batch, 'kb',
			'simData.cPickle'),
			generation_window=(C.ANALYSIS_FIRST_GEN, C.ANALYSIS_LAST_GEN),
			n=int(len(kcat_rows) + len(gfp_rows)),
			notes='Per variant (variant, relative growth, doubling time,'
				' ppGpp uM). kcat batch: %s. GFP batch: %s.' % (
				'; '.join('%d, %.3f, %.1f, %.1f' % tuple(r) for r in
				kcat_rows), '; '.join('%d, %.3f, %.1f, %.1f' % tuple(r)
				for r in gfp_rows)))


if __name__ == '__main__':
	main()
