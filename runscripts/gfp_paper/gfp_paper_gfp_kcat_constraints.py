#!/usr/bin/env python
"""
Growth with GFP, with and without a kcat bound, and the difference.

# Figure 4F

Arms are found by reading each variant's kcat condition from its own sim_data
(gfp_paper_data.kcat_conditions), so the bounded and unbounded GFP designs
may sit in one batch or in two (for example M2 for the bound, B1 without).
Designs are matched across arms by (promoter, RBS) and lineages by seed.
Both arms are normalised to the unbounded no-GFP control of the same seed,
taken from the same batch when it has one and otherwise from another batch
given.

Top: each arm's design means over surviving lineages against GFP share, with
t-based 95% CIs. Bottom: bound minus no bound, per seed pair, mean with a
t-based 95% CI, at the mean of the two arms' GFP shares. Two versions while
decision D4 is open:
	growth rate     relative growth rate and its difference (plain name)
	doubling time   minutes and its difference (suffix DOUBLING_TIME_PART)

The supporting claim: enforcing the enzyme-supply ceiling barely changes
growth with GFP, so enzyme availability is not what slows these cells.

Usage:
	python runscripts/gfp_paper/gfp_paper_gfp_kcat_constraints.py \\
		<batch> [<batch> ...] [--out DIR] [--cpus N]
"""

import argparse
import os

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

from wholecell.analysis import gfp_paper_constants as C
from wholecell.analysis import gfp_paper_data
from wholecell.analysis import gfp_paper_style as style

PERCENT = 100.0


def lineage_rows(batches, cpus):
	"""Every lineage in the batches, with its batch, arm and design."""
	frames = []
	for batch in batches:
		data = gfp_paper_data.load_batch(batch, cpus)
		variants = sorted(int(v) for v in data.lineages['variant'].unique())
		conditions = gfp_paper_data.kcat_conditions(batch, variants,
			os.path.join(batch, 'plotOut'))
		frame = data.lineages.merge(conditions[['variant', 'bound']],
			on='variant')
		frame['batch'] = batch
		frames.append(frame)
	return pd.concat(frames, ignore_index=True)


def normalise(lineages, column):
	"""Each lineage's value over its seed's unbounded no-GFP control."""
	controls = lineages[~lineages['bound'] & (lineages['promoter'] == 0)]
	same = controls.set_index(['batch', 'seed'])[column]
	any_batch = controls.groupby('seed')[column].mean()
	out = []
	for _, row in lineages.iterrows():
		ref = same.get((row['batch'], row['seed']), np.nan)
		if not np.isfinite(ref):
			ref = any_batch.get(row['seed'], np.nan)
		out.append(row[column] / ref if ref and np.isfinite(ref) else np.nan)
	return np.array(out)


def summarise(lineages, column):
	"""Design means and CIs per arm, and paired differences between arms."""
	alive = lineages[lineages['survived'] & (lineages['promoter'] > 0)]
	arms = alive.groupby(['bound', 'promoter', 'rbs']).agg(
		value=(column, 'mean'), share=('gfp_proteome_share', 'mean'),
		n=(column, 'count')).reset_index()
	arms['ci'] = [gfp_paper_data._ci95(alive[(alive['bound'] == b)
		& (alive['promoter'] == p) & (alive['rbs'] == r)][column].dropna()
		.values) for b, p, r in zip(arms['bound'], arms['promoter'],
		arms['rbs'])]
	failed = lineages[lineages['promoter'] > 0].groupby(['bound', 'promoter',
		'rbs'])['survived'].apply(lambda s: not s.all())
	arms['failed'] = [bool(failed.get((b, p, r), False)) for b, p, r in
		zip(arms['bound'], arms['promoter'], arms['rbs'])]
	deltas = []
	for (p, r), group in alive.groupby(['promoter', 'rbs']):
		on = group[group['bound']].set_index('seed')[column]
		off = group[~group['bound']].set_index('seed')[column]
		paired = (on - off).dropna()
		if paired.empty:
			continue
		deltas.append(dict(promoter=p, rbs=r, delta=float(paired.mean()),
			ci=gfp_paper_data._ci95(paired.values), n=int(paired.size),
			share=float(group['gfp_proteome_share'].mean())))
	return arms, pd.DataFrame(deltas)


def main():
	parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
	parser.add_argument('batches', nargs='+')
	parser.add_argument('--out', default=None,
		help='plotOut directory; default the first batch\'s plotOut')
	parser.add_argument('--cpus', type=int, default=1)
	args = parser.parse_args()
	out = args.out or os.path.join(args.batches[0], 'plotOut')
	lineages = lineage_rows(args.batches, args.cpus)
	if not (lineages['bound'] & (lineages['promoter'] > 0)).any() or not (
			~lineages['bound'] & (lineages['promoter'] > 0)).any():
		raise SystemExit('Need GFP designs both with and without a kcat'
			' bound across the batches given.')
	lineages['relative_growth'] = normalise(lineages, 'growth_rate_per_h')

	for part, column, ylabel, dlabel in (
			(None, 'relative_growth', 'growth rate relative to wild type',
				'Δ'),
			(C.DOUBLING_TIME_PART, 'doubling_time_min', 'doubling time (min)',
				'Δ (min)')):
		arms, deltas = summarise(lineages, column)
		style.apply()
		fig, (top, bottom) = plt.subplots(2, 1, sharex=True,
			figsize=style.panel_size('single', 0.85),
			gridspec_kw={'height_ratios': [3, 1]})
		for bound, color, label in ((False, style.INK_MUTED,
				'no k$_{cat}$ bound'), (True, style.COLORS['enzyme'],
				'k$_{cat}$ bound enforced')):
			arm = arms[arms['bound'] == bound].sort_values('share')
			top.errorbar(PERCENT * arm['share'], arm['value'], yerr=np.where(
				np.isfinite(arm['ci']), arm['ci'], 0), fmt='none',
				ecolor=color, elinewidth=0.6, alpha=0.6, zorder=2)
			style.design_scatter(top, PERCENT * arm['share'], arm['value'],
				color, failed=arm['failed'].values, label=label, zorder=3)
		if part is None:
			top.axhline(1.0, color=style.INK_MUTED, linestyle=(0, (4, 3)),
				linewidth=0.8, zorder=1)
		top.set_ylabel(ylabel)
		style.legend_with_failed(top, bool(arms['failed'].any()),
			loc='upper right' if part is not None else 'lower left')
		if len(deltas):
			bottom.errorbar(PERCENT * deltas['share'], deltas['delta'],
				yerr=np.where(np.isfinite(deltas['ci']), deltas['ci'], 0),
				fmt='o', color=style.INK, markersize=3, elinewidth=0.7,
				capsize=1.5, zorder=3)
		bottom.axhline(0.0, color=style.INK_MUTED, linestyle=(0, (4, 3)),
			linewidth=0.8, zorder=1)
		bottom.set_ylabel(dlabel)
		bottom.set_xlabel('GFP share of proteome (%)')
		bottom.set_xlim(left=0)

		style.save_panel(fig, out, style.panel_name(__file__, part),
			script=__file__, batch=';'.join(args.batches),
			sim_data_path=os.path.join(args.batches[0], 'kb',
			'simData.cPickle'),
			generation_window=(C.ANALYSIS_FIRST_GEN, C.ANALYSIS_LAST_GEN),
			n=int(len(deltas)),
			n_lineages_total=int((lineages['promoter'] > 0).sum()),
			n_lineages_survived=int((lineages['survived']
				& (lineages['promoter'] > 0)).sum()),
			notes='%s. Arms from each variant sim_data. Paired differences'
				' (bound - no bound) per design: %s.' % (ylabel, ', '.join(
				'p%g/r%g %.3f +/- %.3f (n %d)' % (d.promoter, d.rbs, d.delta,
				d.ci, d.n) for d in deltas.itertuples())))


if __name__ == '__main__':
	main()
