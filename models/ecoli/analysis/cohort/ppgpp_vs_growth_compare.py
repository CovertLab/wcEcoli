"""
Overlay ppGpp vs relative growth rate for two (or more) variant sweeps, to test
whether a growth slowdown is metabolically driven (stringent response -> high
ppGpp) or not, controlling for the AMOUNT of slowdown.

For each run, for each variant, average ppGpp and instantaneous growth rate over
the settled window (last WINDOW_GENS generations, all seeds), then normalize
growth to that run's wildtype (the lowest variant index, 0).  Plot ppGpp
(y) vs growth / wildtype-growth (x), one series per run.

Reading it: at a MATCHED amount of slowdown (same x), the run whose points sit
HIGHER in ppGpp has the more stringent / supply-limited slowdown.  Use to compare
the wildtype kcat-tightening sweep (where metabolism binds) against the GFP
trl-eff sweep (ribosome-limited): if the wildtype curve rises above the GFP curve
as growth falls, that isolates "metabolic limitation" from "generic slowdown".

This is a standalone two-run tool (the single-run analysis framework can't
overlay runs) and is variant-agnostic -- it needs only mu and ppGpp per variant,
not the estimator / trl-eff layout.  Reads recorded listeners; no re-simulation.

Usage:
	python models/ecoli/analysis/cohort/ppgpp_vs_growth_compare.py \\
		--run "wildtype kcat sweep=<run_dir_1>" \\
		--run "GFP trl-eff sweep=<run_dir_2>" \\
		[--output_dir <dir>] [--window_gens 8]
"""

import argparse
import csv
import os
import sys

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from wholecell.analysis.analysis_tools import read_stacked_columns


SEC_PER_HOUR = 3600.
RUN_COLORS = ['#d62728', '#1f77b4', '#2ca02c', '#9467bd', '#ff7f0e']


def _safe_mean(cells, table, column, scale=1.0):
	if len(cells) == 0:
		return np.nan
	try:
		data = read_stacked_columns(
			cells, table, column, remove_first=True, ignore_exception=True)
		return float(np.nanmean(data)) * scale
	except Exception:
		return np.nan


def _collect(run_dir, window_gens):
	"""Return list of (mu_rel, ppgpp, mu, variant) for a run, sorted by mu_rel."""
	ap = AnalysisPaths(run_dir, variant_plot=True)
	n_gen = ap.n_generation
	window = np.arange(max(n_gen - window_gens, 0), n_gen)
	variants = sorted(ap.get_variants())

	mu, ppgpp = {}, {}
	for vi in variants:
		cells = ap.get_cells(
			variant=[vi], generation=window, only_successful=True)
		mu[vi] = _safe_mean(
			cells, 'Mass', 'instantaneous_growth_rate', scale=SEC_PER_HOUR)
		ppgpp[vi] = _safe_mean(cells, 'GrowthLimits', 'ppgpp_conc')

	# Wildtype = lowest variant index (0 in both sweeps).
	mu_wt = mu.get(variants[0], np.nan)
	if not np.isfinite(mu_wt) or mu_wt == 0:
		# Fall back to the fastest-growing variant as the reference.
		finite = [m for m in mu.values() if np.isfinite(m)]
		mu_wt = max(finite) if finite else np.nan

	rows = []
	for vi in variants:
		if np.isfinite(mu[vi]) and np.isfinite(ppgpp[vi]) and np.isfinite(mu_wt):
			rows.append((mu[vi] / mu_wt, ppgpp[vi], mu[vi], vi))
	rows.sort(key=lambda r: r[0])
	return rows


def main(argv=None):
	parser = argparse.ArgumentParser(description=__doc__,
		formatter_class=argparse.RawDescriptionHelpFormatter)
	parser.add_argument('--run', action='append', required=True,
		metavar='LABEL=PATH',
		help='A run as "label=path"; repeat for each run to overlay.')
	parser.add_argument('--output_dir', default='.',
		help='Directory for the output PDF/CSV (default: cwd).')
	parser.add_argument('--window_gens', type=int, default=8,
		help='Number of trailing generations to average (default 8).')
	args = parser.parse_args(argv)

	runs = []
	for spec in args.run:
		if '=' not in spec:
			print(f'Error: --run must be "label=path", got: {spec}',
				file=sys.stderr)
			return 1
		label, path = spec.split('=', 1)
		path = os.path.abspath(os.path.expandvars(os.path.expanduser(path)))
		if not os.path.isdir(path):
			print(f'Error: run dir not found: {path}', file=sys.stderr)
			return 1
		runs.append((label.strip(), path))

	fig, ax = plt.subplots(figsize=(8, 6))
	csv_rows = []
	for i, (label, path) in enumerate(runs):
		rows = _collect(path, args.window_gens)
		if not rows:
			print(f'No usable data for {label} ({path}); skipping.')
			continue
		color = RUN_COLORS[i % len(RUN_COLORS)]
		x = [r[0] for r in rows]
		y = [r[1] for r in rows]
		ax.plot(x, y, '-o', color=color, lw=1.5, ms=5, label=label)
		for mu_rel, pg, mu_abs, vi in rows:
			csv_rows.append({
				'run': label, 'variant': vi,
				'growth_rate_per_h': f'{mu_abs:.6g}',
				'growth_rel_to_wt': f'{mu_rel:.6g}',
				'ppgpp_uM': f'{pg:.6g}'})
		print(f'{label}: {len(rows)} variants, '
			f'growth {min(x):.2f}-{max(x):.2f}x WT, '
			f'ppGpp {min(y):.1f}-{max(y):.1f} uM')

	ax.set_xlabel('Growth rate / wildtype growth rate')
	ax.set_ylabel('ppGpp (uM)')
	ax.set_title('ppGpp vs growth slowdown\n'
		'(higher ppGpp at matched slowdown = more stringent / '
		'metabolism-limited)')
	ax.invert_xaxis()  # more slowdown (smaller ratio) to the right
	ax.legend(fontsize=9)
	fig.tight_layout()

	out_pdf = os.path.join(args.output_dir, 'ppgpp_vs_growth_compare.pdf')
	fig.savefig(out_pdf)
	plt.close(fig)
	print(f'Wrote {out_pdf}')

	if csv_rows:
		out_csv = os.path.join(args.output_dir, 'ppgpp_vs_growth_compare.csv')
		with open(out_csv, 'w', newline='', encoding='utf-8') as fh:
			writer = csv.DictWriter(fh, fieldnames=[
				'run', 'variant', 'growth_rate_per_h', 'growth_rel_to_wt',
				'ppgpp_uM'])
			writer.writeheader()
			writer.writerows(csv_rows)
		print(f'Wrote {out_csv}')
	return 0


if __name__ == '__main__':
	sys.exit(main())
