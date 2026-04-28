"""
Audit per-category mass doubling across a cohort.

For every numeric column written by the Mass listener, compute the
end/start ratio per cell and aggregate across seeds and generations.
A healthy WT cell should roughly double every category each generation;
categories whose pooled median ratio falls outside DOUBLING_TOL of 2.0
are flagged as "not doubling".

Outputs (in plotOutDir):
	- mass_category_doubling.csv: per-(generation, category) statistics
	  plus a pooled-across-generations summary.
	- mass_category_doubling.{pdf,png}: boxplot of pooled ratios per
	  category (top) and median-ratio-vs-generation line plot (bottom).
"""

import csv
import os

import numpy as np
from matplotlib import pyplot as plt

from wholecell.io.tablereader import TableReader
from wholecell.analysis.analysis_tools import exportFigure
from models.ecoli.analysis import cohortAnalysisPlot


DOUBLING_TARGET = 2.0
DOUBLING_TOL = 0.10  # |median - 2| / 2 must be <= this to pass
SKIP_INITIAL_GENS = 4  # Burn-in: ignore the first N generations.

# Mass-listener columns that are not per-cell mass quantities.
NON_MASS_COLUMNS = {
	"time",
	"simulationStep",
	"growth",
	"instantaneous_growth_rate",
	"cellVolume",
	"processMassDifferences",
	}


class Plot(cohortAnalysisPlot.CohortAnalysisPlot):
	def do_plot(self, variantDir, plotOutDir, plotOutFileName, simDataFile,
			validationDataFile, metadata):
		all_cells = self.ap.get_cells()
		if len(all_cells) == 0:
			print("No cells found -- skipping mass_category_doubling.")
			return

		categories = self._discover_categories(all_cells[0])
		if not categories:
			print("No mass categories discovered -- skipping.")
			return

		n_gen = self.ap.n_generation
		gen_indices = list(range(SKIP_INITIAL_GENS, n_gen))
		if not gen_indices:
			print(f"Cohort has {n_gen} generations; need more than "
				f"SKIP_INITIAL_GENS={SKIP_INITIAL_GENS} to analyze.")
			return
		# ratios[gen][category] -> list of end/start ratios across seeds
		ratios = {g: {c: [] for c in categories} for g in gen_indices}

		for gen_idx in gen_indices:
			gen_cells = self.ap.get_cells(generation=[gen_idx])
			for sim_dir in gen_cells:
				sim_out = os.path.join(sim_dir, "simOut")
				try:
					mass = TableReader(os.path.join(sim_out, "Mass"))
				except Exception as exc:
					print(f"  skipping {sim_dir}: {exc}")
					continue

				for cat in categories:
					try:
						col = mass.readColumn(cat)
					except Exception:
						continue
					if col.ndim != 1 or col.size < 2:
						continue
					m0 = col[0]
					mF = col[-1]
					if not np.isfinite(m0) or not np.isfinite(mF) or m0 <= 0:
						continue
					ratios[gen_idx][cat].append(mF / m0)

		self._write_csv(plotOutDir, categories, ratios, gen_indices)
		self._make_plot(plotOutDir, plotOutFileName, metadata, categories,
			ratios, gen_indices)
		self._print_fail_list(categories, ratios, gen_indices)

	def _discover_categories(self, sim_dir):
		sim_out = os.path.join(sim_dir, "simOut")
		mass = TableReader(os.path.join(sim_out, "Mass"))
		cols = sorted(mass.columnNames())
		return [c for c in cols if c not in NON_MASS_COLUMNS]

	def _pooled(self, ratios, cat, gen_indices):
		out = []
		for g in gen_indices:
			out.extend(ratios[g][cat])
		return np.array(out, dtype=float)

	def _stats(self, values):
		if values.size == 0:
			return dict(n=0, median=np.nan, mean=np.nan, std=np.nan,
				p05=np.nan, p95=np.nan)
		return dict(
			n=int(values.size),
			median=float(np.median(values)),
			mean=float(np.mean(values)),
			std=float(np.std(values)),
			p05=float(np.percentile(values, 5)),
			p95=float(np.percentile(values, 95)),
			)

	def _doubling_ok(self, median):
		if not np.isfinite(median):
			return False
		return abs(median - DOUBLING_TARGET) / DOUBLING_TARGET <= DOUBLING_TOL

	def _write_csv(self, plotOutDir, categories, ratios, gen_indices):
		csv_path = os.path.join(plotOutDir, "mass_category_doubling.csv")
		with open(csv_path, "w", newline="") as fh:
			writer = csv.writer(fh)
			writer.writerow(["scope", "generation", "category", "n_cells",
				"median_ratio", "mean_ratio", "std_ratio", "p05", "p95",
				"doubling_ok"])

			for g in gen_indices:
				for cat in categories:
					values = np.array(ratios[g][cat], dtype=float)
					s = self._stats(values)
					writer.writerow(["per_gen", g, cat, s["n"],
						s["median"], s["mean"], s["std"], s["p05"], s["p95"],
						self._doubling_ok(s["median"])])

			writer.writerow([])
			scope = f"pooled_gens_{gen_indices[0]}_to_{gen_indices[-1]}"
			for cat in categories:
				values = self._pooled(ratios, cat, gen_indices)
				s = self._stats(values)
				writer.writerow([scope, "all", cat, s["n"],
					s["median"], s["mean"], s["std"], s["p05"], s["p95"],
					self._doubling_ok(s["median"])])

	def _make_plot(self, plotOutDir, plotOutFileName, metadata, categories,
			ratios, gen_indices):
		pooled = {cat: self._pooled(ratios, cat, gen_indices)
			for cat in categories}
		drawable = [c for c in categories if pooled[c].size > 0]
		if not drawable:
			print("No drawable categories -- skipping plot.")
			return

		# Sort by median so outliers stand out.
		drawable.sort(key=lambda c: np.median(pooled[c]))
		data = [pooled[c] for c in drawable]

		fig, (ax_box, ax_gen) = plt.subplots(2, 1, figsize=(
			max(8, 0.45 * len(drawable)), 10))

		ax_box.boxplot(data, showfliers=False)
		ax_box.set_xticks(range(1, len(drawable) + 1))
		ax_box.set_xticklabels(drawable, rotation=60, ha="right")
		ax_box.axhline(DOUBLING_TARGET, color="k", linestyle="--",
			linewidth=1, label="target = 2")
		ax_box.axhline(DOUBLING_TARGET * (1 + DOUBLING_TOL), color="r",
			linestyle=":", linewidth=0.8)
		ax_box.axhline(DOUBLING_TARGET * (1 - DOUBLING_TOL), color="r",
			linestyle=":", linewidth=0.8,
			label=f"+/- {int(DOUBLING_TOL * 100)}%")
		ax_box.set_ylabel("end / start mass ratio (per cell)")
		ax_box.set_title(
			"Per-category mass doubling, pooled across seeds "
			f"& generations {gen_indices[0]}-{gen_indices[-1]}")
		ax_box.legend(loc="upper left", fontsize=8)

		gens = np.array(gen_indices)
		for cat in drawable:
			medians = []
			for g in gens:
				vals = np.array(ratios[g][cat], dtype=float)
				medians.append(np.median(vals) if vals.size else np.nan)
			ax_gen.plot(gens, medians, marker="o", label=cat, linewidth=1)
		ax_gen.axhline(DOUBLING_TARGET, color="k", linestyle="--",
			linewidth=1)
		ax_gen.set_xlabel("generation")
		ax_gen.set_ylabel("median end/start ratio")
		ax_gen.set_title("Median ratio vs. generation")
		ax_gen.legend(loc="center left", bbox_to_anchor=(1.0, 0.5),
			fontsize=7, ncol=1)

		plt.tight_layout()
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close("all")

	def _print_fail_list(self, categories, ratios, gen_indices):
		fails = []
		for cat in categories:
			values = self._pooled(ratios, cat, gen_indices)
			if values.size == 0:
				continue
			median = float(np.median(values))
			if not self._doubling_ok(median):
				fails.append((cat, median, values.size))

		if not fails:
			print("mass_category_doubling: all categories within "
				f"+/- {int(DOUBLING_TOL * 100)}% of 2x.")
			return

		print("mass_category_doubling: categories outside tolerance "
			f"(target={DOUBLING_TARGET}, tol={DOUBLING_TOL}):")
		for cat, median, n in sorted(fails, key=lambda r: abs(r[1] - 2.0),
				reverse=True):
			print(f"  {cat}: median={median:.3f} (n={n})")


if __name__ == "__main__":
	Plot().cli()
