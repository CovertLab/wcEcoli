"""
Amino-acid supply / demand / concentration / per-AA charging diagnosis for the
new_gene_trl_eff_v2_kcat_multiplier_sweep (GFP burden crossed with a kcat
multiplier, two estimators).

The 3-D (estimator x trl-eff x multiplier) counterpart of the AA analyses on the
2-D GFP and wildtype sweeps.  The AA render core is 2-D (group x x-axis), so this
slices the cube the way the sweep's sibling analyses do (_comparison /
_growth_limitation / _mass): x = kcat multiplier (the tightening axis expected to
make metabolism bite), one line / heatmap panel per GFP translation-efficiency
level, and a separate figure set per estimator.

For each estimator (gap, smoothed_max) it emits the full 5-figure AA set, suffixed
by estimator short name:
  <name>_<est>.pdf                     overview (supply, demand, supply/demand,
                                       charging, AA conc, ppGpp)
  <name>_<est>_supply_composition.pdf  synthesis / import / export
  <name>_<est>_charging_by_aa.pdf      per-AA charged-tRNA fraction heatmap
  <name>_<est>_conc_by_aa.pdf          per-AA free AA concentration, fold vs WT
  <name>_<est>_supplydemand_by_aa.pdf  per-AA supply/demand, fold vs WT

Question: as the bound tightens (multiplier down) under GFP load, does the
supply-limited signature finally appear (supply below demand, PHE charging
collapsing, ppGpp rising), and at what multiplier?  This maps the boundary between
the mechanical (ribosome-dilution) and controller-mediated (stringent) regimes,
AA by AA.  The multiplier = 1.0 column reproduces the GFP-only response; the
trl_eff = 0.0 line reproduces that estimator's pure multiplier sweep.

Reuses the sweep-agnostic collect/render core; reads only recorded listeners --
no re-simulation.
"""

import numpy as np
from matplotlib import pyplot as plt

from models.ecoli.analysis import variantAnalysisPlot
from models.ecoli.analysis.variant.new_gene_trl_eff_v2_estimator_sweep_aa_supply_demand import (
	AALayout,
	collect_aa_metrics,
	render_all,
)
from models.ecoli.sim.variants.new_gene_trl_eff_v2_kcat_multiplier_sweep import (
	ESTIMATORS,
	ESTIMATOR_SHORT,
	MULTIPLIERS,
	TRL_EFF_VALUES,
	is_wildtype,
	variant_to_cell,
)


# Match the sibling multiplier-sweep analyses: colour lines by trl-eff level.
TRL_COLORS = plt.cm.viridis(np.linspace(0, 0.9, len(TRL_EFF_VALUES)))


class Plot(variantAnalysisPlot.VariantAnalysisPlot):
	def do_plot(self, inputDir, plotOutDir, plotOutFileName, simDataFile,
				validationDataFile, metadata):
		variant_indexes = self.ap.get_variants()
		n_total_gens = self.ap.n_generation
		ignore_first_n_gens = max(n_total_gens - 8, 0)
		window = np.arange(ignore_first_n_gens, n_total_gens)

		# Collect once over all variants (both estimator blocks + the shared
		# wildtype reuse the same per-variant metrics).
		metrics = collect_aa_metrics(self.ap, variant_indexes, window)
		wt = next((vi for vi in variant_indexes if is_wildtype(vi)), None)

		groups = [(f'trl {te:g}', TRL_COLORS[ti])
			for ti, te in enumerate(TRL_EFF_VALUES)]

		# Facet by estimator: one full figure set per estimator block.
		for est_idx, estimator in enumerate(ESTIMATORS):
			assign = {}
			for vi in variant_indexes:
				if is_wildtype(vi):
					continue
				e_idx, trl_idx, mult_idx = variant_to_cell(vi)
				if e_idx != est_idx:
					continue
				assign[vi] = (trl_idx, mult_idx)
			short = ESTIMATOR_SHORT[estimator]
			layout = AALayout(
				groups=groups, assign=assign,
				x_values=np.array(MULTIPLIERS), xlabel='kcat multiplier',
				invert_x=True, wt=wt,
				title_prefix=(
					'new_gene_trl_eff_v2_kcat_multiplier_sweep '
					f'[{short}] (line/panel per trl-eff)'))
			render_all(plotOutDir, f'{plotOutFileName}_{short}', metadata,
				metrics, layout)


if __name__ == '__main__':
	Plot().cli()
