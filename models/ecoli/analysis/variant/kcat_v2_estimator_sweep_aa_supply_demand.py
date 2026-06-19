"""
Amino-acid supply / demand / concentration / per-AA charging diagnosis for the
kcat_v2_estimator_sweep (wildtype, gap / smoothed_max / drop_top_20 estimators x a
kcat multiplier grid, no GFP).

The wildtype-multiplier counterpart of
new_gene_trl_eff_v2_estimator_sweep_aa_supply_demand.py.  There, GFP ribosome
competition drove the burden and metabolism stayed slack; here the kcat bound is
tightened directly, so this is the regime in which the AA supply *should*
eventually fall short of demand.  This analysis is the positive control for the
supply-limited signature -- as the multiplier tightens we expect AA supply to drop
below demand, free AA pools and charged-tRNA fraction to fall, and ppGpp to rise.

Reuses the sweep-agnostic collect/render core from the GFP AA analysis; only the
layout differs (x = kcat multiplier, one series per estimator).  Reads only
recorded listeners -- no re-simulation.
"""

import numpy as np

from models.ecoli.analysis import variantAnalysisPlot
from models.ecoli.analysis.variant.kcat_v2_estimator_sweep_growth_limitation import (
	ESTIMATOR_BASE_COLORS)
from models.ecoli.analysis.variant.new_gene_trl_eff_v2_estimator_sweep_aa_supply_demand import (
	AALayout,
	collect_aa_metrics,
	render_all,
)
from models.ecoli.sim.variants.kcat_v2_estimator_sweep import (
	ESTIMATOR_LABELS,
	ESTIMATOR_SHORT_NAMES,
	MULTIPLIERS,
	is_wildtype,
	variant_to_estimator,
)


class Plot(variantAnalysisPlot.VariantAnalysisPlot):
	def do_plot(self, inputDir, plotOutDir, plotOutFileName, simDataFile,
				validationDataFile, metadata):
		variant_indexes = self.ap.get_variants()
		n_total_gens = self.ap.n_generation
		ignore_first_n_gens = max(n_total_gens - 8, 0)
		window = np.arange(ignore_first_n_gens, n_total_gens)

		metrics = collect_aa_metrics(self.ap, variant_indexes, window)

		wt = next((vi for vi in variant_indexes if is_wildtype(vi)), None)
		assign = {}
		for vi in variant_indexes:
			if is_wildtype(vi):
				continue
			est_idx, mult_idx = variant_to_estimator(vi)
			assign[vi] = (est_idx, mult_idx)
		groups = [(ESTIMATOR_SHORT_NAMES[i], ESTIMATOR_BASE_COLORS[label])
			for i, label in enumerate(ESTIMATOR_LABELS)]
		layout = AALayout(
			groups=groups, assign=assign, x_values=np.array(MULTIPLIERS),
			xlabel='kcat multiplier', invert_x=True, wt=wt,
			title_prefix='kcat_v2_estimator_sweep')

		render_all(plotOutDir, plotOutFileName, metadata, metrics, layout)


if __name__ == '__main__':
	Plot().cli()
