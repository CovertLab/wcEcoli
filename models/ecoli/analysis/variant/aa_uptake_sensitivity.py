"""
Compare cell cycle times and growth rates across variants.  Useful as validation
with the add_one_aa variant and can also be used to compare variants in
remove_one_aa variant.
"""

import pickle

from matplotlib import colors, pyplot as plt
import numpy as np
from scipy.stats import pearsonr

from models.ecoli.analysis import variantAnalysisPlot
from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from models.ecoli.sim.variants.aa_uptake_sensitivity import get_aa_index, get_factor
from wholecell.analysis.analysis_tools import exportFigure, read_stacked_columns
from wholecell.utils import units


CONTROL_LABEL = 'L-SELENOCYSTEINE'  # control because SEL is already included for uptake in minimal media
GLC_ID = 'GLC[p]'
FLUX_UNITS = units.mmol / units.g / units.h
MASS_UNITS = units.fg
GROWTH_UNITS = MASS_UNITS / units.s
AXIS_LIMITS = [0.5, 1.5]


def plot_validation(mean, std, labels, val_rates, val_std, val_aa_ids, label, text_highlight=None):
	# Normalize simulation data by the control condition
	rate_mapping = {label: rate for label, rate in zip(labels, mean)}
	std_mapping = {label: std for label, std in zip(labels, std)}
	wcm_control = rate_mapping.get(CONTROL_LABEL, 1)
	wcm_normalized_growth_rates = np.array([
		rate_mapping.get(aa, 0) / wcm_control
		for aa in val_aa_ids
		])
	wcm_normalized_std = np.array([
		std_mapping.get(aa, 0) / wcm_control
		for aa in val_aa_ids
		])

	# Statistics
	r, p = pearsonr(val_rates, wcm_normalized_growth_rates)
	n = len(val_rates)

	plt.errorbar(val_rates, wcm_normalized_growth_rates,
		xerr=val_std, yerr=wcm_normalized_std, fmt='o', alpha=0.5,
		label=f'{label} r={r:.2f} (p={p:.2g}, n={n})')

	if text_highlight:
		for aa, x, y in zip(val_aa_ids, val_rates, wcm_normalized_growth_rates):
			color = 'r' if text_highlight.get(aa, False) else 'k'
			if y < AXIS_LIMITS[0]:
				y = AXIS_LIMITS[0]
			elif y > AXIS_LIMITS[1]:
				y = AXIS_LIMITS[1]
			plt.text(x, 0.01 + y, aa, ha='center', fontsize=6, color=color)


class Plot(variantAnalysisPlot.VariantAnalysisPlot):
	def do_plot(self, inputDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		ap = AnalysisPaths(inputDir, variant_plot=True)
		variants = ap.get_variants()

		with open(simDataFile, 'rb') as f:
			sim_data = pickle.load(f)
		with open(validationDataFile, 'rb') as f:
			validation_data = pickle.load(f)

		aa_ids = sim_data.molecule_groups.amino_acids
		aa_idx = {aa[:-3]: i for i, aa in enumerate(aa_ids)}
		aa_labels = ['Control' if aa == CONTROL_LABEL + '[c]' else aa[:-3] for aa in aa_ids]

		growth_rates = {}
		control_rates = []
		for variant in variants:
			cells = ap.get_cells(variant=[variant])
			growth_rate = read_stacked_columns(cells, 'Mass', 'instantaneous_growth_rate',
				remove_first=True, ignore_exception=True).mean()
			aa_id = aa_ids[get_aa_index(variant)][:-3]
			if aa_id == CONTROL_LABEL:
				control_rates.append(growth_rate)
			growth_rates[(aa_id, get_factor(variant))] = growth_rate
		control_rate = np.mean(control_rates) if control_rates else 1

		# Load validation growth rates
		all_aa_ids = {aa[:-3] for aa in aa_ids}
		val_control = validation_data.amino_acid_growth_rates['minimal']['mean']
		val_aa_ids = []
		val_normalized_growth_rates = []
		# val_normalized_std = []
		for media, rates in validation_data.amino_acid_growth_rates.items():
			aa_id = media.split('_')[-1]
			if aa_id in all_aa_ids:
				val_aa_ids.append(aa_id)
				val_normalized_growth_rates.append(units.strip_empty_units(rates['mean'] / val_control))
				# val_normalized_std.append(units.strip_empty_units(rates['std'] / val_control))
		# val_normalized_growth_rates = np.array(val_normalized_growth_rates)
		# val_normalized_std = np.array(val_normalized_std)
		val_x = [aa_idx[aa] for aa in val_aa_ids]

		# Create plots
		plt.figure()

		x = []
		y = []
		c = []
		edges = []
		for (aa, factor), rate in growth_rates.items():
			color = 'k'
			x.append(aa_idx[aa])
			y.append(rate / control_rate)
			if factor == 0:
				color = 'gray'
				c.append(1)
			else:
				c.append(factor)
			edges.append(color)

		scatter = plt.scatter(x, y, c=c, cmap='RdBu', edgecolors=edges, alpha=0.8, norm=colors.LogNorm())
		plt.scatter(val_x, val_normalized_growth_rates, marker='_', c='k')
		plt.colorbar(scatter)

		plt.xticks(range(len(aa_labels)), aa_labels, rotation=45, fontsize=6, ha='right')

		# Plot formatting
		# ax.tick_params(axis='x', labelsize=6)
		# ax.axhline(1, linestyle='--', linewidth=0.5, color='k', alpha=0.5)
		# ax.axvline(1, linestyle='--', linewidth=0.5, color='k', alpha=0.5)
		plt.xlabel('Validation growth rate\n(Normalized to minimal media)')
		plt.ylabel('Simulation data\n(Normalized to minimal media)')

		# Plot y=x diagonal
		# x_min, x_max = plt.xlim()
		# y_min, y_max = plt.ylim()
		# min_rate = min(x_min, y_min)
		# max_rate = max(x_max, y_max)
		# plt.plot([min_rate, max_rate], [min_rate, max_rate], '--k')

		# # Limit axes to reasonable range
		# plt.xlim(AXIS_LIMITS)
		# plt.ylim(AXIS_LIMITS)

		plt.tight_layout()
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close('all')


if __name__ == "__main__":
	Plot().cli()
