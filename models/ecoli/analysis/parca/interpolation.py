"""
Plots for interpolation functions.
"""

import os
import pickle

from matplotlib import pyplot as plt
from matplotlib import gridspec
import numpy as np
from scipy import stats

from models.ecoli.analysis import parcaAnalysisPlot
from wholecell.analysis.analysis_tools import exportFigure
from wholecell.utils import constants, units


def get_raw(data, x_col, y_col, factor=1):
	xs = []
	ys = []
	for row in data:
		xs.append(row[x_col].asNumber(units.min))
		y = row[y_col] * factor
		if units.hasUnit(y):
			y = y.asNumber()
		ys.append(y)
	return xs, ys


class Plot(parcaAnalysisPlot.ParcaAnalysisPlot):
	def do_plot(self, input_dir, plot_out_dir, plot_out_filename, sim_data_file, validation_data_file, metadata):
		with open(os.path.join(input_dir, constants.SERIALIZED_RAW_DATA), 'rb') as f:
			raw_data = pickle.load(f)
		with open(sim_data_file, 'rb') as f:
			sim_data = pickle.load(f)
		growth = sim_data.growth_rate_parameters
		mass = sim_data.mass

		# Mapping of functions that perform interpolation to raw data
		interpolation_functions = {
			(growth.get_fraction_active_ribosome, None):
				get_raw(raw_data.growth_rate_dependent_parameters, 'doublingTime',
					'fractionActiveRibosome'),
			(growth.get_fraction_active_rnap, None):
				get_raw(raw_data.growth_rate_dependent_parameters, 'doublingTime',
					'fractionActiveRnap'),
			(growth.get_ppGpp_conc, None):
				get_raw(raw_data.growth_rate_dependent_parameters, 'doublingTime',
					'ppGpp_conc', factor=growth._per_dry_mass_to_per_volume),
			(growth.get_ribosome_elongation_rate, None):
				get_raw(raw_data.growth_rate_dependent_parameters, 'doublingTime',
					'ribosomeElongationRate'),
			(growth.get_rnap_elongation_rate, None):
				get_raw(raw_data.growth_rate_dependent_parameters, 'doublingTime',
					'rnaPolymeraseElongationRate'),
			(mass.get_dna_critical_mass, None): None,
			(mass.get_avg_cell_dry_mass, None):
				get_raw(raw_data.dry_mass_composition, 'doublingTime',
					'averageDryMass'),
			}

		# Interpolation functions that return values in a dictionary
		interpolation_functions.update({
			(mass.get_mass_fractions, fraction):
				get_raw(raw_data.dry_mass_composition, 'doublingTime',
					'{}MassFraction'.format(fraction))
			for fraction in mass.get_mass_fractions(45 * units.min)
			})
		interpolation_functions.update({
			(mass.get_component_masses, fraction): None
			for fraction in mass.get_component_masses(45 * units.min)
			})

		# TODO: handle getTrnaDistribution and all 86 outputs from 'molar_ratio_to_16SrRNA'
		# mass.getTrnaDistribution(45.*units.min)['molar_ratio_to_16SrRNA']

		# Doubling times to show on plot (extended range including all conditions)
		doubling_times = np.unique([
			dt.asNumber(units.min)
			for dt in sim_data.condition_to_doubling_time.values()
			])
		doubling_time_range = np.arange(0.5 * doubling_times.min(), 1.2 * doubling_times.max())

		# Create Plot
		plt.figure(figsize=(20, 20))
		n_plots = len(interpolation_functions)
		cols = 5
		gs = gridspec.GridSpec(int(np.ceil(n_plots / cols)), cols)

		funs = {
			'none': lambda x: x,
			'sqrt': lambda x: np.sqrt(x),
			'exp': lambda x: np.exp(x),
			'log': lambda x: np.log(x),
			'log2': lambda x: np.log(x**2),
			'logsqrt': lambda x: np.log(np.sqrt(x)),
			'2': lambda x: x**2,
			'3': lambda x: x**3,
			'1/sqrt': lambda x: 1/np.sqrt(x),
			'1/x': lambda x: 1/x,
			'1/x2': lambda x: 1/x**2,
			}
		inverse = {
			'none': lambda x: x,
			'sqrt': lambda x: x**2,
			'exp': lambda x: np.log(x),
			'log': lambda x: np.exp(x),
			'log2': lambda x: np.sqrt(np.exp(x)),
			'logsqrt': lambda x: np.exp(x)**2,
			'2': lambda x: np.sqrt(x),
			'3': lambda x: x**(1/3),
			'1/sqrt': lambda x: (1/x)**2,
			'1/x': lambda x: 1/x,
			'1/x2': lambda x: np.sqrt(1/x),
			}
		for i, ((fun, key), data) in enumerate(interpolation_functions.items()):
			ax = plt.subplot(gs[i // cols, i % cols])

			# Get interpolation values and handle units
			values = []
			unit = None
			for dt in units.min * doubling_time_range:
				# noinspection PyBroadException
				try:
					value = fun(dt)
					if key:
						value = value[key]

					if units.hasUnit(value):
						unit = str(units.getUnit(value)).strip('1 []')
						value = value.asNumber()
				except Exception as e:
					value = np.nan
				values.append(value)
			y_interp = np.array(values)

			# Create y labels
			label = fun.__name__ + ('\n{}'.format(key) if key else '')
			if unit:
				label += '\n({})'.format(unit)

			# Plot data
			ax.plot(doubling_time_range, y_interp)
			if data:
				slope = 1
				intercept = 0
				rval = 0
				print(fun.__name__)
				x = np.array(data[0])
				y = np.array(data[1])
				for xname, fx in funs.items():
					for yname, fy in funs.items():
						result = stats.linregress(fx(x), fy(y))
						if np.abs(result.rvalue) > rval:
							rval = np.abs(result.rvalue)
							xtransform = xname
							ytransform = yname
							slope = result.slope
							intercept = result.intercept
						if np.abs(result.rvalue) > 0.998:
							print('\t{} {}: {:.3f} {:.1e}'.format(xname, yname, result.rvalue, result.pvalue))
							ax.plot(doubling_time_range, inverse[yname](funs[xname](doubling_time_range) * result.slope + result.intercept), alpha=0.3)
				ax.plot(x, y, 'or')
				print(xtransform, ytransform)
				ax.plot(doubling_time_range, inverse[ytransform](funs[xtransform](doubling_time_range) * slope + intercept), '--')

			for dt in doubling_times:
				ax.axvline(dt, linestyle='--', color='k', linewidth=0.5)

			# Formatting
			ax.spines['top'].set_visible(False)
			ax.spines['right'].set_visible(False)
			ax.set_ylabel(label)
			ax.set_xlim([np.min(doubling_time_range), np.max(doubling_time_range)])

		## Save figure
		plt.tight_layout()
		exportFigure(plt, plot_out_dir, plot_out_filename, metadata)
		plt.close('all')


if __name__ == "__main__":
	Plot().cli()
