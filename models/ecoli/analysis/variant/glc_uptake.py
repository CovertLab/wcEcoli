"""
Template for variant analysis plots

@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 8/2/18
"""

from __future__ import absolute_import
from __future__ import division

import cPickle
from matplotlib import pyplot as plt
import numpy as np
import os

from models.ecoli.analysis import variantAnalysisPlot
from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from wholecell.analysis.analysis_tools import exportFigure
from wholecell.io.tablereader import TableReader
from wholecell.utils import constants, filepath, units


GLUCOSE_ID = 'GLC[p]'
FLUX_UNITS = units.mmol / units.g / units.h
MASS_UNITS = units.fg
GROWTH_UNITS = units.fg / units.s


class Plot(variantAnalysisPlot.VariantAnalysisPlot):
	def do_plot(self, inputDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		if not os.path.isdir(inputDir):
			raise Exception, 'inputDir does not currently exist as a directory'

		filepath.makedirs(plotOutDir)

		ap = AnalysisPaths(inputDir, variant_plot=True)
		variants = ap.get_variants()
		# TODO: select variants of interest
		n_variants = len(variants)

		with open(os.path.join(inputDir, 'kb', constants.SERIALIZED_FIT1_FILENAME), 'rb') as f:
			sim_data = cPickle.load(f)

		all_yields = []
		time_series = {}
		for variant in variants:
			yields = []
			trace = {}

			for sim_dir in ap.get_cells(variant=[variant]):
				sim_out_dir = os.path.join(sim_dir, 'simOut')

				# Listeners used
				fba_reader = TableReader(os.path.join(sim_out_dir, 'FBAResults'))
				main_reader = TableReader(os.path.join(sim_out_dir, 'Main'))
				mass_reader = TableReader(os.path.join(sim_out_dir, 'Mass'))

				# Load data
				time = main_reader.readColumn('time') - main_reader.readAttribute('initialTime')
				time_step_sec = main_reader.readColumn('timeStepSec')

				external_fluxes = fba_reader.readColumn('externalExchangeFluxes')
				external_molecules = fba_reader.readAttribute('externalMoleculeIDs')

				dry_mass = MASS_UNITS * mass_reader.readColumn('dryMass')
				growth = GROWTH_UNITS * mass_reader.readColumn('growth') / time_step_sec

				glc_idx = external_molecules.index(GLUCOSE_ID)
				glc_flux = FLUX_UNITS * external_fluxes[:, glc_idx]
				glc_mw = sim_data.getter.getMass([GLUCOSE_ID])[0]
				glc_mass_flux = glc_flux * glc_mw * dry_mass
				glc_mass_yield = growth / -glc_mass_flux

				yields += [np.nanmean(glc_mass_yield.asNumber())]
				trace['time'] = time
				trace['yield'] = glc_mass_yield.asNumber()

			all_yields += [yields]
			time_series[variant] = trace

		all_yields = np.vstack(all_yields).T
		means = all_yields.mean(axis=0)
		std = all_yields.std(axis=0)
		n_cells = all_yields.shape[0]

		plt.figure()

		plt.errorbar(range(n_variants), means, yerr=std, fmt='o')
		# plt.plot(range(n_cells), all_yields, 'x')

		# for v, data in time_series.items():
		# 	plt.plot(data['time'], data['yield'])

		exportFigure(plt, plotOutDir, plotOutFileName, metadata)

		plt.close('all')


if __name__ == '__main__':
	Plot().cli()
