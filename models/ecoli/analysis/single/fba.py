"""
Troubleshooting metabolism at certain time steps.

@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 6/1/20
"""

from __future__ import absolute_import, division, print_function

import cPickle
import os

from matplotlib import pyplot as plt
import numpy as np

from models.ecoli.analysis import singleAnalysisPlot
from models.ecoli.processes.metabolism import CONC_UNITS, CONVERSION_UNITS, FluxBalanceAnalysisModel, GDCW_BASIS
from wholecell.analysis.analysis_tools import exportFigure
from wholecell.io.tablereader import TableReader
from wholecell.utils import units


class Plot(singleAnalysisPlot.SingleAnalysisPlot):
	def do_plot(self, simOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		with open(simDataFile, 'rb') as f:
			sim_data = cPickle.load(f)
		exchange_molecules = np.array(sim_data.external_state.all_external_exchange_molecules)

		# Listeners used
		fba_reader = TableReader(os.path.join(simOutDir, 'FBAResults'))
		kinetics_reader = TableReader(os.path.join(simOutDir, 'EnzymeKinetics'))
		main_reader = TableReader(os.path.join(simOutDir, 'Main'))

		# Load data
		catalyst_ids = fba_reader.readAttribute('catalyst_ids')
		update_molecules = fba_reader.readAttribute('conc_update_molecules')
		media_ids = fba_reader.readColumn('media_id')[1:]
		all_conc_updates = fba_reader.readColumn('conc_updates')[1:, :]
		all_catalyst_counts = fba_reader.readColumn('catalyst_counts')[1:, :]
		all_translation_gtp = fba_reader.readColumn('translation_gtp')[1:]
		coefficients = CONVERSION_UNITS * fba_reader.readColumn('coefficient')[1:]
		all_unconstrained_molecules = fba_reader.readColumn('unconstrained_molecules')[1:, :]
		all_constrained_molecules = fba_reader.readColumn('constrained_molecules')[1:, :]
		all_uptake_constraints = fba_reader.readColumn('uptake_constraints')[1:, :]

		metabolite_ids = kinetics_reader.readAttribute('metaboliteNames')
		all_metabolite_counts = kinetics_reader.readColumn('metaboliteCountsInit')[1:, :]
		all_counts_to_molar = CONC_UNITS * kinetics_reader.readColumn('countsToMolar')[1:]

		time_steps = units.s * main_reader.readColumn('timeStepSec')[1:]

		# Model
		model = FluxBalanceAnalysisModel(sim_data)

		for (media_id, conc_updates, catalyst_counts, translation_gtp,
				coefficient, unconstrained_molecules, constrained_molecules,
				uptake_constraints, metabolite_counts, counts_to_molar,
				time_step) in zip(media_ids, all_conc_updates,
				all_catalyst_counts, all_translation_gtp, coefficients,
				all_unconstrained_molecules, all_constrained_molecules,
				all_uptake_constraints, all_metabolite_counts,
				all_counts_to_molar, time_steps):
			## Calculations
			conc_updates = dict(zip(update_molecules, conc_updates))
			catalyst_dict = dict(zip(catalyst_ids, catalyst_counts))
			metabolite_dict = dict(zip(metabolite_ids, metabolite_counts))
			unconstrained = set(exchange_molecules[unconstrained_molecules])
			constrained = {
				mol: uptake * GDCW_BASIS
				for mol, uptake, present
				in zip(exchange_molecules, uptake_constraints, constrained_molecules)
				if present
				}
			kinetic_enzyme_counts = np.array([catalyst_dict[e]
				for e in model.kinetic_constraint_enzymes])
			kinetic_substrate_counts = np.array([metabolite_dict[s]
				for s in model.kinetic_constraint_substrates])

			## Set molecule availability (internal and external)
			model.set_molecule_levels(metabolite_counts, counts_to_molar,
				coefficient, media_id, unconstrained, constrained, conc_updates)

			## Set reaction limits for maintenance and catalysts present
			model.set_reaction_bounds(catalyst_counts, counts_to_molar,
				coefficient, translation_gtp)

			## Constrain reactions based on targets
			targets, upper_targets, lower_targets = model.set_reaction_targets(kinetic_enzyme_counts,
				kinetic_substrate_counts, counts_to_molar, time_step)

		plt.figure()

		### Create Plot ###

		plt.tight_layout()
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close('all')


if __name__ == '__main__':
	Plot().cli()
