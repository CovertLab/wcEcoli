"""
Violin plots to compare predicted k_cats of the new list of disabled constraints with the baseline list.

@organization: Covert Lab, Department of Bioengineering, Stanford University
"""

from __future__ import absolute_import, division, print_function

import cPickle
from matplotlib import pyplot as plt
import numpy as np
import os
import re
import csv

from models.ecoli.analysis import variantAnalysisPlot
from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from wholecell.analysis.analysis_tools import exportFigure
from reconstruction.spreadsheets import JsonReader
from wholecell.io.tablereader import TableReader
from wholecell.utils import filepath
from wholecell.utils import units

from models.ecoli.sim.variants.kinetic_constraints_factorial_experiments import get_disabled_constraints

from models.ecoli.processes.metabolism import COUNTS_UNITS, VOLUME_UNITS, TIME_UNITS

# additional disabled constraints that are to be compared to baseline
ADDITIONAL_DISABLED_CONSTRAINTS = [
	'SUCCINATE-DEHYDROGENASE-UBIQUINONE-RXN-SUC/UBIQUINONE-8//FUM/CPD-9956.31.',
	'NADH-DEHYDROG-A-RXN-NADH/UBIQUINONE-8/PROTON//NAD/CPD-9956/PROTON.46. (reverse)',
	'INORGPYROPHOSPHAT-RXN[CCO-CYTOSOL]-PPI/WATER//Pi/PROTON.34.',
	'GLUTATHIONE-REDUCT-NADPH-RXN',
	]

REACTIONS = [
	'ISOCITDEH-RXN',
	'SUCCINATE-DEHYDROGENASE-UBIQUINONE-RXN-SUC/UBIQUINONE-8//FUM/CPD-9956.31.',
	'R601-RXN-FUM/REDUCED-MENAQUINONE//SUC/CPD-9728.38.',
	'NADH-DEHYDROG-A-RXN-NADH/UBIQUINONE-8/PROTON//NAD/CPD-9956/PROTON.46. (reverse)',
	'PSERTRANSAM-RXN',
	'GLUTATHIONE-REDUCT-NADPH-RXN',
	'GLYOXYLATE-REDUCTASE-NADP+-RXN__CPLX0-235',
	'INORGPYROPHOSPHAT-RXN[CCO-CYTOSOL]-PPI/WATER//Pi/PROTON.34.',
	'XANPRIBOSYLTRAN-RXN',
	]

OLD_MEASUREMENTS = {
	'ISOCITDEH-RXN': {
		'measurements': [106.3, 88.1], 'temps': [40, 40]},
	'SUCCINATE-DEHYDROGENASE-UBIQUINONE-RXN-SUC/UBIQUINONE-8//FUM/CPD-9956.31.': {
		'measurements': [78, 110, 24, 85], 'temps': [30, 30, 30, 30]},
	'R601-RXN-FUM/REDUCED-MENAQUINONE//SUC/CPD-9728.38.': {
		'measurements': [1128, 177, 250, 230], 'temps': [38, 30, 30, 30]},
	'NADH-DEHYDROG-A-RXN-NADH/UBIQUINONE-8/PROTON//NAD/CPD-9956/PROTON.46. (reverse)': {
		'measurements': [26], 'temps': [30]},
	'PSERTRANSAM-RXN': {
		'measurements': [1.75], 'temps': [37]},
	'GLUTATHIONE-REDUCT-NADPH-RXN': {
		'measurements': [733.3], 'temps': [30]},
	'GLYOXYLATE-REDUCTASE-NADP+-RXN__CPLX0-235': {
		'measurements': [203], 'temps': [25]},
	'INORGPYROPHOSPHAT-RXN[CCO-CYTOSOL]-PPI/WATER//Pi/PROTON.34.': {
		'measurements': [187, 390, 390], 'temps': [25, 25, 25]},
	'XANPRIBOSYLTRAN-RXN': {
		'measurements': [150], 'temps': [25]},
	}

NEW_MEASUREMENTS = {
	'ISOCITDEH-RXN': {
		'measurements': [], 'temps': []},
	'SUCCINATE-DEHYDROGENASE-UBIQUINONE-RXN-SUC/UBIQUINONE-8//FUM/CPD-9956.31.': {
		'measurements': [], 'temps': []},
	'R601-RXN-FUM/REDUCED-MENAQUINONE//SUC/CPD-9728.38.': {
		'measurements': [], 'temps': []},
	'NADH-DEHYDROG-A-RXN-NADH/UBIQUINONE-8/PROTON//NAD/CPD-9956/PROTON.46. (reverse)': {
		'measurements': [600], 'temps': [30]},
	'PSERTRANSAM-RXN': {
		'measurements': [], 'temps': []},
	'GLUTATHIONE-REDUCT-NADPH-RXN': {
		'measurements': [], 'temps': []},
	'GLYOXYLATE-REDUCTASE-NADP+-RXN__CPLX0-235': {
		'measurements': [], 'temps': []},
	'INORGPYROPHOSPHAT-RXN[CCO-CYTOSOL]-PPI/WATER//Pi/PROTON.34.': {
		'measurements': [42], 'temps': [25]},
	'XANPRIBOSYLTRAN-RXN': {
		'measurements': [], 'temps': []},
	}

CSV_DIALECT = csv.excel_tab
REACTIONS_FILE = os.path.join("reconstruction", "ecoli", "flat", "reactions.tsv")

# ignore data from first time steps
START_TIME_STEP = 2

def set_ticks(ax, labels):
	ax.spines['top'].set_visible(False)
	ax.spines['right'].set_visible(False)
	ax.get_xaxis().set_tick_params(direction='out')
	ax.xaxis.set_ticks_position('bottom')
	ax.set_xticks(np.arange(1, len(labels) + 1))
	ax.set_xticklabels(labels)
	ax.set_xlim(0.25, len(labels) + 0.75)


class Plot(variantAnalysisPlot.VariantAnalysisPlot):
	def do_plot(self, inputDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		if not os.path.isdir(inputDir):
			raise Exception, 'inputDir does not currently exist as a directory'

		filepath.makedirs(plotOutDir)

		ap = AnalysisPaths(inputDir, variant_plot=True)
		variants = ap.get_variants()

		# scan all variants to find variant indexes for comparison
		old_variant = None
		new_variant = None
		for v, variant in enumerate(variants):
			disable_constraints, additional_disabled = get_disabled_constraints(variant)
			if not additional_disabled:
				old_variant = variant
			if set(ADDITIONAL_DISABLED_CONSTRAINTS) == set(additional_disabled):
				new_variant = variant
		compared_variants = [new_variant, old_variant]

		# make dict of reactions
		all_reactions = {}
		with open(REACTIONS_FILE, 'rU') as csvfile:
			reader = JsonReader(csvfile, dialect=CSV_DIALECT)
			for row in reader:
				reaction_id = row["reaction id"]
				stoichiometry = row["stoichiometry"]
				reversible = row["is reversible"]
				catalyzed = row["catalyzed by"]
				all_reactions[reaction_id] = {
					"stoichiometry": stoichiometry,
					"is reversible": reversible,
					"catalyzed by": catalyzed,
				}

		# make map of reaction names to longer reaction names used by wcEcoli
		reaction_id_map = {}
		for reaction_id in REACTIONS:
			for reaction_id2, specs in all_reactions.iteritems():
				if reaction_id2 in reaction_id:
					reaction_id_map[reaction_id] = reaction_id2

		# get enzyme ids for all reactions, used to get concentrations
		reaction_enzymes = {}
		for reaction_id in REACTIONS:
			for reaction_id2, specs in all_reactions.iteritems():
				if reaction_id2 in reaction_id:
					reaction_enzymes[reaction_id] = all_reactions[reaction_id2]["catalyzed by"]

		# initialize dictionaries for fluxes and concentrations
		reaction_fluxes = {variant: {reaction_id: [] for reaction_id in REACTIONS}
			for variant in compared_variants}
		enzyme_concentrations = {variant: {} for variant in compared_variants}

		for variant in compared_variants:
			with open(ap.get_variant_kb(variant), 'rb') as f:
				sim_data = cPickle.load(f)

			cellDensity = sim_data.constants.cellDensity
			nAvogadro = sim_data.constants.nAvogadro

			for sim_dir in ap.get_cells(variant=[variant]):
				simOutDir = os.path.join(sim_dir, "simOut")

				# Listeners used
				try:
					massListener = TableReader(os.path.join(simOutDir, "Mass"))
					fbaResults = TableReader(os.path.join(simOutDir, "FBAResults"))
					bulkMolecules = TableReader(os.path.join(simOutDir, "BulkMolecules"))
				except Exception as e:
					print(e)
					continue

				# read from mass listener
				cellMass = massListener.readColumn("cellMass")[START_TIME_STEP:]  # units.fg
				cell_volume = np.array([mass * units.fg / cellDensity for mass in cellMass])
				counts_to_millimolar = np.array([1 / (vol.asNumber(VOLUME_UNITS) * nAvogadro.asNumber(1/COUNTS_UNITS))
					for vol in cell_volume])

				# read from FBA listener
				reactionIDs = np.array(fbaResults.readAttribute("reactionIDs"))
				reactionFluxes = (COUNTS_UNITS / TIME_UNITS) * fbaResults.readColumn("reactionFluxes")[START_TIME_STEP:,:]
				reaction_flux_dict = dict(zip(reactionIDs, reactionFluxes.asNumber(COUNTS_UNITS / TIME_UNITS).T))

				# append values to reaction_fluxes
				for reaction_id in REACTIONS:
					reaction_fluxes[variant][reaction_id].extend(list(reaction_flux_dict[reaction_id]))

				# get concentrations of enzymes
				enzymes = [item + compartment for sublist in reaction_enzymes.values()
					for item in sublist for compartment in ['[i]', '[c]', '[p]']]

				# read from bulkMolecules reader
				molecule_ids = bulkMolecules.readAttribute("objectNames")
				enzyme_ids = np.array([enzymeId for enzymeId in enzymes if enzymeId in molecule_ids])
				enzyme_indexes = np.array([molecule_ids.index(enzymeId)
					for enzymeId in enzymes if enzymeId in molecule_ids], np.int)
				molecule_counts = bulkMolecules.readColumn("counts")[START_TIME_STEP:, enzyme_indexes]

				# convert to millimolar concentrations
				concentrations = counts_to_millimolar * molecule_counts.T

				# put concentrations into a dict
				enzyme_concentrations_dict = dict(zip(enzyme_ids, concentrations))

				for enzyme_id, conc_time_series in enzyme_concentrations_dict.iteritems():
					enzyme_id_no_location = re.sub("[\(\[].*?[\)\]]", "", enzyme_id)
					if enzyme_id_no_location in enzyme_concentrations[variant]:
						enzyme_concentrations[variant][enzyme_id_no_location].extend(list(conc_time_series))
					else:
						enzyme_concentrations[variant][enzyme_id_no_location] = list(conc_time_series)

		### Make figure ###
		cols = 1
		rows = len(REACTIONS)
		fig = plt.figure(figsize=(cols * 5, rows * 3))

		# go through each reaction to show predicted k_cat distribution for the
		# new and old variant, and experimental measurements
		for reaction_idx, reaction_id in enumerate(REACTIONS):

			# reaction_id_2 can be used by all_reactions
			reaction_id_2 = reaction_id_map[reaction_id]
			reaction_specs = all_reactions[reaction_id_2]

			enzyme_id = reaction_specs['catalyzed by']

			if 'CPLX0-235' in enzyme_id:
				enzyme_id = ['CPLX0-235']  # remove 'G6539-MONOMER' from GLYOXYLATE-REDUCTASE-NADP+-RXN

			# old measurements
			reaction_measurements = OLD_MEASUREMENTS[reaction_id]
			measurements = reaction_measurements['measurements']
			temps = reaction_measurements['temps']
			adjusted_measurements = np.array([2**((37.-t)/10.)*m for (m, t) in zip(measurements, temps)])

			# new measurements
			reaction_measurements = NEW_MEASUREMENTS[reaction_id]
			measurements = reaction_measurements['measurements']
			temps = reaction_measurements['temps']
			new_adjusted_measurements = np.array([2**((37.-t)/10.)*m for (m, t) in zip(measurements, temps)])

			if reaction_id == 'GLUTATHIONE-REDUCT-NADPH-RXN':
				saturated_fraction = 0.238
				# get effective kcat
				new_adjusted_measurements = adjusted_measurements * saturated_fraction

			# Initialize subplots
			ax = plt.subplot(rows, cols, reaction_idx+1)

			# calculate the reaction's k_cat distribution for each compared variant
			k_cat_distribution = {}
			for variant in compared_variants:
				## Get data
				rxn_fluxes = np.array(reaction_fluxes[variant][reaction_id]) # mmol / L / s
				enzyme_concs = np.array(enzyme_concentrations[variant][enzyme_id[0]])  # mmol / L

				# calculate k_cats
				k_cats = rxn_fluxes / enzyme_concs
				k_cat_distribution[variant] = k_cats

			data = [k_cat_distribution[variant] for variant in compared_variants]

			# plot
			violin_pos = [3, 1]	# position of violin plot [old, new]
			measure_pos = 2  	# position of measurements
			ax.violinplot(data, violin_pos, widths=1.0, showmeans=True, showextrema=True, showmedians=False)
			ax.scatter(np.full_like(adjusted_measurements, measure_pos), adjusted_measurements, marker='*', color='Black')
			ax.scatter(np.full_like(new_adjusted_measurements, measure_pos+0.05), new_adjusted_measurements, marker='*', color='Red')

			# format
			rxn_id_length = 25
			text_reaction_id = ('reaction: %s' % reaction_id[:rxn_id_length])
			labels = ['old', 'measured', 'new']
			ax.set_title(text_reaction_id, fontsize=12)
			ax.set_ylabel('$k_{cat}$', fontsize=12)
			set_ticks(ax, labels)
			ax.set_yscale('log')

		### Create Plot ###
		plt.tight_layout()
		plt.subplots_adjust(hspace=1.0)
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close('all')


if __name__ == "__main__":
	Plot().cli()