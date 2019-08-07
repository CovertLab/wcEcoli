"""
Violin plots to compare predicted k_cats of the new list of disabled constraints with the baseline list.

@organization: Covert Lab, Department of Bioengineering, Stanford University
"""

from __future__ import absolute_import, division, print_function

import cPickle
import os

from matplotlib import pyplot as plt
import numpy as np

from models.ecoli.analysis import variantAnalysisPlot
from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from models.ecoli.processes.metabolism import COUNTS_UNITS, VOLUME_UNITS, TIME_UNITS
from models.ecoli.sim.variants.kinetic_constraints_factorial_experiments import get_disabled_constraints
from wholecell.analysis.analysis_tools import exportFigure
from wholecell.io.tablereader import TableReader
from wholecell.utils import constants, filepath, units


# additional disabled constraints that are to be compared to baseline
ADDITIONAL_DISABLED_CONSTRAINTS = [
	'SUCCINATE-DEHYDROGENASE-UBIQUINONE-RXN-SUC/UBIQUINONE-8//FUM/CPD-9956.31.',
	'NADH-DEHYDROG-A-RXN-NADH/UBIQUINONE-8/PROTON//NAD/CPD-9956/PROTON.46. (reverse)',
	'INORGPYROPHOSPHAT-RXN[CCO-CYTOSOL]-PPI/WATER//Pi/PROTON.34.',
	'GLUTATHIONE-REDUCT-NADPH-RXN',
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
	'RXN-11832': {
		'measurements': [103], 'temps': [30]},
	'CYTDEAM-RXN': {
		'measurements': [185, 165, 49.68 ,132 ,45], 'temps': [25, 37, 25, 30, 25]},
	}

NEW_MEASUREMENTS = {
	'NADH-DEHYDROG-A-RXN-NADH/UBIQUINONE-8/PROTON//NAD/CPD-9956/PROTON.46. (reverse)': {
		'measurements': [600], 'temps': [30]},
	'INORGPYROPHOSPHAT-RXN[CCO-CYTOSOL]-PPI/WATER//Pi/PROTON.34.': {
		'measurements': [42], 'temps': [25]},
	}

REACTIONS = sorted(OLD_MEASUREMENTS.keys())

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
			if additional_disabled is None:
				old_variant = variant
			elif set(ADDITIONAL_DISABLED_CONSTRAINTS) == set(additional_disabled):
				new_variant = variant

		# if the baseline variant or the new variant are missing, stop plotting
		if (old_variant is None) or (new_variant is None):
			print('Variant simulations missing!')
			return

		compared_variants = [old_variant, new_variant]

		# Load sim_data
		with open(os.path.join(inputDir, 'kb', constants.SERIALIZED_FIT1_FILENAME), 'rb') as f:
			sim_data = cPickle.load(f)

		# get reactions from sim_data
		reactionStoich = sim_data.process.metabolism.reactionStoich
		reactionCatalysts = sim_data.process.metabolism.reactionCatalysts

		reaction_enzymes = {}
		for reaction_id in REACTIONS:
			for reaction_id2, specs in reactionStoich.iteritems():
				if reaction_id2 in reaction_id:
					reaction_enzymes[reaction_id] = reactionCatalysts[reaction_id2]
		enzymes = [mol_id for rxn in reaction_enzymes.values() for mol_id in rxn]

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

				# read from bulkMolecules reader
				molecule_ids = bulkMolecules.readAttribute("objectNames")
				enzyme_ids = np.array([enzymeId for enzymeId in enzymes if enzymeId in molecule_ids])
				enzyme_indexes = np.array([molecule_ids.index(enzymeId)
					for enzymeId in enzymes if enzymeId in molecule_ids], np.int)
				molecule_counts = bulkMolecules.readColumn("counts")[START_TIME_STEP:, enzyme_indexes]

				# convert to millimolar concentrations
				concentrations = counts_to_millimolar * molecule_counts.T

				# add concentration timeseries to enzyme_concentrations
				for enzyme_id, conc_time_series in zip(enzyme_ids, concentrations):
					if enzyme_id in enzyme_concentrations[variant]:
						enzyme_concentrations[variant][enzyme_id].extend(list(conc_time_series))
					else:
						enzyme_concentrations[variant][enzyme_id] = list(conc_time_series)

		### Make figure ###
		cols = 1
		rows = len(REACTIONS)
		fig = plt.figure(figsize=(cols * 5, rows * 3))

		# go through each reaction to show predicted k_cat distribution for the
		# new and old variant, and experimental measurements
		for reaction_idx, reaction_id in enumerate(REACTIONS):
			enzyme_id = reaction_enzymes[reaction_id]

			# old measurements
			reaction_measurements = OLD_MEASUREMENTS[reaction_id]
			measurements = reaction_measurements['measurements']
			temps = reaction_measurements['temps']
			adjusted_measurements = np.array([2**((37.-t)/10.)*m for (m, t) in zip(measurements, temps)])

			# new measurements
			reaction_measurements = NEW_MEASUREMENTS.get(reaction_id, {})
			measurements = reaction_measurements.get('measurements', [])
			temps = reaction_measurements.get('temps', [])
			new_adjusted_measurements = np.array([2**((37.-t)/10.)*m for (m, t) in zip(measurements, temps)])

			# get effective kcat for GLUTATHIONE-REDUCT
			if reaction_id == 'GLUTATHIONE-REDUCT-NADPH-RXN':
				# saturated_fraction calculated from Smirnova, et al. (2005). "Effects of cystine and
				# hydrogen peroxideon glutathione status and expression of	antioxidant	genes in Escherichia coli"
				# Oxidized glutathione (GSSG in table 2) gives ~19 uM concentration (with 0.3 dry fraction and 1.1 g/mL density)
				# With 61 uM Km for this reaction, that gives a saturated fraction of 0.238
				saturated_fraction = 0.238
				new_adjusted_measurements = adjusted_measurements * saturated_fraction

			# Initialize subplots
			ax = plt.subplot(rows, cols, reaction_idx+1)

			# calculate the reaction's k_cat distribution for each compared variant
			k_cat_distribution = {}
			for variant in compared_variants:
				## Get data
				rxn_fluxes = np.array(reaction_fluxes[variant][reaction_id]) 			# mmol / L / s
				enzyme_concs = np.array(enzyme_concentrations[variant][enzyme_id[0]])  	# mmol / L

				# calculate k_cats, remove zeros, save to this variant's distribution
				k_cats = rxn_fluxes / enzyme_concs
				k_cats = k_cats[k_cats > 0]
				k_cat_distribution[variant] = k_cats

			data = [k_cat_distribution[old_variant], k_cat_distribution[new_variant]]

			# plot
			violin_pos = [1, 3]	# position of violin plots [old, new]
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