"""
Template for variant analysis plots
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 8/2/18
"""

from __future__ import absolute_import
from __future__ import division

import cPickle
from matplotlib import pyplot as plt
from scipy.stats import pearsonr
import numpy as np
import os
import re
import csv

from models.ecoli.analysis import variantAnalysisPlot
from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from wholecell.analysis.analysis_tools import exportFigure
from wholecell.io.tablereader import TableReader
from wholecell.utils import filepath
from wholecell.utils import units

from models.ecoli.processes.metabolism import COUNTS_UNITS, VOLUME_UNITS, TIME_UNITS, MASS_UNITS

# COUNTS_UNITS = units.mmol
# VOLUME_UNITS = units.L
# MASS_UNITS = units.g
# TIME_UNITS = units.h

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

EXCHANGES = ['GLC[p]']

# ignore data from metabolism burnin period
BURN_IN_TIME = 1

class Plot(variantAnalysisPlot.VariantAnalysisPlot):
	def do_plot(self, inputDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		if not os.path.isdir(inputDir):
			raise Exception, 'inputDir does not currently exist as a directory'

		filepath.makedirs(plotOutDir)

		ap = AnalysisPaths(inputDir, variant_plot=True)
		variants = ap.get_variants()

		# reactions/exchanges of interest
		glc_id = 'GLC[p]'
		succ_id = 'SUCCINATE-DEHYDROGENASE-UBIQUINONE-RXN-SUC/UBIQUINONE-8//FUM/CPD-9956.31.'
		iso_id = 'ISOCITDEH-RXN'
		nadh_id = 'NADH-DEHYDROG-A-RXN-NADH/UBIQUINONE-8/PROTON//NAD/CPD-9956/PROTON.46. (reverse)'

		validation_data = cPickle.load(open(validationDataFile, "rb"))
		toyaReactions = validation_data.reactionFlux.toya2010fluxes["reactionID"]
		toyaFluxes = validation_data.reactionFlux.toya2010fluxes["reactionFlux"]
		toyaStdev = validation_data.reactionFlux.toya2010fluxes["reactionFluxStdev"]
		toyaFluxesDict = dict(zip(toyaReactions, toyaFluxes))
		toyaStdevDict = dict(zip(toyaReactions, toyaStdev))

		save_output_file = 'kinetic_constraint_comparison.csv'
		with open(save_output_file, 'w') as csvfile:
			fieldnames = [
				'variant',
				'disabled reactions',
				'succ avg sim flux (mmol/gDCW/h)',
				'succ target avg (mmol/gDCW/h)',
				'succ distance from toya (mmol/gDCW/h)',
				'nadh avg sim flux (mmol/gDCW/h)',
				'nadh target avg (mmol/gDCW/h)',
				'glc avg sim flux (mmol/gDCW/h)',
				'Pearson R Toya',
				]
			writer = csv.DictWriter(csvfile, delimiter=',',	fieldnames=fieldnames)
			writer.writeheader()

		for variant in variants:
			with open(ap.get_variant_kb(variant), 'rb') as f:
				sim_data = cPickle.load(f)

			cellDensity = sim_data.constants.cellDensity

			additional_disabled = sim_data.process.metabolism.additional_disabled

			# initialize kinetic flux comparison
			targetFluxList = []
			actualFluxList = []
			exchange_fluxes = {entry: [] for entry in EXCHANGES}
			reaction_fluxes = {entry: [] for entry in REACTIONS}

			modelFluxes = {}
			toyaOrder = []
			for rxn in toyaReactions:
				modelFluxes[rxn] = []
				toyaOrder.append(rxn)

			for sim_dir in ap.get_cells(variant=[variant]):
				simOutDir = os.path.join(sim_dir, "simOut")

				# Listeners used
				try:
					massListener = TableReader(os.path.join(simOutDir, "Mass"))
					cellMass = massListener.readColumn("cellMass")
					dryMass = massListener.readColumn("dryMass")
					massListener.close()
				except Exception as e:
					print(e)
					continue

				coefficient = dryMass / cellMass * cellDensity.asNumber(MASS_UNITS / VOLUME_UNITS)

				mainListener = TableReader(os.path.join(simOutDir, "Main"))
				time = mainListener.readColumn("time")

				# skip if no data
				if cellMass.shape is ():
					continue
				burnIn = time > BURN_IN_TIME
				burnIn[0] = False

				## Read from FBA listener
				fbaResults = TableReader(os.path.join(simOutDir, "FBAResults"))
				reactionIDs = np.array(fbaResults.readAttribute("reactionIDs"))
				reactionFluxes = (COUNTS_UNITS / MASS_UNITS / TIME_UNITS) * (fbaResults.readColumn("reactionFluxes").T / coefficient).T
				exFlux = fbaResults.readColumn("externalExchangeFluxes")
				exMolec = fbaResults.readAttribute("externalMoleculeIDs")
				fbaResults.close()

				# make dicts for exchange fluxes, targets, and reactions.
				exchange_flux_dict = dict(zip(exMolec, exFlux.T))
				reaction_flux_dict = dict(zip(reactionIDs, reactionFluxes.asNumber(units.mmol / units.g / units.h).T))

				## Append values for relevant reactions.
				# append to exchanges
				for entry in EXCHANGES:
					exchange_fluxes[entry].extend(list(exchange_flux_dict[entry]))
				# append to reaction fluxes
				for entry in REACTIONS:
					reaction_fluxes[entry].extend(list(reaction_flux_dict[entry]))

				## get all Toya reactions, and corresponding simulated fluxes.
				for toyaReaction in toyaReactions:
					fluxTimeCourse = []

					for rxn in reactionIDs:
						if re.findall(toyaReaction, rxn):
							reverse = 1
							if re.findall("(reverse)", rxn):
								reverse = -1
							if len(fluxTimeCourse):
								fluxTimeCourse += reverse * reactionFluxes[:, np.where(reactionIDs == rxn)]
							else:
								fluxTimeCourse = reverse * reactionFluxes[:, np.where(reactionIDs == rxn)]

					if len(fluxTimeCourse):
						modelFluxes[toyaReaction].append(np.mean(fluxTimeCourse).asNumber(units.mmol / units.g / units.h))

				## Read from EnzymeKinetics listener
				enzymeKineticsReader = TableReader(os.path.join(simOutDir, "EnzymeKinetics"))
				targetFluxes = (COUNTS_UNITS / MASS_UNITS / TIME_UNITS) * (enzymeKineticsReader.readColumn("targetFluxes").T / coefficient).T
				actualFluxes = (COUNTS_UNITS / MASS_UNITS / TIME_UNITS) * (enzymeKineticsReader.readColumn("actualFluxes").T / coefficient).T
				reactionConstraint = enzymeKineticsReader.readColumn("reactionConstraint")
				constrainedReactions = np.array(enzymeKineticsReader.readAttribute("constrainedReactions"))
				enzymeKineticsReader.close()

				targetFluxes = targetFluxes.asNumber(units.mmol / units.g / units.h)
				actualFluxes = actualFluxes.asNumber(units.mmol / units.g / units.h)
				targetAve = np.nanmean(targetFluxes[burnIn, :], axis=0)
				actualAve = np.nanmean(actualFluxes[burnIn, :], axis=0)

				# save relevant values
				if len(targetFluxList) == 0:
					targetFluxList = np.array([targetAve])
					actualFluxList = np.array([actualAve])
					reactionConstraintList = np.array(reactionConstraint[burnIn, :])
				else:
					targetFluxList = np.concatenate((targetFluxList, np.array([targetAve])), axis=0)
					actualFluxList = np.concatenate((actualFluxList, np.array([actualAve])), axis=0)
					reactionConstraintList = np.concatenate((reactionConstraintList, np.array(reactionConstraint[burnIn, :])), axis=0)

			## Get averages from FBA
			fba_exchange_flux_dict_avg = {entry: -1 * np.mean(values) for entry, values in exchange_fluxes.iteritems() if entry in EXCHANGES}
			fba_reaction_flux_dict_avg = {entry: np.mean(values) for entry, values in reaction_fluxes.iteritems() if entry in REACTIONS}

			## Flux comparison with Toya
			toyaVsReactionAve = []
			rxn_order = []
			for rxn, toyaFlux in toyaFluxesDict.iteritems():
				rxn_order.append(rxn)
				if rxn in modelFluxes:
					toyaVsReactionAve.append((np.mean(modelFluxes[rxn]), toyaFlux.asNumber(units.mmol / units.g / units.h), np.std(modelFluxes[rxn]), toyaStdevDict[rxn].asNumber(units.mmol / units.g / units.h)))

			outlier_indicies = np.zeros(len(toyaReactions), bool)
			outlier_indicies[rxn_order.index(succ_id)] = True
			outlier_indicies[rxn_order.index(iso_id)] = True

			toyaVsReactionAve = np.array(toyaVsReactionAve)
			rWithAll = pearsonr(toyaVsReactionAve[:,0], toyaVsReactionAve[:,1])

			succ_avg_sim_flux = toyaVsReactionAve[rxn_order.index(succ_id), 0]
			succ_toya_flux = toyaVsReactionAve[rxn_order.index(succ_id), 1]
			succ_difference = succ_avg_sim_flux - succ_toya_flux

			## Constrained Enzyme Kinetics
			constrainedReactions = list(constrainedReactions)
			constrain_rxn_indicies = np.zeros(len(constrainedReactions), bool)
			for rxn_id in REACTIONS:
				if rxn_id in constrainedReactions:
					constrain_rxn_indicies[
						constrainedReactions.index(rxn_id)] = True
					print(rxn_id + ' constrained')

			# determine average values across all cells
			targetAve = np.nanmean(targetFluxList, axis=0)
			actualAve = np.nanmean(actualFluxList, axis=0)

			# put in dicts
			target_flux = {}
			actual_flux = {}
			distance = {}
			for index, bool_value in enumerate(constrain_rxn_indicies):
				if bool_value:
					target_flux[constrainedReactions[index]] = targetAve[index]
					actual_flux[constrainedReactions[index]] = actualAve[index]
					distance[constrainedReactions[index]] = actualAve[index] - targetAve[index]

			# SUCC
			succ_avg_sim_flux = fba_reaction_flux_dict_avg[succ_id]
			if succ_id in target_flux:
				succ_target_avg = target_flux[succ_id]
			else:
				succ_target_avg = 0

			# NADH
			nadh_avg_sim_flux = fba_reaction_flux_dict_avg[nadh_id]
			if nadh_id in target_flux:
				nadh_target_avg = target_flux[nadh_id]
			else:
				nadh_target_avg = 0

			glc_avg_sim_flux = fba_exchange_flux_dict_avg[glc_id]

			# list of 4-character strings, for each disabled reaction
			additional_disabled_reduced = [reaction[:4] for reaction in additional_disabled]
			print('{}: var {} -- sim succ {}, target succ {}, succ dist toya {}, sim nadh {}, target nadh {}, avg glc flux {}, Pearson R Toya {},'.format(
				variant,
				str(additional_disabled_reduced),
				succ_avg_sim_flux,
				succ_target_avg,
				succ_difference,
				nadh_avg_sim_flux,
				nadh_target_avg,
				glc_avg_sim_flux,
				rWithAll,
				))

			# save reaction to csv
			with open(save_output_file, 'a') as csvfile:
				append_line = [
					variant,
					str(additional_disabled_reduced),
					succ_avg_sim_flux,
					succ_target_avg,
					succ_difference,
					nadh_avg_sim_flux,
					nadh_target_avg,
					glc_avg_sim_flux,
					rWithAll,
					]

				writer = csv.writer(csvfile, delimiter=',')
				writer.writerow(append_line)



if __name__ == "__main__":
	Plot().cli()