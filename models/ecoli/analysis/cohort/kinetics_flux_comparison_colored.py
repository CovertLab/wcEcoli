"""
Compare fluxes in simulation to target fluxes
"""

import io
import os
import pickle
import ast
import numpy as np
from matplotlib import pyplot as plt
from scipy.stats import pearsonr
import pandas as pd
from wholecell.io.tablereader import TableReader
from wholecell.io import tsv
from wholecell.utils import units
from wholecell.utils.sparkline import whitePadSparklineAxis

from models.ecoli.processes.metabolism import COUNTS_UNITS, VOLUME_UNITS, TIME_UNITS, MASS_UNITS
from wholecell.analysis.analysis_tools import exportFigure
from models.ecoli.analysis import cohortAnalysisPlot

# have the relevant tsv file paths:
os.chdir(os.path.expanduser('~/wcEcoli/'))
metabolic_reactions_tsv_path = "reconstruction/ecoli/flat/metabolic_reactions.tsv"
complexation_reactions_tsv_path = "reconstruction/ecoli/flat/complexation_reactions.tsv"


# ignore data from metabolism burnin period
BURN_IN_TIME = 1

IMPORTANT_MONOMERS = ["EG10022-MONOMER"]

class Plot(cohortAnalysisPlot.CohortAnalysisPlot):
	_suppress_numpy_warnings = True

	def find_complexes(self, stoich_list):
		# check if there happens to be a complex within the list of stoich:
		complex_appearances = []
		for substrate in stoich_list:
			if substrate in self.complex_IDs:
				complex_appearances.append(substrate)

		return complex_appearances

	def find_monomer_roots(self, stoich_list):
		# find the monomer roots of each complex:
		original_stoich_list = stoich_list
		# determine if any substrates are complexes too:
		substrate_complexes = self.find_complexes(original_stoich_list)
		substrate_complexes = list(substrate_complexes)
		hi = 5
		# remove the complexes that have been found:
		new_stoich_list = original_stoich_list
		for complex in substrate_complexes:
			new_stoich_list.remove(complex)

		monomers_caught = []
		substrate_complexes_caught = []
		substrate_complexes_updated = substrate_complexes
		while len(substrate_complexes_updated) != 0:
			new_substrates = []
			for complex in substrate_complexes_updated:
				rxn = self.complex_IDs_to_reactions.get(complex)
				# get the
				idx = self.complexation_reactions_df.index[
					self.complexation_reactions_df['id'] == rxn].tolist()
				if len(idx) != 1:
					print(
						f"More than one index for {rxn}: {idx}")  # hopefully do not get any of these
					monomers = "None"

				# find the substrates associated with the reaction:
				stoich_series = self.complexation_reactions_df["stoichiometry"][idx]
				stoich = ast.literal_eval(stoich_series.iloc[0])

				# get a list of the dictionary keys
				monomers_list = list(stoich.keys())
				# remove the complex itself from the keys:
				complex_ID = self.complex_reactions_to_IDs.get(rxn)
				monomers_list.remove(complex_ID)

				# put the new substrate list in the mix:
				new_substrates.append(monomers_list)

				# remove the current complex from the list of substrates:
				substrate_complexes_updated.remove(complex)
				substrate_complexes_caught.append(complex)

				# check if the new monomers are complexes too:
				new_complexes = self.find_complexes(new_substrates)
				if new_complexes != 0:
					substrate_complexes_updated.append(new_complexes)
					# remove the complexes form the list:
					monomers_list.remove(new_complexes)
					monomers_caught.append(monomers_list)

				else:
					monomers_caught.append(monomers_list)

		final_stoich = new_stoich_list.append(monomers_caught)

		return final_stoich







	def load_data(self, simDataFile):
		# code  from processes/metabolism.py  and debug/metabolism.py
		# load in the data:
		allDir = self.ap.get_cells()
		sim_data = self.read_pickle_file(simDataFile)
		simOutDir = os.path.join(allDir[0], "simOut") # only need one directory to access the things needed here
		metabolism = sim_data.process.metabolism

		# Data structures to compute reaction bounds based on enzyme presence/absence
		#self.catalyst_ids = metabolism.catalyst_ids # this is length 1553, which might be more than needed
		#self.reactions_with_catalyst = metabolism.reactions_with_catalyst
		hi = 5
		# technically, I belive there are around 8k catalysts. However, I think I only want the ones involved in the 415 relevant reactions:
		enzymeKineticsReader = TableReader(os.path.join(simOutDir, "EnzymeKinetics"))
		self.constrainedReactions = np.array(enzymeKineticsReader.readAttribute("constrainedReactions")) # length: 415
		#enzymeIDs = enzymeKineticsReader.readAttribute("enzymeIDs") # lenth: 307 (maybe multiple reactions use the same enzyme??)
		enzymeKineticsReader.close()

		# create a complexation reaction to complex dictionary:
		# get the complex ids:
		complex_IDs = sim_data.process.complexation.ids_complexes
		self.complex_IDs = [i[:-3] for i in complex_IDs]
		complexation_reactions = sim_data.process.complexation.ids_reactions
		#complexation_reactions = [i[:-3] for i in complexation_reactions]
		self.complex_IDs_to_reactions = dict(zip(self.complex_IDs, complexation_reactions))
		self.complex_reactions_to_IDs = dict(zip(complexation_reactions, self.complex_IDs))

		# generate a dictionary that maps the complex to the monomers inside it:
		complex_to_monomers = {}
		self.complexation_reactions_df = pd.read_csv(complexation_reactions_tsv_path, sep='\t', comment="#")
		for rxn in complexation_reactions:
			complex_ID = self.complex_reactions_to_IDs.get(rxn)
			# find the row where the complex_ID shows up:
			idx = self.complexation_reactions_df.index[self.complexation_reactions_df['id'] == rxn].tolist()
			hi = 5
			if len(idx) != 1:
				print(f"More than one index for {rxn}: {idx}") # hopefully do not get any of these
				monomers = "None"
			else:
				# find the substrates associated with the reaction:
				stoich_series = self.complexation_reactions_df["stoichiometry"][idx]
				stoich = ast.literal_eval(stoich_series.iloc[0])
				# get a list of the dictionary keys
				monomers_list = list(stoich.keys())
				# remove the complex itself from the keys:
				hi = 5
				if complex_ID in monomers_list:
					monomers_list.remove(complex_ID)
				# check if there are other monomers in the list:
				hi = 4
				substrate_complexes = self.find_complexes(monomers_list)
				hi = 5
				if substrate_complexes != 0:
					monomer_roots = self.find_monomer_roots(monomers_list)
					monomers = monomer_roots

				else:
					monomers = monomers_list
			# append to the dictionary:
			complex_to_monomers[complex_ID] = monomers

		# ok so now I have the massive dictionary, time to find the catalysts:

		metabolic_reaction_to_catalysts = {}
		self.metabolic_reactions_df = pd.read_csv(metabolic_reactions_tsv_path, sep='\t', comment="#")
		for rxn in self.constrainedReactions:
			idx = self.metabolic_reactions_df.index[
				self.metabolic_reactions_df['id'] == rxn].tolist()
			if len(idx) != 1:
				print(f"More than one index for metabloic reaction {rxn}: {idx}") # hopefully do not get any of these
				monomers = "None"
			else:
				# find the substrates associated with the reaction:
				stoich_series = self.metabolic_reactions_df["catalyzed_by"][idx]
				stoich = ast.literal_eval(stoich_series.iloc[0])
				# get a list of the dictionary keys
				monomers_list = stoich.keys()
				# check if there are other complexes in the list:
				substrate_complexes = self.find_complexes(monomers_list)
				if substrate_complexes != 0:
					monomer_roots = self.find_monomer_roots(monomers_list)
					monomers = monomer_roots
				else:
					monomers = monomers_list
			# append to the dictionary:
			metabolic_reaction_to_catalysts[rxn] = monomers

		# finally, search the monomers in the in the input list and see if any show up in metabolic reactions plotted:
		relevant_reactions = {}
		for monomer in IMPORTANT_MONOMERS:
			found_reaction = any(
				monomer in v if isinstance(v, list) else monomer == v
				for v in metabolic_reaction_to_catalysts.values())
			relevant_reactions[monomer] = found_reaction

		# return the reactions with relevancy:
		return(relevant_reactions)




	def do_plot(self, variantDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		# Get all cells
		allDir = self.ap.get_cells()

		sim_data = self.read_pickle_file(simDataFile)

		targetFluxList = []
		actualFluxList = []

		for simDir in allDir:
			simOutDir = os.path.join(simDir, "simOut")

			mainListener = TableReader(os.path.join(simOutDir, "Main"))
			time = mainListener.readColumn("time")
			mainListener.close()
			burnIn = time > BURN_IN_TIME
			burnIn[0] = False

			massListener = TableReader(os.path.join(simOutDir, "Mass"))
			cellMass = massListener.readColumn("cellMass")
			dryMass = massListener.readColumn("dryMass")
			massListener.close()

			coefficient = dryMass / cellMass * sim_data.constants.cell_density.asNumber(MASS_UNITS / VOLUME_UNITS)

			# read constraint data
			enzymeKineticsReader = TableReader(os.path.join(simOutDir, "EnzymeKinetics"))
			allTargetFluxes = (COUNTS_UNITS / MASS_UNITS / TIME_UNITS) * (enzymeKineticsReader.readColumn("targetFluxes").T / coefficient).T
			allActualFluxes = (COUNTS_UNITS / MASS_UNITS / TIME_UNITS) * (enzymeKineticsReader.readColumn("actualFluxes").T / coefficient).T
			kineticsConstrainedReactions = np.array(enzymeKineticsReader.readAttribute("kineticsConstrainedReactions"))
			self.reaction_IDs =  np.array(enzymeKineticsReader.readAttribute("constrainedReactions"))
			constraint_is_kcat_only = np.array(enzymeKineticsReader.readAttribute('constraint_is_kcat_only'))

			allTargetFluxes = allTargetFluxes.asNumber(units.mmol / units.g / units.h)
			allActualFluxes = allActualFluxes.asNumber(units.mmol / units.g / units.h)

			allTargetAve = np.nanmean(allTargetFluxes[burnIn, :], axis = 0)
			allActualAve = np.nanmean(allActualFluxes[burnIn, :], axis = 0)

			if len(targetFluxList) == 0:
				targetFluxList = np.array([allTargetAve])
				actualFluxList = np.array([allActualAve])
			else:
				targetFluxList = np.concatenate((targetFluxList, np.array([allTargetAve])), axis = 0)
				actualFluxList = np.concatenate((actualFluxList, np.array([allActualAve])), axis = 0)

		n_kinetic_constrained_reactions = len(kineticsConstrainedReactions)

		# determine average across all cells
		allTargetAve = np.nanmean(targetFluxList, axis = 0)
		allActualAve = np.nanmean(actualFluxList, axis = 0)

		# see the constrained reactions:
		reaction_IDS = self.reaction_IDs
		# bring in the reactions of interest:
		relevant_rxns = self.load_data(simDataFile)

		# boundary target fluxes
		boundaryTargetAve = allTargetAve[n_kinetic_constrained_reactions:]
		boundaryActualAve = allActualAve[n_kinetic_constrained_reactions:]

		# kinetic target fluxes
		targetAve = allTargetAve[:n_kinetic_constrained_reactions]
		actualAve = allActualAve[:n_kinetic_constrained_reactions]

		# categorize reactions that use constraints with only kcat, Km and kcat, or switch between both types of constraints
		kcatOnlyReactions = constraint_is_kcat_only
		kmAndKcatReactions = ~constraint_is_kcat_only

		# categorize how well the actual flux matches the target flux
		thresholds = [2, 10] # TODO: might want to play around with this
		categorization = np.zeros(n_kinetic_constrained_reactions)
		for i, threshold in enumerate(thresholds):
			categorization[actualAve / targetAve < 1. / threshold] = i + 1
			categorization[actualAve / targetAve > threshold] = i + 1
		categorization[actualAve == 0] = -2
		categorization[actualAve == targetAve] = -1

		# write data for each reaction to a file
		csvFile = io.open(os.path.join(plotOutDir, plotOutFileName + ".tsv"), "wb")
		output = tsv.writer(csvFile)
		output.writerow(["Km and kcat", "Target", "Actual", "Category"])
		for reaction, target, flux, category in zip(kineticsConstrainedReactions[kmAndKcatReactions], targetAve[kmAndKcatReactions], actualAve[kmAndKcatReactions], categorization[kmAndKcatReactions]):
			output.writerow([reaction, target, flux, category])

		output.writerow(["kcat only"])
		for reaction, target, flux, category in zip(kineticsConstrainedReactions[kcatOnlyReactions], targetAve[kcatOnlyReactions], actualAve[kcatOnlyReactions], categorization[kcatOnlyReactions]):
			output.writerow([reaction, target, flux, category])

		csvFile.close() # todo: pressumably this is 415 length (for all the reactions)

		# add small number to allow plotting of 0 flux on log scale
		targetAve += 1e-6
		actualAve += 1e-6

		pearsonAll = pearsonr(np.log10(targetAve), np.log10(actualAve))
		pearsonNoZeros = pearsonr(np.log10(targetAve[(categorization != -2)]), np.log10(actualAve[(categorization != -2)]))

		# plot data
		plt.figure(figsize = (4, 4))
		ax = plt.axes()
		plt.plot([-6, 4], [-6, 4], 'k', linewidth = 0.75)
		plt.plot([-5, 4], [-6, 3], 'k', linewidth = 0.5)
		plt.plot([-6, 3], [-5, 4], 'k', linewidth = 0.5)
		plt.plot(np.log10(targetAve), np.log10(actualAve), 'o', color = "black", markersize = 8, alpha = 0.15, zorder=1, markeredgewidth = 0.0)
		plt.plot(np.log10(boundaryTargetAve), np.log10(boundaryActualAve), "ob", color="red", markeredgewidth=0.25, alpha=0.9, label='boundary fluxes')
		plt.xlabel("Log10(Target Flux [mmol/g/hr])")
		plt.ylabel("Log10(Actual Flux [mmol/g/hr])")
		plt.title("PCC = %.3f, p = %s\n(%.3f, p = %s without points at zero)" % (pearsonAll[0], pearsonAll[1], pearsonNoZeros[0], pearsonNoZeros[1]))
		plt.minorticks_off()
		whitePadSparklineAxis(ax)
		xlim = ax.get_xlim()
		ylim = ax.get_ylim()
		ax.set_ylim(ylim[0] - 0.5, ylim[1])
		ax.set_xlim(xlim[0] - 0.5, xlim[1])
		ax.set_yticks(list(range(-6, int(ylim[1]) + 1, 2)))
		ax.set_xticks(list(range(-6, int(xlim[1]) + 1, 2)))
		ax.legend()

		exportFigure(plt, plotOutDir, plotOutFileName)

		ax.set_xlabel("")
		ax.set_ylabel("")
		ax.set_title("")
		ax.set_xticklabels([])
		ax.set_yticklabels([])

		exportFigure(plt, plotOutDir, plotOutFileName + "_stripped", metadata)
		plt.close("all")


if __name__ == "__main__":
	Plot().cli()
