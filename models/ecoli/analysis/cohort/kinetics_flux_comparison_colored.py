"""
Compare fluxes in simulation to target fluxes
"""

import io
import os
import pickle
import json
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
metabolism_kinetics_tsv_path = "reconstruction/ecoli/flat/metabolism_kinetics.tsv"
WCM_metabolism_gene_xlsx_path = "out/data_tables/WCM_gene_implementation_table_metabolism_06042025.xlsx"

# ignore data from metabolism burnin period
BURN_IN_TIME = 1

IMPORTANT_MONOMERS = ["EG10022-MONOMER", "G6980-MONOMER"]

class Plot(cohortAnalysisPlot.CohortAnalysisPlot):
	_suppress_numpy_warnings = True

	def find_complexes(self, stoich_list):
		if not hasattr(self, 'full_complex_ids_to_rxns'):
			self.find_full_complexation_reaction_dicts()

		complex_IDs_full = list(self.full_complex_ids_to_rxns.keys())

		# check if there happens to be a complex within the list of stoich:
		complex_appearances = []
		for substrate in stoich_list:
			# todo: I do not think self.complex_IDs has all the complexes in it!
			if substrate in complex_IDs_full:
				complex_appearances.append(substrate)

		return complex_appearances


	def find_full_complexation_reaction_dicts(self):
		# since not all complexes are in the complexation reactions data list...
		reaction_IDs = self.complexation_reactions_df['id']

		# find the true complex name in the mix:
		self.full_complex_rxns_to_ids = {}
		self.full_complex_ids_to_rxns = {}
		for rxn in reaction_IDs:
			idx = self.complexation_reactions_df.index[
				self.complexation_reactions_df['id'] == rxn].tolist()
			if len(idx) != 1:
				print(
					f"More than one index for {rxn}: {idx}")  # hopefully do not get any of these
				monomers = "None"
			# next, find the stoich row:
			# find the substrates associated with the reaction:
			#stoich_series = self.complexation_reactions_df["stoichiometry"][idx]
			#stoich = ast.literal_eval(stoich_series.iloc[0])
			hi = 5

			stoich_string = self.complexation_reactions_df["stoichiometry"].iloc[idx[0]]
			hi = 5
			stoich = json.loads(stoich_string) # trying this to handle the null values

			# find the key that has a negative value
			complex_ID = ''
			for key, num in stoich.items():
				if num > 0:
					complex_ID = key
					break


			# append the dictionary entries:
			hi = 5
			self.full_complex_rxns_to_ids[rxn] = complex_ID
			self.full_complex_ids_to_rxns[complex_ID] = rxn


	def find_complex_ID_from_reaction(self, reaction):
		# create the dictionaries if they do not already exist
		if not hasattr(self, 'full_complex_rxns_to_ids'):
			self.find_full_complexation_reaction_dicts()

		# find the complex ID:
		return self.full_complex_rxns_to_ids.get(reaction)

	def find_reaction_ID_from_complex(self, complex):
		if not hasattr(self, 'full_complex_ids_to_rxns'):
			self.find_full_complexation_reaction_dicts()

		# find the complex ID:
		return self.full_complex_ids_to_rxns.get(complex)

	def find_monomer_roots(self, stoich_list):
		# find the monomer roots of each complex:
		original_stoich_list = stoich_list

		# determine if any substrates are complexes too:
		substrate_complexes = self.find_complexes(original_stoich_list)
		substrate_complexes = list(substrate_complexes)
		hi = 5
		# remove the complexes that have been found:
		new_stoich_list = original_stoich_list.copy()

		for complex in substrate_complexes:
			new_stoich_list.remove(complex)

		monomers_caught = []
		substrate_complexes_caught = []
		substrate_complexes_updated = substrate_complexes
		while len(substrate_complexes_updated) != 0:
			new_substrates = []
			for complex in substrate_complexes_updated:
				print(len(substrate_complexes_updated), substrate_complexes_updated)
				print("current complex:", complex)
				#rxn = self.complex_IDs_to_reactions.get(complex) # this does not work actually ughhhhh
				rxn = self.find_reaction_ID_from_complex(complex)


				# get the
				idx = self.complexation_reactions_df.index[
					self.complexation_reactions_df['id'] == rxn].tolist()
				if len(idx) != 1:
					print(
						f"More than one index for {rxn}: {idx}")  # hopefully do not get any of these
					monomers = "None"

				# find the substrates associated with the reaction:
				stoich_string = self.complexation_reactions_df["stoichiometry"].iloc[idx[0]]
				hi = 5
				stoich = json.loads(stoich_string)

				# get a list of the dictionary keys
				monomers_list = list(stoich.keys())
				# remove the complex itself from the keys:
				monomers_list.remove(complex)

				# put the new substrate list in the mix:
				new_substrates.append(monomers_list)

				# remove the current complex from the list of substrates:
				substrate_complexes_updated.remove(complex)
				substrate_complexes_caught.append(complex)

				# check if the new monomers are complexes too:
				new_complexes = self.find_complexes(monomers_list)
				if len(new_complexes) != 0:
					if complex == 'ATPD-CPLX':
						hi = 5
					# append the complexes within this complex to the list if they exist! they technically should be handled!

					# remove the complexes form the list:
					print("complex with a complex inside:", complex)
					print(monomers_list)
					hi = 5
					# a couple complexes within complexes, so those need to be taken out one by one
					for comp in new_complexes:
						monomers_list.remove(comp)
						substrate_complexes_updated.append(comp)
					hi = 5
					# append left over monomers to monomer list (I do it this way so the item handleing at the bottom is consistent)
					for monomer in monomers_list:
						monomers_caught.append([monomer])
					hi = 6
				else:
					monomers_caught.append(monomers_list)

		hi= 5
		final_stoich = []
		for monomer in new_stoich_list:
			hi = 5
			final_stoich.append([monomer])
		for monomer in monomers_caught:
			final_stoich.append(monomer)
		hi = 5
		final_stoich_flattened = [item[0] for item in final_stoich]
		return final_stoich_flattened

	hi = 1

	def load_data(self, simDataFile):
		# code  from processes/metabolism.py  and debug/metabolism.py
		# load in the data:
		allDir = self.ap.get_cells()
		sim_data = self.read_pickle_file(simDataFile)
		simOutDir = os.path.join(allDir[0], "simOut") # only need one directory to access the things needed here

		# todo: check what these are:
		# Data structures to compute reaction bounds based on enzyme presence/absence
		#self.catalyst_ids = metabolism.catalyst_ids # this is length 1553, which might be more than needed
		#self.reactions_with_catalyst = metabolism.reactions_with_catalyst
		hi = 5

		# load in the complex IDs and reaction IDs
		complex_IDs = sim_data.process.complexation.ids_complexes # todo: note: I do not think this has ALL complexes in it, like E20 and E10 are not in it?
		self.complex_IDs = [i[:-3] for i in complex_IDs] # todo: do i even use these?
		complexation_reactions = sim_data.process.complexation.ids_reactions

		# create a dictionary mapping the complexation reaction IDs to their complex IDs (and vice versa)
		# note: do so using this function becuase the above two variables are not a 1:1 match despite being the same length
		self.complexation_reactions_df = pd.read_csv(complexation_reactions_tsv_path, sep='\t',
													 comment="#")
		self.find_full_complexation_reaction_dicts()

		# technically, I belive there are around 8k catalysts. However, I think I only want the ones involved in the 415 relevant reactions:
		enzymeKineticsReader = TableReader(os.path.join(simOutDir, "EnzymeKinetics"))
		self.constrainedReactions = np.array(enzymeKineticsReader.readAttribute("constrainedReactions")) # length: 415
		#enzymeIDs = enzymeKineticsReader.readAttribute("enzymeIDs") # lenth: 307 (maybe multiple reactions use the same enzyme??)
		enzymeKineticsReader.close()

		# create a complexation reaction to complex dictionary:
		# get the complex ids:
		# complex_IDs = sim_data.process.complexation.ids_complexes
		# self.complex_IDs = [i[:-3] for i in complex_IDs] # TODO: this is not in the same order as complexes_IDs!!!!! CANNOT USE THIS!
		# complexation_reactions = sim_data.process.complexation.ids_reactions
		# #complexation_reactions = [i[:-3] for i in complexation_reactions]
		# self.complex_IDs_to_reactions = dict(zip(self.complex_IDs, complexation_reactions)) # todo: cannot use this dictionary
		# self.complex_reactions_to_IDs = dict(zip(complexation_reactions, self.complex_IDs))
		hi = 5
		# generate a dictionary that maps the complex to the monomers inside it:
		complex_to_monomers_dict = {}
		# todo: why is complexation reactions only 1109 long? should I be doing them all? or are these the only relevant ones in the code that are used in the sims?
		for rxn in complexation_reactions:
			complex_ID = self.full_complex_rxns_to_ids.get(rxn) #
			# find the row where the complex_ID shows up:
			idx = self.complexation_reactions_df.index[self.complexation_reactions_df['id'] == rxn].tolist()
			hi = 5
			if len(idx) != 1:
				print(f"More than one index or no indexes for {rxn}: {idx}") # hopefully do not get any of these
				monomers = "None"
			else:
				# find the substrates associated with the reaction:
				#stoich_series = self.complexation_reactions_df["stoichiometry"][idx]
				#stoich = ast.literal_eval(stoich_series.iloc[0])
				# todo: trying to do the json version instead incase the null is causing fail here too:
				stoich_string = self.complexation_reactions_df["stoichiometry"].iloc[idx[0]]
				hi = 5
				stoich = json.loads(stoich_string)  # trying this to handle the null values
				# get a list of the dictionary keys
				monomers_list = list(stoich.keys())
				# remove the complex itself from the keys:
				hi = 5
				print("monomers_list BEFORE:", monomers_list)
				if complex_ID in monomers_list:
					monomers_list.remove(complex_ID)
				# check if there are other monomers in the list:
				hi = 4

				substrate_complexes = self.find_complexes(monomers_list)
				hi = 5
				# there are some complexes that are not in the self.complex_reactions_to_IDs dictionary, so this goes to another funciton that can break down everything
				if len(substrate_complexes) != 0:
					hi = 5
					monomer_roots = self.find_monomer_roots(monomers_list)
					monomers = monomer_roots
					print("monomers_list AFTER:", monomers)

				else:
					monomers = monomers_list
					print("monomers_list AFTER:", monomers)

			# append to the dictionary:
			hi = 5
			complex_to_monomers_dict[complex_ID] = monomers

		# ok so now I have the massive dictionary, time to find the catalysts:
		# todo: consider just using the above funciton for each of the complexes? but I think i did this so the output would come out consistantly formatted ie all were ["protein21", "protein2","protein1"] without extra brakets inside
		hi = 5
		metabolic_reaction_to_catalysts = {}
		self.metabolic_reactions_df = pd.read_csv(metabolic_reactions_tsv_path, sep='\t', comment="#")
		reaction_IDs = self.metabolic_reactions_df["id"]
		# todo: decide if I should only look through self.constrainedReactions here
		for rxn in reaction_IDs:
			idx = self.metabolic_reactions_df.index[
				self.metabolic_reactions_df['id'] == rxn].tolist()
			if len(idx) != 1:
				if len(idx) == 0:
					print(f"no catalysts for {rxn}")
				if len(idx) > 1:
					print(f"More than one index for metabloic reaction {rxn}: {idx}") # hopefully do not get any of these
				monomers = "None"
			else:
				# find the enzyme associated with the reaction:
				enzyme_string = self.metabolic_reactions_df["catalyzed_by"].iloc[idx[0]]
				hi = 5
				enzyme_list = json.loads(enzyme_string)
				# check if there are other complexes in the list:
				substrate_complexes = self.find_complexes(enzyme_list)
				if len(substrate_complexes) != 0:
					hi = 5
					monomer_roots = self.find_monomer_roots(enzyme_list)
					monomers = monomer_roots
				else:
					monomers = enzyme_list
			# append to the dictionary:
			metabolic_reaction_to_catalysts[rxn] = monomers

		#todo: create a gene ID to protein ID dictionary, so that I can map the WCM gene implementation dashboard data into this

		# finally, search the monomers in the in the input list and see if any show up in metabolic reactions plotted:
		hi = 5
		relevant_reactions = {}
		for monomer in IMPORTANT_MONOMERS:
			matches = []
			for reaction, catalysts in metabolic_reaction_to_catalysts.items():
				hi = 5
				if isinstance(catalysts, list) and monomer in catalysts:

					matches.append(reaction)
			relevant_reactions[monomer] = matches

		# now check the metabolism_kinetics.tsv file for matches: (note, there are repeated reactions in this dataframe! just use the first one)
		self.metabolism_kinetics_df = pd.read_csv(metabolism_kinetics_tsv_path, sep='\t',
												  comment="#")

		kinetic_reaction_IDs = self.metabolism_kinetics_df["reactionID"]
		kinetic_reaction_enzymes = self.metabolism_kinetics_df["enzymeID"]
		kinetic_reaction_substrates = self.metabolism_kinetics_df["substrateIDs"]

		# map these!
		kinetic_reaction_to_catalysts = {}
		for rxn in kinetic_reaction_IDs:
			idx = self.metabolism_kinetics_df.index[
				self.metabolism_kinetics_df['reactionID'] == rxn].tolist()
			if len(idx) != 1:
				if len(idx) == 0:
					print(f"no catalysts for {rxn}")
					monomers = "None"
				if len(idx) > 1:
					print(f"More than one index for kinetic reaction {rxn}: {idx}")
					if rxn in kinetic_reaction_to_catalysts.keys():
						print("Skipping repeat appearance of:",rxn)
						continue
					else:
						idx = [idx[0]]
						print(f"Using first appearance of reaction {rxn} at index {idx}")
						# find the enzymes associated with the reaction (in this file, there seems to only be one per):
						enzyme_string = self.metabolism_kinetics_df["enzymeID"].iloc[idx[0]]
						hi = 5
						print(rxn, enzyme_string)
						enzyme_list = [enzyme_string]
						# check if there are other complexes in the list:
						enzyme_complexes = self.find_complexes(enzyme_list)
						if len(enzyme_complexes) != 0:
							hi = 5
							monomer_roots = self.find_monomer_roots(enzyme_list)
							monomers = monomer_roots
						else:
							monomers = enzyme_list

			else:
				# find the enzymes associated with the reaction (in this file, there seems to only be one per):
				enzyme_string = self.metabolism_kinetics_df["enzymeID"].iloc[idx[0]]
				hi = 5
				print(rxn, enzyme_string)
				enzyme_list = [enzyme_string]
				# check if there are other complexes in the list:
				enzyme_complexes = self.find_complexes(enzyme_list)
				if len(enzyme_complexes) != 0:
					hi = 5
					monomer_roots = self.find_monomer_roots(enzyme_list)
					monomers = monomer_roots
				else:
					monomers = enzyme_list
			# append to the dictionary:
			kinetic_reaction_to_catalysts[rxn] = monomers

		# find relevant kinetics reaction matches:
		relevant_kinetics_reactions = {}
		for monomer in IMPORTANT_MONOMERS:
			matches = []
			for reaction, catalysts in metabolic_reaction_to_catalysts.items():
				hi = 5
				if isinstance(catalysts, list) and monomer in catalysts:
					matches.append(reaction)
			relevant_kinetics_reactions[monomer] = matches

		# check if any matches between the substrates in the kinetics file appear:
		hi = 5
		all_monomer_IDs = [i[:-3] for i in sim_data.process.translation.monomer_data["id"]]

		kinetic_reactions_with_monomer_substrates = {}
		for rxn in kinetic_reaction_IDs:
			idx = self.metabolism_kinetics_df.index[
				self.metabolism_kinetics_df['reactionID'] == rxn].tolist()
			if len(idx) != 1:
				if len(idx) == 0:
					print(f"no catalysts for {rxn}")
					continue
				if len(idx) > 1:
					print(f"More than one index for kinetic reaction {rxn}: {idx}")
					if rxn in kinetic_reaction_to_catalysts.keys():
						print("Skipping repeat appearance of:", rxn)
						continue
					else:
						idx = [idx[0]]
						print(f"Using first appearance of reaction {rxn} at index {idx}")
						# check if the substrate is a monomer or complex comprised of monomers:
						substrate_string = self.metabolism_kinetics_df["substrateIDs"].iloc[idx[0]]
						substrate_list = json.loads(substrate_string)
						monomers_caught = []
						for substrate in substrate_list:
							if substrate in all_monomer_IDs:
								monomer_ID = all_monomer_IDs[
									np.where(all_monomer_IDs == substrate)]
								monomers_caught.append(monomer_ID)
							# check if the substrate is a complex:
							substrate_complexes = self.find_complexes(substrate)
							if len(substrate_complexes) != 0:
								monomer_roots = self.find_monomer_roots(substrate_complexes)
								for monomer in monomer_roots:
									if monomer in all_monomer_IDs:
										monomers_caught.append(monomer)
						if len(monomers_caught) != 0:
							kinetic_reactions_with_monomer_substrates[rxn] = monomers_caught

			else:
				# check if the substrate is a monomer or complex comprised of monomers:
				substrate_string = self.metabolism_kinetics_df["substrateIDs"].iloc[idx[0]]
				substrate_list = json.loads(substrate_string)
				monomers_caught = []
				for substrate in substrate_list:
					if substrate in all_monomer_IDs:
						print(substrate)
						hi = 89
						#monomer_idx = all_monomer_IDs.index(substrate)
						#monomer_ID = all_monomer_IDs[monomer_idx] # doing it this way in case only a slice matches
						monomers_caught.append(substrate)
					# check if the substrate is a complex:
					substrate_complexes = self.find_complexes(substrate)
					if len(substrate_complexes) != 0:
						monomer_roots = self.find_monomer_roots(substrate_complexes)
						for monomer in monomer_roots:
							if monomer in all_monomer_IDs:
								monomers_caught.append(monomer)
				if len(monomers_caught) != 0:
					kinetic_reactions_with_monomer_substrates[rxn] = monomers_caught


		hi = 54


		# check if any of the metabolism genes implemented are here:
		self.WCM_gene_implementation_df = pd.read_excel(WCM_metabolism_gene_xlsx_path)

		hi = 65


		# return the reactions with relevancy:
		return relevant_reactions




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

		# boundary target fluxes
		boundaryTargetAve = allTargetAve[n_kinetic_constrained_reactions:]
		boundaryActualAve = allActualAve[n_kinetic_constrained_reactions:]

		# kinetic target fluxes
		targetAve = allTargetAve[:n_kinetic_constrained_reactions]
		actualAve = allActualAve[:n_kinetic_constrained_reactions]

		hello, yellow = self.load_data(simDataFile)
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
