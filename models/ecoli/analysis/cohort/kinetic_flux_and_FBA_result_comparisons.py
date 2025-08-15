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
from wholecell.analysis.analysis_tools import (exportFigure,
	read_stacked_bulk_molecules, read_stacked_columns) # todo: delete later if needed


# ignore data from metabolism burnin period
BURN_IN_TIME = 1
IGNORE_FIRST_N_GENERATIONS = 0
IMPORTANT_MONOMERS = ["ISOCIT-LYASE-MONOMER",  'EG11111-MONOMER', "G6980-MONOMER", "BASS-MONOMER", 'G6986-MONOMER[c]']

# todo: make it so that the important monomners are marked and plotted automatically if they exist in either of the results plotted
# have the relevant tsv file paths:
os.chdir(os.path.expanduser('~/wcEcoli/'))
metabolic_reactions_tsv_path = "reconstruction/ecoli/flat/metabolic_reactions.tsv"
complexation_reactions_tsv_path = "reconstruction/ecoli/flat/complexation_reactions.tsv"
metabolism_kinetics_tsv_path = "reconstruction/ecoli/flat/metabolism_kinetics.tsv"
WCM_metabolism_gene_xlsx_path = "out/data_tables/WCM_gene_implementation_table_metabolism_06042025.xlsx"
ecocyc_proteins_to_genes_xlsx_path = "out/data_tables/ecocyc_proteins_to_genes_05302025.xlsx"

class Plot(cohortAnalysisPlot.CohortAnalysisPlot):
	_suppress_numpy_warnings = True

	def get_monomer_name_with_compartment(self, simDataFile, monomer):
		"""
		Obtain the compartment of the monomer from the simDataFile.
		"""
		sim_data = self.read_pickle_file(simDataFile)
		if not monomer.endswith("]"):
			if sim_data.getter.is_valid_molecule(monomer):
				edited_monomer_name = monomer + "[" + sim_data.getter.get_compartment(monomer)[
					0] + "]"
			else:
				print(f"Monomer {monomer} not found in simDataFile.")
				edited_monomer_name = None
		else:
			# pressumably it already has the compartment tag, so just return it:
			print(f"Monomer {monomer} already has a compartment tag.")
			edited_monomer_name = monomer

		return edited_monomer_name


	def get_relevant_monomers(self, simDataFile):
		"""
		Obtain the catalyst monomers from the simDataFile.
		"""
		sim_data = self.read_pickle_file(simDataFile)

		# attach the compartment tag to the monomers if it does not have one:
		important_monomers = []
		for monomer in IMPORTANT_MONOMERS:
			edited_monomer_name = self.get_monomer_name_with_compartment(simDataFile, monomer)
			if edited_monomer_name is None:
				print(f"Monomer {monomer} not found in simDataFile, skipping.")
				continue
			else:
				important_monomers.append(edited_monomer_name)


		reaction_catalysts_dict = sim_data.process.metabolism.reaction_catalysts

		# should I make a specific reaction to catalyst monomer dict? or a catalyst to monomer dict, and then map that to the reaction?
		# well I want to be able to save the specific reaction catalyst as complex (monomer name) if the monomer only is part of that reaction based on its monomer identity
		reaction_to_catalyst_monomers_dict = {}
		for reaction in reaction_catalysts_dict.keys():
			catalyst_monomers = []
			for catalyst in reaction_catalysts_dict[reaction]:
				# append the compartment tag to the catalyst if it does not have one:
				if not catalyst.endswith("]"):
					catalyst = self.get_monomer_name_with_compartment(simDataFile, catalyst)
				# check if the catalyst is an equilibrium complex:
				if catalyst in sim_data.process.equilibrium.ids_complexes:
					print(f"Found equilibrium complex in {reaction}: {catalyst}")
					equilibrium_monomers = sim_data.process.equilibrium.get_monomers(catalyst)['subunitIds']
					for monomer in equilibrium_monomers:
						catalyst_monomers.append(monomer)
				# check if the catalyst is a complexation complex:
				complexation_monomers = sim_data.process.complexation.get_monomers(catalyst)['subunitIds'] # technically this funciton works without using the if statement used for equilibrium complexes, but without it I can auto add the monomer only names to the list. or should I have an else statement for that to speed up the code and make this one an elif?
				for monomer in complexation_monomers:
					catalyst_monomers.append(monomer)
			# remove duplicates:
			catalyst_monomers = list(set(catalyst_monomers))
			# append to the dictionary:
			reaction_to_catalyst_monomers_dict[reaction] = catalyst_monomers

		# for each monomer listed at the start of the code, check if it exists in the reaction to catalyst monomers dictionary:
		relevant_reactions_as_catalysts = {}
		name = ''
		for monomer in important_monomers:
			reaction_to_name_info = {}
			for reaction, catalysts in reaction_to_catalyst_monomers_dict.items():
				if isinstance(catalysts, list) and monomer in catalysts:
					# get the catalyst(s) for the reaction:
					reaction_catalysts = reaction_catalysts_dict[reaction]
					# if the monomer is directly a catalyst, then add it directly to the reaction dictionary:
					if monomer in reaction_catalysts:
						name = monomer
						print(f"Found monomer {monomer} in reaction {reaction} as a catalyst.")
					else:
						# if the monomer is included as a catalyst due it it being in a complex that contains the catalytic monomer (or the catalyst is the complex itself), then add it as a complex:
						for catalyst in reaction_catalysts:
							# attach the compartment tag to the catalyst if it does not have one:
							if not catalyst.endswith("]"):
								catalyst = self.get_monomer_name_with_compartment(simDataFile,
																			  catalyst)
							# check if it is in an equilibrium complex:
							if catalyst in sim_data.process.equilibrium.ids_complexes:
								equilibrium_monomers = sim_data.process.equilibrium.get_monomers(
									catalyst)['subunitIds']
								if monomer in equilibrium_monomers:
									name = monomer + f" (via {catalyst})"
							# check if it is in a complexation complex:
							complexation_monomers = sim_data.process.complexation.get_monomers(
								catalyst)['subunitIds']
							if monomer in complexation_monomers:
								name = monomer + f" (via {catalyst})"

					# apppend the reaction and the name of it to the dictionary:
					if name != '':
						reaction_to_name_info[reaction] = name

				# append the reactions to the dictionary:
				if reaction_to_name_info != {}:
					relevant_reactions_as_catalysts[monomer] = reaction_to_name_info


		# next, search the reaction_stoich dictionary to see if any of the reactions have monomers in them:
		reaction_stoich_dict = sim_data.process.metabolism.reaction_stoich
		reaction_to_substrate_monomers_dict = {}
		for reaction in reaction_stoich_dict.keys():
			substrate_monomers = []
			for substrate in reaction_stoich_dict[reaction]:
				# append the compartment tag to the substrate if it does not have one:
				if not substrate.endswith("]"):
					substrate = self.get_monomer_name_with_compartment(simDataFile, substrate)
				# check if the substrate is a monomer:
				if substrate in sim_data.process.translation.monomer_data['id']:
					substrate_monomers.append(substrate)
				# check if the substrate is an equilibrium complex:
				elif substrate in sim_data.process.equilibrium.ids_complexes:
					print(f"Found equilibrium complex in {reaction}: {substrate}")
					equilibrium_monomers = sim_data.process.equilibrium.get_monomers(substrate)['subunitIds']
					for monomer in equilibrium_monomers:
						substrate_monomers.append(monomer)
				# check if the catalyst is a complexation complex:
				elif substrate in sim_data.process.complexation.ids_complexes:
					complexation_monomers = sim_data.process.complexation.get_monomers(substrate)['subunitIds'] # technically this funciton works without using the if statement used for equilibrium complexes, but without it I can auto add the monomer only names to the list. or should I have an else statement for that to speed up the code and make this one an elif?
					for monomer in complexation_monomers:
						substrate_monomers.append(monomer)
				else:
					# if the substrate is not a monomer or a complex, then it is not relevant:
					continue

			# remove duplicates:
			substrate_monomers = list(set(substrate_monomers))

			# only append the reaction to the dictionary if there are any monomers in it:
			if substrate_monomers == []:
				continue
			else:
				reaction_to_substrate_monomers_dict[reaction] = substrate_monomers

		# check if there are any reactions that have the monomers in them as substrates:
		relevant_reactions_as_substrates = {}
		name = ''
		for monomer in important_monomers:
			reaction_to_name_info = {}
			for reaction, substrates in reaction_to_substrate_monomers_dict.items():
				if isinstance(substrates, list) and monomer in substrates:
					# get the substrate(s) for the reaction:
					reaction_substrates = reaction_stoich_dict[reaction]
					# if the monomer is directly a substrate, then add it directly to the reaction dictionary:
					if monomer in reaction_substrates:
						name = monomer
					else:
						# if the monomer is included as a substrate due it it being in a complex that contains the substrate monomer (or the substrate is the complex itself), then add it as a complex:
						for substrate in reaction_substrates:
							# attach the compartment tag to the substrate if it does not have one:
							if not substrate.endswith("]"):
								substrate = self.get_monomer_name_with_compartment(simDataFile,
																				  substrate)
							# check if it is in an equilibrium complex:
							if substrate in sim_data.process.equilibrium.ids_complexes:
								equilibrium_monomers = sim_data.process.equilibrium.get_monomers(
									substrate)['subunitIds']
								if monomer in equilibrium_monomers:
									name = monomer + f" (via {substrate})"
							# check if it is in a complexation complex:
							complexation_monomers = sim_data.process.complexation.get_monomers(
								substrate)['subunitIds']
							if monomer in complexation_monomers:
								name = monomer + f" (via {substrate})"

					# apppend the reaction and the name of it to the dictionary:
					if name != '':
						reaction_to_name_info[reaction] = name

				# append the reactions to the dictionary:
				if reaction_to_name_info != {}:
					relevant_reactions_as_substrates[monomer] = reaction_to_name_info

		return relevant_reactions_as_catalysts, relevant_reactions_as_substrates


	# version of funciton 1 WITH THE NAMES FROM THE LISTENER
	def FBA_plots_averaging_method_1a(self, simDataFile, validationDataFile, plotOutFileName, plotOutDir, metadata):
		"""
		average over all time steps in the simulations. So some cells are weighted more than others if they have a longer simulation time.
		"""
		sim_data = self.read_pickle_file(simDataFile)
		simOutDir = os.path.join(self.all_cells[0], "simOut")
		mmol_per_g_per_h = units.mmol / units.g / units.h


		# get the reaction IDs:
		fbaResults = TableReader(os.path.join(simOutDir, "FBAResults"))
		base_reaction_ids = fbaResults.readAttribute("base_reaction_ids") # can also access in sim data, but not sure if the mapping is consistent but the lengths are the same, so maybe? # todo: ask riley

		# get the FBA results from the simDataFile:
		self.rxn_fluxes = read_stacked_columns(self.all_cells, "FBAResults", "reactionFluxes")
		self.base_rxn_fluxes = read_stacked_columns(self.all_cells, "FBAResults", "base_reaction_fluxes")

		# obtain the simulation fluxes:
		cellMass = read_stacked_columns(self.all_cells, "Mass", "cellMass")
		dryMass = read_stacked_columns(self.all_cells, "Mass", "dryMass")
		cellDensity = sim_data.constants.cell_density
		coefficient = dryMass / cellMass * cellDensity.asNumber(MASS_UNITS / VOLUME_UNITS)

		# get the reaction fluxes:
		base_reaction_fluxes = (COUNTS_UNITS / MASS_UNITS / TIME_UNITS) * (
			(self.base_rxn_fluxes).T / coefficient.squeeze()).T
		base_reaction_id_to_index = {
			rxn_id: i for (i, rxn_id) in enumerate(base_reaction_ids)
		}

		# get the Toya fluxes:
		validation_data = self.read_pickle_file(validationDataFile)
		toyaReactions = validation_data.reactionFlux.toya2010fluxes["reactionID"]
		toyaFluxes = validation_data.reactionFlux.toya2010fluxes["reactionFlux"]
		toyaStdev = validation_data.reactionFlux.toya2010fluxes["reactionFluxStdev"]
		toyaFluxesDict = dict(zip(toyaReactions, toyaFluxes))
		toyaStdevDict = dict(zip(toyaReactions, toyaStdev))


		# average the fluxes over all time steps if there is a matching reaction in the model:
		fluxTimeCourse_map = {}
		modelFluxes = {}
		for toyaReaction in toyaReactions:
			if toyaReaction in base_reaction_id_to_index:
				rxn_index = base_reaction_id_to_index[toyaReaction]
				fluxTimeCourse = base_reaction_fluxes[:, rxn_index]
				fluxTimeCourse_map[toyaReaction] = fluxTimeCourse.asNumber(mmol_per_g_per_h)
				modelFluxes[toyaReaction] = np.mean(fluxTimeCourse).asNumber(mmol_per_g_per_h)

		# make a dictionary with the reaction fluxes and their standard deviations (note that the extra np.mean() of the flux does not change anything here as this averaging method already averaged over all the cells present by averaging over the full time)
		toyaVsReactionAve = [] # todo: is this the best way to do the STD for the model?
		for rxn, toyaFlux in toyaFluxesDict.items():
			if rxn in modelFluxes:
				toyaVsReactionAve.append(
					(np.mean(modelFluxes[rxn]), # mean WCM flux
					 toyaFlux.asNumber(mmol_per_g_per_h), # Toya flux
					 np.std(fluxTimeCourse_map[rxn]), # WCM flux STD over the entire simulation duration
					 toyaStdevDict[rxn].asNumber(mmol_per_g_per_h))) # Toya STD

		toyaVsReactionAve = np.array(toyaVsReactionAve)

		# Calulate r for all data points as well as just a select few
		rWithAll = pearsonr(toyaVsReactionAve[:, 0], toyaVsReactionAve[:, 1])
		# removes outliers if there are any:
		idx = np.abs(toyaVsReactionAve[:, 0]) < 5 * np.abs(toyaVsReactionAve[:, 1])
		rWithoutOutliers = pearsonr(toyaVsReactionAve[idx, 0], toyaVsReactionAve[idx, 1])

		# plot the results
		plt.figure(figsize=(7, 7), dpi=100)
		ax = plt.axes()
		plt.title(
			f"Central Carbon Metabolism Flux Comparison ", fontsize=12, pad=20)
		plt.suptitle(f"(averaged over {len(self.all_cells)} cells, each weighted by cell cycle length)", fontsize=8, y=.90)
		plt.errorbar(toyaVsReactionAve[:, 1], toyaVsReactionAve[:, 0],
					 xerr=toyaVsReactionAve[:, 3], yerr=toyaVsReactionAve[:, 2], fmt=".",
					 ecolor="k", alpha=0.5, linewidth=0.5)
		ylim = plt.ylim()
		plt.plot([ylim[0], ylim[1]], [ylim[0], ylim[1]], color="k")
		plt.plot(toyaVsReactionAve[:, 1], toyaVsReactionAve[:, 0], "ob", markeredgewidth=0.1,
				 alpha=0.9)
		plt.xlabel("Toya 2010 Reaction Flux [mmol/g/hr]")
		plt.ylabel(f"Mean WCM Reaction Flux [mmol/g/hr] \n(Simulation ID: {self.sim_id})")

		# noinspection PyTypeChecker
		ax.set_xlim([-20, 30])
		xlim = ax.get_xlim()
		ylim = ax.get_ylim()
		#ax.set_yticks(list(range(int(ylim[0]), int(ylim[1]) + 1, 10)))
		#ax.set_xticks(list(range(int(xlim[0]), int(xlim[1]) + 1, 10)))
		ax.text(.95, 0.075, "Pearson R = %.4f, p = %s\n(%.4f, %s without outliers)" % (
				rWithAll[0], rWithAll[1], rWithoutOutliers[0], rWithoutOutliers[1]), fontsize=6, transform=ax.transAxes,
				ha='right', color="grey")

		plotOutName = plotOutFileName + "_FBA_avgeraging_method_1a_" + self.sim_id + ".png"
		exportFigure(plt, plotOutDir, plotOutName, metadata)
		plt.close("all")


		hi = 5

	# function 2
	def FBA_plots_averaging_method_1b(self, simDataFile, validationDataFile, plotOutFileName, plotOutDir, metadata):
		"""
		average over all time steps in the simulations,
		"""
		sim_data = self.read_pickle_file(simDataFile)
		simOutDir = os.path.join(self.all_cells[0], "simOut")
		mmol_per_g_per_h = units.mmol / units.g / units.h


		# get the reaction IDs:
		#base_reaction_IDs = fbaResults.readAttribute("base_reaction_ids") # can also access in sim data, but not sure if the mapping is consistent but the lengths are the same, so maybe? # todo: ask riley
		base_reaction_ids = sim_data.process.metabolism.base_reaction_ids

		# other ways to get the reaction ids: (already said the base one)
		#all_reaciton_ids = list(sim_data.process.metabolism.reaction_id_to_base_reaction_id.keys())
		#reaction_stoich_ids = list(sim_data.process.metabolism.reaction_stoich.keys())

		# see what reaction IDs are in the all reaction ids of the reaction_stoich_ids:
		#reaction_stoich_ids_set = set(reaction_stoich_ids).intersection(set(all_reaciton_ids))
		# see which do not match:
		#reaction_stoich_ids_not_in_all = set(all_reaciton_ids).difference(set(reaction_stoich_ids)) # I believe the 21 that do not match have to do with aa

		# get the FBA results from the simDataFile:
		self.rxn_fluxes = read_stacked_columns(self.all_cells, "FBAResults", "reactionFluxes")
		self.base_rxn_fluxes = read_stacked_columns(self.all_cells, "FBAResults",
													"base_reaction_fluxes")

		# obtain the simulation fluxes:
		cellMass = read_stacked_columns(self.all_cells, "Mass", "cellMass")
		dryMass = read_stacked_columns(self.all_cells, "Mass", "dryMass")
		cellDensity = sim_data.constants.cell_density
		coefficient = dryMass / cellMass * cellDensity.asNumber(MASS_UNITS / VOLUME_UNITS)

		# get the reaction fluxes:
		base_reaction_fluxes = (COUNTS_UNITS / MASS_UNITS / TIME_UNITS) * (
				(self.base_rxn_fluxes).T / coefficient.squeeze()).T
		base_reaction_id_to_index = {
			rxn_id: i for (i, rxn_id) in enumerate(base_reaction_ids)
		} # TODO: either this base_reaction_ids or the one in 1a. choose! even tho they are the same dont use both

		# get the Toya fluxes:
		validation_data = self.read_pickle_file(validationDataFile)
		toyaReactions = validation_data.reactionFlux.toya2010fluxes["reactionID"]
		toyaFluxes = validation_data.reactionFlux.toya2010fluxes["reactionFlux"]
		toyaStdev = validation_data.reactionFlux.toya2010fluxes["reactionFluxStdev"]
		toyaFluxesDict = dict(zip(toyaReactions, toyaFluxes))
		toyaStdevDict = dict(zip(toyaReactions, toyaStdev))

		# match the toya fluxes to the base reaction fluxes:
		modelFluxes = {}
		toyaOrder = [] # todo: can I take this out I do not use it from what I can tell?
		for rxn in toyaReactions: # todo: should this have .keys()?
			modelFluxes[rxn] = []
			toyaOrder.append(rxn)
		hey = 1
		# average the fluxes over all time steps if there is a matching reaction in the model:
		fluxTimeCourse_map = {}
		for toyaReaction in toyaReactions:
			if toyaReaction in base_reaction_id_to_index: # todo: does this have to have .keys()?
				rxn_index = base_reaction_id_to_index[toyaReaction]
				fluxTimeCourse = base_reaction_fluxes[:, rxn_index]
				fluxTimeCourse_map[toyaReaction] = fluxTimeCourse.asNumber(mmol_per_g_per_h)
				modelFluxes[toyaReaction].append(
					np.mean(fluxTimeCourse).asNumber(mmol_per_g_per_h))
				# todo: check that this is indeed appending to the dictionary as I would like it to be, bc maybe its best to justmake the modelFluxes dict at this loop and not the one before it
		hey = 3
		# make a dictionary with the reaction fluxes and their standard deviations (note that the extra np.mean() of the flux does not change anything)
		# TODO: ask riley if I should just be doing average by the cell, as this would actually yeild a standard deviation for each cell (written the way it is in the original file). I can edit it though to do over all time points, but idk if this is accurate
		toyaVsReactionAve = [] # will contain the mean WCM flux, mean Toya flux, WCM STD
		for rxn, toyaFlux in toyaFluxesDict.items():
			if rxn in modelFluxes:
				toyaVsReactionAve.append(
					(np.mean(modelFluxes[rxn]),
					 toyaFlux.asNumber(mmol_per_g_per_h),
					 np.std(fluxTimeCourse_map[rxn]),
					 toyaStdevDict[rxn].asNumber(mmol_per_g_per_h)))

		toyaVsReactionAve = np.array(toyaVsReactionAve)


		# build the plot
		idx = np.abs(toyaVsReactionAve[:, 0]) < 5 * np.abs(toyaVsReactionAve[:, 1])
		rWithAll = pearsonr(toyaVsReactionAve[:, 0], toyaVsReactionAve[:, 1])
		rWithoutOutliers = pearsonr(toyaVsReactionAve[idx, 0], toyaVsReactionAve[idx, 1])
		hey = 2
		# plot the results
		plt.figure(figsize=(7, 7), dpi=100)
		ax = plt.axes()
		plt.title(
			f"Central Carbon Metabolism Flux Comparison",
			fontsize=12)
		plt.suptitle(
			f"(averaged over {len(self.all_cells)} cells, each weighted by total time spanned)",
			fontsize=8)
		plt.errorbar(toyaVsReactionAve[:, 1], toyaVsReactionAve[:, 0],
					 xerr=toyaVsReactionAve[:, 3], yerr=toyaVsReactionAve[:, 2], fmt=".",
					 ecolor="k", alpha=0.5, linewidth=0.5)
		ylim = plt.ylim()
		plt.plot([ylim[0], ylim[1]], [ylim[0], ylim[1]], color="k")
		plt.plot(toyaVsReactionAve[:, 1], toyaVsReactionAve[:, 0], "ob", markeredgewidth=0.1,
				 alpha=0.9)
		plt.xlabel("Toya 2010 Reaction Flux [mmol/g/hr]")
		plt.ylabel(f"Mean WCM (Simulation ID: {self.sim_id}) Reaction Flux [mmol/g/hr]")

		# noinspection PyTypeChecker
		ax.set_xlim([-20, 30])
		xlim = ax.get_xlim()
		ylim = ax.get_ylim()
		ax.set_yticks(list(range(int(ylim[0]), int(ylim[1]) + 1, 10)))
		ax.set_xticks(list(range(int(xlim[0]), int(xlim[1]) + 1, 10)))
		ax.text(0.5, -0.08, "Pearson R = %.4f, p = %s\n(%.4f, %s without outliers)" % (
			rWithAll[0], rWithAll[1], rWithoutOutliers[0], rWithoutOutliers[1]), fontsize=8,
				transform=ax.transAxes,
				ha='center', va='top')

		plotOutName = plotOutFileName + "_FBA_avgeraging_method_1b_" + self.sim_id + ".png"
		exportFigure(plt, plotOutDir, plotOutName, metadata)
		plt.close("all")
		hi = 5


	# enzyme method 1
	def test_order_plot(self, simDataFile):
		"""
		Plot the order of the reactions for each monomer in the monomer_to_reactions_dict.
		"""
		sim_data = self.read_pickle_file(simDataFile)
		simOutDir = os.path.join(self.all_cells[0], "simOut")

		# get the base reaction IDs:
		fbaResults = TableReader(os.path.join(simOutDir, "FBAResults"))
		base_reaction_IDs = fbaResults.readAttribute("base_reaction_ids") # can also access in sim data, but not sure if the mapping is consistent but the lengths are the same, so maybe? # todo: ask riley
		base_reaction_ids = sim_data.process.metabolism.base_reaction_ids
		if base_reaction_IDs == base_reaction_ids:
			print("base reaction ids are the same and in the same order.")

		# other ways to get the reaction ids: (already said the base one)
		all_reaciton_ids = list(sim_data.process.metabolism.reaction_id_to_base_reaction_id.keys())
		reaction_stoich_ids = list(sim_data.process.metabolism.reaction_stoich.keys())
		# see what reaction IDs are in the all reaction ids of the reaction_stoich_ids:
		reaction_stoich_ids_set = set(reaction_stoich_ids).intersection(set(all_reaciton_ids))
		# see which do not match:
		reaction_stoich_ids_not_in_all = set(all_reaciton_ids).difference(set(reaction_stoich_ids)) # I believe the 21 that do not match have to do with aa

		# reaction ids:
		reaction_ids = reaction_stoich_ids
		reaction_IDs = fbaResults.readAttribute("reactionIDs")
		if reaction_ids == reaction_IDs:
			print("reaction ids are the same and in the same order.") # Seems like this one is not printed, so pressumably I cannot trust this order!

		# get the enzyme IDs:
		# no this one I do not trust to get from sim data.
		hi = 6



	# todo: just plot all reactions by the other version of all reactions



	def get_reactions(self, monomer_to_reactions_dict, reactions_list_to_search):

		monomers_found_in_reaction_list = {}
		for monomer in monomer_to_reactions_dict.keys():
			# check if the monomer is in the reactions list to search:
			reactions_for_monomer = []
			if monomer_to_reactions_dict[monomer] != []:
				for reaction in monomer_to_reactions_dict[monomer]:
					if reaction in reactions_list_to_search:
						reactions_for_monomer.append(reaction)
				if reactions_for_monomer != []:
					monomers_found_in_reaction_list[monomer] = reactions_for_monomer

		return monomers_found_in_reaction_list



	def FBA_plot(self, sim_data, plotOutDir, plotOutFileName, validationDataFile, metadata):
		# taken from models/ecoli/analysis/cohort/centralCarbonMetabolismScatter.py
		# todo: I think this was where I meant to weight differently?
		# Get all cells
		allDir = self.ap.get_cells() # todo: decide if I want to discard some of the initial generations

		validation_data = self.read_pickle_file(validationDataFile)
		toyaReactions = validation_data.reactionFlux.toya2010fluxes["reactionID"]
		toyaFluxes = validation_data.reactionFlux.toya2010fluxes["reactionFlux"]
		toyaStdev = validation_data.reactionFlux.toya2010fluxes["reactionFluxStdev"]
		toyaFluxesDict = dict(zip(toyaReactions, toyaFluxes))
		toyaStdevDict = dict(zip(toyaReactions, toyaStdev))

		cellDensity = sim_data.constants.cell_density

		modelFluxes = {}
		toyaOrder = []
		for rxn in toyaReactions:
			modelFluxes[rxn] = []
			toyaOrder.append(rxn)

		mmol_per_g_per_h = units.mmol / units.g / units.h
		# I bleieve this one is different averaging wise becuase it should be doing STD and flux averaging by the generation
		for simDir in allDir:
			simOutDir = os.path.join(simDir, "simOut")

			massListener = TableReader(os.path.join(simOutDir, "Mass"))
			cellMass = massListener.readColumn("cellMass")
			dryMass = massListener.readColumn("dryMass")
			coefficient = dryMass / cellMass * cellDensity.asNumber(MASS_UNITS / VOLUME_UNITS)

			fbaResults = TableReader(os.path.join(simOutDir, "FBAResults"))
			base_reaction_ids = fbaResults.readAttribute("base_reaction_ids")
			base_reaction_fluxes = (COUNTS_UNITS / MASS_UNITS / TIME_UNITS) * (
						fbaResults.readColumn("base_reaction_fluxes").T / coefficient).T
			base_reaction_id_to_index = {
				rxn_id: i for (i, rxn_id) in enumerate(base_reaction_ids)
			}

			reaction_fluxes = fbaResults.readColumn("reactionFluxes")
			base_reaction_fluxes1 = fbaResults.readColumn("base_reaction_fluxes")

			reaction_IDs = fbaResults.readAttribute("reactionIDs")
			kineticTargetFluxNames = fbaResults.readAttribute("kineticTargetFluxNames")
			overlap_v_k = set(toyaReactions).intersection(set(kineticTargetFluxNames))
			overlap_v_r = set(toyaReactions).intersection(set(reaction_IDs))
			overlap_k_r = set(kineticTargetFluxNames).intersection(set(reaction_IDs))
			overlap_r_b = set(base_reaction_ids).intersection(set(reaction_IDs))
			overlap_k_b = set(base_reaction_ids).intersection(set(kineticTargetFluxNames))
			overlap_r_b_unique = set(set(base_reaction_ids).intersection(set(reaction_IDs)))
			overlap_b_v = set(toyaReactions).intersection(set(base_reaction_ids)) # all 23
			no_overlap_v_r = set(toyaReactions).difference(set(reaction_IDs))
			no_overlap_v_r = list(no_overlap_v_r)
			rxn_id_to_base_rxn_id = sim_data.process.metabolism.reaction_id_to_base_reaction_id
			target_value = no_overlap_v_r[0]
			matching_keys = [k for k, v in rxn_id_to_base_rxn_id.items() if v == target_value]
			unmatched_rxns = {}
			for rxn in no_overlap_v_r:
				matching_keys = [k for k, v in rxn_id_to_base_rxn_id.items() if v == rxn]
				unmatched_rxns[rxn] = matching_keys

			for toyaReaction in toyaReactions:
				if toyaReaction in base_reaction_id_to_index:
					rxn_index = base_reaction_id_to_index[toyaReaction]
					fluxTimeCourse = base_reaction_fluxes[:, rxn_index]
					modelFluxes[toyaReaction].append(
						np.mean(fluxTimeCourse).asNumber(mmol_per_g_per_h))
		hi = 7
		toyaVsReactionAve = []
		for rxn, toyaFlux in toyaFluxesDict.items():
			if rxn in modelFluxes:
				toyaVsReactionAve.append(
					(np.mean(modelFluxes[rxn]),
					 toyaFlux.asNumber(mmol_per_g_per_h),
					 np.std(modelFluxes[rxn]), toyaStdevDict[rxn].asNumber(mmol_per_g_per_h)))
		hi = 5
		toyaVsReactionAve = np.array(toyaVsReactionAve)
		idx = np.abs(toyaVsReactionAve[:, 0]) < 5 * np.abs(toyaVsReactionAve[:, 1])
		rWithAll = pearsonr(toyaVsReactionAve[:, 0], toyaVsReactionAve[:, 1])
		rWithoutOutliers = pearsonr(toyaVsReactionAve[idx, 0], toyaVsReactionAve[idx, 1])
		hi = 5
		plt.figure(figsize=(3.5, 3.5))
		ax = plt.axes()
		plt.title(
			"Central Carbon Metabolism Flux, Pearson R = %.4f, p = %s\n(%.4f, %s without outliers)" % (
			rWithAll[0], rWithAll[1], rWithoutOutliers[0], rWithoutOutliers[1]), fontsize=6)
		plt.errorbar(toyaVsReactionAve[:, 1], toyaVsReactionAve[:, 0],
					 xerr=toyaVsReactionAve[:, 3], yerr=toyaVsReactionAve[:, 2], fmt=".",
					 ecolor="k", alpha=0.5, linewidth=0.5)
		ylim = plt.ylim()
		plt.plot([ylim[0], ylim[1]], [ylim[0], ylim[1]], color="k")
		plt.plot(toyaVsReactionAve[:, 1], toyaVsReactionAve[:, 0], "ob", markeredgewidth=0.1,
				 alpha=0.9)
		plt.xlabel("Toya 2010 Reaction Flux [mmol/g/hr]")
		plt.ylabel("Mean WCM Reaction Flux [mmol/g/hr]")
		whitePadSparklineAxis(ax)

		# noinspection PyTypeChecker
		ax.set_xlim([-20, 30])
		xlim = ax.get_xlim()
		ylim = ax.get_ylim()
		ax.set_yticks(list(range(int(ylim[0]), int(ylim[1]) + 1, 10)))
		ax.set_xticks(list(range(int(xlim[0]), int(xlim[1]) + 1, 10)))

		plotOutName = plotOutFileName + "_FBA"
		exportFigure(plt, plotOutDir, plotOutName, metadata)
		plt.close("all")


	def enzyme_kinetics_plot(self, simDataFile, plotOutDir, plotOutFileName, metadata):
		# adapted from models/ecoli/analysis/cohort/kinetics_flux_comparison.py

		# obtain the simulation data:
		sim_data = self.read_pickle_file(simDataFile)
		simOutDir = os.path.join(self.all_cells[0], "simOut")

		# burn in time:
		cell_burnIn_times = []
		for cell in self.all_cells:
			mainListener = TableReader(os.path.join(os.path.join(cell, "simOut"), "Main"))
			time = mainListener.readColumn("time")
			mainListener.close()
			burnIn = time > BURN_IN_TIME
			burnIn[0] = False
			cell_burnIn_times.append(burnIn)

		# concatenate the burn in times:
		burnIn = np.concatenate(cell_burnIn_times)

		# get the coefficent for the fluxes: # todo: check that this is an ok way to find the coefficent and everything
		cellMass = read_stacked_columns(self.all_cells, "Mass", "cellMass")
		dryMass = read_stacked_columns(self.all_cells, "Mass", "dryMass")
		cellDensity = sim_data.constants.cell_density
		coefficient = dryMass / cellMass * cellDensity.asNumber(MASS_UNITS / VOLUME_UNITS)

		# obtain enzyme kinetics data:
		enzymeKineticsReader = TableReader(os.path.join(simOutDir, "EnzymeKinetics"))
		kineticsConstrainedReactions = np.array(
			enzymeKineticsReader.readAttribute("kineticsConstrainedReactions"))
		constraint_is_kcat_only = np.array(
			enzymeKineticsReader.readAttribute('constraint_is_kcat_only')) # todo: is this needed?

		# make a dictionary of the reaction IDs and their indices:
		kinetic_rxn_id_to_idx = {rxn: i for i, rxn in enumerate(kineticsConstrainedReactions)}



		hi = 5

		# get the enzyme kinetics fluxes using read stacked columns:
		allTargetFluxes = (COUNTS_UNITS / MASS_UNITS / TIME_UNITS) * (
				read_stacked_columns(self.all_cells, "EnzymeKinetics", "targetFluxes").T / coefficient.squeeze()).T
		allActualFluxes = (COUNTS_UNITS / MASS_UNITS / TIME_UNITS) * (
					read_stacked_columns(self.all_cells, "EnzymeKinetics", "actualFluxes").T / coefficient.squeeze()).T

		allTargetFluxes = allTargetFluxes.asNumber(units.mmol / units.g / units.h)
		allActualFluxes = allActualFluxes.asNumber(units.mmol / units.g / units.h)

		allTargetAve = np.nanmean(allTargetFluxes[burnIn, :], axis=0)
		allActualAve = np.nanmean(allActualFluxes[burnIn, :], axis=0)

		# todo: what is this even doing?
		n_kinetic_constrained_reactions = len(kineticsConstrainedReactions)

		# boundary target fluxes
		# TODO: what are these doing?
		boundaryTargetAve = allTargetAve[n_kinetic_constrained_reactions:]
		boundaryActualAve = allActualAve[n_kinetic_constrained_reactions:]

		# kinetic target fluxes
		targetAve = allTargetAve[:n_kinetic_constrained_reactions]
		actualAve = allActualAve[:n_kinetic_constrained_reactions]
		hi = 5

		kcatOnlyReactions = constraint_is_kcat_only
		kmAndKcatReactions = ~constraint_is_kcat_only

		# categorize how well the actual flux matches the target flux
		thresholds = [2, 10]  # TODO: might want to play around with this
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
		for reaction, target, flux, category in zip(
				kineticsConstrainedReactions[kmAndKcatReactions], targetAve[kmAndKcatReactions],
				actualAve[kmAndKcatReactions], categorization[kmAndKcatReactions]):
			output.writerow([reaction, target, flux, category])

		output.writerow(["kcat only"])
		for reaction, target, flux, category in zip(
				kineticsConstrainedReactions[kcatOnlyReactions], targetAve[kcatOnlyReactions],
				actualAve[kcatOnlyReactions], categorization[kcatOnlyReactions]):
			output.writerow([reaction, target, flux, category])
		csvFile.close()

		# add small number to allow plotting of 0 flux on log scale
		targetAve += 1e-6
		actualAve += 1e-6

		pearsonAll = pearsonr(np.log10(targetAve), np.log10(actualAve))
		pearsonNoZeros = pearsonr(np.log10(targetAve[(categorization != -2)]),
								  np.log10(actualAve[(categorization != -2)]))

		# plot data
		plt.figure(figsize=(7, 7), dpi=100)
		ax = plt.axes()
		plt.plot([-6, 4], [-6, 4], 'k', linewidth=0.75)
		plt.plot([-5, 4], [-6, 3], 'k', linewidth=0.5)
		plt.plot([-6, 3], [-5, 4], 'k', linewidth=0.5)
		plt.plot(np.log10(targetAve), np.log10(actualAve), 'o', color="black", markersize=8,
				 alpha=0.15, zorder=1, markeredgewidth=0.0)
		plt.plot(np.log10(boundaryTargetAve), np.log10(boundaryActualAve), "ob", color="red",
				 markeredgewidth=0.25, alpha=0.9, label='boundary fluxes')
		plt.title(f"Enzyme Kinetics Flux Comparison\n(averaged over {len(self.all_cells)} cells, each weighted by simulation time spanned, \nSimulation ID: {self.sim_id})")
		plt.xlabel("Log10(Target Flux [mmol/g/hr])")
		plt.ylabel("Log10(Actual Flux [mmol/g/hr])")
		ax.text(0.5, -0.8,"PCC = %.3f, p = %s\n(%.3f, p = %s without points at zero)" % (
		pearsonAll[0], pearsonAll[1], pearsonNoZeros[0], pearsonNoZeros[1]))
		plt.minorticks_off()
		whitePadSparklineAxis(ax)
		xlim = ax.get_xlim()
		ylim = ax.get_ylim()
		ax.set_ylim(ylim[0] - 0.5, ylim[1])
		ax.set_xlim(xlim[0] - 0.5, xlim[1])
		ax.set_yticks(list(range(-6, int(ylim[1]) + 1, 2)))
		ax.set_xticks(list(range(-6, int(xlim[1]) + 1, 2)))



		if self.relevant_reactions_as_catalysts is not None:
			catalyst_monomers_found_in_kinetic_reactions = self.get_reactions(
				self.relevant_reactions_as_catalysts, kineticsConstrainedReactions)

			if catalyst_monomers_found_in_kinetic_reactions is not None:
				for monomer in catalyst_monomers_found_in_kinetic_reactions.keys():
					for reaction in catalyst_monomers_found_in_kinetic_reactions[monomer]:
						rxn_idx = kinetic_rxn_id_to_idx[reaction]
						plt.plot(np.log10(targetAve[rxn_idx]), np.log10(actualAve[rxn_idx]),'o', color="blue", markersize=8, alpha=0.9, markeredgewidth=0.25,label='catalyst monomers' if rxn_idx == 0 else "")
						# add the monomer names to the plot
						# todo: do I want the reaction name or the monomer name to be the text name?
						ax.text(np.log10(targetAve[rxn_idx])+0.5, np.log10(actualAve[rxn_idx]), reaction ,ha='center', va='bottom', fontsize=8, rotation=0, color="blue", zorder=2)

		# make a table at the bottom of the graph with the important reactions


		# add a legend and save the figure
		plt.legend(loc='upper left', fontsize=8, frameon=False)
		plotOutFileName = plotOutFileName + "_enzyme_kinetics_" + self.sim_id + ".png"
		exportFigure(plt, plotOutDir, plotOutFileName)



				# for monomer, reactions in catalyst_monomers_found_in_kinetic_reactions.items():
				# 	for reaction in reactions:
				# 		words = catalyst_monomers_found_in_kinetic_reactions[monomer]
				# 		plt.annotate(words, xy=(np.log10(targetAve[reaction]), np.log10(actualAve[reaction])),
				# 					 xytext=(5, 5), textcoords='offset points', fontsize=8, color='blue')
				# 		#ax.text(x + .2, y, name, ha='center', va='bottom', fontsize=8,rotation=0, )



		# ax.legend()
		#
		# plotOutFileName = plotOutFileName + "_enzyme_kinetics_" + self.sim_id + ".png"
		# exportFigure(plt, plotOutDir, plotOutFileName)









	def do_plot(self, variantDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		# Get simulation information:
		self.sim_id = metadata["description"]
		self.n_total_gens = self.ap.n_generation
		self.generations_to_plot = np.arange(IGNORE_FIRST_N_GENERATIONS, self.n_total_gens)
		self.all_cells = self.ap.get_cells(generation=self.generations_to_plot, only_successful=True)

		allDir = self.ap.get_cells() # todo: delete later if needed

		# get simulation data:
		sim_data = self.read_pickle_file(simDataFile)

		targetFluxList = []
		actualFluxList = []

		for simDir in allDir:
			simOutDir = os.path.join(simDir, "simOut")

			mainListener = TableReader(os.path.join(simOutDir, "Main"))
			time = mainListener.readColumn("time")
			mainListener.close()
			burnIn = time > BURN_IN_TIME
			hi = 5
			burnIn[0] = False
			hi =5

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
		allTargetAve = np.nanmean(targetFluxList, axis = 0) # todo: ask if this is how it should be averaged? or if it should be averaged over all time points?
		allActualAve = np.nanmean(actualFluxList, axis = 0)

		# boundary target fluxes
		boundaryTargetAve = allTargetAve[n_kinetic_constrained_reactions:]
		boundaryActualAve = allActualAve[n_kinetic_constrained_reactions:]

		# kinetic target fluxes
		targetAve = allTargetAve[:n_kinetic_constrained_reactions]
		actualAve = allActualAve[:n_kinetic_constrained_reactions]

		self.relevant_reactions_as_catalysts, self.relevant_reactions_as_substrates = self.get_relevant_monomers(simDataFile)

		# see if the monomers show up in any of the kinetic constrained reactions:
		#catalyst_monomers_found_in_kinetic_reactions = self.get_reactions(self.relevant_reactions_as_catalysts, kineticsConstrainedReactions)
		#substrate_monomers_found_in_kinetic_reactions = self.get_reactions(self.relevant_reactions_as_substrates, kineticsConstrainedReactions)




		# plot the FBA results
		self.FBA_plot(sim_data, plotOutDir, plotOutFileName, validationDataFile, metadata)
		self.test_order_plot(simDataFile) # todo: since reaction ids are not always in the same order, which should I use?
		self.FBA_plots_averaging_method_1a(simDataFile, validationDataFile, plotOutFileName, plotOutDir, metadata)
		self.FBA_plots_averaging_method_1b(simDataFile, validationDataFile, plotOutFileName, plotOutDir, metadata)

		self.enzyme_kinetics_plot(simDataFile, plotOutDir, plotOutFileName, metadata)

		hi = 5
		self.centralCarbonMetabolismScatterPlot(plotOutDir, plotOutFileName, simDataFile,
										   validationDataFile, metadata)

		hi = 5

		# check if any of the relevant reactions are in the fba reactions:


		# categorize reactions that use constraints with only kcat, Km and kcat, or switch between both types of constraints
		kcatOnlyReactions = constraint_is_kcat_only
		kmAndKcatReactions = ~constraint_is_kcat_only
		hi = 5
		# categorize how well the actual flux matches the target flux
		thresholds = [2, 10]
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

		csvFile.close()

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
