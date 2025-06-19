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
IMPORTANT_MONOMERS = ["ISOCIT-LYASE-MONOMER", "G6980-MONOMER", "BASS-MONOMER", 'G6986-MONOMER[c]']

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
					reaction_to_name_info[reaction] = name

				# append the reactions to the dictionary:
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
					reaction_to_name_info[reaction] = name

				# append the reactions to the dictionary:
				relevant_reactions_as_substrates[monomer] = reaction_to_name_info


		return relevant_reactions_as_catalysts, relevant_reactions_as_substrates



	def obtain_FBA_results(self, simDataFile):
		"""
		Obtain FBA results from the simDataFile.
		"""
		sim_data = self.read_pickle_file(simDataFile)
		fba = sim_data.process.metabolism.fba
		fba_reaction_ids = fba.getReactionIDs()
		base_reaction_ids = sim_data.process.metabolism.base_reaction_ids
		fba_reaction_ids_to_base_reaction_ids = sim_data.process.metabolism.reaction_id_to_base_reaction_id

		return fba_reaction_ids, base_reaction_ids, fba_reaction_ids_to_base_reaction_ids

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

		relevant_reactions_as_catalysts, relevant_reactions_as_substrates = self.get_relevant_monomers(simDataFile)

		# see if the monomers show up in any of the kinetic constrained reactions:
		catalyst_monomers_found_in_kinetic_reactions = self.get_reactions(relevant_reactions_as_catalysts, kineticsConstrainedReactions)
		substrate_monomers_found_in_kinetic_reactions = self.get_reactions(relevant_reactions_as_substrates, kineticsConstrainedReactions)

		hi = 5

		# check if any of the relevant reactions are in the fba reactions:


		# categorize reactions that use constraints with only kcat, Km and kcat, or switch between both types of constraints
		kcatOnlyReactions = constraint_is_kcat_only
		kmAndKcatReactions = ~constraint_is_kcat_only

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
