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

# ignore data from metabolism burnin period
BURN_IN_TIME = 1
IMPORTANT_MONOMERS = ["ISOCIT-LYASE-MONOMER", "G6980-MONOMER", "BASS-MONOMER"]

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


	def get_catalyst_monomers(self, simDataFile):
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


	hi = 6
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

	def generate_WCM_implementation_df(self):
		# first, map the gene names to protein names according to Ecocyc:
		ecocyc_proteins_to_genes_df = pd.read_excel(ecocyc_proteins_to_genes_xlsx_path)
		gene_to_monomer_dict = {}
		for row in range(len(ecocyc_proteins_to_genes_df["Proteins"])):
			gene = ecocyc_proteins_to_genes_df["Genes"][row]
			protein = ecocyc_proteins_to_genes_df["Proteins"][row]
			gene_to_monomer_dict[gene] = protein

		# Add the proteins to the data WCM gene implementation data table:
		WCM_gene_implementation_df = pd.read_excel(WCM_metabolism_gene_xlsx_path)
		WCM_gene_implementation_df["Monomer ID"] = None
		gene_IDs = WCM_gene_implementation_df["Gene ID (EcoCyc)"]
		for gene_ID in gene_IDs:
			idx = WCM_gene_implementation_df.index[
				WCM_gene_implementation_df['Gene ID (EcoCyc)'] == gene_ID].tolist()
			hi = 6
			# TODO: left off here
			protein_ID = gene_to_monomer_dict.get(gene_ID)
			WCM_gene_implementation_df.loc[idx, "Monomer ID"] = protein_ID


		self.WCM_metabolic_protein_implementation_df = WCM_gene_implementation_df[["Monomer ID", "Gene name", "Macklin et al. (2020)", "Latest version (20220602)"]]

	hi = 5
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
		#self.constrainedReactions_in_load_data = np.array(enzymeKineticsReader.readAttribute("constrainedReactions")) # length: 415
		enzymeIDs = enzymeKineticsReader.readAttribute("enzymeIDs") # lenth: 307 (maybe multiple reactions use the same enzyme??)
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
		hi = 6
		enzymeIDs = [i[:-3] for i in enzymeIDs]
		self.shared_kinetic_enzymes = set(enzymeIDs) & set(kinetic_reaction_enzymes) # all the enzymeIDs in the reader exist in the kinetics reaction enzymes list!
		self.non_shared_kinetic_enzymes = set(enzymeIDs) ^ set(kinetic_reaction_enzymes)

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

		# check if any of the metabolism genes implemented are here:
		self.generate_WCM_implementation_df()
		# todo: go to the break point in the function and see if you can see the 415 substrate ids!

		# return the reactions with relevancy:
		return relevant_reactions, kinetic_reaction_IDs, kinetic_reaction_to_catalysts, kinetic_reaction_enzymes,reaction_IDs




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

		get_catalyst_monomers_dict = self.get_catalyst_monomers(simDataFile)
		hi = 5

		# determine if any of the reactions plotted contain an important monomer:
		relevant_reactions, kinetic_reaction_IDs, kinetic_reaction_to_catalysts, kinetic_reaction_enzymes, reaction_IDs = self.load_data(simDataFile)

		# check if any of the relevant reactions overlap with the kinetic reactions:
		kinetic_reactions_with_important_monomers = {}
		for key in relevant_reactions.keys():
			if relevant_reactions[key] != []:
				for reaction in relevant_reactions[key]:
					# check if the reaction is in the kinetic reactions:
					if reaction in kineticsConstrainedReactions:
						kinetic_reactions_with_important_monomers[reaction] = key # make it equal to the monomer

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
