"""
Save the average new gene protein counts and doubling time,
This plot is intended to be run on simulations where the new gene
option was enabled.
"""

import numpy as np
from Bio.Data.IUPACData import avg_ambiguous_dna_weights
from fontTools.merge.util import avg_int
from matplotlib import pyplot as plt
from unum.units import fg

import pickle
import os

from models.ecoli.analysis import variantAnalysisPlot
from wholecell.analysis.analysis_tools import (
	exportFigure, read_stacked_columns, read_stacked_bulk_molecules,
	stacked_cell_identification)
from wholecell.io.tablereader import TableReader

# Remove first N gens from plot
IGNORE_FIRST_N_GENS = 16

INTERACTIVE_PLOT = True

class Plot(variantAnalysisPlot.VariantAnalysisPlot):
	def do_plot(self, inputDir, plotOutDir, plotOutFileName, simDataFile,
				validationDataFile, metadata):

		variants = self.ap.get_variants()

		# Determine new gene ids
		with open(simDataFile, 'rb') as f:
			sim_data = pickle.load(f)
		mRNA_sim_data = sim_data.process.transcription.cistron_data.struct_array
		monomer_sim_data = sim_data.process.translation.monomer_data.struct_array
		new_gene_mRNA_ids = mRNA_sim_data[mRNA_sim_data['is_new_gene']]['id'].tolist()
		mRNA_monomer_id_dict = dict(zip(
			monomer_sim_data['cistron_id'], monomer_sim_data['id']))
		monomer_mRNA_id_dict = dict(zip(
			monomer_sim_data['id'], monomer_sim_data['cistron_id']))
		new_gene_monomer_ids = [mRNA_monomer_id_dict.get(mRNA_id)
								for mRNA_id in new_gene_mRNA_ids]
		if len(new_gene_mRNA_ids) == 0:
			print("This plot is intended to be run on simulations where the"
				  " new gene option was enabled, but no new gene mRNAs were "
				  "found.")
			return
		if len(new_gene_monomer_ids) == 0:
			print("This plot is intended to be run on simulations where the "
				  "new gene option was enabled, but no new gene proteins "
				  "were found.")
			return
		assert len(new_gene_monomer_ids) == len(new_gene_mRNA_ids), \
			'number of new gene monomers and mRNAs should be equal'

		# Data extraction
		plot_variant_mask = np.full(len(variants), True)
		doubling_times = np.zeros(len(variants))
		avg_ng_monomer = np.zeros(len(variants))
		avg_ng_monomer_proteome_mass_fraction = np.zeros(len(variants))
		avg_active_ribosome_counts = np.zeros(len(variants))
		avg_inactive_ribosome_counts = np.zeros(len(variants))
		avg_total_ribosome_counts = np.zeros(len(variants))
		avg_total_ribosome_conc = np.zeros(len(variants))
		avg_active_rnap_counts = np.zeros(len(variants))
		avg_inactive_rnap_counts = np.zeros(len(variants))
		avg_total_rnap_counts = np.zeros(len(variants))
		avg_total_rnap_conc = np.zeros(len(variants))
		avg_dry_mass = np.zeros(len(variants))
		avg_cell_mass = np.zeros(len(variants))
		avg_mRNA_mass = np.zeros(len(variants))
		avg_rRNA_mass = np.zeros(len(variants))
		avg_tRNA_mass = np.zeros(len(variants))
		avg_protein_mass = np.zeros(len(variants))
		avg_dna_mass = np.zeros(len(variants))
		avg_water_mass = np.zeros(len(variants))
		avg_small_molecule_mass = np.zeros(len(variants))
		avg_membrane_mass = np.zeros(len(variants))
		percent_completion = np.zeros(len(variants))
		avg_ribosomal_protein_counts = {}
		avg_rnap_subunit_counts = {}
		avg_rRNA_counts = {}

		avg_s30_limiting_protein_counts = np.zeros(len(variants))
		avg_s30_16s_rRNA_total_counts = np.zeros(len(variants))
		avg_s50_limiting_protein_counts = np.zeros(len(variants))
		avg_s50_23s_rRNA_total_counts = np.zeros(len(variants))
		avg_s50_5s_rRNA_total_counts = np.zeros(len(variants))
		avg_s30_total_counts = np.zeros(len(variants))
		avg_s50_total_counts = np.zeros(len(variants))

		avg_s30_protein_counts = {}
		avg_s50_protein_counts = {}

		n_total_gens = self.ap.n_generation

		min_variant = min(variants)

		for i, variant in enumerate(variants):

			all_cells = self.ap.get_cells(
				variant=[variant],
				generation= np.arange(IGNORE_FIRST_N_GENS, n_total_gens),
				only_successful=True)

			if len(all_cells) == 0:
				plot_variant_mask[i] = False
				continue

			# Doubling times
			dt = read_stacked_columns(
				all_cells, 'Main', 'time',
				fun=lambda x: (x[-1] - x[0]) / 60.).squeeze()
			doubling_times[i] = np.mean(dt)
			if variant == min_variant:
				sim_dir = all_cells[0]
				simOutDir = os.path.join(sim_dir, 'simOut')

				# Extract protein indexes for each new gene
				monomer_counts_reader = TableReader(os.path.join(
					simOutDir, "MonomerCounts"))
				monomer_idx_dict = {monomer: i for i, monomer in
					enumerate(monomer_counts_reader.readAttribute('monomerIds'))}
				new_gene_monomer_indexes = [monomer_idx_dict.get(monomer_id)
					for monomer_id in new_gene_monomer_ids]

				mRNA_cistron_counts_reader = TableReader(os.path.join(
					simOutDir, "RNACounts"))
				mRNA_cistron_ids = mRNA_cistron_counts_reader.readAttribute(
					"mRNA_cistron_ids")
				mRNA_cistron_idx_dict = {
					mRNA_cistron_id: i for i, mRNA_cistron_id in
					enumerate(mRNA_cistron_ids)}

				uniqueMoleculeCounts = TableReader(
					os.path.join(simOutDir, "UniqueMoleculeCounts"))
				ribosome_index = uniqueMoleculeCounts.readAttribute(
					"uniqueMoleculeIds").index('active_ribosome')

				uniqueMoleculeCounts = TableReader(
					os.path.join(simOutDir, "UniqueMoleculeCounts"))
				active_rnap_index = uniqueMoleculeCounts.readAttribute(
					"uniqueMoleculeIds").index('active_RNAP')

				# Ribosome and RNAP components
				uncharged_tRNA_ids = sim_data.process.transcription.uncharged_trna_names
				charged_tRNA_ids = sim_data.process.transcription.charged_trna_names
				tRNA_cistron_ids = [tRNA_id[:-3] for tRNA_id in uncharged_tRNA_ids]
				rRNA_ids = [
					sim_data.molecule_groups.s30_16s_rRNA[0],
					sim_data.molecule_groups.s50_23s_rRNA[0],
					sim_data.molecule_groups.s50_5s_rRNA[0]]
				rRNA_cistron_ids = [rRNA_id[:-3] for rRNA_id in rRNA_ids]
				ribosomal_subunit_ids = [
					sim_data.molecule_ids.s30_full_complex,
					sim_data.molecule_ids.s50_full_complex]

				rnap_subunit_monomer_ids = sim_data.molecule_groups.RNAP_subunits
				rnap_subunit_monomer_indexes = [
					monomer_idx_dict.get(monomer_id) for monomer_id in rnap_subunit_monomer_ids
				]
				rnap_subunit_mRNA_cistron_ids = [
					monomer_mRNA_id_dict.get(monomer_id) for monomer_id in rnap_subunit_monomer_ids
				]
				rnap_subunit_mRNA_cistron_indexes = [
					mRNA_cistron_idx_dict.get(cistron_id) for cistron_id in rnap_subunit_mRNA_cistron_ids
				]

				ribosomal_protein_monomer_ids = sim_data.molecule_groups.ribosomal_proteins
				ribosomal_protein_monomer_indexes = [
					monomer_idx_dict.get(monomer_id) for monomer_id in ribosomal_protein_monomer_ids
				]
				ribosomal_protein_mRNA_cistron_ids = [
					monomer_mRNA_id_dict.get(monomer_id) for monomer_id in ribosomal_protein_monomer_ids
				]
				ribosomal_protein_mRNA_cistron_indexes = [
					mRNA_cistron_idx_dict.get(cistron_id) for cistron_id in ribosomal_protein_mRNA_cistron_ids
				]


				# Set up data structures
				for id in ribosomal_protein_monomer_ids:
					avg_ribosomal_protein_counts[id + "_mRNA"] = np.zeros(len(variants))
					avg_ribosomal_protein_counts[id + "_monomer"] = np.zeros(len(variants))
				for id in rnap_subunit_monomer_ids:
					avg_rnap_subunit_counts[id + "_mRNA"] = np.zeros(len(variants))
					avg_rnap_subunit_counts[id + "_monomer"] = np.zeros(len(variants))
				for id in rRNA_ids:
					avg_rRNA_counts[id + "_rRNA"] = np.zeros(len(variants))

			counts_to_molar = read_stacked_columns(
				all_cells, 'EnzymeKinetics', 'countsToMolar',
				remove_first=True, ignore_exception=True)

			avg_new_gene_monomer_counts = read_stacked_columns(all_cells,
				'MonomerCounts', 'monomerCounts',
				remove_first=True)[:,new_gene_monomer_indexes]
			avg_ng_monomer[i] = np.mean(avg_new_gene_monomer_counts)

			active_ribosome_counts = read_stacked_columns(
				all_cells, 'UniqueMoleculeCounts',
				'uniqueMoleculeCounts',
				remove_first=True)[:,ribosome_index]
			avg_active_ribosome_counts[i] = np.mean(active_ribosome_counts)

			complex_id_30s = [sim_data.molecule_ids.s30_full_complex]
			complex_id_50s = [sim_data.molecule_ids.s50_full_complex]
			(complex_counts_30s, complex_counts_50s) = read_stacked_bulk_molecules(
				all_cells, (complex_id_30s, complex_id_50s),
				remove_first=True)
			inactive_ribosome_counts = np.minimum(
				complex_counts_30s, complex_counts_50s)
			avg_inactive_ribosome_counts[i] = np.mean(inactive_ribosome_counts)
			avg_total_ribosome_counts[i] = (
					avg_active_ribosome_counts[i] + avg_inactive_ribosome_counts[i])
			avg_total_ribosome_conc[i] = np.mean((active_ribosome_counts + inactive_ribosome_counts) * counts_to_molar.squeeze())

			active_rnap_counts = read_stacked_columns(
				all_cells, 'UniqueMoleculeCounts',
				'uniqueMoleculeCounts',
				remove_first=True)[:,active_rnap_index]

			rnap_id = [sim_data.molecule_ids.full_RNAP]
			(rnapCountsBulk,) = read_stacked_bulk_molecules(
				all_cells, (rnap_id,), remove_first=True)
			avg_active_rnap_counts[i] = np.mean(active_rnap_counts)
			avg_inactive_rnap_counts[i] = np.mean(rnapCountsBulk)
			avg_total_rnap_counts[i] = (
					avg_active_rnap_counts[i] + avg_inactive_rnap_counts[i])
			avg_total_rnap_conc[i] = np.mean((active_rnap_counts + rnapCountsBulk) * counts_to_molar.squeeze())

			# Cell masses
			avg_dry_mass[i] = np.mean(read_stacked_columns(
				all_cells, 'Mass', 'dryMass',
				remove_first=True, ignore_exception=True))
			avg_cell_mass[i] = np.mean(read_stacked_columns(
				all_cells, 'Mass', 'cellMass',
				remove_first=True, ignore_exception=True))
			avg_mRNA_mass[i] = np.mean(read_stacked_columns(
				all_cells, 'Mass', 'mRnaMass',
				remove_first=True, ignore_exception=True))
			avg_rRNA_mass[i] = np.mean(read_stacked_columns(
				all_cells, 'Mass', 'rRnaMass',
				remove_first=True, ignore_exception=True))
			avg_tRNA_mass[i] = np.mean(read_stacked_columns(
				all_cells, 'Mass', 'tRnaMass',
				remove_first=True, ignore_exception=True))
			avg_protein_mass[i] = np.mean(read_stacked_columns(
				all_cells, 'Mass', 'proteinMass',
				remove_first=True, ignore_exception=True))
			avg_dna_mass[i] = np.mean(read_stacked_columns(
				all_cells, 'Mass', 'dnaMass',
				remove_first=True, ignore_exception=True))
			avg_water_mass[i] = np.mean(read_stacked_columns(
				all_cells, 'Mass', 'waterMass',
				remove_first=True, ignore_exception=True))
			avg_small_molecule_mass[i] = np.mean(read_stacked_columns(
				all_cells, 'Mass', 'smallMoleculeMass',
				remove_first=True, ignore_exception=True))
			avg_membrane_mass[i] = np.mean(read_stacked_columns(
				all_cells, 'Mass', 'membrane_mass',
				remove_first=True, ignore_exception=True))

			ng_mass = (
				sim_data.getter.get_mass(
				new_gene_monomer_ids[0])/sim_data.constants.n_avogadro).asNumber(fg)
			avg_ng_monomer_proteome_mass_fraction[i] = (
				avg_ng_monomer[i] * ng_mass) / avg_protein_mass[i]

			# count number of seeds that reached full 24 gens
			# of growth
			all_cells_last_gen = self.ap.get_cells(
				variant=[variant],
				generation=np.arange(n_total_gens - 1, n_total_gens),
				only_successful=True)

			all_cells_first_gen = self.ap.get_cells(
				variant=[variant],
				generation=np.arange(0, 1),
				only_successful=True)

			percent_completion[i] = len(all_cells_last_gen) / len(all_cells_first_gen)

			# Get rRNA counts
			(rRNA_counts, ribosomal_subunit_counts) = read_stacked_bulk_molecules(
				all_cells,
				(rRNA_ids, ribosomal_subunit_ids),
				remove_first=True, ignore_exception=True)
			full_ribosome_counts = read_stacked_columns(
				all_cells, 'UniqueMoleculeCounts', 'uniqueMoleculeCounts',
				remove_first=True, ignore_exception=True)[:, ribosome_index]
			rRNA_counts[:, 0] += ribosomal_subunit_counts[:, 0]
			rRNA_counts[:, 1:] += ribosomal_subunit_counts[:, 1:]
			rRNA_counts += full_ribosome_counts[:, None]
			for j in range(len(rRNA_ids)):
				avg_rRNA_counts[rRNA_ids[j] + "_rRNA"][i] = np.mean(rRNA_counts[:, j])

			# Get ribosomal protein counts
			ribosomal_protein_counts = read_stacked_columns(
				all_cells, 'MonomerCounts', 'monomerCounts',
				remove_first=True)[:, ribosomal_protein_monomer_indexes]
			ribosomal_protein_cistron_counts = read_stacked_columns(
				all_cells, 'RNACounts', 'mRNA_cistron_counts',
				remove_first=True)[:, ribosomal_protein_mRNA_cistron_indexes]
			for j, id in enumerate(ribosomal_protein_monomer_ids):
				avg_ribosomal_protein_counts[id + "_mRNA"][i] = np.mean(
					ribosomal_protein_cistron_counts[:, j])
				avg_ribosomal_protein_counts[id + "_monomer"][i] = np.mean(
					ribosomal_protein_counts[:, j])

			# Get RNAP subunit counts
			rnap_subunit_counts = read_stacked_columns(
				all_cells, 'MonomerCounts', 'monomerCounts',
				remove_first=True)[:, rnap_subunit_monomer_indexes]
			rnap_subunit_cistron_counts = read_stacked_columns(
				all_cells, 'RNACounts', 'mRNA_cistron_counts',
				remove_first=True)[:, rnap_subunit_mRNA_cistron_indexes]
			for j, id in enumerate(rnap_subunit_monomer_ids):
				avg_rnap_subunit_counts[id + "_mRNA"][i] = np.mean(
					rnap_subunit_cistron_counts[:, j])
				avg_rnap_subunit_counts[id + "_monomer"][i] = np.mean(
					rnap_subunit_counts[:, j])

			# From ribosome_components.py
			# Load IDs of ribosome components from sim_data
			s30_protein_ids = sim_data.molecule_groups.s30_proteins
			s30_16s_rRNA_ids = sim_data.molecule_groups.s30_16s_rRNA
			s30_full_complex_id = [sim_data.molecule_ids.s30_full_complex]
			s50_protein_ids = sim_data.molecule_groups.s50_proteins
			s50_23s_rRNA_ids = sim_data.molecule_groups.s50_23s_rRNA
			s50_5s_rRNA_ids = sim_data.molecule_groups.s50_5s_rRNA
			s50_full_complex_id = [sim_data.molecule_ids.s50_full_complex]

			# Get complexation stoichiometries of ribosomal proteins
			complexation = sim_data.process.complexation
			s30_monomers = complexation.get_monomers(s30_full_complex_id[0])
			s50_monomers = complexation.get_monomers(s50_full_complex_id[0])
			s30_subunit_id_to_stoich = {
				subunit_id: stoich for (subunit_id, stoich)
				in zip(s30_monomers['subunitIds'], s30_monomers['subunitStoich'])
				}
			s50_subunit_id_to_stoich = {
				subunit_id: stoich for (subunit_id, stoich)
				in zip(s50_monomers['subunitIds'], s50_monomers['subunitStoich'])
				}
			s30_protein_stoich = np.array([
				s30_subunit_id_to_stoich[subunit_id]
				for subunit_id in s30_protein_ids
				])
			s50_protein_stoich = np.array([
				s50_subunit_id_to_stoich[subunit_id]
				for subunit_id in s50_protein_ids
				])

			# Read free counts of rRNAs
			time = read_stacked_columns(all_cells, 'Main', 'time')
			(s30_16s_rRNA_counts, s50_23s_rRNA_counts, s50_5s_rRNA_counts
			) = read_stacked_bulk_molecules(
				all_cells,
				(s30_16s_rRNA_ids, s50_23s_rRNA_ids, s50_5s_rRNA_ids))

			# Read counts of free 30S and 50S subunits
			(s30_full_complex_counts, s50_full_complex_counts
			) = read_stacked_bulk_molecules(
				all_cells, (s30_full_complex_id, s50_full_complex_id))

			# Read counts of active ribosomes
			simOutDir = os.path.join(all_cells[0], 'simOut')
			unique_molecule_counts_reader = TableReader(
				os.path.join(simOutDir, 'UniqueMoleculeCounts'))
			active_ribosome_index = unique_molecule_counts_reader.readAttribute(
				'uniqueMoleculeIds').index('active_ribosome')
			active_ribosome_counts = read_stacked_columns(
				all_cells, 'UniqueMoleculeCounts', 'uniqueMoleculeCounts',
				)[:, active_ribosome_index]

			# Read counts of free and complexed ribosomal proteins (this includes
			# ribosomal proteins in inactive ribosomal subunits and active
			# ribosomes)
			if variant == min_variant:
				simOutDir = os.path.join(all_cells[0], 'simOut')
				monomer_counts_reader = TableReader(
					os.path.join(simOutDir, 'MonomerCounts'))
				monomer_ids = monomer_counts_reader.readAttribute('monomerIds')
				monomer_id_to_index = {
					monomer_id: i for (i, monomer_id) in enumerate(monomer_ids)
					}
				s30_protein_indexes = np.array([
					monomer_id_to_index[protein_id] for protein_id in s30_protein_ids
					])
				s50_protein_indexes = np.array([
					monomer_id_to_index[protein_id] for protein_id in s50_protein_ids
					])

				# Set up data structures
				for j in range(len(s30_protein_ids)):
					avg_s30_protein_counts['s30_' + s30_protein_ids[j]] = np.zeros(len(variants))
				for j in range(len(s50_protein_ids)):
					avg_s50_protein_counts['s50_' + s50_protein_ids[j]] = np.zeros(len(variants))

			monomer_counts = read_stacked_columns(
				all_cells, 'MonomerCounts', 'monomerCounts',
				)
			# Divide by complexation stoichiometric constant
			s30_protein_counts = monomer_counts[:, s30_protein_indexes] / s30_protein_stoich
			s50_protein_counts = monomer_counts[:, s50_protein_indexes] / s50_protein_stoich

			# Calculate total counts of all components
			s30_limiting_protein_counts = s30_protein_counts.min(axis=1)
			s30_16s_rRNA_total_counts = (
					s30_16s_rRNA_counts.sum(axis=1)
					+ s30_full_complex_counts + active_ribosome_counts)
			s50_limiting_protein_counts = s50_protein_counts.min(axis=1)
			s50_23s_rRNA_total_counts = (
					s50_23s_rRNA_counts.sum(axis=1)
					+ s50_full_complex_counts + active_ribosome_counts)
			s50_5s_rRNA_total_counts = (
					s50_5s_rRNA_counts.sum(axis=1)
					+ s50_full_complex_counts + active_ribosome_counts)

			# Calculate total counts of both subunits
			s30_total_counts = s30_full_complex_counts + active_ribosome_counts
			s50_total_counts = s50_full_complex_counts + active_ribosome_counts

			avg_s30_limiting_protein_counts[i] = np.mean(s30_limiting_protein_counts)
			avg_s30_16s_rRNA_total_counts[i] = np.mean(s30_16s_rRNA_total_counts)
			avg_s50_limiting_protein_counts[i] = np.mean(s50_limiting_protein_counts)
			avg_s50_23s_rRNA_total_counts[i] = np.mean(s50_23s_rRNA_total_counts)
			avg_s50_5s_rRNA_total_counts[i] = np.mean(s50_5s_rRNA_total_counts)
			avg_s30_total_counts[i] = np.mean(s30_total_counts)
			avg_s50_total_counts[i] = np.mean(s50_total_counts)

			for j, id in enumerate(s30_protein_ids):
				avg_s30_protein_counts['s30_' + id][i] = np.mean(s30_protein_counts[:, j])
			for j, id in enumerate(s50_protein_ids):
				avg_s50_protein_counts['s50_' + id][i] = np.mean(s50_protein_counts[:, j])

		# Save ribosome and RNAP component counts
		ribosome_and_rnap_components = np.array(variants)
		ribosome_and_rnap_components_header = np.array(['Variants'])
		for id in ribosomal_protein_monomer_ids:
			ribosome_and_rnap_components = np.vstack((ribosome_and_rnap_components,
				avg_ribosomal_protein_counts[id + "_mRNA"]))
			ribosome_and_rnap_components = np.vstack((ribosome_and_rnap_components,
				avg_ribosomal_protein_counts[id + "_monomer"]))
			ribosome_and_rnap_components_header = np.append(
				ribosome_and_rnap_components_header, "Avg mRNA Counts (Ribosomal): " + id)
			ribosome_and_rnap_components_header = np.append(
				ribosome_and_rnap_components_header, "Avg Monomer Counts (Ribosomal): " + id)
		for id in rnap_subunit_monomer_ids:
			ribosome_and_rnap_components = np.vstack((ribosome_and_rnap_components,
				avg_rnap_subunit_counts[id + "_mRNA"]))
			ribosome_and_rnap_components = np.vstack((ribosome_and_rnap_components,
				avg_rnap_subunit_counts[id + "_monomer"]))
			ribosome_and_rnap_components_header = np.append(
				ribosome_and_rnap_components_header, "Avg mRNA Counts (RNAP): " + id)
			ribosome_and_rnap_components_header = np.append(
				ribosome_and_rnap_components_header, "Avg Monomer Counts (RNAP): " + id)
		for id in rRNA_ids:
			ribosome_and_rnap_components = np.vstack((ribosome_and_rnap_components,
				avg_rRNA_counts[id + "_rRNA"]))
			ribosome_and_rnap_components_header = np.append(
				ribosome_and_rnap_components_header, "Avg rRNA Counts: " + id)
		ribosome_and_rnap_components_header = np.reshape(ribosome_and_rnap_components_header, (-1, 1))
		all_avg_component_counts = np.vstack((
			ribosome_and_rnap_components_header.T, ribosome_and_rnap_components.T))
		np.savetxt(os.path.join(plotOutDir, plotOutFileName + "_ribosome_and_rnap_components.csv"),
				   all_avg_component_counts, delimiter=",", fmt="%s")

		# Save ribosome component counts by subunit
		ribosome_components_by_subunit = np.array(variants)
		ribosome_components_by_subunit = np.vstack((
			ribosome_components_by_subunit,
			avg_s30_total_counts,
			avg_s30_limiting_protein_counts,
			avg_s30_16s_rRNA_total_counts))
		ribosome_components_by_subunit_header = np.array([
			'Variants', 'Avg Total Counts (Ribosomal 30s)',
			'Avg Limiting Protein Counts (Ribosomal 30s)',
			'Avg 16s rRNA Counts (Ribosomal 30s)'])
		for id in s30_protein_ids:
			ribosome_components_by_subunit = np.vstack((ribosome_components_by_subunit,
				avg_s30_protein_counts['s30_' + id]))
			ribosome_components_by_subunit_header = np.append(
				ribosome_components_by_subunit_header,
				"Avg Monomer Counts (Ribosomal 30s): " + id)

		ribosome_components_by_subunit = np.vstack((
			ribosome_components_by_subunit,
			avg_s50_total_counts,
			avg_s50_limiting_protein_counts,
			avg_s50_23s_rRNA_total_counts, avg_s50_5s_rRNA_total_counts))
		ribosome_components_by_subunit_header = np.append(
			ribosome_components_by_subunit_header,
			('Avg Total Counts (Ribosomal 50s)',
			'Avg Limiting Protein Counts (Ribosomal 50s)',
			'Avg 23s rRNA Counts (Ribosomal 50s)',
			'Avg 5s rRNA Counts (Ribosomal 50s)'))
		for id in s50_protein_ids:
			ribosome_components_by_subunit = np.vstack((ribosome_components_by_subunit,
				avg_s50_protein_counts['s50_' + id]))
			ribosome_components_by_subunit_header = np.append(
				ribosome_components_by_subunit_header,
				"Avg Monomer Counts (Ribosomal 50s): " + id)
		all_ribosome_components_by_subunit = np.vstack((
			ribosome_components_by_subunit_header.T, ribosome_components_by_subunit.T))
		np.savetxt(os.path.join(plotOutDir, plotOutFileName + "_ribosome_components_by_subunit.csv"),
				   all_ribosome_components_by_subunit, delimiter=",", fmt="%s")

		# Save variant index, doubling time, and new gene protein counts
		# to file
		all_avg_monomer_counts = np.vstack((
			variants, doubling_times,
			avg_ng_monomer, avg_ng_monomer_proteome_mass_fraction,
			avg_active_ribosome_counts, avg_inactive_ribosome_counts,
			avg_total_ribosome_counts, avg_total_ribosome_conc,
			avg_active_rnap_counts, avg_inactive_rnap_counts,
			avg_total_rnap_counts, avg_total_rnap_conc,
			avg_dry_mass, avg_cell_mass, avg_mRNA_mass,
			avg_rRNA_mass, avg_tRNA_mass, avg_protein_mass,
			avg_dna_mass, avg_water_mass, avg_small_molecule_mass,
			avg_membrane_mass, percent_completion,
		)).T
		header = np.array([
			'Variant Index', 'Doubling Time (min)',
			'Avg New Gene Protein Counts',
			'Average New Gene Proteome Mass Fraction',
			'Avg Active Ribosome Counts',
			'Avg Inactive Ribosome Counts', 'Avg Total Ribosome Counts',
			'Avg Total Ribosome Concentration (uM)',
			'Avg Active RNA Polymerase Counts',
			'Avg Inactive RNA Polymerase Counts',
		  	'Avg Total RNA Polymerase Counts',
			'Avg Total RNA Polymerase Concentration (uM)',
			'Avg Dry Mass (fg)', 'Avg Cell Mass (fg)',
			'Avg mRNA Mass (fg)', 'Avg rRNA Mass (fg)',
			'Avg tRNA Mass (fg)', 'Avg Protein Mass (fg)',
			'Avg DNA Mass (fg)', 'Avg Water Mass (fg)',
			'Avg Small Molecule Mass (fg)', 'Avg Membrane Mass (fg)',
			'Percent Completion (fraction of seeds that reached all gens)',
		])
		all_avg_monomer_counts = np.vstack((header, all_avg_monomer_counts))
		np.savetxt(os.path.join(plotOutDir, plotOutFileName + ".csv"),
			all_avg_monomer_counts, delimiter=",", fmt="%s")

		plt.figure()
		plt.xlabel("New Gene Protein Counts")
		plt.ylabel("Doubling Time")

		plt.scatter(
			avg_ng_monomer[plot_variant_mask],
			doubling_times[plot_variant_mask])

		plt.tight_layout()
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close('all')

		if INTERACTIVE_PLOT:
			import plotly.express as px

			fig = px.scatter(
				x=avg_ng_monomer[plot_variant_mask],
				y=doubling_times[plot_variant_mask],
				hover_name=np.array(variants)[plot_variant_mask],
				labels={'x': 'New Gene Protein Counts', 'y': 'Doubling Time'},
				hover_data={
					'Variant Index': np.array(variants)[plot_variant_mask]})

			fig.write_html(os.path.join(
				plotOutDir, plotOutFileName + ".html"))
			# fig.show()

if __name__ == "__main__":
	Plot().cli()