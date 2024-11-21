"""
Generates tables of data to share with EcoCyc for display on the "modeling" tab.

TODO:
	other values
		weighted average for counts (time step weighted and cell cycle progress weighted)
		max/min
"""

import csv
import json
import os
import pickle

import numpy as np
from scipy.stats import pearsonr

from models.ecoli.analysis import cohortAnalysisPlot
from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from models.ecoli.processes.metabolism import (
	COUNTS_UNITS, VOLUME_UNITS, TIME_UNITS, MASS_UNITS)
from wholecell.analysis.analysis_tools import (read_stacked_bulk_molecules,
	read_stacked_columns)
from wholecell.io.tablereader import TableReader
from wholecell.utils import units

IGNORE_FIRST_N_GENS = 2

# TODO (ggsun): Add this to sim_data somewhere?
# Maps media names used in model to IDs used in EcoCyc
MEDIA_NAME_TO_ID = {
	'minimal': 'MIX0-57',
	'minimal_minus_oxygen': 'MIX0-57-ANAEROBIC',
	'minimal_plus_amino_acids': 'MIX0-850',
	'minimal_acetate': 'MIX0-58',
	'minimal_succinate': 'MIX0-844',
	}

# We do not use this in our model, hardcoding it here to route active ribosome
# counts data to the correct EcoCyc page
ECOCYC_FULL_RIBOSOME_ID = 'CPLX0-3964'

def save_file(out_dir, filename, columns, values):
	output_file = os.path.join(out_dir, filename)
	print(f'Saving data to {output_file}')
	with open(output_file, 'w') as f:
		writer = csv.writer(f, delimiter='\t')

		# Header for columns
		writer.writerow(['# Column descriptions:'])

		for col, desc in columns.items():
			writer.writerow([f'# {col}', desc])
		writer.writerow(list(columns.keys()))

		# Data rows
		for i in np.arange(len(values[0])):
			writer.writerow([v[i] for v in values])


class Plot(cohortAnalysisPlot.CohortAnalysisPlot):
	def do_plot(self, variantDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		with open(simDataFile, 'rb') as f:
			sim_data = pickle.load(f)
		with open(validationDataFile, 'rb') as f:
			validation_data = pickle.load(f)

		media_name = sim_data.conditions[sim_data.condition]['nutrients']
		media_id = MEDIA_NAME_TO_ID.get(media_name, media_name)

		ap = AnalysisPaths(variantDir, cohort_plot=True)

		# Ignore first N generations
		cell_paths = ap.get_cells(
			generation=np.arange(IGNORE_FIRST_N_GENS, ap.n_generation),
			only_successful=True)

		if len(cell_paths) == 0:
			print('Skipping analysis -- not enough simulations run.')
			return

		# Load tables and attributes for mRNAs
		RNA_synth_prob_reader = TableReader(
			os.path.join(cell_paths[0], 'simOut', 'RnaSynthProb'))
		RNA_reader = TableReader(
			os.path.join(cell_paths[0], 'simOut', 'RNACounts'))
		mass_reader = TableReader(
			os.path.join(cell_paths[0], 'simOut', 'Mass'))

		mRNA_ids = RNA_reader.readAttribute('mRNA_cistron_ids')
		mass_unit =	mass_reader.readAttribute('cellDry_units')
		assert mass_unit == 'fg'
		gene_ids_rna_synth_prob = RNA_synth_prob_reader.readAttribute("gene_ids")

		# Load tables and attributes for tRNAs and rRNAs
		unique_molecule_counts_reader = TableReader(
			os.path.join(cell_paths[0], 'simOut', 'UniqueMoleculeCounts'))

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
		ribosome_index = unique_molecule_counts_reader.readAttribute('uniqueMoleculeIds').index('active_ribosome')

		# # Read columns
		# # remove_first=True because countsToMolar is 0 at first time step
		# gene_copy_numbers = read_stacked_columns(
		# 	cell_paths, 'RnaSynthProb', 'gene_copy_number',
		# 	remove_first=True, ignore_exception=True)
		# mRNA_counts = read_stacked_columns(
		# 	cell_paths, 'RNACounts', 'mRNA_cistron_counts',
		# 	remove_first=True, ignore_exception=True)
		counts_to_molar = read_stacked_columns(
			cell_paths, 'EnzymeKinetics', 'countsToMolar',
			remove_first=True, ignore_exception=True)
		# dry_masses = read_stacked_columns(
		# 	cell_paths, 'Mass', 'dryMass',
		# 	remove_first=True, ignore_exception=True)
		#
		# (uncharged_tRNA_counts, charged_tRNA_counts, rRNA_counts, ribosomal_subunit_counts) = read_stacked_bulk_molecules(
		# 	cell_paths,
		# 	(uncharged_tRNA_ids, charged_tRNA_ids, rRNA_ids, ribosomal_subunit_ids),
		# 	remove_first=True, ignore_exception=True)
		# full_ribosome_counts = read_stacked_columns(
		# 	cell_paths, 'UniqueMoleculeCounts', 'uniqueMoleculeCounts',
		# 	remove_first=True, ignore_exception=True)[:, ribosome_index]
		#
		# # Add up to total counts of tRNAs and rRNAs
		# tRNA_counts = uncharged_tRNA_counts + charged_tRNA_counts
		# rRNA_counts[:, 0] += ribosomal_subunit_counts[:, 0]
		# rRNA_counts[:, 1:] += ribosomal_subunit_counts[:, 1:]
		# rRNA_counts += full_ribosome_counts[:, None]
		#
		# # Concatenate arrays
		# n_mRNA = len(mRNA_ids)
		# n_tRNA = len(tRNA_cistron_ids)
		# n_rRNA = len(rRNA_cistron_ids)
		# rna_ids = np.concatenate((mRNA_ids, tRNA_cistron_ids, rRNA_cistron_ids))
		# rna_counts = np.hstack((mRNA_counts, tRNA_counts, rRNA_counts))
		#
		# cistron_id_to_mw = {
		# 	cistron_id: cistron_mw for (cistron_id, cistron_mw)
		# 	in zip(
		# 		sim_data.process.transcription.cistron_data['id'],
		# 		sim_data.process.transcription.cistron_data['mw'].asNumber(
		# 			units.fg / units.count))
		# 	}
		# rna_mw = np.array(
		# 	[cistron_id_to_mw[cistron_id] for cistron_id in rna_ids])
		#
		# # Calculate derived RNA values
		# def normalize_within_each_type(values, n_mRNA, n_tRNA, n_rRNA):
		# 	assert len(values) == n_mRNA + n_tRNA + n_rRNA
		# 	normalized_values = np.zeros_like(values)
		#
		# 	normalized_values[:n_mRNA] = values[:n_mRNA] / values[:n_mRNA].sum()
		# 	normalized_values[n_mRNA:(n_mRNA + n_tRNA)] = values[n_mRNA:(n_mRNA + n_tRNA)] / values[n_mRNA:(n_mRNA + n_tRNA)].sum()
		# 	normalized_values[-n_rRNA:] = values[-n_rRNA:] / values[-n_rRNA:].sum()
		#
		# 	return normalized_values
		#
		# gene_copy_numbers_avg = gene_copy_numbers.mean(axis=0)
		# gene_copy_numbers_std = gene_copy_numbers.std(axis=0)
		# rna_counts_avg = rna_counts.mean(axis=0)
		# rna_counts_std = rna_counts.std(axis=0)
		# rna_conc = rna_counts * counts_to_molar
		# rna_conc_avg = rna_conc.mean(axis=0)
		# rna_conc_std = rna_conc.std(axis=0)
		# rna_counts_relative_to_total_rna_counts = rna_counts_avg / rna_counts_avg.sum()
		# rna_counts_relative_to_total_rna_type_counts = normalize_within_each_type(
		# 	rna_counts_avg, n_mRNA, n_tRNA, n_rRNA)
		# rna_masses_avg = rna_counts_avg * rna_mw
		# rna_masses_relative_to_total_rna_mass = rna_masses_avg / rna_masses_avg.sum()
		# rna_masses_relative_to_total_rna_type_mass = normalize_within_each_type(
		# 	rna_masses_avg, n_mRNA, n_tRNA, n_rRNA)
		# rna_masses_relative_to_total_dcw = rna_masses_avg / dry_masses.mean()
		#
		# # Save RNA data in table
		# cistron_id_to_gene_id = {
		# 	cistron['id']: cistron['gene_id']
		# 	for cistron in sim_data.process.transcription.cistron_data
		# 	}
		# gene_ids = [cistron_id_to_gene_id[x] for x in rna_ids]
		#
		# gene_id_to_index = {
		# 	gene_id: i for i,gene_id in enumerate(gene_ids_rna_synth_prob)
		# 	}
		# reordering_indexes = np.array([
		# 	gene_id_to_index[gene_id] for gene_id in gene_ids])
		# assert np.all(
		# 	np.array(gene_ids_rna_synth_prob)[reordering_indexes] == gene_ids)
		# gene_copy_numbers_avg = gene_copy_numbers_avg[reordering_indexes]
		# gene_copy_numbers_std = gene_copy_numbers_std[reordering_indexes]
		#
		# columns = {
		# 	'id': 'Object ID, according to EcoCyc',
		# 	'gene-copy-number-avg': 'A floating point number',
		# 	'gene-copy-number-std': 'A floating point number',
		# 	'rna-count-avg': 'A floating point number',
		# 	'rna-count-std': 'A floating point number',
		# 	'rna-concentration-avg': 'A floating point number in mM units',
		# 	'rna-concentration-std': 'A floating point number in mM units',
		# 	'relative-rna-count-to-total-rna-counts': 'A floating point number',
		# 	'relative-rna-count-to-total-rna-type-counts': 'A floating point number',
		# 	'relative-rna-mass-to-total-rna-mass': 'A floating point number',
		# 	'relative-rna-mass-to-total-rna-type-mass': 'A floating point number',
		# 	'relative-rna-mass-to-total-cell-dry-mass': 'A floating point number',
		# 	}
		# values = [
		# 	gene_ids, gene_copy_numbers_avg, gene_copy_numbers_std,
		# 	rna_counts_avg, rna_counts_std, rna_conc_avg,
		# 	rna_conc_std, rna_counts_relative_to_total_rna_counts,
		# 	rna_counts_relative_to_total_rna_type_counts,
		# 	rna_masses_relative_to_total_rna_mass,
		# 	rna_masses_relative_to_total_rna_type_mass,
		# 	rna_masses_relative_to_total_dcw,
		# 	]
		#
		# save_file(
		# 	plotOutDir, f'wcm_rnas_{media_id}.tsv', columns, values)
		#

		#
		# # Load tables and attributes for proteins
		monomer_reader = TableReader(
			os.path.join(cell_paths[0], 'simOut', 'MonomerCounts'))
		monomer_ids = monomer_reader.readAttribute('monomerIds')
		monomer_mw = sim_data.getter.get_masses(monomer_ids).asNumber(units.fg / units.count)
		#
		# Read columns
		# remove_first=True because countsToMolar is 0 at first time step
		monomer_counts = read_stacked_columns(
			cell_paths, 'MonomerCounts', 'monomerCounts',
			remove_first=True, ignore_exception=True)
		#
		# # Calculate derived protein values
		# monomer_counts_avg = monomer_counts.mean(axis=0)
		# monomer_counts_std = monomer_counts.std(axis=0)
		# monomer_conc = monomer_counts * counts_to_molar
		# monomer_conc_avg = monomer_conc.mean(axis=0)
		# monomer_conc_std = monomer_conc.std(axis=0)
		# monomer_counts_relative_to_total_monomer_counts = monomer_counts_avg / monomer_counts_avg.sum()
		# monomer_mass_avg = monomer_counts_avg * monomer_mw
		# monomer_mass_relative_to_total_monomer_mass = monomer_mass_avg / monomer_mass_avg.sum()
		# monomer_mass_relative_to_total_dcw = monomer_mass_avg / dry_masses.mean()
		#
		# # Save monomer data in table
		# monomer_ecocyc_ids = [monomer[:-3] for monomer in monomer_ids]  # strip [*]
		#
		# columns = {
		# 	'id': 'Object ID, according to EcoCyc',
		# 	'protein-count-avg': 'A floating point number',
		# 	'protein-count-std': 'A floating point number',
		# 	'protein-concentration-avg': 'A floating point number in mM units',
		# 	'protein-concentration-std': 'A floating point number in mM units',
		# 	'relative-protein-count-to-protein-rna-counts': 'A floating point number',
		# 	'relative-protein-mass-to-total-protein-mass': 'A floating point number',
		# 	'relative-protein-mass-to-total-cell-dry-mass': 'A floating point number',
		# 	}
		# values = [
		# 	monomer_ecocyc_ids, monomer_counts_avg, monomer_counts_std,
		# 	monomer_conc_avg, monomer_conc_std,
		# 	monomer_counts_relative_to_total_monomer_counts,
		# 	monomer_mass_relative_to_total_monomer_mass,
		# 	monomer_mass_relative_to_total_dcw,
		# 	]
		#
		# # Add validation data if sims used minimal glucose media
		# if media_name == 'minimal':
		# 	protein_id_to_schmidt_counts = {
		# 		item[0]: item[1] for item in validation_data.protein.schmidt2015Data
		# 		}
		# 	protein_counts_val = np.array([
		# 		protein_id_to_schmidt_counts.get(protein_id, np.nan) for protein_id in monomer_ids
		# 		])
		#
		# 	columns['validation-count'] = 'A floating point number'
		# 	values.append(protein_counts_val)
		#
		# 	protein_val_exists = np.logical_not(np.isnan(protein_counts_val))
		# 	r, _ = pearsonr(
		# 		monomer_counts_avg[protein_val_exists],
		# 		protein_counts_val[protein_val_exists])
		#
		# 	ecocyc_metadata['protein_validation_r_squared'] = r ** 2
		#
		# save_file(
		# 	plotOutDir, f'wcm_monomers_{media_id}.tsv', columns, values)
		#
		# # Load attributes for complexes
		# complex_ids = sim_data.process.complexation.ids_complexes
		# complex_ecocyc_ids = [complex_id[:-3] for complex_id in complex_ids]  # strip [*]
		# complex_mw = sim_data.getter.get_masses(complex_ids).asNumber(units.fg / units.count)
		#
		# # Read columns
		# # remove_first=True because countsToMolar is 0 at first time step
		# (complex_counts, ) = read_stacked_bulk_molecules(
		# 	cell_paths, (complex_ids, ), remove_first=True,
		# 	ignore_exception=True)
		#
		# # Add active RNA polymerase counts to inactive RNA polymerase counts
		# active_RNAP_index = unique_molecule_counts_reader.readAttribute(
		# 	"uniqueMoleculeIds").index('active_RNAP')
		# active_RNAP_counts = read_stacked_columns(
		# 	cell_paths, 'UniqueMoleculeCounts',
		# 	'uniqueMoleculeCounts', remove_first=True,
		# 	ignore_exception=True)[:, active_RNAP_index]
		# inactive_RNAP_id = sim_data.molecule_ids.full_RNAP
		# complex_table_RNAP_index = complex_ids.index(inactive_RNAP_id)
		# complex_counts[:, complex_table_RNAP_index] += active_RNAP_counts
		#
		# # Count active ribosomes as full ribosome complexes
		# complex_ids.append(ECOCYC_FULL_RIBOSOME_ID)
		# complex_ecocyc_ids.append(ECOCYC_FULL_RIBOSOME_ID)
		# full_ribosome_mw = (sim_data.getter.get_mass(
		# 	sim_data.molecule_ids.s50_full_complex) + sim_data.getter.get_mass(
		# 	sim_data.molecule_ids.s30_full_complex)).asNumber(units.fg / units.count)
		# complex_mw = np.append(complex_mw, full_ribosome_mw)
		# complex_counts = np.hstack((complex_counts, full_ribosome_counts[:, None]))
		#
		# # Calculate derived protein values
		# complex_counts_avg = complex_counts.mean(axis=0)
		# complex_counts_std = complex_counts.std(axis=0)
		# complex_conc = complex_counts * counts_to_molar
		# complex_conc_avg = complex_conc.mean(axis=0)
		# complex_conc_std = complex_conc.std(axis=0)
		# complex_mass_avg = complex_counts_avg * complex_mw
		# complex_mass_relative_to_total_protein_mass = complex_mass_avg / monomer_mass_avg.sum()
		# complex_mass_relative_to_total_dcw = complex_mass_avg / dry_masses.mean()
		#
		# # Save complex data in table
		# columns = {
		# 	'id': 'Object ID, according to EcoCyc',
		# 	'complex-count-avg': 'A floating point number',
		# 	'complex-count-std': 'A floating point number',
		# 	'complex-concentration-avg': 'A floating point number in mM units',
		# 	'complex-concentration-std': 'A floating point number in mM units',
		# 	'relative-complex-mass-to-total-protein-mass': 'A floating point number',
		# 	'relative-complex-mass-to-total-cell-dry-mass': 'A floating point number',
		# 	}
		# values = [
		# 	complex_ecocyc_ids, complex_counts_avg, complex_counts_std,
		# 	complex_conc_avg, complex_conc_std,
		# 	complex_mass_relative_to_total_protein_mass,
		# 	complex_mass_relative_to_total_dcw,
		# 	]
		#
		# save_file(
		# 	plotOutDir, f'wcm_complexes_{media_id}.tsv', columns, values)
		#
		# # Load attributes for metabolic fluxes
		# cell_density = sim_data.constants.cell_density
		# reaction_ids = sim_data.process.metabolism.base_reaction_ids
		#
		# # Read columns
		# cell_mass = read_stacked_columns(
		# 	cell_paths, 'Mass', 'cellMass', ignore_exception=True)
		# dry_mass = read_stacked_columns(
		# 	cell_paths, 'Mass', 'dryMass', ignore_exception=True)
		# conversion_coeffs = (
		# 	dry_mass / cell_mass
		# 	* cell_density.asNumber(MASS_UNITS / VOLUME_UNITS)
		# )
		#
		# # Calculate flux in units of mmol/g DCW/h
		# fluxes = (
		# 	(COUNTS_UNITS / MASS_UNITS / TIME_UNITS)
		# 	* (read_stacked_columns(cell_paths, 'FBAResults', 'base_reaction_fluxes', ignore_exception=True) / conversion_coeffs)
		# 	).asNumber(units.mmol / units.g / units.h)
		#
		# # Calculate derived flux values
		# fluxes_avg = fluxes.mean(axis=0)
		# fluxes_std = fluxes.std(axis=0)
		#
		# columns = {
		# 	'id': 'Object ID, according to EcoCyc',
		# 	'flux-avg': 'A floating point number in mmol/g DCW/h units',
		# 	'flux-std': 'A floating point number in mmol/g DCW/h units',
		# 	}
		# values = [
		# 	reaction_ids, fluxes_avg, fluxes_std,
		# 	]
		#
		# save_file(
		# 	plotOutDir, f'wcm_metabolic_reactions_{media_id}.tsv', columns, values)

		# Cell level properties for WCM overview page
		# Load tables and attributes
		ribosomal_proteins = sim_data.molecule_groups.ribosomal_proteins
		mw_ribosomal_proteins = (
			sim_data.getter.get_masses(ribosomal_proteins).asNumber(units.fg / units.count))
		ribosomal_protein_indexes = [
			monomer_ids.index(ribosomal_protein) for ribosomal_protein in ribosomal_proteins]
		# complexation = sim_data.process.complexation
		# equilibrium = sim_data.process.equilibrium
		# two_component_system = sim_data.process.two_component_system
		# enzyme_ids = sim_data.process.metabolism.catalyst_ids
		#
		# enzyme_proteins = set()
		# nonclassified_enzymes = set()
		# for enzyme in enzyme_ids:
		# 	subunit_ids = []
		# 	if enzyme in complexation.complex_names:
		# 		subunit_ids = complexation.get_monomers(enzyme)['subunitIds']
		# 	elif enzyme in equilibrium.ids_complexes:
		# 		subunit_ids = equilibrium.get_monomers(enzyme)['subunitIds']
		# 	elif enzyme in two_component_system.complex_to_monomer.keys():
		# 		subunit_ids = two_component_system.get_monomers(enzyme)['subunitIds']
		# 	elif enzyme in monomer_ids:
		# 		continue
		# 	else:
		# 		nonclassified_enzymes.add(enzyme)
		# 		if enzyme == "CPLX0-7761[c]" or enzyme == "G6260-MONOMER[m]" or enzyme == "MONOMER-48[c]": ### TODO: Investigate why this is happening
		# 			continue
		# 		raise ValueError(f'Enzyme {enzyme} not found in any complex generator.')
		#
		# 	for subunit in subunit_ids:
		# 		enzyme_proteins.add(subunit)
		#
		# enzyme_proteins = sorted(enzyme_proteins)
		#
		# ### TODO: DELETE THIS
		# enzyme_proteins_old = sorted({ # Map enzyme complexes to monomers
		# 	subunit
		# 	for enzyme in enzyme_ids
		# 	for subunit in complexation.get_monomers(enzyme)['subunitIds']
		# })
		#
		# print(nonclassified_enzymes)
		#
		# nonclassified_monomers = set()
		# for protein in enzyme_proteins:
		# 	if protein not in monomer_ids:
		# 		nonclassified_monomers.add(protein)
		#
		# import ipdb
		# ipdb.set_trace()
		#
		# enzyme_proteins_indexes = [
		# 	monomer_ids.index(enzyme_protein) for enzyme_protein in enzyme_proteins]
		# mw_enzyme_proteins = (
		# 	sim_data.getter.get_masses(enzyme_proteins).asNumber(units.fg / units.count))
		#
		#
		#
		# import ipdb
		# ipdb.set_trace()

		# Read columns
		cell_mass = read_stacked_columns(
			cell_paths, 'Mass', 'cellMass', fun=lambda x: x[0],
			ignore_exception=True)
		dry_mass = read_stacked_columns(
			cell_paths, 'Mass', 'dryMass', fun=lambda x: x[0],
			ignore_exception=True)
		protein_mass = read_stacked_columns(
			cell_paths, 'Mass', 'proteinMass', fun=lambda x: x[0],
			ignore_exception=True)
		initial_rna_mass = read_stacked_columns(
			cell_paths, 'Mass', 'rnaMass', fun=lambda x: x[0],
			ignore_exception=True)
		doubling_time = read_stacked_columns(
			cell_paths, 'Main', 'time', fun=lambda x: (x[-1] - x[0]) / 60.,
			ignore_exception=True)
		time_step_length = read_stacked_columns(
			cell_paths, 'Main', 'timeStepSec', ignore_exception=True)
		cell_volume = read_stacked_columns(
			cell_paths, 'Mass', 'cellVolume', fun=lambda x: x[0],
			ignore_exception=True)
		membrane_mass = read_stacked_columns(
			cell_paths, 'Mass', 'membrane_mass', fun=lambda x: x[0],
			ignore_exception=True)
		dna_mass = read_stacked_columns(
			cell_paths, 'Mass', 'dnaMass', fun=lambda x: x[0],
			ignore_exception=True)
		small_molecule_mass = read_stacked_columns(
			cell_paths, 'Mass', 'smallMoleculeMass', fun=lambda x: x[0],
			ignore_exception=True)

		rna_mass = read_stacked_columns(
			cell_paths, 'Mass', 'rnaMass', ignore_exception=True)
		mrna_mass = read_stacked_columns(
			cell_paths, 'Mass', 'mRnaMass', ignore_exception=True)
		rrna_mass = read_stacked_columns(
			cell_paths, 'Mass', 'rRnaMass', ignore_exception=True)
		trna_mass = read_stacked_columns(
			cell_paths, 'Mass', 'tRnaMass', ignore_exception=True)

		protein_mass_over_time = read_stacked_columns(
			cell_paths, 'Mass', 'proteinMass',
			remove_first=True, ignore_exception=True).squeeze()
		ribosomal_monomer_counts = monomer_counts[:, ribosomal_protein_indexes]
		ribosomal_monomer_mass = ribosomal_monomer_counts * mw_ribosomal_proteins
		total_ribosomal_monomer_mass = ribosomal_monomer_mass.sum(axis = 1)
		total_ribosomal_monomer_mass_fraction = total_ribosomal_monomer_mass / protein_mass_over_time

		# enzyme_monomer_counts = monomer_counts[:, enzyme_proteins_indexes]

		# Calculate derived cell values
		cell_mass_avg = cell_mass.mean(axis = 0)
		cell_mass_std = cell_mass.std(axis = 0)
		dry_mass_avg = dry_mass.mean(axis = 0)
		dry_mass_std = dry_mass.std(axis = 0)
		protein_mass_avg = protein_mass.mean(axis = 0)
		protein_mass_std = protein_mass.std(axis = 0)
		initial_rna_mass_avg = initial_rna_mass.mean(axis = 0)
		initial_rna_mass_std = initial_rna_mass.std(axis = 0)
		doubling_time_avg = doubling_time.mean(axis = 0)
		doubling_time_std = doubling_time.std(axis = 0)
		time_step_length_avg = time_step_length.mean(axis = 0)
		time_step_length_std = time_step_length.std(axis = 0)
		cell_volume_avg = cell_volume.mean(axis = 0)
		cell_volume_std = cell_volume.std(axis = 0)
		membrane_mass_avg = membrane_mass.mean(axis = 0)
		membrane_mass_std = membrane_mass.std(axis = 0)
		dna_mass_avg = dna_mass.mean(axis = 0)
		dna_mass_std = dna_mass.std(axis = 0)
		small_molecule_mass_avg = small_molecule_mass.mean(axis = 0)
		small_molecule_mass_std = small_molecule_mass.std(axis = 0)

		mrna_fraction = mrna_mass / rna_mass
		mrna_fraction_avg = mrna_fraction.mean(axis = 0)
		mrna_fraction_std = mrna_fraction.std(axis = 0)

		rrna_fraction = rrna_mass / rna_mass
		rrna_fraction_avg = rrna_fraction.mean(axis = 0)
		rrna_fraction_std = rrna_fraction.std(axis = 0)

		trna_fraction = trna_mass / rna_mass
		trna_fraction_avg = trna_fraction.mean(axis = 0)
		trna_fraction_std = trna_fraction.std(axis = 0)

		ribosomal_protein_mass_fraction_avg = total_ribosomal_monomer_mass_fraction.mean()
		ribosomal_protein_mass_fraction_std = total_ribosomal_monomer_mass_fraction.std()

		# Build dictionary for metadata
		ecocyc_metadata = {
			'git_hash': metadata['git_hash'], ### TODO: ADD A DOC
			'git_hash_doc': 'A string',
			'n_ignored_generations': IGNORE_FIRST_N_GENS,
			'n_ignored_generations_doc': 'An integer',
			'n_total_generations': metadata['total_gens'],
			'n_total_generations_doc': 'An integer',
			'n_seeds': metadata['total_init_sims'],
			'n_seeds_doc': 'An integer',
			'n_cells': len(cell_paths),
			'n_cells_doc': 'An integer',
			'n_timesteps': len(counts_to_molar),
			'n_timesteps_doc': 'An integer',
			'cell_mass_avg': cell_mass_avg.item(),
			'cell_mass_avg_doc': 'A floating point number in fg units',
			'cell_mass_std': cell_mass_std.item(),
			'cell_mass_std_doc': 'A floating point number in fg units',
			'dry_mass_avg': dry_mass_avg.item(),
			'dry_mass_avg_doc': 'A floating point number in fg units',
			'dry_mass_std': dry_mass_std.item(),
			'dry_mass_std_doc': 'A floating point number in fg units',
			'protein_mass_avg': protein_mass_avg.item(),
			'protein_mass_avg_doc': 'A floating point number in fg units',
			'protein_mass_std': protein_mass_std.item(),
			'protein_mass_std_doc': 'A floating point number in fg units',
			'rna_mass_avg': initial_rna_mass_avg.item(),
			'rna_mass_avg_doc': 'A floating point number in fg units',
			'rna_mass_std': initial_rna_mass_std.item(),
			'rna_mass_std_doc': 'A floating point number in fg units',
			'doubling_time_avg': doubling_time_avg.item(),
			'doubling_time_avg_doc': 'A floating point number in min units',
			'doubling_time_std': doubling_time_std.item(),
			'doubling_time_std_doc': 'A floating point number in min units',
			'time_step_length_avg': time_step_length_avg.item(),
			'time_step_length_avg_doc': 'A floating point number in sec units',
			'time_step_length_std': time_step_length_std.item(),
			'time_step_length_std_doc': 'A floating point number in sec units',
			'mrna_mass_to_total_rna_mass_avg': mrna_fraction_avg.item(),
			'mrna_mass_to_total_rna_mass_avg_doc': 'A floating point number',
			'mrna_mass_to_total_rna_mass_std': mrna_fraction_std.item(),
			'mrna_mass_to_total_rna_mass_std_doc': 'A floating point number',
			'rrna_mass_to_total_rna_mass_avg': rrna_fraction_avg.item(),
			'rrna_mass_to_total_rna_mass_avg_doc': 'A floating point number',
			'rrna_mass_to_total_rna_mass_std': rrna_fraction_std.item(),
			'rrna_mass_to_total_rna_mass_std_doc': 'A floating point number',
			'trna_mass_to_total_rna_mass_avg': trna_fraction_avg.item(),
			'trna_mass_to_total_rna_mass_avg_doc': 'A floating point number',
			'trna_mass_to_total_rna_mass_std': trna_fraction_std.item(),
			'trna_mass_to_total_rna_mass_std_doc': 'A floating point number',
			'cell_volume_avg': cell_volume_avg.item(),
			'cell_volume_avg_doc': 'A floating point number in fL units',
			'cell_volume_std': cell_volume_std.item(),
			'cell_volume_std_doc': 'A floating point number in fL units',
			'membrane_mass_avg': membrane_mass_avg.item(),
			'membrane_mass_avg_doc': 'A floating point number in fg units',
			'membrane_mass_std': membrane_mass_std.item(),
			'membrane_mass_std_doc': 'A floating point number in fg units',
			'dna_mass_avg': dna_mass_avg.item(),
			'dna_mass_avg_doc': 'A floating point number in fg units',
			'dna_mass_std': dna_mass_std.item(),
			'dna_mass_std_doc': 'A floating point number in fg units',
			'small_molecule_mass_avg': small_molecule_mass_avg.item(),
			'small_molecule_mass_avg_doc': 'A floating point number in fg units',
			'small_molecule_mass_std': small_molecule_mass_std.item(),
			'small_molecule_mass_std_doc': 'A floating point number in fg units',
			'ribosomal_protein_mass_fraction_avg': ribosomal_protein_mass_fraction_avg.item(),
			'ribosomal_protein_mass_fraction_avg_doc': 'A floating point number',
			'ribosomal_protein_mass_fraction_std': ribosomal_protein_mass_fraction_std.item(),
			'ribosomal_protein_mass_fraction_std_doc': 'A floating point number',
			}

		metadata_file = os.path.join(plotOutDir, f'wcm_metadata_{media_id}.json')
		with open(metadata_file, 'w') as f:
			print(f'Saving data to {metadata_file}')
			json.dump(ecocyc_metadata, f, indent=4)

if __name__ == '__main__':
	Plot().cli()
