"""
Save the average RNA counts and mass portions for each variant index as a
column in a csv file.
"""

import os
import pickle

import numpy as np

from models.ecoli.analysis import variantAnalysisPlot
from wholecell.analysis.analysis_tools import read_stacked_columns, read_stacked_bulk_molecules
from wholecell.io.tablereader import TableReader
from wholecell.utils import units

# Remove first N gens from plot
IGNORE_FIRST_N_GENS = 16

class Plot(variantAnalysisPlot.VariantAnalysisPlot):
	def do_plot(self, inputDir, plotOutDir, plotOutFileName, simDataFile,
				validationDataFile, metadata):

		with open(simDataFile, 'rb') as f:
			sim_data = pickle.load(f)
		with open(validationDataFile, 'rb') as f:
			validation_data = pickle.load(f)

		n_total_gens = self.ap.n_generation

		selected_variant_indexes = self.ap.get_variants()

		n_variants = len(selected_variant_indexes)

		# Initialize everything
		all_cells = self.ap.get_cells(
			variant=[selected_variant_indexes[0]],
			generation=np.arange(IGNORE_FIRST_N_GENS, n_total_gens),
			only_successful=True)

		# Load tables and attributes for RNA cistrons
		RNA_reader = TableReader(
			os.path.join(all_cells[0], 'simOut', 'RNACounts'))
		RNA_cistron_ids = RNA_reader.readAttribute('mRNA_cistron_ids')

		# Load tables and attributes for tRNAs and rRNAs
		unique_molecule_counts_reader = TableReader(
			os.path.join(all_cells[0], 'simOut', 'UniqueMoleculeCounts'))

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

		# Determine which RNAs correspond to essential genes
		cistron_data = sim_data.process.transcription.cistron_data
		cistron_id_to_gene_id = {
			cistron['id']: cistron['gene_id'] for cistron in cistron_data}
		gene_ids = [
			cistron_id_to_gene_id[cistron_id] for cistron_id in RNA_cistron_ids]

		essential_genes = validation_data.essential_genes.essential_genes
		is_essential_RNA = np.array(
			["essential_RNA" if gene_id in essential_genes else "non_essential_RNA" for gene_id in gene_ids])

		is_essential_RNA = np.reshape(is_essential_RNA, (len(is_essential_RNA), 1))
		RNA_cistron_ids = np.reshape(RNA_cistron_ids, (len(RNA_cistron_ids), 1))

		n_mRNA = len(RNA_cistron_ids)
		n_tRNA = len(tRNA_cistron_ids)
		n_rRNA = len(rRNA_cistron_ids)

		# Initialize a pandas dataframe to store the data
		all_avg_RNA_counts = np.zeros((n_mRNA + n_tRNA + n_rRNA, n_variants))
		all_avg_RNA_counts_portion = np.zeros((n_mRNA + n_tRNA + n_rRNA, n_variants))
		all_avg_RNA_transcriptome_mass = np.zeros((n_mRNA + n_tRNA + n_rRNA, n_variants))
		all_avg_RNA_dry_cell_mass = np.zeros((n_mRNA + n_tRNA + n_rRNA, n_variants))

		# Loop through variant indexes
		for i, variant_index in enumerate(selected_variant_indexes):
			# Get all cells (within the generation range) of this variant index
			all_cells = self.ap.get_cells(
				variant=[variant_index],
				generation=np.arange(IGNORE_FIRST_N_GENS, n_total_gens),
				only_successful=True)

			if len(all_cells) == 0:
				continue

			# Read columns
			# remove_first=True because countsToMolar is 0 at first time step
			RNA_cistron_counts = read_stacked_columns(
				all_cells, 'RNACounts', 'mRNA_cistron_counts',
				remove_first=True, ignore_exception=True)

			dry_masses = read_stacked_columns(
				all_cells, 'Mass', 'dryMass',
				remove_first=True, ignore_exception=True)

			(uncharged_tRNA_counts, charged_tRNA_counts, rRNA_counts, ribosomal_subunit_counts) = read_stacked_bulk_molecules(
				all_cells,
				(uncharged_tRNA_ids, charged_tRNA_ids, rRNA_ids, ribosomal_subunit_ids),
				remove_first=True, ignore_exception=True)
			full_ribosome_counts = read_stacked_columns(
				all_cells, 'UniqueMoleculeCounts', 'uniqueMoleculeCounts',
				remove_first=True, ignore_exception=True)[:, ribosome_index]

			# Add up to total counts of tRNAs and rRNAs
			tRNA_counts = uncharged_tRNA_counts + charged_tRNA_counts
			rRNA_counts[:, 0] += ribosomal_subunit_counts[:, 0]
			rRNA_counts[:, 1:] += ribosomal_subunit_counts[:, 1:]
			rRNA_counts += full_ribosome_counts[:, None]

			# Concatenate arrays

			rna_ids = np.concatenate((
				RNA_cistron_ids.squeeze(), np.array(tRNA_cistron_ids),
				np.array(rRNA_cistron_ids)))
			rna_counts = np.hstack((RNA_cistron_counts, tRNA_counts, rRNA_counts))
			tRNA_labels = np.array(['tRNA'] * n_tRNA)
			rRNA_labels = np.array(['rRNA'] * n_rRNA)
			is_essential = np.concatenate((is_essential_RNA.squeeze(), tRNA_labels, rRNA_labels))

			cistron_id_to_mw = {
				cistron_id: cistron_mw for (cistron_id, cistron_mw)
				in zip(
					sim_data.process.transcription.cistron_data['id'],
					sim_data.process.transcription.cistron_data['mw'].asNumber(
						units.fg / units.count))
				}
			rna_mw = np.array(
				[cistron_id_to_mw[cistron_id] for cistron_id in rna_ids])

			rna_counts_avg = rna_counts.mean(axis=0)
			rna_counts_relative_to_total_rna_counts = rna_counts_avg / rna_counts_avg.sum()
			rna_masses_avg = rna_counts_avg * rna_mw
			rna_masses_relative_to_total_rna_mass = rna_masses_avg / rna_masses_avg.sum()
			rna_masses_relative_to_total_dcw = rna_masses_avg / dry_masses.mean()

			all_avg_RNA_counts[:, i] = np.copy(rna_counts_avg)
			all_avg_RNA_counts_portion[:, i] = np.copy(rna_counts_relative_to_total_rna_counts)
			all_avg_RNA_transcriptome_mass[:, i] = np.copy(rna_masses_relative_to_total_rna_mass)
			all_avg_RNA_dry_cell_mass[:, i] = np.copy(rna_masses_relative_to_total_dcw)

		# Save data
		rna_ids = np.reshape(rna_ids, (len(rna_ids), 1))
		is_essential = np.reshape(is_essential, (len(is_essential), 1))

		all_labeled_avg_RNA_counts = np.hstack((
			rna_ids, is_essential, all_avg_RNA_counts))
		header = np.array(
			['RNA_cistron_ids'] + ['is_essential'] + selected_variant_indexes)
		all_labeled_avg_RNA_cistron_counts = np.vstack((
			header, all_labeled_avg_RNA_counts))
		np.savetxt(
			os.path.join(plotOutDir, 'all_labeled_avg_RNA_counts.csv'),
			all_labeled_avg_RNA_cistron_counts, delimiter=',', fmt='%s')

		all_labeled_avg_RNA_counts_portion = np.hstack((
			rna_ids, is_essential, all_avg_RNA_counts_portion))
		header = np.array(
			['RNA_cistron_ids'] + ['is_essential'] + selected_variant_indexes)
		all_labeled_avg_RNA_counts_portion = np.vstack((
			header, all_labeled_avg_RNA_counts_portion))
		np.savetxt(
			os.path.join(plotOutDir, 'all_labeled_avg_RNA_counts_portion.csv'),
			all_labeled_avg_RNA_counts_portion, delimiter=',', fmt='%s')

		all_labeled_avg_RNA_transcriptome_mass = np.hstack((
			rna_ids, is_essential, all_avg_RNA_transcriptome_mass))
		header = np.array(
			['RNA_cistron_ids'] + ['is_essential'] + selected_variant_indexes)
		all_labeled_avg_RNA_transcriptome_mass = np.vstack((
			header, all_labeled_avg_RNA_transcriptome_mass))
		np.savetxt(
			os.path.join(plotOutDir, 'all_labeled_avg_RNA_transcriptome_mass.csv'),
			all_labeled_avg_RNA_transcriptome_mass, delimiter=',', fmt='%s')

		all_labeled_avg_RNA_dry_cell_mass = np.hstack((
			rna_ids, is_essential, all_avg_RNA_dry_cell_mass))
		header = np.array(
			['RNA_cistron_ids'] + ['is_essential'] + selected_variant_indexes)
		all_labeled_avg_RNA_dry_cell_mass = np.vstack((
			header, all_labeled_avg_RNA_dry_cell_mass))
		np.savetxt(
			os.path.join(plotOutDir, 'all_labeled_avg_RNA_dry_cell_mass.csv'),
			all_labeled_avg_RNA_dry_cell_mass, delimiter=',', fmt='%s')

if __name__ == "__main__":
	Plot().cli()
