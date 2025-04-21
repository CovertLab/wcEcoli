"""
Save the average new gene protein counts and doubling time,
This plot is intended to be run on simulations where the new gene
option was enabled.
"""

import numpy as np
from matplotlib import pyplot as plt

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
		mRNA_monomer_id_dict = dict(zip(monomer_sim_data['cistron_id'],
										monomer_sim_data['id']))
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

				uniqueMoleculeCounts = TableReader(
					os.path.join(simOutDir, "UniqueMoleculeCounts"))
				ribosome_index = uniqueMoleculeCounts.readAttribute(
					"uniqueMoleculeIds").index('active_ribosome')

				uniqueMoleculeCounts = TableReader(
					os.path.join(simOutDir, "UniqueMoleculeCounts"))
				active_rnap_index = uniqueMoleculeCounts.readAttribute(
					"uniqueMoleculeIds").index('active_RNAP')

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

		# Save variant index, doubling time, and new gene protein counts
		# to file
		all_avg_monomer_counts = np.vstack((variants, doubling_times,
			avg_ng_monomer, avg_active_ribosome_counts, avg_inactive_ribosome_counts,
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
			'Avg New Gene Protein Counts', 'Avg Active Ribosome Counts',
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
			fig.show()

if __name__ == "__main__":
	Plot().cli()