"""
Plot mRNA and protein counts for genes across multiple generations
"""

import pickle
import os

from matplotlib import pyplot as plt
# noinspection PyUnresolvedReferences
import numpy as np

from models.ecoli.analysis import multigenAnalysisPlot
from wholecell.analysis.analysis_tools import (exportFigure,
	read_stacked_bulk_molecules, read_stacked_columns)
from wholecell.io.tablereader import TableReader

# Replace with the proteins you would like to visualize here:
interest_proteins = np.array([
	# 'ACRD-MONOMER[i]',
	# 'CYNX-MONOMER[i]',
	# 'B0270-MONOMER[i]',
	# 'G7634-MONOMER[i]',
	#'EG11854-MONOMER[c]',
	#'G6606-MONOMER[c]',
	'MONOMER0-2678[c]',
	'EG10037-MONOMER[c]',
	'PD00519[c]',
])

class Plot(multigenAnalysisPlot.MultigenAnalysisPlot):
	def do_plot(self, seedOutDir, plotOutDir, plotOutFileName, simDataFile,
				validationDataFile, metadata):
		with open(simDataFile, 'rb') as f:
			sim_data = pickle.load(f)
		with open(validationDataFile, 'rb') as f:
			validation_data = pickle.load(f)

		cell_paths = self.ap.get_cells()
		variant = self.ap.get_variants()
		seed = self.ap.get_seeds()
		sim_dir = cell_paths[0]
		simOutDir = os.path.join(sim_dir, 'simOut')

		# Determine new gene ids
		with open(simDataFile, 'rb') as f:
			sim_data = pickle.load(f)
		mRNA_cistron_sim_data = sim_data.process.transcription.cistron_data.struct_array
		monomer_sim_data = sim_data.process.translation.monomer_data.struct_array

		# extract info about the protein(s) from the monomer data:
		monomer_data_idxs = []
		for protein in interest_proteins:
			monomer_idx = np.where(monomer_sim_data['id'] == protein)
			monomer_idx = monomer_idx[0][0]
			monomer_data_idxs.append(monomer_idx)
		ip_monomer_data = monomer_sim_data[monomer_data_idxs]
		ip_monomer_ids = ip_monomer_data['id']
		ip_cistron_ids = ip_monomer_data['cistron_id']

		# extract info about the protein(s) from the mRNA/cistron data:
		mRNA_data_idxs = []
		for cistron in ip_cistron_ids:
			mRNA_idx = np.where(mRNA_cistron_sim_data['id'] == cistron)
			mRNA_idx = mRNA_idx[0][0]
			mRNA_data_idxs.append(mRNA_idx)
		ip_mRNA_data = mRNA_cistron_sim_data[mRNA_data_idxs]
		ip_gene_ids = np.array(ip_mRNA_data['gene_id'])
		cistron_ids = ip_cistron_ids

		cistron_monomer_id_dict = dict(zip(monomer_sim_data['cistron_id'],
										monomer_sim_data['id']))
		cistron_monomer_ids = [cistron_monomer_id_dict.get(mRNA_id)
								for mRNA_id in cistron_ids]

		# Extract mRNA indexes for each gene of interest
		mRNA_counts_reader = TableReader(os.path.join(simOutDir,
													  'RNACounts'))
		mRNA_idx_dict = {rna: i for i, rna in enumerate(
			mRNA_counts_reader.readAttribute('mRNA_cistron_ids'))}
		new_gene_mRNA_indexes = [mRNA_idx_dict.get(mRNA_id) for mRNA_id in
								 cistron_ids]

		# Load data
		time = read_stacked_columns(cell_paths, 'Main', 'time')
		(ip_monomer_counts,) = read_stacked_bulk_molecules(
			cell_paths, cistron_monomer_ids)
		all_mRNA_stacked_counts = read_stacked_columns(
			cell_paths, 'RNACounts', 'mRNA_cistron_counts')
		ip_mRNA_counts = all_mRNA_stacked_counts[:, new_gene_mRNA_indexes]

		# Plotting
		plt.figure(figsize = (8.5, 11))

		# section out the generations:
		# Get doubling times from cells with this variant index
		dt = read_stacked_columns(
			cell_paths, 'Main', 'time',
			fun=lambda x: (x[-1] - x[0]) / 60.).squeeze()
		dts = np.zeros(len(dt))
		for i in range(len(dt)):
			if i == 0:
				gt = dt[i]
				dts[i] = gt
			else:
				gt = dt[i] + dts[i-1]
				dts[i] = gt

		# Protein Counts
		plt.subplot(2, 1, 1)
		for x in dts:
			plt.axvline(x=x, color='#bcbd22', linestyle='--', linewidth=2)
		if len(cistron_monomer_ids) == 1:
			plt.plot(time / 60., ip_monomer_counts,
					 label = cistron_monomer_ids[0])
		else:
			for m in range(len(cistron_monomer_ids)):
				plt.plot(time / 60., ip_monomer_counts[:,m],
						 label = cistron_monomer_ids[m])

		plt.xlabel("Time (min)")
		plt.ylabel("Protein Counts")
		plt.title(f"Protein Counts for Proteins of Interest in variant"
				  f" {variant}, seed {seed}")
		plt.legend()

		# mRNA Counts
		plt.subplot(2, 1, 2)
		for x in dts:
			plt.axvline(x=x, color='#bcbd22', linestyle='--', linewidth=2)
		if len(cistron_ids) == 1:
			plt.plot(time / 60., ip_mRNA_counts,
					 label=cistron_ids[0])
		else:
			for r in range(len(cistron_ids)):
				plt.plot(time / 60., ip_mRNA_counts[:,r],
						 label = cistron_ids[r])
		plt.xlabel("Time (min)")
		plt.ylabel("Cistron Counts")
		plt.title(f"mRNA Counts for Proteins of Interest in variant {variant},"
				  f" seed {seed}")
		plt.legend()

		plt.subplots_adjust(hspace = 0.5, top = 0.95, bottom = 0.05)
		exportFigure(plt, plotOutDir, plotOutFileName + '_variant_' +
					 str(variant) + '_seed_' + str(seed), metadata)
		plt.close("all")

if __name__ == '__main__':
	Plot().cli()
