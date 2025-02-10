"""
Plot mRNA and free protein counts for genes across multiple generations.
"""

import pickle
import os
from matplotlib import pyplot as plt
import numpy as np
from models.ecoli.analysis import multigenAnalysisPlot
from wholecell.analysis.analysis_tools import (exportFigure,
	read_stacked_bulk_molecules, read_stacked_columns)
from wholecell.io.tablereader import TableReader

""" USER INPUTS """

# Add the proteins to be plotted here (there can be multiple):
interest_proteins = np.array([
	# 'EG11854-MONOMER[c]',
	'NG-GFP-MONOMER[c]',
])

""" END USER INPUTS """

class Plot(multigenAnalysisPlot.MultigenAnalysisPlot):

	def do_plot(self, seedOutDir, plotOutDir, plotOutFileName, simDataFile,
				validationDataFile, metadata):

		# Obtain the cell path data for the default seed and variant:
		cell_paths = self.ap.get_cells()
		variant = self.ap.get_variants()
		seed = self.ap.get_seeds()
		sim_dir = cell_paths[seed[0]]
		simOutDir = os.path.join(sim_dir, 'simOut')

		# Determine gene ids for the proteins of interest:
		with open(simDataFile, 'rb') as f:
			sim_data = pickle.load(f)
		monomer_sim_data = sim_data.process.translation.monomer_data.struct_array

		# extract info about the protein(s) from the monomer data:
		monomer_data_idxs = []
		for protein in interest_proteins:
			monomer_idx = np.where(monomer_sim_data['id'] == protein)
			monomer_idx = monomer_idx[0][0]
			monomer_data_idxs.append(monomer_idx)
		# monomer_data includes helpful info: monomer_data.dtype.names
		monomer_data = monomer_sim_data[monomer_data_idxs]
		monomer_ids = monomer_data['id']
		cistron_ids = monomer_data['cistron_id']

		# create a dictionary to map cistron ids to monomer ids
		cistron_monomer_id_dict = dict(zip(monomer_sim_data['cistron_id'],
										monomer_sim_data['id']))
		cistron_monomer_ids = [cistron_monomer_id_dict.get(mRNA_id)
								for mRNA_id in cistron_ids]

		# Extract mRNA indexes for each gene of interest
		mRNA_counts_reader = TableReader(os.path.join(simOutDir, 'RNACounts'))
		mRNA_idx_dict = {rna: i for i, rna in enumerate(
			mRNA_counts_reader.readAttribute('mRNA_cistron_ids'))}
		mRNA_indexes = [mRNA_idx_dict.get(mRNA_id) for mRNA_id in
						cistron_ids]

		# Load data
		time = read_stacked_columns(
			cell_paths, 'Main', 'time', ignore_exception=True)
		(monomer_counts,) = read_stacked_bulk_molecules(
			cell_paths, monomer_ids, ignore_exception=True)
		mRNA_counts = read_stacked_columns(
			cell_paths, 'RNACounts', 'mRNA_cistron_counts',
			ignore_exception=True)[:, mRNA_indexes]

		# Generate plots:
		plt.figure(figsize = (8.5, 11))

		# Extract doubling times from cells within this variant index
		dt = read_stacked_columns(
			cell_paths, 'Main', 'time',
			fun=lambda x: (x[-1] - x[0]) / 60.).squeeze()
		# Extract the end time of each generation
		dts = np.zeros(len(dt))
		for i in range(len(dt)):
			if i == 0:
				gt = dt[i]
				dts[i] = gt
			else:
				gt = dt[i] + dts[i-1]
				dts[i] = gt

		# Plot Protein Counts
		plt.subplot(2, 1, 1)
		for x in dts:
			plt.axvline(x=x, color='#bcbd22', linestyle='--', linewidth=2)
		if len(cistron_monomer_ids) == 1:
			plt.plot(time / 60., monomer_counts,
					 label = cistron_monomer_ids[0], alpha=0.5)
		else:
			for m in range(len(cistron_monomer_ids)):
				plt.plot(time / 60., monomer_counts[:,m],
						 label = cistron_monomer_ids[m], alpha=0.5)

		plt.xlabel("Time (min)")
		plt.ylabel("Protein Counts")
		plt.title(f"Free monomer counts for proteins of interest in variant"
				  f" {variant[0]}, seed {seed[0]}")
		plt.legend()

		# Plot mRNA Counts
		plt.subplot(2, 1, 2)
		for x in dts:
			plt.axvline(x=x, color='#bcbd22', linestyle='--', linewidth=2)
		if len(cistron_ids) == 1:
			plt.plot(time / 60., mRNA_counts,
					 label=cistron_ids[0], alpha=0.5)
		else:
			for r in range(len(cistron_ids)):
				plt.plot(time / 60., mRNA_counts[:,r],
						 label = cistron_ids[r], alpha=0.5)

		plt.xlabel("Time (min)")
		plt.ylabel("Cistron Counts")
		plt.title(f"mRNA counts for proteins of interest in variant "
				  f"{variant[0]}, seed {seed[0]}")
		plt.legend()

		# export plot
		plt.subplots_adjust(hspace = 0.5, top = 0.95, bottom = 0.05)
		exportFigure(plt, plotOutDir, plotOutFileName + '_multigenPlot_variant_'
					 + str(variant[0]) + '_seed_' + str(seed[0]) + '_geneIDs_' +
					 '_geneIDs_' + str(interest_proteins), metadata)
		plt.close("all")

if __name__ == '__main__':
	Plot().cli()
