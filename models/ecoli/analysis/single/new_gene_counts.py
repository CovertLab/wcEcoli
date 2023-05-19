"""
Plot mRNA and protein counts for new genes
"""

import os
import pickle

from matplotlib import pyplot as plt

from wholecell.io.tablereader import TableReader
from wholecell.analysis.analysis_tools import exportFigure
from models.ecoli.analysis import singleAnalysisPlot


class Plot(singleAnalysisPlot.SingleAnalysisPlot):
	def do_plot(self, simOutDir, plotOutDir, plotOutFileName, simDataFile,
				validationDataFile, metadata):

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
				  "were "
				  "found.")
			return
		assert len(new_gene_monomer_ids) == len(new_gene_mRNA_ids),\
			'number of new gene monomers and mRNAs should be equal'

		# Extract mRNA indexes for each new gene
		mRNA_counts_reader = TableReader(os.path.join(simOutDir,
													  'RNACounts'))
		mRNA_idx_dict = {rna[:-3]: i for i, rna in enumerate(
			mRNA_counts_reader.readAttribute('mRNA_ids'))}
		new_gene_mRNA_indexes = [mRNA_idx_dict.get(mRNA_id) for mRNA_id in
								 new_gene_mRNA_ids]

		# Extract protein indexes for each new gene
		monomer_counts_reader = TableReader(
			os.path.join(simOutDir, "MonomerCounts"))
		monomer_idx_dict = {monomer: i for i, monomer in
							enumerate(monomer_counts_reader.readAttribute(
								'monomerIds'))}
		new_gene_monomer_indexes = [monomer_idx_dict.get(monomer_id) for
									monomer_id in new_gene_monomer_ids]

		# Extract mRNA and monomer counts for each new gene
		mRNA_counts = mRNA_counts_reader.readColumn('mRNA_counts')
		monomer_counts = monomer_counts_reader.readColumn('monomerCounts')
		new_gene_mRNA_counts = mRNA_counts[:, new_gene_mRNA_indexes]
		new_gene_monomer_counts = monomer_counts[:, new_gene_monomer_indexes]

		main_reader = TableReader(os.path.join(simOutDir, "Main"))
		initialTime = main_reader.readAttribute("initialTime")
		time = main_reader.readColumn("time") - initialTime

		# Plotting
		plt.figure(figsize = (8.5, 11))

		# Protein Counts
		plt.subplot(2, 1, 1)
		for m in range(len(new_gene_monomer_ids)):
			plt.plot(time / 60., new_gene_monomer_counts[:,m],
					 label = new_gene_monomer_ids[m])
		plt.xlabel("Time (min)")
		plt.ylabel("Protein Counts")
		plt.title("New Gene Protein Counts")
		plt.legend()

		# mRNA Counts
		plt.subplot(2, 1, 2)
		for r in range(len(new_gene_mRNA_ids)):
			plt.plot(time / 60., new_gene_mRNA_counts[:,r],
					 label = new_gene_mRNA_ids[r])
		plt.xlabel("Time (min)")
		plt.ylabel("mRNA Counts")
		plt.title("New Gene mRNA Counts")
		plt.legend()

		plt.subplots_adjust(hspace = 0.5, top = 0.95, bottom = 0.05)
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close("all")


if __name__ == "__main__":
	Plot().cli()

