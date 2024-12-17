"""
Template for cohort analysis plots
"""

import pickle
import os

from matplotlib import pyplot as plt
from matplotlib import cm
import numpy as np
import csv

from models.ecoli.analysis import cohortAnalysisPlot
from wholecell.analysis.analysis_tools import (exportFigure,
	read_bulk_molecule_counts, read_stacked_bulk_molecules, read_stacked_columns)
from wholecell.io.tablereader import TableReader
from wholecell.containers.bulk_objects_container import BulkObjectsContainer

complexToMonomer = {
	"CPLX0-7620[c]": "PD00260[c]", # CPLX0-7620's monomer is EG10359-MONOMER, which is ID'ed as PD00260 (proteins.tsv)
	"CPLX0-8801[c]": "G6420-MONOMER[c]",
	"CPLX0-7677[c]": "EG11969-MONOMER[c]",
	"CPLX0-7702[c]": "G6263-MONOMER[c]",
	"CPLX0-7701[c]": "G6420-MONOMER[c]",
	}

monomerToTranslationMonomer = {
	"MONOMER0-1781[c]": "EG11519-MONOMER[c]", # MONOMER0-1781 is a complex, EG11519 is its monomer
	"EG10359-MONOMER[c]": "PD00260[c]", # EG10359-MONOMER is not the ID of fur monomer, it's PD00260 (proteins.tsv)
	}

IGNORE_FIRST_N_GENS = 2

class Plot(cohortAnalysisPlot.CohortAnalysisPlot):
	def do_plot(self, variantDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		with open(simDataFile, 'rb') as f:
			sim_data = pickle.load(f)
		with open(validationDataFile, 'rb') as f:
			validation_data = pickle.load(f)


		if self.ap.n_generation <= IGNORE_FIRST_N_GENS:
			print('Skipping analysis - not enough generations run.')
			return


		# Get list of cistron IDs from sim_data
		cistron_data = sim_data.process.transcription.cistron_data
		cistron_ids = cistron_data['id']

		# Filter list for cistron IDs with associated protein ids
		cistron_id_to_protein_id = {
			protein['cistron_id']: protein['id']
			for protein in sim_data.process.translation.monomer_data
		}
		mRNA_cistron_ids = [
			cistron_id for cistron_id in cistron_ids
			if cistron_id in cistron_id_to_protein_id]

		# Get IDs of associated monomers and genes
		monomer_ids = [
			cistron_id_to_protein_id.get(cistron_id, None)
			for cistron_id in mRNA_cistron_ids]
		cistron_id_to_gene_id = {
			cistron['id']: cistron['gene_id'] for cistron in cistron_data
		}
		gene_ids = [
			cistron_id_to_gene_id[cistron_id]
			for cistron_id in mRNA_cistron_ids]


		# Ignore data from predefined number of generations per seed
		cell_paths = self.ap.get_cells(
			generation=np.arange(IGNORE_FIRST_N_GENS, self.ap.n_generation),
			only_successful=True)

		simOutDir = os.path.join(cell_paths[0], 'simOut')
		rna_counts_reader = TableReader(os.path.join(simOutDir, 'RNACounts'))
		mRNA_cistron_ids_rna_counts_table = rna_counts_reader.readAttribute(
			'mRNA_cistron_ids')

		# Get indexes of mRNA cistrons in this subcolumn
		mRNA_cistron_id_to_index = {
			cistron_id: i for (i, cistron_id)
			in enumerate(mRNA_cistron_ids_rna_counts_table)
		}
		mRNA_cistron_indexes = np.array([
			mRNA_cistron_id_to_index[cistron_id] for cistron_id
			in mRNA_cistron_ids
		])

		# Get subcolumn for monomer IDs in monomer counts table
		monomer_counts_reader = TableReader(
			os.path.join(simOutDir, 'MonomerCounts'))
		monomer_ids_monomer_counts_table = monomer_counts_reader.readAttribute(
			'monomerIds')
		monomer_id_to_index = {
			monomer_id: i for (i, monomer_id)
			in enumerate(monomer_ids_monomer_counts_table)
			}
		monomer_indeces = np.array([
			monomer_id_to_index[monomer_id] for monomer_id in monomer_ids
		])



		# Get all protein ids reqiured
		ids_complexation_complexes = sim_data.process.complexation.ids_complexes  # Only complexes
		ids_equilibrium = sim_data.process.equilibrium.molecule_names  # Complexes of proteins + small molecules, small molecules, protein monomers
		ids_equilibrium_complexes = sim_data.process.equilibrium.ids_complexes  # Only complexes
		ids_protein = sorted(set(ids_complexation_complexes + ids_equilibrium + list(monomer_id_to_index.keys()))) #set of unique complexes and free monomer ids, and small molecules

		bulkContainer = BulkObjectsContainer(ids_protein, dtype=np.float64)

		view_complexation_complexes = bulkContainer.countsView(ids_complexation_complexes)
		view_equilibrium = bulkContainer.countsView(ids_equilibrium)
		view_equilibrium_complexes = bulkContainer.countsView(ids_equilibrium_complexes)
		view_translation = bulkContainer.countsView(monomer_id_to_index.keys())


		# Identify monomers that are not in complexes, in one complex, or found in multiple complexes
		monomersInvolvedInManyComplexes = []
		monomersInvolvedInComplexes = []
		One_to_one_monomer_complex_dict = {}
		for complexId in ids_complexation_complexes:
			subunitIds = sim_data.process.complexation.get_monomers(complexId)["subunitIds"]
			for subunitId in subunitIds:
				if subunitId in monomersInvolvedInComplexes:
					monomersInvolvedInManyComplexes.append(subunitId)
				monomersInvolvedInComplexes.append(subunitId)
				One_to_one_monomer_complex_dict[subunitId] = complexId
		monomersInvolvedInManyComplexes_set = set(monomersInvolvedInManyComplexes)
		monomersInvolvedInManyComplexes_dict = {}
		for x in monomersInvolvedInManyComplexes_set:
			monomersInvolvedInManyComplexes_dict[x] = {}

		all_monomers_set = set(monomer_id_to_index.keys())
		monomersInvolvedInComplexes_set = set(monomersInvolvedInComplexes)
		monomers_in_one_complex_set = monomersInvolvedInComplexes_set.symmetric_difference(monomersInvolvedInManyComplexes_set)
		monomers_not_in_complex_set = all_monomers_set.symmetric_difference(monomersInvolvedInComplexes_set)


		# Get average (over timesteps) counts including all seeds and generations (ie. All cells)
		avgProteinCounts_forAllCells = np.zeros(len(mRNA_cistron_ids), np.float64)
		for i, simDir in enumerate(cell_paths):
			simOutDir = os.path.join(simDir, "simOut")

			# Account for bulk molecules
			bulkMolecules = TableReader(os.path.join(simOutDir, "BulkMolecules"))
			bulkMoleculeCounts = bulkMolecules.readColumn("counts")
			moleculeIds = bulkMolecules.readAttribute("objectNames")
			proteinIndexes = np.array([moleculeIds.index(moleculeId) for moleculeId in ids_protein], int)
			proteinCountsBulk = bulkMoleculeCounts[:, proteinIndexes]

			bulkMolecules.close()

			bulkContainer.countsIs(proteinCountsBulk.mean(axis=0))

			# Unique molecules
			uniqueMoleculeCounts = TableReader(os.path.join(simOutDir, "UniqueMoleculeCounts"))
			ribosomeIndex = uniqueMoleculeCounts.readAttribute("uniqueMoleculeIds").index('active_ribosome')
			rnaPolyIndex = uniqueMoleculeCounts.readAttribute("uniqueMoleculeIds").index('active_RNAP')
			nActiveRibosome = uniqueMoleculeCounts.readColumn("uniqueMoleculeCounts")[:, ribosomeIndex]
			nActiveRnaPoly = uniqueMoleculeCounts.readColumn("uniqueMoleculeCounts")[:, rnaPolyIndex]
			uniqueMoleculeCounts.close()

			# Account for unique molecules
			bulkContainer.countsInc(nActiveRibosome.mean(),
									[sim_data.molecule_ids.s30_full_complex, sim_data.molecule_ids.s50_full_complex])
			bulkContainer.countsInc(nActiveRnaPoly.mean(), [sim_data.molecule_ids.full_RNAP])

			# Account for small-molecule bound complexes
			view_equilibrium.countsInc(
				np.dot(sim_data.process.equilibrium.stoich_matrix_monomers(), view_equilibrium_complexes.counts() * -1))

			# Average counts of monomers
			avgMonomerCounts = view_translation.counts()

			# Get counts of "functional units" (ie. complexed forms)
			avgComplexCounts = view_complexation_complexes.counts()

			for j, complexId in enumerate(ids_complexation_complexes):
				# Map all subsunits to the average counts of the complex (ignores counts of monomers)
				# Some subunits are involved in multiple complexes - these cases are kept track
				subunitIds = sim_data.process.complexation.get_monomers(complexId)["subunitIds"]
				for subunitId in subunitIds:
					if subunitId not in monomer_id_to_index:
						if subunitId in monomerToTranslationMonomer:
							# couple monomers have different ID in ids_translation
							subunitId = monomerToTranslationMonomer[subunitId]
						elif "CPLX" in subunitId:
							# few transcription factors are complexed with ions
							subunitId = complexToMonomer[subunitId]
						elif "RNA" in subunitId:
							continue
					else:
						if subunitId not in monomersInvolvedInManyComplexes_set:
							avgMonomerCounts[monomer_id_to_index[subunitId]] = avgComplexCounts[j]

						else:
							if complexId not in monomersInvolvedInManyComplexes_dict[subunitId]:
								monomersInvolvedInManyComplexes_dict[subunitId][complexId] = 0
							monomersInvolvedInManyComplexes_dict[subunitId][complexId] += avgComplexCounts[j]

			# Store
			avgProteinCounts_forAllCells += avgMonomerCounts

		# Per cell
		avgProteinCounts_perCell = avgProteinCounts_forAllCells / float(len(cell_paths))


		# Plot
		fig = plt.figure(figsize=(12, 4))

		# Create a colormap
		#cmap = cm.get_cmap('Set1')
		#colors = [cmap(i) for i in monomer_indeces]

		def plot_monomer_counts_in_mult_complexes(ax, monomersInvolvedInManyComplexes_dict, monomer_id_to_index):
			for monomer, complexes in monomersInvolvedInManyComplexes_dict.items():
				if monomer in monomer_id_to_index.keys():
					for complexId in complexes:
						avgComplexCount = monomersInvolvedInManyComplexes_dict[monomer][complexId] / float(len(cell_paths))
						if avgComplexCount > 5000:
							ax.scatter(complexId[:-3], avgComplexCount, alpha=0.5, marker=".",  color= 'red')
							ax.annotate(f"{monomer[:-3]}", (complexId[:-3], avgComplexCount))

						else:
							ax.scatter(complexId[:-3], avgComplexCount, alpha=0.5, marker=".", color='grey')

			ax.set_ylabel('Avg Counts')
			ax.set_xlabel('Protein complexes')
			ax.set_xticks([0])
			ax.set_xticklabels([])

		def plot_monomer_counts(ax, monomer_set, monomer_id_to_index, xlabel):
			for monomer in monomer_set:
				if monomer in monomer_id_to_index.keys():
					monomer_index = monomer_id_to_index[monomer]
					if avgProteinCounts_perCell[monomer_index] > 25000:
						ax.scatter(monomer[:-3], avgProteinCounts_perCell[monomer_index], alpha=0.5, marker=".", lw=0.,
								  color='red')
						ax.annotate(f"{monomer[:-3]}", (monomer[:-3], avgProteinCounts_perCell[monomer_index]))
					else:
						ax.scatter(monomer[:-3], avgProteinCounts_perCell[monomer_index], alpha = 0.5, marker = ".", lw = 0., color = 'grey')
			ax.set_ylabel('Avg Counts')
			ax.set_xlabel(xlabel)
			ax.set_xticks([0])
			ax.set_xticklabels([])

		ax1 = plt.subplot(1, 3, 1)
		ax2 = plt.subplot(1, 3, 2)
		ax3 = plt.subplot(1, 3, 3)

		plot_monomer_counts_in_mult_complexes(ax1, monomersInvolvedInManyComplexes_dict, monomer_id_to_index)
		plot_monomer_counts(ax2, monomers_in_one_complex_set, monomer_id_to_index, 'Monomers in one complex')
		plot_monomer_counts(ax3, monomers_not_in_complex_set, monomer_id_to_index, 'Monomers not in complex')


		plt.subplots_adjust(hspace=0.5, wspace=0.5, left=0.1, bottom=0.1, top=0.9, right=0.95)
		# Create a list of patches to represent each category
		#handles = [mpatches.Patch(color=cmap(i), label=monomer) for monomer, i in monomer_id_to_index.items()]
		#fig.legend(handles=handles, bbox_to_anchor=(1, 0.7))

		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close("all")

		def write_table_for_monomers_in_mult_complexes(monomersInvolvedInManyComplexes_dict, monomer_id_to_index):
			# Write data to table
			with open(os.path.join(plotOutDir, plotOutFileName + 'monomers_in_multiple_protein_complexes.tsv'),
					  'w') as f:
				writer = csv.writer(f, delimiter='\t')
				writer.writerow(['monomer_id', 'complex_id', 'avg_count_per_cell'])

				for monomer, complexes in monomersInvolvedInManyComplexes_dict.items():
					if monomer in monomer_id_to_index.keys():
						for complexId in complexes:
							avgComplexCount = monomersInvolvedInManyComplexes_dict[monomer][complexId] / float(len(cell_paths))
							writer.writerow([
								monomer[:-3], complexId[:-3], avgComplexCount
							])

		def write_table_for_monomer_in_one_complex(avgProteinCounts_perCell, monomers_in_one_complex_set, monomer_id_to_index, One_to_one_monomer_complex_dict):
			# Write data to table
			with open(os.path.join(plotOutDir, plotOutFileName + 'Monomers_in_one_complex.tsv'),
					  'w') as f:
				writer = csv.writer(f, delimiter='\t')
				writer.writerow(['monomer_id', 'complex_id', 'avg_count_per_cell'])

				for monomer in monomers_in_one_complex_set:
					if monomer in monomer_id_to_index.keys():
						monomer_index = monomer_id_to_index[monomer]
						complexId = One_to_one_monomer_complex_dict[monomer]
						monomer_counts = avgProteinCounts_perCell[monomer_index]
						writer.writerow([
							monomer[:-3], complexId[:-3], monomer_counts
						])

		def write_table_for_monomer_not_complexed(avgProteinCounts_perCell, monomers_not_in_complex_set, monomer_id_to_index):
			# Write data to table
			with open(os.path.join(plotOutDir, plotOutFileName + 'Monomers_not_in_complex.tsv'),
					  'w') as f:
				writer = csv.writer(f, delimiter='\t')
				writer.writerow(['monomer_id', 'complex_id', 'avg_count_per_cell'])

				for monomer in monomers_not_in_complex_set:
					if monomer in monomer_id_to_index.keys():
						monomer_index = monomer_id_to_index[monomer]
						complexId = None
						monomer_counts = avgProteinCounts_perCell[monomer_index]
						writer.writerow([
							monomer[:-3], complexId, monomer_counts
						])

		write_table_for_monomers_in_mult_complexes(monomersInvolvedInManyComplexes_dict, monomer_id_to_index)
		write_table_for_monomer_in_one_complex(avgProteinCounts_perCell, monomers_in_one_complex_set,
											   monomer_id_to_index, One_to_one_monomer_complex_dict)
		write_table_for_monomer_not_complexed(avgProteinCounts_perCell, monomers_not_in_complex_set,
											  monomer_id_to_index)


if __name__ == '__main__':
	Plot().cli()
