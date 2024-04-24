"""
Calculates the k_eff values of the reaction between the rRNA promoters and RNA
polymerases.
"""

import numpy as np
import os
import pickle

from models.ecoli.analysis import variantAnalysisPlot
from wholecell.analysis.analysis_tools import (
	read_stacked_columns, read_stacked_bulk_molecules,)
from wholecell.io.tablereader import TableReader
from wholecell.utils import units

FONT_SIZE = 9

GENE_ID_TO_RRNA_OPERON_ID = {
	'EG30084': 'rrnA',
	'EG30085': 'rrnB',
	'EG30086': 'rrnC',
	'EG30087': 'rrnD',
	'EG30088': 'rrnE',
	'EG30089': 'rrnG',
    'EG30090': 'rrnH',
	}

class Plot(variantAnalysisPlot.VariantAnalysisPlot):
	def do_plot(self, inputDir, plotOutDir, plotOutFileName, simDataFile,
				validationDataFile, metadata):
		with open(simDataFile, 'rb') as f:
			sim_data = pickle.load(f)

		# Get IDs and constants
		s30_full_complex_id = sim_data.molecule_ids.s30_full_complex
		s50_full_complex_id = sim_data.molecule_ids.s50_full_complex
		inactive_rnap_id = sim_data.molecule_ids.full_RNAP
		n_avogadro = sim_data.constants.n_avogadro

		# Data extraction
		ribosome_concentrations = {}
		doubling_times = {}
		inactive_rnap_concentrations = {}
		rrna_gene_concentrations = {}
		variant_indexes = self.ap.get_variants()

		# Loop through all variant indexes
		for variant_index in variant_indexes:
			# Get all cells (within the generation range) of this variant index
			all_cells = self.ap.get_cells(
				variant=[variant_index],
				only_successful=True)

			if len(all_cells) == 0:
				continue

			# Get index of active ribosomes in unique molecule counts table
			sim_dir = all_cells[0]
			simOutDir = os.path.join(sim_dir, 'simOut')
			unique_molecule_counts_table = TableReader(os.path.join(
				simOutDir, "UniqueMoleculeCounts"))
			active_ribosome_index = unique_molecule_counts_table.readAttribute(
				"uniqueMoleculeIds").index('active_ribosome')

			# Get ribosome counts and volumes from first timestep of every cell
			active_ribosome_counts = read_stacked_columns(
				all_cells, 'UniqueMoleculeCounts','uniqueMoleculeCounts'
				)[:, active_ribosome_index]
			ribosome_subunit_counts = read_stacked_bulk_molecules(
				all_cells, ([s30_full_complex_id, s50_full_complex_id],),)[0]
			cell_volume = (units.L) * 1e-15 * read_stacked_columns(
				all_cells, 'Mass', 'cellVolume').flatten()

			# Calculate concentrations of all ribosomes in uM
			ribosome_conc_this_variant = ((1 / (n_avogadro * cell_volume)) * (
				active_ribosome_counts + ribosome_subunit_counts.min(axis=1)
				)).asNumber(units.mol / units.L) * 1e6

			# Get doubling times of this variant
			dt_this_variant = read_stacked_columns(
				all_cells, 'Main', 'time',
				fun=lambda x: (x[-1] - x[0]) / 60.).squeeze()

			# Get inactive rnap concentrations in uM
			(inactive_rnap_counts,) = read_stacked_bulk_molecules(
				all_cells, ([inactive_rnap_id],))
			inactive_rnap_conc_this_variant = ((1 / (n_avogadro * cell_volume))
				*inactive_rnap_counts).asNumber(units.mol / units.L) * 1e6

			# Get gene_ids attribute from reference cell path
			reference_cell_path = all_cells[0]
			sim_out_dir = os.path.join(reference_cell_path, 'simOut')
			rna_synth_prob_reader = TableReader(
				os.path.join(sim_out_dir, 'RnaSynthProb'))
			gene_ids = rna_synth_prob_reader.readAttribute('gene_ids')

			# Get indexes of 16S genes (first gene in each operon)
			rrna_gene_indexes = np.array([
				gene_ids.index(key) for key in GENE_ID_TO_RRNA_OPERON_ID.keys()
				])

			# Get copy numbers of 16S genes
			rrna_gene_copy_numbers_this_variant = read_stacked_columns(
				all_cells, 'RnaSynthProb', 'gene_copy_number',
				fun=lambda x: x[:, rrna_gene_indexes]).sum(axis=1)
			rrna_gene_conc_this_variant = ((1 / (n_avogadro * cell_volume))
				*rrna_gene_copy_numbers_this_variant).asNumber(units.mol / units.L) * 1e6

			ribosome_concentrations[variant_index] = ribosome_conc_this_variant.mean()
			doubling_times[variant_index] = dt_this_variant.mean()
			inactive_rnap_concentrations[variant_index] = inactive_rnap_conc_this_variant.mean()
			rrna_gene_concentrations[variant_index] = rrna_gene_conc_this_variant.mean()

		variant_indexes_sorted = [k for (k, v) in sorted(doubling_times.items(), key=lambda item: item[1])]

		with open(os.path.join(plotOutDir, plotOutFileName + '.tsv'), 'w') as f:
			f.write(
				'\t'.join([
					'Variant index',
					'Doubling time (min)',
					'Ribosome conc. (uM)',
					'Inactive RNAP conc. (uM)',
					'rRNA promoter conc. (uM)',
					'k_eff'
					]) + '\n'
				)

			for k in variant_indexes_sorted:
				k_eff = np.log(2) * ribosome_concentrations[k] / (
					doubling_times[k] * inactive_rnap_concentrations[k] * rrna_gene_concentrations[k])
				f.write(
					'\t'.join([
						str(k),
						str(doubling_times[k]),
						str(ribosome_concentrations[k]),
						str(inactive_rnap_concentrations[k]),
						str(rrna_gene_concentrations[k]),
						str(k_eff)
						]) + '\n'
					)
if __name__ == "__main__":
	Plot().cli()
