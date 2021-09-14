"""
Template for parca analysis plots
"""

import pickle

from matplotlib import pyplot as plt
import numpy as np

from models.ecoli.analysis import parcaAnalysisPlot
from reconstruction.ecoli.initialization import create_bulk_container
from wholecell.analysis.analysis_tools import exportFigure


def get_monomers(molecules, get_stoich):
	return [monomer for mol in molecules for monomer in get_stoich(mol)['subunitIds']]

def get_container_counts(container, molecules):
	return [container.counts(mol) for mol in molecules]

def get_validation_counts(counts, molecules):
	return [np.array([counts.get(mol, 0) for mol in molecule_group]) for molecule_group in molecules]

def compare_counts(condition1, condition2):
	return [np.sign(c1 - c2) for c1, c2 in zip(condition1, condition2)]

def compare_to_validation(parca, validation):
	matches = [p == v for p, v in zip(parca, validation)]
	return [match.sum() / match.shape[0] for match in matches]


class Plot(parcaAnalysisPlot.ParcaAnalysisPlot):
	def do_plot(self, input_dir, plot_out_dir, plot_out_filename, sim_data_file, validation_data_file, metadata):
		with open(sim_data_file, 'rb') as f:
			sim_data = pickle.load(f)
		with open(validation_data_file, 'rb') as f:
			validation_data = pickle.load(f)

		# sim_data attributes used
		metabolism = sim_data.process.metabolism
		transcription = sim_data.process.transcription
		mol_ids = sim_data.molecule_ids
		get_stoich = sim_data.process.complexation.get_monomers

		# Validation data
		proteomics = validation_data.protein.schmidt2015Data

		# Select molecule groups of interest
		monomer_ids = (
			get_monomers(metabolism.aa_enzymes, get_stoich),
			get_monomers(transcription.synthetase_names, get_stoich),
			get_monomers([mol_ids.RelA, mol_ids.SpoT], get_stoich),
		)
		group_labels = [
			'AA enzymes',
			'Synthetases',
			'ppGpp molecules',
		]

		# Expected bulk containers in different conditions
		rich_container = create_bulk_container(sim_data, condition='with_aa', form_complexes=False)
		basal_container = create_bulk_container(sim_data, form_complexes=False)
		rich_counts = get_container_counts(rich_container, monomer_ids)
		basal_counts = get_container_counts(basal_container, monomer_ids)
		parca_compare = compare_counts(rich_counts, basal_counts)

		# TODO: get mass fraction from validation
		rich_validation = {}
		basal_validation = {}
		for monomer, rich, minimal in zip(proteomics['monomerId'], proteomics['LB_counts'], proteomics['glucoseCounts']):
			rich_validation[monomer] = rich
			basal_validation[monomer] = minimal
		rich_counts_validation = get_validation_counts(rich_validation, monomer_ids)
		basal_counts_validation = get_validation_counts(basal_validation, monomer_ids)
		validation_compare = compare_counts(rich_counts_validation, basal_counts_validation)

		comparison = compare_to_validation(parca_compare, validation_compare)

		plt.figure()
		plt.bar(group_labels, comparison)
		plt.tight_layout()
		exportFigure(plt, plot_out_dir, plot_out_filename, metadata)
		plt.close('all')


if __name__ == "__main__":
	Plot().cli()
