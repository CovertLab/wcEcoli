"""
Template for parca analysis plots

TODO:
	plot positive/negative directionality matches
	show absolute number (matched/unmatched) instead of fraction of matches
	compare different init options (no ppGpp, no attenuation etc)
	use li data or new proteome data
"""

import pickle

from matplotlib import pyplot as plt
import numpy as np

from models.ecoli.analysis import parcaAnalysisPlot
from reconstruction.ecoli.initialization import create_bulk_container
from wholecell.analysis.analysis_tools import exportFigure


def get_mw(mw, molecules):
	return np.array([mw.get(mol) for mol in molecules])

def get_monomers(molecules, get_stoich):
	return [monomer for mol in molecules for monomer in get_stoich(mol)['subunitIds']]

def get_container_counts(container, molecules, ids):
	total_counts = container.counts(ids).sum()
	return [container.counts(mol) / total_counts for mol in molecules]

def get_validation_counts(counts, molecules, ids):
	total_counts = np.sum([counts[m] for m in ids])
	return [np.array([counts.get(mol, 0) / total_counts for mol in molecule_group]) for molecule_group in molecules]

def compare_counts(condition1, condition2):
	return [np.log2(c1 / c2) for c1, c2 in zip(condition1, condition2)]

def compare_to_validation(parca, validation):
	matches = [(np.sign(p) == np.sign(v))[np.isfinite(p) & np.isfinite(v)] for p, v in zip(parca, validation)]
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
		translation = sim_data.process.translation
		mol_ids = sim_data.molecule_ids
		get_stoich = sim_data.process.complexation.get_monomers
		aa_enzymes = [e for e in metabolism.aa_enzymes if e != 'PUTA-CPLXBND[c]']  # PutA is already accounted for as PUTA-CPLX and no easy way to get monomers from both complexes and equilibrium molecules (like PUTA-CPLXBND)

		# Validation data
		rich_validation = {}
		basal_validation = {}
		for p in validation_data.protein.schmidt2015Data:
			monomer = p['monomerId']
			rich_validation[monomer] = p['LB_counts']
			basal_validation[monomer] = p['glucoseCounts']

		common_monomers = [m for m in translation.monomer_data['id'] if rich_validation.get(m, 0) != 0 and basal_validation.get(m, 0) != 0]

		# Select molecule groups of interest
		monomer_ids = (
			get_monomers(aa_enzymes, get_stoich),
			get_monomers(transcription.synthetase_names, get_stoich),
			get_monomers([mol_ids.RelA, mol_ids.SpoT], get_stoich),
		)
		group_labels = [
			'AA enzymes',
			'Synthetases',
			'ppGpp molecules',
		]

		# Expected bulk containers in different conditions
		options = {'ppgpp_regulation': True, 'trna_attenuation': True}  # TODO: iterate on different options
		rich_container = create_bulk_container(sim_data, condition='with_aa', form_complexes=False, **options)
		basal_container = create_bulk_container(sim_data, form_complexes=False, **options)
		rich_counts = get_container_counts(rich_container, monomer_ids, common_monomers)
		basal_counts = get_container_counts(basal_container, monomer_ids, common_monomers)
		parca_compare = compare_counts(rich_counts, basal_counts)

		# Validation mass fractions
		rich_counts_validation = get_validation_counts(rich_validation, monomer_ids, common_monomers)
		basal_counts_validation = get_validation_counts(basal_validation, monomer_ids, common_monomers)
		validation_compare = compare_counts(rich_counts_validation, basal_counts_validation)

		comparison = compare_to_validation(parca_compare, validation_compare)

		_, (bar_ax, scat_ax) = plt.subplots(2, 1, figsize=(5, 10))

		bar_ax.bar(group_labels, comparison)
		bar_ax.set_ylim([0, 1])
		self.remove_border(bar_ax)

		for p, v in zip(parca_compare, validation_compare):
			scat_ax.plot(p, v, 'x')
		scat_ax.axhline(0, color='k', linestyle='--', linewidth=0.5)
		scat_ax.axvline(0, color='k', linestyle='--', linewidth=0.5)
		scat_ax.set_xlabel('Sim log2 fold change\nfrom minimal to rich', fontsize=6)
		scat_ax.set_ylabel('Validation log2 fold change\nfrom minimal to rich', fontsize=8)
		scat_ax.tick_params(labelsize=6)
		self.remove_border(scat_ax)

		plt.tight_layout()
		exportFigure(plt, plot_out_dir, plot_out_filename, metadata)
		plt.close('all')


if __name__ == "__main__":
	Plot().cli()
