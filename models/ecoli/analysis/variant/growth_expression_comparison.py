"""
Template for variant analysis plots
"""

import pickle
import os

from matplotlib import pyplot as plt
import numpy as np

from models.ecoli.analysis import variantAnalysisPlot
from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from wholecell.analysis.analysis_tools import (exportFigure,
	read_bulk_molecule_counts, read_stacked_bulk_molecules, read_stacked_columns)
from wholecell.io.tablereader import TableReader


def read_counts(cell_paths):
	ids = TableReader(os.path.join(cell_paths[0], 'simOut', 'MonomerCounts')).readAttribute('monomerIds')
	counts = read_stacked_columns(cell_paths, 'MonomerCounts', 'monomerCounts').T
	return dict(zip(ids, counts))

def get_mw(mw, molecules):
	return np.array([mw.get(mol) for mol in molecules])

def get_monomers(molecules, get_stoich):
	return [monomer for mol in molecules for monomer in get_stoich(mol)['subunitIds']]

def get_sim_counts(counts, molecules, ids):
	total_counts = np.sum([counts[m] for m in ids], axis=0)
	return [np.array([counts.get(mol, np.zeros_like(total_counts)) / total_counts for mol in molecule_group]).mean(1) for molecule_group in molecules]

def get_validation_counts(counts, molecules, ids):
	total_counts = np.sum([counts[m] for m in ids])
	return [np.array([counts.get(mol, 0) / total_counts for mol in molecule_group]) for molecule_group in molecules]

def compare_counts(condition1, condition2):
	with np.errstate(divide='ignore', invalid='ignore'):
		fc = [np.log2(c1 / c2) for c1, c2 in zip(condition1, condition2)]
	return fc

def compare_to_validation(parca, validation):
	matches = [(np.sign(p) == np.sign(v))[np.isfinite(p) & np.isfinite(v)] for p, v in zip(parca, validation)]
	return [match.sum() / match.shape[0] for match in matches]


class Plot(variantAnalysisPlot.VariantAnalysisPlot):
	def do_plot(self, inputDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		ap = AnalysisPaths(inputDir, variant_plot=True)
		variants = ap.get_variants()

		minimal_sim = None
		rich_sim = None
		for variant in variants:
			with open(ap.get_variant_kb(variant), 'rb') as f:
				variant_sim_data = pickle.load(f)

			cell_paths = ap.get_cells(variant=[variant])
			if variant_sim_data.condition == 'basal':
				minimal_sim = read_counts(cell_paths)
			elif variant_sim_data.condition == 'with_aa':
				rich_sim = read_counts(cell_paths)

		if minimal_sim is None or rich_sim is None:
			print(f'Do not have minimal and rich variant for {plotOutFileName} - skipping...')
			return

		with open(simDataFile, 'rb') as f:
			sim_data = pickle.load(f)
		with open(validationDataFile, 'rb') as f:
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
		rich_counts = get_sim_counts(rich_sim, monomer_ids, common_monomers)
		basal_counts = get_sim_counts(minimal_sim, monomer_ids, common_monomers)
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
		scat_ax.set_xlabel('Sim log2 fold change\nfrom minimal to rich', fontsize=8)
		scat_ax.set_ylabel('Validation log2 fold change\nfrom minimal to rich', fontsize=8)
		scat_ax.tick_params(labelsize=6)
		self.remove_border(scat_ax)

		plt.tight_layout()
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close('all')


if __name__ == "__main__":
	Plot().cli()
