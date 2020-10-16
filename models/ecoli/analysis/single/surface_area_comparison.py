"""
Comparison of surface area calculations
"""

from __future__ import absolute_import, division, print_function

import os

from matplotlib import pyplot as plt
import numpy as np

from models.ecoli.analysis import singleAnalysisPlot
from wholecell.analysis.analysis_tools import exportFigure
from wholecell.analysis.analysis_tools import read_bulk_molecule_counts
from wholecell.io.tablereader import TableReader
from wholecell.analysis.plotting_tools import CMAP_COLORS_255


CMAP_COLORS = [[shade/255. for shade in color] for color in CMAP_COLORS_255]
CMAP_OVER = [0, 1, 0.75]

'''
References:
- width of 0.5 um: Neidhart et al., Physiology of the Bacterial Cell, Chapter 1
	- strain: B/R
- width of 0.73 um: Cluzel et al., Nucleic Acids Research (2008)
	- strain: K-12 Frag1
- surface areas per molecule: Harvard Bionumbers (Properties of the major 
components of outer membrane)
	-sources: Neidhart et al., Escherichia coli and Salmonella: Cellular and 
	Molecular Biology, Chapter 8 ; Smit et al., Journal of Bacteriology (1975)
- assumption that half of all phospholipids produced end up in outer membrane:
Osborn et al., Journal of Biological Chemistry (1972)
- surface area validation point of 6 um^2: Neidhart et al., Physiology of the 
Bacterial Cell, Chapter 1
	- strain: B/R  
- average e. coli time point 44% of the way through cell cycle: Neidhart et al., 
Physiology of the Bacterial Cell, Chapter 1

Note:
- cell volume: this value from the model is derived from the density parameter
	- strain: ML308
'''

width = [0.5, 0.73]		# in um
surface_area_per_molecule = {
	'LPS': 1.42E-6,
	'porins_and_ompA': 9E-6,
	'phospholipids': 4.71E-07,
	'lipoprotein': 7.14E-07,
}

ordered_outer_mem_protein_ids = [
	# phospholipids
	'CPD-12819[c]',
    'CPD-12824[c]',
    'CPD-8260[c]',

	# LPS
    'CPD0-939[c]',

	# porins and ompA
    'EG10669-MONOMER[i]',
	'CPLX0-7533[e]',
	'CPLX0-7534[e]',

	# lipoprotein
	'EG10544-MONOMER[o]',
]

outer_mem_proteins = {
	# phospholipids
	'CPD-12819[c]': 'phosphatidylethanolamine',
    'CPD-12824[c]' : 'cardiolipin',
    'CPD-8260[c]': 'phosphatidylglycerol',

	# LPS
    'CPD0-939[c]': 'LPS',

	# porins and ompA
    'EG10669-MONOMER[i]': 'ompA',
	'CPLX0-7533[e]': 'outer membrane porin C',
	'CPLX0-7534[e]': 'outer membrane porin F',

	# lipoprotein
	'EG10544-MONOMER[o]': 'murein lipoprotein',
}


class Plot(singleAnalysisPlot.SingleAnalysisPlot):
	def do_plot(self, simOutDir, plotOutDir, plotOutFileName, simDataFile,
				validationDataFile, metadata):
		# Listeners used
		main_reader = TableReader(os.path.join(simOutDir, 'Main'))
		mass = TableReader(os.path.join(simOutDir, "Mass"))

		# Load data
		volume = mass.readColumn("cellVolume")
		initial_time = main_reader.readAttribute('initialTime')
		time = (main_reader.readColumn('time') - initial_time) / 60
		(counts,) = read_bulk_molecule_counts(simOutDir, (ordered_outer_mem_protein_ids,))
		counts = counts.astype(float).T


		# Calculate surface area based off of model volume
		radii = np.divide(width, 2)
		surface_area_model = np.zeros((np.shape(radii)[0], np.shape(volume)[0]))
		for i, radius in enumerate(radii):
			surface_area_sphere = 4 * np.pi * radius**2
			length = [(vol - ((4 / 3) * np.pi * radius**3))/(np.pi * radius**2)
					  + 2 * radius for vol in volume]
			surface_area_cylinder = [2 * np.pi * radius * l for l in length]
			surface_area_model[i] = [val + surface_area_sphere
									 for val in surface_area_cylinder]

		# calculate SA based off of molecule counts
		surface_area_LPS = surface_area_per_molecule['LPS'] * \
			counts[ordered_outer_mem_protein_ids.index(
			'CPD0-939[c]')]
		surface_area_porins_and_ompA = \
			surface_area_per_molecule['porins_and_ompA'] * (np.add(np.add(
			counts[ordered_outer_mem_protein_ids.index('CPLX0-7533[e]')],
			counts[ordered_outer_mem_protein_ids.index('CPLX0-7534[e]')]),
			counts[ordered_outer_mem_protein_ids.index('EG10669-MONOMER[i]')]))
		surface_area_phospholipids = \
			surface_area_per_molecule['phospholipids'] * 0.5 * (np.add(np.add(
			counts[ordered_outer_mem_protein_ids.index('CPD-12819[c]')],
			counts[ordered_outer_mem_protein_ids.index('CPD-12824[c]')]),
			counts[ordered_outer_mem_protein_ids.index('CPD-8260[c]')]))
		surface_area_lipoprotein = \
			surface_area_per_molecule['lipoprotein'] * counts[
			ordered_outer_mem_protein_ids.index('EG10544-MONOMER[o]')]

		surface_area_outer_leaflet = np.add(surface_area_LPS,
			surface_area_porins_and_ompA)
		surface_area_inner_leaflet = np.add(np.add(surface_area_phospholipids,
			surface_area_lipoprotein), surface_area_porins_and_ompA)

		plt.figure(figsize=(8.5, 11))
		plt.plot(time, surface_area_model[0], 'r-')
		plt.plot(time, surface_area_model[1], 'r--')
		plt.plot(time, surface_area_outer_leaflet, 'b-')
		plt.plot(time, surface_area_inner_leaflet, 'b--')
		plt.scatter(0.44*time[-1], 6)

		plt.xlabel('time (min)')
		plt.ylabel('surface area (um ^2)')
		plt.title('Comparison of surface area calculations')
		plt.legend(['Model derived SA: width = 0.5', 'Model derived SA: width = 0.73',
					'Molecule count derived SA of outer leaflet',
					'Molecule count derived SA of inner leaflet',
					'Average E. coli SA'])

		plt.tight_layout()
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close('all')


if __name__ == '__main__':
	Plot().cli()
