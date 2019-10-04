"""
Plot the Voronoi diagram of mass fractions

@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 09/27/2019
"""
from __future__ import absolute_import

import os
import cPickle
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

from wholecell.io.tablereader import TableReader
from wholecell.analysis.analysis_tools import exportFigure
from models.ecoli.analysis import singleAnalysisPlot
from wholecell.utils import units
from wholecell.utils import VoronoiPlotMain

VM = VoronoiPlotMain.VoronoiMaster()
CANVAS = np.array([[0, 0], [4, 0], [4, 4], [0, 4]]) #the overall shape of the plot
I_MAX = 75 #the number of iterations for optimizing the Voronoi diagram
ERR_THRES = 1E-6 #the error threshold set to break from optimizing process

class Plot(singleAnalysisPlot.SingleAnalysisPlot):
	def do_plot(self, simOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		if not os.path.isdir(simOutDir):
			raise Exception("simOutDir does not currently exist as a directory")

		if not os.path.exists(plotOutDir):
			os.mkdir(plotOutDir)

		with open(simDataFile, 'rb') as f:
			sim_data = cPickle.load(f)

		# Random seed
		np.random.seed(1)

		# Load data
		mass = TableReader(os.path.join(simOutDir, "Mass"))
		unique_molecule_counts = TableReader(os.path.join(simOutDir, "UniqueMoleculeCounts"))
		bulk_molecules = TableReader(os.path.join(simOutDir, "BulkMolecules"))
		bulk_molecule_counts = bulk_molecules.readColumn("counts")
		bulk_molecule_ids = bulk_molecules.readAttribute("objectNames")

		# Count free rRNAs
		rrna_ids_16s = sim_data.moleculeGroups.s30_16sRRNA
		rrna_ids_23s = sim_data.moleculeGroups.s50_23sRRNA
		rrna_ids_5s = sim_data.moleculeGroups.s50_5sRRNA
		rrna_indexes_16s = np.array(
			[bulk_molecule_ids.index(rRna) for rRna in rrna_ids_16s], np.int64)
		rrna_indexes_23s = np.array(
			[bulk_molecule_ids.index(rRna) for rRna in rrna_ids_23s], np.int64)
		rrna_indexes_5s = np.array(
			[bulk_molecule_ids.index(rRna) for rRna in rrna_ids_5s], np.int64)
		free_rrna_counts_16s = bulk_molecule_counts[:, rrna_indexes_16s]
		free_rrna_counts_23s = bulk_molecule_counts[:, rrna_indexes_23s]
		free_rrna_counts_5s = bulk_molecule_counts[:, rrna_indexes_5s]

		# Get the stoichiometry matrix of 50s & 30s subunits
		ribosome_50s_subunits = sim_data.process.complexation.getMonomers(
			sim_data.moleculeIds.s50_fullComplex)
		ribosome_30s_subunits = sim_data.process.complexation.getMonomers(
			sim_data.moleculeIds.s30_fullComplex)
		ribosome_stoich = np.hstack((ribosome_50s_subunits["subunitStoich"]
			, ribosome_30s_subunits["subunitStoich"]))

		# Count 50s & 30s subunits
		complex_ids = [sim_data.moleculeIds.s50_fullComplex]
		complex_indexes = np.array(
			[bulk_molecule_ids.index(comp) for comp in complex_ids], np.int64)
		[complex_50s_Counts] = bulk_molecule_counts[:, complex_indexes].T
		complex_ids = [sim_data.moleculeIds.s30_fullComplex]
		complex_indexes = np.array(
			[bulk_molecule_ids.index(comp) for comp in complex_ids], np.int64)
		[complex_30s_Counts] = bulk_molecule_counts[:, complex_indexes].T
		n_50s_stoich = np.tensordot(complex_50s_Counts, 
			ribosome_50s_subunits["subunitStoich"], axes = 0)
		n_30s_stoich = np.tensordot(complex_30s_Counts, 
			ribosome_50s_subunits["subunitStoich"], axes = 0)
		n_subunit_ribosome_stoich = np.hstack((n_50s_stoich, n_30s_stoich))

		# Count active ribosomes
		ribosome_subunit_ids = (ribosome_50s_subunits["subunitIds"].tolist() 
			+ ribosome_30s_subunits["subunitIds"].tolist())
		unique_molecule_ids = unique_molecule_counts.readAttribute("objectNames")
		ribosome_idx = unique_molecule_ids.index("activeRibosome")
		n_active_ribosome = unique_molecule_counts.readColumn("uniqueMoleculeCounts")[:, ribosome_idx]
		n_active_ribosome_stoich = np.tensordot(n_active_ribosome, ribosome_stoich, axes = 0)

		# Include the rRNAs already in active ribosomes & ribosome subunits
		for i in range(len(ribosome_subunit_ids)):
			rRna = ribosome_subunit_ids[i]
			if (rRna in rrna_ids_16s):
				rrna_index = rrna_ids_16s.index(rRna)
				free_rrna_counts_16s[:, rrna_index] = (free_rrna_counts_16s[:, rrna_index] 
					+ n_active_ribosome_stoich[:, i])
				free_rrna_counts_16s[:, rrna_index] = (free_rrna_counts_16s[:, rrna_index] 
					+ n_subunit_ribosome_stoich[:, i])
			elif (rRna in rrna_ids_23s):
				rrna_index = rrna_ids_23s.index(rRna)
				free_rrna_counts_23s[:, rrna_index] = (free_rrna_counts_23s[:, rrna_index] 
					+ n_active_ribosome_stoich[:, i])
				free_rrna_counts_23s[:, rrna_index] = (free_rrna_counts_23s[:, rrna_index] 
					+ n_subunit_ribosome_stoich[:, i])
			elif (rRna in rrna_ids_5s):
				rrna_index = rrna_ids_5s.index(rRna)
				free_rrna_counts_5s[:, rrna_index] = (free_rrna_counts_5s[:, rrna_index] 
					+ n_active_ribosome_stoich[:, i])
				free_rrna_counts_5s[:, rrna_index] = (free_rrna_counts_5s[:, rrna_index] 
					+ n_subunit_ribosome_stoich[:, i])

		# Load molecular weights
		rrna_16s_mw = sim_data.getter.getMass(rrna_ids_16s).asNumber(units.g/units.mol)
		rrna_23s_mw = sim_data.getter.getMass(rrna_ids_23s).asNumber(units.g/units.mol)
		rrna_5s_mw = sim_data.getter.getMass(rrna_ids_5s).asNumber(units.g/units.mol)
		nAvogadro = sim_data.constants.nAvogadro.asNumber(1/units.mol)

		# Convert to mass
		rRna_16s = np.dot(free_rrna_counts_16s, rrna_16s_mw)/nAvogadro*(10**15)
		rRna_23s = np.dot(free_rrna_counts_23s, rrna_23s_mw)/nAvogadro*(10**15)
		rRna_5s = np.dot(free_rrna_counts_5s, rrna_5s_mw)/nAvogadro*(10**15)

		# lipids
		lipid_ids = sim_data.moleculeGroups.lipids
		lipid_indexes = np.array([bulk_molecule_ids.index(lipid) for lipid in lipid_ids], np.int64)
		lipid_counts = bulk_molecule_counts[:, lipid_indexes]
		lipids_mw = sim_data.getter.getMass(lipid_ids).asNumber(units.g/units.mol)
		lipid = np.dot(lipid_counts, lipids_mw)/nAvogadro*(10**15)

		# LPS
		lps_id = sim_data.moleculeIds.LPS
		lps_index = np.int64(bulk_molecule_ids.index(lps_id))
		lps_counts = bulk_molecule_counts[:, lps_index]
		lps_mw = sim_data.getter.getMass([lps_id]).asNumber(units.g/units.mol)
		lps = lps_counts*lps_mw/nAvogadro*(10**15)

		# peptidoglycans
		murein_id = sim_data.moleculeIds.murein
		murein_index = np.int64(bulk_molecule_ids.index(murein_id))
		murein_counts = bulk_molecule_counts[:, murein_index]
		murein_mw = sim_data.getter.getMass([murein_id]).asNumber(units.g/units.mol)
		murein = murein_counts*murein_mw/nAvogadro*(10**15)

		# polyamines
		polyamines_ids = sim_data.moleculeGroups.polyamines
		polyamines_indexes = np.array(
			[bulk_molecule_ids.index(polyamine) for polyamine in polyamines_ids],
			np.int64)
		polyamines_counts = bulk_molecule_counts[:, polyamines_indexes]
		polyamines_mw = sim_data.getter.getMass(polyamines_ids).asNumber(units.g/units.mol)
		polyamines = np.dot(polyamines_counts, polyamines_mw)/nAvogadro*(10**15)

		# glycogen
		glycogen_id = sim_data.moleculeIds.glycogen
		glycogen_index = np.int64(bulk_molecule_ids.index(glycogen_id))
		glycogen_counts = bulk_molecule_counts[:, glycogen_index]
		glycogen_mw = sim_data.getter.getMass([glycogen_id]).asNumber(units.g/units.mol)
		glycogen = glycogen_counts*glycogen_mw/nAvogadro*(10**15)

		# other cell components
		protein = mass.readColumn("proteinMass")
		rna = mass.readColumn("rnaMass")
		tRna = mass.readColumn("tRnaMass")
		rRna = mass.readColumn("rRnaMass")
		mRna = mass.readColumn("mRnaMass")
		miscRna = rna - (tRna + rRna + mRna)
		dna = mass.readColumn("dnaMass")
		smallMolecules = mass.readColumn("smallMoleculeMass")
		metabolites = smallMolecules - (lipid + lps + murein + polyamines + glycogen)

		# create dataframe
		mass_final = np.array([protein[-1], 
			dna[-1], mRna[-1], tRna[-1], rRna_16s[-1], rRna_23s[-1], rRna_5s[-1], miscRna[-1], 
			murein[-1], lps[-1], lipid[-1], polyamines[-1], glycogen[-1], metabolites[-1]])
		ids = ['protein',
			'DNA', 'mRNA', 'tRNA', '16srRNA', '23srRNA', '5srRNA', 'miscRNA',
			'peptidoglycan', 'LPS', 'lipid', 'polyamines', 'glycogen', 'metabolites']
		index1 = ['protein',
			'DNA', 'mRNA', 'tRNA', 'rRNA', 'rRNA', 'rRNA', 'miscRNA',
			'peptidoglycan', 'LPS', 'lipid', 'polyamines', 'glycogen', 'metabolites']
		index0 = ['protein']*1 + ['nucleic_acid']*7 + ['metabolites']*6
		df = pd.DataFrame(mass_final.reshape((-1,1)), index = [index0, index1, ids], columns = ['values'])
		df = df.sort_index()

		# create the plot (layered)
		voronoi_0, voronoi_1_all, voronoi_2_all = VM.layered_voronoi(df, CANVAS, I_MAX, ERR_THRES)
		error_all = VM.layered_voronoi_plot(voronoi_0, voronoi_1_all, voronoi_2_all)

		# save the plot
		plt.title("Biomass components")
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close("all")
		
if __name__ == "__main__":
	Plot().cli()
	