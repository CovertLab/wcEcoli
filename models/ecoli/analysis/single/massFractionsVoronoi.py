"""
Plot the Voronoi diagram of mass fractions

@author: Ray Chang
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 09/27/2019
"""
from __future__ import absolute_import

import os
import cPickle
import numpy as np
from matplotlib import pyplot as plt

from wholecell.io.tablereader import TableReader
from wholecell.analysis.analysis_tools import exportFigure
from models.ecoli.analysis import singleAnalysisPlot
from wholecell.utils import units
import sys
sys.path.append('models/ecoli/analysis/single')
from voronoiPlotMain import *

class Plot(singleAnalysisPlot.SingleAnalysisPlot):
	def do_plot(self, simOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		if not os.path.isdir(simOutDir):
			raise Exception, "simOutDir does not currently exist as a directory"
		if not os.path.exists(plotOutDir):
			os.mkdir(plotOutDir)
		with open(simDataFile, 'rb') as f:
			sim_data = cPickle.load(f)
		# Random seed
		np.random.seed(0)
		# Load data
		RibosomeData = TableReader(os.path.join(simOutDir, "RibosomeData"))
		mass = TableReader(os.path.join(simOutDir, "Mass"))
		UniqueMoleculeCounts = TableReader(os.path.join(simOutDir, "UniqueMoleculeCounts"))
		bulkMolecules = TableReader(os.path.join(simOutDir, "BulkMolecules"))
		bulkMoleculeCounts = bulkMolecules.readColumn("counts")
		bulk_molecule_ids = bulkMolecules.readAttribute("objectNames")
		# Count free rRNAs
		rRnaIds_16s = sim_data.moleculeGroups.s30_16sRRNA
		rRnaIds_23s = sim_data.moleculeGroups.s50_23sRRNA
		rRnaIds_5s = sim_data.moleculeGroups.s50_5sRRNA
		rRnaIndexes_16s = np.array([bulk_molecule_ids.index(rRna) for rRna in rRnaIds_16s], np.int)
		rRnaIndexes_23s = np.array([bulk_molecule_ids.index(rRna) for rRna in rRnaIds_23s], np.int)
		rRnaIndexes_5s = np.array([bulk_molecule_ids.index(rRna) for rRna in rRnaIds_5s], np.int)
		freeRRnaCounts_16s = bulkMoleculeCounts[:, rRnaIndexes_16s]
		freeRRnaCounts_23s = bulkMoleculeCounts[:, rRnaIndexes_23s]
		freeRRnaCounts_5s = bulkMoleculeCounts[:, rRnaIndexes_5s]
		# Get the stoichiometry matrix of 50s & 30s subunits
		ribosome_50s_subunits = sim_data.process.complexation.getMonomers(sim_data.moleculeIds.s50_fullComplex)
		ribosome_30s_subunits = sim_data.process.complexation.getMonomers(sim_data.moleculeIds.s30_fullComplex)
		ribosome_stoich = np.hstack((ribosome_50s_subunits["subunitStoich"], ribosome_30s_subunits["subunitStoich"]))
		# Count 50s & 30s subunits
		complexIds = [sim_data.moleculeIds.s50_fullComplex]
		complexIndexes = np.array([bulk_molecule_ids.index(comp) for comp in complexIds], np.int)
		[complex_50s_Counts] = bulkMoleculeCounts[:, complexIndexes].T
		complexIds = [sim_data.moleculeIds.s30_fullComplex]
		complexIndexes = np.array([bulk_molecule_ids.index(comp) for comp in complexIds], np.int)
		[complex_30s_Counts] = bulkMoleculeCounts[:, complexIndexes].T
		n_50s_stoich = np.tensordot(complex_50s_Counts,ribosome_50s_subunits["subunitStoich"],axes=0)
		n_30s_stoich = np.tensordot(complex_30s_Counts,ribosome_50s_subunits["subunitStoich"],axes=0)
		n_subunit_ribosome_stoich = np.hstack((n_50s_stoich,n_30s_stoich))
		# Count active ribosomes
		ribosome_subunit_ids = (ribosome_50s_subunits["subunitIds"].tolist() + ribosome_30s_subunits["subunitIds"].tolist())
		unique_molecule_ids = UniqueMoleculeCounts.readAttribute("objectNames")
		ribosome_idx = unique_molecule_ids.index("activeRibosome")
		n_active_ribosome = UniqueMoleculeCounts.readColumn("uniqueMoleculeCounts")[:, ribosome_idx]
		n_active_ribosome_stoich = np.tensordot(n_active_ribosome,ribosome_stoich,axes=0)
		# Include the rRNAs already in active ribosomes & ribosome subunits
		for i in range(len(ribosome_subunit_ids)):
			rRna = ribosome_subunit_ids[i]
			if (rRna in rRnaIds_16s):
				rRnaIndex = rRnaIds_16s.index(rRna)
				freeRRnaCounts_16s[:,rRnaIndex] = freeRRnaCounts_16s[:,rRnaIndex] + n_active_ribosome_stoich[:,i]
				freeRRnaCounts_16s[:,rRnaIndex] = freeRRnaCounts_16s[:,rRnaIndex] + n_subunit_ribosome_stoich[:,i]
			elif (rRna in rRnaIds_23s):
				rRnaIndex = rRnaIds_23s.index(rRna)
				freeRRnaCounts_23s[:,rRnaIndex] = freeRRnaCounts_23s[:,rRnaIndex] + n_active_ribosome_stoich[:,i]
				freeRRnaCounts_23s[:,rRnaIndex] = freeRRnaCounts_23s[:,rRnaIndex] + n_subunit_ribosome_stoich[:,i]
			elif (rRna in rRnaIds_5s):
				rRnaIndex = rRnaIds_5s.index(rRna)
				freeRRnaCounts_5s[:,rRnaIndex] = freeRRnaCounts_5s[:,rRnaIndex] + n_active_ribosome_stoich[:,i]
				freeRRnaCounts_5s[:,rRnaIndex] = freeRRnaCounts_5s[:,rRnaIndex] + n_subunit_ribosome_stoich[:,i]
		# Load molecular weights
		rRna16sMW = sim_data.getter.getMass(rRnaIds_16s).asNumber(units.g / units.mol)
		rRna23sMW = sim_data.getter.getMass(rRnaIds_23s).asNumber(units.g / units.mol)
		rRna5sMW = sim_data.getter.getMass(rRnaIds_5s).asNumber(units.g / units.mol)
		nAvogadro = sim_data.constants.nAvogadro.asNumber(1 / units.mol)
		# Convert to mass
		rRna_16s = np.dot(freeRRnaCounts_16s, rRna16sMW)/nAvogadro*(10**15)
		rRna_23s = np.dot(freeRRnaCounts_23s, rRna23sMW)/nAvogadro*(10**15)
		rRna_5s = np.dot(freeRRnaCounts_5s, rRna5sMW)/nAvogadro*(10**15)
		# lipids
		lipidIds = sim_data.moleculeGroups.lipids
		lipidIndexes = np.array([bulk_molecule_ids.index(lipid) for lipid in lipidIds], np.int)
		lipidCounts = bulkMoleculeCounts[:, lipidIndexes]
		lipidsMW = sim_data.getter.getMass(lipidIds).asNumber(units.g / units.mol)
		lipid = np.dot(lipidCounts, lipidsMW)/nAvogadro*(10**15)
		# LPS
		lpsIds = sim_data.moleculeGroups.LPS
		lpsIndexes = np.array([bulk_molecule_ids.index(lps) for lps in lpsIds], np.int)
		lpsCounts = bulkMoleculeCounts[:, lpsIndexes]
		lpsMW = sim_data.getter.getMass(lpsIds).asNumber(units.g / units.mol)
		lps = np.dot(lpsCounts, lpsMW)/nAvogadro*(10**15)
		# peptidoglycans
		mureinIds = sim_data.moleculeGroups.murein
		mureinIndexes = np.array([bulk_molecule_ids.index(murein) for murein in mureinIds], np.int)
		mureinCounts = bulkMoleculeCounts[:, mureinIndexes]
		mureinMW = sim_data.getter.getMass(mureinIds).asNumber(units.g / units.mol)
		murein = np.dot(mureinCounts, mureinMW)/nAvogadro*(10**15)
		# polyamines
		polyaminesIds = sim_data.moleculeGroups.polyamines
		polyaminesIndexes = np.array([bulk_molecule_ids.index(polyamine) for polyamine in polyaminesIds], np.int)
		polyaminesCounts = bulkMoleculeCounts[:, polyaminesIndexes]
		polyaminesMW = sim_data.getter.getMass(polyaminesIds).asNumber(units.g / units.mol)
		polyamines = np.dot(polyaminesCounts, polyaminesMW)/nAvogadro*(10**15)
		# glycogen
		glycogenIds = sim_data.moleculeGroups.glycogen
		glycogenIndexes = np.array([bulk_molecule_ids.index(glycogen) for glycogen in glycogenIds], np.int)
		glycogenCounts = bulkMoleculeCounts[:, glycogenIndexes]
		glycogenMW = sim_data.getter.getMass(glycogenIds).asNumber(units.g / units.mol)
		glycogen = np.dot(glycogenCounts, glycogenMW)/nAvogadro*(10**15)
		# other cell components
		protein = mass.readColumn("proteinMass")
		rna = mass.readColumn("rnaMass")
		tRna = mass.readColumn("tRnaMass")
		rRna = mass.readColumn("rRnaMass")
		mRna = mass.readColumn("mRnaMass")
		miscRna = rna - (tRna+rRna+mRna)
		dna = mass.readColumn("dnaMass")
		smallMolecules = mass.readColumn("smallMoleculeMass")
		metabolites = smallMolecules - (lipid+lps+murein+polyamines+glycogen)
		canvas = np.array([[0,0],[4,0],[4,4],[0,4]])
		IDs = np.array(['protein','DNA','mRNA','tRNA','16srRNA','23srRNA','5srRNA','miscRNA','peptidoglycan','LPS','lipid','polyamines','glycogen','metabolites']) 
		massFinal = np.array([protein[-1], dna[-1], mRna[-1], tRna[-1], rRna_16s[-1], rRna_23s[-1], rRna_5s[-1], miscRna[-1], murein[-1], lps[-1], lipid[-1], polyamines[-1], glycogen[-1], metabolites[-1]])
		index1 = np.array(['protein','DNA','mRNA','tRNA','rRNA','rRNA','rRNA','miscRNA','peptidoglycan','LPS','lipid','polyamines','glycogen','metabolites'])
		index0 = np.array(('protein '*1+'NA '*7+'metabolites '*6).split( ))
		df = pd.DataFrame(massFinal.reshape((-1,1)), index = [index0, index1, IDs], columns = ['values'])
		df = df.sort_index()
		i_max = 75
		err_thres = 1E-6
		Voronoi_0, Voronoi_1_all, Voronoi_2_all = Layered_Voronoi(df, canvas, i_max, err_thres)
		error_all = Layered_Voronoi_plot(Voronoi_0, Voronoi_1_all, Voronoi_2_all)
		plt.title("Biomass components")
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close("all")
		
if __name__ == "__main__":
	Plot().cli()