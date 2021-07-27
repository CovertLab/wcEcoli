from __future__ import absolute_import, division, print_function

import os
from six.moves import cPickle

import numpy as np
from matplotlib import pyplot as plt

from wholecell.io.tablereader import TableReader
from wholecell.analysis.analysis_tools import exportFigure
from models.ecoli.analysis import singleAnalysisPlot
from wholecell.analysis.analysis_tools import read_bulk_molecule_counts


class Plot(singleAnalysisPlot.SingleAnalysisPlot):
	def do_plot(self, simOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		
		# Amino acid IDs
		sim_data = cPickle.load(open(simDataFile, "rb"))
		aaIDs = sim_data.molecule_groups.amino_acids
		
		# Amino acid exchanges fluxes
		main_reader = TableReader(os.path.join(simOutDir, "Main"))
		initialTime = main_reader.readAttribute("initialTime")
		time = main_reader.readColumn("time") - initialTime
		fba_results = TableReader(os.path.join(simOutDir, "FBAResults"))
		externalExchangeFluxes = fba_results.readColumn("externalExchangeFluxes")
		externalMoleculeIDs = fba_results.readAttribute("externalMoleculeIDs")

		#Transporter counts
		monomers = TableReader(os.path.join(simOutDir,"BulkMolecules")) 
		monomerNames=monomers.readAttribute('objectNames')
		monomerCounts=monomers.readColumn('counts')

		# kcats to aa mapping
		kcats={aa:kcat for aa,kcat in zip(aaIDs,sim_data.process.metabolism.uptake_kcats_per_aa)}

		# Cell Dry mas at each time step
		mass=TableReader(os.path.join(simOutDir, "Mass"))
		dry_mass=mass.readColumn('dryMass')

		# Map aa to its transporters
		aa_to_transporters={
			"L-ALPHA-ALANINE[c]":["CYCA-MONOMER[i]","EG12713-MONOMER[i]","G7399-MONOMER[i]"],
			"ARG[c]":["YGGA-MONOMER[i]","EG12713-MONOMER[i]","ABC-4-CPLX[i]","CPLX0-7535[i]"],
			"ASN[c]":["ANSP-MONOMER[i]","EG12713-MONOMER[i]"],
			"L-ASPARTATE[c]":["DCUA-MONOMER[i]","EG12713-MONOMER[i]","ABC-13-CPLX[i]","GLTP-MONOMER[i]","YCHM-MONOMER[i]","DCTA-MONOMER[i]"],
			"CYS[c]":[],
			"GLT[c]":["XASA-MONOMER[i]","EG12713-MONOMER[i]","ABC-13-CPLX[i]","GLTS-MONOMER[i]", "GLTP-MONOMER[i]"],
			"GLN[c]":["ABC-12-CPLX[i]", "EG12713-MONOMER[i]"], 
			"GLY[c]":["CYCA-MONOMER[i]","EG12713-MONOMER[i]","CPLX0-7654[i]"], 
			"HIS[c]":["ABC-14-CPLX[i]","EG12713-MONOMER[i]"], 
			"ILE[c]":["BRNQ-MONOMER[i]","EG12713-MONOMER[i]","B4141-MONOMER[i]","ABC-15-CPLX[i]"], 
			"LEU[c]":["BRNQ-MONOMER[i]","EG12713-MONOMER[i]","ABC-15-CPLX[i]","ABC-304-CPLX[i]", "B4141-MONOMER[i]","G6984-MONOMER[i]"], 
			"LYS[c]":["LYSP-MONOMER[i]","EG12713-MONOMER[i]","CADB-MONOMER[i]","G6458-MONOMER[i]","ABC-3-CPLX[i]"], 
			"MET[c]":["METNIQ-METHIONINE-ABC-CPLX[o]","EG12713-MONOMER[i]", "B4141-MONOMER[i]","G6984-MONOMER[i]"], 
			"PHE[c]":["EG12713-MONOMER[i]", "PHEP-MONOMER[i]","AROP-MONOMER[i]", "ABC-15-CPLX[i]","ABC-304-CPLX[i]"], 
			"PRO[c]":["ABC-26-CPLX[i]","PUTP-MONOMER[i]","CPLX0-7642[m]","EG12713-MONOMER[i]"], 
			"SER[c]":["YGJU-MONOMER[m]","EG12713-MONOMER[i]", "SDAC-MONOMER[i]","TDCC-MONOMER[i]"], 
			"THR[c]":["TDCC-MONOMER[i]","EG12713-MONOMER[i]", "RHTB-MONOMER[i]","RHTC-MONOMER[m]","EG12134-MONOMER[i]"], 
			"TRP[c]":["EG12713-MONOMER[i]", "TNAB-MONOMER[i]","AROP-MONOMER[i]","MTR-MONOMER[i]"], 
			"TYR[c]":["EG12713-MONOMER[i]", "TYRP-MONOMER[i]","AROP-MONOMER[i]"], 
			"L-SELENOCYSTEINE[c]":["EG12713-MONOMER[i]"], 
			"VAL[c]":["BRNQ-MONOMER[i]","EG12713-MONOMER[i]","B4141-MONOMER[i]","CPLX0-7684[i]","ABC-15-CPLX[i]"]
		}

		# Plot 1
		rows = 6
		cols = 4
		fig = plt.figure(figsize = (8, 11.5))

		for plotIndex, aa in enumerate(aaIDs):
			ax = plt.subplot(rows, cols, plotIndex + 1)

			# Get actual fluxes
			if not aa.startswith("L-SELENOCYSTEINE"):
				old_aa=aa
				aa = aa[:-3] + "[p]"
			if aa in externalMoleculeIDs:
				aaFlux = externalExchangeFluxes[:, externalMoleculeIDs.index(aa)]
			else:
				aaFlux = np.zeros(len(time))

			#calculate expected rate based on kcats and transporters
			rate_= [kcats[old_aa]]*len(aaFlux)
			t_counts=0
			for i in aa_to_transporters[old_aa]:
				t_counts+=monomerCounts[:, monomerNames.index(i)]
			
			# Right now, code breaks because CYS has 0 transporters
			if 'CYS' not in aa:
				rate_ *= t_counts

			#rate in mol/g/hr
			rate_*= (-3600*1000)/(dry_mass*1e-15*sim_data.constants.n_avogadro.asNumber())

			#Plot, orange is target flux and blue is actual flux
			ax.plot(time / 60., aaFlux, linewidth = 1, label = 'Actual Flux')
			ax.plot(time / 60., rate_,  linewidth = 1, label = 'Expected Flux')
			ax.set_xlabel("Time (min)", fontsize = 6)
			ax.set_ylabel("mol/gDCW/hr", fontsize = 6)
			ax.set_title("%s" % aa, fontsize = 6, y = 1.1)
			ax.tick_params(which = "both", direction = "out", labelsize = 6)

		plt.rc("font", size = 6)
		plt.suptitle("External exchange fluxes of amino acids", fontsize = 10)

		plt.subplots_adjust(hspace = 1, wspace = 1)
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close("all")



		# Plot 2 - Transporters counts
		rows = 8
		cols = 6
		fig = plt.figure(figsize = (8, 11.5))

		mmIDs = ["CYCA-MONOMER[i]","G7399-MONOMER[i]", "YGGA-MONOMER[i]","ABC-4-CPLX[i]",
					"CPLX0-7535[i]", "ANSP-MONOMER[i]", "DCUA-MONOMER[i]","GLTP-MONOMER[i]","YCHM-MONOMER[i]","DCTA-MONOMER[i]",
					"ABC-13-CPLX[i]","EG11902-MONOMER[m]", "EG11639-MONOMER[i]", "G6934-MONOMER[i]", "CPLX0-8152[i]", "ABC-6-CPLX[i]",
					"XASA-MONOMER[i]","GLTS-MONOMER[i]","ABC-12-CPLX[i]","BRNQ-MONOMER[i]","CPLX0-7654[i]", "ABC-14-CPLX[i]", 
					"B4141-MONOMER[i]","LYSP-MONOMER[i]","CADB-MONOMER[i]","G6458-MONOMER[i]","ABC-3-CPLX[i]", "METNIQ-METHIONINE-ABC-CPLX[o]",
					"G6984-MONOMER[i]", "PHEP-MONOMER[i]", "ABC-15-CPLX[i]","ABC-304-CPLX[i]", "ABC-26-CPLX[i]","PUTP-MONOMER[i]",
					"CPLX0-7642[m]", "YGJU-MONOMER[m]","SDAC-MONOMER[i]", "TDCC-MONOMER[i]","RHTB-MONOMER[i]","RHTC-MONOMER[m]",
					"EG12134-MONOMER[i]", "TNAB-MONOMER[i]","MTR-MONOMER[i]", "TYRP-MONOMER[i]","AROP-MONOMER[i]","EG12713-MONOMER[i]",
					"CPLX0-7684[i]"]

		for plotIndex, m in enumerate(mmIDs):
			ax = plt.subplot(rows, cols, plotIndex + 1)

			if m in monomerNames:
				mC = monomerCounts[:, monomerNames.index(m)]
			else:
				mC = np.zeros(len(time))

			ax.plot(time / 60., mC)
			ax.set_xlabel("Time (min)", fontsize = 6)
			ax.set_ylabel("# counts", fontsize = 6)
			ax.set_title("%s" % m, fontsize = 6, y = 1.1)
			ax.tick_params(which = "both", direction = "out", labelsize = 6)

		plt.rc("font", size = 6)
		plt.suptitle("Counts of monomers of interest", fontsize = 10)

		plt.subplots_adjust(hspace = 1, wspace = 1)
		exportFigure(plt, plotOutDir, 'monomers_'+plotOutFileName, metadata)
		plt.close("all")


if __name__ == "__main__":
	Plot().cli()