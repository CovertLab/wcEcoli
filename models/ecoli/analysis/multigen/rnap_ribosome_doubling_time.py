from __future__ import absolute_import
from __future__ import division

import os

import numpy as np
from matplotlib import pyplot as plt

from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from wholecell.io.tablereader import TableReader
from wholecell.analysis.analysis_tools import exportFigure
from wholecell.utils import units
import cPickle
from models.ecoli.analysis import multigenAnalysisPlot


class Plot(multigenAnalysisPlot.MultigenAnalysisPlot):
	def do_plot(self, seedOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		if not os.path.isdir(seedOutDir):
			raise Exception, "seedOutDir does not currently exist as a directory"

		if not os.path.exists(plotOutDir):
			os.mkdir(plotOutDir)

		ap = AnalysisPaths(seedOutDir, multi_gen_plot = True)
		sim_data = cPickle.load(open(simDataFile, "rb"))
		rnap_id = sim_data.moleculeIds.rnapFull
		rnap_subunit_ids = sim_data.moleculeGroups.rnapIds
		ribosome_30s_id = sim_data.moleculeIds.s30_fullComplex
		ribosome_50s_id = sim_data.moleculeIds.s50_fullComplex
		rprotein_ids = sim_data.moleculeGroups.rProteins

		fig, axes_list = plt.subplots(3, 1, figsize=(8.5, 11))
		ax0, ax1, ax2 = axes_list

		for gen, sim_dir in enumerate(ap.get_cells()):
			sim_out_dir = os.path.join(sim_dir, "simOut")

			unique_molecules_reader = TableReader(os.path.join(sim_out_dir, "UniqueMoleculeCounts"))
			rnap_index = unique_molecules_reader.readAttribute("uniqueMoleculeIds").index("activeRnaPoly")
			ribosome_index = unique_molecules_reader.readAttribute("uniqueMoleculeIds").index("activeRibosome")

			bulk_molecules_reader = TableReader(os.path.join(sim_out_dir, "BulkMolecules"))
			molecule_ids = bulk_molecules_reader.readAttribute("objectNames")
			rnap_full_index = molecule_ids.index(rnap_id)
			ribosome_30s_index = molecule_ids.index(ribosome_30s_id)
			ribosome_50s_index = molecule_ids.index(ribosome_50s_id)

			molecule_counts = bulk_molecules_reader.readColumn("counts")


			# RNA synthesis
			# todo

			# RNA polymerase
			n_rnap_active = unique_molecules_reader.readColumn("uniqueMoleculeCounts")[:, rnap_index]
			n_rnap_total = n_rnap_active + molecule_counts[:, rnap_full_index]
			# n_rnap_terminated = TableReader(os.path.join(sim_out_dir, "RnapData"))
			ax0.scatter(gen, n_rnap_total[0], color="tab:blue")
			ax0.scatter(gen, n_rnap_active[0], color="tab:orange")

			# Ribosome
			n_ribosome_active = unique_molecules_reader.readColumn("uniqueMoleculeCounts")[:, ribosome_index]
			n_ribosome_30s = molecule_counts[:, ribosome_30s_index]
			n_ribosome_50s = molecule_counts[:, ribosome_50s_index]
			n_ribosome_total = n_ribosome_active + np.minimum(n_ribosome_30s, n_ribosome_50s)
			ax1.scatter(gen, n_ribosome_total[0], color="tab:blue")
			ax1.scatter(gen, n_ribosome_active[0], color="tab:orange")

			# Close readers
			unique_molecules_reader.close()
			bulk_molecules_reader.close()

			# Doubling time
			time = TableReader(os.path.join(sim_out_dir, "Main")).readColumn("time")
			time_doubling = time[-1] - time[0]
			ax2.scatter(gen, time_doubling / 60.)

		ax0.set_ylabel("RNA polymerase (count)")
		ax0.axhline(760, color="tab:orange")

		ax1.set_ylabel("Ribosome (count)")
		ax1.axhline(12600, color="tab:orange")
		ax2.set_ylabel("Doubling time (min)")
		ax2.set_xlabel("Generation")
		exportFigure(plt, plotOutDir, plotOutFileName,metadata)
		plt.close("all")


if __name__ == "__main__":
	Plot().cli()
