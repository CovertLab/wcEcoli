"""
Plots limited enzyme fluxes, protein counts, and transcription initiation events.
"""

import os

import numpy as np
import matplotlib.pyplot as plt

from wholecell.io.tablereader import TableReader
from models.ecoli.processes.metabolism import COUNTS_UNITS, TIME_UNITS, VOLUME_UNITS
from wholecell.analysis.analysis_tools import exportFigure
from models.ecoli.analysis import multigenAnalysisPlot


class Plot(multigenAnalysisPlot.MultigenAnalysisPlot):
	def do_plot(self, seedOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		enzymeMonomerId = "GLUTCYSLIG-MONOMER[c]"
		enzyme_rna_cistron_id = "EG10418_RNA"
		reactionId = "GLUTCYSLIG-RXN"
		transcriptionFreq = 1.0
		metaboliteId = "GLUTATHIONE[c]"

		# Get all cells
		allDir = self.ap.get_cells()

		simOutDir = os.path.join(allDir[0], "simOut")
		bulkMolecules = TableReader(os.path.join(simOutDir, "BulkMolecules"))
		moleculeIds = bulkMolecules.readAttribute("objectNames")
		enzymeMonomerIndex = moleculeIds.index(enzymeMonomerId)
		metaboliteIndex = moleculeIds.index(metaboliteId)

		RNA_counts_reader = TableReader(
			os.path.join(simOutDir, 'RNACounts'))
		all_mRNA_ids = RNA_counts_reader.readAttribute('mRNA_cistron_ids')
		enzyme_rna_cistron_index = all_mRNA_ids.index(enzyme_rna_cistron_id)

		# TODO (ggsun): Should be tweaked with operons
		rnapDataReader = TableReader(os.path.join(simOutDir, "RnapData"))
		rnap_data_cistron_ids = rnapDataReader.readAttribute('cistron_ids')
		enzyme_RNA_cistron_index_rnap_data = rnap_data_cistron_ids.index(enzyme_rna_cistron_id)

		time = []
		enzymeFluxes = []
		enzymeMonomerCounts = []
		enzymeRnaCounts = []
		enzyme_cistron_init_event = []
		metaboliteCounts = []

		for gen, simDir in enumerate(allDir):
			simOutDir = os.path.join(simDir, "simOut")

			time += TableReader(os.path.join(simOutDir, "Main")).readColumn("time").tolist()

			bulkMolecules = TableReader(os.path.join(simOutDir, "BulkMolecules"))
			moleculeCounts = bulkMolecules.readColumn("counts")
			enzymeMonomerCounts += moleculeCounts[:, enzymeMonomerIndex].tolist()
			metaboliteCounts += moleculeCounts[:, metaboliteIndex].tolist()

			RNA_counts_reader = TableReader(
				os.path.join(simOutDir, 'RNACounts'))
			mRNA_cistron_counts = RNA_counts_reader.readColumn('mRNA_cistron_counts')
			enzymeRnaCounts += mRNA_cistron_counts[:, enzyme_rna_cistron_index].tolist()

			fbaResults = TableReader(os.path.join(simOutDir, "FBAResults"))
			reactionIDs = np.array(fbaResults.readAttribute("reactionIDs"))
			reactionFluxes = np.array(fbaResults.readColumn("reactionFluxes"))
			enzymeFluxes += reactionFluxes[:, np.where(reactionIDs == reactionId)[0][0]].tolist()

			rnapDataReader = TableReader(os.path.join(simOutDir, "RnapData"))
			enzyme_cistron_init_event += rnapDataReader.readColumn("rna_init_event_per_cistron")[:, enzyme_RNA_cistron_index_rnap_data].tolist()

		time = np.array(time)


		# Plot
		plt.figure(figsize = (10, 10))
		rnaInitAxis = plt.subplot(5, 1, 1)
		rnaAxis = plt.subplot(5, 1, 2, sharex = rnaInitAxis)
		monomerAxis = plt.subplot(5, 1, 3, sharex = rnaInitAxis)
		fluxAxis = plt.subplot(5, 1, 4, sharex = rnaInitAxis)
		metAxis = plt.subplot(5, 1, 5, sharex = rnaInitAxis)

		rnaInitAxis.plot(time / 3600., enzyme_cistron_init_event)
		rnaInitAxis.set_title("%s transcription initiation events" % enzyme_rna_cistron_id, fontsize = 10)
		rnaInitAxis.set_ylim([0, rnaInitAxis.get_ylim()[1] * 1.1])
		rnaInitAxis.set_xlim([0, time[-1] / 3600.])

		rnaAxis.plot(time / 3600., enzymeRnaCounts)
		rnaAxis.set_title("%s counts" % enzyme_rna_cistron_id, fontsize = 10)

		monomerAxis.plot(time / 3600., enzymeMonomerCounts)
		monomerAxis.set_title("%s counts" % enzymeMonomerId, fontsize = 10)

		fluxAxis.plot(time / 3600., enzymeFluxes)
		fluxAxis.set_title("%s flux (%s / %s / %s)" % (reactionId, COUNTS_UNITS, VOLUME_UNITS, TIME_UNITS), fontsize = 10)

		metAxis.plot(time / 3600., metaboliteCounts)
		metAxis.set_title("%s counts" % metaboliteId, fontsize = 10)
		metAxis.set_xlabel("Time (hour)\n(%s frequency of at least 1 transcription per generation)" % transcriptionFreq, fontsize = 10)

		plt.subplots_adjust(wspace = 0.4, hspace = 0.4) #, right = 0.83, bottom = 0.05, left = 0.07, top = 0.95)
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close("all")


if __name__ == "__main__":
	Plot().cli()
