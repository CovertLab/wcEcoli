"""
Plots the genomic locations of collisions between replisomes and RNAPs using
polar coordinates.
"""

import os
import pickle

from matplotlib import pyplot as plt
import numpy as np

from models.ecoli.analysis import singleAnalysisPlot
from wholecell.analysis.analysis_tools import exportFigure
from wholecell.io.tablereader import TableReader
from wholecell.utils import units


class Plot(singleAnalysisPlot.SingleAnalysisPlot):
	def do_plot(self, simOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		with open(simDataFile, 'rb') as f:
			sim_data = pickle.load(f)

		# Load replichore lengths
		replichore_lengths = sim_data.process.replication.replichore_lengths

		# Load coordinates of rRNA transcription units
		rna_data = sim_data.process.transcription.rna_data
		rRNA_data = rna_data[rna_data['is_rRNA']]
		rRNA_start_coordinates = rRNA_data['replication_coordinate']
		rRNA_lengths = rRNA_data['length'].asNumber(units.nt)
		rRNA_is_forward = rRNA_data['is_forward']
		rRNA_end_coordinates = (
			rRNA_start_coordinates
			+ 2*(rRNA_is_forward.astype(int) - 0.5) * rRNA_lengths)

		# Listeners used
		rnap_data_reader = TableReader(os.path.join(simOutDir, "RnapData"))

		# Load data
		headon_coordinates = rnap_data_reader.readColumn(
			"headon_collision_coordinates")
		codirectional_coordinates = rnap_data_reader.readColumn(
			"codirectional_collision_coordinates")

		def convert_to_theta(coordinates):
			"""
			Converts genomic coordinates to values of theta in polar
			coordinates.
			"""
			# Flatten and remove NaN's
			coordinates_flat = coordinates.flatten()
			genomic_coordinates = coordinates_flat[~np.isnan(coordinates_flat)]

			# Convert to theta (assuming oriC is at np.pi/2)
			thetas = -np.pi*(genomic_coordinates / replichore_lengths[(genomic_coordinates < 0).astype(int)])

			return thetas

		headon_thetas = convert_to_theta(headon_coordinates)
		codirectional_thetas = convert_to_theta(codirectional_coordinates)
		rRNA_start_thetas = convert_to_theta(rRNA_start_coordinates)
		rRNA_end_thetas = convert_to_theta(rRNA_end_coordinates)

		# Plot
		plt.subplots(subplot_kw={'projection': 'polar'})
		ax = plt.subplot(1, 1, 1)

		# Remove grid and labels
		ax.set_theta_offset(np.pi / 2)  # Put 0 deg at the top instead of the right side
		ax.grid(False)
		ax.set_yticklabels([])
		ax.axis("off")

		# Plot collision locations and genome outline
		ax.vlines(headon_thetas, 1, 1.05, color='C3', lw=0.5)
		ax.vlines(codirectional_thetas, 0.95, 1, color='#555555', lw=0.5)
		ax.plot(np.linspace(0, 2 * np.pi, 1000), np.repeat(1, 1000), color='k')
		ax.vlines(0, 0.98, 1.02, color='k')
		ax.vlines(np.pi, 0.98, 1.02, color='k')

		# Plot locations of rRNA transcripts
		for (s, e) in zip(rRNA_start_thetas, rRNA_end_thetas):
			ax.plot(
				np.linspace(s, e, 10), np.repeat(0.9, 10), color='k')

		# Add labels
		ax.text(0, 1.08, "oriC", ha="center", va="center")
		ax.text(np.pi, 1.08, "terC", ha="center", va="center")

		plt.tight_layout()
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close('all')


if __name__ == "__main__":
	Plot().cli()
