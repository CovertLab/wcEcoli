"""
Plots the locations and structures of rRNA operons for a particular sim_data
object.
"""

import pickle

from matplotlib import pyplot as plt
# noinspection PyUnresolvedReferences
import numpy as np

from models.ecoli.analysis import parcaAnalysisPlot
from wholecell.analysis.analysis_tools import exportFigure
from wholecell.utils import units


TU_ID_TO_RRNA_OPERON_ID = {
	'TU0-1181[c]': 'rrnA',
	'TU0-1182[c]': 'rrnB',
	'TU0-1183[c]': 'rrnC',
	'TU0-1191[c]': 'rrnD',
	'TU0-1186[c]': 'rrnE',
	'TU0-1187[c]': 'rrnG',
	'TU0-1189[c]': 'rrnH',
	}


class Plot(parcaAnalysisPlot.ParcaAnalysisPlot):
	def do_plot(self, input_dir, plot_out_dir, plot_out_filename, sim_data_file, validation_data_file, metadata):
		with open(sim_data_file, 'rb') as f:
			sim_data = pickle.load(f)

		# Load replichore lengths
		replichore_lengths = sim_data.process.replication.replichore_lengths

		# Load attributes of rRNA transcription units
		rna_data = sim_data.process.transcription.rna_data
		rRNA_data = rna_data[rna_data['is_rRNA']]
		rRNA_ids = rRNA_data['id']
		rRNA_id_to_index = {rRNA_id: i for (i, rRNA_id) in enumerate(rRNA_ids)}
		rRNA_start_coordinates = rRNA_data['replication_coordinate']
		rRNA_lengths = rRNA_data['length'].asNumber(units.nt)
		rRNA_is_forward = rRNA_data['is_forward']
		rRNA_end_coordinates = (
			rRNA_start_coordinates
			+ 2 * (rRNA_is_forward.astype(int) - 0.5) * rRNA_lengths)

		# Load attributes of genes that are part of rRNA transcription units
		cistron_data = sim_data.process.transcription.cistron_data
		gene_start_coordinates = []
		gene_end_coordinates = []
		gene_names = []
		gene_is_rRNA = []

		for rRNA_id in rRNA_ids:
			# Extract data for cistrons belonging to this rRNA operon
			rRNA_index = rRNA_id_to_index[rRNA_id]
			cistron_indexes = sim_data.process.transcription.rna_id_to_cistron_indexes(
				rRNA_id)
			gene_ids = cistron_data['gene_id'][cistron_indexes]
			cistron_start_coordinates = cistron_data['replication_coordinate'][cistron_indexes]
			cistron_lengths = cistron_data['length'][cistron_indexes].asNumber(units.nt)
			cistron_is_forward = cistron_data['is_forward'][cistron_indexes]
			cistron_is_rRNA = cistron_data['is_rRNA'][cistron_indexes]
			cistron_end_coordinates = (
				cistron_start_coordinates
				+ 2 * (cistron_is_forward.astype(int) - 0.5) * cistron_lengths)

			# Append data to list
			gene_start_coordinates.append(
				np.abs(cistron_start_coordinates - rRNA_start_coordinates[rRNA_index])
				)
			gene_end_coordinates.append(
				np.abs(cistron_end_coordinates - rRNA_start_coordinates[rRNA_index])
				)
			gene_names.append(
				[sim_data.common_names.get_common_name(x) for x in gene_ids]
				)
			gene_is_rRNA.append(cistron_is_rRNA)

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

		rRNA_start_thetas = convert_to_theta(rRNA_start_coordinates)
		rRNA_end_thetas = convert_to_theta(rRNA_end_coordinates)

		plt.figure(figsize=(9, 4))

		# Plot locations of rRNA operons on circular chromosome
		ax0 = plt.subplot(1, 2, 1, projection='polar')

		# Remove grid and labels
		ax0.set_theta_offset(np.pi / 2)  # Put 0 deg at the top instead of the right side
		ax0.grid(False)
		ax0.set_yticklabels([])
		ax0.set_ylim([0, 1.5])
		ax0.axis("off")

		# Plot collision locations and genome outline
		ax0.plot(
			np.linspace(0, 2 * np.pi, 1000), np.repeat(1, 1000), color='k',
			lw=0.5)
		ax0.vlines(0, 0.98, 1.02, color='k', lw=0.5)
		ax0.vlines(np.pi, 0.98, 1.02, color='k', lw=0.5)

		# Plot locations of rRNA transcripts
		for (s, e) in zip(rRNA_start_thetas, rRNA_end_thetas):
			ax0.plot(
				np.linspace(s, e, 10), np.repeat(1, 10), color='k', lw=2)

		# Add labels
		ax0.text(0, 0.87, "oriC", ha="center", va="center")
		ax0.text(np.pi, 0.9, "terC", ha="center", va="center")

		# Plot locations of genes within each rRNA operon
		ax1 = plt.subplot(1, 2, 2)

		# Remove grid and labels
		ax1.grid(False)
		ax1.set_yticklabels([])
		ax1.set_ylim([-2, 8])
		ax1.axis("off")

		# Plot locations of genes within each operon
		for i, (rRNA_id, operon_id) in enumerate(TU_ID_TO_RRNA_OPERON_ID.items()):
			rRNA_index = rRNA_id_to_index[rRNA_id]

			# Plot lines for template DNA and boundaries of operons
			ax1.hlines(6 - i, -500, 6500, color='k', lw=1)
			ax1.hlines(6 - i, 0, rRNA_lengths[rRNA_index], color='k', lw=2)

			# Plot lines for the boundaries of genes
			for (sc, ec, name, is_rRNA) in zip(
					gene_start_coordinates[rRNA_index],
					gene_end_coordinates[rRNA_index],
					gene_names[rRNA_index],
					gene_is_rRNA[rRNA_index]):
				if is_rRNA:
					ax1.hlines(6 - i, sc, ec, color='k', lw=5)

					# Position adjustments needed to prevent overlapping text
					adj = 0
					if name == 'rrfD':
						adj = -175
					elif name == 'rrfF':
						adj = 175

					# Add labels for gene names
					ax1.text((sc + ec) / 2 + adj, 6.4 - i, name, ha="center", va="center")

				# Use lighter color for non-rRNA genes
				else:
					ax1.hlines(6 - i, sc, ec, color='#bbbbbb', lw=5)

			ax1.text(-1000, 6 - i, operon_id, ha="center", va="center")

		exportFigure(plt, plot_out_dir, plot_out_filename, metadata)
		plt.close('all')


if __name__ == "__main__":
	Plot().cli()
