"""
Template for multigen analysis plots
"""

import pickle
import os
from wholecell.utils import units
from matplotlib import pyplot as plt
# noinspection PyUnresolvedReferences
import numpy as np
import plotly.graph_objects as go
from models.ecoli.analysis import multigenAnalysisPlot
from wholecell.analysis.analysis_tools import (exportFigure,
	read_bulk_molecule_counts, read_stacked_bulk_molecules, read_stacked_columns)
from wholecell.io.tablereader import TableReader
import io
from wholecell.io import tsv
from wholecell.utils.filepath import ROOT_PATH

HIGHLIGHT_IN_RED = []#['EG10863-MONOMER[c]','DETHIOBIOTIN-SYN-MONOMER[c]','DCUR-MONOMER[c]']
HIGHLIGHT_IN_BLUE =[]#['CARBPSYN-SMALL[c]', 'CDPDIGLYSYN-MONOMER[i]','EG10743-MONOMER[c]','GLUTCYSLIG-MONOMER[c]']
HIGHLIGHT_IN_PURPLE = ['ADHP-MONOMER[c]', 'G6988-MONOMER[c]','EG11111-MONOMER[c]', 'PD03867[c]', 'EG50004-MONOMER[c]' ]#['G6890-MONOMER[c]','PD03938[c]','G6737-MONOMER[c]','RPOD-MONOMER[c]','PD02936[c]','RED-THIOREDOXIN2-MONOMER[c]']
# lowest deg rates:  ['PD03867[c]', 'EG50004-MONOMER[c]','ADHP-MONOMER[c]', 'G6988-MONOMER[c]']

# todo: highight complexes as squares


# function to match gene symbols to monomer ids
def get_gene_symbols_for_monomer_ids():
	"""
	Extracts the gene symbols for each monomer id in the model.

	Returns: a dictionary mapping monomer ids to gene symbols.
	Code adapted from convert_to_flat.py.
	"""
	RNAS_FILE = os.path.join(ROOT_PATH, 'reconstruction', 'ecoli',
								 'flat', 'rnas.tsv')
	with (io.open(RNAS_FILE, 'rb') as f):
		reader = tsv.reader(f, delimiter='\t')
		headers = next(reader)
		while headers[0].startswith('#'):
			headers = next(reader)

		# extract relevant information
		gene_symbol_index = headers.index('common_name')
		protein_id_index = headers.index('monomer_ids')
		monomer_ids_to_gene_symbols = {}
		for line in reader:
			gene_symbol = line[gene_symbol_index]
			protein_id = list(
				line[protein_id_index][2:-2].split('", "'))[0]
			monomer_ids_to_gene_symbols[protein_id] = gene_symbol

	return monomer_ids_to_gene_symbols

def get_common_name(protein_id):
	"""
    Obtains the common names for each protein of interest
    Args:
        protein_id: the name of the protein(s) of interest

    Returns: the common name for the protein(s) of interest
    """
	if protein_id == 'NG-GFP-MONOMER[c]':
		return 'GFP'

	else:
		protein = protein_id[:-3]  # subtract the compartment
		common_name = get_gene_symbols_for_monomer_ids()[protein]

	return common_name


class Plot(multigenAnalysisPlot.MultigenAnalysisPlot):
	def do_plot(self, seedOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		with open(simDataFile, 'rb') as f:
			sim_data = pickle.load(f)
		with open(validationDataFile, 'rb') as f:
			validation_data = pickle.load(f)

		cell_paths = self.ap.get_cells()
		sim_dir = cell_paths[0]
		simOutDir = os.path.join(sim_dir, 'simOut')

		# Extract protein information (function from protein_half_lives.py)
		def get_protein_data(sim_data, remove_yibQ):
			protein_ids = sim_data.process.translation.monomer_data['id']
			deg_rate_source = sim_data.process.translation.monomer_data['deg_rate_source']
			degradation_rates = sim_data.process.translation.monomer_data['deg_rate'].asNumber(
				1 / units.s)
			half_lives = np.log(2) / degradation_rates / 60  # in minutes
			if remove_yibQ:
				indices = [i for i, x in enumerate(protein_ids) if x == "EG12298-MONOMER[c]"]
				protein_ids_wo_yibQ = np.delete(protein_ids, indices[0])
				half_lives_wo_yibQ = np.delete(half_lives, indices[0])
				deg_rate_source_wo_yibQ = np.delete(deg_rate_source, indices[0])
				return protein_ids_wo_yibQ, half_lives_wo_yibQ, deg_rate_source_wo_yibQ
			else:
				return protein_ids, half_lives, deg_rate_source

		protein_ids, half_lives, deg_rate_source = (
			get_protein_data(sim_data, remove_yibQ=False))

		# Get the common names for each protein:
		common_names = [get_common_name(protein_id) for protein_id in protein_ids]

		# Extract protein indexes for each new gene
		monomer_counts_reader = TableReader(
			os.path.join(simOutDir, "MonomerCounts"))
		monomer_idx_dict = {monomer: i for i, monomer in
							enumerate(monomer_counts_reader.readAttribute(
								'monomerIds'))}

		# extract the numbers of interest
		monomerIds = monomer_counts_reader.readAttribute(
			"monomerIds")  # this is the one that matches the indexing  I used earlier to construct the listeners!


		# extract the data over all generations:
		# Load data
		time = read_stacked_columns(cell_paths, 'Main', 'time')
		(free_monomer_counts,) = read_stacked_bulk_molecules(
			cell_paths, monomerIds)

		# doubling time function from nora:
		def extract_doubling_times(cell_paths):
			# Load data
			time = read_stacked_columns(cell_paths, 'Main', 'time').squeeze()
			# Determine doubling time
			doubling_times = read_stacked_columns(cell_paths, 'Main', 'time',
												  fun=lambda x: (x[-1] - x[0])).squeeze().astype(
				int)
			end_generation_times = np.cumsum(doubling_times) + time[0]  #
			start_generation_indices = np.searchsorted(time, end_generation_times[:-1],
													   side='left').astype(int)
			start_generation_indices = np.insert(start_generation_indices, 0, 0) + np.arange(
				len(doubling_times))
			end_generation_indices = start_generation_indices + doubling_times
			return time, doubling_times, end_generation_times, start_generation_indices, end_generation_indices

		# extract the doubling times:
		time, doubling_times, end_generation_times, start_generation_indices, end_generation_indices = extract_doubling_times(
			cell_paths)

		# find how many proteins were removed via dilution for each doubling time:
		diluted_counts = np.zeros(((len(doubling_times) - 1), len(monomerIds)))
		diluted_counts_over_time = np.zeros(((len(time)), len(monomerIds)))
		for i in range(
				len(doubling_times) - 1):  # -1 is to account for not doing the last generation
			end_gen = end_generation_indices[i]
			start_gen = start_generation_indices[i + 1]  # the first start is zero, so skip that
			print(end_gen, start_gen)

			# find the protein counts at the end of the generation:
			monomer_counts_at_gen_end = free_monomer_counts[end_gen,:]  # get this for each protein

			# find the protein counts at the start of the next generation:
			monomer_counts_at_gen_start = free_monomer_counts[start_gen, :]

			# find the difference between the two:
			protein_counts_removed = monomer_counts_at_gen_end - monomer_counts_at_gen_start
			diluted_counts[i, :] = protein_counts_removed
			diluted_counts_over_time[end_gen,:] = protein_counts_removed  # put it at the start of the next gen in terms of time


		# compute how many proteins were removed via degradation over the entire sim length:
		degraded_counts = read_stacked_columns(cell_paths, 'MonomerCounts', "protein_deg_ES1_counts") # take from evolve state counter!

		# find the average degradation for each protein over all generations:
		generation_averages = np.zeros((len(cell_paths), len(monomerIds)))
		for i in range(len(end_generation_indices)):
			# retrieve the generation start and end time:
			end_gen = end_generation_indices[i]
			start_gen = start_generation_indices[i]

			# find the protein counts at all time points in the generation:
			monomer_counts_at_gen = degraded_counts[start_gen:end_gen, :]
			# todo: add in the diluted counts here too!

			#
			dilution_counts_at_gen = diluted_counts[i - 1, :]
			# add dilution counts to the last time point of the generation:

			monomer_counts_at_last_time_step = monomer_counts_at_gen[-1,
											   :] + dilution_counts_at_gen
			monomer_counts_at_gen[-1, :] = monomer_counts_at_last_time_step

			# take the average:
			generation_averages[i, :] = np.mean(monomer_counts_at_gen, axis=0)

		# find the average degradation for each protein over all generations:
		avg_degraded_counts = np.mean(generation_averages, axis=0)

		# compute the average loss rate for each protein:
		# todo: affirm this is the correct way to account for dilution
		avg_loss_rate = avg_degraded_counts

		# compute how many counts were added via elongation over the entire sim length:
		elongated_counts = read_stacked_columns(cell_paths, 'MonomerCounts',
												"peptide_elongate_ES1_counts")
		# find the average elongation for each protein over all generations:
		generation_averages = np.zeros((len(cell_paths), len(monomerIds)))
		for i in range(len(end_generation_indices)):
			# retrieve the generation start and end time:
			end_gen = end_generation_indices[i]
			start_gen = start_generation_indices[i]

			# find the protein counts at all time points in the generation:
			monomer_counts_at_gen = elongated_counts[start_gen:end_gen, :]

			# take the average:
			generation_averages[i, :] = np.mean(monomer_counts_at_gen, axis=0)

		# find the average elongation for each protein over all generations:
		avg_elongated_counts = np.mean(generation_averages, axis=0)

		# find the average monomer counts over all generations:
		FMC_generation_averages = np.zeros((len(cell_paths), len(monomerIds)))
		for i in range(len(end_generation_indices)):
			# retrieve the generation start and end time:
			end_gen = end_generation_indices[i]
			start_gen = start_generation_indices[i]

			# find the protein counts at all time points in the generation:
			monomer_counts_at_gen = free_monomer_counts[start_gen:end_gen, :]

			# take the average:
			FMC_generation_averages[i, :] = np.mean(monomer_counts_at_gen, axis=0)

		# find the average monomer counts over all generations:
		avg_FMC = np.mean(FMC_generation_averages, axis=0)

		# determine the minimum value:
		# find the smallest nonzero value in the array:
		min_avg_loss_rate_indices = np.nonzero(avg_loss_rate)
		min_avg_loss_rate = avg_loss_rate[min_avg_loss_rate_indices]
		min_avg_elongated_counts_indices = np.nonzero(avg_elongated_counts)
		min_avg_elongated_counts = avg_elongated_counts[min_avg_elongated_counts_indices]
		min_values = [np.min(min_avg_loss_rate), np.min(min_avg_elongated_counts)]
		min_value = np.min(min_values)
		log_factor = min_value * .1  # add this to avoid negative and zero log values

		# compute the log of the loss and production rates
		log_avg_loss_rate = np.log10(avg_loss_rate + log_factor)
		log_avg_production_rate = np.log10(avg_elongated_counts + log_factor)


		# plot the loss rate vs the production rate:
		# Create figure
		fig = go.Figure()
		common_names = np.array(common_names)

		# Scatter plot for all proteins (grey)
		fig.add_trace(go.Scatter(
			x=log_avg_production_rate,
			y=log_avg_loss_rate,
			mode='markers',
            customdata=np.stack((monomerIds, half_lives, deg_rate_source, avg_FMC, common_names), axis=-1),
            hovertemplate=
            "Monomer ID: %{customdata[0]}<br>" +
            "half life: %{customdata[1]}<br>" +
            "source: %{customdata[2]}<br>" +
			"avgerage free monomer counts: %{customdata[3]}<br>" +
			"common name: %{customdata[4]}<br>" +
            "<extra></extra>",
			marker=dict(size=5, color='lightseagreen', opacity=0.3),
			name="All Proteins"))

		red_name = ''
		if len(HIGHLIGHT_IN_RED) > 0:
			monomerIds = np.array(monomerIds)
			common_names = np.array(common_names)
			# Scatter plot for red proteins
			red_protein_indices = [monomer_idx_dict[protein] for protein in HIGHLIGHT_IN_RED]
			fig.add_trace(go.Scatter(
				x=log_avg_production_rate[red_protein_indices],
				y=log_avg_loss_rate[red_protein_indices],
				mode='markers',
				customdata=np.stack(
					(monomerIds[red_protein_indices], half_lives[red_protein_indices], deg_rate_source[red_protein_indices], avg_FMC[red_protein_indices], common_names[red_protein_indices]), axis=-1),
				hovertemplate=
				"Monomer ID: %{customdata[0]}<br>" +
				"half life: %{customdata[1]}<br>" +
				"source: %{customdata[2]}<br>" +
				"avgerage free monomer counts: %{customdata[3]}<br>" +
				"common name: %{customdata[4]}<br>" +
				"<extra></extra>",
				marker=dict(size=5, color='red', opacity=1),
				name='Red Proteins'))
			# indicate the name should be in the title:
			red_name = str(common_names[red_protein_indices])

		blue_name = ""
		if len(HIGHLIGHT_IN_BLUE) > 0:
			monomerIds = np.array(monomerIds)
			# Scatter plot for blue proteins
			blue_protein_indices = [monomer_idx_dict[protein] for protein in HIGHLIGHT_IN_BLUE]
			fig.add_trace(go.Scatter(
				x=log_avg_production_rate[blue_protein_indices],
				y=log_avg_loss_rate[blue_protein_indices],
				mode='markers',
				customdata=np.stack(
					(monomerIds[blue_protein_indices], half_lives[blue_protein_indices],
					 deg_rate_source[blue_protein_indices], avg_FMC[blue_protein_indices],
					 common_names[blue_protein_indices]), axis=-1),
				hovertemplate=
				"Monomer ID: %{customdata[0]}<br>" +
				"half life: %{customdata[1]}<br>" +
				"source: %{customdata[2]}<br>" +
				"avgerage free monomer counts: %{customdata[3]}<br>" +
				"common name: %{customdata[4]}<br>" +
				"<extra></extra>",
				marker=dict(size=5, color='blue', opacity=1),
				name='Blue Proteins'))
			# indicate the name should be in the title:
			blue_name = str(common_names[blue_protein_indices])

		purple_name = ""
		if len(HIGHLIGHT_IN_PURPLE) > 0:
			monomerIds = np.array(monomerIds)
			# Scatter plot for purple proteins
			purple_protein_indices = [monomer_idx_dict[protein] for protein in HIGHLIGHT_IN_PURPLE]
			fig.add_trace(go.Scatter(
				x=log_avg_production_rate[purple_protein_indices],
				y=log_avg_loss_rate[purple_protein_indices],
				mode='markers',
				customdata=np.stack(
					(monomerIds[purple_protein_indices], half_lives[purple_protein_indices],
					 deg_rate_source[purple_protein_indices], avg_FMC[purple_protein_indices],
					 common_names[purple_protein_indices]), axis=-1),
				hovertemplate=
				"Monomer ID: %{customdata[0]}<br>" +
				"half life: %{customdata[1]}<br>" +
				"source: %{customdata[2]}<br>" +
				"avgerage free monomer counts: %{customdata[3]}<br>" +
				"common name: %{customdata[4]}<br>" +
				"<extra></extra>",
				marker=dict(size=5, color='hotpink', opacity=1),
				name='Purple Proteins'))

			# indicate the name should be in the title:
			purple_name = str(common_names[purple_protein_indices])


		# Generate line data
		x = np.linspace(-4, 3, 100)

		# y = x line (black)
		fig.add_trace(go.Scatter(
			x=x, y=x, mode='lines', line=dict(color='black', width=2,),
			name='y = x'
		))

		# y = x + 1 line (green, dashed)
		fig.add_trace(go.Scatter(
			x=x, y=x + 1, mode='lines',
			line=dict(color='green', width=2, dash='dash'), name='y = x + 1'
		))

		# y = x - 1 line (green, dashed)
		fig.add_trace(go.Scatter(
			x=x, y=x - 1, mode='lines',
			line=dict(color='green', width=2, dash='dash'), name='y = x - 1'
		))

		# Layout settings
		fig.update_layout(
			title="Log10 Average Loss Rate vs Log10 Average Production Rate",
			xaxis_title="Log10 Average Production Rate",
			yaxis_title="Log10 Average Loss Rate",
			width=700, height=700,
			showlegend=True,)

		#save the plot:
		plot_name = plotOutFileName + "_red_" + red_name + "_blue_" + blue_name + "_purple_" + purple_name + ".html"
		fig.write_html(os.path.join(plotOutDir, plot_name))
		#exportFigure(plt, plotOutDir, plotOutFileName, metadata)


if __name__ == '__main__':
	Plot().cli()
