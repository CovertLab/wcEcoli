
import pickle
import os
import pandas as pd
from matplotlib import pyplot as plt
# noinspection PyUnresolvedReferences
from scipy.ndimage import uniform_filter1d
from wholecell.utils import units
import numpy as np

from models.ecoli.analysis import cohortAnalysisPlot
from wholecell.analysis.analysis_tools import (exportFigure,
	read_bulk_molecule_counts, read_stacked_bulk_molecules, read_stacked_columns)
from wholecell.io.tablereader import TableReader

# TODO: try sampling instead!!
# TODO: consider calculating the average as the average per generation, and then STD as the +/- across generations for each seed....


SKIP_INITIAL_GENERATIONS = 2

PROTEINS_OF_INTEREST = ['UHPA-MONOMER']
	# lon and hslv and hslu: ["EG10542-MONOMER[c]", "EG11676-MONOMER", "EG11881-MONOMER"]
	# lon substrates: ['G6890-MONOMER[c]', 'PD03938[c]', 'G6737-MONOMER[c]', 'RPOD-MONOMER[c]', 'PD02936[c]', 'RED-THIOREDOXIN2-MONOMER[c]']

class Plot(cohortAnalysisPlot.CohortAnalysisPlot):
	def check_validity_and_get_compartment(self, sim_data, molecule_list):
		"""
        Validate molecule IDs and add compartment tags.

        Args:
            sim_data: Simulation data object
            molecule_list: List of molecule IDs (with or without compartment tags)

        Returns:
            list: Valid molecule IDs with compartment tags
        """
		revised_molecule_list = []
		for molecule in molecule_list:
			if "[" in molecule:
				molecule = molecule[:-3]  # Remove compartment
			if sim_data.getter.is_valid_molecule(molecule):
				revised_name = molecule + sim_data.getter.get_compartment_tag(molecule)
				revised_molecule_list.append(revised_name)
			else:
				print(f"{molecule} is not a valid molecule in the simulation.")

		return revised_molecule_list

	def _generate_summary_statistics_old(self, cell_paths, protein_ids):
		"""Generate summary statistics for all proteins across all seeds/timepoints"""

		# Read all data at once
		free_counts = read_stacked_bulk_molecules(cell_paths, (protein_ids,))[
			0]  # (n_timepoints, n_proteins)
		total_counts = read_stacked_columns(cell_paths, 'MonomerCounts',
											'monomerCounts')  # (n_timepoints, n_proteins)
		counts_to_molar = read_stacked_columns(cell_paths, 'EnzymeKinetics', 'countsToMolar',
											   ignore_exception=True)

		# Handle shape
		if counts_to_molar.ndim > 1:
			counts_to_molar = counts_to_molar.squeeze()
		counts_to_molar = counts_to_molar / 1000  # Convert to M

		# Compute concentrations (broadcasting counts_to_molar across proteins)
		free_conc = free_counts * counts_to_molar[:, np.newaxis]  # (n_timepoints, n_proteins)
		total_conc = total_counts * counts_to_molar[:, np.newaxis]

		# Compute statistics
		summary_df = pd.DataFrame({
			'Protein ID': protein_ids,
			'Avg Free Counts': free_counts.mean(axis=0),
			'Std Free Counts': free_counts.std(axis=0),
			'Avg Total Counts': total_counts.mean(axis=0),
			'Std Total Counts': total_counts.std(axis=0),
			'Avg Free Conc (M)': free_conc.mean(axis=0),
			'Std Free Conc (M)': free_conc.std(axis=0),
			'Avg Total Conc (M)': total_conc.mean(axis=0),
			'Std Total Conc (M)': total_conc.std(axis=0),
		})

		return summary_df

	def _save_seed_metadata(self, seeds, generations, plotOutDir):
		"""Save time, cell volume, dry mass, and counts_to_molar for each seed"""

		seed_df_dir = os.path.join(plotOutDir, 'seed_summaries')
		os.makedirs(seed_df_dir, exist_ok=True)

		for seed in seeds:
			cell_paths = self.ap.get_cells(
				seed=[seed],
				generation=range(SKIP_INITIAL_GENERATIONS, len(generations) - 1)
			)

			time = read_stacked_columns(cell_paths, 'Main', 'time').flatten()
			cellVolume = read_stacked_columns(cell_paths, 'Mass', 'cellVolume').flatten()
			dryMass = read_stacked_columns(cell_paths, 'Mass', 'dryMass').flatten()
			counts_to_molar = read_stacked_columns(cell_paths, 'EnzymeKinetics',
												   'countsToMolar',
												   ignore_exception=True)
			if counts_to_molar.ndim > 1:
				counts_to_molar = counts_to_molar.squeeze()
			counts_to_molar = (counts_to_molar / 1000).flatten()  # Convert to M

			seed_df = pd.DataFrame({
				'time': time,
				'cellVolume': cellVolume,
				'dryMass': dryMass,
				'counts_to_molar': counts_to_molar,
			})

			seed_df.to_csv(os.path.join(seed_df_dir, f'seed_{seed}_summary.csv'), index=False)

	def _analyze_proteins_of_interest(self, proteins_to_analyze, all_protein_ids, seeds,
									  total_generations, plotOutDir):
		"""Analyze specific proteins of interest in detail"""

		# Store for use in other methods
		self.total_generations = total_generations

		# Create mapping from protein ID to index
		protein_id_to_idx = {prot_id: idx for idx, prot_id in enumerate(all_protein_ids)}

		for prot_id in proteins_to_analyze:
			if prot_id not in protein_id_to_idx:
				print(f"Warning: {prot_id} not found in protein IDs, skipping...")
				continue

			prot_idx = protein_id_to_idx[prot_id]
			print(f"  Processing {prot_id}...")

			# Collect data for this protein across all seeds
			seed_dataframes = {}
			seed_gen_stats = {}  # NEW: Store generation-averaged stats

			for seed in seeds:
				cell_paths = self.ap.get_cells(
					seed=[seed],
					generation=range(SKIP_INITIAL_GENERATIONS, total_generations - 1)
				)

				# Initialize lists to collect data across all generations
				all_time = []
				all_cellVolume = []
				all_dryMass = []
				all_counts_to_molar = []
				all_free_counts = []
				all_total_counts = []
				all_free_conc = []
				all_total_conc = []

				# NEW: Lists to store per-generation averages
				gen_free_counts = []
				gen_total_counts = []
				gen_free_conc = []
				gen_total_conc = []

				# Read data for each generation separately
				for cell_path in cell_paths:
					time = read_stacked_columns([cell_path], 'Main', 'time').flatten()
					cellVolume = read_stacked_columns([cell_path], 'Mass', 'cellVolume').flatten()
					dryMass = read_stacked_columns([cell_path], 'Mass', 'dryMass').flatten()

					# Get counts for this specific protein
					free_counts = read_stacked_bulk_molecules([cell_path], [prot_id])[0]
					total_counts = read_stacked_columns([cell_path], 'MonomerCounts',
														'monomerCounts')[:, prot_idx]

					# Get counts_to_molar
					counts_to_molar = read_stacked_columns([cell_path], 'EnzymeKinetics',
														   'countsToMolar',
														   ignore_exception=True)
					if counts_to_molar.ndim > 1:
						counts_to_molar = counts_to_molar.squeeze()
					counts_to_molar = (counts_to_molar / 1000).flatten()  # Convert to M

					# Compute concentrations
					free_conc = free_counts * counts_to_molar
					total_conc = total_counts * counts_to_molar

					# Append to overall lists
					all_time.extend(time)
					all_cellVolume.extend(cellVolume)
					all_dryMass.extend(dryMass)
					all_counts_to_molar.extend(counts_to_molar)
					all_free_counts.extend(free_counts)
					all_total_counts.extend(total_counts)
					all_free_conc.extend(free_conc)
					all_total_conc.extend(total_conc)

					# NEW: Store generation averages
					gen_free_counts.append(np.mean(free_counts))
					gen_total_counts.append(np.mean(total_counts))
					gen_free_conc.append(np.mean(free_conc))
					gen_total_conc.append(np.mean(total_conc))

				# Create dataframe with all data
				df = pd.DataFrame({
					'time': all_time,
					'cellVolume': all_cellVolume,
					'dryMass': all_dryMass,
					'counts_to_molar': all_counts_to_molar,
					'free_counts': all_free_counts,
					'total_counts': all_total_counts,
					'free_conc': all_free_conc,
					'total_conc': all_total_conc,
				})

				seed_dataframes[seed] = df

				# NEW: Compute and store generation-averaged statistics
				seed_gen_stats[seed] = {
					'free_counts_gen_avg': np.mean(gen_free_counts),
					'free_counts_gen_std': np.std(gen_free_counts),
					'total_counts_gen_avg': np.mean(gen_total_counts),
					'total_counts_gen_std': np.std(gen_total_counts),
					'free_conc_gen_avg': np.mean(gen_free_conc),
					'free_conc_gen_std': np.std(gen_free_conc),
					'total_conc_gen_avg': np.mean(gen_total_conc),
					'total_conc_gen_std': np.std(gen_total_conc),
				}

				# Save dataframe
				safe_prot_id = prot_id.replace('[', '_').replace(']', '_').replace('/', '_')
				df_dir = os.path.join(plotOutDir, 'protein_dataframes', safe_prot_id)
				os.makedirs(df_dir, exist_ok=True)
				df.to_csv(os.path.join(df_dir, f'seed_{seed}.csv'), index=False)

			# Create plots for this protein (now passing gen_stats)
			self._plot_protein_analysis(prot_id, prot_idx, seed_dataframes, seed_gen_stats,
										seeds, plotOutDir)

	def _plot_protein_analysis(self, prot_id, prot_idx, seed_dataframes, seed_gen_stats,
							   seeds, plotOutDir):
		"""Create comprehensive plots for a single protein"""

		# Get half-life and common name for this protein
		half_life = self.half_lives[prot_idx]
		common_name = self.common_names[prot_idx]

		fig = plt.figure(figsize=(20, 12))

		# Create 2x2 subplot layout
		ax1 = plt.subplot(2, 2, 1)  # Free counts
		ax2 = plt.subplot(2, 2, 2)  # Free concentration
		ax3 = plt.subplot(2, 2, 3)  # Total counts
		ax4 = plt.subplot(2, 2, 4)  # Total concentration

		colors = plt.cm.tab10(np.linspace(0, 1, len(seeds)))

		# Plot individual seeds
		for idx, seed in enumerate(seeds):
			df = seed_dataframes[seed]
			stats = seed_gen_stats[seed]
			time = df['time'].values / 60  # Convert to minutes

			# Get generation-averaged statistics from pre-computed values
			free_counts_gen_avg = stats['free_counts_gen_avg']
			free_counts_gen_std = stats['free_counts_gen_std']
			total_counts_gen_avg = stats['total_counts_gen_avg']
			total_counts_gen_std = stats['total_counts_gen_std']
			free_conc_gen_avg = stats['free_conc_gen_avg']
			free_conc_gen_std = stats['free_conc_gen_std']
			total_conc_gen_avg = stats['total_conc_gen_avg']
			total_conc_gen_std = stats['total_conc_gen_std']

			# Calculate simple time-averaged statistics
			free_counts_time_avg = df['free_counts'].mean()
			free_counts_time_std = df['free_counts'].std()
			total_counts_time_avg = df['total_counts'].mean()
			total_counts_time_std = df['total_counts'].std()
			free_conc_time_avg = df['free_conc'].mean()
			free_conc_time_std = df['free_conc'].std()
			total_conc_time_avg = df['total_conc'].mean()
			total_conc_time_std = df['total_conc'].std()

			# Plot free counts with both stats in legend
			ax1.plot(time, df['free_counts'], alpha=0.6, color=colors[idx],
					 label=f'Seed {seed}:\n  Time: {free_counts_time_avg:.0f}±{free_counts_time_std:.0f}\n  Gen: {free_counts_gen_avg:.0f}±{free_counts_gen_std:.0f}')

			# Plot free concentration (µM)
			ax2.plot(time, df['free_conc'] * 1e6, alpha=0.6, color=colors[idx],
					 label=f'Seed {seed}:\n  Time: {free_conc_time_avg * 1e6:.2f}±{free_conc_time_std * 1e6:.2f} µM\n  Gen: {free_conc_gen_avg * 1e6:.2f}±{free_conc_gen_std * 1e6:.2f} µM')

			# Plot total counts
			ax3.plot(time, df['total_counts'], alpha=0.6, color=colors[idx],
					 label=f'Seed {seed}:\n  Time: {total_counts_time_avg:.0f}±{total_counts_time_std:.0f}\n  Gen: {total_counts_gen_avg:.0f}±{total_counts_gen_std:.0f}')

			# Plot total concentration (µM)
			ax4.plot(time, df['total_conc'] * 1e6, alpha=0.6, color=colors[idx],
					 label=f'Seed {seed}:\n  Time: {total_conc_time_avg * 1e6:.2f}±{total_conc_time_std * 1e6:.2f} µM\n  Gen: {total_conc_gen_avg * 1e6:.2f}±{total_conc_gen_std * 1e6:.2f} µM')

		# Formatting
		for ax, title, ylabel in zip(
				[ax1, ax2, ax3, ax4],
				['Free Counts', 'Free Concentration', 'Total Counts', 'Total Concentration'],
				['Counts', 'Concentration (µM)', 'Counts', 'Concentration (µM)']
		):
			ax.set_xlabel('Time (minutes)', fontsize=11)
			ax.set_ylabel(ylabel, fontsize=11)
			ax.set_title(title, fontsize=12, fontweight='bold')
			ax.legend(loc='best', fontsize=7)  # Reduced font size slightly for multi-line legends
			ax.grid(True, alpha=0.3)

		plt.suptitle(
			f'Protein Analysis: {prot_id}\n{common_name} | Half-life: {half_life:.1f} min | generations {SKIP_INITIAL_GENERATIONS -1}-{self.total_generations -1} plotted | {len(seeds) * (self.total_generations - SKIP_INITIAL_GENERATIONS)} cells total',
			fontsize=14, fontweight='bold')
		plt.tight_layout()

		# Save figure
		safe_prot_id = prot_id.replace('[', '_').replace(']', '_').replace('/', '_')
		fig_path = os.path.join(plotOutDir, 'protein_plots')
		os.makedirs(fig_path, exist_ok=True)
		plt.savefig(os.path.join(fig_path, f'{safe_prot_id}_analysis.png'), dpi=300,
					bbox_inches='tight')
		plt.close()

	def _compute_all_proteins_generation_stats(self, protein_ids, seeds):
		"""Compute generation-averaged statistics for all proteins across all seeds"""

		n_proteins = len(protein_ids)

		# Initialize arrays to store generation averages for each protein
		free_counts_gen_means = np.zeros((len(seeds), n_proteins))
		total_counts_gen_means = np.zeros((len(seeds), n_proteins))
		free_conc_gen_means = np.zeros((len(seeds), n_proteins))
		total_conc_gen_means = np.zeros((len(seeds), n_proteins))

		# Loop through seeds
		for seed_idx, seed in enumerate(seeds):
			print(f"  Processing seed {seed} for generation-averaged stats...")

			# Get cells for this seed
			cell_paths = self.ap.get_cells(
				seed=[seed],
				generation=range(SKIP_INITIAL_GENERATIONS, self.total_generations - 1)
			)

			# Initialize lists to store per-generation means
			gen_free_counts = []
			gen_total_counts = []
			gen_free_conc = []
			gen_total_conc = []

			# Process each generation
			for cell_path in cell_paths:
				# Read data for this generation
				free_counts_gen = read_stacked_bulk_molecules([cell_path], (protein_ids,))[0]
				total_counts_gen = read_stacked_columns([cell_path], 'MonomerCounts',
														'monomerCounts')
				counts_to_molar_gen = read_stacked_columns([cell_path], 'EnzymeKinetics',
														   'countsToMolar',
														   ignore_exception=True)

				if counts_to_molar_gen.ndim > 1:
					counts_to_molar_gen = counts_to_molar_gen.squeeze()
				counts_to_molar_gen = counts_to_molar_gen / 1000  # Convert to M

				# Compute concentrations
				free_conc_gen = free_counts_gen * counts_to_molar_gen[:, np.newaxis]
				total_conc_gen = total_counts_gen * counts_to_molar_gen[:, np.newaxis]

				# Store mean for this generation (average over time within generation)
				gen_free_counts.append(free_counts_gen.mean(axis=0))
				gen_total_counts.append(total_counts_gen.mean(axis=0))
				gen_free_conc.append(free_conc_gen.mean(axis=0))
				gen_total_conc.append(total_conc_gen.mean(axis=0))

			# Convert to arrays and compute mean across generations for this seed
			free_counts_gen_means[seed_idx, :] = np.mean(gen_free_counts, axis=0)
			total_counts_gen_means[seed_idx, :] = np.mean(gen_total_counts, axis=0)
			free_conc_gen_means[seed_idx, :] = np.mean(gen_free_conc, axis=0)
			total_conc_gen_means[seed_idx, :] = np.mean(gen_total_conc, axis=0)

		# Now compute statistics across seeds
		return {
			'Avg Free Counts (Gen)': free_counts_gen_means.mean(axis=0),
			'Std Free Counts (Gen)': free_counts_gen_means.std(axis=0),
			'Avg Total Counts (Gen)': total_counts_gen_means.mean(axis=0),
			'Std Total Counts (Gen)': total_counts_gen_means.std(axis=0),
			'Avg Free Conc (M) (Gen)': free_conc_gen_means.mean(axis=0),
			'Std Free Conc (M) (Gen)': free_conc_gen_means.std(axis=0),
			'Avg Total Conc (M) (Gen)': total_conc_gen_means.mean(axis=0),
			'Std Total Conc (M) (Gen)': total_conc_gen_means.std(axis=0),
		}


	def _plot_averaged_data(self, seed_dataframes, seeds, ax1, ax2, ax3, ax4):
		"""Calculate and plot averaged data across seeds with smoothed standard deviation"""

		# Find common time range across all seeds
		min_time = max([seed_dataframes[seed]['time'].min() for seed in seeds])
		max_time = min([seed_dataframes[seed]['time'].max() for seed in seeds])

		# Create common time grid
		n_points = 1000
		common_time = np.linspace(min_time, max_time, n_points)

		# Interpolate all seeds to common time grid
		free_counts_interp = []
		free_conc_interp = []
		total_counts_interp = []
		total_conc_interp = []

		for seed in seeds:
			df = seed_dataframes[seed]
			time = df['time'].values

			free_counts_interp.append(np.interp(common_time, time, df['free_counts'].values))
			free_conc_interp.append(np.interp(common_time, time, df['free_conc'].values))
			total_counts_interp.append(np.interp(common_time, time, df['total_counts'].values))
			total_conc_interp.append(np.interp(common_time, time, df['total_conc'].values))

		# Calculate mean and std
		free_counts_mean = np.mean(free_counts_interp, axis=0)
		free_counts_std = np.std(free_counts_interp, axis=0)
		free_conc_mean = np.mean(free_conc_interp, axis=0)
		free_conc_std = np.std(free_conc_interp, axis=0)
		total_counts_mean = np.mean(total_counts_interp, axis=0)
		total_counts_std = np.std(total_counts_interp, axis=0)
		total_conc_mean = np.mean(total_conc_interp, axis=0)
		total_conc_std = np.std(total_conc_interp, axis=0)

		# Smooth the standard deviations
		window_size = min(50, len(common_time) // 20)  # Adaptive window size
		if window_size > 1:
			free_counts_std_smooth = uniform_filter1d(free_counts_std, size=window_size)
			free_conc_std_smooth = uniform_filter1d(free_conc_std, size=window_size)
			total_counts_std_smooth = uniform_filter1d(total_counts_std, size=window_size)
			total_conc_std_smooth = uniform_filter1d(total_conc_std, size=window_size)
		else:
			free_counts_std_smooth = free_counts_std
			free_conc_std_smooth = free_conc_std
			total_counts_std_smooth = total_counts_std
			total_conc_std_smooth = total_conc_std

		# Convert time to mins
		common_time_hours = common_time / 60

		# Plot averages with shaded standard deviation
		# Free counts
		ax1.plot(common_time_hours, free_counts_mean, 'k-', linewidth=2.5,
				 label=f'Average: {free_counts_mean.mean():.0f}±{free_counts_std.mean():.0f}',
				 zorder=10)
		ax1.fill_between(common_time_hours,
						 free_counts_mean - free_counts_std_smooth,
						 free_counts_mean + free_counts_std_smooth,
						 color='black', alpha=0.2, zorder=5)

		# Free concentration (convert to µM)
		ax2.plot(common_time_hours, free_conc_mean * 1e6, 'k-', linewidth=2.5,
				 label=f'Average: {free_conc_mean.mean() * 1e6:.2f}±{free_conc_std.mean() * 1e6:.2f} µM',
				 zorder=10)
		ax2.fill_between(common_time_hours,
						 (free_conc_mean - free_conc_std_smooth) * 1e6,
						 (free_conc_mean + free_conc_std_smooth) * 1e6,
						 color='black', alpha=0.2, zorder=5)

		# Total counts
		ax3.plot(common_time_hours, total_counts_mean, 'k-', linewidth=2.5,
				 label=f'Average: {total_counts_mean.mean():.0f}±{total_counts_std.mean():.0f}',
				 zorder=10)
		ax3.fill_between(common_time_hours,
						 total_counts_mean - total_counts_std_smooth,
						 total_counts_mean + total_counts_std_smooth,
						 color='black', alpha=0.2, zorder=5)

		# Total concentration (convert to µM)
		ax4.plot(common_time_hours, total_conc_mean * 1e6, 'k-', linewidth=2.5,
				 label=f'Average: {total_conc_mean.mean() * 1e6:.2f}±{total_conc_std.mean() * 1e6:.2f} µM',
				 zorder=10)
		ax4.fill_between(common_time_hours,
						 (total_conc_mean - total_conc_std_smooth) * 1e6,
						 (total_conc_mean + total_conc_std_smooth) * 1e6,
						 color='black', alpha=0.2, zorder=5)

	hi = 5

	def _generate_summary_statistics(self, cell_paths, protein_ids, seeds):
		"""Generate summary statistics for all proteins across all seeds/timepoints"""

		# Read all data at once
		free_counts = read_stacked_bulk_molecules(cell_paths, (protein_ids,))[0]
		total_counts = read_stacked_columns(cell_paths, 'MonomerCounts', 'monomerCounts')
		counts_to_molar = read_stacked_columns(cell_paths, 'EnzymeKinetics', 'countsToMolar',
											   ignore_exception=True)

		# Handle shape
		if counts_to_molar.ndim > 1:
			counts_to_molar = counts_to_molar.squeeze()
		counts_to_molar = counts_to_molar / 1000  # Convert to M

		# Compute concentrations
		free_conc = free_counts * counts_to_molar[:, np.newaxis]
		total_conc = total_counts * counts_to_molar[:, np.newaxis]

		# Compute simple statistics (time-averaged)
		summary_df = pd.DataFrame({
			'Protein ID': protein_ids,
			'Avg Free Counts (Time)': free_counts.mean(axis=0),
			'Std Free Counts (Time)': free_counts.std(axis=0),
			'Avg Total Counts (Time)': total_counts.mean(axis=0),
			'Std Total Counts (Time)': total_counts.std(axis=0),
			'Avg Free Conc (M) (Time)': free_conc.mean(axis=0),
			'Std Free Conc (M) (Time)': free_conc.std(axis=0),
			'Avg Total Conc (M) (Time)': total_conc.mean(axis=0),
			'Std Total Conc (M) (Time)': total_conc.std(axis=0),
		})

		# Compute generation-averaged statistics for all proteins
		print("Computing generation-averaged statistics for all proteins...")
		gen_stats = self._compute_all_proteins_generation_stats(protein_ids, seeds)

		# Add to summary dataframe
		for col_name, values in gen_stats.items():
			summary_df[col_name] = values

		return summary_df

	hi = 5



	def do_plot(self, variantDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		with open(simDataFile, 'rb') as f:
			sim_data = pickle.load(f)

		# Sim information:
		sim_id = metadata['description']
		generations = metadata.get('generations', None)
		seeds = self.ap.get_seeds()
		all_cells = self.ap.get_cells(generation=range(SKIP_INITIAL_GENERATIONS, generations))

		# Get all cell paths excluding initial generations if specified:
		cell_paths = self.ap.get_cells(
			generation=range(SKIP_INITIAL_GENERATIONS, len(all_cells) - 1))

		# Extract the protein IDs:
		protein_ids = sim_data.process.translation.monomer_data['id']
		# Extract monomer information:
		degradation_rates = sim_data.process.translation.monomer_data['deg_rate'].asNumber(
			1 / units.s)
		self.half_lives = np.log(2) / degradation_rates / 60
		cistron_ids = sim_data.process.translation.monomer_data['cistron_id']
		self.common_names = [sim_data.common_names.get_common_name(id) for id in cistron_ids]

		# Validate and get proteins to analyze
		if len(PROTEINS_OF_INTEREST) > 0:
			proteins_to_analyze = self.check_validity_and_get_compartment(sim_data,
																		  PROTEINS_OF_INTEREST)
			print(f"Analyzing {len(proteins_to_analyze)} proteins of interest")
		else:
			print(
				"WARNING: No proteins specified. Generating summary for all proteins but no individual plots.")
			proteins_to_analyze = []

		print("Generating summary statistics for all proteins...")


		# ===== PROTEIN-SPECIFIC ANALYSIS (only for proteins of interest) =====
		if len(proteins_to_analyze) > 0:
			print(f"Analyzing {len(proteins_to_analyze)} proteins of interest...")
			self._analyze_proteins_of_interest(
				proteins_to_analyze,
				protein_ids,
				seeds,
				generations,
				plotOutDir
			)

		# Create a simple summary plot
		plt.figure(figsize=(12, 8))
		summary_text = f'Analysis Complete\n'
		summary_text += f'Total proteins: {len(protein_ids)}\n'
		summary_text += f'Proteins analyzed in detail: {len(proteins_to_analyze)}\n'
		summary_text += f'Seeds: {len(seeds)}\n'
		summary_text += f'See:\n  - protein_summary.csv for all protein statistics\n'
		summary_text += f'  - seed_summaries/ for per-seed metadata\n'
		if len(proteins_to_analyze) > 0:
			summary_text += f'  - protein_dataframes/ for detailed time-series data\n'
			summary_text += f'  - protein_plots/ for visualizations'

		plt.text(0.5, 0.5, summary_text, ha='center', va='center', fontsize=12, family='monospace')
		plt.axis('off')
		plt.tight_layout()
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close('all')

		print("Analysis complete!")



		hi = 5
		# ===== SUMMARY STATISTICS ACROSS ALL SEEDS AND TIMEPOINTS =====
		# Save the summary dataframe last since it takes a while:
		# summary_df = self._generate_summary_statistics(cell_paths, protein_ids)
		# summary_df.to_csv(os.path.join(plotOutDir, 'protein_summary.csv'), index=False)
		# print(f"Saved summary for {len(protein_ids)} proteins")

		self.total_generations = generations  # Store for use in other methods
		summary_df = self._generate_summary_statistics(cell_paths, protein_ids, seeds)
		summary_df.to_csv(os.path.join(plotOutDir, 'protein_summary_new.csv'), index=False)
		print(f"Saved summary for {len(protein_ids)} proteins")

		# ===== SEED-SPECIFIC METADATA =====
		print("Saving seed-specific metadata...")
		self._save_seed_metadata(seeds, all_cells, plotOutDir)
		hi = 6



		# for protein in proteins_to_plot, loop through each seed and save a df with the time, cell volume, dry mass, counts_to_molar, free counts, total counts, free concentration, and total concentration for that protein for that seed.
		# also, for each protein, loop through the seeds and plot the free counts of all seeds on the same plot with the average over the matching timepoints under it (please have STD ploted smoothed around it), and do the same for total counts, free concs, and total concs.





if __name__ == '__main__':
	Plot().cli()
