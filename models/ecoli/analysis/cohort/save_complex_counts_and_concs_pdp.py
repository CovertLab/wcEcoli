import pickle
import os
import pandas as pd
from matplotlib import pyplot as plt
from scipy.ndimage import uniform_filter1d
from wholecell.utils import units
import numpy as np

from models.ecoli.analysis import cohortAnalysisPlot
from wholecell.analysis.analysis_tools import (exportFigure,
    read_bulk_molecule_counts, read_stacked_bulk_molecules, read_stacked_columns)
from wholecell.io.tablereader import TableReader

SKIP_INITIAL_GENERATIONS = 2

COMPLEXES_OF_INTEREST = ['CPLX0-125', ]
# lon complex: ["CPLX0-2881"]
    # hslV ["CPLX0-1163", "CPLX0-1161", "CPLX0-1162"]

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
        hi = 5
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

    def _analyze_complexes_of_interest(self, complexes_to_analyze, seeds,
                                       total_generations, plotOutDir):
        """Analyze specific complexes of interest in detail"""

        # Store for use in other methods
        self.total_generations = total_generations

        for complex_id in complexes_to_analyze:
            print(f"  Processing {complex_id}...")

            # Collect data for this complex across all seeds
            seed_dataframes = {}
            seed_gen_stats = {}  # Store generation-averaged stats

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
                all_counts = []
                all_conc = []

                # Lists to store per-generation averages
                gen_counts = []
                gen_conc = []

                # Read data for each generation separately
                for cell_path in cell_paths:
                    time = read_stacked_columns([cell_path], 'Main', 'time').flatten()
                    cellVolume = read_stacked_columns([cell_path], 'Mass', 'cellVolume').flatten()
                    dryMass = read_stacked_columns([cell_path], 'Mass', 'dryMass').flatten()

                    # Get counts for this specific complex
                    counts = read_stacked_bulk_molecules([cell_path], [complex_id])[0]

                    # Get counts_to_molar
                    counts_to_molar = read_stacked_columns([cell_path], 'EnzymeKinetics',
                                                           'countsToMolar',
                                                           ignore_exception=True)
                    if counts_to_molar.ndim > 1:
                        counts_to_molar = counts_to_molar.squeeze()
                    counts_to_molar = (counts_to_molar / 1000).flatten()  # Convert to M

                    # Compute concentrations
                    conc = counts * counts_to_molar

                    # Append to overall lists
                    all_time.extend(time)
                    all_cellVolume.extend(cellVolume)
                    all_dryMass.extend(dryMass)
                    all_counts_to_molar.extend(counts_to_molar)
                    all_counts.extend(counts)
                    all_conc.extend(conc)

                    # Store generation averages
                    gen_counts.append(np.mean(counts))
                    gen_conc.append(np.mean(conc))

                # Create dataframe with all data
                df = pd.DataFrame({
                    'time': all_time,
                    'cellVolume': all_cellVolume,
                    'dryMass': all_dryMass,
                    'counts_to_molar': all_counts_to_molar,
                    'counts': all_counts,
                    'conc': all_conc,
                })

                seed_dataframes[seed] = df

                # Compute and store generation-averaged statistics
                seed_gen_stats[seed] = {
                    'counts_gen_avg': np.mean(gen_counts),
                    'counts_gen_std': np.std(gen_counts),
                    'conc_gen_avg': np.mean(gen_conc),
                    'conc_gen_std': np.std(gen_conc),
                }

                # Save dataframe
                safe_complex_id = complex_id.replace('[', '_').replace(']', '_').replace('/', '_')
                df_dir = os.path.join(plotOutDir, 'complex_dataframes', safe_complex_id)
                os.makedirs(df_dir, exist_ok=True)
                df.to_csv(os.path.join(df_dir, f'seed_{seed}.csv'), index=False)

            # Create plots for this complex
            self._plot_complex_analysis(complex_id, seed_dataframes, seed_gen_stats,
                                        seeds, plotOutDir)

    def _plot_complex_analysis(self, complex_id, seed_dataframes, seed_gen_stats,
                               seeds, plotOutDir):
        """Create comprehensive plots for a single complex"""

        fig = plt.figure(figsize=(20, 6))

        # Create 1x2 subplot layout (only counts and concentration)
        ax1 = plt.subplot(1, 2, 1)  # Counts
        ax2 = plt.subplot(1, 2, 2)  # Concentration

        colors = plt.cm.tab10(np.linspace(0, 1, len(seeds)))

        # Plot individual seeds
        for idx, seed in enumerate(seeds):
            df = seed_dataframes[seed]
            stats = seed_gen_stats[seed]
            time = df['time'].values / 60  # Convert to minutes

            # Get generation-averaged statistics from pre-computed values
            counts_gen_avg = stats['counts_gen_avg']
            counts_gen_std = stats['counts_gen_std']
            conc_gen_avg = stats['conc_gen_avg']
            conc_gen_std = stats['conc_gen_std']

            # Calculate simple time-averaged statistics
            counts_time_avg = df['counts'].mean()
            counts_time_std = df['counts'].std()
            conc_time_avg = df['conc'].mean()
            conc_time_std = df['conc'].std()

            # Plot counts with both stats in legend
            ax1.plot(time, df['counts'], alpha=0.6, color=colors[idx],
                     label=f'Seed {seed}:\n  Time: {counts_time_avg:.0f}±{counts_time_std:.0f}\n  Gen: {counts_gen_avg:.0f}±{counts_gen_std:.0f}')

            # Plot concentration (µM)
            ax2.plot(time, df['conc'] * 1e6, alpha=0.6, color=colors[idx],
                     label=f'Seed {seed}:\n  Time: {conc_time_avg * 1e6:.2f}±{conc_time_std * 1e6:.2f} µM\n  Gen: {conc_gen_avg * 1e6:.2f}±{conc_gen_std * 1e6:.2f} µM')

        # Formatting
        for ax, title, ylabel in zip(
                [ax1, ax2],
                ['Counts', 'Concentration'],
                ['Counts', 'Concentration (µM)']
        ):
            ax.set_xlabel('Time (minutes)', fontsize=11)
            ax.set_ylabel(ylabel, fontsize=11)
            ax.set_title(title, fontsize=12, fontweight='bold')
            ax.legend(loc='best', fontsize=7)
            ax.grid(True, alpha=0.3)

        plt.suptitle(
            f'Complex Analysis: {complex_id}\nGenerations {SKIP_INITIAL_GENERATIONS}-{self.total_generations - 1} plotted | {len(seeds) * (self.total_generations - SKIP_INITIAL_GENERATIONS)} cells total',
            fontsize=14, fontweight='bold')
        plt.tight_layout()

        # Save figure
        safe_complex_id = complex_id.replace('[', '_').replace(']', '_').replace('/', '_')
        fig_path = os.path.join(plotOutDir, 'complex_plots')
        os.makedirs(fig_path, exist_ok=True)
        plt.savefig(os.path.join(fig_path, f'{safe_complex_id}_analysis.png'), dpi=300,
                    bbox_inches='tight')
        plt.close()

    def _compute_all_complexes_generation_stats(self, complex_ids, seeds):
        """Compute generation-averaged statistics for all complexes across all seeds"""

        n_complexes = len(complex_ids)

        # Initialize arrays to store generation averages for each complex
        counts_gen_means = np.zeros((len(seeds), n_complexes))
        conc_gen_means = np.zeros((len(seeds), n_complexes))

        # Loop through seeds
        for seed_idx, seed in enumerate(seeds):
            print(f"  Processing seed {seed} for generation-averaged stats...")

            # Get cells for this seed
            cell_paths = self.ap.get_cells(
                seed=[seed],
                generation=range(SKIP_INITIAL_GENERATIONS, self.total_generations - 1)
            )

            # Initialize lists to store per-generation means
            gen_counts = []
            gen_conc = []

            # Process each generation
            for cell_path in cell_paths:
                # Read data for this generation
                # FIXED: Remove the extra tuple wrapping
                counts_gen = read_stacked_bulk_molecules([cell_path], complex_ids)[0]
                counts_to_molar_gen = read_stacked_columns([cell_path], 'EnzymeKinetics',
                                                           'countsToMolar',
                                                           ignore_exception=True)

                if counts_to_molar_gen.ndim > 1:
                    counts_to_molar_gen = counts_to_molar_gen.squeeze()
                counts_to_molar_gen = counts_to_molar_gen / 1000  # Convert to M

                # Compute concentrations
                conc_gen = counts_gen * counts_to_molar_gen[:, np.newaxis]

                # Store mean for this generation (average over time within generation)
                gen_counts.append(counts_gen.mean(axis=0))
                gen_conc.append(conc_gen.mean(axis=0))

            # Convert to arrays and compute mean across generations for this seed
            counts_gen_means[seed_idx, :] = np.mean(gen_counts, axis=0)
            conc_gen_means[seed_idx, :] = np.mean(gen_conc, axis=0)

        # Now compute statistics across seeds
        return {
            'Avg Counts (Gen)': counts_gen_means.mean(axis=0),
            'Std Counts (Gen)': counts_gen_means.std(axis=0),
            'Avg Conc (M) (Gen)': conc_gen_means.mean(axis=0),
            'Std Conc (M) (Gen)': conc_gen_means.std(axis=0),
        }

    def _generate_summary_statistics(self, cell_paths, complex_ids, seeds):
        """Generate summary statistics for all complexes across all seeds/timepoints"""

        # Read all data at once
        # FIXED: Remove the extra tuple wrapping
        counts = read_stacked_bulk_molecules(cell_paths, complex_ids)[0]
        counts_to_molar = read_stacked_columns(cell_paths, 'EnzymeKinetics', 'countsToMolar',
                                               ignore_exception=True)

        # Handle shape
        if counts_to_molar.ndim > 1:
            counts_to_molar = counts_to_molar.squeeze()
        counts_to_molar = counts_to_molar / 1000  # Convert to M

        # Compute concentrations
        conc = counts * counts_to_molar[:, np.newaxis]

        # Compute simple statistics (time-averaged)
        summary_df = pd.DataFrame({
            'Complex ID': complex_ids,
            'Avg Counts (Time)': counts.mean(axis=0),
            'Std Counts (Time)': counts.std(axis=0),
            'Avg Conc (M) (Time)': conc.mean(axis=0),
            'Std Conc (M) (Time)': conc.std(axis=0),
        })

        # Compute generation-averaged statistics for all complexes
        print("Computing generation-averaged statistics for all complexes...")
        gen_stats = self._compute_all_complexes_generation_stats(complex_ids, seeds)

        # Add to summary dataframe
        for col_name, values in gen_stats.items():
            summary_df[col_name] = values

        return summary_df

    def do_plot(self, variantDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile,
                metadata):
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

        # Extract the complex IDs:
        complex_ids = sim_data.process.complexation.ids_complexes
        print(f"Found {len(complex_ids)} complexes in simulation data")

        # Validate and get complexes to analyze
        if len(COMPLEXES_OF_INTEREST) > 0:
            complexes_to_analyze = self.check_validity_and_get_compartment(sim_data,
                                                                           COMPLEXES_OF_INTEREST)
            print(f"Analyzing {len(complexes_to_analyze)} complexes of interest")
        else:
            print(
                "WARNING: No complexes specified. Generating summary for all complexes but no individual plots.")
            complexes_to_analyze = []

        # ===== COMPLEX-SPECIFIC ANALYSIS (only for complexes of interest) =====
        if len(complexes_to_analyze) > 0:
            print(f"Analyzing {len(complexes_to_analyze)} complexes of interest...")
            self._analyze_complexes_of_interest(
                complexes_to_analyze,
                seeds,
                generations,
                plotOutDir
            )

        # Create a simple summary plot
        plt.figure(figsize=(12, 8))
        summary_text = f'Analysis Complete\n'
        summary_text += f'Total complexes: {len(complex_ids)}\n'
        summary_text += f'Complexes analyzed in detail: {len(complexes_to_analyze)}\n'
        summary_text += f'Seeds: {len(seeds)}\n'
        summary_text += f'See:\n  - complex_summary.csv for all complex statistics\n'
        summary_text += f'  - seed_summaries/ for per-seed metadata\n'
        if len(complexes_to_analyze) > 0:
            summary_text += f'  - complex_dataframes/ for detailed time-series data\n'
            summary_text += f'  - complex_plots/ for visualizations'

        plt.text(0.5, 0.5, summary_text, ha='center', va='center', fontsize=12, family='monospace')
        plt.axis('off')
        plt.tight_layout()
        exportFigure(plt, plotOutDir, plotOutFileName, metadata)
        plt.close('all')

        print("Analysis complete!")

        # ===== SUMMARY STATISTICS ACROSS ALL SEEDS AND TIMEPOINTS =====
        # Save the summary dataframe last since it takes a while:
        print("Generating summary statistics for all complexes...")
        self.total_generations = generations  # Store for use in other methods
        summary_df = self._generate_summary_statistics(cell_paths, complex_ids, seeds)
        summary_df.to_csv(os.path.join(plotOutDir, 'complex_summary.csv'), index=False)
        print(f"Saved summary for {len(complex_ids)} complexes")

        # ===== SEED-SPECIFIC METADATA =====
        #print("Saving seed-specific metadata...")
        #self._save_seed_metadata(seeds, all_cells, plotOutDir)

if __name__ == '__main__':
    Plot().cli()
