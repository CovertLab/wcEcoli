### TODO: Merge into new_gene_translation_efficiency_heatmaps.py (kept
# separate now for convenience of running on large simulation set)


"""
Plot one value per index via heatmap for
new_gene_expression_and_translation_efficiency variant.

Plots:
- Average cell volume
- Average cell mass
- Average dry mass
- Average cell protein mass
- Average cell mRNA mass

TODO:
- Accomodate more than one new gene
"""

import numpy as np
from matplotlib import pyplot as plt

from models.ecoli.analysis import variantAnalysisPlot
from wholecell.analysis.analysis_tools import exportFigure, \
    read_stacked_columns, stacked_cell_threshold_mask

import os.path

# 1 to exclude cells that took full MAX_CELL_LENGTH, 0 otherwise
exclude_timeout_cells = 1

"""
1 to plot early (before MIN_LATE_CELL_INDEX), and late generations in
addition to all generations
"""
exclude_early_gens = 1

FONT_SIZE = 9
MAX_VARIANT = 55  # do not include any variant >= this index
MAX_CELL_INDEX = 16  # do not include any generation >= this index

"""
Count number of sims that reach this generation (remember index 7 
corresponds to generation 8)
"""
COUNT_INDEX = 15

"""
generations before this may not be representative of dynamics 
due to how they are initialized
"""
MIN_LATE_CELL_INDEX = 4

MAX_CELL_LENGTH = 36000
if (exclude_timeout_cells == 0):
    MAX_CELL_LENGTH += 1000000


class Plot(variantAnalysisPlot.VariantAnalysisPlot):
    ### TODO: move to analysis_tools
    def heatmap(self, ax, mask, data, completion_data, xlabel, ylabel, xlabels,
                ylabels, title, textsize="medium"):
        im = ax.imshow(data, cmap="GnBu")
        ax.set_xticks(np.arange(len(xlabels)))
        ax.set_xticklabels(xlabels)
        ax.set_yticks(np.arange(len(
            ylabels)))
        ax.set_yticklabels(ylabels)
        plt.setp(ax.get_xticklabels(), rotation=45, ha="right",
                 rotation_mode="anchor")
        for i in range(len(ylabels)):
            for j in range(len(xlabels)):
                if mask[i, j]:
                    col = "k"
                    if completion_data[i, j] < 0.9:
                        col = "r"
                    text = ax.text(j, i, data[i, j],
                                   ha="center", va="center", color=col,
                                   fontsize=textsize)
        ax.set_xlabel(xlabel, fontsize=FONT_SIZE)
        ax.set_ylabel(ylabel, fontsize=FONT_SIZE)
        ax.set_title(title)

    def do_plot(self, inputDir, plotOutDir, plotOutFileName, simDataFile,
                validationDataFile, metadata):

        print("Running analysis script with exclude_timeout_cells=",
              exclude_timeout_cells, " and exclude_early_gens=",
              exclude_early_gens)

        # Map variant indices to expression factors and translation efficiency 
        # values
        if 'new_gene_expression_factors' not in metadata or \
                'new_gene_translation_efficiency_values' not in metadata:
            print("This plot is intended to be run on simulations where the"
                  " new gene expression-translation efficiency variant was "
                  "enabled, but no parameters for this variant were found.")

        new_gene_expression_factors = metadata[
            'new_gene_expression_factors']
        new_gene_translation_efficiency_values = metadata[
            'new_gene_translation_efficiency_values']

        separator = len(new_gene_translation_efficiency_values)

        variants = self.ap.get_variants()
        variant_index_to_values = {}
        variant_index_to_list_indices = {}
        variant_mask = np.zeros((  # Track whether we ran this sim
            len(new_gene_translation_efficiency_values),
            len(new_gene_expression_factors)), dtype=bool)

        for index in variants:
            if index >= MAX_VARIANT:
                continue

            if index == 0:
                expression_list_index = 0
                trl_eff_list_index = len(
                    new_gene_translation_efficiency_values) - 1
                expression_variant_index = 0
                # Note: this value should not matter since gene is knocked out
                trl_eff_value = 0
            else:
                trl_eff_list_index = index % separator
                if trl_eff_list_index == 0:
                    expression_list_index = index // separator
                else:
                    expression_list_index = index // separator + 1

                expression_variant_index = new_gene_expression_factors[
                    expression_list_index]
                trl_eff_value = new_gene_translation_efficiency_values[
                    trl_eff_list_index]
            variant_index_to_values[index] = np.array([
                expression_variant_index, trl_eff_value])
            variant_index_to_list_indices[index] = np.array([
                expression_list_index, trl_eff_list_index])
            variant_mask[trl_eff_list_index, expression_list_index] = True

        # Create data structures that we will use for the heatmaps
        doubling_times_heatmap = np.zeros((3,
            len(new_gene_translation_efficiency_values),
            len(new_gene_expression_factors))) - 1
        completed_gens_heatmap = np.zeros((1,
            len(new_gene_translation_efficiency_values),
            len(new_gene_expression_factors)))
        # TODO: Expand to Accomodate Multiple New Genes
        avg_volume_heatmap = np.zeros((3,
            len(new_gene_translation_efficiency_values),
            len(new_gene_expression_factors))) - 1
        avg_mass_heatmap = np.zeros((3,
            len(new_gene_translation_efficiency_values),
            len(new_gene_expression_factors))) - 1
        avg_dry_mass_heatmap = np.zeros((3,
             len(new_gene_translation_efficiency_values),
             len(new_gene_expression_factors))) - 1
        avg_mRNA_mass_heatmap = np.zeros((3,
             len(new_gene_translation_efficiency_values),
             len(new_gene_expression_factors))) - 1
        avg_monomer_mass_heatmap = np.zeros((3,
             len(new_gene_translation_efficiency_values),
             len(new_gene_expression_factors))) - 1

        # Data extraction
        print("---Data Extraction---")
        doubling_times = {}
        reached_count_gen = {}
        generations = {}

        variants = self.ap.get_variants()
        min_variant = min(variants)
        for variant in variants:

            if variant >= MAX_VARIANT:
                continue

            print("Variant: ", variant)
            all_cells = self.ap.get_cells(variant=[variant],
                                          only_successful=True)
            if len(all_cells) == 0:
                continue

            exclude_timeout_cell_mask = stacked_cell_threshold_mask(
                all_cells, 'Main', 'time', MAX_CELL_LENGTH,
                fun=lambda x: (x[-1] - x[0]) / 60.).squeeze()
            all_cells_gens = np.array([int(os.path.basename(os.path.dirname(
                cell_path))[-6:]) for cell_path in all_cells])[
                exclude_timeout_cell_mask]
            generations[variant] = all_cells_gens

            # Doubling times
            dt = read_stacked_columns(all_cells, 'Main', 'time',
                                      fun=lambda x: (x[-1] - x[0]) / 60.)
            doubling_times[variant] = dt[exclude_timeout_cell_mask]

            # Count the number of simulations that reach gen COUNT_INDEX + 1
            num_count_gen = len(self.ap.get_cells(variant=[variant],
                                                  generation=[COUNT_INDEX],
                                                  only_successful=True))
            num_zero_gen = len(self.ap.get_cells(variant=[variant],
                                                 generation=[0],
                                                 only_successful=True))
            reached_count_gen[variant] = num_count_gen / num_zero_gen

            # Total volume
            avg_volume = read_stacked_columns(
                all_cells, 'Mass', 'cellVolume', fun=lambda x: np.mean(x))
            avg_volume = avg_volume[exclude_timeout_cell_mask]

            # Total mass
            avg_mass = read_stacked_columns(
                all_cells, 'Mass', 'cellMass', fun=lambda x: np.mean(x))
            avg_mass = avg_mass[exclude_timeout_cell_mask]

            # Total dry mass
            avg_dry_mass = read_stacked_columns(
                all_cells, 'Mass', 'dryMass', fun=lambda x: np.mean(x))
            avg_dry_mass = avg_dry_mass[exclude_timeout_cell_mask]

            # Total mRNA mass
            avg_mRNA_mass = read_stacked_columns(
                all_cells, 'Mass', 'mRnaMass', fun=lambda x: np.mean(x))
            avg_mRNA_mass = avg_mRNA_mass[exclude_timeout_cell_mask]

            # Total protein mass
            avg_monomer_mass = read_stacked_columns(all_cells, 'Mass',
                'proteinMass', fun=lambda x: np.mean(x))
            avg_monomer_mass = avg_monomer_mass[exclude_timeout_cell_mask]

            # Add values to heatmap data structures
            exp_index, trl_eff_index = variant_index_to_list_indices[variant]
            doubling_times_heatmap[0, trl_eff_index, exp_index] = round(
                np.mean(doubling_times[variant]))
            completed_gens_heatmap[0, trl_eff_index, exp_index] = \
                round(reached_count_gen[variant], 2)
            avg_volume_heatmap[0, trl_eff_index, exp_index] = round(np.mean(
                avg_volume), 2)
            avg_mass_heatmap[0, trl_eff_index, exp_index] = round(
                np.mean(avg_mass), 2)
            avg_dry_mass_heatmap[0, trl_eff_index, exp_index] = round(
                np.mean(avg_dry_mass), 2)
            avg_mRNA_mass_heatmap[0, trl_eff_index, exp_index] = round(
                np.mean(avg_mRNA_mass), 2)
            avg_monomer_mass_heatmap[0, trl_eff_index, exp_index] = round(
                np.mean(avg_monomer_mass), 2)

            if exclude_early_gens == 1:
                # Add early gen values to the heatmap structure
                early_cell_mask = generations[variant] < MIN_LATE_CELL_INDEX
                if len(early_cell_mask) == 1:
                    early_cell_mask = early_cell_mask[0]

                avg_volume_heatmap[1, trl_eff_index, exp_index] = round(
                    np.mean(avg_volume[early_cell_mask]), 2)
                avg_mass_heatmap[1, trl_eff_index, exp_index] = round(
                    np.mean(avg_mass[early_cell_mask]), 2)
                avg_dry_mass_heatmap[1, trl_eff_index, exp_index] = round(
                    np.mean(avg_dry_mass[early_cell_mask]), 2)
                avg_mRNA_mass_heatmap[1, trl_eff_index, exp_index] = round(
                    np.mean(avg_mRNA_mass[early_cell_mask]), 2)
                avg_monomer_mass_heatmap[1, trl_eff_index, exp_index] = round(
                    np.mean(avg_monomer_mass[early_cell_mask]), 2)

                # Add late gen values to the heatmap structure
                late_cell_mask = np.logical_and((generations[variant] >=
                     MIN_LATE_CELL_INDEX), (generations[variant] < MAX_CELL_INDEX))
                if len(late_cell_mask) == 1:
                    late_cell_mask = late_cell_mask[0]
                if sum(late_cell_mask) != 0:
                    avg_volume_heatmap[2, trl_eff_index, exp_index] = round(
                        np.mean(avg_volume[late_cell_mask]), 2)
                    avg_mass_heatmap[2, trl_eff_index, exp_index] = round(
                        np.mean(avg_mass[late_cell_mask]), 2)
                    avg_dry_mass_heatmap[2, trl_eff_index, exp_index] = round(
                        np.mean(avg_dry_mass[late_cell_mask]), 2)
                    avg_mRNA_mass_heatmap[2, trl_eff_index, exp_index] = round(
                        np.mean(avg_mRNA_mass[late_cell_mask]), 2)
                    avg_monomer_mass_heatmap[2, trl_eff_index, exp_index] = round(
                        np.mean(avg_monomer_mass[late_cell_mask]), 2)

        # Plotting
        print("---Plotting---")
        plot_descr = ["_all_gens"]
        if exclude_early_gens == 1:
            plot_descr += ["_early_gens", "_late_gens"]

        for j in range(len(plot_descr)):
            # Volume
            fig, ax = plt.subplots(1, 1, figsize=(10, 5))
            self.heatmap(ax, variant_mask,
                         avg_volume_heatmap[j, :, :],
                         completed_gens_heatmap[0, :, :],
                         "Expression Variant",
                         "Translation Efficiency Value (Normalized)",
                         new_gene_expression_factors,
                         new_gene_translation_efficiency_values,
                         "Cell Volume")
            fig.tight_layout()
            plt.show()
            exportFigure(plt, plotOutDir,
                         'volume_heatmap' + plot_descr[j])

            # Mass
            fig, ax = plt.subplots(1, 1, figsize=(10, 5))
            self.heatmap(ax, variant_mask,
                         avg_mass_heatmap[j, :, :],
                         completed_gens_heatmap[0, :, :],
                         "Expression Variant",
                         "Translation Efficiency Value (Normalized)",
                         new_gene_expression_factors,
                         new_gene_translation_efficiency_values,
                         "Cell Mass", "x-small")
            fig.tight_layout()
            plt.show()
            exportFigure(plt, plotOutDir,
                         'mass_heatmap' + plot_descr[j])

            # Dry Mass
            fig, ax = plt.subplots(1, 1, figsize=(10, 5))
            self.heatmap(ax, variant_mask,
                         avg_dry_mass_heatmap[j, :, :],
                         completed_gens_heatmap[0, :, :],
                         "Expression Variant",
                         "Translation Efficiency Value (Normalized)",
                         new_gene_expression_factors,
                         new_gene_translation_efficiency_values,
                         "Dry Mass", "x-small")
            fig.tight_layout()
            plt.show()
            exportFigure(plt, plotOutDir,
                         'dry_mass_heatmap' + plot_descr[j])

            # mRNA Mass
            fig, ax = plt.subplots(1, 1, figsize=(10, 5))
            self.heatmap(ax, variant_mask,
                         avg_mRNA_mass_heatmap[j, :, :],
                         completed_gens_heatmap[0, :, :],
                         "Expression Variant",
                         "Translation Efficiency Value (Normalized)",
                         new_gene_expression_factors,
                         new_gene_translation_efficiency_values,
                         "mRNA Mass")
            fig.tight_layout()
            plt.show()
            exportFigure(plt, plotOutDir,
                         'mRNA_mass_heatmap' + plot_descr[j])

            # Protein Mass
            fig, ax = plt.subplots(1, 1, figsize=(10, 5))
            self.heatmap(ax, variant_mask,
                         avg_monomer_mass_heatmap[j, :, :],
                         completed_gens_heatmap[0, :, :],
                         "Expression Variant",
                         "Translation Efficiency Value (Normalized)",
                         new_gene_expression_factors,
                         new_gene_translation_efficiency_values,
                         "Protein Mass", "x-small")
            fig.tight_layout()
            plt.show()
            exportFigure(plt, plotOutDir,
                         'monomer_mass_heatmap' + plot_descr[j])

            plt.close('all')


if __name__ == "__main__":
    Plot().cli()