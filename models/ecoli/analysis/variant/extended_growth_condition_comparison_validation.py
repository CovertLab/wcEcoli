"""
Analysis plot for extended growth condition comparison validation.
"""
from __future__ import absolute_import, division, print_function


import os

import numpy as np
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from six.moves import cPickle, range

from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from wholecell.io.tablereader import TableReader
from wholecell.utils import units
from unum import Unum

from wholecell.utils.sparkline import whitePadSparklineAxis
from wholecell.analysis.analysis_tools import exportFigure
from models.ecoli.analysis import variantAnalysisPlot


FONT_SIZE=9
LEGEND_FONT_SIZE = 6


class Plot(variantAnalysisPlot.VariantAnalysisPlot):
    def do_plot(self, inputDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
        """
        Creates a validation plot for the variation of various macromolecular parameters with growth rate.

        Validation data derived from Dennis & Bremmer (2021).
        """
        ap = AnalysisPaths(inputDir, variant_plot=True)

        validation_data = cPickle.load(open(validationDataFile, "rb"))
        val_doubling_time_units = validation_data.macromolecular_growth_rate_modulation.doubling_time
        val_PRD_per_mass = validation_data.macromolecular_growth_rate_modulation.PRD_per_mass
        val_dry_mass_units = validation_data.macromolecular_growth_rate_modulation.mass_per_cell
        val_protein_mass_per_cell_units = validation_data.macromolecular_growth_rate_modulation.protein_per_cell_ug
        val_RNA_mass_per_cell_units = validation_data.macromolecular_growth_rate_modulation.RNA_per_cell_ug
        val_DNA_mass_per_cell_units = validation_data.macromolecular_growth_rate_modulation.DNA_per_cell_ug
        val_PRD_per_cell_units = validation_data.macromolecular_growth_rate_modulation.PRD_per_cell
        val_total_RNA_stable_fraction = validation_data.macromolecular_growth_rate_modulation.total_RNA_stable_fraction
        val_stable_RNA_tRNA_fraction = validation_data.macromolecular_growth_rate_modulation.stable_RNA_tRNA_fraction

        val_doubling_time = np.array([])
        val_dry_mass = np.array([])
        val_protein_mass_per_cell = np.array([])
        val_RNA_mass_per_cell = np.array([])
        val_DNA_mass_per_cell = np.array([])
        val_PRD_per_cell = np.array([])
        for i in val_doubling_time_units:
            val_doubling_time = np.append(val_doubling_time, i.asNumber(units.min))
        for i in val_dry_mass_units:
            val_dry_mass = np.append(val_dry_mass, i.asNumber(units.fg))
        for i in val_protein_mass_per_cell_units:
            val_protein_mass_per_cell = np.append(val_protein_mass_per_cell, i.asNumber(units.fg))
        for i in val_RNA_mass_per_cell_units:
            val_RNA_mass_per_cell = np.append(val_RNA_mass_per_cell, i.asNumber(units.fg))
        for i in val_DNA_mass_per_cell_units:
            val_DNA_mass_per_cell = np.append(val_DNA_mass_per_cell, i.asNumber(units.fg))
        for i in val_PRD_per_cell_units:
            val_PRD_per_cell = np.append(val_PRD_per_cell, i.asNumber(units.fg))

        sim_doubling_time = np.zeros(ap.n_variant)
        sim_PRD_per_mass = np.zeros(ap.n_variant)
        sim_dry_mass = np.zeros(ap.n_variant)
        sim_protein_mass = np.zeros(ap.n_variant)
        sim_RNA_mass = np.zeros(ap.n_variant)
        sim_DNA_mass = np.zeros(ap.n_variant)
        sim_PRD_per_cell = np.zeros(ap.n_variant)
        sim_total_RNA_stable_fraction = np.zeros(ap.n_variant)
        sim_stable_RNA_tRNA_fraction = np.zeros(ap.n_variant)

        sim_PRD_per_mass_std = np.zeros(ap.n_variant)
        sim_dry_mass_std = np.zeros(ap.n_variant)
        sim_protein_mass_std = np.zeros(ap.n_variant)
        sim_RNA_mass_std = np.zeros(ap.n_variant)
        sim_DNA_mass_std = np.zeros(ap.n_variant)
        sim_PRD_per_cell_std = np.zeros(ap.n_variant)
        sim_total_RNA_stable_fraction_std = np.zeros(ap.n_variant)
        sim_stable_RNA_tRNA_fraction_std = np.zeros(ap.n_variant)

        variants = ap.get_variants()

        for varIdx in range(ap.n_variant):
            variant = variants[varIdx]
            cells = ap.get_cells(variant=[variant])
            try:
                sim_data = cPickle.load(open(ap.get_variant_kb(variant), 'rb'))
            except Exception as e:
                print("Couldn't load sim_data object. Exiting.", e)
                return

            doublingTime = np.zeros(len(cells), dtype=Unum)
            proteinMass = np.zeros(len(cells))
            mRnaMass = np.zeros(len(cells))
            tRnaMass = np.zeros(len(cells))
            rRnaMass = np.zeros(len(cells))
            rnaMass = np.zeros(len(cells))
            dnaMass = np.zeros(len(cells))
            dryMass = np.zeros(len(cells))
            prdPerMass = np.zeros(len(cells))
            prdPerCell = np.zeros(len(cells))
            totalRnaStableFraction = np.zeros(len(cells))
            stableRnatRnaFraction = np.zeros(len(cells))

            for idx, simDir in enumerate(cells):
                simOutDir = os.path.join(simDir, "simOut")
                try:
                    main_reader = TableReader(os.path.join(simOutDir, "Main"))
                    time = main_reader.readColumn("time")
                except Exception as e:
                    print('Error with data for %s: %s' % (simDir, e))
                    continue

                doublingTime[idx] = (time[-1] - time[0]) * units.s
                mass = TableReader(os.path.join(simOutDir, "Mass"))
                proteinMass[idx] = mass.readColumn("proteinMass").mean()
                mRnaMass[idx] = mass.readColumn("mRnaMass").mean()
                tRnaMass[idx] = mass.readColumn("tRnaMass").mean()
                rRnaMass[idx] = mass.readColumn("rRnaMass").mean()
                rnaMass[idx] = mass.readColumn("rnaMass").mean()
                dnaMass[idx] = mass.readColumn("dnaMass").mean()
                dryMass[idx] = mass.readColumn("dryMass").mean()
                prdPerMass[idx] = (proteinMass[idx] + rnaMass[idx] + dnaMass[idx]) / dryMass[idx]
                prdPerCell[idx] = proteinMass[idx] + rnaMass[idx] + dnaMass[idx]
                totalRnaStableFraction[idx] = (tRnaMass[idx] + rRnaMass[idx]) / rnaMass[idx]
                stableRnatRnaFraction[idx] = tRnaMass[idx] / (tRnaMass[idx] + rRnaMass[idx])

            sim_doubling_time[varIdx] = np.mean(doublingTime).asNumber(units.min)
            sim_PRD_per_mass[varIdx] = prdPerMass.mean()
            sim_dry_mass[varIdx] = np.mean(dryMass)
            sim_protein_mass[varIdx] = np.mean(proteinMass)
            sim_RNA_mass[varIdx] = np.mean(rnaMass)
            sim_DNA_mass[varIdx] = np.mean(dnaMass)
            sim_PRD_per_cell[varIdx] = np.mean(prdPerCell)
            sim_total_RNA_stable_fraction[varIdx] = totalRnaStableFraction.mean()
            sim_stable_RNA_tRNA_fraction[varIdx] = stableRnatRnaFraction.mean()

            sim_PRD_per_mass_std[varIdx] = prdPerMass.std()
            sim_dry_mass_std[varIdx] = dryMass.std()
            sim_protein_mass_std[varIdx] = proteinMass.std()
            sim_RNA_mass_std[varIdx] = rnaMass.std()
            sim_DNA_mass_std[varIdx] = dnaMass.std()
            sim_PRD_per_cell_std[varIdx] = prdPerCell.std()
            sim_total_RNA_stable_fraction_std[varIdx] = totalRnaStableFraction.std()
            sim_stable_RNA_tRNA_fraction_std[varIdx] = stableRnatRnaFraction.std()


            with PdfPages(os.path.join(plotOutDir, plotOutFileName)) as pdf:
                fig = plt.figure()
                fig.set_figwidth(5)
                fig.set_figheight(5)

                ax0 = plt.subplot2grid((2, 2), (0, 0))
                ax1 = plt.subplot2grid((2, 2), (0, 1))
                ax2 = plt.subplot2grid((2, 2), (1, 0))
                ax3 = plt.subplot2grid((2, 2), (1, 1))
                lines = {'linestyle': 'dashed'}
                plt.rc('lines', **lines)

                ax0.errorbar(sim_doubling_time[np.argsort(sim_doubling_time)[::-1]],
                             sim_PRD_per_mass[np.argsort(sim_doubling_time)[::-1]],
                             yerr=sim_PRD_per_mass_std[np.argsort(sim_doubling_time)[::-1]], label="Simulation",
                             color="black", fmt='', marker='o', markersize=2, linewidth=0.5)
                ax0.errorbar(val_doubling_time[np.argsort(val_doubling_time)[::-1]],
                             np.array(val_PRD_per_mass)[np.argsort(val_doubling_time)[::-1]],
                             yerr=0,
                             label="Bremer & Dennis 2021", color="blue", marker='o', markersize=2,
                             linewidth=0.5)  # , markeredgecolor="blue")
                ax0.set_title("PRD per mass", fontsize=FONT_SIZE)
                ax0.set_ylim([0, 1])

                ax1.errorbar(sim_doubling_time[np.argsort(sim_doubling_time)[::-1]],
                             sim_dry_mass[np.argsort(sim_doubling_time)[::-1]],
                             yerr=sim_dry_mass_std[np.argsort(sim_doubling_time)[::-1]], label="Simulation",
                             color="black", fmt='', marker='o', markersize=2, linewidth=0.5)
                ax1.errorbar(val_doubling_time[np.argsort(val_doubling_time)[::-1]],
                             np.array(val_dry_mass)[np.argsort(val_doubling_time)[::-1]],
                             yerr=0,
                             label="Bremer & Dennis 2021", color="blue", marker='o', markersize=2,
                             linewidth=0.5)  # , markeredgecolor="blue")
                ax1.set_title("Dry mass per cell (fg)", fontsize=FONT_SIZE)

                # In the paper, the stable fraction of total RNA is calculated by estimating the fraction of total RNA that is mRNA, then subtracting from 1.
                # In our model, unstable RNAs also include the miscellaneous RNAs, not accounted for in the mRNA count. So we use the more accurate measure of tRNA + rRNA.
                # The difference between the two methods is quite small, at most 0.01.
                ax2.errorbar(sim_doubling_time[np.argsort(sim_doubling_time)[::-1]],
                             sim_total_RNA_stable_fraction[np.argsort(sim_doubling_time)[::-1]],
                             yerr=sim_total_RNA_stable_fraction_std[np.argsort(sim_doubling_time)[::-1]],
                             label="Simulation",
                             color="black", fmt='', marker='o', markersize=2, linewidth=0.5)
                ax2.errorbar(val_doubling_time[np.argsort(val_doubling_time)[::-1]],
                             np.array(val_total_RNA_stable_fraction)[np.argsort(val_doubling_time)[::-1]],
                             yerr=0,
                             label="Bremer & Dennis 2021", color="blue", marker='o', markersize=2,
                             linewidth=0.5)  # , markeredgecolor="blue")
                ax2.set_title("Stable Fraction of Total RNA", fontsize=FONT_SIZE)
                ax2.set_ylim([0.8, 1])

                ax3.errorbar(sim_doubling_time[np.argsort(sim_doubling_time)[::-1]],
                             sim_stable_RNA_tRNA_fraction[np.argsort(sim_doubling_time)[::-1]],
                             yerr=sim_stable_RNA_tRNA_fraction_std[np.argsort(sim_doubling_time)[::-1]],
                             label="Simulation",
                             color="black", fmt='', marker='o', markersize=2, linewidth=0.5)
                ax3.errorbar(val_doubling_time[np.argsort(val_doubling_time)[::-1]],
                             np.array(val_stable_RNA_tRNA_fraction)[np.argsort(val_doubling_time)[::-1]],
                             yerr=0,
                             label="Bremer & Dennis 2021", color="blue", marker='o', markersize=2,
                             linewidth=0.5)  # , markeredgecolor="blue")
                ax3.set_title("tRNA Fraction of Stable RNA", fontsize=FONT_SIZE)
                ax3.set_ylim([0, 0.2])

                axes_list = [ax0, ax1, ax2, ax3]
                for axes in axes_list:
                    axes.legend(fontsize=LEGEND_FONT_SIZE)
                    axes.set_xlabel("Doubling time (min)", fontsize=FONT_SIZE)

                plt.tight_layout()
                pdf.savefig()
                plt.close()

                ax4 = plt.subplot2grid((2, 2), (0, 0))
                ax5 = plt.subplot2grid((2, 2), (0, 1))
                ax6 = plt.subplot2grid((2, 2), (1, 0))
                ax7 = plt.subplot2grid((2, 2), (1, 1))

                ax4.errorbar(sim_doubling_time[np.argsort(sim_doubling_time)[::-1]],
                             sim_protein_mass[np.argsort(sim_doubling_time)[::-1]],
                             yerr=sim_protein_mass_std[np.argsort(sim_doubling_time)[::-1]], label="Simulation",
                             color="black", fmt='', marker='o', markersize=2, linewidth=0.5)
                ax4.errorbar(val_doubling_time[np.argsort(val_doubling_time)[::-1]],
                             np.array(val_protein_mass_per_cell)[np.argsort(val_doubling_time)[::-1]],
                             yerr=0,
                             label="Bremer & Dennis 2021", color="blue", marker='o', markersize=2,
                             linewidth=0.5)  # , markeredgecolor="blue")
                ax4.set_title("Protein mass per cell (fg)", fontsize=FONT_SIZE)

                ax5.errorbar(sim_doubling_time[np.argsort(sim_doubling_time)[::-1]],
                             sim_RNA_mass[np.argsort(sim_doubling_time)[::-1]],
                             yerr=sim_RNA_mass_std[np.argsort(sim_doubling_time)[::-1]], label="Simulation",
                             color="black", fmt='', marker='o', markersize=2, linewidth=0.5)
                ax5.errorbar(val_doubling_time[np.argsort(val_doubling_time)[::-1]],
                             np.array(val_RNA_mass_per_cell)[np.argsort(val_doubling_time)[::-1]],
                             yerr=0,
                             label="Bremer & Dennis 2021", color="blue", marker='o', markersize=2,
                             linewidth=0.5)  # , markeredgecolor="blue")
                ax5.set_title("RNA mass per cell (fg)", fontsize=FONT_SIZE)

                ax6.errorbar(sim_doubling_time[np.argsort(sim_doubling_time)[::-1]],
                             sim_DNA_mass[np.argsort(sim_doubling_time)[::-1]],
                             yerr=sim_DNA_mass_std[np.argsort(sim_doubling_time)[::-1]], label="Simulation",
                             color="black", fmt='', marker='o', markersize=2, linewidth=0.5)
                ax6.errorbar(val_doubling_time[np.argsort(val_doubling_time)[::-1]],
                             np.array(val_DNA_mass_per_cell)[np.argsort(val_doubling_time)[::-1]],
                             yerr=0,
                             label="Bremer & Dennis 2021", color="blue", marker='o', markersize=2,
                             linewidth=0.5)  # , markeredgecolor="blue")
                ax6.set_title("DNA mass per cell (fg)", fontsize=FONT_SIZE)

                ax7.errorbar(sim_doubling_time[np.argsort(sim_doubling_time)[::-1]],
                             sim_PRD_per_cell[np.argsort(sim_doubling_time)[::-1]],
                             yerr=sim_PRD_per_cell_std[np.argsort(sim_doubling_time)[::-1]], label="Simulation",
                             color="black", fmt='', marker='o', markersize=2, linewidth=0.5)
                ax7.errorbar(val_doubling_time[np.argsort(val_doubling_time)[::-1]],
                             np.array(val_PRD_per_cell)[np.argsort(val_doubling_time)[::-1]],
                             yerr=0,
                             label="Bremer & Dennis 2021", color="blue", marker='o', markersize=2,
                             linewidth=0.5)  # , markeredgecolor="blue")
                ax7.set_title("PRD per cell (fg)", fontsize=FONT_SIZE)

                axes_list = [ax4, ax5, ax6, ax7]
                for axes in axes_list:
                    axes.legend(fontsize=LEGEND_FONT_SIZE)
                    axes.set_xlabel("Doubling time (min)", fontsize=FONT_SIZE)

                plt.tight_layout()
                pdf.savefig()
                plt.close()


if __name__ == "__main__":
    Plot().cli()