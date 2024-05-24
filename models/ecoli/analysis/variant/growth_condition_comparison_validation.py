"""
Compare various cell macromolecular properties across different growth rates
with validation data from Dennis & Bremmer (2021)
https://journals.asm.org/doi/epub/10.1128/ecosal.5.2.3
"""

import os

import numpy as np
from matplotlib import pyplot as plt
import pickle

from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from wholecell.io.tablereader import TableReader
from wholecell.utils import units

from models.ecoli.analysis import variantAnalysisPlot
from wholecell.analysis.analysis_tools import (exportFigure, read_stacked_columns)

FIG_SIZE = (30, 40)

NCOLS = 6
NROWS = 8

LABEL_FONT_SIZE = 9
LEGEND_FONT_SIZE = 6

SIM_RAW_SCATTER_STYLE = dict(
    s = 0.2,
    c = 'black',
    marker = 'o'
)

SIM_MEAN_LINE_STYLE = dict(
    color = "black",
    marker = 'o',
    markersize = 1,
    linewidth = 0.5,
    linestyle = 'dashed'
)

VAL_LINE_STYLE = dict(
    color = "blue",
    marker = 'o',
    markersize = 1,
    linewidth = 0.5,
    linestyle = 'dashed'
)

# Plotted units of data
mass_units = units.fg
time_units = units.min
mole_units = units.mol

# Average E. coli monomer masses as cited in Dennis & Bremmer (2021)
AVG_AA_MASS = 108.0 * units.g / units.mol
AVG_RNA_MASS = 324.0 * units.g / units.mol

# Placeholder for variables where we don't set a y_lim
null_ylim = (-1, -1)

class Plot(variantAnalysisPlot.VariantAnalysisPlot):
    def plot_lines(self, fig, position, sim_raw_x, sim_raw_y, sim_mean_x, sim_mean_y, sim_std, val_x, val_y, title, y_lim=None):
        # Sort plotted variables in order of increasing doubling time
        sim_raw_y = sim_raw_y[np.argsort(sim_raw_x)][::-1]
        sim_raw_x = sim_raw_x[np.argsort(sim_raw_x)][::-1]
        sim_mean_y = sim_mean_y[np.argsort(sim_mean_x)][::-1]
        sim_std = sim_std[np.argsort(sim_mean_x)][::-1]
        sim_mean_x = sim_mean_x[np.argsort(sim_mean_x)][::-1]
        val_y = val_y[np.argsort(val_x)][::-1]
        val_x = val_x[np.argsort(val_x)][::-1]

        # Make plot
        ax = plt.subplot2grid((NROWS, NCOLS), position, rowspan=1, colspan=1, fig=fig)
        ax.scatter(sim_raw_x, sim_raw_y, **SIM_RAW_SCATTER_STYLE)
        ax.errorbar(sim_mean_x, sim_mean_y, yerr=sim_std, label="Simulation",
                    **SIM_MEAN_LINE_STYLE)
        ax.errorbar(val_x, val_y, yerr=0, label="Bremer & Dennis 2021",
                    **VAL_LINE_STYLE)
        if not np.array_equiv(y_lim, np.asarray(null_ylim)):
            ax.set_ylim(y_lim)
        ax.set_xlabel("Doubling time (min)", fontsize=LABEL_FONT_SIZE)
        ax.set_title(title, fontsize=LABEL_FONT_SIZE)
        ax.legend(fontsize=LEGEND_FONT_SIZE)

    def do_plot(self, inputDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
        ap = AnalysisPaths(inputDir, variant_plot=True)

        # Tuples for all plots, consisting of tuple name, units to plot in,
        # title of plot, position (row, column) of plot, and top y-axis limit.
        plotted_variables = np.array([
            ("doubling_time", time_units, None, (0, 0), None),
            ("protein_per_mass", units.aa / mass_units, "Protein amino acids per dry mass", (0, 0), null_ylim),
            ("RNA_per_mass", units.nt / mass_units, "RNA nucleotides per dry mass", (0, 1), null_ylim),
            ("DNA_per_mass", 1 / mass_units, "Genome equivalents per dry mass", (0, 2), null_ylim),
            ("protein_per_genome", units.aa, "Protein amino acids per genome equivalent", (0, 3), null_ylim),
            ("RNA_per_genome", units.nt, "RNA nucleotides per genome equivalent", (0, 4), null_ylim),

            ("origins_per_cell", None, "Origins of replication per cell", (1, 0), null_ylim),
            ("termini_per_cell", None, "Chromosome termini per cell", (1, 1), null_ylim),
            ("replication_forks_per_cell", None, "Replication forks per cell", (1, 2), null_ylim),
            ("origins_per_genome", None, "Origins per genome equivalent", (1, 3), null_ylim),
            ("protein_per_origin", units.aa, "Protein amino acids per origin", (1, 4), null_ylim),

            ("DNA_per_cell_ug", mass_units, "DNA mass per cell (fg)", (2, 0), null_ylim),
            ("mass_per_cell", mass_units, "Dry mass per cell (fg)", (2, 1), null_ylim),
            ("protein_per_cell_ug", mass_units, "Protein mass per cell (fg)", (2, 2), null_ylim),
            ("RNA_per_cell_ug", mass_units, "RNA mass per cell (fg)", (2, 3), null_ylim),
            ("PRD_per_mass", None, "Protein + RNA + DNA mass per dry mass", (2, 4), (0, 1)),
            ("PRD_per_cell", mass_units, "Protein + RNA + DNA mass per cell (fg)", (2, 5), null_ylim),

            ("total_RNA_stable_fraction", None, "Stable Fraction of Total RNA", (3, 0), (0.8, 1)),
            ("stable_RNA_tRNA_fraction", None, "tRNA Fraction of Stable RNA", (3, 1), (0, 0.2)),
            ("ribosomes_per_cell", None, "Ribosomes per cell", (3, 2), null_ylim),
            ("r_prot_per_total_protein", None, "Ribosomal protein per total protein (amino acids)", (3, 3), (0, 0.5)),
            ("tRNA_per_cell", None, "tRNA counts per cell", (3, 4), null_ylim),
            ("peptide_elongation_rate", units.aa / units.s, "Ribosome elongation rate per second (aa/s)", (3, 5), null_ylim),

            ("RNAP_per_total_protein", None, "Protein fraction of RNAPs", (4, 0), (0, 0.02)),
            ("RNAP_per_cell", None, "Total RNAP cores per cell", (4, 1), null_ylim),
            ("RNAP_per_ribosome", None, "RNA polymerase cores per ribosome", (4, 2), (0, 0.4)),

            ("RNA_synth_stable_fraction", None, "Fraction of RNA synthesis that is stable RNAs", (5, 0), (0, 1)),
            ("stable_RNA_synthesis_per_cell", units.nt / units.s, "Cell-wide stable RNA synthesis rate (nt/sec)", (5, 1), null_ylim),
            ("mRNA_synthesis_per_cell", units.nt / units.s, "Cell-wide mRNA synthesis rate (nt/sec)", (5, 2), null_ylim),
            ("active_RNAP_synthesizing_stable_fraction", None, "Fraction of active RNAP synthesizing stable RNA", (5, 3), (0, 1)),
            ("RNAP_active_fraction", None, "Active fraction of RNAPs", (5, 4), (0, 1)),
            ("active_RNAP_per_cell", None, "Active RNAPs per cell", (5, 5), null_ylim),

            ("rrn_genes_per_cell", None, "Rrn genes per cell", (6, 0), null_ylim),
            ("rrn_genes_per_genome", None, "Rrn genes per genome", (6, 1), null_ylim),
            ("rrn_gene_initiation_rate", 1 / units.s, "Ribosomal gene initiation rate (inits/gene/s)", (6, 2), null_ylim),

            ("ppGpp_concn_per_mass", mole_units / mass_units, "ppGpp concentration per dry mass (mol/fg)", (7, 0), null_ylim),
            ("ppGpp_concn_per_protein", mole_units / units.aa, "ppGpp concentration per total protein (mol/amino acid)", (7, 1), null_ylim),
            ],
            dtype=[('name', 'U100'), ('units', units.Unum), ('plot_name', 'U100'), ('position', np.int32, (2,)), ('y_lim', np.float64, (2,))])

        # Extract validation data and strip the units while numerically
        # converting to the plotted units
        validation_data = self.read_pickle_file(validationDataFile)
        db_table = validation_data.macromolecular_growth_rate_modulation
        val_table = dict()
        for idx in range(len(plotted_variables)):
            var_name = plotted_variables[idx]['name']
            var_units = plotted_variables[idx]['units']
            val_table[var_name] = getattr(db_table, var_name).asNumber(var_units)

        # Initialize ndarrays for simulation values per cell, as well as means
        # and standard deviations within each variant
        sim_table = dict()
        sim_mean_table = dict()
        sim_std_table = dict()
        for var_name in plotted_variables['name']:
            sim_table[var_name] = np.array([])
            sim_mean_table[var_name] = np.zeros(ap.n_variant)
            sim_std_table[var_name] = np.zeros(ap.n_variant)

        variants = ap.get_variants()

        # Read universal table attributes and sim data from a test cell for use
        test_cell = ap.get_cells(variant=[variants[0]])[0]
        simOutDir = os.path.join(test_cell, "simOut")
        try:
            sim_data = self.read_pickle_file(self.ap.get_variant_kb(variants[0]))
        except Exception as e:
            print("Couldn't load sim_data object. Exiting.", e)
            return
        try:
            mass = TableReader(os.path.join(simOutDir, "Mass"))
        except Exception as e:
            print('Error with data for %s: %s' % (test_cell, e))

        # Get units that the Mass table uses
        protein_mass_units = getattr(units,
                                     mass.readAttribute("protein_units"))
        nucleic_acid_mass_units = getattr(units,
                                          mass.readAttribute("rna_units"))
        dry_mass_units = getattr(units, mass.readAttribute("cellDry_units"))

        # Create mask for rrn genes for counting rrn gene copy number
        # (note: even though the 'gene_ids' attribute is defined using
        # cistron_data in the RnaSynthProb listener, the cistron_data['is_rRNA']
        # doesn't work because it contains less elements than 'gene_ids' for
        # some reason)
        rna_reader = TableReader(os.path.join(simOutDir, "RnaSynthProb"))
        gene_ids = rna_reader.readAttribute('gene_ids')
        gene_data = sim_data.process.replication.gene_data
        gene_id_to_rna_id = {gene['name']: gene['cistron_id'] + '[c]' for gene in gene_data}
        molecule_groups = sim_data.molecule_groups
        rrna_ids = molecule_groups.s30_16s_rRNA + molecule_groups.s50_23s_rRNA + molecule_groups.s50_5s_rRNA
        is_rrn = [gene_id_to_rna_id[gene] in rrna_ids for gene in gene_ids]

        # Get masks and RNA lengths to calculate rrn initiation rate and
        # RNA synthesis rates
        rna_data = sim_data.process.transcription.rna_data
        is_rRNA = rna_data['is_rRNA']
        is_stable = [rna['is_rRNA'] or rna['is_tRNA'] for rna in rna_data]
        is_mRNA = rna_data['is_mRNA']
        stable_rna_lengths = rna_data[is_stable]['length'].asNumber()
        mRNA_lengths = rna_data[is_mRNA]['length'].asNumber()

        # Get indices/masks of various molecules
        uniqueMolecules = TableReader(os.path.join(simOutDir, "UniqueMoleculeCounts"))
        active_rnap_index = uniqueMolecules.readAttribute("uniqueMoleculeIds").index("active_RNAP")
        active_ribosome_index = uniqueMolecules.readAttribute("uniqueMoleculeIds").index(
            'active_ribosome')
        full_chromosome_index = uniqueMolecules.readAttribute("uniqueMoleculeIds").index(
            'full_chromosome')
        bulkMolecules = TableReader(os.path.join(simOutDir, "BulkMolecules"))
        molecule_ids = sim_data.molecule_ids
        rnap_id = molecule_ids.full_RNAP
        rnap_index = bulkMolecules.readAttribute('objectNames').index(rnap_id)
        inactive_ribosome_ids = [molecule_ids.s50_full_complex, molecule_ids.s30_full_complex]
        inactive_ribosome_mask = [molecule in inactive_ribosome_ids for molecule in
            bulkMolecules.readAttribute('objectNames')]
        uncharged_tRNA_ids = list(rna_data['id'][rna_data['is_tRNA']])
        charged_tRNA_ids = list(sim_data.process.transcription.charged_trna_names)
        bulk_tRNA_ids = uncharged_tRNA_ids + charged_tRNA_ids
        bulk_tRNA_mask = [molecule in bulk_tRNA_ids for molecule in
            bulkMolecules.readAttribute('objectNames')]
        ppGpp_id = molecule_ids.ppGpp
        ppGpp_index = bulkMolecules.readAttribute('objectNames').index(ppGpp_id)
        n_avogadro = sim_data.constants.n_avogadro
        counts_to_moles = (1 / n_avogadro).asNumber(mole_units)

        # Get r-protein mask and r-prot amino acid counts to calculate r-protein
        # per total protein
        ribosomal_protein_ids = molecule_groups.ribosomal_proteins
        monomer_data = sim_data.process.translation.monomer_data
        is_r_prot = [protein['id'] in ribosomal_protein_ids for protein in monomer_data]
        r_prot_aa_counts = np.sum(monomer_data[is_r_prot]['aa_counts'].asNumber(), axis=1)

        # Calculate mass of relevant molecules
        geq_mass = (sim_data.getter.get_mass(molecule_ids.full_chromosome) / n_avogadro).asNumber(
            mass_units)
        # Bremmer & Dennis (2021) only consider beta and beta' subunits when
        # calculating RNAP-related data
        rnap_core_subunits_ids = ['RPOB-MONOMER[c]', 'RPOC-MONOMER[c]']
        rnap_mass = (np.sum(sim_data.getter.get_masses(rnap_core_subunits_ids)) / n_avogadro).asNumber(mass_units)
        ribosome_subunit_masses = (sim_data.getter.get_masses(inactive_ribosome_ids) / n_avogadro).asNumber(mass_units)
        full_ribosome_mass = np.sum(ribosome_subunit_masses)
        ribosome_subunit_mass_fractions = ribosome_subunit_masses / full_ribosome_mass
        avg_aa_mass = (AVG_AA_MASS / n_avogadro).asNumber(mass_units)
        avg_rna_mass = (AVG_RNA_MASS / n_avogadro).asNumber(mass_units)

        # For each variant, extract the simulated values from each cell and
        # store together as an array. These arrays are compiled in a dictionary,
        # raw_sim_table.
        # TODO: implement function to take time average and accounting for
        #  cell age distribution instead of step averages
        # We do averages on a per cell basis, which may give differences from
        # population averages in ratio measurements
        for varIdx in range(ap.n_variant):
            variant = variants[varIdx]
            cells = ap.get_cells(variant=[variant])
            raw_sim_table = dict()
            raw_sim_table['doubling_time_sec'] = read_stacked_columns(
                cells, 'Main', 'time', fun=lambda x: x[-1] - x[0])
            raw_sim_table['doubling_time'] = raw_sim_table['doubling_time_sec'] * (units.s).asNumber(time_units)
            raw_sim_table['protein_per_cell_ug'] = read_stacked_columns(
                cells, 'Mass', 'proteinMass', fun=lambda x: x.mean(),
                remove_first=True) * protein_mass_units.asNumber(mass_units)
            raw_sim_table['RNA_per_cell_ug'] = read_stacked_columns(
                cells, 'Mass', 'rnaMass', fun=lambda x: x.mean(),
                remove_first=True) * nucleic_acid_mass_units.asNumber(mass_units)
            raw_sim_table['DNA_per_cell_ug'] = read_stacked_columns(
                cells, 'Mass', 'dnaMass', fun=lambda x: x.mean(),
                remove_first=True) * nucleic_acid_mass_units.asNumber(mass_units)
            raw_sim_table['mass_per_cell'] = read_stacked_columns(
                cells, 'Mass', 'dryMass', fun=lambda x: x.mean(),
                remove_first=True) * dry_mass_units.asNumber(mass_units)
            raw_sim_table['PRD_per_mass'] = (raw_sim_table['protein_per_cell_ug'] + raw_sim_table['RNA_per_cell_ug'] + raw_sim_table['DNA_per_cell_ug']) / raw_sim_table[
                'mass_per_cell']
            raw_sim_table['PRD_per_cell'] = raw_sim_table['protein_per_cell_ug'] + raw_sim_table['RNA_per_cell_ug'] + raw_sim_table['DNA_per_cell_ug']
            raw_sim_table['mRNA_mass'] = read_stacked_columns(
                cells, 'Mass', 'mRnaMass', fun=lambda x: x.mean(),
                remove_first=True) * nucleic_acid_mass_units.asNumber(mass_units)
            raw_sim_table['rRNA_mass'] = read_stacked_columns(
                cells, 'Mass', 'rRnaMass', fun=lambda x: x.mean(),
                remove_first=True) * nucleic_acid_mass_units.asNumber(mass_units)
            raw_sim_table['tRNA_mass'] = read_stacked_columns(
                cells, 'Mass', 'tRnaMass', fun=lambda x: x.mean(),
                remove_first=True) * nucleic_acid_mass_units.asNumber(mass_units)

            raw_sim_table['DNA_per_cell_geq'] = raw_sim_table['DNA_per_cell_ug'] / geq_mass
            raw_sim_table['DNA_per_mass'] = raw_sim_table['DNA_per_cell_geq'] / raw_sim_table['mass_per_cell']
            raw_sim_table['protein_per_cell_aa'] = raw_sim_table['protein_per_cell_ug'] / avg_aa_mass
            raw_sim_table['RNA_per_cell_nt'] = raw_sim_table['RNA_per_cell_ug'] / avg_rna_mass
            raw_sim_table['protein_per_mass'] = raw_sim_table['protein_per_cell_aa'] / raw_sim_table['mass_per_cell']
            raw_sim_table['RNA_per_mass'] = raw_sim_table['RNA_per_cell_nt'] / raw_sim_table['mass_per_cell']
            raw_sim_table['DNA_per_mass'] = raw_sim_table['DNA_per_cell_geq'] / raw_sim_table['mass_per_cell']
            raw_sim_table['protein_per_genome'] = raw_sim_table['protein_per_cell_aa'] / raw_sim_table['DNA_per_cell_geq']
            raw_sim_table['RNA_per_genome'] = raw_sim_table['RNA_per_cell_nt'] / raw_sim_table['DNA_per_cell_geq']
            raw_sim_table['origins_per_cell'] = read_stacked_columns(cells, 'ReplicationData', 'numberOfOric', fun=lambda x: x.mean(), remove_first=True)
            raw_sim_table['replication_forks_per_cell'] = read_stacked_columns(
                cells, 'ReplicationData', 'fork_coordinates',
                fun=lambda x: (np.logical_not(np.isnan(x)).sum(axis=1)).mean(), remove_first=True)
            raw_sim_table['termini_per_cell'] = read_stacked_columns(
                cells, 'UniqueMoleculeCounts', 'uniqueMoleculeCounts',
                fun=lambda x: x[:, full_chromosome_index].mean(), remove_first=True)
            raw_sim_table['origins_per_genome'] = raw_sim_table['origins_per_cell'] / raw_sim_table['DNA_per_cell_geq']
            raw_sim_table['protein_per_origin'] = raw_sim_table['protein_per_cell_aa'] / raw_sim_table['origins_per_cell']

            # In Dennis & Bremmer (2021), the stable fraction of total RNA
            # is calculated by estimating the fraction of total RNA that is
            # mRNA, then subtracting from 1. In the model, unstable RNAs
            # also include miscellaneous RNAs not accounted for in the mRNA
            # count. We thus use the more "accurate" measure of tRNA + rRNA,
            # though the difference between the two methods is quite small,
            # at most 0.01.
            raw_sim_table['total_RNA_stable_fraction'] = (raw_sim_table['tRNA_mass'] + raw_sim_table['rRNA_mass']) / raw_sim_table['RNA_per_cell_ug']
            raw_sim_table['stable_RNA_tRNA_fraction'] = raw_sim_table['tRNA_mass'] / (raw_sim_table['tRNA_mass'] + raw_sim_table['rRNA_mass'])
            # We divide by 3 because there are three copies of the rrn operon,
            # but Dennis & Bremmer (2021) calculate the copy numbers of a single
            # operon
            raw_sim_table['rrn_genes_per_cell'] = read_stacked_columns(
                cells, 'RnaSynthProb', 'gene_copy_number',
                fun=lambda x: np.sum(x[:, is_rrn], axis=1).mean(), remove_first=True) / 3
            raw_sim_table['rrn_genes_per_genome'] = raw_sim_table['rrn_genes_per_cell'] / raw_sim_table['DNA_per_cell_geq']
            raw_sim_table['rrn_gene_initiation_rate'] = read_stacked_columns(
                cells, 'RnapData', 'rnaInitEvent',
                fun=lambda x: np.sum(x[:, is_rRNA]), remove_first=True) / (
                raw_sim_table['rrn_genes_per_cell'] * raw_sim_table['doubling_time_sec'])
            raw_sim_table['active_RNAP_per_cell'] = read_stacked_columns(
                cells, 'UniqueMoleculeCounts', 'uniqueMoleculeCounts',
                fun=lambda x: x[:, active_rnap_index].mean(), remove_first=True)
            raw_sim_table['RNAP_per_cell'] = read_stacked_columns(
                cells, 'BulkMolecules', 'counts',
                fun=lambda x: x[:, rnap_index].mean(), remove_first=True) + raw_sim_table['active_RNAP_per_cell']
            raw_sim_table['RNAP_active_fraction'] = raw_sim_table['active_RNAP_per_cell'] / raw_sim_table['RNAP_per_cell']
            # Dennis & Bremmer (2021) calculate via amino acid ratios, but we
            # calculate via mass ratios
            raw_sim_table['RNAP_per_total_protein'] = raw_sim_table['RNAP_per_cell'] * rnap_mass / raw_sim_table['protein_per_cell_ug']
            raw_sim_table['active_ribosomes_per_cell'] = read_stacked_columns(
                cells, 'UniqueMoleculeCounts', 'uniqueMoleculeCounts',
                fun=lambda x: x[:, active_ribosome_index].mean(), remove_first=True)
            # Bremmer & Dennis (2021) use the total ribosomal rna mass to
            # calculate total ribosomes, so we use the mass average amount of
            # subunits instead of taking the min across subunits
            raw_sim_table['inactive_ribosomes_per_cell'] = read_stacked_columns(
                cells, 'BulkMolecules', 'counts',
                fun=lambda x: np.mean(x[:, inactive_ribosome_mask],
                axis=0) @ ribosome_subunit_mass_fractions, remove_first=True)
            raw_sim_table['ribosomes_per_cell'] = raw_sim_table['active_ribosomes_per_cell'] + raw_sim_table['inactive_ribosomes_per_cell']
            raw_sim_table['r_prot_aa_counts'] = read_stacked_columns(
                cells, 'MonomerCounts', 'monomerCounts',
                fun=lambda x: np.mean(x[:, is_r_prot], axis=0) @ r_prot_aa_counts, remove_first=True)
            raw_sim_table['r_prot_per_total_protein'] = raw_sim_table['r_prot_aa_counts'] / raw_sim_table['protein_per_cell_aa']
            raw_sim_table['RNAP_per_ribosome'] = raw_sim_table['RNAP_per_cell'] / raw_sim_table['ribosomes_per_cell']
            raw_sim_table['tRNA_per_cell'] = read_stacked_columns(
                cells, 'BulkMolecules', 'counts',
                fun=lambda x: np.mean(np.sum(x[:, bulk_tRNA_mask], axis=1)),
                remove_first=True)
            raw_sim_table['peptide_elongation_rate'] = read_stacked_columns(
                cells, 'RibosomeData', 'effectiveElongationRate',
                fun=lambda x: x.mean(), remove_first=True)

            raw_sim_table['active_RNAP_synthesizing_stable_fraction'] = read_stacked_columns(
                cells, 'RnapData', 'active_rnap_on_stable_RNA_indexes',
                fun=lambda x: np.mean(np.array([np.count_nonzero(~np.isnan(indices)) for indices in x]))) / raw_sim_table['active_RNAP_per_cell']
            raw_sim_table['stable_RNA_synthesis_per_cell'] = read_stacked_columns(
                cells, "TranscriptElongationListener", "countRnaSynthesized",
                fun=lambda x: np.sum(x[:, is_stable], axis=0) @ stable_rna_lengths, remove_first=True) / raw_sim_table['doubling_time_sec']
            raw_sim_table['mRNA_synthesis_per_cell'] = read_stacked_columns(
                cells, "TranscriptElongationListener", "countRnaSynthesized",
                fun=lambda x: np.sum(x[:, is_mRNA], axis=0) @ mRNA_lengths, remove_first=True) / raw_sim_table['doubling_time_sec']
            raw_sim_table['total_RNA_synth'] = read_stacked_columns(
                cells, "RnapData", "actualElongations", fun=lambda x: x.sum(),
                remove_first=True) / raw_sim_table['doubling_time_sec']
            raw_sim_table['RNA_synth_stable_fraction'] = raw_sim_table['stable_RNA_synthesis_per_cell'] / raw_sim_table['total_RNA_synth']
            raw_sim_table['ppGpp_counts'] = read_stacked_columns(
                cells, 'BulkMolecules', 'counts',
                fun=lambda x: np.mean(x[:, ppGpp_index]), remove_first=True)
            raw_sim_table['ppGpp_concn_per_mass'] = raw_sim_table['ppGpp_counts'] * counts_to_moles / raw_sim_table['mass_per_cell']
            raw_sim_table['ppGpp_concn_per_protein'] = raw_sim_table['ppGpp_counts'] * counts_to_moles / (raw_sim_table['protein_per_cell_ug'] / avg_aa_mass)

            # Calculate mean and standard deviation for each value across
            # the cells of the current variant
            for var_name in plotted_variables['name']:
                sim_table[var_name] = np.append(sim_table[var_name], raw_sim_table[var_name])
                sim_mean_table[var_name][varIdx] = np.mean(raw_sim_table[var_name])
                sim_std_table[var_name][varIdx] = raw_sim_table[var_name].std()

        # Make plots
        fig = plt.figure(figsize=FIG_SIZE)
        for (var_name, var_units, title, position, y_lim) in plotted_variables:
            if var_name != 'doubling_time':
                self.plot_lines(fig, position, sim_table['doubling_time'], sim_table[var_name], sim_mean_table['doubling_time'], sim_mean_table[var_name],
                                sim_std_table[var_name], val_table['doubling_time'], val_table[var_name], title, y_lim)
        plt.tight_layout()

        exportFigure(plt, plotOutDir, plotOutFileName, metadata)

if __name__ == "__main__":
    Plot().cli()
