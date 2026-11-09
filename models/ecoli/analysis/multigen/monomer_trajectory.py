"""
# TODO (mia): for the mRNA count plot, add each TU as a separate line, and
have each's avg count + STD in the legend with it. have its basal synth prob
 and avg. synth prob in either the legend or the title too. for the mRNA itself,
include the avg. full, partial, and total count breakdown. you could do the
same on the plot for each TU too since we track that, but might be too much
info on one. could make each TU its own color and have each state a different
linestyle? but total might overalapp all the time with one of those two idk.


Multigen analysis plot for tracking monomer counts and their trajectories over time.

Multigeneration analysis for analyzing monomer counts and their trajectories
over time. The output contains four panels to help track the dynamics
of monomer usage and complex formation for each monomer of interest:

1. Monomer counts over time, decomposed by state: free, in each parent complex
   (complexation, equilibrium, two component system (TCS) complexes), and in active
   ribosomes, RNAPs, replisomes, and DNA-bound transcription factors (TFs).
2. Events affecting the free monomer pool each timestep (+ adds, - removes):
   synthesis (MonomerCounts.monomersElongated) and degradation
   (MonomerCounts.monomersDegraded), complexation/equilibrium/TCS reaction
   events (TCS estimated via NNLS from the protein rows of
   TwoComponentSystems.moleculeChanges, with the number of timesteps the fit
   misses, if any, shown in the legend), ribosome/RNAP/replisome
   assembly-disassembly when the monomer enters those directly, and
   TF binding/unbinding (RnaSynthProb.nActualBound/nActualUnbound).
   Each legend entry gives the line's average |events| per timestep
   ("avg. X/step"), leaving out the first steps of each generation (which may
   be clipped). Equilibrium lines are net events per step, so they go both
   ways.
3. Its mRNA cistron counts, with the cistron's mRNA half-life and the fitted
   half-lives of the transcription units containing it.
4. The monomer's proteome mass fraction.

Note about replisomes and mechanistic_replisome:
The active_replisome count from UniqueMoleculeCounts is the number of
replisomes (one per replication fork) actually doing replication. But when
mechanistic_replisome is False (the default), initiation happens when cell mass
per origin reaches a critical mass (this is a fixed rule that is not actually
DnaA-regulated currently in the model). The simulation never checks whether
enough replisome subunits actually exist (it assumes they are not limiting),
and never removes any from the free pool (those replisomes carry no protein
mass either). Only with mechanistic_replisome True does each initiation remove
3 of each trimer subunit and 1 of each monomer subunit per replisome from the
bulk pool and physically into the unique molecule space (they are then returned
at termination). See chromosome_replication.py for more information on all this.

The MonomerCounts listener only counts subunits as in active replisomes when
mechanistic_replisome is True. This plot reads the option from the simulation
metadata and when it is off:
  - plots active_replisome * stoich as "# of monomers chromosome replication
    assumes are available to make active replisomes" (dotted) without adding
    it to the complexed or total counts,
  - omits replisome assembly from the free-pool events panel.
With the option on, replisome assembly and disassembly are both plotted.
"""
import pickle
import os
from matplotlib import pyplot as plt
import numpy as np
from models.ecoli.analysis import multigenAnalysisPlot
from models.ecoli.analysis.complex_membership_helpers import (
    build_subunit_categories, build_tcs_reaction_event_estimates,
    check_validity_and_get_compartment, ComplexNetwork,
    format_parent_path_sentence, get_subunit_category_names,
    get_tcs_small_molecules, is_a_replisome_subunit, is_a_ribosome_subunit,
    is_a_TF, is_an_RNAP_subunit, replisome_subunit_stoich, tf_index)
from wholecell.analysis.analysis_tools import (count_replisomes_formed_and_removed,
    exportFigure, extract_doubling_times, get_first_daughter_lineage,
    read_stacked_bulk_molecules, read_stacked_columns)
from wholecell.analysis.plotting_tools import format_event_lines
from wholecell.io.tablereader import TableReader
from wholecell.utils import units

# USER INPUTS:

# List of proteins to be plotted (with or without compartment tags):
PLOT_MONOMERS = ["PD00423[c]", "EG10504-MONOMER"]

# Specify the # of timesteps at the start of each generation where lines can be
# truncated in the events panel to avoid having large spikes throw off the
# y axis. The first generation starts with no equilibrium or TCS complexes
# (initialization only forms complexation complexes), so those reactions form
# a significant number of molecules in its first steps after initialization to
# hit steady state. Later generations inherit a count level closer to steady
# state, so their first steps only spike a few times above the usual, but can
# still be large relative to the average activity of these reactions (see
# the truncating logic in the code for more details):
N_INITIAL_STEPS_TO_CLIP = 5

# Colors for the dotted parent-complex lines in panel 1 (tab10 without its red
# and purple, which are used for the total complexed and replisome lines):
PARENT_COMPLEX_COLORS = ['tab:blue', 'tab:orange', 'tab:green', 'tab:brown',
                         'tab:pink', 'tab:olive', 'tab:cyan', 'tab:gray']

# Colors for the reaction event lines in panel 2 (note: synthesis and degradation
# have their own fixed colors, lightseagreen and red, and are left out here):
EVENT_REACTION_COLORS = ['cornflowerblue', 'tab:orange', 'tab:purple', 'tab:olive',
                         'tab:pink', 'tab:brown', 'tab:gray', 'tab:green']
# END USER INPUTS

def classify_monomer(monomer_id, subunit_categories):
    """
    Describe all the groups in build_subunit_categories() a monomer belongs to,
    with its transcription factor type(s) when it is part of a TF.
    """
    categories = []
    for name in get_subunit_category_names(monomer_id, subunit_categories):
        if name == 'Transcription factor':
            tf_types = '/'.join(sorted(subunit_categories['tf_monomer_types'][monomer_id]))
            name = f'{tf_types} transcription factor'
        elif name == 'Two component system':
            name = 'two component system'
        categories.append(name)

    if not categories:
        return 'Monomer only'
    return ' + '.join(sorted(categories)) + ' subunit'


class Plot(multigenAnalysisPlot.MultigenAnalysisPlot):
    def check_monomer_validity(self, sim_data, molecule_list):
        """
        Verify that molecules are valid monomers in the simulation.

        Args:
            molecule_list: List of molecule IDs to check
        Returns:
            list: valid monomer IDs
        """
        # First, check validity of molecule IDs and get compartment tags:
        revised_molecule_list = check_validity_and_get_compartment(
                                    sim_data, molecule_list)

        # Next, check that the molecules are valid monomers:
        all_monomers = sim_data.process.translation.monomer_data['id'].tolist()
        valid_monomers = []
        for molecule in revised_molecule_list:
            if molecule in all_monomers:
                valid_monomers.append(molecule)
            else:
                print(f"{molecule} is not a valid monomer.")

        return valid_monomers

    def do_plot(self, seedOutDir, plotOutDir, plotOutFileName, simDataFile,
                validationDataFile, metadata):
        with open(simDataFile, 'rb') as f:
            sim_data = pickle.load(f)
        # Only the first-daughter lineage is plotted:
        cell_paths = get_first_daughter_lineage(self.ap.get_cells())

        # Check if replication actually sequestered replisome subunits (see module
        # docstring). If this is missing from the metadata, it means the
        # default was used (which is typically set to False as of 10/21/2026):
        mechanistic_replisome = metadata.get('mechanistic_replisome', False)

        # Monomer-level reconstruction properties:
        self.monomerIDs = sim_data.process.translation.monomer_data['id']
        monomer_idx_dict = {m: i for i, m in enumerate(self.monomerIDs)}
        deg_rate = sim_data.process.translation.monomer_data['deg_rate'].asNumber(
            1 / units.s)
        with np.errstate(divide='ignore'):
            half_lives = np.log(2) / deg_rate / 60  # minutes
        cistron_ids = sim_data.process.translation.monomer_data['cistron_id']
        common_names = [sim_data.common_names.get_common_name(c) for c in cistron_ids]
        mw = sim_data.process.translation.monomer_data['mw']

        # Extract each cistron's HL (used by ParCa to fit the HL of each TU
        # containing the cistron) and each TU's HL:
        transcription = sim_data.process.transcription
        cistron_row = {c: i for i, c in enumerate(transcription.cistron_data['id'])}
        with np.errstate(divide='ignore'):
            cistron_half_lives = (np.log(2)
                / transcription.cistron_data['deg_rate'].asNumber(1 / units.s) / 60)
            tu_half_lives = (np.log(2)
                / transcription.rna_data['deg_rate'].asNumber(1 / units.s) / 60)
        tu_ids = transcription.rna_data['id']

        def mrna_half_life_str(cistron_id):
            """
            Get a string describing the mRNA half-life (and its respective TU
            HLs) for a given cistron ID.
            """
            if cistron_id not in cistron_row:
                return "mRNA half-life: n/a"
            i = cistron_row[cistron_id]
            tus = transcription.cistron_tu_mapping_matrix.getrow(i).indices
            tu_str = ", ".join(
                f"{tu_ids[j].split('[')[0]} {tu_half_lives[j]:.1f}" for j in tus)
            return (f"mRNA half-life: cistron {cistron_half_lives[i]:.1f} min | "
                    f"TU(s): {tu_str} min")

        complexation_rxn_ids = list(sim_data.process.complexation.ids_reactions)
        equilibrium_rxn_ids = list(sim_data.process.equilibrium.rxn_ids)
        tcs_rxn_ids = list(sim_data.process.two_component_system.rxn_ids)
        tcs_stoich = sim_data.process.two_component_system.stoich_matrix()
        tcs_molecule_names = list(sim_data.process.two_component_system.molecule_names)
        tcs_small_molecules = get_tcs_small_molecules(sim_data)

        # Classify all proteins by the groups they are subunits of:
        subunit_categories = build_subunit_categories(sim_data, monomer_idx_dict)
        self.monomer_classifications = {
            protein_id: classify_monomer(protein_id, subunit_categories)
            for protein_id in self.monomerIDs}

        # Build the molecule to parent complex dictionaries of the three types
        # of complexes:
        self.complex_network = ComplexNetwork(sim_data)

        # Check that all the valid inputs are indeed valid monomers:
        valid_monomers = self.check_monomer_validity(sim_data, PLOT_MONOMERS)
        if not valid_monomers:
            print(f"No valid monomers in {PLOT_MONOMERS}; skipping.")
            return

        # Extract relevant monomer tracking information from listeners:
        free_monomer_counts = read_stacked_columns(
            cell_paths, 'MonomerCounts', "freeMonomerCounts")
        total_monomer_counts = read_stacked_columns(
            cell_paths, 'MonomerCounts', "monomerCounts")
        monomers_degraded = read_stacked_columns(
            cell_paths, 'MonomerCounts', "monomersDegraded")
        monomers_synthesized = read_stacked_columns(
            cell_paths, 'MonomerCounts', "monomersElongated")
        cistron_counts = read_stacked_columns(
            cell_paths, 'RNACounts', 'mRNA_cistron_counts')
        mRNA_cistron_ids = transcription.cistron_data['id'][
            transcription.cistron_data['is_mRNA']]
        total_protein_mass = np.squeeze(
            read_stacked_columns(cell_paths, "Mass", "proteinMass"))
        tf_bound = read_stacked_columns(cell_paths, 'RnaSynthProb', "nActualBound")
        tf_unbound = read_stacked_columns(cell_paths, 'RnaSynthProb', "nActualUnbound")

        # Active ribosome, RNAP and replisome counts:
        unique_molecule_ids = TableReader(os.path.join(
            cell_paths[0], 'simOut', 'UniqueMoleculeCounts')).readAttribute(
            'uniqueMoleculeIds')
        unique_molecule_counts = read_stacked_columns(
            cell_paths, 'UniqueMoleculeCounts', 'uniqueMoleculeCounts')
        active_ribosome = unique_molecule_counts[
            :, unique_molecule_ids.index('active_ribosome')]
        active_rnap = unique_molecule_counts[
            :, unique_molecule_ids.index('active_RNAP')]
        active_replisome = unique_molecule_counts[
            :, unique_molecule_ids.index('active_replisome')]

        # Extract time and generation boundaries:
        (time, _, end_generation_times, generation_start_idxs, _) = (
            extract_doubling_times(cell_paths))

        # Complexation/equilibrium events are read only if a monomer of interest
        # needs them, and TCS reaction events are NNLS-estimated once since
        # they are not recorded during the simulation:
        listener_cache = {}

        def complexation_events():
            if 'complexation' not in listener_cache:
                listener_cache['complexation'] = read_stacked_columns(
                    cell_paths, 'ComplexationListener', "complexationEvents")
            return listener_cache['complexation']

        def equilibrium_events():
            if 'equilibrium' not in listener_cache:
                listener_cache['equilibrium'] = read_stacked_columns(
                    cell_paths, 'EquilibriumListener', "complexationEvents")
            return listener_cache['equilibrium']

        def tcs_events():
            if 'tcs' not in listener_cache:
                listener_cache['tcs'] = build_tcs_reaction_event_estimates(
                    sim_data, read_stacked_columns(
                        cell_paths, 'TwoComponentSystems', 'moleculeChanges'))
            return listener_cache['tcs']

        def tcs_label(rxn_idx, label):
            # Misfits on this reaction's own proteins (shared ATP/ADP/Pi/H2O/H+
            # rows are not fit; see build_tcs_reaction_event_estimates):
            _, tcs_misfit = tcs_events()
            rows = [r for r in np.nonzero(tcs_stoich[:, rxn_idx])[0]
                    if tcs_molecule_names[r] not in tcs_small_molecules]
            n_misfit = int(tcs_misfit[:, rows].any(axis=1).sum())
            if n_misfit:
                label += f"\n[fit misses at {n_misfit}/{len(time)} steps]"
            return label

        # Generate a separate trajectory plot for each monomer:
        for monomer_id in valid_monomers:
            monomer_idx = monomer_idx_dict[monomer_id]
            monomer_classification = self.monomer_classifications[monomer_id]
            half_life = half_lives[monomer_idx]
            common_name = common_names[monomer_idx]
            molecular_weight = mw[monomer_idx]

            print(f"\n{'=' * 60}\nAnalyzing monomer: {monomer_id}\n{'=' * 60}")
            print(f"Classification: {monomer_classification} | Common name: {common_name}")
            print(f"Half-life: {half_life:.2f} min, {half_life / 60:.2f} hrs")

            # Trace downstream complexes the monomer contributes to, and the
            # unique paths to each complex:
            downstream_complexes, _ = self.complex_network.analyze_molecule_parents(
                monomer_id, molecule_label="Monomer")
            complexes_by_id = {
                complex_id: [
                    {'path_data': path_data,
                     'path_string': format_parent_path_sentence(
                         "Monomer", complex_id, path_data)}
                    for path_data in path_list]
                for complex_id, path_list in downstream_complexes.items()}
            unique_complex_path_pairs = [
                (complex_id, p)
                for complex_id, paths in complexes_by_id.items()
                for p in paths]
            print(f"Unique complex-path combinations: {len(unique_complex_path_pairs)}")

            # Determine the free counts of each parent complex:
            complex_free_counts = {}
            for complex_id in complexes_by_id:
                complex_free_counts[complex_id] = (
                    read_stacked_bulk_molecules(cell_paths, [complex_id]))[0]

            # Extract free and total counts for the monomer:
            monomer_total_counts = total_monomer_counts[:, monomer_idx]
            monomer_free_counts = free_monomer_counts[:, monomer_idx]

            # Generate the plot:
            fig, (ax1, ax2, ax3, ax4) = plt.subplots(
                4, 1, sharex=True, figsize=(10, 12))

            # PLOT 1: Monomer counts over time
            ax1.plot(time, monomer_total_counts, color='black',
                     label='Total monomer counts', linewidth=.5, alpha=1)
            ax1.plot(time, monomer_free_counts, color='lightseagreen',
                     label='Free monomer counts', linewidth=.5, alpha=0.9)
            # Total complexed: created here so it is listed third and drawn under
            # the individual complex lines, but its values are filled in once every
            # contribution has been summed below:
            (total_complexed_line,) = ax1.plot(
                time, np.zeros(len(time)), color='red',
                label='Total complexed monomer counts\n(sum of all complex contributions)',
                linewidth=0.6, alpha=0.6, linestyle='--')

            # Track complexed counts to verify total = free + complexed, and
            # which unique molecule types were already counted via complexes
            # (to avoid double counting):
            total_complexed_counts = np.zeros_like(monomer_total_counts, dtype=float)
            already_counted_in_replisomes = False
            already_counted_as_tf = False

            if unique_complex_path_pairs:
                for i, (complex_id, path_info) in enumerate(unique_complex_path_pairs):
                    total_stoich = abs(path_info['path_data']['total_stoich'])
                    contrib = complex_free_counts[complex_id] * total_stoich
                    total_complexed_counts += contrib

                    n_paths = len(complexes_by_id[complex_id])
                    label = path_info['path_string']
                    if n_paths > 1:
                        path_num = complexes_by_id[complex_id].index(path_info) + 1
                        label += f"\n[Path {path_num} of {n_paths}]"
                    ax1.plot(time, contrib,
                             color=PARENT_COMPLEX_COLORS[i % len(PARENT_COMPLEX_COLORS)],
                             alpha=0.9, label=label, linewidth=0.75, linestyle=":")

                    # Complex that is itself a DNA-bound TF:
                    if is_a_TF(sim_data, complex_id) != []:
                        tf_idx = tf_index(sim_data, complex_id)
                        contrib = tf_bound[:, tf_idx] * total_stoich
                        total_complexed_counts += contrib
                        ax1.plot(time, contrib, color='lightcoral',
                                 label=f'Monomer in bound {complex_id} TFs\n'
                                       f'({total_stoich:g} per)',
                                 linewidth=0.75, alpha=0.75, linestyle="-.")
                        already_counted_as_tf = True

                    if is_a_ribosome_subunit(sim_data, complex_id) != []:
                        contrib = active_ribosome * total_stoich
                        total_complexed_counts += contrib
                        ax1.plot(time, contrib, color='gold',
                                 label=f'Monomer in ribosomes via {complex_id}\n'
                                       f'({total_stoich:g} per)',
                                 linewidth=0.75, alpha=0.75, linestyle="-.")

                    if is_an_RNAP_subunit(sim_data, complex_id) != []:
                        contrib = active_rnap * total_stoich
                        total_complexed_counts += contrib
                        ax1.plot(time, contrib, color='orange',
                                 label=f'Monomer in active RNAPs via {complex_id}\n'
                                       f'({total_stoich:g} per)',
                                 linewidth=0.75, alpha=0.75, linestyle="-.")

                    if is_a_replisome_subunit(sim_data, complex_id) != []:
                        rep_stoich = replisome_subunit_stoich(sim_data, complex_id)
                        contrib = active_replisome * total_stoich * rep_stoich
                        if mechanistic_replisome:
                            total_complexed_counts += contrib
                            label = (f'Monomer in active replisomes via {complex_id}\n'
                                     f'({total_stoich * rep_stoich:g} per)')
                        else:
                            # Not real protein: not added to the totals
                            label = ('# of monomers chromosome replication\n'
                                     'assumes are available to make active\n'
                                     'replisomes, but mechanistic_replisome is\n'
                                     'off so they are never actually removed\n'
                                     f'from bulk ({total_stoich * rep_stoich:g} per, '
                                     f'via {complex_id})')
                        ax1.plot(time, contrib, color='purple', label=label,
                                 linewidth=0.75, alpha=0.9,
                                 linestyle="-." if mechanistic_replisome else ":")
                        already_counted_in_replisomes = True

            # Monomers directly in active replisomes or DNA-bound TFs (only if
            # not via a complex). Ribosome and RNAP subunits are not checked
            # here as monomers only enter active ribosomes via the 30S or 50S
            # subunit and active RNAPs via the inactive RNAP complex
            # (APORNAP-CPLX), which are already counted in the complex loop
            # above:
            if (is_a_replisome_subunit(sim_data, monomer_id) != []
                    and not already_counted_in_replisomes):
                rep_stoich = replisome_subunit_stoich(sim_data, monomer_id)
                contrib = active_replisome * rep_stoich
                if mechanistic_replisome:
                    total_complexed_counts += contrib
                    label = f'Monomer directly in active replisomes ({rep_stoich:g} per)'
                else:
                    # if mechanistic replisome is not on, then the replisomes
                    # recorded do not actually pull monomers from the free pool
                    # since they are assumed to be in excess, so they are not
                    # counted in the total:
                    label = ('# of monomers chromosome replication\n'
                             'assumes are available to make active\n'
                             'replisomes, but mechanistic_replisome is\n'
                             'off so they are never actually removed\n'
                             f'from bulk ({rep_stoich:g} per)')
                ax1.plot(time, contrib, color='purple', label=label,
                         linewidth=0.75, alpha=0.9,
                         linestyle="-." if mechanistic_replisome else ":")
            if is_a_TF(sim_data, monomer_id) != [] and not already_counted_as_tf:
                contrib = tf_bound[:, tf_index(sim_data, monomer_id)]
                total_complexed_counts += contrib
                ax1.plot(time, contrib, color='lightcoral',
                         label='Monomer directly in bound TFs (1 per)',
                         linewidth=0.75, alpha=0.75, linestyle="-.")

            total_complexed_line.set_ydata(total_complexed_counts)

            # Verify mass balance: total should equal free + complexed:
            calculated_total = monomer_free_counts + total_complexed_counts
            discrepancy = monomer_total_counts - calculated_total
            max_discrepancy = float(np.max(np.abs(discrepancy))) if len(time) else 0.
            print(f"Max mass-balance discrepancy = {max_discrepancy:.1f}")
            if max_discrepancy > 1:
                # Drawn on top and dotted so it shows even where it overlaps
                # another line (i.e. the free counts when nothing is complexed):
                ax1.plot(time, calculated_total, color='deeppink',
                         label=f'Calculated total (free + complexed)\n'
                               f'Max discrepancy: {max_discrepancy:.1f}',
                         linewidth=0.75, alpha=0.9, linestyle=':', zorder=5)

            stats_str = (
                f"Mean counts - Total: {np.mean(monomer_total_counts):.1f} | "
                f"Free: {np.mean(monomer_free_counts):.1f} | "
                f"Complexed: {np.mean(total_complexed_counts):.1f} | "
                f"Median counts - Total: {np.median(monomer_total_counts):.1f} | "
                f"Free: {np.median(monomer_free_counts):.1f} | "
                f"Complexed: {np.median(total_complexed_counts):.1f}")
            ax1.legend(fontsize=6, loc="center left", bbox_to_anchor=(1, 0.5))
            ax1.set_ylabel("Monomer Counts", fontsize=10)
            ax1.set_title(stats_str, fontsize=7)

            fig.suptitle(
                f"{monomer_classification} | "
                f"Half-life: {half_life:.1f} min, {half_life / 60:.2f} hrs\n"
                f"Sim ID: {metadata['description']} | "
                f"variant {metadata.get('variant', '?')} {metadata.get('variant_index', '')} | "
                f"seed {int(metadata.get('seed', 0))} | {len(cell_paths)} gen(s)",
                fontsize=9, x=.4, y=0.976, color='grey')
            # Bold first line, drawn separately (one title can't mix weights)
            # and placed just above the suptitle (inside the figure, since
            # exportFigure does not expand the canvas to fit outside text):
            fig.text(0.4, 0.998, f"Monomer count trajectory of {monomer_id} ({common_name})",
                     fontsize=10, fontweight='bold', ha='center', va='top')

            # PLOT 2: Events affecting the free monomer pool

            # Reaction event lines in panel 2 colors:
            ax2.set_prop_cycle(color=EVENT_REACTION_COLORS)

            ax2.plot(time, monomers_synthesized[:, monomer_idx], color='lightseagreen',
                     label='Synthesis events',
                     linewidth=0.75, alpha=0.75)
            ax2.plot(time, -monomers_degraded[:, monomer_idx], color='red',
                     label='Degradation events',
                     linewidth=0.75, alpha=0.75)

            # Extract the list of reactions that will be plotted:
            rxns_plotted = []
            # Reactions that directly consume the free monomer (first path step):
            for parent_complex, path_info in unique_complex_path_pairs:
                path = path_info['path_data']['path']
                if path[0][0] != monomer_id:
                    continue
                rxn_id = path[0][3]
                if rxn_id in rxns_plotted:
                    continue
                rxns_plotted.append(rxn_id)
                rxn_stoich = abs(path[0][1])

                if rxn_id in complexation_rxn_ids:
                    rxn_events = complexation_events()[
                        :, complexation_rxn_ids.index(rxn_id)]
                    ax2.plot(time, -rxn_events * rxn_stoich, alpha=0.6,
                             label=f'{rxn_id}\n(consumes {rxn_stoich:g} free {monomer_id})',
                             linewidth=0.6)
                elif rxn_id in equilibrium_rxn_ids:
                    rxn_events = equilibrium_events()[
                        :, equilibrium_rxn_ids.index(rxn_id)]
                    # Equilibrium events are net per step, so the same line
                    # consumes (complex formation) or releases (dissociation):
                    ax2.plot(time, -rxn_events * rxn_stoich, alpha=0.6,
                             label=f'{rxn_id} net events\n(consumes/releases {rxn_stoich:g} free '
                                   f'{monomer_id})',
                             linewidth=0.6)
                elif rxn_id in tcs_rxn_ids:
                    rxn_idx = tcs_rxn_ids.index(rxn_id)
                    ax2.plot(time, -tcs_events()[0][:, rxn_idx] * rxn_stoich, alpha=0.6,
                             label=tcs_label(rxn_idx, f'{rxn_id}\n'
                                 f'(estimated, consumes {rxn_stoich:g} free {monomer_id})'),
                             linewidth=0.6)

                    # The same TCS complex can also form via the other
                    # phosphotransfer reaction type:
                    if "POS" in rxn_id:
                        immediate_product = path[1][0] if len(path) > 1 else parent_complex
                        for extra_rxn in self.complex_network.find_extra_tcs_reactions(
                                monomer_id, immediate_product, rxns_plotted):
                            rxns_plotted.append(extra_rxn)
                            extra_idx = tcs_rxn_ids.index(extra_rxn)
                            ax2.plot(time, -tcs_events()[0][:, extra_idx] * rxn_stoich,
                                     alpha=0.6,
                                     label=tcs_label(extra_idx, f'{extra_rxn}\n'
                                         f'(estimated, consumes {rxn_stoich:g} free {monomer_id})'),
                                     linewidth=0.6)

            # Equilibrium reverse reactions where the monomer is a direct product
            # (equilibrium events are net, so dissociation = negative events):
            if monomer_id in self.complex_network.eq_m2pc_dict:
                for eq_complex_id, info in self.complex_network.eq_m2pc_dict[monomer_id].items():
                    rxn_id = info[0]['reaction_id']
                    if rxn_id in rxns_plotted:
                        continue
                    rxns_plotted.append(rxn_id)
                    release_stoich = abs(info[1]['stoichiometry'])
                    rxn_events = equilibrium_events()[
                        :, equilibrium_rxn_ids.index(rxn_id)]
                    dissociation_events = np.minimum(rxn_events, 0) * -1
                    if np.any(dissociation_events > 0):
                        ax2.plot(time, dissociation_events * release_stoich, alpha=0.6,
                                 label=f'{rxn_id} (reverse)\n'
                                       f'(releases {release_stoich:g} free {monomer_id})',
                                 linewidth=0.6)

            # TCS reactions that directly produce the free monomer (positive
            # stoichiometry) are dephosphorylations for response regulators,
            # phosphotransfers for histidine kinases:
            if monomer_id in tcs_molecule_names:
                monomer_row = tcs_stoich[tcs_molecule_names.index(monomer_id)]
                for rxn_idx in np.where(monomer_row > 0)[0]:
                    rxn_id = tcs_rxn_ids[rxn_idx]
                    if rxn_id in rxns_plotted:
                        continue
                    rxns_plotted.append(rxn_id)
                    produced = monomer_row[rxn_idx]
                    ax2.plot(time, tcs_events()[0][:, rxn_idx] * produced, alpha=0.6,
                             label=tcs_label(rxn_idx, f'{rxn_id}\n'
                                 f'(estimated, produces {produced:g} free {monomer_id})'),
                             linewidth=0.6)

            # No need to check for direct ribosome membership as monomers only
            # become active ribosomes once they have already become the 30S or
            # 50S subunit, so ribosome assembly would never affect the free
            # monomer pool. Same goes for RNAPs, as monomers only become active
            # RNAPs once they have already become the inactive
            # RNAP complex (APORNAP-CPLX).

            # Replisomes: only use if the monomer enters replisomes directly, and only
            # when mechanistic_replisome is on (otherwise replisome formation never
            # removes subunits from the free pool, see module docstring):
            if (is_a_replisome_subunit(sim_data, monomer_id) != []
                    and mechanistic_replisome and not any(
                        is_a_replisome_subunit(sim_data, c) != []
                        for c, _ in unique_complex_path_pairs)):
                rep_stoich = replisome_subunit_stoich(sim_data, monomer_id)
                # Replisomes formed in a step consume subunits (assembly) and
                # replisomes removed return them to the free pool (termination),
                # counted separately (see count_replisomes_formed_and_removed).
                # Changes across division are not events, so those rows are 0:
                replisomes_formed, replisomes_removed = (
                    np.concatenate(counts) for counts in zip(*(
                        count_replisomes_formed_and_removed(cell_path)
                        for cell_path in cell_paths)))
                ax2.plot(time, replisomes_removed * rep_stoich,
                         color='purple',
                         label=f'Replisome disassembly\n'
                               f'(releases {rep_stoich:g} free {monomer_id})',
                         linewidth=0.75, alpha=0.75)
                ax2.plot(time, -replisomes_formed * rep_stoich,
                         color='purple',
                         label=f'Replisome assembly\n'
                               f'(consumes {rep_stoich:g} free {monomer_id})',
                         linewidth=0.75, alpha=0.75)

            # TFs: binding / unbinding only if the free monomer itself is the TF
            # (not via a TF complex):
            if is_a_TF(sim_data, monomer_id) != [] and not any(
                    is_a_TF(sim_data, c) != [] for c, _ in unique_complex_path_pairs):
                tf_idx = tf_index(sim_data, monomer_id)
                # nActualBound is TF occupancy at the end of each step, not a
                # count of binding events. That works as "binding events" here
                # because TfUnbinding releases every bound TF at the start of
                # each step, so every TF bound at the end of a step was bound
                # during it. Two things to note about TF binding:
                # 1. The first row of each generation is the initial state,
                #    recorded before any process runs, so occupancy there is
                #    inherited from the parent cell (with no matching unbinding
                #    row), not new binding. Plotting it would show a fake drop
                #    of the TF's whole bound amount at every generation start,
                #    so those rows are set to 0.
                # 2. When a replication fork passes a promoter,
                #    ChromosomeStructure returns its bound TFs to the free pool
                #    after TfBinding ran, so occupancy = bindings - those
                #    released TFs. The binding trace therefore reads low by a
                #    few TFs at those (rare) steps, but the panel total still
                #    matches the real free-pool change.
                tf_binding = tf_bound[:, tf_idx].astype(float)
                tf_binding[generation_start_idxs] = 0
                ax2.plot(time, -tf_binding, color='lightcoral',
                         label=f'TF binding events\n(consumes free {monomer_id};\n'
                               'net of TFs released when replication\n'
                               'passes a promoter)',
                         linewidth=0.75, alpha=0.75)
                ax2.plot(time, tf_unbound[:, tf_idx], color='pink',
                         label=f'TF unbinding events\n(releases free {monomer_id})',
                         linewidth=0.75, alpha=0.75)

            # Order and label the events by activity, and clip any large spikes
            # in the first steps of each generation so they do not dwarf the
            # rest of the panel. The first generation starts with no equilibrium
            # or TCS complexes (initialization only forms complexation
            # complexes), so those reactions form nearly all of them at once in
            # its first steps. In later generations, division splits each
            # molecule between the daughters independently (binomially), so
            # complexes and their free subunits start slightly off steady state
            # and the reactions correct this over the first few steps:
            ordered_event_lines = format_event_lines(
                ax2, time, generation_start_idxs, N_INITIAL_STEPS_TO_CLIP)

            ax2.axhline(y=0, color='black', linestyle='-', linewidth=0.5, alpha=0.5)
            ax2.legend(handles=ordered_event_lines, fontsize=6, loc="center left",
                       bbox_to_anchor=(1, 0.5))
            ax2.set_ylabel("Events\n(+ adds to free pool,\n- removes from free pool)",
                           fontsize=9)

            # PLOT 3: mRNA cistron counts over time
            monomer_cistron_id = cistron_ids[monomer_idx]
            cistron_index = np.where(mRNA_cistron_ids == monomer_cistron_id)[0]
            if len(cistron_index) > 0:
                monomer_cistron_counts = cistron_counts[:, cistron_index[0]]
                ax3.plot(time, monomer_cistron_counts, color='mediumvioletred',
                         label='mRNA counts', linewidth=1, alpha=0.9)
                ax3.set_title(
                    f"mRNA ID: {monomer_cistron_id} | "
                    f"{mrna_half_life_str(monomer_cistron_id)}\n"
                    f"Mean: {np.mean(monomer_cistron_counts):.1f} | "
                    f"Median: {np.median(monomer_cistron_counts):.1f} | "
                    f"Range: [{np.min(monomer_cistron_counts):.0f}, "
                    f"{np.max(monomer_cistron_counts):.0f}]",
                    fontsize=7)
                ax3.legend(fontsize=6, loc="center left", bbox_to_anchor=(1, 0.5))
            else:
                ax3.text(0.5, 0.5,
                         f'Cistron {monomer_cistron_id} not found in mRNA_cistron_counts',
                         ha='center', va='center', transform=ax3.transAxes, fontsize=8)
            ax3.set_ylabel("mRNA Cistron\nCounts", fontsize=9)

            # PLOT 4: Proteome mass fraction over time
            monomer_mass = monomer_total_counts * molecular_weight.asNumber(
                units.fg / units.count)
            with np.errstate(divide='ignore', invalid='ignore'):
                proteome_mass_fraction = np.where(
                    total_protein_mass > 0, monomer_mass / total_protein_mass, np.nan)
            ax4.plot(time, proteome_mass_fraction * 100, color='mediumslateblue',
                     label='Proteome fraction',
                     linewidth=1, alpha=0.9)
            ax4.set_title(
                f"Proteome mass fraction\n"
                f"MW: {molecular_weight.asNumber(units.g / units.mol):,.0f} g/mol | "
                f"Mean: {np.nanmean(proteome_mass_fraction) * 100:.4f}% | "
                f"Median: {np.nanmedian(proteome_mass_fraction) * 100:.4f}% | "
                f"Range: [{np.nanmin(proteome_mass_fraction) * 100:.4f}%, "
                f"{np.nanmax(proteome_mass_fraction) * 100:.4f}%]",
                fontsize=7)
            ax4.legend(fontsize=6, loc="center left", bbox_to_anchor=(1, 0.5))
            ax4.set_ylabel("Proteome Mass\nFraction (%)", fontsize=9)
            ax4.set_xlabel("Time (s)", fontsize=10)

            # Plot generation boundaries on every panel:
            for ax in (ax1, ax2, ax3, ax4):
                for end_time in end_generation_times:
                    ax.axvline(x=end_time, linestyle='--', color="yellowgreen",
                               alpha=0.4, linewidth=0.75)

            plt.tight_layout()
            # Leave room at the top for the two title blocks:
            plt.subplots_adjust(hspace=0.3, top=0.93)

            file_name = f"{plotOutFileName}_{monomer_id}_{common_name}"
            exportFigure(plt, plotOutDir, file_name, metadata, dpi=300)
            plt.close()
            print(f"Plot saved: {file_name}")


if __name__ == '__main__':
    Plot().cli()
