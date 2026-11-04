"""
This plot visualizes the counts of a specified complex over time, along with:
- Free monomer subunit availability
- Complexation events that produce or consume the complex. Each legend entry
  gives the line's average |events| per timestep ("avg. X/step"), leaving out
  the first steps of each generation (which may be truncated so they do not
  throw off the y-axis from the typical average range). Equilibrium lines
  are net events per step, so they can go both ways.
- Counts within parent complexes and unique molecules (ribosomes, RNAPs, etc.)

Supported complex types:
- Complexation complexes (irreversible assembly)
- Equilibrium complexes (reversible ligand binding)
- Two-component system complexes (phosphotransfer/phosporylation/dephosporylation reactions)

In total, there are 1131 unique complexation complexes, 39 unique equilibrium
complexes, and 18 unique two component system complexes.
"""
import html
import pickle
import os
import re
from matplotlib import pyplot as plt
import numpy as np
from models.ecoli.analysis import multigenAnalysisPlot
from models.ecoli.analysis.complex_membership_helpers import (
    build_tcs_reaction_event_estimates, check_validity_and_get_compartment,
    ComplexNetwork, format_parent_path_sentence, is_a_replisome_subunit,
    is_a_ribosome_subunit, is_a_TF, is_an_RNAP_subunit, replisome_subunit_stoich,
    tf_index)
from wholecell.analysis.analysis_tools import (count_replisomes_formed_and_removed,
    exportFigure, extract_doubling_times, get_first_daughter_lineage,
    read_stacked_bulk_molecules, read_stacked_columns)
from wholecell.analysis.plotting_tools import format_event_lines
from wholecell.io.tablereader import TableReader

# USER INPUTS:

# List complexes you would like to plot here (with or without compartment tags):
PLOT_COMPLEXES = ["APORNAP-CPLX"]

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

# END USER INPUTS

def get_clean_common_name(sim_data, molecule_id):
    """
    Get a molecule's common name for use in plot text, without HTML entities
    or tags (i.e. "&beta;", "<i>sn</i>").

    Args:
        sim_data: simulation data
        molecule_id: ID of the molecule (with or without compartment tag)
    Returns:
        str or None: the common name, or None if there is none (some
        equilibrium and TCS complexes only have a note about the reaction that
        makes them in place of a name, e.g. "Reaction required to produce 1CS
        TF ...", so those are left out too)
    """
    base_id = molecule_id.split('[')[0]
    name = sim_data.common_names.get_common_name(base_id)
    name = re.sub(r'<[^>]+>', '', html.unescape(name))
    if name == base_id or name.startswith('Reaction required'):
        return None
    return name

class Plot(multigenAnalysisPlot.MultigenAnalysisPlot):
    def check_complex_validity(self, sim_data, molecule_list):
        """
        Verify that molecules are valid complexes and determine their types.

        Args:
            molecule_list: List of molecule IDs to check
        Returns:
            tuple: (valid_complexes, complex_types)
                - valid_complexes: List of valid complex IDs
                - complex_types: Dict mapping complex IDs to their types
        """
        # First, check validity of molecule IDs and get compartment tags:
        revised_molecule_list = check_validity_and_get_compartment(
                                    sim_data, molecule_list)

        # Next, check that the molecules are valid complexes and determine their types:
        network = self.complex_network
        all_complexes = network.complex_ids + network.eq_complex_ids + network.tcs_complex_ids
        valid_complexes = []
        complex_types = {}

        for molecule in revised_molecule_list:
            if molecule in all_complexes:
                valid_complexes.append(molecule)
                if molecule in network.complex_ids:
                    complex_types[molecule] = "complexation"
                elif molecule in network.eq_complex_ids:
                    complex_types[molecule] = "equilibrium"
                else:
                    complex_types[molecule] = "two component system"
            else:
                print(f"{molecule} is not a valid complex. Only complexation, "
                      f"equilibrium, and two component system complexes can "
                      f"be plotted.")

        return valid_complexes, complex_types

    def determine_complex_makeup(self, complex_id, complex_type):
        """
        Depending on the type of complex (complexation, equilibrium, or two
        component system), determine what combination of monomer subunits the
        complex is made from. If there is only one type of monomer subunit in
        the complex, it is considered homogeneous. If there are multiple types
        of monomer subunits, it is considered heterogeneous.

        NOTE: this classification methodology is only dependent on monomer
        subunit make up and not other molcule types (i.e. if an equilibrium
        complex is made from a ligand binding to a single monomer, the complex
        will still be labeled as homogeneous because it is made from only one
        type of monomer subunit).
        Args:
            complex_id: ID of the complex being analyzed
            complex_type: complex classification: complexation, eq, or TCS.
        Returns:
            str: "homogeneous" or "heterogeneous" for equilibrium and
            complexation complexes, and something from the "pairings" dictionary
            in the build_tcs_dictionaries() function for TCS complexes.

        """
        # Determine the makeup of the complex based on its type:
        if complex_type == 'complexation':
            makeup = self.complex_network.c2ct[complex_id]
        elif complex_type == 'equilibrium':
            makeup = self.complex_network.eq_c2ct[complex_id]
        elif complex_type == 'two component system':
            makeup = self.complex_network.tcs_c2ct[complex_id]
        else:
            makeup = "unknown"

        return makeup

    def find_base_monomers_with_paths(self, complex_id):
        """
        Find all base monomers in a complex and trace the assembly path from monomer to complex.

        This function recursively decomposes a complex into its constituent monomers,
        tracking the stoichiometry and reaction steps at each level of assembly.

        Args:
            complex_id: The ID of the complex to analyze
        Returns:
            dict: {monomer_id: [{
                'path': [(molecule_id, stoich, type, reaction_id), ...],  # Ordered monomer → complex
                'total_stoich': int  # Total copies of this monomer per complex
            }]}
        """

        def get_direct_subunits(molecule_id, molecule_type):
            """
            Get the direct subunits (reactants) that form a given molecule.

            Args:
                molecule_id: The molecule whose subunits to find
                molecule_type: Type of molecule ('complex', 'eq_complex', or 'two component system')

            Returns:
                dict: {subunit_id: (stoichiometry, reaction_id)}
                    Where stoichiometry is positive (representing reactant consumption)
            """
            # Select the appropriate molecule-to-parent-complex dictionary:
            if molecule_type == 'complex':
                m2pc = self.complex_network.c_m2pc_dict  # Complexation process
            elif molecule_type == 'eq_complex':
                m2pc = self.complex_network.eq_m2pc_dict  # Equilibrium process
            elif molecule_type == 'two component system':
                m2pc = self.complex_network.tcs_m2pc_dict  # Two-component system process
            else:
                return {}

            # Find all molecules that have this molecule_id as a parent complex
            subunits = {}
            for mol, parents in m2pc.items():
                if molecule_id in parents:
                    # Extract reaction and stoichiometry info
                    reaction_id = parents[molecule_id][0]['reaction_id']
                    stoich = parents[molecule_id][1]['stoichiometry']
                    subunits[mol] = (np.negative(stoich), reaction_id)

            return subunits

        def explore_molecule(molecule_id, current_path, current_stoich, parent_stoich,
                             parent_reaction):
            """
            Recursively decomposes a molecule to find all constituent monomers
            by working through all assembly steps.

            Args:
                molecule_id: Current molecule being explored
                current_path: Path built so far, as list of (molecule_id, stoich, type, reaction_id)
                             This path goes from complex → subcomplex → ... (will be reversed at end)
                current_stoich: Cumulative stoichiometry - how many copies of this molecule
                               are needed per target complex
                parent_stoich: Stoichiometry of this molecule in its immediate
                            parent (the molecule one level up in the path)
                parent_reaction: Reaction ID that produces the parent from this molecule
            """
            mol_type = self.complex_network.determine_molecule_type(molecule_id)

            # Skip molecules that are not monomers or complexes, e.g. rRNAs
            # in the 30S/50S, metabolite ligands in equilibrium complexes:
            if mol_type == 'unknown':
                print(f"  → Skipping {molecule_id} (not a monomer or complex)")
                return

            # BASE CASE: the molecule is a fundamental monomer
            if mol_type == 'monomer':
                # Initialize storage of this monomer if 1st time encountering it
                if molecule_id not in monomer_paths:
                    monomer_paths[molecule_id] = []

                # Complete the path by adding the monomer itself
                # Use parent_stoich (not current_stoich) to show how many of this
                # monomer go into its immediate parent
                complete_path = current_path + [
                    (molecule_id, parent_stoich, mol_type, parent_reaction)]

                # Store the path in reverse order (monomer → intermediate → ... → complex)
                monomer_paths[molecule_id].append({
                    'path': complete_path[::-1],  # Reverse to go monomer → complex
                    'total_stoich': current_stoich  # Total copies needed per final complex
                })
                return

            # RECURSIVE CASE: the molecule is a complex itself, keep decomposing
            subunits = get_direct_subunits(molecule_id, mol_type)

            for subunit_id, (stoich, reaction_id) in subunits.items():
                # Add current molecule to the path with its parent stoichiometry
                # (how many of this molecule go into the level above it)
                new_path = current_path + [(molecule_id, parent_stoich, mol_type, parent_reaction)]

                # Update cumulative stoichiometry (need X of this molecule,
                # and each needs Y of the subunit, then need X*Y of the subunit total)
                new_stoich = current_stoich * stoich

                # Recursively explore this subunit and pass stoich as
                # parent_stoich for the next level:
                explore_molecule(subunit_id, new_path, new_stoich, stoich, reaction_id)

        # Initialize storage for all monomer paths that will be found:
        monomer_paths = {}

        # Determine what type of complex is being analyzed:
        complex_type = self.complex_network.determine_molecule_type(complex_id)

        if complex_type == 'unknown':
            print(f"Warning: {complex_id} is not a valid complex")
            return {}

        # Get the direct subunits of the target complex and start exploring each
        subunits = get_direct_subunits(complex_id, complex_type)
        for subunit_id, (stoich, reaction_id) in subunits.items():
            # Start with empty path (to be built upon):
            initial_path = []
            # Start with stoich as both current and parent stoich (will need
            # 'stoich' copies of this subunit per complex):
            explore_molecule(subunit_id, initial_path, stoich, stoich, reaction_id)

        return monomer_paths

    def format_monomer_path_sentence(self, monomer_id, complex_id, path_info):
        """
        Converts the path data structure into a readable sentence showing the
        assembly pathway and stoichiometry required at each step.

        Args:
            monomer_id: The base monomer ID
            complex_id: The target complex ID containing the monomer
            path_info: Dictionary with:
                - 'path': List of (molecule_id, stoich, mol_type, reaction_id) tuples
                         ordered from monomer → intermediate(s) → complex
                - 'total_stoich': Total copies of monomer per final complex
        Returns:
            str: Formatted sentence explaining the monomer counts within the
            complex and the assembly pathway.
        """
        path = path_info['path']
        total_stoich = path_info['total_stoich']

        # Build the assembly pathway string with stoichiometries
        path_parts = []

        # Add all molecules in the path with their stoichiometries
        # Each entry shows how many copies are needed at that assembly step
        for molecule_id, stoich, mol_type, reaction_id in path:
            path_parts.append(f"{stoich:g} {molecule_id}")

        # Add the final target complex (always 1 copy at the end of the path)
        path_parts.append(f"1 {complex_id}")

        # Format the path with line breaks for readability (break into chunks
        # of 2 molecules per line to avoid very long lines)
        path_string_parts = []
        for i in range(0, len(path_parts), 2):
            chunk = " → ".join(path_parts[i:i + 2])
            path_string_parts.append(chunk)
        path_string = "\n→ ".join(path_string_parts)

        # Construct the final sentence:
        sentence = (f"{monomer_id} counts within {complex_id}, "
                    f"{total_stoich:g} per\n(via {path_string})")

        return sentence

    def analyze_complex(self, complex_id):
        """
        Analyze a complex to identify all constituent monomers and their
        assembly pathways.

        Args:
            complex_id: The ID of the complex to analyze
        Returns:
            tuple: (monomer_paths, descriptions)
                - monomer_paths: Dict mapping monomer IDs to their path information
                    {monomer_id: [{
                        'path': [(molecule_id, stoich, type, reaction_id), ...],
                        'total_stoich': int
                    }]}
                - descriptions: List of formatted strings describing each pathway
        """
        print(f"\nAnalyzing complex: {complex_id}\n")

        # Find all monomers and their assembly pathways within the complex:
        monomer_paths = self.find_base_monomers_with_paths(complex_id)

        # Iterate through each monomer found in the complex:
        descriptions = []
        for monomer_id, path_list in monomer_paths.items():
            for path_info in path_list:
                sentence = self.format_monomer_path_sentence(
                    monomer_id, complex_id, path_info
                )
                descriptions.append(sentence)

        for desc in descriptions:
            print(desc)

        return monomer_paths, descriptions

    # Build the complexation parent and downstream molecule to complex dictionaries:
    def get_ribosome_subunit_counts(self, cell_paths):
        """
        Extract ribosome subunit counts from simulation data.
        Args:
            cell_paths: paths to the cell simulations to extract data from
        Returns:
            tuple: (total_active_ribosome_counts, total_ribosome_initiation_events,
                    total_ribosome_termination_events)
        """
        # Active ribosome counts need to be extracted from the Unique Molecules listener:
        ribosome_counts = []
        ribosome_initiation_events = []
        ribosome_termination_events = []
        ribosome_ID = 'active_ribosome'
        # Loop through the cells to extract the active ribosome counts @ each time point:
        for sim_dir in cell_paths:
            simOutDir = os.path.join(sim_dir, 'simOut')
            unique_molecule_counts_reader = TableReader(
                os.path.join(simOutDir, "UniqueMoleculeCounts"))
            ribosome_idx = unique_molecule_counts_reader.readAttribute(
                "uniqueMoleculeIds").index(ribosome_ID)
            active_ribosome_counts = unique_molecule_counts_reader.readColumn(
                "uniqueMoleculeCounts")[:, ribosome_idx]
            ribosome_counts.append(active_ribosome_counts)

            ribosome_data_reader = TableReader(os.path.join(simOutDir, "RibosomeData"))
            ribosomes_activated = ribosome_data_reader.readColumn("didInitialize")
            ribosomes_terminated = ribosome_data_reader.readColumn("didTerminate")

            # Add all the counts to the respective lists:
            ribosome_initiation_events.append(ribosomes_activated)
            ribosome_termination_events.append(ribosomes_terminated)

        # Concatenate the ribosome counts and events from each cell into one:
        total_active_ribosome_counts = np.concatenate(ribosome_counts)
        total_ribosome_initiation_events = np.concatenate(ribosome_initiation_events)
        total_ribosome_termination_events = np.concatenate(ribosome_termination_events)

        return (total_active_ribosome_counts,
                total_ribosome_initiation_events, total_ribosome_termination_events)

    def get_rnap_subunit_counts(self, metadata, cell_paths):
        """
        Extract RNAP subunit counts and events from simulation data.

        Args:
            metadata: simulation metadata
            cell_paths: cell paths to extract data from
        Returns:
            tuple: (total_active_rnap_counts, total_activation_events,
                    total_inactivation_events)
        """
        # Active replisome counts need to be extracted from the Unique Molecules listener:
        rnap_counts = []
        activation_events = []
        attenuation_events = []
        termination_events = []
        stalling_events = []
        rnap_ID = 'active_RNAP'
        r_s_e = metadata['recycle_stalled_elongation']

        # Loop through the cells to extract the RNAP counts and events:
        for sim_dir in cell_paths:
            simOutDir = os.path.join(sim_dir, 'simOut')
            unique_molecule_counts_reader = TableReader(
                os.path.join(simOutDir, "UniqueMoleculeCounts"))
            rnap_idx = unique_molecule_counts_reader.readAttribute(
                "uniqueMoleculeIds").index(rnap_ID)
            active_rnap_counts = unique_molecule_counts_reader.readColumn(
                "uniqueMoleculeCounts")[:, rnap_idx]
            rnap_counts.append(active_rnap_counts)

            rnap_data_reader = TableReader(os.path.join(simOutDir, "RnapData"))
            rnap_activated = rnap_data_reader.readColumn("didInitialize")

            # Active RNAP can be returned back to inactive via attenuation,
            # termination, or stalling:
            transcription_elongation_reader = TableReader(
                os.path.join(simOutDir, "TranscriptElongationListener"))
            rnap_attenuated = (
                transcription_elongation_reader.readColumn("counts_attenuated")).sum(
                axis=1)
            rnap_terminated = rnap_data_reader.readColumn("didTerminate")

            # Add all the counts to the respective lists:
            activation_events.append(rnap_activated)
            attenuation_events.append(rnap_attenuated)
            termination_events.append(rnap_terminated)

            # If recycle_stalled_elongation is turned on, count the stalled too:
            if r_s_e:
                rnap_stalled = rnap_data_reader.readColumn("didStall")
            else:
                rnap_stalled = np.zeros(len(rnap_activated))
            stalling_events.append(rnap_stalled)

        # Concatenate the RNAP counts and events from each cell into one array:
        total_active_rnap_counts = np.concatenate(rnap_counts)
        activation_events = np.concatenate(activation_events)
        attenuation_events = np.concatenate(attenuation_events)
        termination_events = np.concatenate(termination_events)
        stalling_events = np.concatenate(stalling_events)

        total_activation_events = activation_events
        total_inactivation_events = (
                attenuation_events + termination_events + stalling_events)

        return (total_active_rnap_counts,
                total_activation_events, total_inactivation_events)


    def get_replisome_subunit_counts(self, sim_data, cell_paths, molecule):
        """
        Extract replisome subunit counts and events from simulation data.
        Args:
            metadata: simulation metadata
            cell_paths: cell paths to extract data from
                molecule: replisome subunit molecule ID
        Returns:
            tuple: (replisome_subunit_stoich, replisome_counts, replisome_events)
        """
        # Active replisome counts need to be extracted from the Unique Molecules listener:
        replisome_counts = []
        replisome_events = []
        replisome_ID = 'active_replisome'
        # loop through the cell paths to extract the active replisome counts at each time point:
        for sim_dir in cell_paths:
            simOutDir = os.path.join(sim_dir, 'simOut')
            unique_molecule_counts_reader = TableReader(
                os.path.join(simOutDir, "UniqueMoleculeCounts"))
            replisome_idx = unique_molecule_counts_reader.readAttribute(
                "uniqueMoleculeIds").index(replisome_ID)
            active_replisome_counts = unique_molecule_counts_reader.readColumn(
                "uniqueMoleculeCounts")[:, replisome_idx]
            replisome_counts.append(active_replisome_counts)

            # Replisome assembly events are the replisomes formed in each step
            # (see count_replisomes_formed_and_removed), so inherited
            # replisomes are not counted as assembly at each generation start
            # and terminations in the same step can't cancel them out:
            replisomes_formed, _ = count_replisomes_formed_and_removed(sim_dir)
            replisome_events.append(replisomes_formed)

        # Concatenate the replisome counts and events from each simulation into one array:
        replisome_counts = np.concatenate(replisome_counts)
        replisome_events = np.concatenate(replisome_events)

        return (replisome_subunit_stoich(sim_data, molecule),
                replisome_counts, replisome_events)


    def do_plot(self, seedOutDir, plotOutDir, plotOutFileName, simDataFile,
                validationDataFile, metadata):
        with open(simDataFile, 'rb') as f:
            sim_data = pickle.load(f)
        # Only the first-daughter lineage is plotted:
        cell_paths = get_first_daughter_lineage(self.ap.get_cells())
        sim_dir = cell_paths[0]
        simOutDir = os.path.join(sim_dir, 'simOut')

        # Check if replication sequestered replisome subunits. With mechanistic
        # replisome off (the default setting), chromosome replication never
        # removes replisome subunits from the bulk pool, so active replisomes
        # contain no real subunits. If it is on, then subunits are pulled from
        # bullk to generate replilsomes in the unique molecule space.
        mechanistic_replisome = metadata.get('mechanistic_replisome', False)

        # Extract protein indexes for each monomer:
        monomer_counts_reader = TableReader(
            os.path.join(simOutDir, "MonomerCounts"))
        monomer_idx_dict = {monomer: i for i, monomer in
                            enumerate(monomer_counts_reader.readAttribute(
                                'monomerIds'))}

        # Map each monomer to its gene's common name (e.g. dnaE) for the
        # subunit availability legend:
        monomer_data = sim_data.process.translation.monomer_data
        monomer_to_gene_name = {
            monomer: get_clean_common_name(sim_data, cistron)
            for monomer, cistron in zip(monomer_data['id'], monomer_data['cistron_id'])}

        # Extract the free monomer counts:
        free_monomer_counts = read_stacked_columns(
            cell_paths, 'MonomerCounts', "freeMonomerCounts")

        # Extract the doubling times:
        (time, doubling_times, end_generation_times, start_generation_indices,
         end_generation_indices) = extract_doubling_times(cell_paths)

        # Build the molecule to parent complex dictionaries of the three types
        # of complexes:
        self.complex_network = ComplexNetwork(sim_data)

        # Estimate the TCS reaction events:
        self.tcs_reaction_events, _ = build_tcs_reaction_event_estimates(
            sim_data, read_stacked_columns(
                cell_paths, 'TwoComponentSystems', 'moleculeChanges'))

        # Check that all the valid inputs are indeed valid and complexes:
        valid_complexes, complex_type_dict = (
            self.check_complex_validity(sim_data, PLOT_COMPLEXES))

        # Generate a plot for each complex:
        for complex_id in valid_complexes:
            # Extract the upstream monomers and sentences describing their paths:
            upstream_monomers, _ = self.analyze_complex(complex_id)

            # Keep each distinct path from a monomer into the complex only. As
            # of 10/24/2026, CPLX0-3803 is the only complex with more than one
            # path for a monomer: EG11413-MONOMER and EG11414-MONOMER are in it
            # both via CPLX0-7910 directly, and via CPLX0-3801 (and each route
            # contributes its own copies of the monomers). Currently CPLX0-3803
            # never actually forms because its DnaX subunit (EG10245-MONOMER)
            # is never made (since translation.py only translates the first
            # protein listed for each cistron and dnaX lists the frameshift
            # product first):
            unique_upstream_paths = {}
            for monomer, path_list in upstream_monomers.items():
                unique_paths = []
                for path_info in path_list:
                    if not any(path_info['path'] == unique_path['path']
                               and path_info['total_stoich'] == unique_path['total_stoich']
                               for unique_path in unique_paths):
                        unique_paths.append(path_info)
                unique_upstream_paths[monomer] = unique_paths
            upstream_path_pairs = [
                (monomer, path_info)
                for monomer, path_list in unique_upstream_paths.items()
                for path_info in path_list]

            # Extract all downstream complexes the complex is a subunit of:
            downstream_complexes, _ = self.complex_network.analyze_molecule_parents(complex_id)

            unique_downstream_paths = downstream_complexes
            parent_path_pairs = [
                (parent_complex, path_info)
                for parent_complex, path_list in unique_downstream_paths.items()
                for path_info in path_list]

            # Copies of the complex per parent complex, summed over all
            # distinct paths:
            parent_total_stoich = {
                parent_complex: sum(path_info['total_stoich'] for path_info in path_list)
                for parent_complex, path_list in unique_downstream_paths.items()}

            # Determine the complex type (complexation vs. equilibrium vs. TCS):
            complex_type = complex_type_dict[complex_id]

            # Determine the complex makeup (homogeneous vs. heterogeneous):
            complex_makeup = self.determine_complex_makeup(complex_id, complex_type)

            # Extract the counts of the complex from bulk molecule counts:
            complex_counts = (
                read_stacked_bulk_molecules(cell_paths, [complex_id]))[0]

            # Track the theoretical total counts of the complex (i.e., free
            # counts + counts within other complexes + counts within unique molecules)
            total_counts = complex_counts.copy().astype(float)

            # Generate the plots:
            if is_a_ribosome_subunit(sim_data, complex_id) == []:
                # Normal plot size for non-ribosome related complexes:
                fig, (ax1, ax2, ax3) = (
                    plt.subplots(3, 1, sharex=True, figsize=(10, 6)))
            else:
                # Make the plot larger for ribosome subunits since there are
                # so many monomer subunits to plot for the ribosome subunits:
                fig, (ax1, ax2, ax3) = (
                    plt.subplots(3, 1, sharex=True, figsize=(30, 20)))


            # First, add the free counts of the complex over time to the first plot:
            ax1.plot(time, complex_counts, color='lightseagreen',
                     label=f'{complex_id} free complex counts', linewidth=.75, alpha=0.75)

            # Also plot the counts of each base monomer within the complex
            # (plot one line per distinct path, labeled with that path):
            for monomer, path_list in unique_upstream_paths.items():
                for path_num, path_info in enumerate(path_list, 1):
                    path = self.format_monomer_path_sentence(
                        monomer, complex_id, path_info)
                    if len(path_list) > 1:
                        path += f'\n[Path {path_num} of {len(path_list)}]'
                    monomer_counts = complex_counts * path_info['total_stoich']

                    # Plot the counts within the complex itself:
                    ax1.plot(time, monomer_counts, alpha=1,
                             label=path,
                             linewidth=1, linestyle=":")

                # Second, plot the free monomer subunits on the available free
                # monomer subunit counts plot:
                protein_idx = monomer_idx_dict[monomer]
                protein_FMC = free_monomer_counts[:, protein_idx]
                gene_name = monomer_to_gene_name.get(monomer)
                ax2.plot(time, protein_FMC, alpha=0.75,
                         label=f'{monomer}' + (f' ({gene_name})' if gene_name else ''),
                         linewidth=0.75)

            # Add the parent complexes to the complex counts plot as well:
            # (one line per distinct path, labeled with that path):
            complexes = list(unique_downstream_paths.keys())
            for parent_complex, path_list in unique_downstream_paths.items():
                parent_complex_free_counts = (
                    read_stacked_bulk_molecules(cell_paths, [parent_complex]))[0]
                for path_num, path_info in enumerate(path_list, 1):
                    total_stoich = path_info['total_stoich']
                    path = format_parent_path_sentence(
                        complex_id, parent_complex, path_info)
                    if len(path_list) > 1:
                        path += f'\n[Path {path_num} of {len(path_list)}]'
                    parent_complex_counts = parent_complex_free_counts * total_stoich
                    ax1.plot(time, parent_complex_counts, alpha=0.5,
                             label=path,
                             linewidth=0.75, linestyle="-.")
                    total_counts += parent_complex_counts

            # Third, add a plot of the complexation events over time:
            rxns_plotted = []
            consumption_rxns_plotted = []

            # Define the differnt complex reactions for easy reference below:
            complexation_rxn_ids = sim_data.process.complexation.ids_reactions
            equilibrium_rxn_ids = sim_data.process.equilibrium.rxn_ids
            tcs_rxn_ids = sim_data.process.two_component_system.rxn_ids

            # Plot the reactions that generate the complex:
            for monomer, path_info in upstream_path_pairs:
                # Find the reaction that makes the complex from this monomer
                # (the last step of each distinct path):
                rxn_id = path_info['path'][-1][3]
                if rxn_id in rxns_plotted:
                    continue # skip if the reacion has been plotted
                else:
                    rxns_plotted.append(rxn_id)
                rxn_stoich = 1 # always 1 for complex formation
                name = f'{rxn_id} \n(generates {rxn_stoich:g} {complex_id})'

                # Plot the reaction events based on the type of complex:
                if complex_type == "complexation":
                    rxn_idx = complexation_rxn_ids.index(rxn_id)
                    rxn_events = read_stacked_columns(cell_paths,
                        'ComplexationListener', "complexationEvents")[:, rxn_idx]
                    ax3.plot(time, rxn_events, alpha=0.5, label=name, linewidth=0.75)
                elif complex_type == "equilibrium":
                    rxn_idx = equilibrium_rxn_ids.index(rxn_id)
                    rxn_events = read_stacked_columns(cell_paths,
                        'EquilibriumListener', "complexationEvents")[:, rxn_idx]
                    # Equilibrium events are net per step, so the same line
                    # generates (formation) or dissociates the complex:
                    ax3.plot(time, rxn_events, alpha=0.5, linewidth=0.75,
                             label=f'{rxn_id} net events\n(generates/dissociates '
                                   f'{rxn_stoich:g} {complex_id})')
                else: # TCS reactions
                    rxn_idx = tcs_rxn_ids.index(rxn_id)
                    rxn_events = self.tcs_reaction_events[:, rxn_idx]
                    ax3.plot(time, rxn_events, alpha=0.5, label=
                            f'{rxn_id} \n estimated events '
                            f'(generates 1 {complex_id})', linewidth=0.75)

                    # Also check if the complex is made from a different TCS
                    # reaction by the same base monomer:
                    if "POS" in rxn_id:
                        reaction_IDs = (self.complex_network.find_extra_tcs_reactions(
                            monomer, complex_id, rxns_plotted))
                        for reaction_ID in reaction_IDs:
                            print(f"Plotting extra reaction {reaction_ID} for "
                                  f"{complex_id} formation from {monomer}")
                             # Plot the reaction events based on the type of complex:
                            rxns_plotted.append(reaction_ID)
                            rxn_idx = tcs_rxn_ids.index(reaction_ID)
                            rxn_events = self.tcs_reaction_events[:, rxn_idx]
                            ax3.plot(time, rxn_events, alpha=0.5, label=
                            f'{reaction_ID} \n estimated events '
                            f'(generates 1 {complex_id})', linewidth=0.75)

                    # Lastly, if there is a TCS reaction in the mix, check if
                    # the complex itself is consumed by any TCS reactions:
                    if complex_id in self.complex_network.tcs_c2m_via_dprs_dict.keys():
                        subunit_complex_info = self.complex_network.tcs_c2m_via_dprs_dict[complex_id]
                        for subunit_complex in subunit_complex_info:
                            # Extract the reverse reaction reactant's information:
                            subunit_complex_dict = subunit_complex_info[subunit_complex]
                            rxn_id = subunit_complex_dict[0]['reaction_id']
                            if rxn_id in rxns_plotted or rxn_id in consumption_rxns_plotted:
                                continue
                            # Otherwise, plot the dephosphorylation reaction:
                            print(f"Plotting dephosphorylation reaction {rxn_id}"
                                  f" that consumes {complex_id}")
                            consumption_rxns_plotted.append(rxn_id)
                            rxn_idx = tcs_rxn_ids.index(rxn_id)
                            rxn_events = self.tcs_reaction_events[:, rxn_idx]
                            ax3.plot(time, np.negative(rxn_events), alpha=0.5,
                                     label=f'{rxn_id} \n estimated events '
                                    f'(consumes 1 {complex_id})', linewidth=0.75)


            # Next, Find the reactions that consume the complex as a subunit
            # of a larger complex, and plot those as well:
            # (the first step of each distinct path):
            for parent_complex, path_info in parent_path_pairs:
                rxn_id = path_info['path'][0][3]
                if rxn_id in consumption_rxns_plotted or rxn_id in rxns_plotted:
                    continue # skip if the reacion has been plotted
                else:
                    consumption_rxns_plotted.append(rxn_id)
                rxn_stoich = path_info['path'][0][1]
                name = f'{rxn_id} \n(consumes {rxn_stoich:g} {complex_id})'

                if rxn_id in complexation_rxn_ids:
                    rxn_idx = complexation_rxn_ids.index(rxn_id)
                    rxn_events = read_stacked_columns(cell_paths,
                        'ComplexationListener', "complexationEvents")[:, rxn_idx]
                    ax3.plot(time, np.negative(rxn_events),
                             alpha=0.5, label=name, linewidth=0.75)

                elif rxn_id in equilibrium_rxn_ids:
                    rxn_idx = equilibrium_rxn_ids.index(rxn_id)
                    rxn_events = read_stacked_columns(cell_paths,
                        'EquilibriumListener', "complexationEvents")[:, rxn_idx]
                    # Equilibrium events are net per step, so the same line
                    # consumes (parent formation) or releases (dissociation):
                    ax3.plot(time, np.negative(rxn_events), alpha=0.5,
                             linewidth=0.75,
                             label=f'{rxn_id} net events\n(consumes/releases '
                                   f'{rxn_stoich:g} {complex_id})')

                else:
                    rxn_idx = tcs_rxn_ids.index(rxn_id)
                    rxn_events = self.tcs_reaction_events[:, rxn_idx]
                    ax3.plot(time, np.negative(rxn_events),
                             alpha=0.5, label=f'{rxn_id} \n estimated events '
                                 f'(consumes 1 {complex_id})', linewidth=0.75)

                    # Also check if the parent is made from a different TCS
                    # reaction by the same base complex:
                    if "POS" in rxn_id:
                        reaction_IDs = self.complex_network.find_extra_tcs_reactions(
                            complex_id, parent_complex, consumption_rxns_plotted)
                        for reaction_ID in reaction_IDs:
                            # Plot the reaction events based on the type of complex:
                            consumption_rxns_plotted.append(reaction_ID)
                            rxn_idx = tcs_rxn_ids.index(reaction_ID)
                            rxn_events = self.tcs_reaction_events[:, rxn_idx]
                            ax3.plot(time, np.negative(rxn_events), alpha=0.5, label=
                            f'{reaction_ID} \n estimated events (consumes 1 {complex_id})', linewidth=0.75)

                    # Lastly, check if the parent dephosphorylates to the complex:
                    if parent_complex in self.complex_network.tcs_c2m_via_dprs_dict.keys():
                        subunit_complex_info = self.complex_network.tcs_c2m_via_dprs_dict[parent_complex]
                        if complex_id in subunit_complex_info.keys():
                            subunit_complex_dict = subunit_complex_info[complex_id]
                            rxn_id = subunit_complex_dict[0]['reaction_id']
                            if rxn_id in rxns_plotted or rxn_id in consumption_rxns_plotted:
                                continue
                            # Otherwise, plot the dephosphorylation reaction:
                            rxns_plotted.append(rxn_id)
                            rxn_idx = tcs_rxn_ids.index(rxn_id)
                            rxn_events = self.tcs_reaction_events[:, rxn_idx]
                            ax3.plot(time, rxn_events, alpha=0.5, label=
                                f'{rxn_id}\n estimated events (generates 1 {complex_id})', linewidth=0.75)

            # Finally, check if the complex (or its parent complexes) are
            # involved in any special events that should be plotted:

            # Check if the complex is a TF:
            if is_a_TF(sim_data, complex_id) != []:
                # If the complex is a TF, plot its bound counts and binding/unbinding events:
                tf_idx = tf_index(sim_data, complex_id)
                tf_counts = read_stacked_columns(cell_paths,'RnaSynthProb',
                                                "nActualBound")[:,tf_idx]
                tfs_unbound = read_stacked_columns(cell_paths,'RnaSynthProb',
                                                "nActualUnbound")[:,tf_idx]
                ax1.plot(time, tf_counts, color='lightcoral',
                         label=f'transcription unit bound\n {complex_id} '
                               f'transcription factors',
                         linewidth=0.75, alpha=0.75)
                # nActualBound is TF occupancy at the end of each step, which
                # works as binding events because TfUnbinding releases every
                # bound TF at the start of each step. The exception is the
                # first row/timestep of each generation, which is the initial state
                # inherited from the parent cell (with no matching unbinding
                # row), so those rows are set to 0 to avoid a fake drop at
                # every generation start:
                tf_binding = tf_counts.astype(float)
                tf_binding[start_generation_indices] = 0
                ax3.plot(time, (tf_binding * -1), color='lightcoral',
                         label=f'{complex_id} Transcription factor\nbinding events',
                         linewidth=0.75, alpha=0.75)
                ax3.plot(time, tfs_unbound, color='purple',
                         label=f'{complex_id} Transcription factor\nunbinding events',
                         linewidth=0.75, alpha=0.75)
                total_counts += tf_counts

            # Check if the complex is a replisome subunit:
            if is_a_replisome_subunit(sim_data, complex_id) != []:
                # If the complex is a replisome subunit, plot its counts within
                # active replisomes and the replisome assembly events that consume it:
                replisome_stoich, replisome_counts, replisome_events = (
                    self.get_replisome_subunit_counts(sim_data, cell_paths, complex_id))
                subunit_counts_within_replisomes = (
                        replisome_counts * replisome_stoich)
                if mechanistic_replisome:
                    ax1.plot(time, subunit_counts_within_replisomes, color='gold',
                             label=f'{complex_id} counts within active replisomes, '
                                   f'{replisome_stoich:g} per',
                             linewidth=0.75, alpha=0.75, linestyle="--")
                    ax3.plot(time, np.negative(replisome_events), color='gold',
                             label=f'Active replisome assembly events\n(consumes '
                                   f'{replisome_stoich:g} {complex_id})',
                             linewidth=0.75, alpha=0.75)
                    # Add the replisome subunit counts to the total counts of the complex:
                    total_counts += (subunit_counts_within_replisomes)
                else:
                    # Not real protein: not added to the total counts, and
                    # replisome assembly never consumes the free complex
                    ax1.plot(time, subunit_counts_within_replisomes, color='gold',
                             label=f'# of {complex_id} chromosome replication\n'
                                   f'assumes are available to make active\n'
                                   f'replisomes, but mechanistic_replisome is\n'
                                   f'off so they are never actually removed\n'
                                   f'from bulk ({replisome_stoich:g} per)',
                             linewidth=0.75, alpha=0.75, linestyle=":")

            # Check if the complex is an inactive RNAP (1:1 ratio):
            if is_an_RNAP_subunit(sim_data, complex_id) != []:
                rnap_counts, rnap_activation_events, rnap_deactivation_events = (
                    self.get_rnap_subunit_counts(metadata, cell_paths))
                ax1.plot(time, rnap_counts, color='gold',
                         label=f'{complex_id} counts within active RNAPs, 1 per',
                         linewidth=0.75, alpha=0.75, linestyle="--")
                ax3.plot(time, rnap_activation_events * -1, color='gold',
                         label=f'RNAP activation events\n(consumes 1 {complex_id})',
                         linewidth=0.75, alpha=0.75)
                ax3.plot(time, rnap_deactivation_events, color='orange',
                         label=f'RNAP inactivation events\n(produces 1 {complex_id})',
                         linewidth=0.75, alpha=0.75)
                total_counts += rnap_counts

            # Check if the complex is a ribosome subunit (1:1 ratio):
            if is_a_ribosome_subunit(sim_data, complex_id) != []:
                ribosome_counts, ribosome_initiation_events, ribosome_termination_events = (
                    self.get_ribosome_subunit_counts(cell_paths))
                ax1.plot(time, ribosome_counts, color='gold',
                         label=f'{complex_id} counts within active ribosomes, 1 per',
                         linewidth=0.75, alpha=0.75, linestyle="--")
                ax3.plot(time, ribosome_initiation_events * -1, color='gold',
                         label=f'Active ribosome assembly events\n(consumes 1 {complex_id})',
                         linewidth=0.75, alpha=0.75)
                ax3.plot(time, ribosome_termination_events, color='orange',
                         label=f'Active ribosome disassembly events\n(produces 1 {complex_id})',
                         linewidth=0.75, alpha=0.75)
                total_counts += ribosome_counts

            # Parent complex checks:
            for i in range(len(complexes)):
                parent_complex = complexes[i]
                if is_a_replisome_subunit(sim_data, parent_complex) != []:
                    replisome_stoich, replisome_counts, _ = (
                        self.get_replisome_subunit_counts(
                            sim_data, cell_paths, parent_complex))

                    # Get the stoich of the parent complex to the complex of interest:
                    total_stoich = parent_total_stoich[parent_complex]

                    # Multiply the replisome subunit counts by the stoich of
                    # the parent complex and the stoich of the parent complex in the replisome:
                    subunit_counts = replisome_counts * total_stoich * replisome_stoich
                    if mechanistic_replisome:
                        ax1.plot(time, subunit_counts, color='gold',
                                 label=f'{complex_id} counts within active replisomes,'
                                       f'\n{replisome_stoich*total_stoich:g} per '
                                       f'(via {parent_complex})',
                                 linewidth=0.75, alpha=0.75, linestyle="--")

                        # Add the replisome counts to the total counts of the complex:
                        total_counts += subunit_counts
                    else:
                        # monomers are not actually used to create replisomes
                        # and mechanistic replisome assumes they are unlimited:
                        ax1.plot(time, subunit_counts, color='gold',
                                 label=f'# of {complex_id} chromosome replication\n'
                                       f'assumes are available to make active\n'
                                       f'replisomes, but mechanistic_replisome is\n'
                                       f'off so they are never actually removed\n'
                                       f'from bulk ({replisome_stoich*total_stoich:g} '
                                       f'per, via {parent_complex})',
                                 linewidth=0.75, alpha=0.75, linestyle=":")

                if is_a_TF(sim_data, parent_complex) != []:
                    tf_idx = tf_index(sim_data, parent_complex)
                    tf_counts = read_stacked_columns(cell_paths,
                        'RnaSynthProb',"nActualBound")[:,tf_idx]

                    # Get the stoich of the parent complex to the complex of interest:
                    total_stoich = parent_total_stoich[parent_complex]

                    # Multiply the TF counts by the stoich of the complex to
                    # the parent complex:
                    subunit_counts = tf_counts * total_stoich
                    ax1.plot(time, subunit_counts, color='lightcoral',
                             label=f'bound {complex_id} transcription'
                                   f' factors,\n {total_stoich:g} per (via {parent_complex})',
                             linewidth=0.75, alpha=0.75, linestyle="--")

                    # Add the TF counts to the total counts of the complex:
                    total_counts += subunit_counts

                # No need to check RNAPs, they only have monomer subunits.

                if is_a_ribosome_subunit(sim_data, parent_complex) != []:
                    ribosome_counts, _, _ = (
                        self.get_ribosome_subunit_counts(cell_paths))

                    # Get the stoich of the parent complex to the complex of interest:
                    total_stoich = parent_total_stoich[parent_complex]

                    # Multiply the ribosome subunit counts by the stoich of
                    # the parent complex and the stoich of the parent complex in the ribosome:
                    subunit_counts = ribosome_counts * total_stoich

                    ax1.plot(time, subunit_counts, color='gold',
                             label=f'{complex_id} counts within active ribosomes,'
                                   f'\n{total_stoich:g} per (via {parent_complex})',
                             linewidth=0.75, alpha=0.75, linestyle="--")

                    # Add the ribosome counts to the total counts of the complex:
                    total_counts += subunit_counts


            # Plot the total theoretical complex counts behind the other lines
            # (so a free count that flips every step doesn't get covered).
            # If the total equals the free count at every step, keep the
            # legend entry but don't draw it on top of the free line:
            if np.array_equal(total_counts, complex_counts):
                ax1.plot([], [], color='black',
                         label=f'Theoretical total {complex_id} counts\n(not '
                               f'drawn here since free = total '
                               f'complex counts at every step)',
                         linewidth=.25, alpha=1)
            else:
                ax1.plot(time, total_counts, color='black',
                         label=f'Theoretical total {complex_id} counts\n(free complex '
                               f'counts + counts within parent complexes)',
                         linewidth=.25, alpha=1, zorder=0)

            # Set the axes and legends :
            ax1.legend(fontsize=5, loc="center left", bbox_to_anchor=(1, 0.5))
            ax1.set_ylabel("Complex Counts")

            # Create stats above the counts panel ("in other complexes" is
            # everything in the theoretical total that is not free):
            in_other_complexes = total_counts - complex_counts
            ax1.set_title(
                f"Mean counts - Total: {np.mean(total_counts):.1f} | "
                f"Free: {np.mean(complex_counts):.1f} | "
                f"In other complexes: {np.mean(in_other_complexes):.1f} | "
                f"Median counts - Total: {np.median(total_counts):.1f} | "
                f"Free: {np.median(complex_counts):.1f} | "
                f"In other complexes: {np.median(in_other_complexes):.1f}",
                fontsize=7)

            # Format the title a bold line drawn separately above grey subtitle
            # lines. The titles for the TCS role labels are built from
            # build_tcs_dictionaries (and this assumes only the phosphorylated
            # forms are TCS complexes themselves):
            tcs_role_names = {
                "PHOSPHO-RR": "TCS phosphorylated response regulator",
                "PHOSPHO-HK": "TCS phosphorylated histidine kinase",
                "PHOSPHO-HK-LIGAND": "TCS phosphorylated ligand-bound histidine kinase",
            }
            if complex_makeup in tcs_role_names:
                classification = tcs_role_names[complex_makeup]
            else:
                classification = f"{complex_makeup} {complex_type} complex"
            tf_type = sim_data.process.transcription_regulation.tf_to_tf_type.get(
                complex_id.split('[')[0])
            if tf_type is not None:
                classification += f" | {tf_type} transcription factor"
            common_name = get_clean_common_name(sim_data, complex_id)
            if common_name:
                classification += f" | {common_name}"
            fig_height = fig.get_figheight()
            fig.text(0.4, 1 - 0.024 / fig_height,
                     f"Count trajectory for {complex_id}",
                     fontsize=10, fontweight='bold', ha='center', va='top')
            fig.suptitle(
                f"{classification[0].upper() + classification[1:]}\n"
                f"Sim ID: {metadata['description']} | "
                f"variant {metadata.get('variant', '?')} {metadata.get('variant_index', '')} | "
                f"seed {int(metadata.get('seed', 0))} | {len(cell_paths)} gen(s)",
                fontsize=9, x=0.4, y=1 - 0.288 / fig_height, color='grey')
            ax2.set_ylabel("Free Monomer \nSubunit Availability")
            ax2.legend(fontsize=5, loc="center left", bbox_to_anchor=(1, 0.5))

            # In the first generation, equilibrium (1CS ligand binding) and TCS
            # reactions do not execute during initialization to create
            # complexes for the initial state the way complexation complexes
            # do, so there are usually large spikes for these molecules at the
            # start. Later generations can also have smaller spikes as these
            # reactions relax after division (i.e. a +/- pair at the 3rd and
            # 4th steps). To avoid obscuring the y axis on the events graph
            # from these changes, these large spikes are truncated and labeled
            # instead:
            ordered_event_lines = format_event_lines(
                ax3, time, start_generation_indices, N_INITIAL_STEPS_TO_CLIP)
            ax3.axhline(y=0, color='black', linestyle='-', linewidth=0.5, alpha=0.5)

            ax3.set_xlabel("Time (s)")
            ax3.set_ylabel("Free Complex\nReaction Events")
            ax3.legend(handles=ordered_event_lines, fontsize=5, loc="center left",
                       bbox_to_anchor=(1, 0.5))

            # Add vertical lines for the end of each generation:
            for i in range(len(end_generation_times)):
                dt = end_generation_times[i]
                a = 0.5
                ax1.axvline(x=dt, linestyle='--', color="yellowgreen", alpha=a)
                ax2.axvline(x=dt, linestyle='--', color="yellowgreen", alpha=a)
                ax3.axvline(x=dt, linestyle='--', color="yellowgreen", alpha=a)

            if is_a_ribosome_subunit(sim_data, complex_id) == []:
                plt.tight_layout()
            # Leave room at the top for the two title blocks:
            plt.subplots_adjust(top=1 - 0.95 / fig_height)

            # Save the plot:
            file_name = plotOutFileName + "_" + complex_id
            exportFigure(plt, plotOutDir, file_name, metadata, dpi=300)
            plt.close()


if __name__ == '__main__':
    Plot().cli()
