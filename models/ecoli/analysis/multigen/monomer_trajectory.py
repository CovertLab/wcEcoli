"""
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
import scipy.optimize
from models.ecoli.analysis import multigenAnalysisPlot
from wholecell.analysis.analysis_tools import (exportFigure,
         read_stacked_bulk_molecules, read_stacked_columns)
from wholecell.io.tablereader import TableReader
from wholecell.utils import units

# USER INPUTS:

# List of proteins to be plotted (with or without compartment tags):
PLOT_MONOMERS = ["MONOMER0-4539", "MONOMER0-4538", "EG10504-MONOMER",
                 "BAES-MONOMER[i]", "BAER-MONOMER[c]"]

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

def count_replisomes_formed_and_removed(cell_path):
    """
    Number of replisomes formed and removed at each timestep of one cell, from
    the ReplicationData listener's fork_unique_index (the unique index of each
    active replisome, NaN-padded). Comparing which replisomes exist at
    successive timesteps counts initiation and termination separately, even
    when both happen in the same timestep, where the net change in the
    active_replisome count would cancel them out. The first timestep has no
    events, since the cell inherited those replisomes.
    """
    fork_indices = TableReader(
        os.path.join(cell_path, 'simOut', 'ReplicationData')).readColumn(
        'fork_unique_index', squeeze=False)
    forks_by_step = [set(row[~np.isnan(row)].astype(np.int64))
                     for row in fork_indices]
    formed = np.zeros(len(forks_by_step))
    removed = np.zeros(len(forks_by_step))
    for i in range(1, len(forks_by_step)):
        formed[i] = len(forks_by_step[i] - forks_by_step[i - 1])
        removed[i] = len(forks_by_step[i - 1] - forks_by_step[i])
    return formed, removed


class Plot(multigenAnalysisPlot.MultigenAnalysisPlot):
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
                # If the molecule is not valid print a warning message:
                print(f"{molecule} is not a valid molecule in the simulation.")

        return revised_molecule_list

    def generate_subunit_categories(self, sim_data):
        """
        Generate dictionaries to classify monomer subunits into categories
        based on their participation in different complexes
        (macromolecule/complexation, equilibrium, TCS) and actively modeled
        processes in the simulation (TFs, ribosomes, replisomes).
        """
        protein_ids = sim_data.process.translation.monomer_data['id']
        # Create monomer dict for indexing
        monomer_dict = {mol: i for i, mol in enumerate(protein_ids)}

        def get_monomer_indexes(keys):
            indices = []
            for x in keys:
                if x in monomer_dict:
                    indices.append(monomer_dict[x])
            return np.array(indices)

        # COMPLEXATION COMPLEXES:
        complexation_molecule_ids = sim_data.process.complexation.molecule_names
        complexation_complex_ids = sim_data.process.complexation.ids_complexes
        complexation_smm = sim_data.process.complexation.stoich_matrix_monomers()
        complexes_to_monomers = {}
        complexation_monomers = {}

        for id in complexation_molecule_ids:
            if id in monomer_dict:
                complexation_monomers[id] = monomer_dict[id]
            elif id in complexation_complex_ids:
                # Find constituent monomers of this complex
                complex_idx = complexation_complex_ids.index(id)
                monomer_idxs = np.where(complexation_smm[:, complex_idx] < 0)[0]
                monomers = []
                for idx in monomer_idxs:
                    monomer_id = complexation_molecule_ids[idx]
                    if monomer_id in monomer_dict:
                        complexation_monomers[monomer_id] = monomer_dict[monomer_id]
                        monomers.append(monomer_id)
                complexes_to_monomers[id] = monomers

        # EQUILIBRIUM COMPLEXES:
        equilibrium_molecule_ids = sim_data.process.equilibrium.molecule_names
        equilibrium_complex_ids = sim_data.process.equilibrium.ids_complexes
        eq_smm = sim_data.process.equilibrium.stoich_matrix_monomers()
        eq_complexes_to_monomers = {}
        equilibrium_monomers = {}

        for id in equilibrium_molecule_ids:
            if id in monomer_dict:
                equilibrium_monomers[id] = monomer_dict[id]
            elif id in complexation_complex_ids:
                # Use complexation monomers
                if id in complexes_to_monomers:
                    for monomer_id in complexes_to_monomers[id]:
                        if monomer_id in monomer_dict:
                            equilibrium_monomers[monomer_id] = monomer_dict[monomer_id]
            elif id in equilibrium_complex_ids:
                # Find constituent molecules of this equilibrium complex
                eq_complex_idx = equilibrium_complex_ids.index(id)
                molecule_idxs = np.where(eq_smm[:, eq_complex_idx] < 0)[0]
                monomers = []
                for idx in molecule_idxs:
                    molecule_name = equilibrium_molecule_ids[idx]
                    if molecule_name in monomer_dict:
                        monomers.append(molecule_name)
                        equilibrium_monomers[molecule_name] = monomer_dict[molecule_name]
                    elif molecule_name in complexation_complex_ids:
                        # Expand complexation complexes to monomers
                        if molecule_name in complexes_to_monomers:
                            for monomer_id in complexes_to_monomers[molecule_name]:
                                if monomer_id in monomer_dict:
                                    monomers.append(monomer_id)
                                    equilibrium_monomers[monomer_id] = monomer_dict[monomer_id]
                eq_complexes_to_monomers[id] = monomers

        # TCS COMPLEXES:
        two_component_system_molecule_ids = (
            list(sim_data.process.two_component_system.modified_molecules))
        two_component_system_complex_ids = (
            list(sim_data.process.two_component_system.complex_to_monomer.keys()))
        tcs_complex_to_monomer_dict = sim_data.process.two_component_system.complex_to_monomer
        tcs_complexes_to_monomers = {}
        two_component_system_monomers = {}

        for id in two_component_system_molecule_ids:
            if id in protein_ids:
                two_component_system_monomers[id] = get_monomer_indexes([id])[0]
            elif id in complexation_complex_ids:
                monomer_ids = complexes_to_monomers[id]
                tcs_complexes_to_monomers[id] = monomer_ids
            elif id in equilibrium_complex_ids:
                monomer_ids = eq_complexes_to_monomers[id]
                tcs_complexes_to_monomers[id] = monomer_ids
            elif id in two_component_system_complex_ids:
                monomer_ids = tcs_complex_to_monomer_dict[id]
                tcs_complexes_to_monomers[id] = monomer_ids.keys()

        # RIBOSOMES:
        ribosome_50s_subunits = sim_data.process.complexation.get_monomers(
            sim_data.molecule_ids.s50_full_complex)
        ribosome_30s_subunits = sim_data.process.complexation.get_monomers(
            sim_data.molecule_ids.s30_full_complex)
        ribosome_subunits = (ribosome_50s_subunits["subunitIds"].tolist() +
                             ribosome_30s_subunits["subunitIds"].tolist())
        ribosome_monomers = {}
        for subunit in ribosome_subunits:
            if subunit in protein_ids:
                ribosome_monomers[subunit] = get_monomer_indexes([subunit])[0]

        # RNAPs:
        rnap_subunits = sim_data.process.complexation.get_monomers(
            sim_data.molecule_ids.full_RNAP)
        rnap_subunit_ids = rnap_subunits["subunitIds"].tolist()
        rnap_monomers = {}
        for subunit in rnap_subunit_ids:
            if subunit in protein_ids:
                rnap_monomers[subunit] = get_monomer_indexes([subunit])[0]
            if subunit in complexation_complex_ids:
                monomer_ids = complexes_to_monomers[subunit]
                for id in monomer_ids:
                    if id in protein_ids:
                        rnap_monomers[id] = get_monomer_indexes([id])[0]

        # REPLISOMES:
        replisome_trimer_subunits = sim_data.molecule_groups.replisome_trimer_subunits
        replisome_monomer_subunits = sim_data.molecule_groups.replisome_monomer_subunits
        replisome_subunit_ids = replisome_trimer_subunits + replisome_monomer_subunits
        replisome_monomers = {}
        for subunit in replisome_subunit_ids:
            if subunit in protein_ids:
                replisome_monomers[subunit] = get_monomer_indexes([subunit])[0]
            elif subunit in complexation_complex_ids:
                monomer_ids = complexes_to_monomers[subunit]
                for id in monomer_ids:
                    if id in protein_ids:
                        replisome_monomers[id] = get_monomer_indexes([id])[0]

        # TF MONOMERS
        tfs = sim_data.process.transcription_regulation.tf_ids
        tf_to_tf_type = sim_data.process.transcription_regulation.tf_to_tf_type
        tf_subunit_ids = [tf_id + f'[{sim_data.getter.get_compartment(tf_id)[0]}]'
                          for tf_id in tfs]
        tf_monomers = {}
        # TF type(s) (0CS, 1CS or 2CS) of the TFs each monomer is part of:
        tf_monomer_types = {}
        for tf_id, subunit in zip(tfs, tf_subunit_ids):
            if subunit in protein_ids:
                monomer_ids = [subunit]
            elif subunit in complexation_complex_ids:
                monomer_ids = complexes_to_monomers[subunit]
            elif subunit in equilibrium_complex_ids:
                monomer_ids = eq_complexes_to_monomers[subunit]
            elif subunit in two_component_system_complex_ids:
                monomer_ids = tcs_complexes_to_monomers[subunit]
            else:
                monomer_ids = []
            for id in monomer_ids:
                if id in protein_ids:
                    tf_monomers[id] = get_monomer_indexes([id])[0]
                    tf_monomer_types.setdefault(id, set()).add(tf_to_tf_type[tf_id])

        # create self variables:
        self.complexation_monomers = complexation_monomers
        self.equilibrium_monomers = equilibrium_monomers
        self.two_component_system_monomers = two_component_system_monomers
        self.ribosome_monomers = ribosome_monomers
        self.rnap_monomers = rnap_monomers
        self.replisome_monomers = replisome_monomers
        self.tf_monomers = tf_monomers
        self.tf_monomer_types = tf_monomer_types


    def classify_molecule_types_multi(self, monomer_id, monomer_idx):
        """
        Classify a monomer into all categories it belongs to.
        """
        categories = []

        if monomer_id in self.ribosome_monomers and self.ribosome_monomers[monomer_id] == monomer_idx:
            categories.append('Ribosome')
        if monomer_id in self.rnap_monomers and self.rnap_monomers[monomer_id] == monomer_idx:
            categories.append('RNAP')
        if monomer_id in self.replisome_monomers and self.replisome_monomers[monomer_id] == monomer_idx:
            categories.append('Replisome')
        if monomer_id in self.tf_monomers and self.tf_monomers[monomer_id] == monomer_idx:
            tf_types = '/'.join(sorted(self.tf_monomer_types[monomer_id]))
            categories.append(f'{tf_types} transcription factor')
        if monomer_id in self.complexation_monomers and self.complexation_monomers[
            monomer_id] == monomer_idx:
            categories.append('Complexation complex')
        if monomer_id in self.equilibrium_monomers and self.equilibrium_monomers[
            monomer_id] == monomer_idx:
            categories.append('Equilibrium complex')
        if monomer_id in self.two_component_system_monomers and self.two_component_system_monomers[
            monomer_id] == monomer_idx:
            categories.append('two component system')

        if not categories:
            return 'Monomer only'
        elif len(categories) == 1:
            return categories[0] + ' subunit'
        else:
            # Create compound category name
            return ' + '.join(sorted(categories)) + ' subunit'


    def check_monomer_validity(self, sim_data, molecule_list):
        """
        Verify that molecules are valid monomers in the simulation.

        Args:
            molecule_list: List of molecule IDs to check
        Returns:
            tuple: (valid_monomers, monomer_type)
                - valid_monomers: List of valid monomer IDs
                - monomer_type: Dict mapping monomer IDs to their subunit types
        """
        # First, check validity of molecule IDs and get compartment tags:
        revised_molecule_list = self.check_validity_and_get_compartment(
                                    sim_data, molecule_list)

        # Next, check that the molecules are valid complexes and determine their types:
        all_monomers = sim_data.process.translation.monomer_data['id'].tolist()
        valid_monomers = []
        monomer_type = {}

        for molecule in revised_molecule_list:
            if molecule in all_monomers:
                valid_monomers.append(molecule)
                molecule_idx = all_monomers.index(molecule)
                monomer_type[molecule] = (
                    self.classify_molecule_types_multi(molecule, molecule_idx))
            else:
                print(f"{molecule} is not a valid monomer.")

        return valid_monomers, monomer_type

    def extract_doubling_times(self, cell_paths):
        """
        Extract time data and generation boundaries from simulation output.

        Args:
            cell_paths: List of paths to cell simulation directories
        Returns:
            tuple: (time, doubling_times, end_generation_times,
                    start_generation_indices, end_generation_indices)
                - time: Array of simulation timepoints
                - doubling_times: Duration of each generation
                - end_generation_times: Cumulative time at end of each generation
                - start_generation_indices: Index where each generation starts
                - end_generation_indices: Index where each generation ends
        """
        time = read_stacked_columns(cell_paths, 'Main', 'time').squeeze()

        # Generation boundaries come from each generation's own row count and
        # recorded timestamps, so they stay exact when timesteps are not 1 s
        # (the sim shortens a timestep whenever a process requires it):
        rows_per_generation = read_stacked_columns(
            cell_paths, 'Main', 'time',
            fun=lambda x: np.array([[x.shape[0]]])).reshape(-1).astype(int)
        start_generation_indices = np.concatenate(
            ([0], np.cumsum(rows_per_generation)[:-1]))
        end_generation_indices = start_generation_indices + rows_per_generation - 1
        end_generation_times = time[end_generation_indices]
        doubling_times = end_generation_times - time[start_generation_indices]

        return time, doubling_times, end_generation_times, start_generation_indices, end_generation_indices

    def find_downstream_complexes_with_paths(self, molecule_id):
        """
        Find all parent complexes that contain the molecule and trace assembly
        paths.

        Note: some complexes may have intermediate complexes between the
        base monomers and the final "grandparent" complexes.
        """
        parent_complexes = {}
        visited_paths = set()  # Path states already explored (to skip repeats)

        def get_parent_complexes(molecule_id):
            """
            Get all complexes that directly contain this molecule as a subunit.
            """
            parents = {}

            # Check complexation complexes (irreversible assembly)
            if molecule_id in self.c_m2pc_dict:
                for parent_id in self.c_m2pc_dict[molecule_id].keys():
                    stoich = self.c_m2pc_dict[molecule_id][parent_id][1]['stoichiometry']
                    reaction_id = self.c_m2pc_dict[molecule_id][parent_id][0]['reaction_id']
                    parents[parent_id] = ('complex', np.negative(stoich), reaction_id)

            # Check equilibrium complexes (reversible ligand binding)
            if molecule_id in self.eq_m2pc_dict:
                for parent_id in self.eq_m2pc_dict[molecule_id].keys():
                    stoich = self.eq_m2pc_dict[molecule_id][parent_id][1]['stoichiometry']
                    reaction_id = self.eq_m2pc_dict[molecule_id][parent_id][0]['reaction_id']
                    parents[parent_id] = ('eq_complex', np.negative(stoich), reaction_id)

            # Check two-component system complexes (phosphotransfer reactions)
            if molecule_id in self.tcs_m2pc_dict:
                for parent_id in self.tcs_m2pc_dict[molecule_id].keys():
                    stoich = self.tcs_m2pc_dict[molecule_id][parent_id][1]['stoichiometry']
                    reaction_id = self.tcs_m2pc_dict[molecule_id][parent_id][0]['reaction_id']
                    parents[parent_id] = ('two component system', np.negative(stoich), reaction_id)

            return parents

        def determine_molecule_type(mol_id):
            if mol_id in self.monomerIDs:
                return 'monomer'
            elif mol_id in self.complexIDs:
                return 'complex'
            elif mol_id in self.eqComplexIDs:
                return 'eq_complex'
            elif mol_id in self.tcsComplexIDs:
                return 'two component system'
            else:
                return 'unknown'

        def explore_parents(current_molecule, current_path, current_stoich, is_start=False):
            """
            Recursively explore parent complexes, building paths upward in the hierarchy.
            """
            mol_type = determine_molecule_type(current_molecule)

            # Skip a path state (molecule + path so far + stoich) already explored:
            path_signature = (current_molecule, tuple((m[0], m[1]) for m in current_path),
                              current_stoich)
            if path_signature in visited_paths:
                # Already explored this exact path state
                return
            visited_paths.add(path_signature)

            # Get all complexes that directly use this molecule as a subunit
            parents = get_parent_complexes(current_molecule)

            # BASE CASE: No parents found (the molecule is a top-level complex)
            if not parents:
                # Only save if this isn't the starting molecule
                if not is_start:
                    if current_molecule not in parent_complexes:
                        parent_complexes[current_molecule] = []

                    # Complete the path by adding the terminal complex
                    complete_path = current_path + [(current_molecule, 1, mol_type, None)]

                    # Only record a path that is not already recorded for this complex:
                    path_exists = False
                    for existing_entry in parent_complexes[current_molecule]:
                        if (existing_entry['path'] == complete_path and
                                existing_entry['total_stoich'] == current_stoich):
                            path_exists = True
                            break

                    if not path_exists:
                        parent_complexes[current_molecule].append({
                            'path': complete_path,
                            'total_stoich': current_stoich
                        })
                return

            # RECURSIVE CASE: Explore each parent of the molecule:
            for parent_id, (parent_type, stoich, reaction_id) in parents.items():
                if parent_id not in parent_complexes:
                    parent_complexes[parent_id] = []

                # Build complete path to this parent
                path_to_parent = current_path + [
                    (current_molecule, stoich, mol_type, reaction_id),
                    (parent_id, 1, parent_type, None)]

                # Only record a path that is not already recorded for this complex:
                path_exists = False
                for existing_entry in parent_complexes[parent_id]:
                    if (existing_entry['path'] == path_to_parent and
                            existing_entry['total_stoich'] == current_stoich * stoich):
                        path_exists = True
                        break

                if not path_exists:
                    # Save this parent with its complete path
                    parent_complexes[parent_id].append({
                        'path': path_to_parent,
                        'total_stoich': current_stoich * stoich})

                # Continue exploring to see if this parent is itself a subunit of other complexes:
                new_path = current_path + [(current_molecule, stoich, mol_type, reaction_id)]
                new_stoich = current_stoich * stoich

                # Recursively explore the parent:
                explore_parents(parent_id, new_path, new_stoich, is_start=False)

        # Verify the starting molecule is valid:
        mol_type = determine_molecule_type(molecule_id)

        if mol_type == 'unknown':
            print(f"Warning: {molecule_id} is not a valid molecule")
            return {}

        # Start exploration from the given molecule
        explore_parents(molecule_id, [], 1, is_start=True)

        return parent_complexes

    def format_parent_path_sentence(self, molecule_id, parent_complex_id, path_info):
        """
        Format a sentence describing how a molecule appears in a parent complex.

        Args:
            molecule_id: The starting molecule (monomer or complex)
            parent_complex_id: The parent complex
            path_info: Dictionary with 'path' and 'total_stoich'
        Returns:
            str: Formatted sentence
        """
        path = path_info['path']
        total_stoich = path_info['total_stoich']
        path_parts = []

        # Add all steps in the path with their stoichiometries:
        for mol_id, stoich, mol_type, reaction_id in path:
            path_parts.append(f"{abs(stoich):g} {mol_id}")

        # Join with line breaks every 2 arrows to prevent very long lines:
        path_string_parts = []
        for i in range(0, len(path_parts), 2):
            chunk = " → ".join(path_parts[i:i + 2])
            path_string_parts.append(chunk)
        path_string = " →\n".join(path_string_parts)

        sentence = (f"Monomer counts within {parent_complex_id}, "
                    f"\n{abs(total_stoich):g} per (via {path_string})")

        return sentence

    def analyze_molecule_parents(self, molecule_id):
        """
        Main function to analyze what parent complexes contain a molecule.

        Args:
            molecule_id: The molecule (monomer or complex) to analyze
        Returns:
            dict: Dictionary of parent complexes with their paths
        """
        print(f"\nAnalyzing parent complexes for: {molecule_id}\n")
        parent_complexes = self.find_downstream_complexes_with_paths(molecule_id)

        # Filter out cases where the molecule is its own parent (no actual parents):
        filtered_complexes = {k: v for k, v in parent_complexes.items() if k != molecule_id}

        descriptions = []
        for parent_id, path_list in filtered_complexes.items():
            for path_info in path_list:
                sentence = self.format_parent_path_sentence(
                    molecule_id, parent_id, path_info
                )
                descriptions.append(sentence)

        if not descriptions:
            print(f"No parent complexes found for {molecule_id}")
        else:
            print(f"Found {len(filtered_complexes)} parent complex(es):\n")
            for desc in descriptions:
                print(desc)
                print()

        return filtered_complexes, descriptions

    # Build the complexation parent and downstream molecule to complex dictionaries:
    def build_complexation_dictionaries(self, sim_data):
        """
        Build a dictionary mapping each molecule to the complexation complexes
        it is a direct subunit of.

        Args:
            sim_data: Simulation data object
        Returns:
            dict: {subunit_id: {complex_id: [
                  {'reaction_id': str},
                  {'stoichiometry': int}]}}  # Negative value (consumed)
        """
        molecule_names = sim_data.process.complexation.molecule_names
        subunit_names = sim_data.process.complexation.subunit_names
        rxn_ids = sim_data.process.complexation.ids_reactions
        sm = sim_data.process.complexation.stoich_matrix()

        molecules_to_parent_complexes_dict = {}
        for subunit in subunit_names:
            # Find the matrix index where this subunit is as a molecule:
            subunit_index = molecule_names.index(subunit)

            # Find the indices for reactions the subunit is a reactant in:
            reaction_indices = np.where(sm[subunit_index, :] < 0)[0]

            # For each reaction index, find the complex that is formed
            # (NOTE: each reaction forms one complexation complex):
            parent_complexes = {}
            for rxn_idx in reaction_indices:
                complex_index = np.where(sm[:, rxn_idx] > 0)[0]
                complex_name = molecule_names[complex_index[0]]
                # NOTE: this is the subunit's stoichiometry!
                parent_complexes[complex_name] = [
                    {'reaction_id': rxn_ids[rxn_idx]},
                    {'stoichiometry': sm[subunit_index, rxn_idx]}]

            molecules_to_parent_complexes_dict[subunit] = parent_complexes

        return molecules_to_parent_complexes_dict

    def build_equilibrium_dictionaries(self, sim_data):
        """
        Build a dictionary mapping each molecule to the equilibrium complexes
        it is a direct subunit of. Note: complexation complexes can be subunits
        of equilibrium complexes.

        Args:
            sim_data: Simulation data object
        Returns:
            dict: {subunit_id: {complex_id: [
                  {'reaction_id': str},
                  {'stoichiometry': int}]}}  # Negative value (consumed)
        """
        molecule_names = sim_data.process.equilibrium.molecule_names
        rxn_ids = sim_data.process.equilibrium.rxn_ids
        sm = sim_data.process.equilibrium.stoich_matrix()

        molecules_to_parent_complexes_dict = {}
        for subunit_index, subunit in enumerate(molecule_names):
            # Find the indices for reactions the subunit is a reactant in:
            reaction_indices = np.where(sm[subunit_index, :] < 0)[0]

            if len(reaction_indices) == 0:
                # skip molecules that do not form any complexes
                continue

            # For each reaction index, find the complex that is formed (each
            # reaction forms one equilibrium complex):
            parent_complexes = {}
            for rxn_idx in reaction_indices:
                complex_index = np.where(sm[:, rxn_idx] > 0)[0]
                complex_name = molecule_names[complex_index[0]]
                # NOTE: this is the subunit's stoichiometry!
                parent_complexes[complex_name] = [
                    {'reaction_id': rxn_ids[rxn_idx]},
                    {'stoichiometry': sm[subunit_index, rxn_idx]}]

            molecules_to_parent_complexes_dict[subunit] = parent_complexes

        return molecules_to_parent_complexes_dict

    def get_tcs_small_molecules(self, sim_data):
        """
        Get the small molecules shared by all two component system (TCS)
        reactions (ATP, ADP, Pi, water and H+), i.e. every TCS molecule whose
        sim_data molecule type is not one of the TCS protein types (response
        regulators, histidine kinases and ligand-bound histidine kinases, plus
        their phosphorylated forms).

        Args:
            sim_data: Simulation data object
        Returns:
            set: IDs of the TCS small molecules
        """
        tcs_protein_types = {"RR", "PHOSPHO-RR", "HK", "PHOSPHO-HK",
            "HK-LIGAND", "PHOSPHO-HK-LIGAND"}
        tcs = sim_data.process.two_component_system
        return {name for name, molecule_type
            in zip(tcs.molecule_names, tcs.molecule_types)
            if molecule_type not in tcs_protein_types}

    def build_tcs_dictionaries(self, sim_data):
        """
        Build dictionaries mapping each unphosphorylated molecule to the
        phosphorylated two component system (TCS) complexes it forms. Note:
        complexation and equilibrium complexes can be subunits of TCS complexes.

        For positively oriented TCS systems, a subunit can form the same TCS
        complex via BOTH a ligand-bound HK phosphotransfer reaction AND a
        ligand-free HK phosphotransfer reaction, so the two reaction types are
        kept in separate dictionaries (one dictionary keyed by complex would
        keep only one of the two reactions).

        Args:
            sim_data: Simulation data object
        Returns:
            tuple: (molecules_to_parent_complexes_via_lbptrs,
                    molecules_to_parent_complexes_via_ptrs)
            - each is {subunit_id: {complex_id: [
                  {'reaction_id': str},
                  {'stoichiometry': int}]}}  # Negative value (consumed)
              for ligand-bound (lbptrs) and all other (ptrs) reactions
        """
        molecule_names = list(sim_data.process.two_component_system.molecule_names)
        complex_ids = list(sim_data.process.two_component_system.complex_to_monomer.keys())
        molecule_types = sim_data.process.two_component_system.molecule_types
        rxn_ids = sim_data.process.two_component_system.rxn_ids
        sm = sim_data.process.two_component_system.stoich_matrix()
        tcs_small_molecules = self.get_tcs_small_molecules(sim_data)
        pairings = {"RR": "PHOSPHO-RR", "PHOSPHO-RR": "RR", "HK": "PHOSPHO-HK",
                    "PHOSPHO-HK": "HK", "HK-LIGAND": "PHOSPHO-HK-LIGAND",
                    "PHOSPHO-HK-LIGAND": "HK-LIGAND"}

        molecules_to_parent_complexes_via_lbptrs = {}
        molecules_to_parent_complexes_via_ptrs = {}
        for subunit_index, subunit in enumerate(molecule_names):
            # If the molecule is a TCS complex or not a complex or monomer, skip it:
            if subunit in tcs_small_molecules or subunit in complex_ids:
                continue
            complement = pairings[molecule_types[subunit_index]]

            # Find the indices for reactions the subunit is a reactant in:
            reaction_indices = np.where(sm[subunit_index, :] < 0)[0]
            if len(reaction_indices) == 0:
                continue

            parent_complexes_lbptrs = {}
            parent_complexes_ptrs = {}
            for rxn_idx in reaction_indices:
                # Determine if this is a ligand-bound or ligand-free HK phosphotransfer reaction:
                rxn_id = rxn_ids[rxn_idx]
                if "BOUND" in rxn_id:
                    parent_dict_to_use = parent_complexes_lbptrs
                else:
                    parent_dict_to_use = parent_complexes_ptrs

                # TCS reactions can have multiple products, so keep only the
                # complement of this reactant subunit:
                for product_index in np.where(sm[:, rxn_idx] > 0)[0]:
                    if molecule_types[product_index] == complement:
                        # NOTE: this is the subunit's stoichiometry:
                        parent_dict_to_use[molecule_names[product_index]] = [
                            {'reaction_id': rxn_id},
                            {'stoichiometry': sm[subunit_index, rxn_idx]}]

            molecules_to_parent_complexes_via_lbptrs[subunit] = parent_complexes_lbptrs
            molecules_to_parent_complexes_via_ptrs[subunit] = parent_complexes_ptrs

        return (molecules_to_parent_complexes_via_lbptrs,
                molecules_to_parent_complexes_via_ptrs)

    def build_tcs_reaction_event_estimates(self, sim_data, cell_paths):
        """
        Estimate reaction fluxes from TCS molecule changes using non-negative
        least squares.

        NOTE: there could technically be many possible reaction flux solutions
        that produce a fit similar to molecule_changes, but non-negative least
        squares (NNLS) tends to give a sparse solution that uses few reactions
        to hit the molecule changes. This is a reasonable assumption to make
        given the small number of reactions and given that molecule_changes are
        already weighted by the rates of reaction.

        If this does not appear to be working well because the reaction rates
        change or something, it might be best to try a different method.

        Args:
            sim_data: Simulation data object
            cell_paths: paths to the cell simulations whose
                TwoComponentSystems moleculeChanges (the per-step molecule
                count changes from the ODE integration) are read
        Returns:
            tuple: ((n_timesteps, n_reactions) estimated reaction counts,
                (n_timesteps, n_molecules) bool mask of molecules whose change
                the estimate misses by more than 0.5)
        """
        # Extract the molecule changes for each TCS molecule across the cell paths:
        all_molecule_changes = read_stacked_columns(
            cell_paths, 'TwoComponentSystems', 'moleculeChanges')
        sm = sim_data.process.two_component_system.stoich_matrix()
        molecule_names = sim_data.process.two_component_system.molecule_names
        tcs_small_molecules = self.get_tcs_small_molecules(sim_data)
        n_timesteps = all_molecule_changes.shape[0]
        n_reactions = sm.shape[1]

        # Estimate the reaction fluxes for each step using non-negative least
        # squares on the protein rows only. The emitted ATP/ADP/Pi/H2O/H+
        # changes are not the sum of the individual reactions because the TCS
        # process sets ATP use to the net phosphorylation summed over all
        # systems and never records Pi, so whenever a protein is
        # dephosphorylated they understate what the reactions used. Matching
        # them would need negative reaction counts, which NNLS can't give, so
        # including those rows would pull the estimates too low. The fit is
        # still checked per protein:
        protein_rows = np.array(
            [name not in tcs_small_molecules for name in molecule_names])
        sm_proteins = sm[protein_rows]
        estimated_reaction_fluxes = np.zeros((n_timesteps, n_reactions))
        for i in range(n_timesteps):
            reaction_fluxes, _ = scipy.optimize.nnls(
                sm_proteins, all_molecule_changes[i, protein_rows])
            estimated_reaction_fluxes[i] = np.round(reaction_fluxes)
        misfit = np.abs(
            estimated_reaction_fluxes @ sm.T - all_molecule_changes) > 0.5
        misfit[:, ~protein_rows] = False
        n_bad = misfit.any(axis=1).sum()
        if n_bad:
            print(f"TCS reaction estimate misfits some molecule at "
                  f"{n_bad}/{n_timesteps} timesteps (per-reaction counts in legend)")

        return estimated_reaction_fluxes, misfit

    def find_extra_tcs_reactions(self, molecule, complex, rxns_plotted):
        """
        Some TCS complexes can be formed from the same reactant by two
        different phosphotransfer reactions (ligand-free and ligand-bound).
        This function checks whether the molecule also forms the parent TCS
        complex via the other reaction type, and if so, returns that reaction
        so its events are plotted too (if they have not already been plotted).

        This matters for the events panel, which shows every reaction that
        consumes the free monomer: if one of the two reactions were skipped,
        the monomer's TCS events would be undercounted.

        Args:
            molecule: ID of the monomer (the reactant)
            complex: ID of the parent TCS complex the molecule forms
            rxns_plotted: reaction IDs already plotted
        Returns:
            list: IDs of the other reactions forming the parent TCS complex
            from the molecule that have not been plotted yet
        """
        reaction_IDs = []

        # Check if the reaction events for the TCS complex via the ligand-free
        # phosphotransfer reactions have already been plotted:
        if molecule in self.tcs_m2pc_via_ptrs_dict.keys():
            alt_parent_complex_info = self.tcs_m2pc_via_ptrs_dict[molecule]
            if complex in alt_parent_complex_info.keys():
                alt_rxn_id = alt_parent_complex_info[complex][0]['reaction_id']
                if alt_rxn_id not in rxns_plotted:
                    reaction_IDs.append(alt_rxn_id)

        # Check if the reaction events for the TCS complex via the ligand-bound
        # phosphotransfer reactions have already been plotted:
        if molecule in self.tcs_m2pc_via_lbptrs_dict.keys():
            alt_parent_complex_info = self.tcs_m2pc_via_lbptrs_dict[molecule]
            if complex in alt_parent_complex_info.keys():
                alt_rxn_id = alt_parent_complex_info[complex][0]['reaction_id']
                if alt_rxn_id not in rxns_plotted:
                    reaction_IDs.append(alt_rxn_id)

        return reaction_IDs

    def is_a_TF(self, sim_data, molecule):
        """
        Determine whether a molecule is a TF.
        Args:
            sim_data: simulation data object
            molecule: molecule ID to check
        Returns:
            list: TF ID that match the molecule (if applicable),
            but returns an empty list if no match is found.

        """
        # Get the TF IDs:
        tfs = sim_data.process.transcription_regulation.tf_ids
        tf_ids = [tf_id + f'[{sim_data.getter.get_compartment(tf_id)[0]}]'
                          for tf_id in tfs]

        # Determine if the molecule is a TF:
        TF_id_match = [tf_id for tf_id in tf_ids if molecule == tf_id]

        if TF_id_match != []:
            print(f"{molecule} is a TF.")

        return TF_id_match

    def tf_index(self, sim_data, molecule):
        """
        Get the index of a TF in the list of TFs.
        Args:
            sim_data: simulation data object
            molecule: molecule ID to check

        Returns:
            int: index of the TF in the list of TFs
        """
        # Get the TF IDs:
        tfs = sim_data.process.transcription_regulation.tf_ids
        tf_ids = [tf_id + f'[{sim_data.getter.get_compartment(tf_id)[0]}]'
                          for tf_id in tfs]

        # Find the index of the TF in the TF list:
        tf_index = tf_ids.index(molecule)

        return tf_index

    def is_a_ribosome_subunit(self, sim_data, molecule):
        """
        Determine whether a molecule is a ribosome subunit.
        Args:
            sim_data: simulation data object
            molecule: molecule ID to check
        Returns:
            list: ribosome subunit ID that matches the molecule (if applicable)
        """
        # Get the 30s ribosome subunit IDs:
        ribosome30s_subunits = sim_data.process.complexation.get_monomers(
            sim_data.molecule_ids.s30_full_complex)
        ribosome30s_subunit_ids = ribosome30s_subunits["subunitIds"].tolist()

        # Get the 50s ribosome subunit IDs:
        ribosome50s_subunits = sim_data.process.complexation.get_monomers(
            sim_data.molecule_ids.s50_full_complex)
        ribosome50s_subunit_ids = ribosome50s_subunits["subunitIds"].tolist()

        ribosome_subunits = (ribosome30s_subunit_ids + ribosome50s_subunit_ids
                             + [sim_data.molecule_ids.s30_full_complex]
                             + [sim_data.molecule_ids.s50_full_complex])

        # Determine if the molecule is a ribosome subunit:
        ribosome_subunit_id_match = [subunit_id for subunit_id in
                                     ribosome_subunits if molecule == subunit_id]

        if ribosome_subunit_id_match != []:
            print(f"{molecule} is a ribosome subunit.")

        return ribosome_subunit_id_match

    def is_an_RNAP_subunit(self, sim_data, molecule):
        """
        Determine whether a molecule is an RNAP subunit.
        Args:
            sim_data: simulation data object
            molecule: molecule ID to check
        Returns:
            list: RNAP subunit ID that match the molecule (if applicable)
        """
        # Get the RNAP subunit IDs:
        rnap_subunits = sim_data.process.complexation.get_monomers(
            sim_data.molecule_ids.full_RNAP)
        rnap_subunit_ids = rnap_subunits["subunitIds"].tolist()

        # Determine if the molecule is an RNAP subunit:
        rnap_subunit_id_match = [subunit_id for subunit_id in
                                 rnap_subunit_ids if molecule == subunit_id]

        if rnap_subunit_id_match != []:
            print(f"{molecule} is an RNAP subunit.")

        # Check if the molecule is the inactive RNAP complex:
        if molecule == "APORNAP-CPLX[c]":
            rnap_subunit_id_match += ["APORNAP-CPLX[c]"]
            print(f"{molecule} is the inactive RNAP complex.")

        return rnap_subunit_id_match

    def is_a_replisome_subunit(self, sim_data, molecule):
        """
        Determine whether a molecule is a replisome subunit.
        Args:
            sim_data: simulation data object
            molecule: molecule ID to check
        Returns:
            list: replisome subunit ID that match the molecule (if applicable)
        """
        # Get the replisome subunit IDs:
        replisome_trimer_subunits = sim_data.molecule_groups.replisome_trimer_subunits
        replisome_monomer_subunits = sim_data.molecule_groups.replisome_monomer_subunits
        replisome_subunit_ids = replisome_trimer_subunits + replisome_monomer_subunits

        # Determine if the molecule is a replisome subunit:
        replisome_subunit_id_match = [subunit_id for subunit_id in
                replisome_subunit_ids if molecule == subunit_id]

        if replisome_subunit_id_match != []:
            print(f"{molecule} is a replisome subunit.")

        return replisome_subunit_id_match

    def replisome_subunit_stoich(self, sim_data, molecule):
        """
        Get the number of copies of a replisome subunit in one replisome.
        Args:
            sim_data: simulation data object
            molecule: replisome subunit molecule ID
        Returns:
            float: copies of the subunit per active replisome
        """
        replisome_trimer_subunits = sim_data.molecule_groups.replisome_trimer_subunits
        replisome_monomer_subunits = sim_data.molecule_groups.replisome_monomer_subunits
        replisome_subunit_ids = replisome_trimer_subunits + replisome_monomer_subunits
        replisome_stoich = np.hstack(
            (3 * np.ones(len(replisome_trimer_subunits)),
             np.ones(len(replisome_monomer_subunits))))
        return replisome_stoich[replisome_subunit_ids.index(molecule)]

    def do_plot(self, seedOutDir, plotOutDir, plotOutFileName, simDataFile,
                validationDataFile, metadata):
        with open(simDataFile, 'rb') as f:
            sim_data = pickle.load(f)
        # Only the first-daughter lineage (daughter cell 000000 of every
        # generation) is plotted. With SINGLE_DAUGHTERS=0, each generation also
        # has sibling cells, and stacking those would join unrelated cells into
        # one trajectory with time jumping backwards:
        all_cell_paths = self.ap.get_cells()
        cell_paths = np.array([
            cell_path for cell_path in all_cell_paths
            if int(os.path.basename(os.path.normpath(cell_path))) == 0])
        if len(cell_paths) < len(all_cell_paths):
            print(f"NOTE: plotting only the first-daughter lineage "
                  f"({len(cell_paths)} of {len(all_cell_paths)} cells), so "
                  f"sibling cells are left out.")

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
        tcs_small_molecules = self.get_tcs_small_molecules(sim_data)

        # Generate subunit categories and classify all proteins:
        self.generate_subunit_categories(sim_data)
        self.monomer_classifications = {}
        for i, protein_id in enumerate(self.monomerIDs):
            self.monomer_classifications[protein_id] = (
                self.classify_molecule_types_multi(protein_id, i))

        # Extract the IDs of the three types of complexes:
        self.complexIDs = sim_data.process.complexation.ids_complexes
        self.eqComplexIDs = sim_data.process.equilibrium.ids_complexes
        self.tcsComplexIDs = list(
            sim_data.process.two_component_system.complex_to_monomer.keys())

        # Generate the molecule to parent complexes dictionaries:
        self.c_m2pc_dict = self.build_complexation_dictionaries(sim_data)
        self.eq_m2pc_dict = self.build_equilibrium_dictionaries(sim_data)
        self.tcs_m2pc_via_lbptrs_dict, self.tcs_m2pc_via_ptrs_dict = (
            self.build_tcs_dictionaries(sim_data))

        # Define a TCS molecules to parent complexes dictionary:
        self.tcs_m2pc_dict = {}
        for molecule in set(list(self.tcs_m2pc_via_lbptrs_dict.keys()) + list(
                self.tcs_m2pc_via_ptrs_dict.keys())):
            self.tcs_m2pc_dict[molecule] = {}
            if molecule in self.tcs_m2pc_via_lbptrs_dict:
                self.tcs_m2pc_dict[molecule].update(
                    self.tcs_m2pc_via_lbptrs_dict[molecule])
            if molecule in self.tcs_m2pc_via_ptrs_dict:
                self.tcs_m2pc_dict[molecule].update(self.tcs_m2pc_via_ptrs_dict[molecule])

        # Check that all the valid inputs are indeed valid monomers:
        valid_monomers, _ = self.check_monomer_validity(sim_data, PLOT_MONOMERS)
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
        (time, _, end_generation_times, _, _) = self.extract_doubling_times(cell_paths)
        cell_lengths = np.atleast_1d(read_stacked_columns(
            cell_paths, 'Main', 'time', fun=lambda x: len(x)).squeeze())
        generation_start_idxs = np.r_[0, np.cumsum(cell_lengths)[:-1]].astype(int)

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
                listener_cache['tcs'] = self.build_tcs_reaction_event_estimates(
                    sim_data, cell_paths)
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
            downstream_complexes, _ = self.analyze_molecule_parents(monomer_id)
            complexes_by_id = {}
            for complex_id, path_list in downstream_complexes.items():
                unique_paths = []
                for path_data in path_list:
                    if any(path_data['path'] == up['path_data']['path']
                           and path_data['total_stoich'] == up['path_data']['total_stoich']
                           for up in unique_paths):
                        continue
                    unique_paths.append({
                        'path_data': path_data,
                        'path_string': self.format_parent_path_sentence(
                            monomer_id, complex_id, path_data)})
                complexes_by_id[complex_id] = unique_paths
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
                    if self.is_a_TF(sim_data, complex_id) != []:
                        tf_idx = self.tf_index(sim_data, complex_id)
                        contrib = tf_bound[:, tf_idx] * total_stoich
                        total_complexed_counts += contrib
                        ax1.plot(time, contrib, color='lightcoral',
                                 label=f'Monomer in bound {complex_id} TFs\n'
                                       f'({total_stoich:g} per)',
                                 linewidth=0.75, alpha=0.75, linestyle="-.")
                        already_counted_as_tf = True

                    if self.is_a_ribosome_subunit(sim_data, complex_id) != []:
                        contrib = active_ribosome * total_stoich
                        total_complexed_counts += contrib
                        ax1.plot(time, contrib, color='gold',
                                 label=f'Monomer in ribosomes via {complex_id}\n'
                                       f'({total_stoich:g} per)',
                                 linewidth=0.75, alpha=0.75, linestyle="-.")

                    if self.is_an_RNAP_subunit(sim_data, complex_id) != []:
                        contrib = active_rnap * total_stoich
                        total_complexed_counts += contrib
                        ax1.plot(time, contrib, color='orange',
                                 label=f'Monomer in active RNAPs via {complex_id}\n'
                                       f'({total_stoich:g} per)',
                                 linewidth=0.75, alpha=0.75, linestyle="-.")

                    if self.is_a_replisome_subunit(sim_data, complex_id) != []:
                        rep_stoich = self.replisome_subunit_stoich(sim_data, complex_id)
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
            if (self.is_a_replisome_subunit(sim_data, monomer_id) != []
                    and not already_counted_in_replisomes):
                rep_stoich = self.replisome_subunit_stoich(sim_data, monomer_id)
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
            if self.is_a_TF(sim_data, monomer_id) != [] and not already_counted_as_tf:
                contrib = tf_bound[:, self.tf_index(sim_data, monomer_id)]
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
                        for extra_rxn in self.find_extra_tcs_reactions(
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
            if monomer_id in self.eq_m2pc_dict:
                for eq_complex_id, info in self.eq_m2pc_dict[monomer_id].items():
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
            if (self.is_a_replisome_subunit(sim_data, monomer_id) != []
                    and mechanistic_replisome and not any(
                        self.is_a_replisome_subunit(sim_data, c) != []
                        for c, _ in unique_complex_path_pairs)):
                rep_stoich = self.replisome_subunit_stoich(sim_data, monomer_id)
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
            if self.is_a_TF(sim_data, monomer_id) != [] and not any(
                    self.is_a_TF(sim_data, c) != [] for c, _ in unique_complex_path_pairs):
                tf_idx = self.tf_index(sim_data, monomer_id)
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

            # Order the events by activity so the most active are listed first in
            # the legend and drawn underneath the rest, and the least active
            # are drawn on top so sparse events stay visible. Each label gets
            # the line's average |events| per timestep (excluding the initial
            # steps that hit the truncation limit below):
            event_lines = [line for line in ax2.get_lines()
                           if len(line.get_ydata()) == len(time)]
            # Rows that may be cut off/truncated in the first steps of each generation:
            initial_rows = np.unique(np.minimum(
                generation_start_idxs[:, None] + np.arange(N_INITIAL_STEPS_TO_CLIP),
                len(time) - 1))
            average_mask = np.ones(len(time), bool)
            average_mask[initial_rows] = False
            event_activity = {}
            for line in event_lines:
                y = np.abs(np.asarray(line.get_ydata(), float))[average_mask]
                activity = y.mean() if y.size else 0.
                event_activity[line] = activity
                average_str = f'{activity:.0f}' if activity >= 100 else f'{activity:.2g}'
                first_line, newline, rest = line.get_label().partition('\n')
                if rest.startswith('(estimated, '):
                    rest = rest.replace(
                        '(estimated, ', f'(estimated avg. {average_str}/step, ', 1)
                    line.set_label(f'{first_line}{newline}{rest}')
                elif rest.startswith('('):
                    rest = rest.replace('(', f'(avg. {average_str}/step, ', 1)
                    line.set_label(f'{first_line}{newline}{rest}')
                else:
                    line.set_label(f'{line.get_label()} (avg. {average_str}/step)')
            ordered_event_lines = sorted(event_lines, key=lambda line: -event_activity[line])
            for rank, line in enumerate(ordered_event_lines):
                line.set_zorder(2 + rank / len(ordered_event_lines))
                # plot using thin lines so dense event traces overlap less and
                # keep their own color instead of blending into each other:
                line.set_linewidth(0.4)
                line.set_alpha(0.85)

            # The first steps of each generation can spike well above the rest
            # of the panel. The first generation starts with no equilibrium or
            # TCS complexes (initialization only forms complexation complexes),
            # so those reactions form nearly all of them at once in its first
            # steps. In later generations, division splits each molecule
            # between the daughters independently (binomially), so complexes
            # and their free subunits start slightly off steady state and the
            # reactions correct this over the first few steps. To keep these
            # spikes from dwarfing the rest of the panel, if any value in the
            # first N_INITIAL_STEPS_TO_CLIP steps of a generation is more than
            # 1.5x beyond the y range of all other steps (plus 10% padding),
            # the y axis is limited to that range and each value outside it is
            # labeled at the axis edge in its line's color:
            event_lines = [line for line in ax2.get_lines()
                           if len(line.get_ydata()) == len(time)]
            # Clipping needs both initial rows and rows outside them (none when
            # N_INITIAL_STEPS_TO_CLIP is 0 or the generations are that short):
            if event_lines and 0 < len(initial_rows) < len(time):
                rest_mask = np.ones(len(time), bool)
                rest_mask[initial_rows] = False
                all_rest = np.concatenate(
                    [np.asarray(line.get_ydata(), float)[rest_mask]
                     for line in event_lines])
                y_low, y_high = min(all_rest.min(), 0), max(all_rest.max(), 0)
                pad = 0.1 * max(y_high - y_low, 1)
                y_low, y_high = y_low - pad, y_high + pad
                # Largest and smallest event of each line at the start of
                # each generation:
                peaks = []
                for line in event_lines:
                    y_all = np.asarray(line.get_ydata(), float)
                    for gen_idx, start_idx in enumerate(generation_start_idxs):
                        rows = initial_rows[(initial_rows >= start_idx)
                                            & (initial_rows < start_idx + N_INITIAL_STEPS_TO_CLIP)]
                        for i in {rows[np.argmax(y_all[rows])], rows[np.argmin(y_all[rows])]}:
                            peaks.append((y_all[i], time[i], line.get_color(), gen_idx))
                # Only truncate if some value dwarfs the rest of the panel, then
                # label every value that falls outside the new y range:
                if any(value < 1.5 * y_low or value > 1.5 * y_high
                       for value, _, _, _ in peaks):
                    ax2.set_ylim(y_low, y_high)
                    n_stacked = {}
                    for value, x, color, gen_idx in sorted(peaks, key=lambda c: -abs(c[0])):
                        if y_low <= value <= y_high:
                            continue
                        below = value < y_low
                        # Stack labels from the same generation start:
                        n = n_stacked.get((gen_idx, below), 0)
                        n_stacked[(gen_idx, below)] = n + 1
                        ax2.annotate(f'{value:.0f}', xy=(x, y_low if below else y_high),
                                     xytext=(4, (3 + 7 * n) * (1 if below else -1)),
                                     textcoords='offset points', color=color, fontsize=6,
                                     ha='left', va='bottom' if below else 'top')

            # No line is plotted along 0 where an event does not happen so the
            # plot doesn't get too crowded. Every event spike is still drawn up
            # from 0 though:
            for line in ordered_event_lines:
                y = np.asarray(line.get_ydata(), float).copy()
                keep = y != 0
                keep[1:] |= keep[:-1].copy()
                keep[:-1] |= (y[1:] != 0)
                y[~keep] = np.nan
                line.set_ydata(y)

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
