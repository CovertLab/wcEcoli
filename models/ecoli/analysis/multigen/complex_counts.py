"""
This plot visualizes the counts of a specified complex over time, along with:
- Free monomer subunit availability
- Complexation events that produce or consume the complex
- Counts within parent complexes and unique molecules (ribosomes, RNAPs, etc.)

Supported complex types:
- Complexation complexes (irreversible assembly)
- Equilibrium complexes (reversible ligand binding)
- Two-component system complexes (phosphotransfer/phosporylation/dephosporylation reactions)

In total, there are 1131 unique complexation complexes, 39 unique equilibrium
complexes, and 18 unique two component system complexes.
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

# List complexes you would like to plot here (with or without compartment tags):
PLOT_COMPLEXES = ["APORNAP-CPLX"]

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
                print(f"{molecule} is not a valid molecule in the simulation.")

        return revised_molecule_list

    def check_complex_validity(self, molecule_list):
        """
        Verify that molecules are valid complexes and determine their types.

        Args:
            molecule_list: List of molecule IDs to check
        Returns:
            tuple: (valid_complexes, complex_types)
                - valid_complexes: List of valid complex IDs
                - complex_types: Dict mapping complex IDs to their types
        """
        all_complexes = self.complexIDs + self.eqComplexIDs + self.tcsComplexIDs
        valid_complexes = []
        complex_types = {}

        for molecule in molecule_list:
            if molecule in all_complexes:
                valid_complexes.append(molecule)
                if molecule in self.complexIDs:
                    complex_types[molecule] = "complexation"
                elif molecule in self.eqComplexIDs:
                    complex_types[molecule] = "equilibrium"
                else:
                    complex_types[molecule] = "two component system"
            else:
                print(f"{molecule} is not a valid complex. Only complexation, "
                      f"equilibrium, and two component system complexes can be plotted.")

        return valid_complexes, complex_types

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

        doubling_times = read_stacked_columns(
            cell_paths, 'Main', 'time',
            fun=lambda x: (x[-1] - x[0])).squeeze().astype(int)

        end_generation_times = np.cumsum(doubling_times) + time[0]
        start_generation_indices = np.searchsorted(
            time, end_generation_times[:-1], side='left').astype(int)
        start_generation_indices = np.insert(start_generation_indices, 0, 0) + np.arange(
            len(doubling_times))
        end_generation_indices = start_generation_indices + doubling_times

        return time, doubling_times, end_generation_times, start_generation_indices, end_generation_indices

    def determine_complex_makeup(self, complex, complex_type):
        # Determine the makeup of the complex based on its type:
        if complex_type == 'complexation':
            makeup = self.c2ct[complex]
        elif complex_type == 'equilibrium':
            makeup = self.eq_c2ct[complex]
        elif complex_type == 'two component system':
            makeup = self.tcs_c2ct[complex]
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
                m2pc = self.c_m2pc_dict  # Complexation process
            elif molecule_type == 'eq_complex':
                m2pc = self.eq_m2pc_dict  # Equilibrium process
            elif molecule_type == 'two component system':
                m2pc = self.tcs_m2pc_dict  # Two-component system process
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

        def determine_molecule_type(molecule_id):
            """
            Classify a molecule as monomer, complex, equilibrium complex,
            or two-component system.

            Args:
                molecule_id: The molecule to classify
            Returns:
                str: 'monomer', 'complex', 'eq_complex', 'two component system', or 'unknown'
            """
            if molecule_id in self.monomerIDs:
                return 'monomer'
            elif molecule_id in self.complexIDs:
                return 'complex'
            elif molecule_id in self.eqComplexIDs:
                return 'eq_complex'
            elif molecule_id in self.tcsComplexIDs:
                return 'two component system'
            else:
                return 'unknown'

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
            mol_type = determine_molecule_type(molecule_id)

            # Skip molecules that aren't in our database
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
        complex_type = determine_molecule_type(complex_id)

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
            path_parts.append(f"{stoich} {molecule_id}")

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
                    f"{total_stoich} per\n(via {path_string})")

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

    def find_downstream_complexes_with_paths(self, molecule_id):
        """
        Find all parent complexes that contain the molecule and trace assembly
        paths. Note: some complexes may have intermediate complexes between the
        base monomers and the final "grandparent" complex.

        Args:
            molecule_id: The molecule (monomer or complex) to search for
        Returns:
            dict: {parent_complex_id: [{
                'path': [(molecule_id, stoich, type, reaction_id), ...],
                        # Ordered from starting molecule → intermediate → parent complex
                'total_stoich': int  # Total copies of starting molecule per parent complex
            }]}
        """
        parent_complexes = {}

        def get_parent_complexes(molecule_id):
            """
            Get all complexes that directly contain this molecule as a subunit.

            Args:
                molecule_id: The molecule to find parents for
            Returns:
                dict: {parent_id: (parent_type, stoichiometry, reaction_id)}
                    Where stoichiometry is positive (representing consumption of molecule_id)
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

            Args:
                current_molecule: Current molecule being explored
                current_path: Path built so far, as list of (molecule_id, stoich, type, reaction_id)
                             Ordered from starting molecule → current molecule
                current_stoich: Cumulative stoichiometry - how many copies of the starting
                               molecule are needed per current molecule
                is_start: True if this is the starting molecule (to avoid saving it as its own parent)
            """
            mol_type = determine_molecule_type(current_molecule)

            # Skip molecules not in our database
            if mol_type == 'unknown':
                print(f"  → Skipping {current_molecule} (not a monomer or complex)")
                return

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

                # Save this parent with its complete path
                parent_complexes[parent_id].append({
                    'path': path_to_parent,
                    'total_stoich': current_stoich * stoich})

                # Continue exploring to see if this parent is itself a subunit
                # of other complexes:
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
            path_parts.append(f"{abs(stoich)} {mol_id}")

        # Join with line breaks every 2 arrows to prevent very long lines:
        path_string_parts = []
        for i in range(0, len(path_parts), 2):
            chunk = " → ".join(path_parts[i:i + 2])
            path_string_parts.append(chunk)
        path_string = " →\n".join(path_string_parts)

        sentence = (f"{molecule_id} counts within {parent_complex_id}, "
                    f"\n{abs(total_stoich)} per (via {path_string})")

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
        Build dictionaries mapping molecules to complexation complexes they form.

        Args:
            sim_data: Simulation data object
        Returns:
            tuple: (molecules_to_parent_complexes_dict,
                    molecules_to_all_downstream_complexes_dict,
                    complex_to_complex_type)

            - molecules_to_parent_complexes_dict: {subunit_id: {complex_id: [
                  {'reaction_id': str},
                  {'stoichiometry': int},  # Negative value (consumed)
                  {'stoich_unknown': bool},
                  {'complex_type': 'homogeneous' or 'heterogeneous'}]}}
            - molecules_to_all_downstream_complexes_dict: Same structure but
            maps from complex to base monomers (bypassing intermediate complexes)
            - complex_to_complex_type: {complex_id: 'homogeneous' or 'heterogeneous'}
        """
        molecule_names = sim_data.process.complexation.molecule_names
        subunit_names = sim_data.process.complexation.subunit_names
        rxn_ids = sim_data.process.complexation.ids_reactions
        rsu = sim_data.process.complexation.reaction_stoichiometry_unknown
        sm = sim_data.process.complexation.stoich_matrix()
        smm = sim_data.process.complexation.stoich_matrix_monomers()

        molecules_to_parent_complexes_dict = {}
        molecules_to_all_downstream_complexes_dict = {}
        complex_to_complex_type = {}

        # Generate dictionary mapping molecules to direct parent complexes they form:
        for subunit in subunit_names:
            # Find the matrix index where this subunit is as a molecule:
            subunit_index = molecule_names.index(subunit)

            # Find the indices for reactions the subunit is a reactant in:
            reaction_indices = np.where(sm[subunit_index,:] < 0)[0]

            # For each reaction index, find the complex(es) that is(are) formed:
            parent_complexes = {}
            for rxn_idx in reaction_indices:
                # Initialize variables to hold info about the complex:
                complex_information = []
                stoich = {} # NOTE: this will be the subunit's stoichiometry!
                stoich_known = {}
                complex_type = {}
                reaction_name = {}

                # Find the complex formed in this reaction
                # (NOTE: each reaction forms one complexation complex):
                complex_index = np.where(sm[:, rxn_idx] > 0)[0]
                complex_name = molecule_names[complex_index[0]]
                cplx_idx = self.complexIDs.index(complex_name)

                # Find the number of unique subunits in this complex reaction:
                num_unique_subunits = len(np.where(smm[:, cplx_idx] < 0)[0])
                if num_unique_subunits > 1:
                    cplx_type = 'heterogeneous'
                else:
                    cplx_type = 'homogeneous'
                # Save the complex type in a separate dict for easy access later:
                complex_to_complex_type[complex_name] = cplx_type

                # Add complex information to lists
                reaction_name['reaction_id'] = rxn_ids[rxn_idx]
                stoich['stoichiometry'] = sm[subunit_index, rxn_idx]
                stoich_known['stoich_unknown'] = rsu[rxn_idx]
                complex_type['complex_type'] = cplx_type
                complex_information.append(reaction_name)
                complex_information.append(stoich)
                complex_information.append(stoich_known)
                complex_information.append(complex_type)

                # Append the complex name and stoich as a dictionary entry
                parent_complexes[complex_name] = complex_information

            molecules_to_parent_complexes_dict[subunit] = parent_complexes

        # Make a dictionary mapping between molecules to all downstream complexes
        # they form, both directly and indirectly (via another complex):
        for subunit in subunit_names:
            # Find the index in the matrix where this subunit is as a reactant:
            subunit_index = molecule_names.index(subunit)

            # Find the indices of complexes containing the subunit:
            complex_indices = np.where(smm[subunit_index, :] < 0)[0]

            downstream_complexes = {}
            for complex_idx in complex_indices:
                # Find the complex's name:
                complex_name = self.complexIDs[complex_idx]

                # Obtain the index of the complex within molecule_names:
                cplx_idx = molecule_names.index(complex_name)

                # Use the stoichMatrix() to find the reaction index:
                reaction_idx = np.where(sm[cplx_idx, :] > 0)[0]

                # Initialize variables to hold info about the complex:
                downstream_complex_information = []
                stoich = {}
                stoich_known = {}
                complex_type = {}
                reaction_name = {}

                # Add complex information to lists:
                reaction_name['reaction_id'] = rxn_ids[reaction_idx[0]]
                stoich['stoichiometry'] = smm[subunit_index, complex_idx]
                stoich_known['stoich_unknown'] = rsu[reaction_idx[0]]
                complex_type['complex_type'] = complex_to_complex_type[complex_name]
                downstream_complex_information.append(reaction_name)
                downstream_complex_information.append(stoich)
                downstream_complex_information.append(stoich_known)
                downstream_complex_information.append(complex_type)

                # Append the complex name and stoich as a dictionary entry
                downstream_complexes[complex_name] = downstream_complex_information

            molecules_to_all_downstream_complexes_dict[subunit] = downstream_complexes

        return (molecules_to_parent_complexes_dict,
                molecules_to_all_downstream_complexes_dict,
                complex_to_complex_type)


    def build_equilibrium_dictionaries(self, sim_data):
        """
        Build dictionaries mapping molecules to equilibirum complexation
        complexes they form. Note: complexation complexes can be subunits of
        equilibrium complexes.

        Args:
            sim_data: Simulation data object
        Returns:
            tuple: (molecules_to_parent_complexes_dict,
                    molecules_to_all_downstream_complexes_dict,
                    complex_to_complex_type)

            - molecules_to_parent_complexes_dict: {subunit_id: {complex_id: [
                  {'reaction_id': str},
                  {'stoichiometry': int},  # Negative value (consumed)
                  {'stoich_unknown': bool},
                  {'complex_type': 'homogeneous' or 'heterogeneous'}]}}
            - molecules_to_all_downstream_complexes_dict: Same structure but
            maps from complex to base monomers (bypassing intermediate complexes)
            - complex_to_complex_type: {complex_id: 'homogeneous' or 'heterogeneous'}
        """
        molecule_names = sim_data.process.equilibrium.molecule_names
        rxn_ids = sim_data.process.equilibrium.rxn_ids
        rsu = sim_data.process.equilibrium.reaction_stoichiometry_unknown
        sm = sim_data.process.equilibrium.stoich_matrix()
        smm = sim_data.process.equilibrium.stoich_matrix_monomers()

        c_molecule_names = sim_data.process.complexation.molecule_names
        c_smm = sim_data.process.complexation.stoich_matrix_monomers()

        molecules_to_parent_complexes_dict = {}
        molecules_to_all_downstream_complexes_dict = {}
        complex_to_complex_type = {}

        # Generate dictionary mapping equilibrium molecules to the direct parent
        # complexes they form:
        for subunit in molecule_names:
            # Find the matrix index where this subunit is as a molecule:
            subunit_index = molecule_names.index(subunit)

            # Find the indices for reactions the subunit is a reactant in:
            reaction_indices = np.where(sm[subunit_index,:] < 0)[0]

            if len(reaction_indices) == 0:
                # skip molecules that do not form any complexes
                continue

            # For each reaction index, find the complex(es) that is(are) formed:
            parent_complexes = {}
            for rxn_idx in reaction_indices:
                # Initialize variables to hold info about the complex:
                complex_information = []
                stoich = {} # NOTE: this will be the subunit's stoichiometry!
                stoich_known = {}
                complex_type = {}
                reaction_name = {}

                # Find the complex formed in the reaction (each reaction forms
                # one equilibrium complex):
                complex_index = np.where(sm[:, rxn_idx] > 0)[0]
                complex_name = molecule_names[complex_index[0]]
                cplx_idx = self.eqComplexIDs.index(complex_name)

                # Find the number of unique monomer subunits in the complex
                # (i.e. not counting metabolite ligands):
                unique_subunits_in_complex = np.where(smm[:, cplx_idx] < 0)[0]
                subunit_list = []
                for idx in unique_subunits_in_complex:
                    subunit_name = molecule_names[idx]
                    if subunit_name in self.monomerIDs:
                        if subunit_name not in subunit_list:
                            subunit_list.append(subunit_name)
                    # Check if the reactant is an eq complex:
                    if subunit_name in self.eqComplexIDs:
                        # Find the index of the eq complex in the molecule names list:
                        eq_complex_index = molecule_names.index(subunit_name)
                        # Append the index to the unique subunits list (as it
                        # could also be a complexation complex):
                        if eq_complex_index not in unique_subunits_in_complex:
                            unique_subunits_in_complex = (
                                np.append(unique_subunits_in_complex, eq_complex_index))
                    # Check if the reactant is a complexation complex:
                    elif subunit_name in self.complexIDs:
                        # See how many unique monomer subunits the complex has:
                        complexation_complex_index = self.complexIDs.index(subunit_name)
                        complexation_unique_subunits = (
                            np.where(c_smm[:, complexation_complex_index] < 0))[0]
                        # Add the unique monomer subunits of the complexation
                        # complex to the list of unique subunits in the eq complex:
                        for idx in complexation_unique_subunits:
                            c_su_name = c_molecule_names[idx]
                            if c_su_name not in subunit_list:
                                subunit_list.append(c_su_name)
                    # If the reactant is a metabolite ligand, skip it:
                    else:
                        continue

                num_unique_subunits = len(subunit_list)
                if num_unique_subunits > 1:
                    cplx_type = 'heterogeneous'
                else:
                    cplx_type = 'homogeneous'
                complex_to_complex_type[complex_name] = cplx_type

                # Add info about the complex to the parent complex info lists:
                reaction_name['reaction_id'] = rxn_ids[rxn_idx]
                stoich['stoichiometry'] = sm[subunit_index, rxn_idx]
                stoich_known['stoich_unknown'] = rsu[rxn_idx]
                complex_type['complex_type'] = cplx_type
                complex_information.append(reaction_name)
                complex_information.append(stoich)
                complex_information.append(stoich_known)
                complex_information.append(complex_type)

                # Append the complex name and stoich as a dictionary entry
                parent_complexes[complex_name] = complex_information

            molecules_to_parent_complexes_dict[subunit] = parent_complexes

        # Make a dictionary mapping between molecules to all downstream complexes
        # they form, both directly and indirectly (via another complex):
        for subunit in molecule_names:
            # Find the index of this subunit in the matrix:
            subunit_index = molecule_names.index(subunit)

            # Find the indices of complexes containing the subunit:
            complex_indices = np.where(smm[subunit_index, :] < 0)[0]

            if len(complex_indices) == 0:
                # skip molecules that form no downsteam complexes
                continue

            downstream_complexes = {}
            for complex_idx in complex_indices:
                # Find the complex's name:
                complex_name = self.eqComplexIDs[complex_idx]

                # Obtain the index of the complex within molecule_names:
                cplx_idx = molecule_names.index(complex_name)

                # Find the reaction index:
                reaction_idx = np.where(sm[cplx_idx, :] > 0)[0]

                # Initialize variables to hold info about the complex:
                downstream_complex_information = []
                stoich = {}
                stoich_known = {}
                complex_type = {}
                reaction_name = {}

                # Add complex information to lists:
                reaction_name['reaction_id'] = rxn_ids[reaction_idx[0]]
                stoich['stoichiometry'] = smm[subunit_index, complex_idx]
                stoich_known['stoich_unknown'] = rsu[reaction_idx[0]]
                complex_type['complex_type'] = complex_to_complex_type[complex_name]
                downstream_complex_information.append(reaction_name)
                downstream_complex_information.append(stoich)
                downstream_complex_information.append(stoich_known)
                downstream_complex_information.append(complex_type)

                # Append the complex name and stoich as a dictionary entry
                downstream_complexes[complex_name] = downstream_complex_information

            molecules_to_all_downstream_complexes_dict[subunit] = downstream_complexes

        return (molecules_to_parent_complexes_dict,
                molecules_to_all_downstream_complexes_dict,
                complex_to_complex_type)

    def build_tcs_dictionaries(self, sim_data):
        """
        Build dictionaries mapping molecules to two component system (TCS) complexes
        they form. Note: complexation and equilibrium complexes can be subunits
        of TCS complexes.

        Args:
            sim_data: Simulation data object
        Returns:
            tuple: (molecules_to_parent_complexes_via_lbptrs,
                molecules_to_parent_complexes_via_ptrs,
                complexes_to_reactants_via_dephosphorylation_reactions,
                molecules_to_all_downstream_complexes_dict,
                complex_to_complex_type)
            - read the descriptions in the function for details on the structure
                of molecules_to_parent_complexes_via_lbptrs,
                molecules_to_parent_complexes_via_ptrs,
                complexes_to_reactants_via_dephosphorylation_reactions
            - molecules_to_all_downstream_complexes_dict: Same structure but
                maps from complex to base monomers (bypassing intermediate complexes)
            - complex_to_complex_type: {complex_id: 'PHOSPHO-HK' or 'PHOSPHO-RR'}
        """
        molecule_names = list(sim_data.process.two_component_system.molecule_names)
        modified_molecules = list(sim_data.process.two_component_system.modified_molecules)
        complex_ids = list(sim_data.process.two_component_system.complex_to_monomer.keys())
        molecule_types = sim_data.process.two_component_system.molecule_types
        rxn_ids = sim_data.process.two_component_system.rxn_ids
        sm = sim_data.process.two_component_system.stoich_matrix()
        smm = sim_data.process.two_component_system.stoich_matrix_monomers()
        molecules_to_skip = ["ATP[c]", "ADP[c]", "Pi[c]", "WATER[c]", "PROTON[c]"]
        pairings = {"RR": "PHOSPHO-RR", "PHOSPHO-RR": "RR", "HK": "PHOSPHO-HK",
                    "PHOSPHO-HK": "HK", "HK-LIGAND": "PHOSPHO-HK-LIGAND",
                    "PHOSPHO-HK-LIGAND": "HK-LIGAND"}

        molecules_to_all_downstream_complexes_dict = {}
        complex_to_complex_type = {}

        # NOTE: for positively oriented TCS systems, a subunit can form a TCS
        # complex via BOTH a POS-LIGAND-BOUND-HK-PHOSPHOTRANSFER reaction AND
        # a POS-HK-PHOSPHOTRANSFER reaction! Thus, two dictionaries need to be
        # generated to avoid overwriting previous dictionary entries!

        # Dictionary that maps non-phosphorylated subunits to the TCS complexes
        # they form via ligand-bound HK phosphotransfer reactions:
        molecules_to_parent_complexes_via_lbptrs = {}

        # Dictionary that maps non-phosphorylated subunits to the TCS complexes they
        # form via ligand-free HK phosphotransfer reactions:
        molecules_to_parent_complexes_via_ptrs = {}

        # Also need to model dephosphorylation reactions, where the
        # phosphorylated complex (HK or RR) is the reactant:
        complexes_to_reactants_via_dephosphorylation_reactions = {}

        # Generate mapping from molecules to direct parent complexes they form:
        for subunit in molecule_names:
            # If the molecule is a TCS complex or not a complex or monomer, skip it:
            if subunit in molecules_to_skip or subunit in complex_ids:
                continue

            # Find the matrix index where this subunit is as a molecule:
            subunit_index = molecule_names.index(subunit)
            reactant_type = molecule_types[molecule_names.index(subunit)]
            compliment = pairings[reactant_type]

            # Find the reaction index within the stoich matrix:
            reaction_indices = np.where(sm[subunit_index, :] < 0)[0]

            # For each reaction index, find the product(s) that is(are) formed:
            parent_complexes_lbptrs = {}
            parent_complexes_ptrs = {}
            for rxn_idx in reaction_indices:
                # Initialize data structures to hold product information
                product_information = []
                stoich = {} # NOTE: this will be the subunit's stoichiometry!
                stoich_known = {}
                product_type = {}
                reaction_name = {}

                # Determine if this is a ligand-bound or ligand-free HK phosphotransfer reaction:
                rxn_id = rxn_ids[rxn_idx]
                if "BOUND" in rxn_id:
                    parent_dict_to_use = parent_complexes_lbptrs
                else:
                    parent_dict_to_use = parent_complexes_ptrs

                # Find the products formed in this reaction:
                product_indices = np.where(sm[:, rxn_idx] > 0)[0]

                # TCS reactions can have multiple products, so they need to be filtered:
                for product_index in product_indices:
                    product_name = molecule_names[product_index]

                    # Since there are multiple products, determine which product
                    # is the compliment of the specific reactant subunit:
                    pdt_type = molecule_types[molecule_names.index(product_name)]
                    if pdt_type != compliment:
                        continue
                    else:
                        # Save the complex type in a separate dict for easy access later:
                        complex_to_complex_type[product_name] = pdt_type

                        # Add product's information to lists:
                        reaction_name['reaction_id'] = rxn_id
                        stoich['stoichiometry'] = sm[subunit_index, rxn_idx]
                        stoich_known['stoich_unknown'] = 'not applicable'
                        product_type['complex_type'] = pdt_type
                        product_information.append(reaction_name)
                        product_information.append(stoich)
                        product_information.append(stoich_known)
                        product_information.append(product_type)

                        # Append the product name and stoich as a dictionary entry
                        parent_dict_to_use[product_name] = product_information

                molecules_to_parent_complexes_via_lbptrs[subunit] = parent_complexes_lbptrs
                molecules_to_parent_complexes_via_ptrs[subunit] = parent_complexes_ptrs

        # Make a modified version that specifically handles dephosphorylation
        # reactions where phosphorylated TCS complexes are reactants:
        for subunit in complex_ids:
            # Find the matrix index where this subunit is as a molecule:
            subunit_index = molecule_names.index(subunit)

            # Find the reaction index within the stoich matrix:
            reaction_indices = np.where(sm[subunit_index, :] < 0)[0]

            # Find the type of the reactant complex:
            reactant_type = molecule_types[molecule_names.index(subunit)]
            compliment = pairings[reactant_type]

            # For each reaction index, find the complex(es) that is(are) formed:
            products = {}
            for rxn_idx in reaction_indices:
                # Initialize data structures to hold product information
                product_information = []
                stoich = {}  # NOTE: this will be the complexes's stoichiometry!
                stoich_known = {}
                product_type = {} # NOTE: this is really the subunit type
                reaction_name = {}

                # Find the products formed in this reaction:
                product_indices = np.where(sm[:, rxn_idx] > 0)[0]

                # TCS reactions can have multiple products, so they need to be filtered:
                for product_index in product_indices:
                    product_name = molecule_names[product_index]

                    # Since there are multiple products, determine which product
                    # is the compliment of the specific reactant subunit:
                    pdt_type = molecule_types[molecule_names.index(product_name)]
                    if pdt_type != compliment:
                        continue
                    else:
                        # Add product's information to lists:
                        reaction_name['reaction_id'] = rxn_ids[rxn_idx]
                        stoich['stoichiometry'] = sm[subunit_index, rxn_idx]
                        stoich_known['stoich_unknown'] = 'not applicable'
                        product_type['complex_type'] = \
                            (molecule_types[molecule_names.index(product_name)])
                        product_information.append(reaction_name)
                        product_information.append(stoich)
                        product_information.append(stoich_known)
                        product_information.append(product_type)

                        # Append the product name and stoich as a dictionary entry
                        products[product_name] = product_information

                complexes_to_reactants_via_dephosphorylation_reactions[subunit] = products

        # Make a dictionary mapping relevant reactants to the phosphorylated
        # TCS products they form:
        for subunit in modified_molecules:
            # If the molecule is not a complex or monomer, skip it:
            if subunit in molecules_to_skip:
                continue

            # Find the matrix index where this subunit is as a molecule:
            subunit_index = modified_molecules.index(subunit)

            # Find the indices of complexes containing the subunit:
            complex_indices = np.where(smm[subunit_index, :] < 0)[0]

            if len(complex_indices) == 0:
                # skip molecules that are not base subunits of a TCS complex
                continue

            downstream_complexes = {}
            for complex_idx in complex_indices:
                # Find the complex's name:
                complex_name = complex_ids[complex_idx]

                # Obtain the index of the complex within self.molecule_names:
                cplx_idx_molecule = molecule_names.index(complex_name)

                # Use the stoichMatrix() to find the reaction index:
                reaction_idx = np.where(sm[cplx_idx_molecule, :] > 0)[0]

                # Initialize data structures to hold info about the complex:
                downstream_complex_information = []
                stoich = {}
                stoich_known = {}
                complex_type = {}
                reaction_name = {}

                # Add complex information to lists:
                reaction_name['reaction_id'] = rxn_ids[reaction_idx[0]]
                stoich['stoichiometry'] = smm[subunit_index, complex_idx]
                stoich_known['stoich_unknown'] = "not applicable"
                complex_type['complex_type'] = (
                    molecule_types)[molecule_names.index(complex_name)]
                downstream_complex_information.append(reaction_name)
                downstream_complex_information.append(stoich)
                downstream_complex_information.append(stoich_known)
                downstream_complex_information.append(complex_type)

                # Append the complex name and stoich as a dictionary entry
                downstream_complexes[complex_name] = downstream_complex_information

            molecules_to_all_downstream_complexes_dict[subunit] = downstream_complexes

        return (molecules_to_parent_complexes_via_lbptrs,
                molecules_to_parent_complexes_via_ptrs,
                complexes_to_reactants_via_dephosphorylation_reactions,
                molecules_to_all_downstream_complexes_dict,
                complex_to_complex_type)

    def build_tcs_reaction_event_estimates(self, sim_data, cell_paths):
        """
        Estimate reaction fluxes from TCS molecule changes using non-negative
        least squares.

        NOTE: there could technically be many possible reaction flux solutions
        that could product a fit similar to molecule_changes, but the nnl method
        will optimize for the sparset solution involving the least number
        of reactions total to hit the desired molecule changes. This is a
        reasonable assumption to make given the small number of reactions and
        given that molecule_changes are already weighted by the rates of reaction.

        If this does not appear to be working well because the reaction rates
        change or something, it might be best to try a different method.

        Args:
            molecule_changes: (n_molecules,) array of molecule count changes
            from the ODE integration
        Returns:
            reaction_fluxes: (n_reactions,) array of ESTIMATED reaction counts
        """
        # Extract the molecule changes for each TCS molecule across the cell paths:
        all_molecule_changes = read_stacked_columns(
            cell_paths, 'TwoComponentSystems', 'moleculeChanges')
        sm = sim_data.process.two_component_system.stoich_matrix()
        n_timesteps = all_molecule_changes.shape[0]
        n_reactions = sm.shape[1]

        # Estimate the reaction fluxes for each step using non-negative least squares:
        estimated_reaction_fluxes = np.zeros((n_timesteps, n_reactions))

        for i in range(n_timesteps):
            reaction_fluxes, residual = scipy.optimize.nnls(sm, all_molecule_changes[i])
            estimated_reaction_fluxes[i] = np.round(reaction_fluxes).astype(int)

            # Print a warning if the fit was not good::
            reconstructed = sm.dot(estimated_reaction_fluxes[i])
            max_error = np.max(np.abs(reconstructed - all_molecule_changes[i]))

            if max_error > 0.5:
                print(f"Warning: reaction flux reconstruction error (t={i}) = {max_error}")
                print(f"Residual from NNLS: {residual}")

        return estimated_reaction_fluxes

    def find_extra_tcs_reactions(self, molecule, complex, rxns_plotted):
        """
        Since some of TCS complexes can be formed by two different types of
        reactions via the same reactant, this function checks whether the TCS
        complex being plotted can also be formed via the other reaction type,
        and if so, plots the events for that reaction as well (if they have not
         already been plotted).

        This is important for accurately interpreting the events plot since it
        shows the total events for all reactions that form the complex, and if
        one of the reactions is skipped, then the events will be undercounted
        and the plot will be misleading.
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

        # Determine if the complex is a TF:
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
                             + ["CPLX0-3956[c]","CPLX0-3955[c]"]
                             + [sim_data.molecule_ids.s30_full_complex]
                             + [sim_data.molecule_ids.s50_full_complex])

        # Determine if the complex is a ribosome subunit:
        ribosome_subunit_id_match = [subunit_id for subunit_id in
                                     ribosome_subunits if molecule == subunit_id]

        if ribosome_subunit_id_match != []:
            print(f"{molecule} is a ribosome subunit.")

        return ribosome_subunit_id_match

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

        # Determine if the complex is an RNAP subunit:
        rnap_subunit_id_match = [subunit_id for subunit_id in
                                 rnap_subunit_ids if molecule == subunit_id]

        if rnap_subunit_id_match != []:
            print(f"{molecule} is an RNAP subunit.")

        # Check if the molecule is the inactive RNAP complex:
        if molecule == "APORNAP-CPLX[c]":
            rnap_subunit_id_match += ["APORNAP-CPLX[c]"]
            print(f"{molecule} is the inactive RNAP complex.")

        return rnap_subunit_id_match

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


    def is_a_replisome_subunit(self, sim_data, molecule):
        """
        Determine whether a molecule is an RNAP subunit.
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

        # Determine if the complex is a replisome subunit:
        replisome_subunit_id_match = [subunit_id for subunit_id in
                replisome_subunit_ids if molecule == subunit_id]

        if replisome_subunit_id_match != []:
            print(f"{molecule} is a replisome subunit.")

        return replisome_subunit_id_match

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

            # Find the replisome complexation events by subtracting the counts
            # of the previous time step from the current time step:
            replisome_complexation_events = np.diff(active_replisome_counts, prepend=0)

            # Delete all negative values (which correspond to replisomes
            # deleting from the simulations):
            replisome_complexation_events[replisome_complexation_events < 0] = 0
            replisome_events.append(replisome_complexation_events)

        # Concatenate the replisome counts and events from each simulation into one array:
        replisome_counts = np.concatenate(replisome_counts)
        replisome_events = np.concatenate(replisome_events)

        # Get the replisome subunit IDs:
        replisome_trimer_subunits = sim_data.molecule_groups.replisome_trimer_subunits
        replisome_monomer_subunits = sim_data.molecule_groups.replisome_monomer_subunits
        replisome_subunit_ids = replisome_trimer_subunits + replisome_monomer_subunits

        # Generate the replisome stoichiometry matrix:
        replisome_stoich = np.hstack(
            (3 * np.ones(len(replisome_trimer_subunits)),
             np.ones(len(replisome_monomer_subunits))))

        # Find the index of the molecule within the replisome subunit list:
        replisome_subunit_index = replisome_subunit_ids.index(molecule)

        # Determine the stoich of the molecule within the replisome complex and
        # the counts of the molecule in active replisomes at each moment:
        replisome_subunit_stoich = replisome_stoich[replisome_subunit_index]

        return replisome_subunit_stoich, replisome_counts, replisome_events


    def do_plot(self, seedOutDir, plotOutDir, plotOutFileName, simDataFile,
                validationDataFile, metadata):
        with open(simDataFile, 'rb') as f:
            sim_data = pickle.load(f)
        cell_paths = self.ap.get_cells()
        sim_dir = cell_paths[0]
        simOutDir = os.path.join(sim_dir, 'simOut')

        # Extract protein indexes for each monomer:
        monomer_counts_reader = TableReader(
            os.path.join(simOutDir, "MonomerCounts"))
        monomer_idx_dict = {monomer: i for i, monomer in
                            enumerate(monomer_counts_reader.readAttribute(
                                'monomerIds'))}

        # Extract the monomer IDs:
        self.monomerIDs = monomer_counts_reader.readAttribute("monomerIds")

        # Extract the free monomer counts:
        free_monomer_counts = read_stacked_columns(
            cell_paths, 'MonomerCounts', "freeMonomerCounts")

        # Extract the doubling times:
        (time, doubling_times, end_generation_times, start_generation_indices,
         end_generation_indices) = self.extract_doubling_times(cell_paths)

        # Extract the IDs of the three types of complexes:
        self.complexIDs = sim_data.process.complexation.ids_complexes
        self.eqComplexIDs = sim_data.process.equilibrium.ids_complexes
        self.tcsComplexIDs = list(sim_data.process.two_component_system.complex_to_monomer.keys())

        # Generatae the complexation molecules to parent complexes
        # and molecules to all downstream complexes dictionaries:
        self.c_m2pc_dict, self.c_m2adc_dict, self.c2ct = (
            self.build_complexation_dictionaries(sim_data))

        # Warning: these will not unpack subunits to monomers if the subunit
        # is a complexation complex!
        self.eq_m2pc_dict, self.eq_m2adc_dict, self.eq_c2ct = (
            self.build_equilibrium_dictionaries(sim_data))

        # Warning: these will not unpack subunits to baseline monomers if the
        # subunit is a complexation or equilibrium complex!
        (self.tcs_m2pc_via_lbptrs_dict, self.tcs_m2pc_via_ptrs_dict,
         self.tcs_c2m_via_dprs_dict, self.tcs_m2adc_dict, self.tcs_c2ct)\
            = self.build_tcs_dictionaries(sim_data)

        # Define a TCS molecules to parent complexes dictionary to be used in
        # other functions, given that intermediate complex stoich is always 1:1:
        self.tcs_m2pc_dict = {}
        for molecule in set(list(self.tcs_m2pc_via_lbptrs_dict.keys()) + list(
                self.tcs_m2pc_via_ptrs_dict.keys())):
            self.tcs_m2pc_dict[molecule] = {}

            if molecule in self.tcs_m2pc_via_lbptrs_dict:
                self.tcs_m2pc_dict[molecule].update(self.tcs_m2pc_via_lbptrs_dict[molecule])

            if molecule in self.tcs_m2pc_via_ptrs_dict:
                self.tcs_m2pc_dict[molecule].update(self.tcs_m2pc_via_ptrs_dict[molecule])

        # Estimate the TCS reaction events:
        self.tcs_reaction_events = (
            self.build_tcs_reaction_event_estimates(sim_data, cell_paths))

        # Make sure thee inputted entries are valid and assign a compartment tag:
        PLOT_COMPLEXES_revised = (
            self.check_validity_and_get_compartment(sim_data, PLOT_COMPLEXES))

        # Check that all the valid inputs are indeed complexes:
        valid_complexes, complex_type_dict = (
            self.check_complex_validity(PLOT_COMPLEXES_revised))

        # Generate a plot for each complex:
        for complex_id in valid_complexes:
            # Extract the upstream monomers and sentences describing their paths:
            upstream_monomers, upstream_monomer_paths = self.analyze_complex(complex_id)

            # Extract all downstream complexes the complex is a subunit of:
            downstream_complexes, downstream_complex_paths = (
                self.analyze_molecule_parents(complex_id))

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
            if self.is_a_ribosome_subunit(sim_data, complex_id) == []:
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

            # Also plot the counts of each base monomer within the complex:
            monomer_ids = list(upstream_monomers.keys())
            for i in range(len(monomer_ids)):
                monomer = monomer_ids[i]
                path = upstream_monomer_paths[i]  # Get corresponding path from the list
                total_stoich = upstream_monomers[monomer][0]['total_stoich']
                monomer_counts = complex_counts * total_stoich

                # Plot the counts within the complex itself:
                ax1.plot(time, monomer_counts, alpha=1,
                         label=path,
                         linewidth=1, linestyle=":")

                # Second, plot the free monomer subunits on the available free
                # monomer subunit counts plot:
                protein_idx = monomer_idx_dict[monomer]
                protein_FMC = free_monomer_counts[:, protein_idx]
                ax2.plot(time, protein_FMC, alpha=0.75,
                         label=f'{monomer}', linewidth=0.75)

            # Add the parent complexes to the complex counts plot as well:
            complexes = list(downstream_complexes.keys())
            for i in range(len(complexes)):
                parent_complex = complexes[i]
                path = downstream_complex_paths[i]
                total_stoich = downstream_complexes[parent_complex][0]['total_stoich']
                parent_complex_free_counts = (
                    read_stacked_bulk_molecules(cell_paths, [parent_complex]))[0]
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
            for i in range(len(monomer_ids)):
                monomer = monomer_ids[i]
                # Find the reaction that makes the complex from this monomer:
                rxn_id = upstream_monomers[monomer][0]['path'][-1][3]
                if rxn_id in rxns_plotted:
                    continue # skip if the reacion has been plotted
                else:
                    rxns_plotted.append(rxn_id)
                rxn_stoich = 1 # always 1 for complex formation
                name = f'{rxn_id} \n(generates {rxn_stoich} {complex_id})'

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
                    ax3.plot(time, rxn_events, alpha=0.5, label=name, linewidth=0.75)
                else: # TCS reactions
                    rxn_idx = tcs_rxn_ids.index(rxn_id)
                    rxn_events = self.tcs_reaction_events[:, rxn_idx]
                    ax3.plot(time, rxn_events, alpha=0.5, label=
                            f'{rxn_id} \n estimated events '
                            f'(generates 1 {complex_id})', linewidth=0.75)

                    # Also check if the complex is made from a different TCS
                    # reaction by the same base monomer:
                    if "POS" in rxn_id:
                        reaction_IDs = (self.find_extra_tcs_reactions(
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
                    if complex_id in self.tcs_c2m_via_dprs_dict.keys():
                        subunit_complex_info = self.tcs_c2m_via_dprs_dict[complex_id]
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
            for i in range(len(complexes)):
                parent_complex = complexes[i]
                rxn_id = downstream_complexes[parent_complex][0]['path'][0][3]
                if rxn_id in consumption_rxns_plotted or rxn_id in rxns_plotted:
                    continue # skip if the reacion has been plotted
                else:
                    consumption_rxns_plotted.append(rxn_id)
                rxn_stoich = downstream_complexes[parent_complex][0]['path'][0][1]
                name = f'{rxn_id} \n(consumes {rxn_stoich} {complex_id})'

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
                    ax3.plot(time, np.negative(rxn_events),
                             alpha=0.5, label=name, linewidth=0.75)

                else:
                    rxn_idx = tcs_rxn_ids.index(rxn_id)
                    rxn_events = self.tcs_reaction_events[:, rxn_idx]
                    ax3.plot(time, np.negative(rxn_events),
                             alpha=0.5, label=f'{rxn_id} \n estimated events '
                                 f'(consumes 1 {complex_id})', linewidth=0.75)

                    # Also check if the parent is made from a different TCS
                    # reaction by the same base complex:
                    if "POS" in rxn_id:
                        reaction_IDs = self.find_extra_tcs_reactions(
                            complex_id, parent_complex, consumption_rxns_plotted)
                        for reaction_ID in reaction_IDs:
                            # Plot the reaction events based on the type of complex:
                            consumption_rxns_plotted.append(reaction_ID)
                            rxn_idx = tcs_rxn_ids.index(reaction_ID)
                            rxn_events = self.tcs_reaction_events[:, rxn_idx]
                            ax3.plot(time, np.negative(rxn_events), alpha=0.5, label=
                            f'{reaction_ID} \n estimated events (consumes 1 {complex_id})', linewidth=0.75)

                    # Lastly, check if the parent dephosphorylates to the complex:
                    if parent_complex in self.tcs_c2m_via_dprs_dict.keys():
                        subunit_complex_info = self.tcs_c2m_via_dprs_dict[parent_complex]
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
            if self.is_a_TF(sim_data, complex_id) != []:
                # If the complex is a TF, plot its bound counts and binding/unbinding events:
                tf_idx = self.tf_index(sim_data, complex_id)
                tf_counts = read_stacked_columns(cell_paths,'RnaSynthProb',
                                                "nActualBound")[:,tf_idx]
                tfs_unbound = read_stacked_columns(cell_paths,'RnaSynthProb',
                                                "nActualUnbound")[:,tf_idx]
                ax1.plot(time, tf_counts, color='lightcoral',
                         label=f'transcription unit bound\n {complex_id} '
                               f'transcription factors',
                         linewidth=0.75, alpha=0.75)
                ax3.plot(time, (tf_counts * -1), color='lightcoral',
                         label=f'{complex_id} Transcription factor\nbinding events',
                         linewidth=0.75, alpha=0.75)
                ax3.plot(time, tfs_unbound, color='purple',
                         label=f'{complex_id} Transcription factor\nunbinding events',
                         linewidth=0.75, alpha=0.75)
                total_counts += tf_counts

            # Check if the complex is a replisome subunit:
            if self.is_a_replisome_subunit(sim_data, complex_id) != []:
                # If the complex is a replisome subunit, plot its counts within
                # active replisomes and the replisome assembly events that consume it:
                replisome_stoich, replisome_counts, replisome_events = (
                    self.get_replisome_subunit_counts(sim_data, cell_paths, complex_id))
                subunit_counts_within_replisomes = (
                        replisome_counts * replisome_stoich)
                ax1.plot(time, subunit_counts_within_replisomes, color='gold',
                         label=f'{complex_id} counts within active replisomes, '
                               f'{replisome_stoich} per',
                         linewidth=0.75, alpha=0.75, linestyle="--")
                ax3.plot(time, np.negative(replisome_events), color='gold',
                         label=f'Active replisome assebly events\n(consumes '
                               f'{replisome_stoich} {complex_id})',
                         linewidth=0.75, alpha=0.75)
                # Add the replisome subunit counts to the total counts of the complex:
                total_counts += (subunit_counts_within_replisomes)

            # Check if the complex is an inactive RNAP (1:1 ratio):
            if self.is_an_RNAP_subunit(sim_data, complex_id) != []:
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
            if self.is_a_ribosome_subunit(sim_data, complex_id) != []:
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
                if self.is_a_replisome_subunit(sim_data, parent_complex) != []:
                    replisome_stoich, replisome_counts, _ = (
                        self.get_replisome_subunit_counts(
                            sim_data, cell_paths, parent_complex))

                    # Get the stoich of the parent complex to the complex of interest:
                    total_stoich = downstream_complexes[parent_complex][0]['total_stoich']

                    # Multiply the replisome subunit counts by the stoich of
                    # the parent complex and the stoich of the parent complex in the replisome:
                    subunit_counts = replisome_counts * total_stoich * replisome_stoich

                    ax1.plot(time, subunit_counts, color='gold',
                             label=f'{complex_id} counts within active replisomes,'
                                   f'\n{replisome_stoich*replisome_stoich} per '
                                   f'(via {parent_complex})',
                             linewidth=0.75, alpha=0.75, linestyle="--")

                    # Add the replisome counts to the total counts of the complex:
                    total_counts += subunit_counts

                if self.is_a_TF(sim_data, parent_complex) != []:
                    tf_idx = self.tf_index(sim_data, parent_complex)
                    tf_counts = read_stacked_columns(cell_paths,
                        'RnaSynthProb',"nActualBound")[:,tf_idx]

                    # Get the stoich of the parent complex to the complex of interest:
                    total_stoich = downstream_complexes[parent_complex][0]['total_stoich']

                    # Multiply the TF counts by the stoich of the complex to
                    # the parent complex:
                    subunit_counts = tf_counts * total_stoich
                    ax1.plot(time, subunit_counts, color='lightcoral',
                             label=f'transcription unit bound {complex_id} transcription'
                                   f' factors,\n {total_stoich} per (via {parent_complex})',
                             linewidth=0.75, alpha=0.75, linestyle="--")

                    # Add the TF counts to the total counts of the complex:
                    total_counts += subunit_counts

                # No need to check RNAPs, they only have monomer subunits.

                if self.is_a_ribosome_subunit(sim_data, parent_complex) != []:
                    ribosome_counts, _, _ = (
                        self.get_ribosome_subunit_counts(cell_paths))

                    # Get the stoich of the parent complex to the complex of interest:
                    total_stoich = downstream_complexes[parent_complex][0]['total_stoich']

                    # Multiply the ribosome subunit counts by the stoich of
                    # the parent complex and the stoich of the parent complex in the ribosome:
                    subunit_counts = ribosome_counts * total_stoich

                    ax1.plot(time, subunit_counts, color='gold',
                             label=f'{complex_id} counts within active ribosomes,'
                                   f'\n{total_stoich} per (via {parent_complex})',
                             linewidth=0.75, alpha=0.75, linestyle="--")

                    # Add the ribosome counts to the total counts of the complex:
                    total_counts += subunit_counts


            # Plot the total theoretical complex counts:
            ax1.plot(time, total_counts, color='black',
                     label=f'Theoretical total {complex_id} counts\n(free complex '
                           f'counts + counts within parent complexes)',
                     linewidth=.25, alpha=1)

            # Set the axes and legends :
            ax1.legend(fontsize=5, loc="center left", bbox_to_anchor=(1, 0.5))
            ax1.set_ylabel("Complex Counts")
            ax1.set_title(
                f"Complex counts for the {complex_makeup} {complex_type} complex "
                f"{complex_id}\n Sim ID: {metadata['description']}", fontsize=6)
            ax2.set_ylabel(f"Free Monomer \nSubunit Availability")
            ax2.legend(fontsize=5, loc="center left", bbox_to_anchor=(1, 0.5))

            ax3.set_xlabel("Time (s)")
            ax3.set_ylabel("Free Complex\nReaction Events")
            ax3.legend(fontsize=5, loc="center left", bbox_to_anchor=(1, 0.5))

            # Add vertical lines for the end of each generation:
            for i in range(len(end_generation_times)):
                dt = end_generation_times[i]
                a = 0.5
                ax1.axvline(x=dt, linestyle='--', color="yellowgreen", alpha=a)
                ax2.axvline(x=dt, linestyle='--', color="yellowgreen", alpha=a)
                ax3.axvline(x=dt, linestyle='--', color="yellowgreen", alpha=a)

            if self.is_a_ribosome_subunit(sim_data, complex_id) == []:
                plt.tight_layout()

            # Save the plot:
            file_name = plotOutFileName + "_" + complex_id
            exportFigure(plt, plotOutDir, file_name, metadata)


if __name__ == '__main__':
    Plot().cli()
