"""
Shared helpers for analysis plots that need to know which complexes and active
molecules (ribosomes, RNAPs, replisomes, DNA-bound transcription factors (TFs))
a molecule belongs to:
- build_subunit_categories() groups monomers by what they are a subunit of.
- ComplexNetwork holds the reaction dictionaries linking each molecule to the
  complexation, equilibrium and two component system (TCS) complexes it forms,
  and traces a molecule's assembly paths through them.
- The is_a_*()/is_an_*() functions check one molecule ID against sim_data.
"""
import numpy as np
import scipy.optimize


def check_validity_and_get_compartment(sim_data, molecule_list):
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


def build_subunit_categories(sim_data, monomer_dict):
    """
    Group monomers by the complexes and active molecules they are subunits of.
    Monomers not in monomer_dict are ignored.

    Args:
        sim_data: Simulation data object
        monomer_dict: {monomer_id: monomer_idx} of the monomers to categorize
    Returns:
        dict with:
        - the group membership maps, each {monomer_id: monomer_idx}:
          'ribosome_monomers', 'rnap_monomers', 'replisome_monomers',
          'tf_monomers', 'complexation_monomers', 'equilibrium_monomers',
          'two_component_system_monomers'
        - 'tf_monomer_types': {monomer_id: set of the types (0CS, 1CS or 2CS)
          of the TFs the monomer is part of}
        - the complex to subunit maps, each {complex_id: [monomer_id, ...]}:
          'complexes_to_monomers', 'eq_complexes_to_monomers',
          'tcs_complexes_to_monomers'
    """
    # COMPLEXATION COMPLEXES
    complexation_molecule_ids = sim_data.process.complexation.molecule_names
    complexation_complex_ids = sim_data.process.complexation.ids_complexes
    complexation_smm = sim_data.process.complexation.stoich_matrix_monomers()
    complexation_monomers = {}
    complexes_to_monomers = {}

    for mol_id in complexation_molecule_ids:
        if mol_id in monomer_dict:
            complexation_monomers[mol_id] = monomer_dict[mol_id]
        elif mol_id in complexation_complex_ids:
            complex_idx = complexation_complex_ids.index(mol_id)
            monomer_idxs = np.where(complexation_smm[:, complex_idx] < 0)[0]
            monomers = []
            for idx in monomer_idxs:
                monomer_id = complexation_molecule_ids[idx]
                if monomer_id in monomer_dict:
                    complexation_monomers[monomer_id] = monomer_dict[monomer_id]
                    monomers.append(monomer_id)
            complexes_to_monomers[mol_id] = monomers

    # EQUILIBRIUM COMPLEXES
    equilibrium_molecule_ids = sim_data.process.equilibrium.molecule_names
    equilibrium_complex_ids = sim_data.process.equilibrium.ids_complexes
    eq_smm = sim_data.process.equilibrium.stoich_matrix_monomers()
    eq_complexes_to_monomers = {}
    equilibrium_monomers = {}

    for mol_id in equilibrium_molecule_ids:
        if mol_id in monomer_dict:
            equilibrium_monomers[mol_id] = monomer_dict[mol_id]
        elif mol_id in complexation_complex_ids:
            if mol_id in complexes_to_monomers:
                for monomer_id in complexes_to_monomers[mol_id]:
                    if monomer_id in monomer_dict:
                        equilibrium_monomers[monomer_id] = monomer_dict[monomer_id]
        elif mol_id in equilibrium_complex_ids:
            eq_complex_idx = equilibrium_complex_ids.index(mol_id)
            molecule_idxs = np.where(eq_smm[:, eq_complex_idx] < 0)[0]
            monomers = []
            for idx in molecule_idxs:
                molecule_name = equilibrium_molecule_ids[idx]
                if molecule_name in monomer_dict:
                    monomers.append(molecule_name)
                    equilibrium_monomers[molecule_name] = monomer_dict[molecule_name]
                elif molecule_name in complexation_complex_ids:
                    if molecule_name in complexes_to_monomers:
                        for monomer_id in complexes_to_monomers[molecule_name]:
                            if monomer_id in monomer_dict:
                                monomers.append(monomer_id)
                                equilibrium_monomers[monomer_id] = monomer_dict[monomer_id]
            eq_complexes_to_monomers[mol_id] = monomers

    # TWO COMPONENT SYSTEM COMPLEXES
    two_component_system_monomers = {}
    tcs_complexes_to_monomers = {}
    two_component_system_complex_ids = []

    try:
        two_component_system_molecule_ids = list(
            sim_data.process.two_component_system.modified_molecules)
        two_component_system_complex_ids = list(
            sim_data.process.two_component_system.complex_to_monomer.keys())
        tcs_complex_to_monomer_dict = sim_data.process.two_component_system.complex_to_monomer

        for mol_id in two_component_system_molecule_ids:
            subunit_ids = []
            if mol_id in monomer_dict:
                two_component_system_monomers[mol_id] = monomer_dict[mol_id]
            elif mol_id in complexation_complex_ids:
                subunit_ids = complexes_to_monomers[mol_id]
                tcs_complexes_to_monomers[mol_id] = subunit_ids
            elif mol_id in equilibrium_complex_ids:
                subunit_ids = eq_complexes_to_monomers[mol_id]
                tcs_complexes_to_monomers[mol_id] = subunit_ids
            elif mol_id in two_component_system_complex_ids:
                subunit_ids = list(tcs_complex_to_monomer_dict[mol_id].keys())
                tcs_complexes_to_monomers[mol_id] = subunit_ids
            for monomer_id in subunit_ids:
                if monomer_id in monomer_dict:
                    two_component_system_monomers[monomer_id] = monomer_dict[monomer_id]
    except AttributeError:
        # Print a warning if the two_component_system.modified_molecules
        # attribute is not found since the listener was added in 03/26
        print("Warning: two_component_system.modified_molecules not found "
              "(likely due to using an older sim version where sim_data version "
              "is out of date with current listeners used in the model).")

    # RIBOSOMES
    ribosome_50s_subunits = sim_data.process.complexation.get_monomers(
        sim_data.molecule_ids.s50_full_complex)
    ribosome_30s_subunits = sim_data.process.complexation.get_monomers(
        sim_data.molecule_ids.s30_full_complex)
    ribosome_subunits = (ribosome_50s_subunits["subunitIds"].tolist() +
                         ribosome_30s_subunits["subunitIds"].tolist())
    ribosome_monomers = {}
    for subunit in ribosome_subunits:
        if subunit in monomer_dict:
            ribosome_monomers[subunit] = monomer_dict[subunit]

    # RNAPs
    rnap_subunits = sim_data.process.complexation.get_monomers(
        sim_data.molecule_ids.full_RNAP)
    rnap_subunit_ids = rnap_subunits["subunitIds"].tolist()
    rnap_monomers = {}
    for subunit in rnap_subunit_ids:
        if subunit in monomer_dict:
            rnap_monomers[subunit] = monomer_dict[subunit]

    # REPLISOMES
    replisome_trimer_subunits = sim_data.molecule_groups.replisome_trimer_subunits
    replisome_monomer_subunits = sim_data.molecule_groups.replisome_monomer_subunits
    replisome_subunit_ids = replisome_trimer_subunits + replisome_monomer_subunits
    replisome_monomers = {}
    for subunit in replisome_subunit_ids:
        if subunit in monomer_dict:
            replisome_monomers[subunit] = monomer_dict[subunit]
        elif subunit in complexation_complex_ids:
            for mol_id in complexes_to_monomers[subunit]:
                if mol_id in monomer_dict:
                    replisome_monomers[mol_id] = monomer_dict[mol_id]

    # TRANSCRIPTION FACTORS
    tfs = sim_data.process.transcription_regulation.tf_ids
    tf_to_tf_type = sim_data.process.transcription_regulation.tf_to_tf_type
    tf_monomers = {}
    tf_monomer_types = {}
    for tf_id in tfs:
        subunit = tf_id + f'[{sim_data.getter.get_compartment(tf_id)[0]}]'
        if subunit in monomer_dict:
            monomer_ids_list = [subunit]
        elif subunit in complexation_complex_ids:
            monomer_ids_list = complexes_to_monomers[subunit]
        elif subunit in equilibrium_complex_ids:
            monomer_ids_list = eq_complexes_to_monomers[subunit]
        elif subunit in two_component_system_complex_ids:
            monomer_ids_list = tcs_complexes_to_monomers[subunit]
        else:
            monomer_ids_list = []
        for mol_id in monomer_ids_list:
            if mol_id in monomer_dict:
                tf_monomers[mol_id] = monomer_dict[mol_id]
                tf_monomer_types.setdefault(mol_id, set()).add(tf_to_tf_type[tf_id])

    return {
        'ribosome_monomers': ribosome_monomers,
        'rnap_monomers': rnap_monomers,
        'replisome_monomers': replisome_monomers,
        'tf_monomers': tf_monomers,
        'tf_monomer_types': tf_monomer_types,
        'complexation_monomers': complexation_monomers,
        'equilibrium_monomers': equilibrium_monomers,
        'two_component_system_monomers': two_component_system_monomers,
        'complexes_to_monomers': complexes_to_monomers,
        'eq_complexes_to_monomers': eq_complexes_to_monomers,
        'tcs_complexes_to_monomers': tcs_complexes_to_monomers,
    }


def get_subunit_category_names(monomer_id, subunit_categories):
    """
    List the groups in build_subunit_categories() a monomer belongs to.

    Args:
        monomer_id: ID of the monomer
        subunit_categories: output of build_subunit_categories()
    Returns:
        list: group names, in a fixed order (empty if the monomer is in none)
    """
    category_names = (
        ('ribosome_monomers', 'Ribosome'),
        ('rnap_monomers', 'RNAP'),
        ('replisome_monomers', 'Replisome'),
        ('tf_monomers', 'Transcription factor'),
        ('complexation_monomers', 'Complexation complex'),
        ('equilibrium_monomers', 'Equilibrium complex'),
        ('two_component_system_monomers', 'Two component system'),
    )
    return [name for key, name in category_names
            if monomer_id in subunit_categories[key]]


def get_tcs_small_molecules(sim_data):
    """
    Get the small molecules shared by all TCS reactions (ATP, ADP, Pi, water and
    H+), i.e. every TCS molecule whose sim_data molecule type is not one of the
    TCS protein types (response regulators, histidine kinases and ligand-bound
    histidine kinases, plus their phosphorylated forms).

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


def build_tcs_reaction_event_estimates(sim_data, all_molecule_changes):
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
        all_molecule_changes: (n_timesteps, n_molecules) TwoComponentSystems
            moleculeChanges (the per-step molecule count changes from the ODE
            integration)
    Returns:
        tuple: ((n_timesteps, n_reactions) estimated reaction counts,
            (n_timesteps, n_molecules) bool mask of molecules whose change
            the estimate misses by more than 0.5)
    """
    sm = sim_data.process.two_component_system.stoich_matrix()
    molecule_names = sim_data.process.two_component_system.molecule_names
    tcs_small_molecules = get_tcs_small_molecules(sim_data)
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
        print(f"Warning: TCS reaction estimate misfits some protein at "
              f"{n_bad}/{n_timesteps} timesteps")

    return estimated_reaction_fluxes, misfit


def build_complexation_dictionaries(sim_data):
    """
    Build the dictionaries describing the complexation complexes each molecule
    is a direct subunit of.

    Args:
        sim_data: Simulation data object
    Returns:
        tuple: (molecules_to_parent_complexes_dict, complex_to_complex_type)
        - molecules_to_parent_complexes_dict: {subunit_id: {complex_id: [
              {'reaction_id': str},
              {'stoichiometry': int}]}}  # Negative value (consumed)
        - complex_to_complex_type: {complex_id: 'homogeneous' or 'heterogeneous'}
    """
    complex_ids = sim_data.process.complexation.ids_complexes
    molecule_names = sim_data.process.complexation.molecule_names
    subunit_names = sim_data.process.complexation.subunit_names
    rxn_ids = sim_data.process.complexation.ids_reactions
    sm = sim_data.process.complexation.stoich_matrix()
    smm = sim_data.process.complexation.stoich_matrix_monomers()

    molecules_to_parent_complexes_dict = {}
    complex_to_complex_type = {}
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
            cplx_idx = complex_ids.index(complex_name)

            # Find the number of unique subunits in this complex reaction:
            num_unique_subunits = len(np.where(smm[:, cplx_idx] < 0)[0])
            complex_to_complex_type[complex_name] = (
                'heterogeneous' if num_unique_subunits > 1 else 'homogeneous')

            # NOTE: this is the subunit's stoichiometry!
            parent_complexes[complex_name] = [
                {'reaction_id': rxn_ids[rxn_idx]},
                {'stoichiometry': sm[subunit_index, rxn_idx]}]

        molecules_to_parent_complexes_dict[subunit] = parent_complexes

    return molecules_to_parent_complexes_dict, complex_to_complex_type


def build_equilibrium_dictionaries(sim_data):
    """
    Build the dictionaries describing the equilibrium complexes each molecule
    is a direct subunit of. Note: complexation complexes can be subunits of
    equilibrium complexes.

    Args:
        sim_data: Simulation data object
    Returns:
        tuple: (molecules_to_parent_complexes_dict, complex_to_complex_type)
        - molecules_to_parent_complexes_dict: {subunit_id: {complex_id: [
              {'reaction_id': str},
              {'stoichiometry': int}]}}  # Negative value (consumed)
        - complex_to_complex_type: {complex_id: 'homogeneous' or 'heterogeneous'}
    """
    monomer_ids = sim_data.process.translation.monomer_data['id']
    complexation_complex_ids = sim_data.process.complexation.ids_complexes
    equilibrium_complex_ids = sim_data.process.equilibrium.ids_complexes
    molecule_names = sim_data.process.equilibrium.molecule_names
    rxn_ids = sim_data.process.equilibrium.rxn_ids
    sm = sim_data.process.equilibrium.stoich_matrix()
    smm = sim_data.process.equilibrium.stoich_matrix_monomers()

    c_molecule_names = sim_data.process.complexation.molecule_names
    c_smm = sim_data.process.complexation.stoich_matrix_monomers()

    molecules_to_parent_complexes_dict = {}
    complex_to_complex_type = {}
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
            cplx_idx = equilibrium_complex_ids.index(complex_name)

            # Find the number of unique monomer subunits in the complex
            # (i.e. not counting metabolite ligands):
            subunit_list = []
            for idx in np.where(smm[:, cplx_idx] < 0)[0]:
                subunit_name = molecule_names[idx]
                if subunit_name in monomer_ids:
                    if subunit_name not in subunit_list:
                        subunit_list.append(subunit_name)
                # Check if the reactant is a complexation complex (metabolite
                # ligands and equilibrium complexes add no subunits here):
                if (subunit_name not in equilibrium_complex_ids
                        and subunit_name in complexation_complex_ids):
                    # Add the unique monomer subunits of the complexation
                    # complex to the list of unique subunits in the eq complex:
                    complexation_complex_index = complexation_complex_ids.index(subunit_name)
                    for c_idx in np.where(c_smm[:, complexation_complex_index] < 0)[0]:
                        c_su_name = c_molecule_names[c_idx]
                        if c_su_name not in subunit_list:
                            subunit_list.append(c_su_name)

            complex_to_complex_type[complex_name] = (
                'heterogeneous' if len(subunit_list) > 1 else 'homogeneous')

            # NOTE: this is the subunit's stoichiometry!
            parent_complexes[complex_name] = [
                {'reaction_id': rxn_ids[rxn_idx]},
                {'stoichiometry': sm[subunit_index, rxn_idx]}]

        molecules_to_parent_complexes_dict[subunit] = parent_complexes

    return molecules_to_parent_complexes_dict, complex_to_complex_type


def build_tcs_dictionaries(sim_data):
    """
    Build the dictionaries describing the phosphorylated TCS complexes each
    unphosphorylated molecule forms, and the molecules each TCS complex
    returns to when it is dephosphorylated. Note: complexation and equilibrium
    complexes can be subunits of TCS complexes.

    For positively oriented TCS systems, a subunit can form the same TCS
    complex via BOTH a ligand-bound HK phosphotransfer reaction AND a
    ligand-free HK phosphotransfer reaction, so the two reaction types are
    kept in separate dictionaries (one dictionary keyed by complex would
    keep only one of the two reactions).

    Args:
        sim_data: Simulation data object
    Returns:
        tuple: (molecules_to_parent_complexes_via_lbptrs,
                molecules_to_parent_complexes_via_ptrs,
                complexes_to_reactants_via_dephosphorylation_reactions,
                complex_to_complex_type)
        - the first three are each {molecule_id: {product_id: [
              {'reaction_id': str},
              {'stoichiometry': int}]}}  # Negative value (consumed)
          for ligand-bound (lbptrs) reactions, all other (ptrs) reactions, and
          reactions where a TCS complex is the reactant
        - complex_to_complex_type: {complex_id: 'PHOSPHO-HK' or 'PHOSPHO-RR'}
    """
    molecule_names = list(sim_data.process.two_component_system.molecule_names)
    complex_ids = list(sim_data.process.two_component_system.complex_to_monomer.keys())
    molecule_types = sim_data.process.two_component_system.molecule_types
    rxn_ids = sim_data.process.two_component_system.rxn_ids
    sm = sim_data.process.two_component_system.stoich_matrix()
    tcs_small_molecules = get_tcs_small_molecules(sim_data)
    pairings = {"RR": "PHOSPHO-RR", "PHOSPHO-RR": "RR", "HK": "PHOSPHO-HK",
                "PHOSPHO-HK": "HK", "HK-LIGAND": "PHOSPHO-HK-LIGAND",
                "PHOSPHO-HK-LIGAND": "HK-LIGAND"}

    molecules_to_parent_complexes_via_lbptrs = {}
    molecules_to_parent_complexes_via_ptrs = {}
    complex_to_complex_type = {}
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
                    product_name = molecule_names[product_index]
                    complex_to_complex_type[product_name] = molecule_types[product_index]
                    # NOTE: this is the subunit's stoichiometry:
                    parent_dict_to_use[product_name] = [
                        {'reaction_id': rxn_id},
                        {'stoichiometry': sm[subunit_index, rxn_idx]}]

        molecules_to_parent_complexes_via_lbptrs[subunit] = parent_complexes_lbptrs
        molecules_to_parent_complexes_via_ptrs[subunit] = parent_complexes_ptrs

    # Dephosphorylation reactions, where the phosphorylated complex (HK or RR)
    # is the reactant:
    complexes_to_reactants_via_dephosphorylation_reactions = {}
    for complex_id in complex_ids:
        complex_index = molecule_names.index(complex_id)
        complement = pairings[molecule_types[complex_index]]

        reaction_indices = np.where(sm[complex_index, :] < 0)[0]
        if len(reaction_indices) == 0:
            continue

        products = {}
        for rxn_idx in reaction_indices:
            for product_index in np.where(sm[:, rxn_idx] > 0)[0]:
                if molecule_types[product_index] == complement:
                    # NOTE: this is the complex's stoichiometry:
                    products[molecule_names[product_index]] = [
                        {'reaction_id': rxn_ids[rxn_idx]},
                        {'stoichiometry': sm[complex_index, rxn_idx]}]
        complexes_to_reactants_via_dephosphorylation_reactions[complex_id] = products

    return (molecules_to_parent_complexes_via_lbptrs,
            molecules_to_parent_complexes_via_ptrs,
            complexes_to_reactants_via_dephosphorylation_reactions,
            complex_to_complex_type)


def format_parent_path_sentence(molecule_label, parent_complex_id, path_info):
    """
    Format a sentence describing how a molecule appears in a parent complex.

    Args:
        molecule_label: What to call the molecule in the sentence (e.g. its ID)
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

    sentence = (f"{molecule_label} counts within {parent_complex_id}, "
                f"\n{abs(total_stoich):g} per (via {path_string})")

    return sentence


class ComplexNetwork:
    """
    The complexation, equilibrium and TCS reactions of one sim_data, as
    dictionaries linking each molecule to the complexes it directly forms:

    - c_m2pc_dict, eq_m2pc_dict, tcs_m2pc_dict: molecule to parent complexes
      for each reaction type (see the build_*_dictionaries functions)
    - tcs_m2pc_via_lbptrs_dict, tcs_m2pc_via_ptrs_dict: tcs_m2pc_dict split by
      ligand-bound and ligand-free phosphotransfer reactions
    - tcs_c2m_via_dprs_dict: TCS complex to the molecules it dephosphorylates to
    - c2ct, eq_c2ct, tcs_c2ct: complex to its makeup or TCS role
    """
    def __init__(self, sim_data):
        self.monomer_ids = sim_data.process.translation.monomer_data['id']
        self.complex_ids = sim_data.process.complexation.ids_complexes
        self.eq_complex_ids = sim_data.process.equilibrium.ids_complexes
        self.tcs_complex_ids = list(
            sim_data.process.two_component_system.complex_to_monomer.keys())

        self.c_m2pc_dict, self.c2ct = build_complexation_dictionaries(sim_data)
        self.eq_m2pc_dict, self.eq_c2ct = build_equilibrium_dictionaries(sim_data)
        (self.tcs_m2pc_via_lbptrs_dict, self.tcs_m2pc_via_ptrs_dict,
         self.tcs_c2m_via_dprs_dict, self.tcs_c2ct) = build_tcs_dictionaries(sim_data)

        # A TCS molecule to parent complexes dictionary for tracing paths, given
        # that intermediate complex stoich is always 1:1:
        self.tcs_m2pc_dict = {}
        for molecule in set(list(self.tcs_m2pc_via_lbptrs_dict.keys()) + list(
                self.tcs_m2pc_via_ptrs_dict.keys())):
            self.tcs_m2pc_dict[molecule] = {}
            if molecule in self.tcs_m2pc_via_lbptrs_dict:
                self.tcs_m2pc_dict[molecule].update(
                    self.tcs_m2pc_via_lbptrs_dict[molecule])
            if molecule in self.tcs_m2pc_via_ptrs_dict:
                self.tcs_m2pc_dict[molecule].update(self.tcs_m2pc_via_ptrs_dict[molecule])

    def determine_molecule_type(self, molecule_id):
        """
        Classify a molecule as monomer, complex, equilibrium complex, or
        two-component system.

        Args:
            molecule_id: The molecule to classify
        Returns:
            str: 'monomer', 'complex', 'eq_complex', 'two component system', or 'unknown'
        """
        if molecule_id in self.monomer_ids:
            return 'monomer'
        elif molecule_id in self.complex_ids:
            return 'complex'
        elif molecule_id in self.eq_complex_ids:
            return 'eq_complex'
        elif molecule_id in self.tcs_complex_ids:
            return 'two component system'
        else:
            return 'unknown'

    def find_downstream_complexes_with_paths(self, molecule_id):
        """
        Find all parent complexes that contain the molecule and trace assembly
        paths. Note: some complexes may have intermediate complexes between the
        base monomers and the final "grandparent" complex. Each distinct path
        to a parent complex is recorded once.

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
        visited_paths = set()  # Path states already explored (to skip repeats)

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

        def record_path(complex_id, path, total_stoich):
            """
            Record a path to a complex if it is not already recorded.
            """
            entries = parent_complexes.setdefault(complex_id, [])
            if not any(entry['path'] == path and entry['total_stoich'] == total_stoich
                       for entry in entries):
                entries.append({'path': path, 'total_stoich': total_stoich})

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
            mol_type = self.determine_molecule_type(current_molecule)

            # Skip a path state (molecule + path so far + stoich) already explored:
            path_signature = (current_molecule, tuple((m[0], m[1]) for m in current_path),
                              current_stoich)
            if path_signature in visited_paths:
                return
            visited_paths.add(path_signature)

            # Get all complexes that directly use this molecule as a subunit
            parents = get_parent_complexes(current_molecule)

            # BASE CASE: No parents found (the molecule is a top-level complex)
            if not parents:
                # Only save if this isn't the starting molecule
                if not is_start:
                    # Complete the path by adding the terminal complex
                    record_path(current_molecule,
                        current_path + [(current_molecule, 1, mol_type, None)],
                        current_stoich)
                return

            # RECURSIVE CASE: Explore each parent of the molecule:
            for parent_id, (parent_type, stoich, reaction_id) in parents.items():
                new_path = current_path + [(current_molecule, stoich, mol_type, reaction_id)]
                new_stoich = current_stoich * stoich

                # Save this parent with its complete path:
                record_path(parent_id, new_path + [(parent_id, 1, parent_type, None)],
                            new_stoich)

                # Continue exploring to see if this parent is itself a subunit of other complexes:
                explore_parents(parent_id, new_path, new_stoich, is_start=False)

        # Verify the starting molecule is valid:
        if self.determine_molecule_type(molecule_id) == 'unknown':
            print(f"Warning: {molecule_id} is not a valid molecule")
            return {}

        # Start exploration from the given molecule
        explore_parents(molecule_id, [], 1, is_start=True)

        return parent_complexes

    def analyze_molecule_parents(self, molecule_id, molecule_label=None):
        """
        Main function to analyze what parent complexes contain a molecule.

        Args:
            molecule_id: The molecule (monomer or complex) to analyze
            molecule_label: What to call the molecule in the path sentences
                (defaults to molecule_id)
        Returns:
            tuple: (filtered_complexes, descriptions)
                - filtered_complexes: Dictionary of parent complexes with their paths
                - descriptions: List of sentences describing each path
        """
        print(f"\nAnalyzing parent complexes for: {molecule_id}\n")
        parent_complexes = self.find_downstream_complexes_with_paths(molecule_id)

        # Filter out cases where the molecule is its own parent (no actual parents):
        filtered_complexes = {k: v for k, v in parent_complexes.items() if k != molecule_id}

        descriptions = []
        for parent_id, path_list in filtered_complexes.items():
            for path_info in path_list:
                sentence = format_parent_path_sentence(
                    molecule_label or molecule_id, parent_id, path_info)
                descriptions.append(sentence)

        if not descriptions:
            print(f"No parent complexes found for {molecule_id}")
        else:
            print(f"Found {len(filtered_complexes)} parent complex(es):\n")
            for desc in descriptions:
                print(desc)
                print()

        return filtered_complexes, descriptions

    def find_extra_tcs_reactions(self, molecule, complex, rxns_plotted):
        """
        Some TCS complexes can be formed from the same reactant by two
        different phosphotransfer reactions (ligand-free and ligand-bound).
        This function checks whether the molecule also forms the parent TCS
        complex via the other reaction type, and if so, returns that reaction
        so its events are plotted too (if they have not already been plotted).

        This matters for events panels, which show every reaction that forms
        or consumes the molecule: if one of the two reactions were skipped, the
        molecule's TCS events would be undercounted.

        Args:
            molecule: ID of the monomer or complex (the reactant)
            complex: ID of the parent TCS complex the molecule forms
            rxns_plotted: reaction IDs already plotted
        Returns:
            list: IDs of the other reactions forming the parent TCS complex
            from the molecule that have not been plotted yet
        """
        reaction_IDs = []

        # Check if the reaction events for the TCS complex via the ligand-free
        # and then the ligand-bound phosphotransfer reactions were plotted:
        for m2pc_dict in (self.tcs_m2pc_via_ptrs_dict, self.tcs_m2pc_via_lbptrs_dict):
            if complex in m2pc_dict.get(molecule, {}):
                alt_rxn_id = m2pc_dict[molecule][complex][0]['reaction_id']
                if alt_rxn_id not in rxns_plotted:
                    reaction_IDs.append(alt_rxn_id)

        return reaction_IDs


def is_a_TF(sim_data, molecule):
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


def tf_index(sim_data, molecule):
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
    return tf_ids.index(molecule)


def is_a_ribosome_subunit(sim_data, molecule):
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


def is_an_RNAP_subunit(sim_data, molecule):
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


def is_a_replisome_subunit(sim_data, molecule):
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


def replisome_subunit_stoich(sim_data, molecule):
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
