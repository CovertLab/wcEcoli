"""
This plot allows one to visualize the counts of a specified complex over time,
along with the counts of the free monomers that make up that complex,
and the complexation events that produce that complex.

# NOTES: as of 03/01/2026:
- 81 complexation complexes are subunits of other complexation complexes
- 16 complexation complexes are subunits of equilibrium complexes
- 4 complexation complexes are subunits of two component system complexes
- 1 equilibrium complex is a subunit of another equilibrium complex
- 4 equilibrium complexes are subunits of two component system complexes

In total, there are 1131 unique complexation complexes, 39 unique equilibrium
complexes, and 18 unique two component system complexes.

"""
import pickle
import os
from scipy.optimize import curve_fit
from matplotlib import pyplot as plt
# noinspection PyUnresolvedReferences
import numpy as np

from models.ecoli.analysis import multigenAnalysisPlot
from wholecell.analysis.analysis_tools import (exportFigure,
    read_bulk_molecule_counts, read_stacked_bulk_molecules, read_stacked_columns)
from wholecell.io.tablereader import TableReader
import wholecell.utils.units as units

PLOT_COMPLEXES = ['PC00027'] #
                  #"MONOMER0-160", "MONOMER0-155", ]


# TODO LIST:
# todo: clean up the extract_doubling_times function
# todo: determine if the is valid molecule function should be moved outside the class
# todo: see if it is possible to build the molecules_to_all_downstream_complexes_dict
#  and molecules_to_parent_complexes_dict in this file using just SM and SMM
# doubling time function from nora (note the normal doubling time extraction is not working):
# TODO: clean up this function:
def extract_doubling_times(cell_paths):
    # Load simulation time span data
    time = read_stacked_columns(cell_paths, 'Main', 'time').squeeze()

    # Determine doubling time
    doubling_times = read_stacked_columns(cell_paths, 'Main', 'time',
                                                  fun=lambda x: (x[-1] - x[0])).squeeze().astype(
                int)
    end_generation_times = np.cumsum(doubling_times) + time[0]  #
    start_generation_indices = np.searchsorted(time, end_generation_times[:-1],
                                                       side='left').astype(int)
    start_generation_indices = np.insert(start_generation_indices, 0, 0) + np.arange(
                len(doubling_times))
    end_generation_indices = start_generation_indices + doubling_times

    return time, doubling_times, end_generation_times, start_generation_indices, end_generation_indices


class Plot(multigenAnalysisPlot.MultigenAnalysisPlot):

    # Get the validity:
    def check_validity_and_get_compartment(self, sim_data, molecule_list):
        revised_molecule_list = []
        for molecule in molecule_list:
            if "[" in molecule:
                molecule = molecule[:-3]  # remove compartment
            if sim_data.getter.is_valid_molecule(molecule):
                revised_name = molecule + sim_data.getter.get_compartment_tag(molecule)
                revised_molecule_list.append(revised_name)

        return revised_molecule_list

    # Function to check if an input is indeed a complex
    def check_complex_validity(self, molecule_list):
         # The the valid complexs are complexation complexes, equilibrium complexes, and two component system complexes:
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
                      f"equiliibrium, and two component system complexes can "
                      f"be plotted here.")
        return valid_complexes, complex_types


    # Build the complexation parent and downstream molecule to complex dictionaries:
    def build_complexation_dictionaries(self, sim_data):
        molecule_names = sim_data.process.complexation.molecule_names
        subunit_names = sim_data.process.complexation.subunit_names # TODO: note: eq does NOT have subunit names!
        rxn_ids = sim_data.process.complexation.ids_reactions
        rsu = sim_data.process.complexation.reaction_stoichiometry_unknown
        smI = sim_data.process.complexation._stoich_matrix_I # TODO: see if just the SM can be used here!
        smJ = sim_data.process.complexation._stoich_matrix_J
        smV = sim_data.process.complexation._stoich_matrix_V
        smm = sim_data.process.complexation.stoich_matrix_monomers()

        molecules_to_parent_complexes_dict = {}
        molecules_to_all_downstream_complexes_dict = {}

        # Generate dictionary mapping molecules to direct parent complexes they form:
        for subunit in subunit_names:
            # Find the matrix index where this subunit is as a molecule:
            subunit_index = molecule_names.index(subunit)

            # Find the indicies of smJ where the value will correspond to the
            # index of reaction in rxn_ids the subunit is a reactant in:
            reaction_indicies = np.where(
                (smI == subunit_index) &
                (smV < 0))[0]

            # For each reaction index, find the complex(es) that is(are) formed:
            parent_complexes = {}
            for reaction_idx in reaction_indicies:
                # Find the value of index that corresponds to the index of the reaction in rxn_ids:
                rxn_idx = smJ[reaction_idx]

                # Initialize variables to hold info about the complex:
                complex_information = []
                stoich = {} # NOTE: this will be the subunit's stoichiometry needed to make the complex!
                stoich_known = {}
                complex_type = {}
                reaction_name = {}

                # Find the complex formed in this reaction
                # (NOTE: each reaction forms one complex):
                complex_index = np.where(
                    (smJ == rxn_idx) &
                    (smV > 0))[0]

                # Find the number of unique subunits in this complex reaction:
                unique_subunits_in_complex = np.where(
                    (smJ == rxn_idx) &
                    (smV < 0))[0]
                num_unique_subunits = len(unique_subunits_in_complex)
                if num_unique_subunits > 1:
                    cplx_type = 'heterogeneous'
                else:
                    cplx_type = 'homogeneous'

                # Add complex information to lists
                complex_name = molecule_names[smI[complex_index][0]]
                reaction_name['reaction_id'] = rxn_ids[rxn_idx]
                stoich['stoichiometry'] = smV[reaction_idx]
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
                reaction_indices = np.where(
                    (smI == cplx_idx) & (smV > 0))[0]

                # Obtain the value that corresponds to the index of the reaction:
                reaction_idx = smJ[reaction_indices]

                # Initialize variables to hold info about the complex:
                downstream_complex_information = []
                stoich = {}
                stoich_known = {}
                complex_type = {}
                reaction_name = {}

                # Find the number of unique subunits in this complex reaction:
                unique_subunits_in_complex = np.where(
                    (smJ == reaction_idx) &
                    (smV < 0))[0]
                num_unique_subunits = len(unique_subunits_in_complex)
                if num_unique_subunits > 1:
                    cplx_type = 'heterogeneous'
                elif num_unique_subunits == 1:
                    cplx_type = 'homogeneous'
                else:
                    cplx_type = 'unknown' # TODO why is this here

                # Add complex information to lists:
                reaction_name['reaction_id'] = rxn_ids[reaction_idx[0]]
                stoich['stoichiometry'] = smm[subunit_index, complex_idx]
                stoich_known['stoich_unknown'] = rsu[
                    reaction_idx[0]]
                complex_type['complex_type'] = cplx_type
                downstream_complex_information.append(reaction_name)
                downstream_complex_information.append(stoich)
                downstream_complex_information.append(stoich_known)
                downstream_complex_information.append(complex_type)

                # Append the complex name and stoich as a dictionary entry
                downstream_complexes[complex_name] = downstream_complex_information

            molecules_to_all_downstream_complexes_dict[subunit] = downstream_complexes

        return molecules_to_parent_complexes_dict, molecules_to_all_downstream_complexes_dict


    def build_equilibrium_dictionaries(self, sim_data):
        molecule_names = sim_data.process.equilibrium.molecule_names
        rxn_ids = sim_data.process.equilibrium.rxn_ids
        rsu = sim_data.process.equilibrium.reaction_stoichiometry_unknown
        smI = sim_data.process.equilibrium._stoichMatrixI  # TODO: see if just the SM can be used here!
        smJ = sim_data.process.equilibrium._stoichMatrixJ
        smV = sim_data.process.equilibrium._stoichMatrixV
        smm = sim_data.process.equilibrium.stoich_matrix_monomers()

        molecules_to_parent_complexes_dict = {}
        molecules_to_all_downstream_complexes_dict = {}
        complex_to_complex_type = {}

        # Generate dictionary mapping equilibrium molecules to the direct parent
        # complexes they form:
        for subunit in molecule_names:
            # Find the matrix index where this subunit is as a molecule:
            subunit_index = molecule_names.index(subunit)

            # Find the indicies of smJ where the value will
            # correspond to the index of the correct reaction in rxn_ids:
            reaction_indicies = np.where(
                (smI == subunit_index) &
                (smV < 0))[0]

            if len(reaction_indicies) == 0:
                # skip molecules that do not form any complexes
                continue

            # For each reaction index, find the complex(es) that is(are) formed:
            parent_complexes = {}
            for reaction_idx in reaction_indicies:
                # Find the value of index that corresponds to the index of the reaction in rxn_ids:
                rxn_idx = smJ[reaction_idx]

                # Initialize variables to hold info about the complex:
                complex_information = []
                stoich = {} # NOTE: this will be the subunit's stoichiometry needed to make the complex!
                stoich_known = {}
                complex_type = {}
                reaction_name = {}

                # Find the complex formed in the reaction (each reaction forms
                # one equilibrium complex):
                complex_index = np.where(
                    (smJ == rxn_idx) &
                    (smV > 0))[0]

                # Find the number of unique monomer subunits in the complex
                # (i.e. not counting metabolite ligands):
                unique_subunits_in_complex = np.where(
                    (smJ == rxn_idx) &
                    (smV < 0))[0]
                subunit_list = []
                for idx in unique_subunits_in_complex:
                    subunit_name = molecule_names[smI[idx]]
                    if subunit_name in self.monomerIDs:
                        if subunit_name not in subunit_list:
                            subunit_list.append(subunit_name)
                    # Check the reactant is a complexation complex:
                    elif subunit_name in self.complexIDs:
                        # See how many unique monomer subunits the complex has:
                        complexation_complex_index = self.complexIDs.index(subunit_name)
                        c_smm = sim_data.process.complexation.stoich_matrix_monomers()
                        complexation_unique_subunits = np.where(c_smm[:, complexation_complex_index] < 0)[0]
                        c_molecule_names = sim_data.process.complexation.molecule_names
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

                # Add info about the complex to the parent complex info lists:
                complex_name = molecule_names[smI[complex_index][0]]
                complex_to_complex_type[complex_name] = cplx_type
                reaction_name['reaction_id'] = rxn_ids[rxn_idx]
                stoich['stoichiometry'] = smV[reaction_idx]
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
            # Find the matrix index where this subunit is as a molecule:
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

                # Use the stoichMatrix() to find the reaction index:
                reaction_indices = np.where(
                    (smI == cplx_idx) & (smV > 0))[0]

                # Obtain the value that corresponds to the index of the reaction:
                reaction_idx = smJ[reaction_indices]

                # Initialize variables to hold info about the complex:
                downstream_complex_information = []
                stoich = {}
                stoich_known = {}
                complex_type = {}
                reaction_name = {}

                # Add complex information to lists:
                reaction_name['reaction_id'] = rxn_ids[reaction_idx[0]]
                stoich['stoichiometry'] = smm[subunit_index, complex_idx]
                stoich_known['stoich_unknown'] = rsu[
                    reaction_idx[0]]
                complex_type['complex_type'] = complex_to_complex_type[complex_name]
                downstream_complex_information.append(reaction_name)
                downstream_complex_information.append(stoich)
                downstream_complex_information.append(stoich_known)
                downstream_complex_information.append(complex_type)

                # Append the complex name and stoich as a dictionary entry
                downstream_complexes[complex_name] = downstream_complex_information

            molecules_to_all_downstream_complexes_dict[subunit] = downstream_complexes

        return molecules_to_parent_complexes_dict, molecules_to_all_downstream_complexes_dict

    def build_tcs_dictionaries(self, sim_data):
        molecule_names = list(sim_data.process.two_component_system.molecule_names)
        modified_molecules = list(sim_data.process.two_component_system.modified_molecules)
        molecule_types = sim_data.process.two_component_system.molecule_types
        rxn_ids = sim_data.process.two_component_system.rxn_ids
        smI = sim_data.process.two_component_system._stoichMatrixI  # TODO: see if just the SM can be used here!
        smJ = sim_data.process.two_component_system._stoichMatrixJ
        smV = sim_data.process.two_component_system._stoichMatrixV
        smm = sim_data.process.two_component_system.stoich_matrix_monomers_TEMP(sim_data)
        complex_ids = list(sim_data.process.two_component_system.complex_to_monomer.keys())
        molecules_to_skip = ["ATP[c]", "ADP[c]", "Pi[c]", "WATER[c]", "PROTON[c]"]
        pairings = {"RR": "PHOSPHO-RR", "PHOSPHO-RR": "RR", "HK": "PHOSPHO-HK",
                    "PHOSPHO-HK": "HK", "HK-LIGAND": "PHOSPHO-HK-LIGAND",
                    "PHOSPHO-HK-LIGAND": "HK-LIGAND"}

        molecules_to_parent_complexes_dict = {}
        molecules_to_all_downstream_complexes_dict = {}

        # Generate mapping from molecules to direct parent complexes they form:
        for subunit in molecule_names:
            # If the molecule is not a complex or monomer, skip it:
            if subunit in molecules_to_skip:
                continue

            # Find the matrix index where this subunit is as a molecule:
            subunit_index = molecule_names.index(subunit)

            reactant_type = molecule_types[molecule_names.index(subunit)]
            compliment = pairings[reactant_type]

            # Find the indicies of smJ where the value will
            # correspond to the index of the correct reaction in rxn_ids:
            reaction_indicies = np.where(
                (smI == subunit_index) &
                (smV < 0))[0]

            # For each reaction index, find the product(s) that is(are) formed:
            parent_complexes = {}
            for reaction_idx in reaction_indicies:
                # Find the value of the reaction index (which will be the index
                # that corresponds to the index of the reaction in self.ids_reactions):
                rxn_idx = smJ[reaction_idx]

                # Initialize data structures to hold product information
                product_information = []
                stoich = {} # NOTE: this will be the subunit's stoichiometry needed to make the complex!
                stoich_known = {}
                product_type = {}
                reaction_name = {}

                # Find the products formed in this reaction:
                product_indices = np.where(
                    (smJ == rxn_idx) &
                    (smV > 0))[0]

                # TCS reactions can have multiple products, so they need to be filtered:
                for product_index in product_indices:
                    product_name = molecule_names[smI[product_index]]

                    # Since there are multiple products, determine which product
                    # is the compliment of the specific reactant subunit:
                    pdt_type = molecule_types[molecule_names.index(product_name)]
                    if pdt_type != compliment:
                        continue
                    else:
                        # Add product's information to lists:
                        reaction_name['reaction_id'] = rxn_ids[rxn_idx]
                        stoich['stoichiometry'] = smV[reaction_idx]
                        stoich_known['stoich_unknown'] = 'not applicable'
                        product_type['complex_type'] = pdt_type
                        product_information.append(reaction_name)
                        product_information.append(stoich)
                        product_information.append(stoich_known)
                        product_information.append(product_type)

                        # Append the product name and stoich as a dictionary entry
                        parent_complexes[product_name] = product_information

                molecules_to_parent_complexes_dict[subunit] = parent_complexes

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
                # skip molecules that are not downstream subunits of any TCS complexes
                continue

            downstream_complexes = {}
            for complex_idx in complex_indices:
                # Find the complex's name:
                complex_name = complex_ids[complex_idx]

                # Obtain the index of the complex within self.molecule_names:
                cplx_idx = modified_molecules.index(complex_name)

                # Use the stoichMatrix() to find the reaction index:
                reaction_indices = np.where(
                    (smI == cplx_idx) & (smV > 0))[0]

                # Obtain the value that corresponds to the index of the reaction:
                reaction_idx = smJ[reaction_indices]

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
                complex_type['complex_type'] = molecule_types[molecule_names.index(complex_name)]
                downstream_complex_information.append(reaction_name)
                downstream_complex_information.append(stoich)
                downstream_complex_information.append(stoich_known)
                downstream_complex_information.append(complex_type)

                # Append the complex name and stoich as a dictionary entry
                downstream_complexes[complex_name] = downstream_complex_information

            molecules_to_all_downstream_complexes_dict[subunit] = downstream_complexes

        return molecules_to_parent_complexes_dict, molecules_to_all_downstream_complexes_dict


    # UNIQUE TYPES: TFs, ribosomes, RNAps, replisomes

    def is_a_TF(self, sim_data, molecule):
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
        # Get the TF IDs:
        tfs = sim_data.process.transcription_regulation.tf_ids
        tf_ids = [tf_id + f'[{sim_data.getter.get_compartment(tf_id)[0]}]'
                          for tf_id in tfs]

        # Find the index of the TF in the TF list:
        tf_index = tf_ids.index(molecule)

        return tf_index

    def is_a_ribosome_subunit(self, sim_data, molecule):
        # Get the ribosome subunit IDs:
        ribosome_50s_subunits = sim_data.process.complexation.get_monomers(
            sim_data.molecule_ids.s50_full_complex)
        ribosome_30s_subunits = sim_data.process.complexation.get_monomers(
            sim_data.molecule_ids.s30_full_complex)
        ribosome_subunit_ids = (ribosome_50s_subunits["subunitIds"].tolist() +
                                ribosome_30s_subunits["subunitIds"].tolist())

        # Determine if the complex is a ribosome subunit:
        ribosome_subunit_id_match = [subunit_id for subunit_id in ribosome_subunit_ids if molecule == subunit_id]

        if ribosome_subunit_id_match != []:
            print(f"{molecule} is a ribosome subunit.")

        return ribosome_subunit_id_match

    def ribosome_subunit_index(self, sim_data, molecule):
        # Get the ribosome subunit IDs:
        ribosome_50s_subunits = sim_data.process.complexation.get_monomers(
            sim_data.molecule_ids.s50_full_complex)
        ribosome_30s_subunits = sim_data.process.complexation.get_monomers(
            sim_data.molecule_ids.s30_full_complex)
        ribosome_subunit_ids = (ribosome_50s_subunits["subunitIds"].tolist() +
                                ribosome_30s_subunits["subunitIds"].tolist())

        # Find the index of the ribosome subunit in the ribosome subunit list:
        ribosome_subunit_index = ribosome_subunit_ids.index(molecule)

        return ribosome_subunit_index

    def is_an_RNAP_subunit(self, sim_data, molecule):
        # Get the RNAP subunit IDs:
        rnap_subunits = sim_data.process.complexation.get_monomers(
            sim_data.molecule_ids.full_RNAP)
        rnap_subunit_ids = rnap_subunits["subunitIds"].tolist()

        # Determine if the complex is an RNAP subunit:
        rnap_subunit_id_match = [subunit_id for subunit_id in rnap_subunit_ids if molecule == subunit_id]

        if rnap_subunit_id_match != []:
            print(f"{molecule} is an RNAP subunit.")

        return rnap_subunit_id_match

    def rnap_subunit_index(self, sim_data, molecule):
        # Get the RNAP subunit IDs:
        rnap_subunits = sim_data.process.complexation.get_monomers(
            sim_data.molecule_ids.full_RNAP)
        rnap_subunit_ids = rnap_subunits["subunitIds"].tolist()

        # Find the index of the RNAP subunit in the RNAP subunit list:
        rnap_subunit_index = rnap_subunit_ids.index(molecule)

        return rnap_subunit_index

    def is_a_replisome_subunit(self, sim_data, molecule):
        # Get the replisome subunit IDs:
        replisome_trimer_subunits = sim_data.molecule_groups.replisome_trimer_subunits
        replisome_monomer_subunits = sim_data.molecule_groups.replisome_monomer_subunits
        replisome_subunit_ids = replisome_trimer_subunits + replisome_monomer_subunits

        # Determine if the complex is a replisome subunit:
        replisome_subunit_id_match = [subunit_id for subunit_id in replisome_subunit_ids if molecule == subunit_id]

        if replisome_subunit_id_match != []:
            print(f"{molecule} is a replisome subunit.")

        return replisome_subunit_id_match

    def replisome_subunit_index(self, sim_data, molecule):
        # Get the replisome subunit IDs:
        replisome_trimer_subunits = sim_data.molecule_groups.replisome_trimer_subunits
        replisome_monomer_subunits = sim_data.molecule_groups.replisome_monomer_subunits
        replisome_subunit_ids = replisome_trimer_subunits + replisome_monomer_subunits

        # Find the index of the replisome subunit in the replisome subunit list:
        replisome_subunit_index = replisome_subunit_ids.index(molecule)

        return replisome_subunit_index









    def do_plot(self, seedOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
        with open(simDataFile, 'rb') as f:
            sim_data = pickle.load(f)
        cell_paths = self.ap.get_cells()
        sim_dir = cell_paths[0]
        simOutDir = os.path.join(sim_dir, 'simOut')

        # extract the doubling times:
        time, doubling_times, end_generation_times, start_generation_indices, end_generation_indices = extract_doubling_times(
            cell_paths)

        # Extract protein indexes for each monomer
        monomer_counts_reader = TableReader(
            os.path.join(simOutDir, "MonomerCounts"))
        monomer_idx_dict = {monomer: i for i, monomer in
                            enumerate(monomer_counts_reader.readAttribute(
                                'monomerIds'))}
        # Extract the monomer IDs:
        self.monomerIDs = monomer_counts_reader.readAttribute("monomerIds")

        # Extract the free monomer counts using the monomer counts listener:
        free_monomer_counts = read_stacked_columns(cell_paths, 'MonomerCounts', "freeMonomerCounts")


        # Extract the three types of complexes:
        # TODO: delete these readers!
        complex_counts_listener = TableReader(
            os.path.join(simOutDir, "ComplexationListener"))
        self.complexIDs = sim_data.process.complexation.ids_complexes
        complex_idx_dict = {complexID: i for i, complexID in
                            enumerate(self.complexIDs)}
        eq_complex_counts_listener = TableReader(
            os.path.join(simOutDir, "EquilibriumListener"))
        self.eqComplexIDs = sim_data.process.equilibrium.ids_complexes
        eq_complex_idx_dict = {complexID: i for i, complexID in
                            enumerate(self.eqComplexIDs)}
        self.tcsComplexIDs = list(sim_data.process.two_component_system.complex_to_monomer.keys())
        tcs_molecule_ids = list(sim_data.process.two_component_system.molecule_names)
        tcs_complex_idx_dict = {complexID: i for i, complexID in
                            enumerate(tcs_molecule_ids)}

        self.c_m2pc_dict, self.c_m2adc_dict = self.build_complexation_dictionaries(sim_data)

        # Warning: these will not unpack subunits to monomers if the subunit is a complexation complex!
        self.eq_m2pc_dict, self.eq_m2adc_dict = self.build_equilibrium_dictionaries(sim_data)

        # Warning: these will not unpack subunits to baseline monomers if the subunit is a complexation or equilibrium complex!
        self.tcs_m2pc_dict, self.tcs_m2adc_dict = self.build_tcs_dictionaries(sim_data)


        hi = 5


        # Make sure the proteins inputted are valid and have a compartment tag:
        PLOT_COMPLEXES_revised = self.check_validity_and_get_compartment(sim_data, PLOT_COMPLEXES)

        # Check that all the valid inputs are indeed complexes:
        valid_complexes, complex_type_dict = self.check_complex_validity(PLOT_COMPLEXES_revised)

        # for each complex, we need to
        # 1. determine what type of complex it is (complexation, equilibrium, or two component system)
        # 2. find the monomers that make up the complex and the stoichiometry of the complex
        # 3. if it is a complexation complex, determine if it is a subunit of another complexatoin complex, or an EQ or 2CS complex
        # 4. if it is an EQ complex, determine if it is a subunit of a 2CS complex. also determine if its subunits are a complexation complex, and if so, find its monomers that way.
        # 5. determine if it is a subunit of a unique molecule complex (i.e. ribosome, RNAP, replisome, or TF)


        # Generate a plot for each complex:
        for complex in valid_complexes:
            all_base_monomers = {}
            reaction_to_idx_dict = {}
            complexation_complex_subunits = []
            equilibrium_complex_subunits = []
            complexation_parent_complexes = []
            equilibirum_parent_complexes = []
            tcs_parent_complexes = []

            # Unique molecule subunit tracking:
            tf_existance = self.is_a_TF(sim_data, complex)
            ribosome_existance = self.is_a_ribosome_subunit(sim_data, complex)
            rnap_existance = self.is_an_RNAP_subunit(sim_data, complex)
            replisome_existance = self.is_a_replisome_subunit(sim_data, complex)


            # Determine the complex type:
            complex_type = complex_type_dict[complex]

            if complex_type == "complexation":
                complex_idx = complex_idx_dict[complex]
                complex_counts = read_stacked_columns(cell_paths, 'ComplexationListener',
                                                      "complexCounts")[:, complex_idx]

                # Check if the complex is a subunit of another complexation complex itself:
                if complex in self.c_m2pc_dict.keys():
                    complexation_parent_complexes.append(self.c_m2pc_dict[complex])

                # Find all constiuent monomers of the complex based on where the
                # complex shows up as a value in molecules_to_all_downstream_complexes dict:
                monomers = {k: v for k, v in self.c_m2adc_dict.items() if
                            complex in v}

                # Find which complexation reactions the complex is involved in:
                complex_reactions = complex_counts_listener.readAttribute("reactionIDs") # see if we can get rid of the listener by directly indexing to the reaction ids in sim data?
                complex_makeup = '' # TODO: check that this is where we want this
                for monomer in monomers:
                    monomer_info = monomers[monomer]
                    monomer_complex_info = monomer_info[complex][0]
                    reaction_id = monomer_complex_info['reaction_id']
                    complex_makeup = monomer_info[complex][3]['complex_type']
                    all_base_monomers[monomer] = monomer_info
                    if reaction_id in complex_reactions:
                        reaction_idx = complex_reactions.index(reaction_id)
                        reaction_to_idx_dict[reaction_id] = reaction_idx

                # Check if the complex is a subunit of an equilibrium complex:
                if complex in self.eq_m2pc_dict.keys():
                    equilibirum_parent_complexes.append(self.eq_m2pc_dict[complex])

                # Check if the complex is a subunit of a TCS complex:
                if complex in self.tcs_m2pc_dict.keys():
                    tcs_parent_complexes.append(self.tcs_m2pc_dict[complex])

            elif complex_type == "equilibrium":
                complex_idx = eq_complex_idx_dict[complex]
                complex_counts = read_stacked_columns(cell_paths, 'EquilibriumListener',
                                                      "complexCounts")[:, complex_idx]

                # Check if the complex is a subunit of another EQ complex itself:
                if complex in self.eq_m2pc_dict.keys():
                    equilibirum_parent_complexes.append(self.eq_m2pc_dict[complex])

                # Find all base molecules of the complex based on where the
                # complex shows up as a value in molecules_to_all_downstream_complexes dict:
                molecules = {k: v for k, v in self.eq_m2adc_dict.items() if
                            complex in v}

                # Find all base monomers of the equilibrium complex:
                subunits = []
                for molecule in molecules:
                    if molecule in self.monomerIDs:
                        all_base_monomers[molecule] = molecules[molecule]
                        subunits.append(molecule)
                    elif molecule in self.complexIDs:
                        complexation_complex_subunits.append(molecule)
                        subunits.append(molecule)
                    else:
                        # if the molecule is a ligand, skip it
                        # (i.e. do not add it to the subunits list)
                        continue

                # Find which EQ reactions this complex is involved in:
                eq_complex_reactions = eq_complex_counts_listener.readAttribute("reactionIDs")
                complex_makeup = ''
                for subunit in subunits:
                    subunit_info = molecules[subunit]
                    subunit_complex_info = subunit_info[complex][0]
                    reaction_id = subunit_complex_info['reaction_id']
                    complex_makeup = subunit_info[complex][3]['complex_type']
                    if reaction_id in eq_complex_reactions:
                        reaction_idx = eq_complex_reactions.index(reaction_id)
                        reaction_to_idx_dict[reaction_id] = reaction_idx

                # Check if the complex is a subunit of a TCS complex:
                if complex in self.tcs_m2pc_dict.keys():
                    tcs_parent_complexes.append(self.tcs_m2pc_dict[complex])

                # TODO: add checks for TFs, ribosome subunits, RNAP subunits, and replisome subunits here as well (and add to the plot if so)

            # TCS complexes:
            else:
                complex_idx = self.tcsComplexIDs.index(complex)
                complex_counts = read_stacked_columns(cell_paths, 'MonomerCounts',
                                                      "twoComponentSystemComplexCounts")[:, complex_idx]

                # Find all constiuent monomers of the complex based on where the
                # complex shows up as a value in molecules_to_all_downstream_complexes dict
                # (NOTE: the tcs_m2adc dict is similar to c_m2adc because it
                # does actually map each complex down to its monomers, whereas
                # the eq_m2adc dict does not, hence why it needed extra unpacking
                # of its subunits to affirm their identity)
                monomers = {k: v for k, v in self.tcs_m2adc_dict.items() if
                            complex in v}

                # Check if any of the direct subunits are complexation or EQ complexes:
                molecules = {k: v for k, v in self.tcs_m2pc_dict.items() if
                             complex in v}
                for molecule in molecules:
                    if molecule in self.complexIDs:
                        complexation_complex_subunits.append(molecule)
                    elif molecule in self.eqComplexIDs:
                        equilibrium_complex_subunits.append(molecule)


                # Find which TCS reactions the complex is involved in:
                tcs_complex_reactions = sim_data.process.two_component_system.rxn_ids
                complex_makeup = ''
                for monomer in monomers:
                    monomer_info = monomers[monomer]
                    monomer_complex_info = monomer_info[complex][0]
                    all_base_monomers[monomer] = monomer_info
                    reaction_id = monomer_complex_info['reaction_id']
                    complex_makeup = monomer_info[complex][3]['complex_type']
                    if reaction_id in tcs_complex_reactions:
                        reaction_idx = tcs_complex_reactions.index(reaction_id)
                        reaction_to_idx_dict[reaction_id] = reaction_idx


                # TODO: check if the complex is a subunit of a unique molecule type!


            # Generate the plots:
            fig, (ax1, ax2, ax3) = plt.subplots(3, 1, sharex=True, figsize=(10, 6))

            # First, make a plot of the complex itself and its subunits (as well as the complexes it forms)
            ax1.plot(time, complex_counts, color='lightseagreen', label=f'{complex} free complex counts', linewidth=.75, alpha=0.75)

            # Plot the counts of each monomer within the complex:

            # OK so if I have an eq complex being plotted, and one of its subunits is a complexation complex, I want it to say MONOMERX counts within COMPLEX (Y per, via Z COMPLEXATION COMPLEX)
            # same with TCS complexes!

            """monomers_to_add = {k: v for k, v in self.c_m2adc_dict.items() if
                               molecule in v}
            for monomer in monomers_to_add:
                if monomer not in all_base_monomers:
                    all_base_monomers[monomer] = monomers_to_add[monomer][molecule][0]"""

            # ACTUALLY, how I want to do this is have the monomer subunits literally just be the monomers from that complexes actual reaction. then use the complexation or equilibirum subunits [] to get the monomers within!

            if complex_type == "complexation":
                molecules_to_parent_complexes = self.c_m2pc_dict
                molecules_to_all_downstream_complexes = self.c_m2adc_dict
            elif complex_type == "equilibrium":
                molecules_to_parent_complexes = self.eq_m2pc_dict
                molecules_to_all_downstream_complexes = self.eq_m2adc_dict
            else:
                molecules_to_parent_complexes = self.tcs_m2pc_dict
                molecules_to_all_downstream_complexes = self.tcs_m2adc_dict

            all_monomer_subunits = list(all_base_monomers.keys())
            for monomer in all_base_monomers.keys():
                # get the counts of the monomer in the complex:
                monomer_info = all_base_monomers[monomer][complex][1]
                monomer_stoich = monomer_info['stoichiometry']  * -1 # use 1 to get the stoichiometry dict
                monomer_complex_counts = complex_counts * monomer_stoich

                # check if the monomer is the result of the complex or a different parent complex:
                monomer_parent_complexes = molecules_to_parent_complexes[monomer]
                if complex in monomer_parent_complexes.keys():
                    # if the complex is a parent complex of the monomer, plot with a different style:
                    ax1.plot(time, monomer_complex_counts, alpha=1, label=f'{monomer} counts within {complex} \n({monomer_stoich} per)', linewidth=1, linestyle=":")
                else:
                    for pc in monomer_parent_complexes.keys():
                        gparent_complexes = molecules_to_parent_complexes[pc]
                        if complex in gparent_complexes.keys():
                            # if the complex is a parent complex of the monomer's parent complex, note that:
                            ax1.plot(time, monomer_complex_counts, alpha=1, label=f'{monomer} counts within {complex} \n({monomer_stoich} per, via {pc})', linewidth=.8, linestyle="-.")
                            break # Do not continue checking the other grandparents!
                        else:
                            # if the complex is not a parent complex of the monomer or its parent complexes:
                            ax1.plot(time, monomer_complex_counts, alpha=1,
                                     label=f'{monomer} counts within\n{complex} ({monomer_stoich} per, via an unknown intermediate)',
                                     linewidth=.8, linestyle="-.")

            # Also check if any of the complex's subunits are complexation complexes that need to be further unpacked:
            if complexation_complex_subunits != []:
                for subunit in complexation_complex_subunits:

                    # Find the counts of the subunit within the complex first:
                    subunit_info = molecules_to_all_downstream_complexes[subunit][complex][1]
                    subunit_stoich = subunit_info['stoichiometry'] * -1

                    # Find where the subunit shows up in the complexation complexes:
                    monomers = {k: v for k, v in self.c_m2adc_dict.items() if
                                       subunit in v}

                    for monomer in monomers:
                        all_monomer_subunits.append(monomer)
                        monomer_info = monomers[monomer][subunit][1]
                        monomer_stoich = monomer_info['stoichiometry'] * -1
                        total_stoich = subunit_stoich * int(monomer_stoich) # TODO check this!
                        monomer_complex_counts = complex_counts * total_stoich
                        ax1.plot(time, monomer_complex_counts, alpha=1,
                                 label=f'{monomer} counts within {complex} \n({total_stoich} per, via {subunit}, {subunit_stoich} per)', linewidth=.8, linestyle=":")

            if equilibrium_complex_subunits != []:
                for subunit in equilibrium_complex_subunits:
                    # Find the counts of the subunit within the complex first:
                    subunit_info = molecules_to_all_downstream_complexes[subunit][complex][0]
                    subunit_stoich = subunit_info['stoichiometry'] * -1

                    # Find where the subunit shows up in the equilibrium complexes:
                    molecules = {k: v for k, v in self.eq_m2adc_dict.items() if
                                       subunit in v}

                    # check if the monomers are within a complexation complex too:
                    temp_complexation_complex_subunits = []
                    monomers = {}
                    for molecule in molecules:
                        if molecule in self.complexIDs:
                            temp_complexation_complex_subunits.append(molecule)
                        if molecule in self.monomerIDs:
                            monomers[molecule] = molecules[molecule]

                    for monomer in monomers:
                        all_monomer_subunits.append(monomer)
                        monomer_info = monomers[monomer][subunit][1]
                        monomer_stoich = monomer_info['stoichiometry'] * -1
                        total_stoich = subunit_stoich * int(monomer_stoich) # TODO check this!
                        monomer_complex_counts = complex_counts * total_stoich
                        ax1.plot(time, monomer_complex_counts, alpha=1,
                                 label=f'{monomer} counts within {complex} \n({total_stoich} per, via {subunit}, {subunit_stoich} per)', linewidth=.8, linestyle=":")

                    # Unpack any complexation complex subunits within the equilibrium complex subunit as well:
                    if temp_complexation_complex_subunits != []:
                        for cplx_subunit in temp_complexation_complex_subunits:
                            # Find the counts of the subunit within the complex first:
                            cplx_subunit_info = molecules_to_all_downstream_complexes[cplx_subunit][subunit][1]
                            cplx_subunit_stoich = cplx_subunit_info['stoichiometry'] * -1

                            # Find where the subunit shows up in the complexation complexes:
                            temp_monomers = {k: v for k, v in self.c_m2adc_dict.items() if
                                               cplx_subunit in v}

                            for monomer in temp_monomers:
                                all_monomer_subunits.append(monomer)
                                monomer_info = temp_monomers[monomer][cplx_subunit][1]
                                monomer_stoich = monomer_info['stoichiometry'] * -1
                                # subunit stoich = how many of the eq subunit are in the complex, complexation subunit stoich = how many of the complexation subunit are in the eq subunit, monomer stoich = how many of the monomer are in the complexation subunit
                                total_stoich = subunit_stoich * int(cplx_subunit_stoich) * int(monomer_stoich)
                                monomer_complex_counts = complex_counts * total_stoich
                                ax1.plot(time, monomer_complex_counts, alpha=1,
                                         label=f'{monomer} counts within {complex} '
                                               f'\n({total_stoich} per, via {subunit}, {subunit_stoich} per)', linewidth=.8, linestyle=":")



            # Next, plot the free monomer counts for each protein that makes up the complex:
            for monomer in all_monomer_subunits:
                protein_idx = monomer_idx_dict[monomer]
                protein_FMC = free_monomer_counts[:, protein_idx]
                ax2.plot(time, protein_FMC, alpha=0.75, label=f'{monomer}', linewidth=0.75)

            ax2.set_ylabel(f"Free Monomer \nSubunit Availability")
            ax2.legend(fontsize=5, loc="center left", bbox_to_anchor=(1, 0.5))

            # Third, add a plot of the complexation events over time:
            for rxn in reaction_to_idx_dict.keys():
                rxn_idx = reaction_to_idx_dict[rxn]
                if complex_type == "complexation":
                    rxn_events = read_stacked_columns(cell_paths, 'ComplexationListener', "complexationEvents")[:, rxn_idx]
                elif complex_type == "equilibrium":
                    # NOTE: equilibrium reactions can go negative (reverse reactions), so the events here may be negative
                    rxn_events = read_stacked_columns(cell_paths, 'EquilibriumListener', "complexationEvents")[:, rxn_idx]
                else: # TCS reactions
                    # Find this by indexing to where the complex is made within the TCS delta molecules counter, since they are all one-to-one reactions:
                    tcs_index = tcs_complex_idx_dict[complex]
                    rxn_events = read_stacked_columns(cell_paths, 'MonomerCounts', "delta2CMolecules")[:, tcs_index]
                ax3.plot(time, rxn_events, alpha=0.75, label=f'{rxn} \n(generates 1 {complex})', linewidth=0.75)


            # If the complex is itself a subunit of a larger complex, plot the
            # counts of the parent complex and the events that generate it as well:
            if complexation_parent_complexes != []:
                parent_complexes = complexation_parent_complexes[0].keys()
                complex_reactions = sim_data.process.complexation.ids_reactions
                for parent_complex in parent_complexes:
                    # Find the counts of complex within the parent complex:
                    parent_complex_idx = complex_idx_dict[parent_complex]
                    parent_complex_counts = read_stacked_columns(cell_paths,
                                                                     'ComplexationListener',
                                                                     "complexCounts")[:,
                                                parent_complex_idx]

                    complex_stoich = self.c_m2pc_dict[complex][parent_complex][1]['stoichiometry']
                    complex_counts_within_parent_complex = parent_complex_counts * complex_stoich * -1

                    # Find the events that generate the parent complex:
                    parent_complex_reaction_ID = self.c_m2pc_dict[complex][parent_complex][0]['reaction_id']
                    reaction_idx = complex_reactions.index(parent_complex_reaction_ID)
                    parent_complex_events = read_stacked_columns(cell_paths,
                                                                     'ComplexationListener',
                                                                     "complexationEvents")[:,
                                                reaction_idx]
                    ax1.plot(time, complex_counts_within_parent_complex, alpha=0.5,
                             label=f'{complex} counts within\n{parent_complex} ({complex_stoich * -1} per)',
                             linewidth=0.75,
                             linestyle="-.")
                    ax3.plot(time, np.negative(parent_complex_events), alpha=0.5,
                             label=f'{parent_complex_reaction_ID} \n(consumes {complex_stoich * -1} {complex} to generate {parent_complex})',
                             linewidth=0.75)

            if equilibirum_parent_complexes != []:
                parent_complexes = equilibirum_parent_complexes[0].keys()
                complex_reactions = eq_complex_counts_listener.readAttribute("reactionIDs")
                for parent_complex in parent_complexes:
                    # Find the counts of complex within the parent complex:
                    parent_complex_idx = eq_complex_idx_dict[parent_complex]
                    parent_complex_counts = read_stacked_columns(cell_paths,
                                                                 'EquilibriumListener',
                                                                 "complexCounts")[:,
                                            parent_complex_idx]

                    complex_stoich = self.eq_m2pc_dict[complex][parent_complex][1]['stoichiometry']
                    complex_counts_within_parent_complex = parent_complex_counts * complex_stoich * -1

                    # Find the events that generate the parent complex:
                    parent_complex_reaction_ID = self.eq_m2pc_dict[complex][parent_complex][0][
                        'reaction_id']
                    reaction_idx = complex_reactions.index(parent_complex_reaction_ID)
                    parent_complex_events = read_stacked_columns(cell_paths,
                                                                 'EquilibriumListener',
                                                                 "complexationEvents")[:,
                                            reaction_idx]
                    ax1.plot(time, complex_counts_within_parent_complex, alpha=0.5,
                             label=f'{complex} counts within\n{parent_complex} ({complex_stoich * -1} per)',
                             linewidth=0.75,
                             linestyle="-.")
                    ax3.plot(time, np.negative(parent_complex_events), alpha=0.5,
                             label=f'{parent_complex_reaction_ID} \n(consumes {complex_stoich * -1} {complex} to generate {parent_complex})',
                             linewidth=0.75)

            if tcs_parent_complexes != []:
                parent_complexes = tcs_parent_complexes[0].keys()
                for parent_complex in parent_complexes:
                    # Find the counts of complex within the parent complex:
                    parent_complex_idx = tcs_complex_idx_dict[parent_complex]
                    parent_complex_counts = read_stacked_columns(cell_paths,
                                                                 'MonomerCounts',
                                                                 "twoComponentSystemComplexCounts")[:,
                                            parent_complex_idx]

                    complex_stoich = self.tcs_m2pc_dict[complex][parent_complex][1]['stoichiometry']
                    complex_counts_within_parent_complex = parent_complex_counts * complex_stoich * -1

                    # Find the events that generate the parent complex:
                    parent_complex_reaction_ID = self.tcs_m2pc_dict[complex][parent_complex][0][
                        'reaction_id']
                    reaction_idx = tcs_molecule_ids.index(parent_complex)
                    parent_complex_events = read_stacked_columns(cell_paths,
                                                                 'MonomerCounts',
                                                                 "delta2CMolecules")[:,
                                            reaction_idx]
                    ax1.plot(time, complex_counts_within_parent_complex, alpha=0.5,
                             label=f'{complex} counts within\n{parent_complex} ({complex_stoich * -1} per)',
                             linewidth=0.75,
                             linestyle="-.")
                    ax3.plot(time, np.negative(parent_complex_events), alpha=0.5,
                             label=f'{parent_complex_reaction_ID} \n(consumes {complex_stoich * -1} {complex} to generate {parent_complex})',
                             linewidth=0.75)

            # Check if the complex is a TF:
            if tf_existance != []:
                tf_idx = self.tf_index(sim_data, complex)
                tf_counts = read_stacked_columns(cell_paths,
                                                                 'RnaSynthProb',
                                                                 "nActualBound")[:,
                                            tf_idx]
                tfs_unbound = read_stacked_columns(cell_paths,
                                                                 'RnaSynthProb',
                                                                 "nActualUnbound")[:,
                                            tf_idx]
                ax1.plot(time, tf_counts, color='lightcoral', label=f'transcription unit bound\n {complex} transcription factors', linewidth=0.75, alpha=0.75)

                ax3.plot(time, (tf_counts * -1), color='lightcoral', label=f'Binding events of transcription factor \n{complex}', linewidth=0.75, alpha=0.75)
                ax3.plot(time, tfs_unbound, color='indianred', label=f'Unbinding events of \n{complex} transcription factor', linewidth=0.75, alpha=0.75)

            # TODO: decide whether or not to delete this based off whether there are any complex subunits of ribosomes.
            if ribosome_existance != []:
                ribosome_subunit_idx = self.ribosome_subunit_index(sim_data, complex) # todo: does this matter? arent they all 1:1?
                ribosome_subunit_counts = 5 # TODO: find ths by summing monomersElongated over all time steps. that is the number of ribosomes that get used!
                ax1.plot(time, ribosome_subunit_counts, color='mediumpurple', label=f'Ribosome subunit \n{complex} counts', linewidth=0.75, alpha=0.75)
                ax3.plot(time, ribosome_subunit_counts * -1, color='mediumpurple', label=f'Ribosome assembly events of subunit \n{complex}', linewidth=0.75, alpha=0.75)


            # TODO: add RNAPs and replisomes!




            # Set the axes and legends (since some might have been added with the parent complexes):
            ax1.legend(fontsize=5, loc="center left", bbox_to_anchor=(1, 0.5))
            ax1.set_ylabel("Complex Counts")
            ax1.set_title(
                f"Complex counts for the {complex_makeup} {complex_type} complex {complex}\n Sim ID: {metadata['description']}")
            ax3.set_xlabel("Time (s)")
            ax3.set_ylabel("Complexation Events")
            ax3.legend(fontsize=5, loc="center left", bbox_to_anchor=(1, 0.5))

            # Add vertical lines for the end of each generation:
            for i in range(len(end_generation_times)):
                dt = end_generation_times[i]
                a = 0.5
                ax1.axvline(x=dt, linestyle='--', color="yellowgreen", alpha=a)
                ax2.axvline(x=dt, linestyle='--', color="yellowgreen", alpha=a)
                ax3.axvline(x=dt, linestyle='--', color="yellowgreen", alpha=a)


            plt.tight_layout()

            # Save the plot:
            file_name = plotOutFileName + "_" +complex
            exportFigure(plt, plotOutDir, file_name, metadata)








if __name__ == '__main__':
    Plot().cli()
