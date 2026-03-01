"""
This plot allows one to visualize the counts of a specified complex over time,
along with the counts of the free monomers that make up that complex,
and the complexation events that produce that complex.
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

PLOT_COMPLEXES = ["MONOMER0-160",] #
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
                    complex_types[molecule] = "tcs"
            else:
                print(f"{molecule} is not a valid complex. Only complexation, "
                      f"equiliibrium, and two component system complexes can "
                      f"be plotted here.")
        return valid_complexes, complex_types


    # Build the complexation parent and downstream molecule to complex dictionaries:
    def build_complexation_dictionaries(self, sim_data):
        # NOTE: 81 complexation complexes are complexes of other complexation complexes.
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

        # Generate dictionary mapping molecules to the direct parent complexes they form:
        for subunit in subunit_names:
            # find the matrix index where this subunit is as a molecule:
            subunit_index = molecule_names.index(subunit)

            # Find the indicies of self._stoich_matrix_J where the value will
            # correspond to the index of the correct reaction in self.ids_reactions
            reaction_indicies = np.where(
                (smI == subunit_index) &
                (smV < 0))[0]

            # For each reaction index, find the complex(es) that is(are) formed:
            parent_complexes = {}
            for reaction_idx in reaction_indicies:
                # find the value of the reaction index (which will be the index
                # that corresponds to the index of the reaction in self.ids_reactions):
                rxn_idx = smJ[reaction_idx]

                # Initialize data structures to hold complex information
                complex_information = []
                stoich = {}
                stoich_known = {}
                complex_type = {}
                reaction_name = {}

                # Find the complex formed in this reaction (each reaction forms one complex):
                complex_index = np.where(
                    (smJ == rxn_idx) &
                    (smV > 0))[0]

                # Find the number of unique subunits in this complex reaction:
                unique_subunits_in_complex = np.where(
                    (smJ == rxn_idx) &
                    (smV < 0))[0]
                num_unique_subunits = len(unique_subunits_in_complex)
                # TODO: for the EQ complex, edit this to not count the metabolite ids. then check if the remaining are in monomerIDs.
                #  if they are not in monomerIDs, they are probably a complexation complex. check if that complexation complex is unique
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

        # Make a dictionary mapping molecules to all downstream complexes they form
        # (both directly and indirectly via another complex):
        for subunit in subunit_names:
            # find the matrix index where this subunit is as a molecule:
            subunit_index = molecule_names.index(subunit)

            # Find the indices of complexes containing the subunit:
            complex_indices = np.where(smm[subunit_index, :] < 0)[0]

            downstream_complexes = {}
            for complex_idx in complex_indices:
                # Find the complex's name:
                complex_name = self.complexIDs[complex_idx] # TODO: double check that complexIDs is the exact same as sim_data.process.complexation.ids_complexes

                # Obtain the index of the complex within self.molecule_names:
                cplx_idx = molecule_names.index(complex_name)

                # Use the stoichMatrix() to find the reaction index:
                reaction_indices = np.where(
                    (smI == cplx_idx) & (smV > 0))[0]

                # Obtain the value that corresponds to the index of the reaction in self.ids_reactions:
                reaction_idx = smJ[reaction_indices]

                # Initialize data structures to hold complex information
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

    # NOTE: complexation complexes can be subunits of other complexation complexes, equilibrium complexes, and two component system complexes
    def is_a_complexation_complex_subunit(self, sim_data, molecule):
        # Get the complexation complex subunit IDs:
        complexation_complex_ids = sim_data.process.complexation.complexation_complex_ids

        # Determine if the complex is a complexation complex subunit:
        complexation_complex_id_match = [complex_id for complex_id in complexation_complex_ids if molecule == complex_id]

        if complexation_complex_id_match != []:
            print(f"{molecule} is a complexation complex subunit.")

        return complexation_complex_id_match








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
        monomerIDs = monomer_counts_reader.readAttribute("monomerIds")

        # Extract the free monomer counts using the monomer counts listener:
        free_monomer_counts = read_stacked_columns(cell_paths, 'MonomerCounts', "freeMonomerCounts")


        # Extract the three types of complexes:
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
        tcs_complex_idx_dict = {complexID: i for i, complexID in
                            enumerate(self.tcsComplexIDs)}

        self.c_m2pc_dict, self.c_m2adc_dict = self.build_complexation_dictionaries(sim_data)



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
            monomers = {}
            complexation_parent_complexes = {}
            equilibirum_parent_complexes = {}
            tcs_parent_complexes = {}

            # Determine the complex type:
            complex_type = complex_type_dict[complex]

            if complex_type == "complexation":
                complex_idx = complex_idx_dict[complex]
                complex_counts = read_stacked_columns(cell_paths, 'ComplexationListener',
                                                      "complexCounts")[:, complex_idx]
                # find the complex and its consituent monomers:
                # TODO: see if it is possible to generate this dictionary in this file!
                molecules_to_all_downstream_complexes = sim_data.process.complexation.molecules_to_all_downstream_complexes_dict
                # determine if the complex itself is a subunit of a larger complex:
                # TODO: figure out how to search for this in any complex type! Maybe make one large dictionary that combines both complexation and equilibrium?
                molecules_to_parent_complexes = sim_data.process.complexation.molecules_to_parent_complexes_dict

                # If the complex is a subunit of another complexation complex,
                # it will show up in the molecules_to_parent_complexes dictionary:
                # TODO: note: for eq and tcs, the complexes will show up in the
                #  keys no matter what bc there are no subunit_names to make the dicts with.
                #  instead do parent_complexes = molecules_to_parent_complexes[complex] and then check if its not empty and add it if so.
                if complex in molecules_to_parent_complexes.keys():
                    complexation_parent_complexes = molecules_to_parent_complexes[complex]


                # Find all constiuent monomers of the complex based on where the
                # complex shows up as a value in MDSC dict:
                monomers = {k: v for k, v in molecules_to_all_downstream_complexes.items() if
                            complex in v}

                # Find which complexation reactions this complex is involved in:
                complex_reactions = complex_counts_listener.readAttribute("reactionIDs")
                complexation_rxn_to_idx_dict = {}
                complex_makeup = '' # TODO: check that this is where we want this
                for monomer in monomers:
                    monomer_info = monomers[monomer]
                    monomer_complex_info = monomer_info[complex][0]
                    reaction_id = monomer_complex_info['reaction_id']
                    complex_makeup = monomer_info[complex][3]['complex_type']
                    if reaction_id in complex_reactions:
                        reaction_idx = complex_reactions.index(reaction_id)
                        complexation_rxn_to_idx_dict[reaction_id] = reaction_idx

                # TODO: check if the complex is part of an equilibrium complex or tcs as well, using both itself AND its parent complex if applicable!

            # BEFORE:


            elif complex in complex_idx_dict: # TODO: need to use keys here, otherwise the complexation complexes part of eq complexes will be recognized.
                complex_type = "complexation"
                complex_idx = complex_idx_dict[complex]
                # find the complex and its consituent monomers:
                molecules_to_all_downstream_complexes = sim_data.process.complexation.molecules_to_all_downstream_complexes_dict
                # determine if the complex itself is a subunit of a larger complex:
                # TODO: figure out how to search for this in any complex type! Maybe make one large dictionary that combines both complexation and equilibrium?
                molecules_to_parent_complexes = sim_data.process.complexation.molecules_to_parent_complexes_dict
                if complex in molecules_to_parent_complexes.keys(): # TODO: double check this does not include the grandparent complexes in it? I feel like maybe it could? and that could cuase issues?
                    parent_complexes = molecules_to_parent_complexes[complex]
                # find where the complex shows up as a value in the dict:
                monomers = {k: v for k, v in molecules_to_all_downstream_complexes.items() if complex in v}
                complex_counts = read_stacked_columns(cell_paths, 'ComplexationListener', "complexCounts")[:, complex_idx]

                # TODO: determine if any complexes are generated by multiple reactions (if so, probalby need to generate a dictionary)
                # Find which reactions this complex is involved in:
                complex_reactions = complex_counts_listener.readAttribute("reactionIDs")
                reaction_to_idx_dict = {}
                complex_makeup = ''
                for monomer in monomers:
                    monomer_info = monomers[monomer]
                    monomer_complex_info = monomer_info[complex][0]
                    reaction_id = monomer_complex_info['reaction_id']
                    complex_makeup = monomer_info[complex][3]['complex_type']
                    if reaction_id in complex_reactions:
                        reaction_idx = complex_reactions.index(reaction_id)
                        reaction_to_idx_dict[reaction_id] = reaction_idx

            elif complex in eq_complex_idx_dict:
                complex_type = "equilibrium"
                complex_idx = eq_complex_idx_dict[complex]
                # find the complex and its consituent monomers:
                molecules_to_all_downstream_complexes = sim_data.process.equilibrium.molecules_to_all_downstream_complexes_dict
                # determine if the complex itself is a subunit of a larger complex: # TODO: figure out how to search for this in any complex type! Maybe make one large dictionary that combines both complexation and equilibrium?
                molecules_to_parent_complexes = sim_data.process.equilibrium.molecules_to_parent_complexes_dict
                if complex in molecules_to_parent_complexes.keys():
                    parent_complexes = molecules_to_parent_complexes[complex]
                # find where the complex shows up as a value in the dict:
                monomers_temp = {k: v for k, v in molecules_to_all_downstream_complexes.items() if
                            complex in v}
                # remove any monomers that are not actually monomers (i.e. ATP):
                for key in monomers_temp.keys():
                    if key in monomerIDs:
                        monomers[key] = monomers_temp[key]

                complex_counts = read_stacked_columns(cell_paths, 'EquilibriumListener',
                                                      "complexCounts")[:, complex_idx]

                # find which reactions this complex is involved in:
                complex_reactions = eq_complex_counts_listener.readAttribute("reactionIDs")
                reaction_to_idx_dict = {}
                complex_makeup = ''
                for monomer in monomers:
                    monomer_info = monomers[monomer]
                    monomer_complex_info = monomer_info[complex][0]
                    reaction_id = monomer_complex_info['reaction_id']
                    complex_makeup = monomer_info[complex][3]['complex_type']
                    if reaction_id in complex_reactions:
                        reaction_idx = complex_reactions.index(reaction_id)
                        reaction_to_idx_dict[reaction_id] = reaction_idx

            else:
                # TODO: determine if any complexes are in both complexation and equilibrium (if so, need to handle that here)?
                raise ValueError(f"Complex {complex} not found in complexation complex or equilibrium complex  listeners.")

            # Generate the plots:
            fig, (ax1, ax2, ax3) = plt.subplots(3, 1, sharex=True, figsize=(10, 6))

            # First, make a plot of the complex itself and its subunits (as well as the complexes it forms)
            ax1.plot(time, complex_counts, color='lightseagreen', label=f'{complex} free complex counts', linewidth=.75, alpha=0.75)

            # plot the counts of each monomer that makes up the complex:
            for monomer in monomers.keys():
                # get the counts of the monomer in the complex:
                monomer_info = monomers[monomer]
                monomer_complex_info = monomer_info[complex][1] # use 1 to get the stoichiometry dict
                monomer_stoich = monomer_complex_info['stoichiometry']
                monomer_complex_counts = complex_counts * monomer_stoich
                # check if the monomer is the result of the complex or a different parent complex:
                monomer_parent_complexes = molecules_to_parent_complexes[monomer]
                if complex in monomer_parent_complexes.keys():
                    # if the complex is a parent complex of the monomer, plot with a different style:
                    ax1.plot(time, monomer_complex_counts, alpha=1, label=f'{monomer} counts within\n{complex} ({monomer_stoich} per)', linewidth=1, linestyle=":")
                else:
                    for pc in monomer_parent_complexes.keys():
                        gparent_complexes = molecules_to_parent_complexes[pc]
                        if complex in gparent_complexes.keys():
                            # if the complex is a parent complex of the monomer's parent complex, note that:
                            ax1.plot(time, monomer_complex_counts, alpha=1, label=f'{monomer} counts within\n{complex} ({monomer_stoich} per, via {pc})', linewidth=.8, linestyle="-.")
                            break
                        else:
                            # if the complex is not a parent complex of the monomer or its parent complexes:
                            ax1.plot(time, monomer_complex_counts, alpha=1,
                                     label=f'{monomer} counts within\n{complex} ({monomer_stoich} per, via an unknown intermediate)',
                                     linewidth=.8, linestyle="-.")

            # Next, plot the free monomer counts for each protein that makes up the complex:
            for monomer in monomers.keys():
                protein_idx = monomer_idx_dict[monomer]
                protein_FMC = free_monomer_counts[:, protein_idx]
                ax2.plot(time, protein_FMC, alpha=0.75, label=f'{monomer}', linewidth=0.75)

            ax2.set_ylabel(f"Free Monomer \nSubunit Counts")
            ax2.legend(fontsize=5, loc="center left", bbox_to_anchor=(1, 0.5))

            # Third, add a plot of the complexation events over time:
            for rxn in reaction_to_idx_dict.keys():
                rxn_idx = reaction_to_idx_dict[rxn]
                if complex_type == "equilibrium":
                    # NOTE: equilibrium reactions can go negative (reverse reactions), so the events here may be negative
                    rxn_events = read_stacked_columns(cell_paths, 'EquilibriumListener', "complexationEvents")[:, rxn_idx]
                else:
                    rxn_events = read_stacked_columns(cell_paths, 'ComplexationListener', "complexationEvents")[:, rxn_idx]
                ax3.plot(time, rxn_events, alpha=0.75, label=f'{rxn} \n(generates 1 {complex})', linewidth=0.75)


            # If the complex is itself a subunit of a larger complex, plot the
            # counts of the parent complex and the events that generate it as well:
            if parent_complexes != {}:
                for parent_complex in parent_complexes:
                    if complex_type == "complexation":
                        # find the counts of complex within the parent complex:
                        parent_complex_idx = complex_idx_dict[parent_complex]
                        parent_complex_counts = read_stacked_columns(cell_paths,
                                                                     'ComplexationListener',
                                                                     "complexCounts")[:,
                                                parent_complex_idx]
                        complex_stoich = parent_complexes[parent_complex][1]['stoichiometry']
                        complex_counts_within_parent_complex = parent_complex_counts * complex_stoich

                        # find the events that generate the parent complex:
                        parent_complex_reaction_ID = parent_complexes[parent_complex][0]['reaction_id']
                        reaction_idx = complex_reactions.index(parent_complex_reaction_ID)
                        parent_complex_events = read_stacked_columns(cell_paths,
                                                                     'ComplexationListener',
                                                                     "complexationEvents")[:,
                                                reaction_idx]
                    else:
                        parent_complex_idx = eq_complex_idx_dict[parent_complex]
                        parent_complex_counts = read_stacked_columns(cell_paths,
                                                                     'EquilibriumListener',
                                                                     "complexCounts")[:,
                                                parent_complex_idx]
                        complex_stoich = parent_complexes[parent_complex][1]['stoichiometry']
                        complex_counts_within_parent_complex = parent_complex_counts * complex_stoich

                        # find the events that generate the parent complex:
                        parent_complex_reaction_ID = parent_complexes[parent_complex][0][
                            'reaction_id']
                        reaction_idx = complex_reactions.index(parent_complex_reaction_ID)
                        parent_complex_events = read_stacked_columns(cell_paths,
                                                                     'ComplexationListener',
                                                                     "complexationEvents")[:,
                                                reaction_idx]


                    ax1.plot(time, complex_counts_within_parent_complex, alpha=0.5,
                             label=f'{complex} counts within\n{parent_complex} ({complex_stoich} per)', linewidth=0.75,
                             linestyle="-.")
                    ax3.plot(time, np.negative(parent_complex_events), alpha=0.5,
                             label=f'{parent_complex_reaction_ID} \n(consumes {complex_stoich} {complex} to generate {parent_complex})',
                             linewidth=0.75)


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
