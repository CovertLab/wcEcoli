#! /usr/bin/env python

"""
Find valid reactions and metabolites in reaction network.
"""

import pickle
import csv
import os
import sys
from typing import Any, Dict, Iterable, List, Optional, Set, Tuple

from reconstruction.ecoli.dataclasses.process.metabolism import Metabolism
from reconstruction.ecoli.dataclasses.state.external_state import ExternalState
from reconstruction.ecoli.simulation_data import SimulationDataEcoli


# Directories
FILE_LOCATION = os.path.dirname(os.path.realpath(__file__))

# Filenames
METABOLITE_FILE = os.path.join(FILE_LOCATION, 'invalid_metabolites.tsv')
REACTION_FILE = os.path.join(FILE_LOCATION, 'invalid_reactions.tsv')
VALID_MONOMERS_FILE = os.path.join(FILE_LOCATION, 'valid_monomers.tsv')
INVALID_MONOMERS_FILE = os.path.join(FILE_LOCATION, 'invalid_monomers.tsv')


def get_boundaries(metabolism, external_state, media=None):
	# type: (Metabolism, ExternalState, Optional[str]) -> Tuple[Set[str], Set[str]]
	"""
	Get source and sink metabolites. Imports are sources. Concentration targets
	and secretions are sinks.

	Args:
		metabolism: sim_data metabolism process class
		external_state: provides the saved_media and exchange_data_from_media
		media: media label for potential import molecules, if None, uses all
			possible import molecules

	Returns:
		sources: set of metabolite IDs with location tag that are sources
		sinks: set of metabolite IDs with location tag that are sinks
	"""

	if media is None:
		sources = set([
			met for media in external_state.saved_media
			for met in external_state.exchange_data_from_media(media)['importExchangeMolecules']
			])
	else:
		sources = set(external_state.exchange_data_from_media(media)['importExchangeMolecules'])
	sinks = set(external_state.secretion_exchange_molecules)
	sinks.update(metabolism.conc_dict)

	return sources, sinks

def get_mappings(reactions, map_reactants):
	# type: (Dict[str, Dict[str, int]], bool) -> Tuple[Dict[str, List[str]], Dict[str, List[str]]]
	"""
	Mappings of metabolites to reactions and reactions to metabolites for easy
	parsing of all reactions.

	Args:
		reactions: sim_data stoichiometry structure
		map_reactants: True if mapping for reactants, False if mapping for products

	Returns:
		met_to_rxn: map each metabolite (reactant or product) to reactions that
			contain it
		rxn_to_met: map each reaction to metabolites that are either reactants
			or products
	"""

	if map_reactants:
		direction = -1
	else:
		direction = 1

	met_to_rxn = {}  # type: Dict[str, List[str]]
	rxn_to_met = {}  # type: Dict[str, List[str]]
	for rxn, stoich in reactions.items():
		for met, factor in stoich.items():
			if factor * direction > 0:
				met_to_rxn[met] = met_to_rxn.get(met, []) + [rxn]
				rxn_to_met[rxn] = rxn_to_met.get(rxn, []) + [met]

	return met_to_rxn, rxn_to_met

def trace_possible_reactions(start, met_to_rxn, rxn_to_met, excluded_rxns):
	# type: (Set[str], Dict[str, List[str]], Dict[str, List[str]], Set[str]) -> Set[str]
	"""
	Iteratively trace metabolites through valid reactions.

	Args:
		start: starting set of metabolites (sources or sinks)
		met_to_rxn: mapping of metabolites to reactions
		rxn_to_met: mapping of reactions to metabolites
		excluded_rxns: invalid reactions that will not link to new metabolites

	Returns:
		mets: set of metabolites that can be traced from start through valid
			reactions (source -> product or reactant -> sink)
	"""

	mets = set()
	rxns = set(excluded_rxns)

	new_mets = set(start)
	while new_mets:
		# Get reactions that have a source metabolite as an input
		possible_rxns = set()
		for met in new_mets:
			possible_rxns.update({r for r in met_to_rxn.get(met, [])})
		new_rxns = possible_rxns.difference(rxns)

		# Get products from reactions with sources
		possible_mets = set()
		for rxn in new_rxns:
			possible_mets.update({m for m in rxn_to_met.get(rxn, [])})

		# Update metabolites for the next round
		mets.update(new_mets)
		new_mets = possible_mets.difference(mets)

	return mets

def prune_reactions(reactions, valid_mets):
	# type: (Dict[str, Dict[str, int]], Set[str]) -> Set[str]
	"""
	Identify reactions that are not possible because mass balance must exist.

	Args:
		reactions: sim_data stoichiometry structure
		valid_mets: metabolites that have a path to a source and sink
			so that flux can flow through

	Returns:
		excluded_rxns: reactions that are not valid because one or more
			metabolites can not have flux
	"""

	excluded_rxns = set()
	for rxn, stoich in reactions.items():
		if not all([m in valid_mets for m in stoich]):
			excluded_rxns.add(rxn)

	return excluded_rxns

def get_monomers(sim_data, enzymes):
	# type: (SimulationDataEcoli, Set[str]) -> Set[str]
	"""
	Break down complexes into monomers.

	Args:
		sim_data: simulation data
		enzymes: enzyme IDs to break to monomers, must include the location tag

	Returns:
		set of monomers from enzymes and monomers that make up any complexes in enzymes
	"""

	monomers = [
		x for x in enzymes
		if x not in sim_data.process.complexation.complex_names
		and x not in sim_data.process.equilibrium.complex_name_to_rxn_idx
		]
	complexes = [
		x for x in enzymes
		if x in sim_data.process.complexation.complex_names
		or x in sim_data.process.equilibrium.complex_name_to_rxn_idx
		]

	assert len(monomers) + len(complexes) == len(enzymes)

	for complex_ in complexes:
		if complex_ in sim_data.process.complexation.complex_names:
			monomers += sim_data.process.complexation.get_monomers(complex_)['subunitIds'].tolist()
		elif complex_ in sim_data.process.equilibrium.complex_name_to_rxn_idx:
			for subunit in sim_data.process.equilibrium.get_monomers(complex_)['subunitIds'].tolist():
				if subunit in sim_data.process.complexation.complex_names:
					monomers += sim_data.process.complexation.get_monomers(subunit)['subunitIds'].tolist()
				elif subunit in sim_data.process.translation.monomer_data['id']:
					monomers += [subunit]
		else:
			raise ValueError(f'Complex ({complex_}) is not in equilibrium or complexation')

	return {m[:-3] for m in monomers}

def save_to_file(path, data, reactants, products):
	# type: (str, Iterable[Any], Dict[str, List[str]], Dict[str, List[str]]) -> None
	"""Saves data to a file."""

	print('Saving to {}'.format(path))
	with open(path, 'w') as f:
		writer = csv.writer(f, delimiter='\t')
		writer.writerow(['Label', 'Reactants', 'Products'])
		for d in sorted(data):
			writer.writerow([d, reactants.get(d), products.get(d)])

def save_list(path, data):
	# type: (str, Iterable[str]) -> None
	"""Saves data to a file."""

	print('Saving to {}'.format(path))
	with open(path, 'w') as f:
		f.write('\n'.join(sorted(data)))


if __name__ == '__main__':
	# Load data to analyze
	if len(sys.argv) != 2:
		raise RuntimeError('Must pass a path to a sim_data object as an arg to the script.')
	sim_data_file = sys.argv[1]
	print('Loading data from {}'.format(sim_data_file))
	with open(sim_data_file, 'rb') as f:
		sim_data = pickle.load(f)
	metabolism = sim_data.process.metabolism
	external_state = sim_data.external_state
	reactions = metabolism.reaction_stoich

	# Extract necessary data
	print('Analyzing reaction network')
	sources, sinks = get_boundaries(metabolism, external_state)
	reactant_to_rxn, rxn_to_reactant = get_mappings(reactions, True)
	product_to_rxn, rxn_to_product = get_mappings(reactions, False)

	# Find reactions that are not possible due to metabolites not being balanced
	excluded_rxns = set()  # type: Set[str]
	rxn_len = -1  # used for break condition checking
	while len(excluded_rxns) != rxn_len:
		# Store old length to find when no new reactions are found
		rxn_len = len(excluded_rxns)

		# Find metabolites that are source -> product or reactant -> sink
		potential_products = trace_possible_reactions(
			sources, reactant_to_rxn, rxn_to_product, excluded_rxns)
		potential_reactants = trace_possible_reactions(
			sinks, product_to_rxn, rxn_to_reactant, excluded_rxns)

		# Only metabolites with path from a source and to a sink can carry flux and be balanced
		valid_mets = potential_products.intersection(potential_reactants)

		# Exclude reactions that have invalid metabolites
		excluded_rxns = prune_reactions(reactions, valid_mets)

	# Analyze results
	all_mets = set(list(reactant_to_rxn) + list(product_to_rxn))
	all_rxns = set(list(rxn_to_reactant) + list(rxn_to_product))
	excluded_mets = all_mets.difference(valid_mets)
	valid_rxns = all_rxns.difference(excluded_rxns)
	all_enzymes = {cat for cat in metabolism.catalyst_ids}
	valid_enzymes = {
		enz
		for rxn in valid_rxns
		for enz in metabolism.reaction_catalysts.get(rxn, [])
		}
	all_monomers = get_monomers(sim_data, all_enzymes)
	valid_monomers = get_monomers(sim_data, valid_enzymes)
	invalid_monomers = all_monomers - valid_monomers

	# Print summary
	print('\t{}/{} metabolites are valid'.format(len(valid_mets), len(all_mets)))
	print('\t{}/{} reactions are valid'.format(len(valid_rxns), len(all_rxns)))
	print('\t{}/{} enzymes are valid'.format(len(valid_enzymes), len(all_enzymes)))
	print('\t{}/{} enzyme monomers are valid'.format(len(valid_monomers), len(all_monomers)))

	# Save output to files
	save_to_file(METABOLITE_FILE, excluded_mets, reactant_to_rxn, product_to_rxn)
	save_to_file(REACTION_FILE, excluded_rxns, rxn_to_reactant, rxn_to_product)
	save_list(VALID_MONOMERS_FILE, valid_monomers)
	save_list(INVALID_MONOMERS_FILE, invalid_monomers)
