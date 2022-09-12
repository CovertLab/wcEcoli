"""
SimulationData moleculeGroups
"""

from __future__ import absolute_import, division, print_function

POLYMERIZED_FRAGMENT_PREFIX = 'polymerized_'

class MoleculeGroups(object):
	"""
	Helper class to extract molecule IDs of "special" groups of molecules. All
	values returned are lists of strings.
	"""

	def __init__(self, raw_data, sim_data):
		self._build_molecule_groups(raw_data, sim_data)

	def _build_molecule_groups(self, raw_data, sim_data):
		aa_ids = list(sim_data.amino_acid_code_to_id_ordered.values())
		ntp_ids = list(sim_data.ntp_code_to_id_ordered.values())
		dntp_ids = list(sim_data.dntp_code_to_id_ordered.values())
		polymerized_aa_ids = [
			POLYMERIZED_FRAGMENT_PREFIX + aa_id for aa_id in aa_ids]
		polymerized_ntp_ids = [
			POLYMERIZED_FRAGMENT_PREFIX + ntp_id for ntp_id in ntp_ids]
		polymerized_dntp_ids = [
			POLYMERIZED_FRAGMENT_PREFIX + dntp_id for dntp_id in dntp_ids]

		# Build list of rRNA IDs from raw data
		s30_16s_rRNA = [
			rna['id'] + '[c]' for rna in raw_data.rnas
			if rna['id'].startswith('RRS')]
		s50_23s_rRNA = [
			rna['id'] + '[c]' for rna in raw_data.rnas
			if rna['id'].startswith('RRL')]
		s50_5s_rRNA = [
			rna['id'] + '[c]' for rna in raw_data.rnas
			if rna['id'].startswith('RRF')]

		# Build list of ribosomal proteins from raw data
		monomer_ids = set([monomer['id'] for monomer in raw_data.proteins])

		complex_id_to_subunit_ids = {}
		for rxn in raw_data.complexation_reactions:
			complex_ids = []
			subunit_ids = []

			for mol, v in rxn['stoichiometry'].items():
				if (v is None or v < 0):
					subunit_ids.append(mol)
				else:
					complex_ids.append(mol)

			assert len(complex_ids) == 1

			complex_id_to_subunit_ids[complex_ids[0]] = subunit_ids

		def find_protein_subunits(molecule_id):
			"""
			Recursive function to find all protein monomers that are subunits
			of a given molecule.
			"""
			subunits = []

			if molecule_id in monomer_ids:
				subunits.append(molecule_id)
			elif molecule_id in complex_id_to_subunit_ids:
				for mol in complex_id_to_subunit_ids[molecule_id]:
					subunits.extend(find_protein_subunits(mol))

			return subunits

		s30_proteins = [
			mol + '[c]' for mol
			in find_protein_subunits(sim_data.molecule_ids.s30_full_complex[:-3])
			]
		s50_proteins = [
			mol + '[c]' for mol
			in find_protein_subunits(sim_data.molecule_ids.s50_full_complex[:-3])
			]

		assert len(set(s30_proteins) & set(s50_proteins)) == 0
		ribosomal_proteins = sorted(s30_proteins + s50_proteins)

		# Build list of RNA polymerase subunits from raw data
		RNAP_subunits = [
			mol + '[c]' for mol
			in find_protein_subunits(sim_data.molecule_ids.full_RNAP[:-3])
			]

		molecule_groups = {
			'amino_acids': aa_ids,
			'ntps': ntp_ids,
			'dntps': dntp_ids,

			'polymerized_amino_acids': polymerized_aa_ids,
			'polymerized_ntps': polymerized_ntp_ids,
			'polymerized_dntps': polymerized_dntp_ids,
			'polymerized_subunits': polymerized_aa_ids + polymerized_ntp_ids + polymerized_dntp_ids,

			's30_proteins':	s30_proteins,
			's30_16s_rRNA': s30_16s_rRNA,

			's50_protein_complexes': ['CPLX0-3956[c]'],
			's50_proteins':	s50_proteins,
			's50_23s_rRNA': s50_23s_rRNA,
			's50_5s_rRNA': s50_5s_rRNA,

			'lipids': ['CPD-8260[c]', 'CPD-12819[c]', 'CPD-12824[c]'],
			'polyamines': ['GAMMA-GLUTAMYL-PUTRESCINE[c]', 'PUTRESCINE[c]',
				'GLUTATHIONYLSPERMIDINE[c]', 'SPERMIDINE[c]',
				'N1-ACETYLSPERMINE[c]', 'SPERMINE[c]'],

			# TODO: 'EG10245-MONOMER[c]' (DNAP III subunit tau) should be added
			# 	to the list of trimer subunits once frame-shifting proteins are
			# 	produced.
			'replisome_trimer_subunits': ['CPLX0-2361[c]', 'CPLX0-3761[c]'],
			'replisome_monomer_subunits': ['CPLX0-3621[c]', 'EG10239-MONOMER[c]',
				'EG11500-MONOMER[c]', 'EG11412-MONOMER[c]'],

			'exoRNases': ['EG11620-MONOMER[c]', 'G7175-MONOMER[c]',
				'EG10858-MONOMER[c]', 'EG10863-MONOMER[c]', 'EG11259-MONOMER[c]',
				'EG11547-MONOMER[c]', 'EG10746-MONOMER[c]', 'G7842-MONOMER[c]',
				'EG10743-MONOMER[c]'],
			'endoRNase_rnas': ['EG10856_RNA', 'EG10857_RNA', 'EG10859_RNA',
				'EG10860_RNA', 'EG10861_RNA', 'EG10862_RNA', 'EG11299_RNA',
				'G7175_RNA', 'G7365_RNA'],
			'exoRNase_rnas': ['EG11620_RNA', 'G7175_RNA', 'EG10858_RNA',
				'EG10863_RNA', 'EG11259_RNA', 'EG11547_RNA', 'EG10746_RNA',
				'G7842_RNA', 'EG10743_RNA'],
			'RNAP_subunits': RNAP_subunits,

			'ribosomal_proteins': ribosomal_proteins,

			'carbon_sources': ['GLC[p]', 'ACET[p]', 'SUC[p]'],
		}

		# Initialize molecule groups for how molecules are split between two
		# daughter cells at cell division (populated later by InternalState)
		molecule_groups['bulk_molecules_binomial_division'] = []
		molecule_groups['bulk_molecules_equal_division'] = []

		molecule_groups['unique_molecules_active_ribosome_division'] = []
		molecule_groups['unique_molecules_RNA_division'] = []
		molecule_groups['unique_molecules_domain_index_division'] = []
		molecule_groups['unique_molecules_chromosomal_segment_division'] = []

		self.__dict__.update(molecule_groups)
