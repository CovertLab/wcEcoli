"""
MonomerCounts Listener
"""

import numpy as np
import wholecell.listeners.listener


class MonomerCounts(wholecell.listeners.listener.Listener):
	"""
	Listener for the full counts of protein monomers, including those that are
	part of a complex.
	"""
	_name = 'MonomerCounts'

	def __init__(self, *args, **kwargs):
		super(MonomerCounts, self).__init__(*args, **kwargs)

	def initialize(self, sim, sim_data):
		super(MonomerCounts, self).initialize(sim, sim_data)

		# Get IDs of all bulk molecules
		self.bulkMolecules = sim.internal_states["BulkMolecules"]
		bulk_molecule_ids = self.bulkMolecules.container.objectNames()

		# Get IDs of molecules involved in complexation and equilibrium
		complexation_molecule_ids = sim_data.process.complexation.molecule_names
		complexation_complex_ids = sim_data.process.complexation.ids_complexes
		equilibrium_molecule_ids = sim_data.process.equilibrium.molecule_names
		equilibrium_complex_ids = sim_data.process.equilibrium.ids_complexes
		self.monomer_ids = sim_data.process.translation.monomer_data["id"].tolist()

		# Get IDs of molecules involved in two component system complexes:
		two_component_system_molecule_ids = (
			list(sim_data.process.two_component_system.modified_molecules))
		two_component_system_complex_ids = (
			list(sim_data.process.two_component_system.complex_to_monomer.keys()))

		# Get IDs of monomers within ribosome subunits for unpacking of active
		# ribosomes (stored within unique molecules):
		ribosome_50s_subunits = sim_data.process.complexation.get_monomers(
			sim_data.molecule_ids.s50_full_complex)
		ribosome_30s_subunits = sim_data.process.complexation.get_monomers(
			sim_data.molecule_ids.s30_full_complex)
		ribosome_subunit_ids = (ribosome_50s_subunits["subunitIds"].tolist() +
			ribosome_30s_subunits["subunitIds"].tolist())

		# Get IDs of RNA polymerase subunits for unpacking of active RNA
		# polymerases (stored within unique molecules):
		rnap_subunits = sim_data.process.complexation.get_monomers(
			sim_data.molecule_ids.full_RNAP)
		rnap_subunit_ids = rnap_subunits["subunitIds"].tolist()

		# Get IDs of replisome subunits for unpacking of active replisomes
		# (stored within unique molecules):
		replisome_trimer_subunits = sim_data.molecule_groups.replisome_trimer_subunits
		replisome_monomer_subunits = sim_data.molecule_groups.replisome_monomer_subunits
		replisome_subunit_ids = replisome_trimer_subunits + replisome_monomer_subunits

		# Get the IDs of transcription factor subunits for unpacking active TFs
		# (stored within unique molecules):
		tfs = sim_data.process.transcription_regulation.tf_ids
		tf_subunit_ids = [tf_id + f'[{sim_data.getter.get_compartment(tf_id)[0]}]'
						  for tf_id in tfs]

		# Get stoichiometric matrices for complexation, equilibrium, two
		# component system complexes and the assembly of active unique molecules
		self.complexation_stoich = sim_data.process.complexation.stoich_matrix_monomers()
		self.equilibrium_stoich = sim_data.process.equilibrium.stoich_matrix_monomers()
		self.two_component_system_stoich = (
			sim_data.process.two_component_system.stoich_matrix_monomers())
		self.ribosome_stoich = np.hstack(
			(ribosome_50s_subunits["subunitStoich"],
			ribosome_30s_subunits["subunitStoich"]))
		self.rnap_stoich = rnap_subunits["subunitStoich"]
		self.replisome_stoich = np.hstack(
			(3*np.ones(len(replisome_trimer_subunits)),
			np.ones(len(replisome_monomer_subunits))))

		# Construct dictionary to quickly find bulk molecule indexes from IDs
		molecule_dict = {mol: i for i, mol in enumerate(bulk_molecule_ids)}

		def get_molecule_indexes(keys):
			return np.array([molecule_dict[x] for x in keys])

		# Get indexes of all relevant bulk molecules
		self.monomer_idx = get_molecule_indexes(self.monomer_ids)
		self.complexation_molecule_idx = get_molecule_indexes(complexation_molecule_ids)
		self.complexation_complex_idx = get_molecule_indexes(complexation_complex_ids)
		self.equilibrium_molecule_idx = get_molecule_indexes(equilibrium_molecule_ids)
		self.equilibrium_complex_idx = get_molecule_indexes(equilibrium_complex_ids)
		self.two_component_system_molecule_idx = (
			get_molecule_indexes(two_component_system_molecule_ids))
		self.two_component_system_complex_idx = (
			get_molecule_indexes(two_component_system_complex_ids))
		self.ribosome_subunit_idx = get_molecule_indexes(ribosome_subunit_ids)
		self.rnap_subunit_idx = get_molecule_indexes(rnap_subunit_ids)
		self.replisome_subunit_idx = get_molecule_indexes(replisome_subunit_ids)
		self.bound_TF_subunit_idx = get_molecule_indexes(list(tf_subunit_ids))

		# Get indexes of all unique molecules that need to be accounted for
		self.uniqueMolecules = sim.internal_states["UniqueMolecules"]
		unique_molecule_ids = self.uniqueMolecules.container.objectNames()
		self.ribosome_idx = unique_molecule_ids.index('active_ribosome')
		self.rnap_idx = unique_molecule_ids.index('active_RNAP')
		self.replisome_idx = unique_molecule_ids.index("active_replisome")


	def allocate(self):
		super(MonomerCounts, self).allocate()

		self.monomerCounts = np.zeros(
			len(self.monomer_ids),
			np.int64
			)

		self.freeMonomerCounts = np.zeros(
			len(self.monomer_ids),
			np.int64
		)

		self.monomersDegraded = np.zeros(
			len(self.monomer_ids),
			np.int64
		)

		self.monomersElongated = np.zeros(
			len(self.monomer_ids),
			np.int64
		)


	def update(self):
		# Get current counts of bulk and unique molecules:
		bulkMoleculeCounts = self.bulkMolecules.container.counts()
		uniqueMoleculeCounts = self.uniqueMolecules.container.counts()

		# Get promoter objects to determine the number of bound TFs:
		promoters = self.uniqueMolecules.container.objectsInCollection('promoter')

		# Get current counts of active unique molecules:
		n_active_ribosome = uniqueMoleculeCounts[self.ribosome_idx]
		n_active_rnap = uniqueMoleculeCounts[self.rnap_idx]
		n_active_replisome = uniqueMoleculeCounts[self.replisome_idx]
		n_bound_TFs = promoters.attr('bound_TF')

		# Calculate the subunit counts using stoich matrices:
		n_ribosome_subunit = n_active_ribosome * self.ribosome_stoich
		n_rnap_subunit = n_active_rnap * self.rnap_stoich
		n_replisome_subunit = n_active_replisome * self.replisome_stoich
		n_bound_TF_subunit = n_bound_TFs.sum(axis=0)

		# Add the counts of all active unique molecule complex subunits to the
		# free counts of each "inactive" subunit in the bulk (note this must
		# happen before the bulk molecule complexes are unpacked, as some of
		# these subunits are bulk molecule complexes):
		bulkMoleculeCounts[self.ribosome_subunit_idx] += n_ribosome_subunit.astype(int)
		bulkMoleculeCounts[self.rnap_subunit_idx] += n_rnap_subunit.astype(int)
		bulkMoleculeCounts[self.replisome_subunit_idx] += n_replisome_subunit.astype(int)
		bulkMoleculeCounts[self.bound_TF_subunit_idx] += n_bound_TF_subunit.astype(int)

		# Account for monomer subunits that make up TCS complexes
		# (note: TCS complex subunits can be equilibrium or complexation
		# complexes, so these must be unpacked before those):
		two_component_monomer_counts = np.dot(self.two_component_system_stoich,
			np.negative(bulkMoleculeCounts[self.two_component_system_complex_idx]))
		bulkMoleculeCounts[self.two_component_system_molecule_idx] += (
			two_component_monomer_counts.astype(int))

		# Account for monomer subunits that make up equilibrium complexes
		# (note: eq complex subunits can be complexation complexes, so these
		# must be unpacked before those):
		equilibrium_monomer_counts = np.dot(self.equilibrium_stoich,
			np.negative(bulkMoleculeCounts[self.equilibrium_complex_idx]))
		bulkMoleculeCounts[self.equilibrium_molecule_idx] += (
			equilibrium_monomer_counts.astype(int))

		# Account for monomer subunits that make up complexation complexes:
		complex_monomer_counts = np.dot(self.complexation_stoich,
			np.negative(bulkMoleculeCounts[self.complexation_complex_idx]))
		bulkMoleculeCounts[self.complexation_molecule_idx] += (
			complex_monomer_counts.astype(int))

		# Update the total and free monomer counts at the start of each time step:
		self.monomerCounts = bulkMoleculeCounts[self.monomer_idx]
		self.freeMonomerCounts = (
			self.bulkMolecules.container.counts())[self.monomer_idx]

	def tableCreate(self, tableWriter):
		subcolumns = {
			'monomerCounts': 'monomerIds',
			'freeMonomerCounts': 'monomerIds',
			'monomersElongated': 'monomerIds',
			'monomersDegraded': 'monomerIds',
			}

		tableWriter.writeAttributes(
			monomerIds = self.monomer_ids,
			subcolumns = subcolumns)

	def tableAppend(self, tableWriter):
		tableWriter.append(
			time = self.time(),
			simulationStep = self.simulationStep(),
			monomerCounts = self.monomerCounts,
			freeMonomerCounts=self.freeMonomerCounts,
			monomersElongated = self.monomersElongated,
			monomersDegraded = self.monomersDegraded,
			)
