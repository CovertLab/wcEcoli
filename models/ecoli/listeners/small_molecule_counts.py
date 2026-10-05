"""
# TODO: find where the missing PPI counter is?
Small Molecule Counts Listener

Records the total counts of all trackable intracellular small molecules at each
timestep, plus extracellular small molecule concentrations in the media (NOTE:
the free counts for all small molecules is accessible via the bulk container).

The tracked set is built from the model's reaction network (see below), so it
spans the cell's small molecules broadly (i.e. metabolites, amino acids,
nucleotides (NTPs/dNTPs), cofactors, and inorganic ions).

reconstruction/ecoli/flat/metabolites.tsv has ~7,600 metabolite compounds,
but most never appear as a countable species in a running simulation. This
listener tracks every small molecule that BOTH (a) exists as a bulk molecule
(has a countable ID) AND (b) participates in an encoded reaction within the
model. The tracked set is built in initialize() as:

    ( homeostatic concentration targets (across all saved media)
      ∪ equilibrium-reaction ligands  (equilibrium.metabolite_set)
      ∪ every species in metabolism.reaction_stoich )
      ∩  bulk molecule ids   (filters by countable species)
      −  macromolecules      (proteins, RNAs, and equilibrium / TCS /
         complexation complexes that get looped in from equilibrium and TCS
         stoich maps and need to be removed as they are tracked elsewhere)

The exact tracked list for a given run is emitted as the
totalSmallMoleculeCounts field metadata (``smallMoleculeIds``).

In the default model media conditions, only a subset of tracked small
molecules are present in nonzero amounts over the course of the simulation;
the rest exist in the reaction network but sit at zero count in that
condition.

Arrays emitted: totalSmallMoleculeCounts and environmentSmallMoleculeConcentrations.
Per-process flux event arrays that are NOT reconstructable from other existing
listeners are also emitted (dntpUsedInReplication, nmpFromRnaDegradation,
ppGppReactionSmallMoleculeDelta, nmpFromRnaMaturation, ppiFromTranscription,
ppiFromChromosomeStructure, aaFromChromosomeStructure).

Total counts are computed by adding back small molecules currently sequestered
in non-bulk locations:

  1. Equilibrium complexes (e.g., TF-ligand (1CS), 2CS-ligand bound complexes):
     Small molecule ligands are bound inside these complexes and are released
     upon dissociation. Unpacked via equilibrium.stoich_matrix_monomers().

  2. Two component system (TCS) complexes (PHOSPHO-HK, PHOSPHO-RR, etc.):
     Each phosphorylated TCS molecule carries exactly 1 Pi[c] covalently.
     ATP is consumed and ADP is released to the free pool during
     phosphorylation, so only Pi needs to be recovered here.

  3. TCS complexes that contain equilibrium complexes as subunits
     (e.g. PHOSPHO-HK-LIGAND contains the HK-LIGAND eq complex, which
     contains a small molecule ligand). When HK-LIGAND is phosphorylated, the
     ligand is no longer counted by the equilibrium unpacking (HK-LIGAND
     count dropped) but is still sequestered in the TCS complex, so it is
     recovered here separately.

  4. TCS and equilibrium complexes can become bound to transcription units
     on DNA, so small molecules are technically in bound transcription factors
     (TFs) as well (which are tracked in the promoters table as the bound_TF
     field).

Small molecules consumed by ribosomes/RNAP/replisomes are spent (not
sequestered), so no unpacking is needed for those. Complexation complexes are
protein-only.

Extracellular small molecules are stored as concentrations (not counts) because
there is no single cell volume to convert with at the environment level.
These use a separate molecule ID list (environment_sm_ids).
"""

import numpy as np
import wholecell.listeners.listener


class SmallMoleculeCounts(wholecell.listeners.listener.Listener):
	"""
	Listener for total counts of intracellular small molecules,
	including small molecules within complexes and
	bound transcription factors (which can be both equilibrium complexes and
	TCS complexes). Also tracks extracellular small molecule concentrations.
	"""
	_name = 'SmallMoleculeCounts'

	def __init__(self, *args, **kwargs):
		super(SmallMoleculeCounts, self).__init__(*args, **kwargs)

	def initialize(self, sim, sim_data):
		super(SmallMoleculeCounts, self).initialize(sim, sim_data)

		# Obtain bulk molecule container and build index mapping:
		self.bulkMolecules = sim.internal_states['BulkMolecules']
		bulk_molecule_ids = self.bulkMolecules.container.objectNames()
		molecule_dict = {mol: i for i, mol in enumerate(bulk_molecule_ids)}
		bulk_id_set = set(bulk_molecule_ids)

		# Build the small molecule ID list: EVERY small molecule that
		# participates in a metabolic reaction AND exists as a bulk molecule,
		# across ALL compartments (includes homeostatic targets, equilibrium
		# ligands, and ppGpp):
		concentration_updates = sim_data.process.metabolism.concentration_updates
		sm_id_set = set()
		sm_id_set.add(sim_data.molecule_ids.ppGpp)
		for media_id in sim_data.external_state.saved_media.keys():
			exchanges = sim_data.external_state.exchange_data_from_media(media_id)
			conc_dict = concentration_updates.concentrations_based_on_nutrients(
				imports=exchanges['importExchangeMolecules'])
			sm_id_set.update(conc_dict.keys())
		sm_id_set.update(sim_data.process.equilibrium.metabolite_set)
		# Every small molecule appearing in any metabolic reaction, all compartments:
		for stoich in sim_data.process.metabolism.reaction_stoich.values():
			sm_id_set.update(stoich.keys())
		# Keep only actual bulk molecules (so the counts are storable):
		sm_id_set &= bulk_id_set

		# Exclude macromolecules: some proteins/complexes appear as participants
		# in metabolic / equilibrium / TCS reactions and need to be filtered out:
		macromolecule_ids = (
			set(sim_data.process.translation.monomer_data['id'])
			| set(sim_data.process.transcription.rna_data['id'])
			| set(sim_data.process.equilibrium.ids_complexes)
			| set(sim_data.process.two_component_system.complex_to_monomer.keys())
			| set(sim_data.process.complexation.ids_complexes))
		sm_id_set -= macromolecule_ids

		self.sm_ids = sorted(sm_id_set)
		sm_to_idx = {m: i for i, m in enumerate(self.sm_ids)}

		def get_indexes(keys):
			return np.array([molecule_dict[x] for x in keys])

		self.sm_idx = get_indexes(self.sm_ids)

		# Obtain equilibrium complex information for unpacking bound small molecules:
		equilibrium = sim_data.process.equilibrium
		eq_molecule_names = equilibrium.molecule_names
		eq_complex_ids = equilibrium.ids_complexes

		self.equilibrium_stoich = equilibrium.stoich_matrix_monomers()
		self.eq_complex_idx = get_indexes(eq_complex_ids)

		# Identify which rows of the equilibrium stoich matrix correspond to
		# tracked small molecules:
		eq_sm_rows = []
		eq_sm_to_tracked_idx = []
		for row_idx, mol_name in enumerate(eq_molecule_names):
			if mol_name in sm_to_idx:
				eq_sm_rows.append(row_idx)
				eq_sm_to_tracked_idx.append(sm_to_idx[mol_name])

		self.eq_sm_rows = np.array(eq_sm_rows)
		self.eq_sm_to_tracked_idx = np.array(eq_sm_to_tracked_idx)

		if len(self.eq_sm_rows) > 0:
			self.eq_sm_stoich = self.equilibrium_stoich[
				self.eq_sm_rows, :]
		else:
			self.eq_sm_stoich = np.zeros(
				(0, len(eq_complex_ids)), dtype=np.float64)

		# Obtain two component system complex information (each phosphorylated
		# TCS molecule carries exactly 1 Pi[c] covalently. ATP is consumed and
		# ADP is released to the free pool during phosphorylation, so only Pi
		# is "hidden" in these molecules):
		tcs = sim_data.process.two_component_system
		phosphorylated_mol_ids = list(
			tcs.independent_to_dependent_molecules.values())

		# Bulk indexes for all phosphorylated TCS molecules
		self.tcs_phosphorylated_idx = get_indexes(phosphorylated_mol_ids)

		# Index of Pi[c] in our small molecule list (may be absent in edge cases)
		pi_id = 'Pi[c]'
		if pi_id in sm_to_idx:
			self.pi_sm_idx = sm_to_idx[pi_id]
			self.track_tcs_pi = True
		else:
			self.pi_sm_idx = -1
			self.track_tcs_pi = False

		# TCS complexes that contain equilibrium-complex subunits with
		# small molecule ligands (e.g. PHOSPHO-HK-LIGAND contains HK-LIGAND
		# which contains a small molecule LIGAND). When HK-LIGAND is
		# phosphorylated to PHOSPHO-HK-LIGAND, the ligand is no longer tracked
		# by the equilibrium complex unpacking (HK-LIGAND count drops) but IS
		# still sequestered in PHOSPHO-HK-LIGAND (so it must be tracked
		# separately):
		equilibrium = sim_data.process.equilibrium
		eq_complex_id_set = set(equilibrium.ids_complexes)
		eq_mol_names = list(equilibrium.molecule_names)
		eq_stoich = equilibrium.stoich_matrix_monomers()

		# Build eq_complex → {sm_id: stoich_per_complex} mapping
		eq_complex_to_sms = {}
		for col_idx, cplx_id in enumerate(equilibrium.ids_complexes):
			for row_idx, mol_name in enumerate(eq_mol_names):
				if (mol_name in equilibrium.metabolite_set
						and eq_stoich[row_idx, col_idx] < 0):
					if cplx_id not in eq_complex_to_sms:
						eq_complex_to_sms[cplx_id] = {}
					eq_complex_to_sms[cplx_id][mol_name] = (
						-eq_stoich[row_idx, col_idx])

		# For each TCS complex, check if any of its reaction reactants are
		# equilibrium complexes containing small molecule ligands:
		tcs_mol_names = list(tcs.molecule_names)
		tcs_stoich = tcs.stoich_matrix()

		# tcs_complex_to_ligands: {tcs_complex_id: {sm_id: stoich}}
		tcs_complex_to_ligands = {}
		for tcs_cplx in tcs.complex_to_monomer.keys():
			if tcs_cplx not in tcs_mol_names:
				continue
			cplx_row = tcs_mol_names.index(tcs_cplx)
			# Find reactions where this TCS complex is produced (+stoich)
			prod_rxn_cols = np.where(tcs_stoich[cplx_row, :] > 0)[0]
			for rxn_col in prod_rxn_cols:
				# Find reactants consumed in that reaction
				reactant_rows = np.where(tcs_stoich[:, rxn_col] < 0)[0]
				for r_row in reactant_rows:
					r_mol = tcs_mol_names[r_row]
					if (r_mol in eq_complex_id_set
							and r_mol in eq_complex_to_sms):
						stoich_eq_per_tcs = -tcs_stoich[r_row, rxn_col]
						for sm_id, stoich_sm in (
								eq_complex_to_sms[r_mol].items()):
							if sm_id not in sm_to_idx:
								continue
							if tcs_cplx not in tcs_complex_to_ligands:
								tcs_complex_to_ligands[tcs_cplx] = {}
							tcs_complex_to_ligands[tcs_cplx][sm_id] = (
								stoich_eq_per_tcs * stoich_sm)

		# Store as parallel arrays for fast update() computation
		self.tcs_ligand_complex_ids = []
		self.tcs_ligand_sm_idxs = []
		self.tcs_ligand_stoichs = []
		for cplx_id, sm_stoichs in tcs_complex_to_ligands.items():
			for sm_id, stoich in sm_stoichs.items():
				if sm_id in sm_to_idx:
					self.tcs_ligand_complex_ids.append(cplx_id)
					self.tcs_ligand_sm_idxs.append(sm_to_idx[sm_id])
					self.tcs_ligand_stoichs.append(stoich)

		if self.tcs_ligand_complex_ids:
			self.tcs_ligand_complex_idx = get_indexes(
				self.tcs_ligand_complex_ids)
			self.tcs_ligand_sm_idxs = np.array(self.tcs_ligand_sm_idxs)
			self.tcs_ligand_stoichs = np.array(self.tcs_ligand_stoichs)
			self.track_tcs_ligands = True
		else:
			self.track_tcs_ligands = False

		# Small molecule ligands inside transcription factors bound to DNA.
		# When a TF that is an equilibrium complex (e.g. a ligand-bound
		# repressor like PdhR-pyruvate) binds a promoter, its complex leaves
		# the bulk pool, so its small molecule ligand drops out of
		# smallMoleculesInEquilibriumComplexes. We add it back here (analogous
		# to how monomer_counts.py adds back bound-TF protein subunits). Bound
		# TF counts come from the 'promoter' unique molecule's 'bound_TF'
		# attribute, whose columns are ordered by tf_ids.
		self.uniqueMolecules = sim.internal_states['UniqueMolecules']
		tf_ids = sim_data.process.transcription_regulation.tf_ids

		# Set of phosphorylated TCS molecules (each carries 1 Pi). Many active
		# TFs are phospho-response-regulators (PHOSPHO-PHOB, PHOSPHO-NARL, ...);
		# when bound to DNA they leave the bulk pool, so their Pi is no longer
		# counted by smallMoleculesInTCSPhosphorylation (which reads bulk). Add
		# it back here via the same bound-TF machinery (mapping the TF column
		# to +1 Pi).
		phospho_tcs_set = set(phosphorylated_mol_ids)

		self.tf_bound_col_idxs = []   # column in the bound_TF array
		self.tf_bound_sm_idxs = []    # our small molecule index
		self.tf_bound_stoichs = []    # small molecule stoich in the TF complex
		for col, tf_id in enumerate(tf_ids):
			tf_mol = tf_id + f'[{sim_data.getter.get_compartment(tf_id)[0]}]'
			# Case 1: TF is an equilibrium complex carrying a small molecule ligand
			if tf_mol in eq_complex_to_sms:
				for sm_id, stoich in eq_complex_to_sms[tf_mol].items():
					if sm_id in sm_to_idx:
						self.tf_bound_col_idxs.append(col)
						self.tf_bound_sm_idxs.append(sm_to_idx[sm_id])
						self.tf_bound_stoichs.append(stoich)
			# Case 2: TF is a phosphorylated TCS response regulator (carries
			# 1 Pi covalently)
			if tf_mol in phospho_tcs_set and self.track_tcs_pi:
				self.tf_bound_col_idxs.append(col)
				self.tf_bound_sm_idxs.append(self.pi_sm_idx)
				self.tf_bound_stoichs.append(1)

		if self.tf_bound_col_idxs:
			self.tf_bound_col_idxs = np.array(self.tf_bound_col_idxs)
			self.tf_bound_sm_idxs = np.array(self.tf_bound_sm_idxs)
			self.tf_bound_stoichs = np.array(self.tf_bound_stoichs)
			self.track_bound_tf_sm = True
		else:
			self.track_bound_tf_sm = False

		# The Environment container stores concentrations (not molecule counts)
		# because there is no single cell volume to convert with (and these use
		# a separate molecule ID list that needs to be accounted for).
		self.environment = sim.external_states['Environment']
		self.environment_sm_ids = list(self.environment._moleculeIDs)

		# Store the molecule ID lists so that changes in small molecules from
		# all processes that use them can be tracked:
		self.aa_ids = list(sim_data.molecule_groups.amino_acids)
		self.ntp_ids = ['ATP[c]', 'CTP[c]', 'GTP[c]', 'UTP[c]']
		self.dntp_ids = list(sim_data.molecule_groups.dntps)
		self.nmp_ids = list(sim_data.molecule_groups.nmps)

		# Protein degradation (releases AAs and consumes water):
		self.prot_deg_sm_ids = (
			list(sim_data.molecule_groups.amino_acids)
			+ [sim_data.molecule_ids.water])

		# tRNA charging in polypeptide elongation produces AMP + PPi and
		# consumes ATP + AAs. This is the main source of AMP in the sim and
		# is NOT captured by FBA deltaMetabolites.
		self.charging_molecule_ids = list(
			sim_data.process.transcription.charging_molecules)

		# ppGpp synthesis (GDPPYPHOSKIN-RXN: ATP + GDP → ppGpp + AMP) and
		# degradation (PPGPPSYN-RXN: ppGpp + H2O → GDP + PPi) run inside
		# polypeptide elongation outside FBA. ppGpp synthesis produces AMP,
		# which is the missing ~2000 AMP/timestep gap in the events balance.
		self.ppgpp_reaction_sm_ids = list(
			sim_data.process.metabolism.ppgpp_reaction_metabolites)

		# RNA maturation (rna_maturation.py) releases NMPs when pre-tRNAs and
		# pre-rRNAs are processed. Uses same NMP ordering as RNA degradation.
		self.rna_maturation_nmp_ids = list(sim_data.molecule_groups.nmps)

		# RNA degradation endo-nuclease cleavage fragment small molecules:
		# polymerized NTPs (fragment bases) + water, PPi, proton.
		# These are released at each endo-cleavage event before exo-nuclease
		# digestion converts fragment bases into free NMPs.
		self.endo_cleavage_sm_ids = (
			list(sim_data.molecule_groups.polymerized_ntps)
			+ [sim_data.molecule_ids.water,
			   sim_data.molecule_ids.ppi,
			   sim_data.molecule_ids.proton])

	def allocate(self):
		super(SmallMoleculeCounts, self).allocate()

		n = len(self.sm_ids)

		# Save free and total counts of all intracellular small molecules at
		# each timestep:
		self.freeSmallMoleculeCounts = np.zeros(n, np.int64)
		self.smallMoleculesInEquilibriumComplexes = np.zeros(n, np.int64)
		self.smallMoleculesInTCSPhosphorylation = np.zeros(n, np.int64)
		# Small molecule ligands sequestered in TCS complexes that contain
		# equilibrium complexes as subunits (e.g. LIGAND in PHOSPHO-HK-LIGAND)
		self.smallMoleculesInTCSComplexes = np.zeros(n, np.int64)
		# Small molecule ligands inside transcription factors bound to DNA
		self.smallMoleculesInBoundTFs = np.zeros(n, np.int64)
		self.totalSmallMoleculeCounts = np.zeros(n, np.int64)

		# Save extracellular small molecule concentrations at each timestep:
		self.environmentSmallMoleculeConcentrations = np.zeros(
			len(self.environment_sm_ids), np.float64)

		# Save changes in small molecule counts at each timestep.
		# NOTE: 7 event columns (aaUsedInTranslation, ntpUsedInTranscription,
		# ppiFromReplication, smallMoleculesFromProteinDegradation,
		# chargingMoleculeDeltaInTranslation, fragmentSmallMoleculesFromEndoCleavage,
		# ppiFromRnaMaturation) were removed because they are exactly
		# reconstructable analysis-side from OTHER existing listeners + sim_data
		# (see small_molecule_events.py). The columns kept below are NOT fully
		# reconstructable from existing listeners.
		self.dntpUsedInReplication = np.zeros(len(self.dntp_ids), np.int64)
		self.nmpFromRnaDegradation = np.zeros(len(self.nmp_ids), np.int64)
		self.ppGppReactionSmallMoleculeDelta = np.zeros(
			len(self.ppgpp_reaction_sm_ids), np.int64)
		self.nmpFromRnaMaturation = np.zeros(
			len(self.rna_maturation_nmp_ids), np.int64)
		# PPi[c] released/consumed coupled to polymerization, tracked per
		# process (these processes consume NTP/dNTP/NMP — already tracked — but
		# the coupled PPi was not). Scalars: net PPi count change for PPi[c].
		self.ppiFromTranscription = 0
		self.ppiFromChromosomeStructure = 0
		# Amino acids released when chromosome structure removes ribosomes
		# (replication-transcription conflicts) and degrades the incomplete
		# nascent polypeptides back to free amino acids.
		self.aaFromChromosomeStructure = np.zeros(len(self.aa_ids), np.int64)

	def update(self):
		bulk_counts = self.bulkMolecules.container.counts()

		# Free small molecules tracked in bulk molecule counts:
		self.freeSmallMoleculeCounts = bulk_counts[self.sm_idx].copy()

		# Extract counts of small molecules bound in equilibrium complexes:
		eq_bound = np.zeros(len(self.sm_ids), np.int64)
		if len(self.eq_sm_rows) > 0:
			eq_complex_counts = bulk_counts[self.eq_complex_idx]
			bound = np.dot(
				self.eq_sm_stoich,
				np.negative(eq_complex_counts))
			eq_bound[self.eq_sm_to_tracked_idx] += bound.astype(np.int64)
		self.smallMoleculesInEquilibriumComplexes = eq_bound

		# Obtain Pi counts in phosphorylated TCS complexes:
		tcs_bound = np.zeros(len(self.sm_ids), np.int64)
		if self.track_tcs_pi:
			n_phosphorylated = int(
				bulk_counts[self.tcs_phosphorylated_idx].sum())
			tcs_bound[self.pi_sm_idx] = n_phosphorylated
		self.smallMoleculesInTCSPhosphorylation = tcs_bound

		# Small molecule ligands in TCS complexes that contain eq-complex
		# subunits (e.g. the LIGAND small molecule in PHOSPHO-HK-LIGAND). These
		# are NOT counted by the equilibrium unpacking (because HK-LIGAND count
		# dropped when it was phosphorylated) so must be tracked separately.
		tcs_complex_bound = np.zeros(len(self.sm_ids), np.int64)
		if self.track_tcs_ligands:
			tcs_cplx_counts = bulk_counts[self.tcs_ligand_complex_idx]
			for i in range(len(self.tcs_ligand_complex_ids)):
				tcs_complex_bound[self.tcs_ligand_sm_idxs[i]] += int(
					tcs_cplx_counts[i] * self.tcs_ligand_stoichs[i])
		self.smallMoleculesInTCSComplexes = tcs_complex_bound

		# Small molecule ligands inside transcription factors bound to DNA.
		# When a TF that is an equilibrium complex binds a promoter, its
		# complex leaves the bulk pool, so its small molecule ligand is no
		# longer counted in smallMoleculesInEquilibriumComplexes. Add it back
		# here using the bound_TF counts from the promoter unique molecules.
		bound_tf_sm = np.zeros(len(self.sm_ids), np.int64)
		if self.track_bound_tf_sm:
			promoters = self.uniqueMolecules.container.objectsInCollection(
				'promoter')
			n_bound_TFs = promoters.attr('bound_TF')
			n_bound_per_tf = n_bound_TFs.sum(axis=0)
			for i in range(len(self.tf_bound_col_idxs)):
				col = self.tf_bound_col_idxs[i]
				bound_tf_sm[self.tf_bound_sm_idxs[i]] += int(
					n_bound_per_tf[col] * self.tf_bound_stoichs[i])
		self.smallMoleculesInBoundTFs = bound_tf_sm

		# Sum up total counts:
		self.totalSmallMoleculeCounts = (
			self.freeSmallMoleculeCounts
			+ self.smallMoleculesInEquilibriumComplexes
			+ self.smallMoleculesInTCSPhosphorylation
			+ self.smallMoleculesInTCSComplexes
			+ self.smallMoleculesInBoundTFs)

		# Obtain extracellular small molecule concentrations:
		self.environmentSmallMoleculeConcentrations = (
			self.environment.container.counts().copy())

	def tableCreate(self, tableWriter):
		subcolumns = {
			'totalSmallMoleculeCounts': 'smallMoleculeIds',
			'environmentSmallMoleculeConcentrations': 'environmentSmallMoleculeIds',
			'dntpUsedInReplication': 'dntpIds',
			'nmpFromRnaDegradation': 'nmpIds',
			'ppGppReactionSmallMoleculeDelta': 'ppGppReactionSmallMoleculeIds',
			'nmpFromRnaMaturation': 'rnaMaturationNmpIds',
			'aaFromChromosomeStructure': 'aaIds',
			}
		tableWriter.writeAttributes(
			smallMoleculeIds=self.sm_ids,
			environmentSmallMoleculeIds=self.environment_sm_ids,
			aaIds=self.aa_ids,
			dntpIds=self.dntp_ids,
			nmpIds=self.nmp_ids,
			ppGppReactionSmallMoleculeIds=self.ppgpp_reaction_sm_ids,
			rnaMaturationNmpIds=self.rna_maturation_nmp_ids,
			subcolumns=subcolumns)

	def tableAppend(self, tableWriter):
		tableWriter.append(
			time=self.time(),
			simulationStep=self.simulationStep(),
			totalSmallMoleculeCounts=self.totalSmallMoleculeCounts,
			environmentSmallMoleculeConcentrations=(
				self.environmentSmallMoleculeConcentrations),
			dntpUsedInReplication=self.dntpUsedInReplication,
			nmpFromRnaDegradation=self.nmpFromRnaDegradation,
			ppGppReactionSmallMoleculeDelta=(
				self.ppGppReactionSmallMoleculeDelta),
			nmpFromRnaMaturation=self.nmpFromRnaMaturation,
			ppiFromTranscription=self.ppiFromTranscription,
			ppiFromChromosomeStructure=self.ppiFromChromosomeStructure,
			aaFromChromosomeStructure=self.aaFromChromosomeStructure,
			)
