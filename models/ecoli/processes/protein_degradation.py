"""
ProteinDegradation

Protein degradation sub-model. Encodes molecular simulation of protein
degradation as a Poisson process.

Every copy of a protein monomer degrades with its own first-order rate
(deg_rate, from monomer_data), whether the copy is free or sits inside a
complexation or equilibrium complex. This matches the ParCa, which applies each
monomer's degradation rate to its total count (free + complexed) when it fits
expression.

- Free monomers: uses n ~ Poisson(k_p * dt * free_p) to decide how many degrade
- Complexes: every instance of complex c is "hit" with rate
  lambda_c = sum_p(s_pc * k_p), where s_pc is the number of copies of monomer p
  in c (nested complexes are flattened down to their base subunits). Each hit
  dissociates one instance and degrades ONE subunit copy, picked with
  probability s_pc * k_p / lambda_c; all other subunits (proteins, RNAs and
  small-molecule ligands) are returned to their free pools. In expectation each
  monomer copy in a complex degrades at k_p, same as a free copy. (Two subunits
  of one instance degrading in the same time step is second order in k_p * dt
  and is ignored.)

Complexes that are never degraded (DEGRADE_COMPLEXES_EXCLUDED):
- ribosomal subunits, RNA polymerase and replisome subunits, plus any complex
  nested inside them, since the active (unique-molecule) forms of these are not
  degraded either;
- every complex that is, or contains, a two-component-system molecule.

Bound transcription factors: this process runs right after TfUnbinding (see
models/ecoli/sim/simulation.py), when every TF has been returned to the bulk
pool and before TfBinding rebinds them, so promoter-bound TFs (e.g. ~90% of
LexA dimers) are degraded at the same rate as free ones.
"""

import numpy as np

import wholecell.processes.process
from wholecell.utils.constants import REQUEST_PRIORITY_DEGRADATION
from wholecell.utils import units

# Set to False to restore the original behavior (only free monomers degrade).
DEGRADE_COMPLEXES = True


class ProteinDegradation(wholecell.processes.process.Process):
	""" ProteinDegradation """

	_name = "ProteinDegradation"

	# Constructor
	def __init__(self):

		# Parameters
		self.proteinLengths = None		# Protein lengths
		self.proteinDegSMatrix = None	# Protein degradation stoichiometry matrix [metabolite x rna]

		# Views
		self.metabolites = None
		self.h2o = None
		self.proteins = None
		self.complexes = None
		self.eq_complexes = None
		self.released_molecules = None

		super(ProteinDegradation, self).__init__()

	# Construct object graph
	def initialize(self, sim, sim_data):
		super(ProteinDegradation, self).initialize(sim, sim_data)

		# Load protein degradation rates (based on N-end rule)
		self.rawDegRate = sim_data.process.translation.monomer_data['deg_rate'].asNumber(1 / units.s)

		# Build metabolite IDs for S matrix
		h2oId = [sim_data.molecule_ids.water]
		metaboliteIds = sim_data.molecule_groups.amino_acids + h2oId
		aaIdxs = np.arange(0, len(sim_data.molecule_groups.amino_acids))
		h2oIdx = metaboliteIds.index(sim_data.molecule_ids.water)

		# Build protein IDs for S matrix
		proteinIds = sim_data.process.translation.monomer_data['id']

		# Load protein length
		self.proteinLengths = sim_data.process.translation.monomer_data['length']

		# Build S matrix
		self.proteinDegSMatrix = np.zeros((len(metaboliteIds), len(proteinIds)), np.int64)
		self.proteinDegSMatrix[aaIdxs, :] = np.transpose(sim_data.process.translation.monomer_data['aa_counts'].asNumber())
		self.proteinDegSMatrix[h2oIdx, :]  = -(np.sum(self.proteinDegSMatrix[aaIdxs, :], axis = 0) - 1)

		# Build Views
		self.metabolites = self.bulkMoleculesView(metaboliteIds)
		self.h2o = self.bulkMoleculeView(sim_data.molecule_ids.water)
		self.proteins = self.bulkMoleculesView(proteinIds)

		complexIds = np.array(sim_data.process.complexation.ids_complexes)
		eqComplexIds = np.array(sim_data.process.equilibrium.ids_complexes)
		self.complexes = self.bulkMoleculesView(complexIds)
		self.eq_complexes = self.bulkMoleculesView(eqComplexIds)

		# Add the IDs for proteins and complexes to self:
		self.protein_IDs = proteinIds
		self.complex_IDs = complexIds
		self.eq_complex_IDs = eqComplexIds

		self._build_complex_degradation(sim_data)

		# Obtain Avogadro's number for conversion of counts to molar concentration:
		self.n_avogadro = sim_data.constants.n_avogadro.asNumber(units.mol ** -1)

		self.bulkMoleculesRequestPriorityIs(REQUEST_PRIORITY_DEGRADATION)

	def _build_complex_degradation(self, sim_data):
		"""
		Builds the flattened subunit matrix for complexation + equilibrium
		complexes, the per-complex hit rates, and the set of excluded complexes.
		"""
		complexation = sim_data.process.complexation
		equilibrium = sim_data.process.equilibrium
		cx_ids = self.complex_IDs.tolist()
		eq_ids = self.eq_complex_IDs.tolist()
		cx_id_set = set(cx_ids)
		eq_id_set = set(eq_ids)

		# Flattened subunits of each complex: {base subunit ID: stoich}.
		# Complexation's stoich_matrix_monomers() is already flattened to base
		# subunits (monomers, RNAs, ligands). Equilibrium's is not: its
		# complexation-complex subunits are expanded here with the complexation
		# matrix.
		cx_molecule_ids = complexation.molecule_names
		cx_smm = complexation.stoich_matrix_monomers()
		cx_subunits = []
		for j in range(len(cx_ids)):
			rows = np.where(cx_smm[:, j] < 0)[0]
			cx_subunits.append(
				{cx_molecule_ids[i]: int(-cx_smm[i, j]) for i in rows})

		eq_molecule_ids = equilibrium.molecule_names
		eq_smm = equilibrium.stoich_matrix_monomers()
		eq_subunits = []
		for j in range(len(eq_ids)):
			subunits = {}
			for i in np.where(eq_smm[:, j] < 0)[0]:
				subunit_id = eq_molecule_ids[i]
				stoich = int(-eq_smm[i, j])
				if subunit_id in cx_id_set:
					nested = cx_subunits[cx_ids.index(subunit_id)]
					for base_id, base_stoich in nested.items():
						subunits[base_id] = subunits.get(base_id, 0) + stoich*base_stoich
				else:
					subunits[subunit_id] = subunits.get(subunit_id, 0) + stoich
			eq_subunits.append(subunits)

		# Every molecule (complex or base subunit) that goes into each complex,
		# following nested complexes all the way down
		def direct_constituents(process):
			molecule_ids = process.molecule_names
			stoich = process.stoich_matrix()
			constituents = {}
			for r in range(stoich.shape[1]):
				for product_row in np.where(stoich[:, r] > 0)[0]:
					constituents.setdefault(molecule_ids[product_row], set()).update(
						molecule_ids[i] for i in np.where(stoich[:, r] < 0)[0])
			return constituents
		direct = direct_constituents(complexation)
		for complex_id, reactants in direct_constituents(equilibrium).items():
			direct.setdefault(complex_id, set()).update(reactants)

		def closure(complex_id, seen=None):
			if seen is None:
				seen = set()
			for subunit_id in direct.get(complex_id, ()):
				if subunit_id not in seen:
					seen.add(subunit_id)
					closure(subunit_id, seen)
			return seen

		# Structural complexes whose active forms are unique molecules
		structural_ids = {
			sim_data.molecule_ids.s30_full_complex,
			sim_data.molecule_ids.s50_full_complex,
			sim_data.molecule_ids.full_RNAP,
			} | set(sim_data.molecule_groups.replisome_trimer_subunits) | set(
			sim_data.molecule_groups.replisome_monomer_subunits)
		structural_ids &= (cx_id_set | eq_id_set)
		for complex_id in list(structural_ids):
			structural_ids |= (closure(complex_id) & (cx_id_set | eq_id_set))

		# Two-component-system proteins and complexes (the TCS molecule list
		# also holds ATP, ADP, Pi, etc., which are not used here)
		tcs = sim_data.process.two_component_system
		protein_id_set = set(self.protein_IDs.tolist())
		tcs_ids = {
			mol_id for mol_id in tcs.modified_molecules
			if mol_id in protein_id_set or mol_id in cx_id_set
			or mol_id in eq_id_set or mol_id.startswith('PHOSPHO-')
			} | set(tcs.complex_to_monomer.keys())

		def is_excluded(complex_id):
			constituents = closure(complex_id) | {complex_id}
			return bool(constituents & (structural_ids | tcs_ids))

		all_complex_ids = cx_ids + eq_ids
		all_subunits = cx_subunits + eq_subunits
		self.n_cx = len(cx_ids)
		self.complex_excluded = np.array(
			[is_excluded(c) for c in all_complex_ids], dtype=bool)
		self.excluded_complex_ids = [
			c for c, ex in zip(all_complex_ids, self.complex_excluded) if ex]

		# Base subunit matrix [base subunit x (complexation + equilibrium
		# complexes)], positive stoichiometry
		base_ids = sorted({s for subunits in all_subunits for s in subunits})
		base_index = {s: i for i, s in enumerate(base_ids)}
		self.subunit_matrix = np.zeros((len(base_ids), len(all_complex_ids)), np.int64)
		for j, subunits in enumerate(all_subunits):
			for subunit_id, stoich in subunits.items():
				self.subunit_matrix[base_index[subunit_id], j] = stoich
		self.released_molecules = self.bulkMoleculesView(base_ids)

		# Map base subunits that are protein monomers to monomer indexes
		protein_index = {p: i for i, p in enumerate(self.protein_IDs)}
		self.base_is_protein = np.array(
			[s in protein_index for s in base_ids], dtype=bool)
		self.base_protein_idx = np.array(
			[protein_index[s] for s in base_ids if s in protein_index], dtype=np.int64)

		# Per-subunit degradation rates within each complex (s_pc * k_p, 1/s),
		# zero for non-protein subunits and for excluded complexes
		base_deg_rate = np.zeros(len(base_ids))
		base_deg_rate[self.base_is_protein] = self.rawDegRate[self.base_protein_idx]
		self.subunit_deg_rates = self.subunit_matrix * base_deg_rate[:, None]
		self.subunit_deg_rates[:, self.complex_excluded] = 0
		self.complex_hit_rates = self.subunit_deg_rates.sum(axis=0)

	def calculateRequest(self):
		# Determine how many free proteins to degrade based on the degradation
		# rates and counts of each protein
		nProteinsToDegrade = np.fmin(
			self.randomState.poisson(self._proteinDegRates() * self.proteins.total_counts()),
			self.proteins.total_counts()
			)

		# Determine how many complexes to dissociate and which subunit of each
		# is degraded
		n_complexes = len(self.complex_hit_rates)
		self.complex_hits = np.zeros(n_complexes, np.int64)
		self.degraded_subunits = {}  # {complex index: degraded base subunits}
		if DEGRADE_COMPLEXES:
			complex_counts = np.concatenate((
				self.complexes.total_counts(), self.eq_complexes.total_counts()))
			self.complex_hits = np.fmin(
				self.randomState.poisson(
					self.complex_hit_rates * self.timeStepSec() * complex_counts),
				complex_counts).astype(np.int64)
			for j in np.where(self.complex_hits > 0)[0]:
				self.degraded_subunits[j] = self._pick_degraded_subunits(
					j, self.complex_hits[j])

		degraded_base = np.zeros(self.subunit_matrix.shape[0], np.int64)
		for picks in self.degraded_subunits.values():
			degraded_base += picks
		nProteinsFromComplexes = self._proteins_from_complexes(degraded_base)

		# Determine the number of hydrolysis reactions
		nReactions = np.dot(
			self.proteinLengths.asNumber(), nProteinsToDegrade + nProteinsFromComplexes)

		# Determine the amount of water required to degrade the selected proteins
		# Assuming one N-1 H2O is required per peptide chain length N
		self.h2o.requestIs(
			nReactions - np.sum(nProteinsToDegrade + nProteinsFromComplexes))
		self.proteins.requestIs(nProteinsToDegrade)
		self.complexes.requestIs(self.complex_hits[:self.n_cx])
		self.eq_complexes.requestIs(self.complex_hits[self.n_cx:])

	def evolveState(self):
		# Proteins and complexes actually allocated to this process
		free_degraded = self.proteins.counts()
		complex_hits = np.concatenate((
			self.complexes.counts(), self.eq_complexes.counts()))

		# If fewer complexes were allocated than requested, redraw which
		# subunits are degraded for those complexes
		degraded_subunits = self.degraded_subunits
		for j in np.where(complex_hits != self.complex_hits)[0]:
			degraded_subunits[j] = self._pick_degraded_subunits(j, complex_hits[j])
		degraded_base = np.zeros(self.subunit_matrix.shape[0], np.int64)
		for picks in degraded_subunits.values():
			degraded_base += picks

		# Subunits of dissociated complexes: one copy per hit is degraded, the
		# rest are returned to their free pools
		released_subunits = self.subunit_matrix.dot(complex_hits) - degraded_base
		complex_degraded = self._proteins_from_complexes(degraded_base)
		complex_released = np.zeros(len(self.protein_IDs), np.int64)
		complex_released[self.base_protein_idx] = released_subunits[self.base_is_protein]

		# Remove degraded free proteins and dissociated complexes (must happen
		# before released subunits are added back, which may share molecules)
		self.proteins.countsIs(0)
		self.complexes.countsIs(0)
		self.eq_complexes.countsIs(0)
		self.released_molecules.countsInc(released_subunits)

		# Degrade selected proteins, release amino acids from those proteins back into the cell,
		# and consume H_2O that is required for the degradation process
		self.metabolites.countsInc(np.dot(
			self.proteinDegSMatrix,
			free_degraded + complex_degraded
			))

		# Record how many monomers were degraded and where they came from:
		self.writeToListener(
			"MonomerCounts", "monomersDegraded", free_degraded + complex_degraded)
		self.writeToListener("MonomerCounts", "monomersDegradedFree", free_degraded)
		self.writeToListener(
			"MonomerCounts", "monomersDegradedFromComplexes", complex_degraded)
		self.writeToListener(
			"MonomerCounts", "monomersProducedViaDissociation", complex_released)

		# Record dissociated complexes:
		cx_dissociated = complex_hits[:self.n_cx]
		eq_dissociated = complex_hits[self.n_cx:]
		self.writeToListener(
			"ComplexationListener", "complexesDissociated", cx_dissociated)
		self.writeToListener(
			"ComplexationListener", "dissociationIterations", int(cx_dissociated.sum()))
		self.writeToListener(
			"EquilibriumListener", "complexesDissociated", eq_dissociated)
		self.writeToListener(
			"EquilibriumListener", "dissociationIterations", int(eq_dissociated.sum()))

	def _pick_degraded_subunits(self, complex_idx, n_hits):
		"""
		Picks which base subunit is degraded in each of n_hits dissociated
		instances of a complex, weighted by s_pc * k_p.
		"""
		rates = self.subunit_deg_rates[:, complex_idx]
		if n_hits == 0 or rates.sum() == 0:
			return np.zeros(len(rates), np.int64)
		return self.randomState.multinomial(n_hits, rates / rates.sum())

	def _proteins_from_complexes(self, degraded_base):
		"""
		Maps degraded base subunit counts onto monomer indexes.
		"""
		degraded = np.zeros(len(self.protein_IDs), np.int64)
		degraded[self.base_protein_idx] = degraded_base[self.base_is_protein]
		return degraded

	def _proteinDegRates(self):
		return self.rawDegRate * self.timeStepSec()
