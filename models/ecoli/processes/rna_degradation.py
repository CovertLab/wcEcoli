"""
Submodel for RNA degradation.

Mathematical formulation:

dr/dt = Kb - kcatEndoRNase * EndoRNase * r/Km / (1 + Sum(r/Km))

	where	r = RNA counts
			Kb = RNA production given a RNAP synthesis rate
			kcatEndoRNase = enzymatic activity for EndoRNases
			Km = Michaelis-Menten constants fitted to recapitulate first-order
			RNA decay:
				kd * r = kcatEndoRNase * EndoRNase * r/Km / (1 + Sum(r/Km))

This sub-model encodes molecular simulation of RNA degradation as two main
steps guided by RNases, "endonucleolytic cleavage" and "exonucleolytic
digestion":

1. Compute total counts of RNA to be degraded (D) and total capacity for
endo-cleavage (C) at each time point
2. Evaluate C and D. If C > D, then define a fraction of active endoRNases 
3. Dissect RNA degraded into different species (mRNA, tRNA, and rRNA) by
accounting endoRNases specificity
4. Update RNA fragments (assumption: fragments are represented as a pool of
nucleotides) created because of endonucleolytic cleavage
5. Compute total capacity of exoRNases and determine fraction of nucleotides
that can be digested
6. Update pool of metabolites (H and H2O) created because of exonucleolytic
digestion
"""

import numpy as np

import wholecell.processes.process
from wholecell.utils import units
from wholecell.utils.constants import REQUEST_PRIORITY_DEGRADATION

class RnaDegradation(wholecell.processes.process.Process):
	""" RnaDegradation """

	_name = "RnaDegradation"

	# Constructor
	def __init__(self):
		super(RnaDegradation, self).__init__()

	# Construct object graph
	def initialize(self, sim, sim_data):
		super(RnaDegradation, self).initialize(sim, sim_data)
		transcription = sim_data.process.transcription

		all_rna_ids = (
			list(transcription.rna_data['id'])
			+ list(transcription.mature_rna_data['id']))
		rna_id_to_index = {rna_id: i for (i, rna_id) in enumerate(all_rna_ids)}
		self.n_total_RNAs = len(all_rna_ids)

		# Load constants
		self.n_avogadro = sim_data.constants.n_avogadro
		self.cell_density = sim_data.constants.cell_density

		# Load RNase kinetic parameters
		endoRNase_ids = sim_data.process.rna_decay.endoRNase_ids
		exoRNase_ids = sim_data.molecule_groups.exoRNases
		self.Kcat_endoRNases = sim_data.process.rna_decay.kcats
		self.kcat_exoRNase = sim_data.constants.kcat_exoRNase

		# Load information about charged tRNAs
		uncharged_trna_names = transcription.uncharged_trna_names
		charged_trna_names = transcription.charged_trna_names
		self.uncharged_trna_indexes = np.array([
			rna_id_to_index[trna_id] for trna_id in uncharged_trna_names
			])

		# Load first-order RNA degradation rate constants (estimated by mRNA
		# half-life data)
		self.rna_deg_rates = (1 / units.s) * np.concatenate((
			transcription.rna_data['deg_rate'].asNumber(1/units.s),
			transcription.mature_rna_data['deg_rate'].asNumber(1/units.s)
			))

		self.is_mRNA = np.concatenate((
			transcription.rna_data['is_mRNA'].astype(np.int64),
			np.zeros(len(transcription.mature_rna_data), np.int64)
			))  # All mature RNAs are not mRNAs
		self.is_rRNA = np.concatenate((
			transcription.rna_data['is_rRNA'].astype(np.int64),
			transcription.mature_rna_data['is_rRNA'].astype(np.int64)
			))
		self.is_tRNA = np.concatenate((
			transcription.rna_data['is_tRNA'].astype(np.int64),
			transcription.mature_rna_data['is_tRNA'].astype(np.int64)
			))

		# Load lengths and nucleotide counts for all degradable RNAs
		self.rna_lengths = np.concatenate((
			transcription.rna_data['length'].asNumber(),
			transcription.mature_rna_data['length'].asNumber()
			))
		nt_counts = np.concatenate((
			transcription.rna_data['counts_ACGU'].asNumber(units.nt),
			transcription.mature_rna_data['counts_ACGU'].asNumber(units.nt)
			))

		# Build stoichiometric matrix
		polymerized_ntp_ids = sim_data.molecule_groups.polymerized_ntps
		end_cleavage_metabolite_ids = polymerized_ntp_ids + [
			sim_data.molecule_ids.water,
			sim_data.molecule_ids.ppi,
			sim_data.molecule_ids.proton
			]
		nmp_idx = np.arange(4)
		water_idx = end_cleavage_metabolite_ids.index(sim_data.molecule_ids.water)
		ppi_idx = end_cleavage_metabolite_ids.index(sim_data.molecule_ids.ppi)
		proton_idx = end_cleavage_metabolite_ids.index(sim_data.molecule_ids.proton)
		self.endo_degradation_stoich_matrix = np.zeros(
			(len(end_cleavage_metabolite_ids), self.n_total_RNAs), np.int64)
		self.endo_degradation_stoich_matrix[nmp_idx, :] = nt_counts.T
		self.endo_degradation_stoich_matrix[water_idx, :] = 0
		self.endo_degradation_stoich_matrix[ppi_idx, :] = 1
		self.endo_degradation_stoich_matrix[proton_idx, :] = 0

		# Build Views
		self.bulk_RNAs = self.bulkMoleculesView(all_rna_ids)
		self.unique_RNAs = self.uniqueMoleculesView('RNA')
		self.water = self.bulkMoleculeView(sim_data.molecule_ids.water)
		self.nmps = self.bulkMoleculesView(sim_data.molecule_groups.nmps)
		self.proton = self.bulkMoleculeView(sim_data.molecule_ids.proton)

		self.fragment_metabolites = self.bulkMoleculesView(end_cleavage_metabolite_ids)
		self.fragment_bases = self.bulkMoleculesView(polymerized_ntp_ids)

		self.endoRNases = self.bulkMoleculesView(endoRNase_ids)
		self.exoRNases = self.bulkMoleculesView(exoRNase_ids)
		self.charged_trnas = self.bulkMoleculesView(charged_trna_names)

		# Set priority
		self.bulkMoleculesRequestPriorityIs(REQUEST_PRIORITY_DEGRADATION)

		# Load Michaelis-Menten constants fitted to recapitulate first-order RNA decay model
		self.Kms = (units.mol / units.L) * np.concatenate((
			transcription.rna_data['Km_endoRNase'].asNumber(units.mol/units.L),
			transcription.mature_rna_data['Km_endoRNase'].asNumber(units.mol/units.L)
			))


	def calculateRequest(self):
		# Compute factor that convert counts into concentration, and vice versa
		cell_mass = self.readFromListener("Mass", "cellMass") * units.fg
		cell_volume = cell_mass / self.cell_density
		counts_to_molar = 1 / (self.n_avogadro * cell_volume)

		# Get total counts of free rRNAs, uncharged and charged tRNAs, and
		# active (translatable) unique mRNAs
		bulk_RNA_counts = self.bulk_RNAs.total_counts().copy()
		bulk_RNA_counts[self.uncharged_trna_indexes] += self.charged_trnas.total_counts()

		TU_index, can_translate, is_full_transcript = self.unique_RNAs.attrs(
			'TU_index', 'can_translate', 'is_full_transcript')
		TU_index_translatable_mRNAs = TU_index[can_translate]
		unique_RNA_counts = np.bincount(
			TU_index_translatable_mRNAs, minlength=self.n_total_RNAs)
		total_RNA_counts = bulk_RNA_counts + unique_RNA_counts

		# Compute RNA concentrations
		rna_conc_molar = counts_to_molar * total_RNA_counts

		# Get counts of endoRNases
		endoRNase_counts = self.endoRNases.total_counts().copy()
		total_kcat_endornase = units.dot(self.Kcat_endoRNases, endoRNase_counts)

		# Calculate the fraction of active endoRNases for each RNA based on
		# Michaelis-Menten kinetics
		frac_endoRNase_saturated = (
			rna_conc_molar / self.Kms / (1 + units.sum(rna_conc_molar / self.Kms))
		).asNumber()

		# Calculate difference in degradation rates from first-order decay
		# and the number of EndoRNases per one molecule of RNA
		total_endornase_counts = np.sum(endoRNase_counts)
		diff_relative_first_order_decay = units.sum(
			units.abs(self.rna_deg_rates * total_RNA_counts -
					  total_kcat_endornase * frac_endoRNase_saturated)
			)
		endornase_per_rna = total_endornase_counts / np.sum(total_RNA_counts)

		self.writeToListener(
			"RnaDegradationListener", "FractionActiveEndoRNases",
			np.sum(frac_endoRNase_saturated)
			)
		self.writeToListener(
			"RnaDegradationListener", "DiffRelativeFirstOrderDecay",
			diff_relative_first_order_decay.asNumber()
			)
		self.writeToListener(
			"RnaDegradationListener", "FractEndoRRnaCounts", endornase_per_rna)

		# Dissect RNAse specificity into mRNA, tRNA, and rRNA
		mrna_specificity = np.dot(frac_endoRNase_saturated, self.is_mRNA)
		trna_specificity = np.dot(frac_endoRNase_saturated, self.is_tRNA)
		rrna_specificity = np.dot(frac_endoRNase_saturated, self.is_rRNA)

		n_total_mrnas_to_degrade = self._calculate_total_n_to_degrade(
			mrna_specificity, total_kcat_endornase)
		n_total_trnas_to_degrade = self._calculate_total_n_to_degrade(
			trna_specificity, total_kcat_endornase)
		n_total_rrnas_to_degrade = self._calculate_total_n_to_degrade(
			rrna_specificity, total_kcat_endornase)

		# Compute RNAse specificity
		rna_specificity = frac_endoRNase_saturated / np.sum(frac_endoRNase_saturated)

		# Boolean variable that tracks existence of each RNA
		rna_exists = (total_RNA_counts > 0).astype(np.int64)

		# Compute degradation probabilities of each RNA: for mRNAs and rRNAs,
		# this is based on the specificity of each RNA. For tRNAs, this is
		# distributed evenly.
		mrna_deg_probs = 1. / np.dot(rna_specificity, self.is_mRNA * rna_exists) * rna_specificity * self.is_mRNA * rna_exists
		rrna_deg_probs = 1. / np.dot(rna_specificity, self.is_rRNA * rna_exists) * rna_specificity * self.is_rRNA * rna_exists
		trna_deg_probs = 1. / np.dot(self.is_tRNA, rna_exists) * self.is_tRNA * rna_exists

		# Mask RNA counts into each class of RNAs
		mrna_counts = total_RNA_counts * self.is_mRNA
		trna_counts = total_RNA_counts * self.is_tRNA
		rrna_counts = total_RNA_counts * self.is_rRNA

		# Determine number of individual RNAs to be degraded for each class
		# of RNA.
		n_mrnas_to_degrade = self._get_rnas_to_degrade(
			n_total_mrnas_to_degrade, mrna_deg_probs, mrna_counts)
		n_trnas_to_degrade = self._get_rnas_to_degrade(
			n_total_trnas_to_degrade, trna_deg_probs, trna_counts)
		n_rrnas_to_degrade = self._get_rnas_to_degrade(
			n_total_rrnas_to_degrade, rrna_deg_probs, rrna_counts)
		n_RNAs_to_degrade = n_mrnas_to_degrade + n_trnas_to_degrade + n_rrnas_to_degrade

		# Bulk RNAs (tRNAs and rRNAs) are degraded immediately. Unique RNAs
		# (mRNAs) are immediately deactivated (becomes unable to bind
		# ribosomes), but not degraded until transcription is finished and the
		# mRNA becomes a full transcript to simplify the transcript elongation
		# process.
		n_bulk_RNAs_to_degrade = n_RNAs_to_degrade.copy()
		n_bulk_RNAs_to_degrade[self.is_mRNA.astype(bool)] = 0
		self.n_unique_RNAs_to_deactivate = n_RNAs_to_degrade.copy()
		self.n_unique_RNAs_to_deactivate[np.logical_not(self.is_mRNA.astype(bool))] = 0

		self.bulk_RNAs.requestIs(n_bulk_RNAs_to_degrade)
		self.unique_RNAs.request_access(self.EDIT_DELETE_ACCESS)
		self.fragment_bases.requestAll()

		# Calculate the amount of water required for total RNA hydrolysis by
		# endo and exonucleases. We first calculate the number of unique RNAs
		# that should be degraded at this timestep.
		self.unique_mRNAs_to_degrade = np.logical_and(
			np.logical_not(can_translate), is_full_transcript)
		self.n_unique_RNAs_to_degrade = np.bincount(
			TU_index[self.unique_mRNAs_to_degrade],
			minlength=self.n_total_RNAs)

		# Assuming complete hydrolysis for now. Note that one additional water
		# molecule is needed for each RNA to hydrolyze the 5' diphosphate.
		water_for_degraded_rnas = np.dot(
			n_bulk_RNAs_to_degrade + self.n_unique_RNAs_to_degrade,
			self.rna_lengths)
		water_for_fragments = self.fragment_bases.total_counts().sum()
		self.water.requestIs(water_for_degraded_rnas + water_for_fragments)
		

	def evolveState(self):
		# Get vector of numbers of RNAs to degrade for each RNA species
		n_degraded_bulk_RNA = self.bulk_RNAs.counts()
		n_degraded_unique_RNA = self.n_unique_RNAs_to_degrade
		n_degraded_RNA = n_degraded_bulk_RNA + n_degraded_unique_RNA

		self.writeToListener(
			"RnaDegradationListener", "count_RNA_degraded", n_degraded_RNA
			)
		self.writeToListener(
			"RnaDegradationListener", "nucleotidesFromDegradation",
			np.dot(n_degraded_RNA, self.rna_lengths)
			)

		# Degrade bulk RNAs
		self.bulk_RNAs.countsIs(0)

		# Deactivate and degrade unique RNAs
		TU_index, can_translate = self.unique_RNAs.attrs(
			'TU_index', 'can_translate')
		n_deactivated_unique_RNA = self.n_unique_RNAs_to_deactivate

		# Deactive unique RNAs
		non_zero_deactivation = (n_deactivated_unique_RNA > 0)

		for index, n_degraded in zip(
				np.arange(n_deactivated_unique_RNA.size)[non_zero_deactivation],
				n_deactivated_unique_RNA[non_zero_deactivation]):
			# Get mask for translatable mRNAs belonging to the degraded species
			mask = np.logical_and(TU_index == index, can_translate)

			# Choose n_degraded indexes randomly to deactivate
			can_translate[self.randomState.choice(
				size=n_degraded, a=np.where(mask)[0], replace=False)] = False

		self.unique_RNAs.attrIs(can_translate=can_translate)

		# Degrade full mRNAs that are inactive
		self.unique_RNAs.delByIndexes(
			np.where(self.unique_mRNAs_to_degrade)[0])

		# Modeling assumption: Once a RNA is cleaved by an endonuclease its
		# resulting nucleotides are lumped together as "polymerized fragments".
		# These fragments can carry over from previous timesteps. We are also
		# assuming that during endonucleolytic cleavage the 5'terminal
		# phosphate is removed. This is modeled as all of the fragments being
		# one long linear chain of "fragment bases".

		# Example:
		# PPi-Base-PO4(-)-Base-PO4(-)-Base-OH
		# ==>
		# Pi-FragmentBase-PO4(-)-FragmentBase-PO4(-)-FragmentBase + PPi
		# Note: Lack of -OH on 3' end of chain
		metabolites_endo_cleavage = np.dot(
			self.endo_degradation_stoich_matrix, n_degraded_RNA)

		# Increase polymerized fragment counts
		self.fragment_metabolites.countsInc(metabolites_endo_cleavage)

		# Check if exonucleolytic digestion can happen 
		if self.fragment_bases.counts().sum() == 0:
			return

		# Calculate exolytic cleavage events

		# Modeling assumption: We model fragments as one long fragment chain of
		# polymerized nucleotides. We are also assuming that there is no
		# sequence specificity or bias towards which nucleotides are
		# hydrolyzed.

		# Example:
		# Pi-FragmentBase-PO4(-)-FragmentBase-PO4(-)-FragmentBase + 3 H2O
		# ==>
		# 3 NMP + 3 H(+)
		# Note: Lack of -OH on 3' end of chain

		n_exoRNases = self.exoRNases.counts()
		n_fragment_bases = self.fragment_bases.counts()
		n_fragment_bases_sum = n_fragment_bases.sum()

		exornase_capacity = n_exoRNases.sum() * self.kcat_exoRNase * (
				units.s * self.timeStepSec())

		if exornase_capacity >= n_fragment_bases_sum:
			self.nmps.countsInc(n_fragment_bases)
			self.water.countDec(n_fragment_bases_sum)
			self.proton.countInc(n_fragment_bases_sum)
			self.fragment_bases.countsIs(0)

			total_fragment_bases_digested = n_fragment_bases_sum

		else:
			fragment_specificity = n_fragment_bases / n_fragment_bases_sum
			possible_bases_to_digest = self.randomState.multinomial(
				exornase_capacity, fragment_specificity)
			n_fragment_bases_digested = n_fragment_bases - np.fmax(
				n_fragment_bases - possible_bases_to_digest, 0)

			total_fragment_bases_digested = n_fragment_bases_digested.sum()

			self.nmps.countsInc(n_fragment_bases_digested)
			self.water.countDec(total_fragment_bases_digested)
			self.proton.countInc(total_fragment_bases_digested)
			self.fragment_bases.countsDec(n_fragment_bases_digested)

		self.writeToListener("RnaDegradationListener",
			"fragmentBasesDigested", total_fragment_bases_digested)


	def _calculate_total_n_to_degrade(self, specificity, total_kcat_endornase):
		"""
		Calculate the total number of RNAs to degrade for a specific class of
		RNAs, based on the specificity of endoRNases on that specific class and
		the total kcat value of the endoRNases.

		Args:
			specificity: Sum of fraction of active endoRNases for all RNAs
			in a given class
			total_kcat_endornase: The summed kcat of all existing endoRNases
		Returns:
			Total number of RNAs to degrade for the given class of RNAs
		"""
		return np.round(
			(specificity * total_kcat_endornase
			 * (units.s * self.timeStepSec())).asNumber()
			)


	def _get_rnas_to_degrade(self, n_total_rnas_to_degrade, rna_deg_probs,
			rna_counts):
		"""
		Distributes the total count of RNAs to degrade for each class of RNAs
		into individual RNAs, based on the given degradation probabilities
		of individual RNAs. The upper bound is set by the current count of the
		specific RNA.

		Args:
			n_total_rnas_to_degrade: Total number of RNAs to degrade for the
			given class of RNAs (integer, scalar)
			rna_deg_probs: Degradation probabilities of each RNA (vector of
			equal length to the total number of different RNAs)
			rna_counts: Current counts of each RNA molecule (vector of equal
			length to the total number of different RNAs)
		Returns:
			Vector of equal length to rna_counts, specifying the number of
			molecules to degrade for each RNA
		"""
		n_rnas_to_degrade = np.zeros_like(rna_counts)
		remaining_rna_counts = rna_counts

		while n_rnas_to_degrade.sum() < n_total_rnas_to_degrade and remaining_rna_counts.sum() != 0:
			n_rnas_to_degrade += np.fmin(
				self.randomState.multinomial(
					n_total_rnas_to_degrade - n_rnas_to_degrade.sum(),
					rna_deg_probs
					),
				remaining_rna_counts
				)
			remaining_rna_counts = rna_counts - n_rnas_to_degrade

		return n_rnas_to_degrade
