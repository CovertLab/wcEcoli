"""
PolypeptideElongation

Translation elongation sub-model.

TODO:
- see the initiation process for more TODOs
"""

from __future__ import annotations

from typing import Any, Callable, Dict, Optional, Set, Tuple

from numba import njit
import numpy as np
from scipy.integrate import solve_ivp

import wholecell.processes.process
from wholecell.utils.polymerize import buildSequences, polymerize, computeMassIncrease
from wholecell.utils.random import stochasticRound
from wholecell.utils import units


CONC_UNITS = units.umol / units.L
REMOVED_FROM_CHARGING = {'L-SELENOCYSTEINE[c]'}


class PolypeptideElongation(wholecell.processes.process.Process):
	""" PolypeptideElongation """

	_name = "PolypeptideElongation"

	def __init__(self):
		super(PolypeptideElongation, self).__init__()

	def initialize(self, sim, sim_data):
		super(PolypeptideElongation, self).initialize(sim, sim_data)

		# Simulation options
		self.aa_supply_in_charging = sim._aa_supply_in_charging
		self.adjust_timestep_for_charging = sim._adjust_timestep_for_charging
		self.mechanistic_translation_supply = sim._mechanistic_translation_supply
		self.mechanistic_aa_transport = sim._mechanistic_aa_transport
		self.ppgpp_regulation = sim._ppgpp_regulation
		self.disable_ppgpp_elongation_inhibition = sim._disable_ppgpp_elongation_inhibition
		self.variable_elongation = sim._variable_elongation_translation
		self.variable_polymerize = self.ppgpp_regulation or self.variable_elongation
		translation_supply = sim._translationSupply
		trna_charging = sim._trna_charging

		constants = sim_data.constants
		translation = sim_data.process.translation
		transcription = sim_data.process.transcription

		self.max_time_step = translation.max_time_step

		# Load parameters
		self.n_avogadro = constants.n_avogadro
		proteinIds = translation.monomer_data['id']
		self.proteinLengths = translation.monomer_data["length"].asNumber()
		self.proteinSequences = translation.translation_sequences
		self.aaWeightsIncorporated = translation.translation_monomer_weights
		self.endWeight = translation.translation_end_weight
		self.make_elongation_rates = translation.make_elongation_rates
		self.next_aa_pad = translation.next_aa_pad

		self.ribosomeElongationRate = float(sim_data.growth_rate_parameters.ribosomeElongationRate.asNumber(units.aa / units.s))

		# Amino acid supply calculations
		self.translation_aa_supply = sim_data.translation_supply_rate
		self.import_threshold = sim_data.external_state.import_constraint_threshold

		# Used for figure in publication
		self.trpAIndex = np.where(proteinIds == "TRYPSYN-APROTEIN[c]")[0][0]

		# Create view onto actively elongating 70S ribosomes
		self.active_ribosomes = self.uniqueMoleculesView('active_ribosome')

		# Create views onto 30S and 50S ribosomal subunits for termination
		self.ribosome30S = self.bulkMoleculeView(sim_data.molecule_ids.s30_full_complex)
		self.ribosome50S = self.bulkMoleculeView(sim_data.molecule_ids.s50_full_complex)

		# Create view onto all proteins
		self.bulkMonomers = self.bulkMoleculesView(proteinIds)

		# Create views onto all polymerization reaction small molecules
		self.aaNames = sim_data.molecule_groups.amino_acids
		self.aas = self.bulkMoleculesView(self.aaNames)

		self.elngRateFactor = 1.

		# Data structures for charging
		self.aa_from_trna = transcription.aa_from_trna

		# Set modeling method
		if trna_charging:
			self.elongation_model = SteadyStateElongationModel(sim_data, self)
		elif translation_supply:
			self.elongation_model = TranslationSupplyElongationModel(sim_data, self)
		else:
			self.elongation_model = BaseElongationModel(sim_data, self)

		# Growth associated maintenance energy requirements for elongations
		self.gtpPerElongation = constants.gtp_per_translation
		## Need to account for ATP hydrolysis for charging that has been
		## removed from measured GAM (ATP -> AMP is 2 hydrolysis reactions)
		## if charging reactions are not explicitly modeled
		if not trna_charging:
			self.gtpPerElongation += 2
		## Variable for metabolism to read to consume required energy
		self.gtp_to_hydrolyze = 0

		self.aa_exchange_rates = CONC_UNITS / units.s * np.zeros(len(self.aaNames))

	def calculateRequest(self):
		# Set ribosome elongation rate based on simulation medium environment and elongation rate factor
		# which is used to create single-cell variability in growth rate
		# The maximum number of amino acids that can be elongated in a single timestep is set to 22 intentionally as the minimum number of padding values
		# on the protein sequence matrix is set to 22. If timesteps longer than 1.0s are used, this feature will lead to errors in the effective ribosome
		# elongation rate.

		current_media_id = self._external_states['Environment'].current_media_id

		# MODEL SPECIFIC: get ribosome elongation rate
		self.ribosomeElongationRate = self.elongation_model.elongation_rate()

		# If there are no active ribosomes, return immediately
		if self.active_ribosomes.total_count() == 0:
			return

		# Build sequences to request appropriate amount of amino acids to
		# polymerize for next timestep
		proteinIndexes, peptideLengths = self.active_ribosomes.attrs(
			'protein_index', 'peptide_length'
			)

		self.elongation_rates = self.make_elongation_rates(
			self.randomState,
			self.ribosomeElongationRate,
			self.timeStepSec(),
			self.variable_elongation)

		sequences = buildSequences(
			self.proteinSequences,
			proteinIndexes,
			peptideLengths,
			self.elongation_rates)

		sequenceHasAA = (sequences != polymerize.PAD_VALUE)
		aasInSequences = np.bincount(sequences[sequenceHasAA], minlength=21)

		# Calculate AA supply for expected doubling of protein
		dryMass = (self.readFromListener("Mass", "dryMass") * units.fg)
		translation_supply_rate = self.translation_aa_supply[current_media_id] * self.elngRateFactor
		mol_aas_supplied = translation_supply_rate * dryMass * self.timeStepSec() * units.s
		self.aa_supply = units.strip_empty_units(mol_aas_supplied * self.n_avogadro)
		self.writeToListener("RibosomeData", "translationSupply", translation_supply_rate.asNumber())

		# MODEL SPECIFIC: Calculate AA request
		fraction_charged, aa_counts_for_translation = self.elongation_model.request(aasInSequences)

		# Write to listeners
		self.writeToListener("GrowthLimits", "fraction_trna_charged", np.dot(fraction_charged, self.aa_from_trna))
		self.writeToListener("GrowthLimits", "aaPoolSize", self.aas.total_counts())
		self.writeToListener("GrowthLimits", "aaRequestSize", aa_counts_for_translation)

		# Request full access to active ribosome molecules
		self.active_ribosomes.request_access(self.EDIT_DELETE_ACCESS)

	def evolveState(self):
		# Set values for metabolism in case of early return
		self.gtp_to_hydrolyze = 0
		self.aa_count_diff = {}

		# Write allocation data to listener
		self.writeToListener("GrowthLimits", "aaAllocated", self.aas.counts())

		# Get number of active ribosomes
		n_active_ribosomes = self.active_ribosomes.total_count()
		self.writeToListener("GrowthLimits", "activeRibosomeAllocated", n_active_ribosomes)

		if n_active_ribosomes == 0:
			return

		# Build amino acids sequences for each ribosome to polymerize
		protein_indexes, peptide_lengths, positions_on_mRNA = self.active_ribosomes.attrs(
			'protein_index', 'peptide_length', 'pos_on_mRNA'
			)

		all_sequences = buildSequences(
			self.proteinSequences,
			protein_indexes,
			peptide_lengths,
			self.elongation_rates + self.next_aa_pad)
		sequences = all_sequences[:, :-self.next_aa_pad].copy()

		if sequences.size == 0:
			return

		# Calculate elongation resource capacity
		aaCountInSequence = np.bincount(sequences[(sequences != polymerize.PAD_VALUE)])
		total_aa_counts = self.aas.counts()

		# MODEL SPECIFIC: Get amino acid counts
		aa_counts_for_translation = self.elongation_model.final_amino_acids(total_aa_counts)

		# Using polymerization algorithm elongate each ribosome up to the limits
		# of amino acids, sequence, and GTP
		result = polymerize(
			sequences,
			aa_counts_for_translation,
			10000000, # Set to a large number, the limit is now taken care of in metabolism
			self.randomState,
			self.elongation_rates[protein_indexes],
			variable_elongation=self.variable_polymerize,
			)

		sequence_elongations = result.sequenceElongation
		aas_used = result.monomerUsages
		nElongations = result.nReactions

		next_amino_acid = all_sequences[np.arange(len(sequence_elongations)), sequence_elongations]
		next_amino_acid_count = np.bincount(next_amino_acid[next_amino_acid != polymerize.PAD_VALUE], minlength=21)

		# Update masses of ribosomes attached to polymerizing polypeptides
		added_protein_mass = computeMassIncrease(
			sequences,
			sequence_elongations,
			self.aaWeightsIncorporated
			)

		updated_lengths = peptide_lengths + sequence_elongations
		updated_positions_on_mRNA = positions_on_mRNA + 3*sequence_elongations

		didInitialize = (
			(sequence_elongations > 0) &
			(peptide_lengths == 0)
			)

		added_protein_mass[didInitialize] += self.endWeight

		# Write current average elongation to listener
		currElongRate = (sequence_elongations.sum() / n_active_ribosomes) / self.timeStepSec()
		self.writeToListener("RibosomeData", "effectiveElongationRate", currElongRate)

		# Update active ribosomes, terminating if necessary
		self.active_ribosomes.attrIs(
			peptide_length=updated_lengths,
			pos_on_mRNA=updated_positions_on_mRNA)
		self.active_ribosomes.add_submass_by_name("protein", added_protein_mass)

		# Ribosomes that reach the end of their sequences are terminated and
		# dissociated into 30S and 50S subunits. The polypeptide that they are polymerizing
		# is converted into a protein in BulkMolecules
		terminalLengths = self.proteinLengths[protein_indexes]

		didTerminate = (updated_lengths == terminalLengths)

		terminatedProteins = np.bincount(
			protein_indexes[didTerminate],
			minlength = self.proteinSequences.shape[0]
			)

		self.active_ribosomes.delByIndexes(np.where(didTerminate)[0])
		self.bulkMonomers.countsInc(terminatedProteins)

		nTerminated = didTerminate.sum()
		nInitialized = didInitialize.sum()

		self.ribosome30S.countInc(nTerminated)
		self.ribosome50S.countInc(nTerminated)

		# MODEL SPECIFIC: evolve
		# TODO: use something other than a class attribute to pass aa diff to metabolism
		net_charged, self.aa_count_diff = self.elongation_model.evolve(
			total_aa_counts, aas_used, next_amino_acid_count, nElongations, nInitialized)

		# GTP hydrolysis is carried out in Metabolism process for growth
		# associated maintenance. This is set here for metabolism to use.
		self.gtp_to_hydrolyze = self.gtpPerElongation * nElongations

		# Write data to listeners
		self.writeToListener("GrowthLimits", "net_charged", net_charged)
		self.writeToListener("GrowthLimits", "aasUsed", aas_used)
		self.writeToListener("GrowthLimits", "aaCountDiff", [self.aa_count_diff.get(id_, 0) for id_ in self.aaNames])

		self.writeToListener("RibosomeData", "aaCountInSequence", aaCountInSequence)
		self.writeToListener("RibosomeData", "aaCounts", aa_counts_for_translation)

		self.writeToListener("RibosomeData", "actualElongations", sequence_elongations.sum())
		self.writeToListener("RibosomeData", "actualElongationHist", np.histogram(sequence_elongations, bins = np.arange(0,23))[0])
		self.writeToListener("RibosomeData", "elongationsNonTerminatingHist", np.histogram(sequence_elongations[~didTerminate], bins=np.arange(0,23))[0])

		self.writeToListener("RibosomeData", "didTerminate", didTerminate.sum())
		self.writeToListener("RibosomeData", "terminationLoss", (terminalLengths - peptide_lengths)[didTerminate].sum())
		self.writeToListener("RibosomeData", "numTrpATerminated", terminatedProteins[self.trpAIndex])

		self.writeToListener("RibosomeData", "processElongationRate", self.ribosomeElongationRate / self.timeStepSec())

	def isTimeStepShortEnough(self, inputTimeStep, timeStepSafetyFraction):
		model_specific = self.elongation_model.isTimeStepShortEnough(inputTimeStep, timeStepSafetyFraction)
		max_time_step = inputTimeStep <= self.max_time_step
		return model_specific and max_time_step


class BaseElongationModel(object):
	"""
	Base Model: Request amino acids according to upcoming sequence, assuming
	max ribosome elongation.
	"""
	def __init__(self, sim_data, process):
		self.process = process
		self.basal_elongation_rate = sim_data.constants.ribosome_elongation_rate_basal.asNumber(units.aa / units.s)
		self.ribosomeElongationRateDict = sim_data.process.translation.ribosomeElongationRateDict
		self.uncharged_trna_names = sim_data.process.transcription.uncharged_trna_names
		self.aaNames = sim_data.molecule_groups.amino_acids
		self.proton = self.process.bulkMoleculeView(sim_data.molecule_ids.proton)
		self.water = self.process.bulkMoleculeView(sim_data.molecule_ids.water)

	def elongation_rate(self):
		current_media_id = self.process._external_states['Environment'].current_media_id
		rate = self.process.elngRateFactor * self.ribosomeElongationRateDict[
			current_media_id].asNumber(units.aa / units.s)
		return np.min([self.basal_elongation_rate, rate])

	def amino_acid_counts(self, aasInSequences):
		return aasInSequences

	def request(self, aasInSequences):
		aa_counts_for_translation = self.amino_acid_counts(aasInSequences)

		self.process.aas.requestIs(aa_counts_for_translation)

		# Not modeling charging so set fraction charged to 0 for all tRNA
		fraction_charged = np.zeros(len(self.aaNames))

		return fraction_charged, aa_counts_for_translation

	def final_amino_acids(self, total_aa_counts):
		return total_aa_counts

	def evolve(self, total_aa_counts, aas_used, next_amino_acid_count, nElongations, nInitialized):
		# Update counts of amino acids and water to reflect polymerization reactions
		self.process.aas.countsDec(aas_used)
		self.water.countInc(nElongations - nInitialized)
		net_charged = np.zeros(len(self.uncharged_trna_names))

		return net_charged, {}

	def isTimeStepShortEnough(self, inputTimeStep, timeStepSafetyFraction):
		return True

class TranslationSupplyElongationModel(BaseElongationModel):
	"""
	Translation Supply Model: Requests minimum of 1) upcoming amino acid
	sequence assuming max ribosome elongation (ie. Base Model) and 2) estimation
	based on doubling the proteome in one cell cycle (does not use ribosome
	elongation, computed in Parca).
	"""
	def __init__(self, sim_data, process):
		super(TranslationSupplyElongationModel, self).__init__(sim_data, process)

	def elongation_rate(self):
		return self.basal_elongation_rate

	def amino_acid_counts(self, aasInSequences):
		return np.fmin(self.process.aa_supply, aasInSequences)  # Check if this is required. It is a better request but there may be fewer elongations.

class SteadyStateElongationModel(TranslationSupplyElongationModel):
	"""
	Steady State Charging Model: Requests amino acids based on the
	Michaelis-Menten competitive inhibition model.
	"""
	def __init__(self, sim_data, process):
		super(SteadyStateElongationModel, self).__init__(sim_data, process)
		constants = sim_data.constants
		transcription = sim_data.process.transcription
		metabolism = sim_data.process.metabolism
		molecule_ids = sim_data.molecule_ids
		molecule_groups = sim_data.molecule_groups

		# Cell parameters
		self.cellDensity = constants.cell_density

		# Names of molecules associated with tRNA charging
		self.charged_trna_names = transcription.charged_trna_names
		self.charging_molecule_names = transcription.charging_molecules
		self.synthetase_names = transcription.synthetase_names

		# Data structures for charging
		self.aa_from_synthetase = transcription.aa_from_synthetase
		self.charging_stoich_matrix = transcription.charging_stoich_matrix()
		self.charging_molecules_not_aa = np.array([
			mol not in set(self.aaNames)
			for mol in self.charging_molecule_names
			])

		# Create views for tRNA charging molecules
		self.uncharged_trna = self.process.bulkMoleculesView(self.uncharged_trna_names)
		self.charged_trna = self.process.bulkMoleculesView(self.charged_trna_names)
		self.charging_molecules = self.process.bulkMoleculesView(self.charging_molecule_names)
		self.synthetases = self.process.bulkMoleculesView(self.synthetase_names)

		# ppGpp synthesis
		self.ppgpp_reaction_metabolites = self.process.bulkMoleculesView(metabolism.ppgpp_reaction_metabolites)
		self.rela = self.process.bulkMoleculeView(molecule_ids.RelA)
		self.spot = self.process.bulkMoleculeView(molecule_ids.SpoT)
		self.ppgpp = self.process.bulkMoleculeView(molecule_ids.ppGpp)
		self.elong_rate_by_ppgpp = sim_data.growth_rate_parameters.get_ribosome_elongation_rate_by_ppgpp

		# Parameters for tRNA charging, ribosome elongation and ppGpp reactions
		self.charging_params = get_charging_params(sim_data,
			variable_elongation=self.process.variable_elongation)
		self.ppgpp_params = get_ppgpp_params(sim_data)

		# Amino acid supply calculations
		self.aa_supply_scaling = metabolism.aa_supply_scaling
		self.aa_environment = self.process.environmentView([aa[:-3] for aa in self.aaNames])

		# Manage unstable charging with too long time step by setting
		# time_step_short_enough to False during updates. Other variables
		# manage when to trigger an adjustment and how quickly the time step
		# increases after being reduced
		self.time_step_short_enough = True
		self.max_time_step = self.process.max_time_step
		self.time_step_increase = 1.01
		self.max_amino_acid_adjustment = 0.05

		self.aa_enzymes = self.process.bulkMoleculesView(metabolism.aa_enzymes)
		self.aa_aas = self.process.bulkMoleculesView(molecule_groups.amino_acids)
		self.amino_acid_synthesis = metabolism.amino_acid_synthesis
		self.amino_acid_import = metabolism.amino_acid_import
		self.amino_acid_export = metabolism.amino_acid_export
		self.get_pathway_enzyme_counts_per_aa = metabolism.get_pathway_enzyme_counts_per_aa

		self.aa_importers = self.process.bulkMoleculesView(metabolism.aa_importer_names)
		self.aa_exporters = self.process.bulkMoleculesView(metabolism.aa_exporter_names)

	def elongation_rate(self):
		if self.process.ppgpp_regulation and not self.process.disable_ppgpp_elongation_inhibition:
			cell_mass = self.process.readFromListener("Mass", "cellMass") * units.fg
			cell_volume = cell_mass / self.cellDensity
			counts_to_molar = 1 / (self.process.n_avogadro * cell_volume)
			ppgpp_conc = self.ppgpp.total_count() * counts_to_molar
			rate = self.elong_rate_by_ppgpp(ppgpp_conc, self.basal_elongation_rate).asNumber(units.aa / units.s)
		else:
			rate = super().elongation_rate()
		return rate

	def request(self, aasInSequences):
		self.max_time_step = min(self.process.max_time_step, self.max_time_step * self.time_step_increase)

		# Conversion from counts to molarity
		cell_mass = self.process.readFromListener("Mass", "cellMass") * units.fg
		dry_mass = self.process.readFromListener("Mass", "dryMass") * units.fg
		cell_volume = cell_mass / self.cellDensity
		self.counts_to_molar = 1 / (self.process.n_avogadro * cell_volume)

		# ppGpp related concentrations
		ppgpp_conc = self.counts_to_molar * self.ppgpp.total_count()
		rela_conc = self.counts_to_molar * self.rela.total_count()
		spot_conc = self.counts_to_molar * self.spot.total_count()

		# Get counts and convert synthetase and tRNA to a per AA basis
		synthetase_counts = np.dot(self.aa_from_synthetase, self.synthetases.total_counts())
		aa_counts = self.process.aas.total_counts()
		uncharged_trna_counts = np.dot(self.process.aa_from_trna, self.uncharged_trna.total_counts())
		charged_trna_counts = np.dot(self.process.aa_from_trna, self.charged_trna.total_counts())
		ribosome_counts = self.process.active_ribosomes.total_count()

		# Get concentration
		f = aasInSequences / aasInSequences.sum()
		synthetase_conc = self.counts_to_molar * synthetase_counts
		aa_conc = self.counts_to_molar * aa_counts
		uncharged_trna_conc = self.counts_to_molar * uncharged_trna_counts
		charged_trna_conc = self.counts_to_molar * charged_trna_counts
		ribosome_conc = self.counts_to_molar * ribosome_counts

		# Calculate amino acid supply
		aa_in_media = self.aa_environment.import_present()
		fwd_enzyme_counts, rev_enzyme_counts = self.get_pathway_enzyme_counts_per_aa(
			self.aa_enzymes.total_counts())
		importer_counts = self.aa_importers.total_counts()
		exporter_counts = self.aa_exporters.total_counts()
		synthesis, fwd_saturation, rev_saturation = self.amino_acid_synthesis(fwd_enzyme_counts, rev_enzyme_counts, aa_conc)
		import_rates = self.amino_acid_import(aa_in_media, dry_mass, aa_conc, importer_counts, self.process.mechanistic_aa_transport)
		export_rates = self.amino_acid_export(exporter_counts, aa_conc, self.process.mechanistic_aa_transport)
		exchange_rates = import_rates - export_rates

		supply_function = get_charging_supply_function(
			self.process.aa_supply_in_charging, self.process.mechanistic_translation_supply,
			self.process.mechanistic_aa_transport, self.amino_acid_synthesis,
			self.amino_acid_import, self.amino_acid_export, self.aa_supply_scaling,
			self.counts_to_molar, self.process.aa_supply, fwd_enzyme_counts, rev_enzyme_counts,
			dry_mass, importer_counts, exporter_counts, aa_in_media,
			)

		self.process.writeToListener('GrowthLimits', 'original_aa_supply', self.process.aa_supply)
		self.process.writeToListener('GrowthLimits', 'aa_in_media', aa_in_media)

		# Calculate steady state tRNA levels and resulting elongation rate
		self.charging_params['max_elong_rate'] = self.elongation_rate()
		fraction_charged, v_rib, synthesis_in_charging, import_in_charging, export_in_charging = calculate_trna_charging(
			synthetase_conc,
			uncharged_trna_conc,
			charged_trna_conc,
			aa_conc,
			ribosome_conc,
			f,
			self.charging_params,
			supply=supply_function,
			limit_v_rib=True,
			time_limit=self.process.timeStepSec())

		# Use the supply calculated from each sub timestep while solving the charging steady state
		if self.process.aa_supply_in_charging:
			conversion = 1 / self.counts_to_molar.asNumber(CONC_UNITS) / self.process.timeStepSec()
			synthesis = conversion * synthesis_in_charging
			import_rates = conversion * import_in_charging
			export_rates = conversion * export_in_charging
			self.process.aa_supply = synthesis + import_rates - export_rates
		# Use the supply calculated from the starting amino acid concentrations only
		else:
			if self.process.mechanistic_translation_supply:
				# Set supply based on mechanistic synthesis and supply
				self.process.aa_supply = self.process.timeStepSec() * (synthesis + exchange_rates)
			else:
				# Adjust aa_supply higher if amino acid concentrations are low
				# Improves stability of charging and mimics amino acid synthesis
				# inhibition and export
				self.process.aa_supply *= self.aa_supply_scaling(aa_conc, aa_in_media)
		self.process.aa_exchange_rates = self.counts_to_molar / units.s * (import_rates - export_rates)

		self.process.writeToListener('GrowthLimits', 'synthetase_conc', synthetase_conc.asNumber(CONC_UNITS))
		self.process.writeToListener('GrowthLimits', 'uncharged_trna_conc', uncharged_trna_conc.asNumber(CONC_UNITS))
		self.process.writeToListener('GrowthLimits', 'charged_trna_conc', charged_trna_conc.asNumber(CONC_UNITS))
		self.process.writeToListener('GrowthLimits', 'aa_conc', aa_conc.asNumber(CONC_UNITS))
		self.process.writeToListener('GrowthLimits', 'ribosome_conc', ribosome_conc.asNumber(CONC_UNITS))
		self.process.writeToListener('GrowthLimits', 'fraction_aa_to_elongate', f)

		aa_counts_for_translation = v_rib * f * self.process.timeStepSec() / self.counts_to_molar.asNumber(CONC_UNITS)

		total_trna = self.charged_trna.total_counts() + self.uncharged_trna.total_counts()
		final_charged_trna = stochasticRound(self.process.randomState, np.dot(fraction_charged, self.process.aa_from_trna * total_trna))

		charged_trna_request = self.charged_trna.total_counts() - final_charged_trna
		charged_trna_request[charged_trna_request < 0] = 0
		uncharged_trna_request = final_charged_trna - self.charged_trna.total_counts()
		uncharged_trna_request[uncharged_trna_request < 0] = 0

		self.aa_counts_for_translation = np.array(aa_counts_for_translation)

		fraction_trna_per_aa = total_trna / np.dot(np.dot(self.process.aa_from_trna, total_trna), self.process.aa_from_trna)
		total_charging_reactions = stochasticRound(self.process.randomState,
				np.dot(aa_counts_for_translation, self.process.aa_from_trna)
				* fraction_trna_per_aa + uncharged_trna_request)

		self.process.writeToListener('GrowthLimits', 'aa_supply', self.process.aa_supply)
		self.process.writeToListener('GrowthLimits', 'aa_synthesis', synthesis * self.process.timeStepSec())
		self.process.writeToListener('GrowthLimits', 'aa_import', import_rates * self.process.timeStepSec())
		self.process.writeToListener('GrowthLimits', 'aa_export', export_rates * self.process.timeStepSec())
		self.process.writeToListener('GrowthLimits', 'aa_supply_enzymes_fwd', fwd_enzyme_counts)
		self.process.writeToListener('GrowthLimits', 'aa_supply_enzymes_rev', rev_enzyme_counts)
		self.process.writeToListener('GrowthLimits', 'aa_importers', importer_counts)
		self.process.writeToListener('GrowthLimits', 'aa_exporters', exporter_counts)
		self.process.writeToListener('GrowthLimits', 'aa_supply_aa_conc', aa_conc.asNumber(units.mmol/units.L))
		self.process.writeToListener('GrowthLimits', 'aa_supply_fraction_fwd', fwd_saturation)
		self.process.writeToListener('GrowthLimits', 'aa_supply_fraction_rev', rev_saturation)

		# Only request molecules that will be consumed in the charging reactions
		aa_from_uncharging = -self.charging_stoich_matrix @ charged_trna_request
		aa_from_uncharging[self.charging_molecules_not_aa] = 0
		requested_molecules = -np.dot(self.charging_stoich_matrix, total_charging_reactions) - aa_from_uncharging
		requested_molecules[requested_molecules < 0] = 0
		self.charging_molecules.requestIs(requested_molecules)

		# Request charged tRNA that will become uncharged
		self.charged_trna.requestIs(charged_trna_request)
		self.uncharged_trna_to_charge = uncharged_trna_request

		# Request water for transfer of AA from tRNA for initial polypeptide.
		# This is severe overestimate assuming the worst case that every
		# elongation is initializing a polypeptide. This excess of water
		# shouldn't matter though.
		self.water.requestIs(aa_counts_for_translation.sum())

		# ppGpp reactions based on charged tRNA
		self.process.writeToListener('GrowthLimits', 'ppgpp_conc', ppgpp_conc.asNumber(CONC_UNITS))
		self.process.writeToListener('GrowthLimits', 'rela_conc', rela_conc.asNumber(CONC_UNITS))
		self.process.writeToListener('GrowthLimits', 'spot_conc', spot_conc.asNumber(CONC_UNITS))
		if self.process.ppgpp_regulation:
			total_trna_conc = self.counts_to_molar * (uncharged_trna_counts + charged_trna_counts)
			updated_charged_trna_conc = total_trna_conc * fraction_charged
			updated_uncharged_trna_conc = total_trna_conc - updated_charged_trna_conc
			delta_metabolites, *_ = ppgpp_metabolite_changes(
				updated_uncharged_trna_conc, updated_charged_trna_conc, ribosome_conc,
				f, rela_conc, spot_conc, ppgpp_conc, self.counts_to_molar, v_rib,
				self.charging_params, self.ppgpp_params, self.process.timeStepSec(),
				request=True, random_state=self.process.randomState,
			)

			request_ppgpp_metabolites = -delta_metabolites
			self.ppgpp_reaction_metabolites.requestIs(request_ppgpp_metabolites)
			self.ppgpp.requestAll()

		return fraction_charged, aa_counts_for_translation

	def final_amino_acids(self, total_aa_counts):
		charged_counts_to_uncharge = self.process.aa_from_trna @ self.charged_trna.counts()
		return np.fmin(total_aa_counts + charged_counts_to_uncharge, self.aa_counts_for_translation)

	def evolve(self, total_aa_counts, aas_used, next_amino_acid_count, nElongations, nInitialized):
		# Get tRNA counts
		uncharged_trna = self.uncharged_trna.counts()
		charged_trna = self.charged_trna.counts()
		total_trna = uncharged_trna + charged_trna

		# Adjust molecules for number of charging reactions that occurred
		## Determine limitations for charging and uncharging reactions
		charged_and_elongated_per_aa = np.fmax(0, (aas_used - self.process.aa_from_trna @ charged_trna))
		aa_for_charging = total_aa_counts - charged_and_elongated_per_aa
		n_aa_charged = np.fmin(aa_for_charging, np.dot(self.process.aa_from_trna, np.fmin(self.uncharged_trna_to_charge, uncharged_trna)))
		n_uncharged_per_aa = aas_used - charged_and_elongated_per_aa

		## Calculate changes in tRNA based on limitations
		n_trna_charged = self.distribution_from_aa(n_aa_charged, uncharged_trna, True)
		n_trna_uncharged = self.distribution_from_aa(n_uncharged_per_aa, charged_trna, True)

		## Determine reactions that are charged and elongated in same time step without changing
		## charged or uncharged counts
		charged_and_elongated = self.distribution_from_aa(charged_and_elongated_per_aa, total_trna)

		## Determine total number of reactions that occur
		total_uncharging_reactions = charged_and_elongated + n_trna_uncharged
		total_charging_reactions = charged_and_elongated + n_trna_charged
		net_charged = total_charging_reactions - total_uncharging_reactions
		self.charging_molecules.countsInc(np.dot(self.charging_stoich_matrix, total_charging_reactions))

		## Account for uncharging of tRNA during elongation
		self.charged_trna.countsDec(total_uncharging_reactions)
		self.uncharged_trna.countsInc(total_uncharging_reactions)

		# Update proton counts to reflect polymerization reactions and transfer of AA from tRNA
		# Peptide bond formation releases a water but transferring AA from tRNA consumes a OH-
		# Net production of H+ for each elongation, consume extra water for each initialization
		# since a peptide bond doesn't form
		self.proton.countInc(nElongations)
		self.water.countDec(nInitialized)

		# Create or degrade ppGpp
		# This should come after all countInc/countDec calls since it shares some molecules with
		# other views and those counts should be updated to get the proper limits on ppGpp reactions
		if self.process.ppgpp_regulation:
			v_rib = (nElongations * self.counts_to_molar).asNumber(CONC_UNITS) / self.process.timeStepSec()
			ribosome_conc = self.counts_to_molar * self.process.active_ribosomes.total_count()
			updated_uncharged_trna_counts = self.uncharged_trna.total_counts() - net_charged
			updated_charged_trna_counts = self.charged_trna.total_counts() + net_charged
			uncharged_trna_conc = self.counts_to_molar * np.dot(
				self.process.aa_from_trna, updated_uncharged_trna_counts)
			charged_trna_conc = self.counts_to_molar * np.dot(
				self.process.aa_from_trna, updated_charged_trna_counts)
			ppgpp_conc = self.counts_to_molar * self.ppgpp.total_count()
			rela_conc = self.counts_to_molar * self.rela.total_count()
			spot_conc = self.counts_to_molar * self.spot.total_count()

			# Need to include the next amino acid the ribosome sees for certain
			# cases where elongation does not occur, otherwise f will be NaN
			aa_at_ribosome = aas_used + next_amino_acid_count
			f = aa_at_ribosome / aa_at_ribosome.sum()
			limits = self.ppgpp_reaction_metabolites.counts()
			delta_metabolites, ppgpp_syn, ppgpp_deg, rela_syn, spot_syn, spot_deg, spot_deg_inhibited = ppgpp_metabolite_changes(
				uncharged_trna_conc, charged_trna_conc,	ribosome_conc, f, rela_conc,
				spot_conc, ppgpp_conc, self.counts_to_molar, v_rib,
				self.charging_params, self.ppgpp_params, self.process.timeStepSec(),
				random_state=self.process.randomState, limits=limits,
				)

			self.process.writeToListener('GrowthLimits', 'rela_syn', rela_syn)
			self.process.writeToListener('GrowthLimits', 'spot_syn', spot_syn)
			self.process.writeToListener('GrowthLimits', 'spot_deg', spot_deg)
			self.process.writeToListener('GrowthLimits', 'spot_deg_inhibited', spot_deg_inhibited)

			self.ppgpp_reaction_metabolites.countsInc(delta_metabolites)

		# Use the difference between (expected AA supply based on expected doubling time
		# and current DCW) and AA used to charge tRNA to update the concentration target
		# in metabolism during the next time step
		aa_used_trna = np.dot(self.process.aa_from_trna, total_charging_reactions)
		aa_diff = self.process.aa_supply - aa_used_trna
		if np.any(np.abs(aa_diff / self.process.aas.total_counts()) > self.max_amino_acid_adjustment):
			self.time_step_short_enough = False

		self.process.writeToListener('GrowthLimits', 'trnaCharged', aa_used_trna)
		return net_charged, {aa: diff for aa, diff in zip(self.aaNames, aa_diff)}

	def distribution_from_aa(self, n_aa, n_trna, limited=False):
		'''
		Distributes counts of amino acids to tRNAs that are associated with each amino acid.
		Uses self.process.aa_from_trna mapping to distribute from amino acids to tRNA based on the
		fraction that each tRNA species makes up for all tRNA species that code for the
		same amino acid.

		Inputs:
			n_aa (array of ints) - counts of each amino acid to distribute to each tRNA
			n_trna (array of ints) - counts of each tRNA to determine the distribution
			limited (bool) - optional, if True, limits the amino acids distributed to
				each tRNA to the number of tRNA that are available (n_trna)

		Returns:
			array of ints - distributed counts for each tRNA
		'''

		# Determine the fraction each tRNA species makes up out of all tRNA of the
		# associated amino acid
		with np.errstate(invalid='ignore'):
			f_trna = n_trna / np.dot(np.dot(self.process.aa_from_trna, n_trna), self.process.aa_from_trna)
		f_trna[~np.isfinite(f_trna)] = 0

		trna_counts = np.zeros(f_trna.shape, np.int64)
		for count, row in zip(n_aa, self.process.aa_from_trna):
			idx = (row == 1)
			frac = f_trna[idx]

			counts = np.floor(frac * count)
			diff = int(count - counts.sum())

			# Add additional counts to get up to counts to distribute
			# Prevent adding over the number of tRNA available if limited
			if diff > 0:
				if limited:
					for _ in range(diff):
						frac[(n_trna[idx] - counts) == 0] = 0
						frac /= frac.sum()  # normalize for multinomial distribution
						adjustment = self.process.randomState.multinomial(1, frac)
						counts += adjustment
				else:
					adjustment = self.process.randomState.multinomial(diff, frac)
					counts += adjustment

			trna_counts[idx] = counts

		return trna_counts

	def isTimeStepShortEnough(self, inputTimeStep, timeStepSafetyFraction):
		short_enough = True

		# Needs to be less than the max time step to prevent oscillatory behavior
		if inputTimeStep > self.max_time_step:
			short_enough = False

		# Decrease the max time step to get more stable charging
		if not self.time_step_short_enough and self.process.adjust_timestep_for_charging:
			self.max_time_step = inputTimeStep / 2
			self.time_step_short_enough = True
			short_enough = False

		return short_enough

def get_ppgpp_params(sim_data) -> Dict[str, Any]:
	"""
	Get parameters required for ppGpp reaction calulations to help
	encapsulate the function so that it does not need to be a class method.

	Args:
		sim_data: SimulationData object

	Returns:
		parameters that are used in ppgpp_metabolite_changes
	"""

	constants = sim_data.constants
	metabolism = sim_data.process.metabolism
	transcription = sim_data.process.transcription

	return dict(
		KD_RelA=transcription.KD_RelA.asNumber(CONC_UNITS),
		k_RelA=constants.k_RelA_ppGpp_synthesis.asNumber(1 / units.s),
		k_SpoT_syn=constants.k_SpoT_ppGpp_synthesis.asNumber(1 / units.s),
		k_SpoT_deg=constants.k_SpoT_ppGpp_degradation.asNumber(1 / (CONC_UNITS * units.s)),
		KI_SpoT=transcription.KI_SpoT.asNumber(CONC_UNITS),
		ppgpp_reaction_stoich=metabolism.ppgpp_reaction_stoich,
		synthesis_index=metabolism.ppgpp_reaction_names.index(metabolism.ppgpp_synthesis_reaction),
		degradation_index=metabolism.ppgpp_reaction_names.index(metabolism.ppgpp_degradation_reaction),
		)

def ppgpp_metabolite_changes(uncharged_trna_conc, charged_trna_conc,
		ribosome_conc, f, rela_conc, spot_conc, ppgpp_conc, counts_to_molar,
		v_rib, charging_params, ppgpp_params, time_step,
		request=False, limits=None, random_state=None):
	'''
	Calculates the changes in metabolite counts based on ppGpp synthesis and
	degradation reactions.

	Args:
		uncharged_trna_conc (np.array[float] with concentration units):
			concentration of uncharged tRNA associated with each amino acid
		charged_trna_conc (np.array[float] with concentration units):
			concentration of charged tRNA associated with each amino acid
		ribosome_conc (float with concentration units): concentration of active ribosomes
		f (np.array[float]): fraction of each amino acid to be incorporated
			to total amino acids incorporated
		rela_conc (float with concentration units): concentration of RelA
		spot_conc (float with concentration units): concentration of SpoT
		ppgpp_conc (float with concentration units): concentration of ppGpp
		counts_to_molar (float with concentration units): conversion factor
			from counts to molarity
		v_rib (float): rate of amino acid incorporation at the ribosome,
			in units of uM/s
		charging_params (Dict[str, Any]): parameters used in charging equations
			- this should be generated by get_charging_params
		ppgpp_params (Dict[str, Any]): parameters used in ppGpp reactions
			- this should be generated by get_ppgpp_params
		time_step (float): length of the current time step
		request (bool): if True, only considers reactant stoichiometry,
			otherwise considers reactants and products. For use in
			calculateRequest. GDP appears as both a reactant and product
			and the request can be off the actual use if not handled in this
			manner.
		limits (np.array[float]): counts of molecules that are available to prevent
			negative total counts as a result of delta_metabolites.
			If None, no limits are placed on molecule changes.
		random_state (np.random.RandomState): random state for the process

	Returns:
		delta_metabolites (np.array[int]): the change in counts of each metabolite
			involved in ppGpp reactions
		n_syn_reactions (int): the number of ppGpp synthesis reactions
		n_deg_reactions (int): the number of ppGpp degradation reactions
		v_rela_syn (np.ndarray[float]): rate of synthesis from RelA per amino
			acid tRNA species
		v_spot_syn (float): rate of synthesis from SpoT
		v_deg (float): rate of degradation from SpoT
		v_deg_inhibited (np.ndarray[float]): rate of degradation from SpoT per
			amino acid tRNA species
	'''

	if random_state is None:
		random_state = np.random.RandomState()

	uncharged_trna_conc = uncharged_trna_conc.asNumber(CONC_UNITS)
	charged_trna_conc = charged_trna_conc.asNumber(CONC_UNITS)
	ribosome_conc = ribosome_conc.asNumber(CONC_UNITS)
	rela_conc = rela_conc.asNumber(CONC_UNITS)
	spot_conc = spot_conc.asNumber(CONC_UNITS)
	ppgpp_conc = ppgpp_conc.asNumber(CONC_UNITS)
	counts_to_micromolar = counts_to_molar.asNumber(CONC_UNITS)

	numerator = 1 + charged_trna_conc / charging_params['krta'] + uncharged_trna_conc / charging_params['krtf']
	saturated_charged = charged_trna_conc / charging_params['krta'] / numerator
	saturated_uncharged = uncharged_trna_conc / charging_params['krtf'] / numerator
	if v_rib == 0:
		ribosome_conc_a_site = f * ribosome_conc
	else:
		ribosome_conc_a_site = f * v_rib / (saturated_charged * charging_params['max_elong_rate'])
	ribosomes_bound_to_uncharged = ribosome_conc_a_site * saturated_uncharged

	# Handle rare cases when tRNA concentrations are 0
	# Can result in inf and nan so assume a fraction of ribosomes
	# bind to the uncharged tRNA if any tRNA are present or 0 if not
	mask = ~np.isfinite(ribosomes_bound_to_uncharged)
	ribosomes_bound_to_uncharged[mask] = ribosome_conc * f[mask] * np.array(
		uncharged_trna_conc[mask] + charged_trna_conc[mask] > 0)

	# Calculate active fraction of RelA
	competitive_inhibition = 1 + ribosomes_bound_to_uncharged / ppgpp_params['KD_RelA']
	inhibition_product = np.prod(competitive_inhibition)
	with np.errstate(divide='ignore'):
		frac_rela = 1 / (ppgpp_params['KD_RelA'] / ribosomes_bound_to_uncharged * inhibition_product / competitive_inhibition + 1)

	# Calculate rates for synthesis and degradation
	v_rela_syn = ppgpp_params['k_RelA'] * rela_conc * frac_rela
	v_spot_syn = ppgpp_params['k_SpoT_syn'] * spot_conc
	v_syn = v_rela_syn.sum() + v_spot_syn
	max_deg = ppgpp_params['k_SpoT_deg'] * spot_conc * ppgpp_conc
	fractions = uncharged_trna_conc / ppgpp_params['KI_SpoT']
	v_deg =  max_deg / (1 + fractions.sum())
	v_deg_inhibited = (max_deg - v_deg) * fractions / fractions.sum()

	# Convert to discrete reactions
	n_syn_reactions = stochasticRound(random_state, v_syn * time_step / counts_to_micromolar)[0]
	n_deg_reactions = stochasticRound(random_state, v_deg * time_step / counts_to_micromolar)[0]

	# Only look at reactant stoichiometry if requesting molecules to use
	if request:
		ppgpp_reaction_stoich = np.zeros_like(ppgpp_params['ppgpp_reaction_stoich'])
		reactants = ppgpp_params['ppgpp_reaction_stoich'] < 0
		ppgpp_reaction_stoich[reactants] = ppgpp_params['ppgpp_reaction_stoich'][reactants]
	else:
		ppgpp_reaction_stoich = ppgpp_params['ppgpp_reaction_stoich']

	# Calculate the change in metabolites and adjust to limits if provided
	# Possible reactions are adjusted down to limits if the change in any
	# metabolites would result in negative counts
	max_iterations = int(n_deg_reactions + n_syn_reactions + 1)
	old_counts = None
	for it in range(max_iterations):
		delta_metabolites = (ppgpp_reaction_stoich[:, ppgpp_params['synthesis_index']] * n_syn_reactions
			+ ppgpp_reaction_stoich[:, ppgpp_params['degradation_index']] * n_deg_reactions)

		if limits is None:
			break
		else:
			final_counts = delta_metabolites + limits

			if np.all(final_counts >= 0) or (old_counts is not None and np.all(final_counts == old_counts)):
				break

			limited_index = np.argmin(final_counts)
			if ppgpp_reaction_stoich[limited_index, ppgpp_params['synthesis_index']] < 0:
				limited = np.ceil(final_counts[limited_index] / ppgpp_reaction_stoich[limited_index, ppgpp_params['synthesis_index']])
				n_syn_reactions -= min(limited, n_syn_reactions)
			if ppgpp_reaction_stoich[limited_index, ppgpp_params['degradation_index']] < 0:
				limited = np.ceil(final_counts[limited_index] / ppgpp_reaction_stoich[limited_index, ppgpp_params['degradation_index']])
				n_deg_reactions -= min(limited, n_deg_reactions)

			old_counts = final_counts
	else:
		raise ValueError('Failed to meet molecule limits with ppGpp reactions.')

	return delta_metabolites, n_syn_reactions, n_deg_reactions, v_rela_syn, v_spot_syn, v_deg, v_deg_inhibited

def get_charging_params(
		sim_data,
		aa_removed_from_charging: Optional[Set[str]] = None,
		variable_elongation: bool = False,
		) -> Dict[str, Any]:
	"""
	Get parameters required for tRNA charging calulations to help
	encapsulate the function so that it does not need to be a class method.

	Args:
		sim_data: SimulationData object
		aa_removed_from_charging: any amino acid IDs that should be ignored
			when calculating charging
		variable_elongation: if True, the max elongation rate is set to be
			higher

	Returns:
		parameters that are used in calculate_trna_charging
	"""

	constants = sim_data.constants
	metabolism = sim_data.process.metabolism
	transcription = sim_data.process.transcription
	if aa_removed_from_charging is None:
		aa_removed_from_charging = REMOVED_FROM_CHARGING
	aa_charging_mask = np.array([
		aa not in aa_removed_from_charging
		for aa in sim_data.molecule_groups.amino_acids
		])
	elongation_max = (constants.ribosome_elongation_rate_max
		if variable_elongation else constants.ribosome_elongation_rate_basal)

	return dict(
		kS=constants.synthetase_charging_rate.asNumber(1 / units.s),
		KMaa=transcription.aa_kms.asNumber(CONC_UNITS),
		KMtf=transcription.trna_kms.asNumber(CONC_UNITS),
		krta=constants.Kdissociation_charged_trna_ribosome.asNumber(CONC_UNITS),
		krtf=constants.Kdissociation_uncharged_trna_ribosome.asNumber(CONC_UNITS),
		max_elong_rate=float(elongation_max.asNumber(units.aa / units.s)),
		charging_mask=aa_charging_mask,
		unit_conversion=metabolism.get_amino_acid_conc_conversion(CONC_UNITS),
		)

def calculate_trna_charging(synthetase_conc, uncharged_trna_conc, charged_trna_conc, aa_conc, ribosome_conc,
		f, params, supply=None, time_limit=1000, limit_v_rib=False, use_disabled_aas=False):
	'''
	Calculates the steady state value of tRNA based on charging and incorporation through polypeptide elongation.
	The fraction of charged/uncharged is also used to determine how quickly the ribosome is elongating.

	Inputs:
		synthetase_conc (array of floats with concentration units) - concentration of synthetases associated
			with each amino acid
		uncharged_trna_conc (array of floats with concentration units) - concentration of uncharged tRNA associated
			with each amino acid
		charged_trna_conc (array of floats with concentration units) - concentration of charged tRNA associated
			with each amino acid
		aa_conc (array of floats with concentration units) - concentration of each amino acid
		ribosome_conc (float with concentration units) - concentration of active ribosomes
		f (array of floats) - fraction of each amino acid to be incorporated to total amino acids incorporated
		params (Dict[str, Any]) - parameters used in charging equations - this should be
			generated by get_charging_params
		supply (Callable) - function to get the rate of amino acid supply (synthesis and import)
			based on amino acid concentrations. If None, amino acid concentrations remain constant
			during charging
		time_limit (float) - time limit to reach steady state
		limit_v_rib (bool) - if True, v_rib is limited to the number of amino acids that are
			available
		use_disabled_aas (bool) - if True, all amino acids will be used for charging calculations,
			if False, some will be excluded as determined in initialize

	Returns:
		new_fraction_charged (array of floats) - fraction of total tRNA that is charged for each
			amino acid species
		v_rib (float) - ribosomal elongation rate in units of uM/s
		total_synthesis (np.ndarray) - the total amount of amino acids synthesized during charging
			in units of CONC_UNITS.  Will be zeros if supply function is not given.
		total_import (np.ndarray) - the total amount of amino acids imported during charging
			in units of CONC_UNITS.  Will be zeros if supply function is not given.
		total_export (np.ndarray) - the total amount of amino acids exported during charging
			in units of CONC_UNITS.  Will be zeros if supply function is not given.
	'''

	def negative_check(trna1, trna2):
		'''
		Check for floating point precision issues that can lead to small
		negative numbers instead of 0. Adjusts both species of tRNA to
		bring concentration of trna1 to 0 and keep the same total concentration.

		Args:
			trna1 (ndarray[float]): concentration of one tRNA species (charged or uncharged)
			trna2 (ndarray[float]): concentration of another tRNA species (charged or uncharged)
		'''

		mask = trna1 < 0
		trna2[mask] = trna1[mask] + trna2[mask]
		trna1[mask] = 0

	def dcdt(t, c):
		'''
		Function for solve_ivp to integrate

		Args:
			c (ndarray[float]): 1D array of concentrations of uncharged and charged tRNAs
				dims: 2 * number of amino acids (uncharged tRNA come first, then charged)
			t (float): time of integration step

		Returns:
			ndarray[float]: dc/dt for tRNA concentrations
				dims: 2 * number of amino acids (uncharged tRNA come first, then charged)
		'''
		v_charging, dtrna, daa = dcdt_jit(t, c, n_aas_masked, n_aas, mask, 
			params['kS'], synthetase_conc, params['KMaa'], params['KMtf'],
			f, params['krta'], params['krtf'], params['max_elong_rate'],
			ribosome_conc, limit_v_rib, aa_rate_limit, v_rib_max)
		if supply is None:
			v_synthesis = np.zeros(n_aas)
			v_import = np.zeros(n_aas)
			v_export = np.zeros(n_aas)
		else:
			aa_conc = c[2*n_aas_masked:2*n_aas_masked+n_aas]
			v_synthesis, v_import, v_export = supply(unit_conversion * aa_conc)
			v_supply = v_synthesis + v_import - v_export
			daa[mask] = v_supply[mask] - v_charging

		return np.hstack((-dtrna, dtrna, daa, v_synthesis, v_import, v_export))

	# Convert inputs for integration
	synthetase_conc = synthetase_conc.asNumber(CONC_UNITS)
	uncharged_trna_conc = uncharged_trna_conc.asNumber(CONC_UNITS)
	charged_trna_conc = charged_trna_conc.asNumber(CONC_UNITS)
	aa_conc = aa_conc.asNumber(CONC_UNITS)
	ribosome_conc = ribosome_conc.asNumber(CONC_UNITS)
	unit_conversion = params['unit_conversion']

	# Remove disabled amino acids from calculations
	n_total_aas = len(aa_conc)
	if use_disabled_aas:
		mask = np.ones(n_total_aas, bool)
	else:
		mask = params['charging_mask']
	synthetase_conc = synthetase_conc[mask]
	original_uncharged_trna_conc = uncharged_trna_conc[mask]
	original_charged_trna_conc = charged_trna_conc[mask]
	original_aa_conc = aa_conc[mask]
	f = f[mask]

	n_aas = len(aa_conc)
	n_aas_masked = len(original_aa_conc)

	# Limits for integration
	aa_rate_limit = original_aa_conc / time_limit
	trna_rate_limit = original_charged_trna_conc / time_limit
	v_rib_max = max(0, ((aa_rate_limit + trna_rate_limit) / f).min())

	# Integrate rates of charging and elongation
	c_init = np.hstack((original_uncharged_trna_conc, original_charged_trna_conc,
		aa_conc, np.zeros(n_aas), np.zeros(n_aas), np.zeros(n_aas)))
	sol = solve_ivp(dcdt, [0, time_limit], c_init, method='BDF')
	c_sol = sol.y.T

	# Determine new values from integration results
	final_uncharged_trna_conc = c_sol[-1, :n_aas_masked]
	final_charged_trna_conc = c_sol[-1, n_aas_masked:2*n_aas_masked]
	total_synthesis = c_sol[-1, 2*n_aas_masked+n_aas:2*n_aas_masked+2*n_aas]
	total_import = c_sol[-1, 2*n_aas_masked+2*n_aas:2*n_aas_masked+3*n_aas]
	total_export = c_sol[-1, 2*n_aas_masked+3*n_aas:2*n_aas_masked+4*n_aas]

	negative_check(final_uncharged_trna_conc, final_charged_trna_conc)
	negative_check(final_charged_trna_conc, final_uncharged_trna_conc)

	fraction_charged = final_charged_trna_conc / (final_uncharged_trna_conc + final_charged_trna_conc)
	numerator_ribosome = 1 + np.sum(f * (params['krta'] / final_charged_trna_conc + final_uncharged_trna_conc / final_charged_trna_conc * params['krta'] / params['krtf']))
	v_rib = params['max_elong_rate'] * ribosome_conc / numerator_ribosome
	if limit_v_rib:
		v_rib_max = max(0, ((original_aa_conc + (original_charged_trna_conc - final_charged_trna_conc)) / time_limit / f).min())
		v_rib = min(v_rib, v_rib_max)

	# Replace SEL fraction charged with average
	new_fraction_charged = np.zeros(n_total_aas)
	new_fraction_charged[mask] = fraction_charged
	new_fraction_charged[~mask] = fraction_charged.mean()

	return new_fraction_charged, v_rib, total_synthesis, total_import, total_export

@njit(error_model='numpy')
def dcdt_jit(t, c, n_aas_masked, n_aas, mask, 
	kS, synthetase_conc, KMaa, KMtf,
	f, krta, krtf, max_elong_rate,
	ribosome_conc, limit_v_rib, aa_rate_limit, v_rib_max
):
	uncharged_trna_conc = c[:n_aas_masked]
	charged_trna_conc = c[n_aas_masked:2*n_aas_masked]
	aa_conc = c[2*n_aas_masked:2*n_aas_masked+n_aas]
	masked_aa_conc = aa_conc[mask]

	v_charging = (kS * synthetase_conc * uncharged_trna_conc * masked_aa_conc / (KMaa[mask] * KMtf[mask])
		/ (1 + uncharged_trna_conc/KMtf[mask] + masked_aa_conc/KMaa[mask] + uncharged_trna_conc*masked_aa_conc/KMtf[mask]/KMaa[mask]))
	numerator_ribosome = 1 + np.sum(f * (krta / charged_trna_conc + uncharged_trna_conc / charged_trna_conc * krta / krtf))
	v_rib = max_elong_rate * ribosome_conc / numerator_ribosome

	# Handle case when f is 0 and charged_trna_conc is 0
	if not np.isfinite(v_rib):
		v_rib = 0

	# Limit v_rib and v_charging to the amount of available amino acids
	if limit_v_rib:
		v_charging = np.fmin(v_charging, aa_rate_limit)
		v_rib = min(v_rib, v_rib_max)

	dtrna = v_charging - v_rib*f
	daa = np.zeros(n_aas)
	
	return v_charging, dtrna, daa

def get_charging_supply_function(
		supply_in_charging: bool,
		mechanistic_supply: bool,
		mechanistic_aa_transport: bool,
		amino_acid_synthesis: Callable,
		amino_acid_import: Callable,
		amino_acid_export: Callable,
		aa_supply_scaling: Callable,
		counts_to_molar: units.Unum,
		aa_supply: np.ndarray,
		fwd_enzyme_counts: np.ndarray,
		rev_enzyme_counts: np.ndarray,
		dry_mass: units.Unum,
		importer_counts: np.ndarray,
		exporter_counts: np.ndarray,
		aa_in_media: np.ndarray,
		) -> Optional[Callable[[np.ndarray], Tuple[np.ndarray, np.ndarray, np.ndarray]]]:
	"""
	Get a function mapping internal amino acid concentrations to the amount of
	amino acid supply expected.

	Args:
		supply_in_charging: True if using the aa_supply_in_charging option
		mechanistic_supply: True if using the mechanistic_translation_supply option
		mechanistic_aa_transport: True if using the mechanistic_aa_transport option
		amino_acid_synthesis: function to provide rates of synthesis for amino
			acids based on the internal state
		amino_acid_import: function to provide import rates for amino
			acids based on the internal and external state
		amino_acid_export: function to provide export rates for amino
			acids based on the internal state
		aa_supply_scaling: function to scale the amino acid supply based
			on the internal state
		counts_to_molar: conversion factor for counts to molar in units of counts/volume
		aa_supply: rate of amino acid supply expected
		fwd_enzyme_counts: counts for enzymes in forward reactions for each amino acid
		rev_enzyme_counts: counts for enzymes in loss reactions for each amino acid
		dry_mass: dry mass of the cell with mass units
		importer_counts: counts for amino acid importers
		exporter_counts: counts for amino acid exporters
		aa_in_media: True for each amino acid that is present in the media

	Returns:
		supply_function: function that provides the amount of supply (synthesis, import, export)
			for each amino acid based on the internal state of the cell
	"""

	# Create functions that are only dependent on amino acid concentrations for more stable
	# charging and amino acid concentrations.  If supply_in_charging is not set, then
	# setting None will maintain constant amino acid concentrations throughout charging.
	supply_function = None
	if supply_in_charging:
		counts_to_molar = counts_to_molar.asNumber(CONC_UNITS)
		zeros = counts_to_molar * np.zeros_like(aa_supply)
		if mechanistic_supply:
			if mechanistic_aa_transport:
				supply_function = lambda aa_conc: (
					counts_to_molar * amino_acid_synthesis(fwd_enzyme_counts, rev_enzyme_counts, aa_conc)[0],
					counts_to_molar * amino_acid_import(aa_in_media, dry_mass, aa_conc, importer_counts, mechanistic_aa_transport),
					counts_to_molar * amino_acid_export(exporter_counts, aa_conc, mechanistic_aa_transport),
					)
			else:
				supply_function = lambda aa_conc: (
					counts_to_molar * amino_acid_synthesis(fwd_enzyme_counts, rev_enzyme_counts, aa_conc)[0],
					counts_to_molar * amino_acid_import(aa_in_media, dry_mass, aa_conc, importer_counts, mechanistic_aa_transport),
					zeros,
					)
		else:
			supply_function = lambda aa_conc: (
				counts_to_molar * aa_supply * aa_supply_scaling(aa_conc, aa_in_media),
				zeros,
				zeros,
				)

	return supply_function
