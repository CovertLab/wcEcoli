"""
Optional ParCa step that sets tRNA charging from measured data and makes
charging, attenuation, amino acid supply and ppGpp consistent with each other
before the simulation starts. Port of the vivarium-ecoli (vEcoli) step of the
same name.

Without this step every synthetase shares one rate constant (kS = 100 /s), so
charging capacity ranges from about 3x to 40x demand between amino acids, and
the ParCa calibrates attenuation as if all tRNA were charged. Genome-wide
measurements put tRNA charging in glucose minimal medium near 60%
(Avcilar-Kucukgoze et al. 2016, which uses E. coli MC4100 (different from
MG1655, especially relA) but at a 43 min doubling time, close to the 44 min
ParCa target; a tRNA-Phe Northern in MG1655 gives about 65%, Bullwinkle and
Ibba 2016). At that charging, the literature RelA constants (Bosdriesz et al.
2015, fit alongside kS near full charging) would drive ppGpp far above its
target, so they are refit to the ParCa's context here.

:py:func:`fit_charging_attenuation_and_ppgpp` does the following:

1. Fixed-point loop:

   a. Fit amino acid supply constants (synthesis, export, uptake) and the
      ppGpp kinetics spread from the current expected bulk counts (these
      depend on attenuation through the readthrough of attenuated operons).
   b. Per-synthetase charging rate constants are set once, after the first
      supply fit, so each synthetase has its measured capacity (Dykeman 2025
      kcat x copies per cell) on the model's own counts. (Optionally, one
      shared kS can be fit to a measured mean charged fraction instead; see
      DEFAULT_ANCHOR_CHARGED_FRACTION for why that is off by default.)
   c. Solve the steady-state tRNA charging and free amino acid concentrations
      for each fitted condition with the simulation's own rate laws
      (:py:func:`models.ecoli.processes.polypeptide_elongation.dcdt_jit` for
      charging and elongation, and the mechanistic synthesis, import and
      export functions for supply). Amino acids released by protein
      degradation are included in supply, matching the simulation, which
      adds them to the amino acid targets.
   d. Recalibrate attenuation on that charged tRNA and repeat until the
      charged fractions stop changing.

2. ppGpp kinetics, with RelA's maximum rate (k_RelA) fixed:

   a. The anchor condition's SpoT synthesis factor is set so SpoT alone holds
      ppGpp at the measured relA-null level (Murray and Bremer, 1996) at that
      condition's charging. The level is scaled with the condition's target
      ppGpp when its doubling time differs from the default basal 44 min.
   b. KD_RelA (one scale on the per amino acid values) is fit so SpoT plus
      RelA holds the anchor condition at its target ppGpp.
   c. Every other fitted condition's SpoT factor is set from its own target
      ppGpp at its predicted charging (Murray and Bremer's principle that
      SpoT synthesis sets the difference between media, so the value is the
      model's prediction).
"""

from typing import Any, Dict, Optional, Tuple

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

from reconstruction.ecoli.initialization import create_expected_bulk_container
from wholecell.utils import units

UM = units.umol / units.L

# Conditions whose charging is solved (and used for attenuation)
DEFAULT_FIT_CONDITIONS = ('basal', 'with_aa')
# Condition whose measured charging anchors kS, RelA and its SpoT factor
ANCHOR_CONDITION = 'basal'
# Measured mean charged fraction in glucose minimal medium (Avcilar-Kucukgoze
# et al. 2016: "oscillates around 50-60%", MC4100, which carries relA1; a
# tRNA-Phe Northern in MG1655 M9 gives about 65%). Not used by default: with
# tRNA KMs (0.16-1.8 uM) far below tRNA levels (1-24 uM per amino acid) the
# charging rate law saturates in uncharged tRNA, so charging is either nearly
# full or collapses when capacity falls below demand. Stable intermediate
# charging needs charged tRNA inhibition (DEFAULT_CHARGED_TRNA_KI_OVER_KM).
MEASURED_MINIMAL_CHARGED_FRACTION = 0.60
DEFAULT_ANCHOR_CHARGED_FRACTION = None  # type: Optional[float]
# Conditions whose growth rates set the per-synthetase rate constants from
# Dykeman 2025 capacities (kcat x copies per cell)
SYNTHETASE_CAPACITY_CONDITIONS = ('basal', 'with_aa')
# Headroom floor: each synthetase's charging capacity is raised to at least
# MEASURED_HEADROOM_MULTIPLIER x its measured margin (Dykeman 2025 S1 Tab E,
# trna_synthetase_headroom.tsv, interpolated to each fitted condition) over the
# model's own demand. Without a floor, ThrRS sits at about 1.2x in rich (on the
# r-protein operon, it rises less with growth rate than Dykeman's copies),
# leaving Thr about 65% uncharged. No floor is applied when the table is missing.
# The measured margin is capacity at saturating substrates. True applies it at
# the model's own amino acid pools (divides by the KM saturation), so
# low-saturation amino acids (Cys, Trp) keep the measured margin in practice;
# False applies it as if pools were saturated (with KI/KM 10, Cys went to 74%
# uncharged and basal and rich sims grew 5-7 min slower).
MEASURED_HEADROOM_AT_MODEL_POOLS = True
# Measured margins replaced by the median of the other synthetases' margins.
# AlaRS (1.14x, the thinnest) is the least constrained value in Dykeman 2025:
# its turnover was lowered to about 10 /s because single-turnover experiments
# put the transfer step at 16.4 /s (S1 Text), and it is the largest kcat
# discrepancy against in vitro data (6.6x), which the main text says needs
# revisiting. At 1.14x, Ala went to 56% uncharged in with_aa sims.
UNRELIABLE_HEADROOM = ('L-ALPHA-ALANINE',)
# One factor on every measured margin. Dykeman's usage is sized for
# post-upshift elongation (margins about 1.09x larger against steady-state
# demand), and macromolecular crowding can raise effective kcat by up to 1.6x
# (Weilandt and Hatzimanikatis, 2019, cited by Dykeman), so about 1.1-1.75 is
# defensible. An offline grid of this factor and the charged-tRNA KI/KM put
# 1.4 (with KI/KM 3) as the smallest value at which no synthetase dropped
# 20% raises its amino acid's uncharged fraction by more than 0.2. In sims
# (4 seeds, gens 3-10), 1.4 left an amino acid over 70% uncharged in 16%
# (basal) and 31% (with_aa) of generations, with ppGpp swings as large as
# the former uniform 2.2x floor; 1.6 cut this to 9% and 3%, halved the
# ppGpp spread and grew at 44.5 and 24.2 min. The factor and KI/KM
# compensate for each other; the charging data alone do not pin either one.
MEASURED_HEADROOM_MULTIPLIER = 1.6
# Inhibition of synthetases by their own charged tRNA, as a multiple of the
# uncharged tRNA KM (1 = charged and uncharged tRNA bind equally; Dykeman
# 2025's mechanistic scheme includes charged-tRNA-bound synthetase states, and
# for class I synthetases release of the charged tRNA is often rate limiting).
# With tRNA levels far above the KMs, this makes the uncharged fraction track
# demand/capacity instead of switching between nearly full charging and
# collapse. None turns it off.
# One shared value for all synthetases. Not fit in ParCa: with the headroom
# floor, KI/KM and MEASURED_HEADROOM_MULTIPLIER trade off against each other,
# and 5-10 fit per amino acid charging compiled by Dykeman 2025 (S1 Tab I)
# about as well as 3. 3 is used with the 1.6x multiplier.
DEFAULT_CHARGED_TRNA_KI_OVER_KM = 3.0  # type: Optional[float]
# relA-null ppGpp in the anchor condition: Murray and Bremer, 1996, glucose
# minimal, 24 pmol/OD460, about 46 uM with 1.9 uM per pmol/OD460
# (Ahn-Horst et al. 2022)
DEFAULT_RELA_NULL_PPGPP_UM = 46.0
# Basal doubling time (min) that the 46 uM relA-null level is paired with (the
# model's default basal condition). When the anchor condition's doubling time
# differs, the relA-null level is scaled by its target ppGpp over the target
# at this doubling time, so SpoT keeps the same share of basal ppGpp (about
# 94%). Holding 46 uM at a 63 min basal doubling time made the fit tighten
# KD_RelA about 11x and let RelA make 30-60% of basal ppGpp.
RELA_NULL_REFERENCE_DOUBLING_TIME_MIN = 44.0
MAX_LOOP_ITERATIONS = 20
# Fraction of the way each pass moves the charged fractions used for
# attenuation toward the new charging solve (1 = no damping). Damping only
# changes the path, not a converged answer.
ATTENUATION_DAMPING = 0.5
# Largest change in any amino acid's charged fraction to call the loop converged
CHARGED_FRACTION_TOLERANCE = 2e-3
# Integration time per round of the steady-state solve (s) and rounds allowed
STEADY_STATE_TIME = 2000.0
MAX_STEADY_STATE_ROUNDS = 10
STEADY_STATE_TOLERANCE = 1e-4
# Search ranges (multiples of the current value) for kS and KD_RelA
KS_RANGE = (1e-3, 10.0)
KD_RELA_RANGE = (1e-3, 1e5)


def ribosome_elongation_rate(charged_trna, uncharged_trna, ribosome_conc, f,
		krta, krtf, max_elong_rate):
	# type: (np.ndarray, np.ndarray, float, np.ndarray, float, float, float) -> float
	"""
	Rate of amino acid incorporation (uM/s) for a given tRNA state, with the
	same competitive binding of charged and uncharged tRNA at the A site as
	the simulation (calculate_trna_charging).
	"""
	numerator = 1 + np.sum(
		f * (krta / charged_trna + uncharged_trna / charged_trna * krta / krtf))
	return max_elong_rate * ribosome_conc / numerator


def ppgpp_rates(uncharged_trna, charged_trna, ribosome_conc, f, v_rib,
		rela_conc, spot_conc, ppgpp_conc, krta, krtf, max_elong_rate,
		KD_RelA, k_RelA, k_SpoT_syn, k_SpoT_deg, KI_SpoT):
	"""
	ppGpp synthesis and degradation rates (uM/s) with the rate law of
	:py:func:`models.ecoli.processes.polypeptide_elongation.ppgpp_metabolite_changes`
	(RelA activated by ribosomes bound to uncharged tRNA, SpoT synthesis, and
	SpoT degradation inhibited by uncharged tRNA).

	Returns:
		v_rela_syn: RelA synthesis rate attributed to each amino acid
		v_spot_syn: SpoT synthesis rate
		v_deg: SpoT degradation rate
		v_deg_inhibited: degradation prevented by each amino acid's uncharged tRNA
	"""
	numerator = 1 + charged_trna / krta + uncharged_trna / krtf
	saturated_charged = charged_trna / krta / numerator
	saturated_uncharged = uncharged_trna / krtf / numerator
	if v_rib == 0:
		ribosome_conc_a_site = f * ribosome_conc
	else:
		ribosome_conc_a_site = f * v_rib / (saturated_charged * max_elong_rate)
	ribosomes_bound_to_uncharged = ribosome_conc_a_site * saturated_uncharged

	# Rare cases with no tRNA give inf or nan: assume ribosomes bind uncharged
	# tRNA in proportion to f if any tRNA of that amino acid is present
	mask = ~np.isfinite(ribosomes_bound_to_uncharged)
	ribosomes_bound_to_uncharged[mask] = (ribosome_conc * f[mask]
		* np.array(uncharged_trna[mask] + charged_trna[mask] > 0))

	# Fraction of RelA active, with competition between amino acids
	competitive_inhibition = 1 + ribosomes_bound_to_uncharged / KD_RelA
	inhibition_product = np.prod(competitive_inhibition)
	with np.errstate(divide='ignore'):
		frac_rela = 1 / (KD_RelA / ribosomes_bound_to_uncharged
			* inhibition_product / competitive_inhibition + 1)

	v_rela_syn = k_RelA * rela_conc * frac_rela
	v_spot_syn = k_SpoT_syn * spot_conc
	max_deg = k_SpoT_deg * spot_conc * ppgpp_conc
	fractions = uncharged_trna / KI_SpoT
	v_deg = max_deg / (1 + fractions.sum())
	v_deg_inhibited = (max_deg - v_deg) * fractions / fractions.sum()

	return v_rela_syn, v_spot_syn, v_deg, v_deg_inhibited


def _bulk_dry_mws(sim_data):
	# type: (Any) -> np.ndarray
	"""Dry molecular weight (g/mol) of each bulk molecule (water excluded)."""
	molecule_mws = sim_data.internal_state.bulk_molecules.bulk_data[
		'mass'].asNumber(units.g / units.mol)
	water_index = sim_data.submass_name_to_index['water']
	return molecule_mws.sum(axis=1) - molecule_mws[:, water_index]


def _dry_mass_of_counts(sim_data, container):
	# type: (Any, Any) -> float
	"""Dry mass (fg) of a bulk container, as the simulation's mass listener sums it."""
	n_avogadro = sim_data.constants.n_avogadro.asNumber(1 / units.mol)
	grams = container.counts() @ _bulk_dry_mws(sim_data) / n_avogadro
	return float(grams * 1e15)


def _aa_in_media(sim_data, nutrients, aa_ids):
	# type: (Any, str, list) -> np.ndarray
	"""True for each amino acid present in the media above the import threshold."""
	media = sim_data.external_state.saved_media[nutrients]
	threshold = sim_data.external_state.import_constraint_threshold
	if units.hasUnit(threshold):
		threshold = threshold.asNumber(units.mmol / units.L)
	present = []
	for aa in aa_ids:
		conc = media.get(aa[:-3], 0)
		if units.hasUnit(conc):
			conc = conc.asNumber(units.mmol / units.L)
		present.append(conc > threshold)
	return np.array(present)


def condition_charging_state(sim_data, condition):
	# type: (Any, str) -> Dict[str, Any]
	"""
	Concentrations (uM), counts and rate parameters for one condition, from the
	expected bulk counts the simulation starts from.
	"""
	transcription = sim_data.process.transcription
	metabolism = sim_data.process.metabolism
	translation = sim_data.process.translation
	constants = sim_data.constants
	growth = sim_data.growth_rate_parameters
	nutrients = sim_data.conditions[condition]['nutrients']
	doubling_time = sim_data.condition_to_doubling_time[condition]

	container = create_expected_bulk_container(sim_data, condition=condition)

	# Volume from the dry mass of the counts present, as in the simulation
	dry_mass = _dry_mass_of_counts(sim_data, container) * units.fg
	volume = dry_mass / sim_data.mass.cell_dry_mass_fraction / constants.cell_density
	to_um = (1 / (constants.n_avogadro * volume)).asNumber(UM)

	trna_counts = (container.counts(transcription.uncharged_trna_names)
		+ container.counts(transcription.charged_trna_names))
	total_trna = transcription.aa_from_trna @ trna_counts * to_um
	synthetase = (transcription.aa_from_synthetase
		@ container.counts(transcription.synthetase_names) * to_um)
	rela = container.count(sim_data.molecule_ids.RelA) * to_um
	spot = container.count(sim_data.molecule_ids.SpoT) * to_um
	n_30s = container.count(sim_data.molecule_ids.s30_full_complex)
	n_50s = container.count(sim_data.molecule_ids.s50_full_complex)
	# Active ribosomes: subunit pairs times the active fraction
	ribosome = (min(n_30s, n_50s)
		* growth.get_fraction_active_ribosome(doubling_time) * to_um)

	aa_ids = sim_data.molecule_groups.amino_acids
	conc = metabolism.concentration_updates.concentrations_based_on_nutrients(
		media_id=nutrients)
	aa_conc = np.array([conc[aa].asNumber(UM) for aa in aa_ids])

	# Amino acids used by translation (dilution plus degradation) and released
	# by protein degradation, per amino acid. The expected container has
	# already formed complexes, so dilution counts every monomer, free or in
	# a complex. Only free monomers are degraded in the simulation, so
	# degradation and recycling use free monomers alone.
	monomers = translation.monomer_data
	free_monomer_counts = container.counts(monomers['id'])
	complexation = sim_data.process.complexation
	in_complexes = dict(zip(complexation.molecule_names,
		-complexation.stoich_matrix_monomers() @ container.counts(complexation.ids_complexes)))
	total_monomer_counts = free_monomer_counts + np.array([
		max(in_complexes.get(monomer_id, 0), 0) for monomer_id in monomers['id']])
	deg_rate = monomers['deg_rate'].asNumber(1 / units.s)
	aa_counts = monomers['aa_counts'].asNumber()
	dilution_rate = np.log(2) / doubling_time.asNumber(units.s)
	aa_use = (total_monomer_counts * dilution_rate + free_monomer_counts * deg_rate) @ aa_counts
	recycled = (free_monomer_counts * deg_rate) @ aa_counts * to_um

	# Same amino acids as the simulation's charging (selenocysteine excluded)
	charging_mask = np.array([aa[:-3] != 'L-SELENOCYSTEINE' for aa in aa_ids], dtype=bool)
	aa_use[~charging_mask] = 0
	f = aa_use / aa_use.sum()

	# Supply inputs, as the simulation builds them
	fwd_enzymes, rev_enzymes = metabolism.get_pathway_enzyme_counts_per_aa(
		container.counts(metabolism.aa_enzymes))
	importers = container.counts(metabolism.aa_importer_names)
	exporters = container.counts(metabolism.aa_exporter_names)
	aa_in_media = _aa_in_media(sim_data, nutrients, aa_ids)

	target_ppgpp = growth.get_ppGpp_conc(doubling_time).asNumber(UM)
	max_elong_rate = growth.get_ribosome_elongation_rate_by_ppgpp(
		target_ppgpp * UM,
		constants.ribosome_elongation_rate_basal.asNumber(units.aa / units.s),
		).asNumber(units.aa / units.s)

	synthetase_kcats = getattr(transcription, 'synthetase_kcats', None)
	charged_trna_kis = getattr(transcription, 'charged_trna_kis', None)
	params = {
		'kS': (synthetase_kcats.asNumber(1 / units.s)
			if synthetase_kcats is not None
			else constants.synthetase_charging_rate.asNumber(1 / units.s)),
		'KMaa': transcription.aa_kms.asNumber(UM),
		'KMtf': transcription.trna_kms.asNumber(UM),
		'KI_charged_trna': (charged_trna_kis.asNumber(UM)
			if charged_trna_kis is not None
			else np.full(len(aa_ids), np.inf)),
		'krta': constants.Kdissociation_charged_trna_ribosome.asNumber(UM),
		'krtf': constants.Kdissociation_uncharged_trna_ribosome.asNumber(UM),
		'max_elong_rate': max_elong_rate,
		'unit_conversion': metabolism.get_amino_acid_conc_conversion(UM),
		'KD_RelA': transcription.KD_RelA.asNumber(UM),
		'k_RelA': constants.k_RelA_ppGpp_synthesis.asNumber(1 / units.s),
		'k_SpoT_syn': constants.k_SpoT_ppGpp_synthesis.asNumber(1 / units.s),
		'k_SpoT_deg': constants.k_SpoT_ppGpp_degradation.asNumber(
			units.L / units.umol / units.s),
		'KI_SpoT': transcription.KI_SpoT.asNumber(UM),
		}

	return {
		'condition': condition,
		'nutrients': nutrients,
		'doubling_time_min': doubling_time.asNumber(units.min),
		'dilution_rate': dilution_rate,
		'target_ppgpp': target_ppgpp,
		'to_um': to_um,
		'dry_mass': dry_mass,
		'total_trna': total_trna,
		'synthetase': synthetase,
		'rela': rela,
		'spot': spot,
		'ribosome': ribosome,
		'aa_conc': aa_conc,
		'recycled': recycled,
		'f': f,
		'charging_mask': charging_mask,
		'fwd_enzymes': fwd_enzymes,
		'rev_enzymes': rev_enzymes,
		'importers': importers,
		'exporters': exporters,
		'aa_in_media': aa_in_media,
		'params': params,
		'aa_ids': [aa[:-3] for aa in aa_ids],
		}


def steady_state_charging(sim_data, state, include_recycled=True):
	# type: (Any, Dict[str, Any], bool) -> Dict[str, Any]
	"""
	Steady-state tRNA charging and free amino acid concentrations with the
	simulation's charging, elongation and supply rate laws, holding enzyme,
	tRNA and ribosome concentrations at their expected values.

	Integrates the same system as
	:py:func:`models.ecoli.processes.polypeptide_elongation.calculate_trna_charging`
	(with supply) in rounds until the state stops changing.
	"""
	# Imported here to keep the ParCa module importable without the process
	from models.ecoli.processes.polypeptide_elongation import dcdt_jit

	metabolism = sim_data.process.metabolism
	p = state['params']
	mask = state['charging_mask']
	n_aas = len(mask)
	n_masked = int(mask.sum())
	to_um = state['to_um']
	synthetase = state['synthetase'][mask]
	synthetase_rate = p['kS']
	if np.ndim(synthetase_rate):
		synthetase_rate = np.asarray(synthetase_rate, dtype=float)[mask]
	f = state['f'][mask]
	recycled = state['recycled'] if include_recycled else np.zeros(n_aas)
	unit_conversion = p['unit_conversion']
	ki_charged = np.asarray(p['KI_charged_trna'], dtype=float)

	def supply(aa_conc):
		aa_conc = unit_conversion * aa_conc
		synthesis = metabolism.amino_acid_synthesis(
			state['fwd_enzymes'], state['rev_enzymes'], aa_conc)[0]
		imported = metabolism.amino_acid_import(state['aa_in_media'],
			state['dry_mass'], aa_conc, state['importers'], True)
		exported = metabolism.amino_acid_export(state['exporters'], aa_conc, True)
		return to_um * synthesis, to_um * imported, to_um * exported

	def dcdt(t, c):
		v_charging, dtrna, _ = dcdt_jit(t, c, n_masked, n_aas, mask,
			synthetase_rate, synthetase, p['KMaa'], p['KMtf'], f,
			p['krta'], p['krtf'], p['max_elong_rate'], state['ribosome'],
			False, np.zeros(n_masked), 0.0, ki_charged)
		aa_conc = c[2 * n_masked:]
		synthesis, imported, exported = supply(aa_conc)
		daa = np.zeros(n_aas)
		daa[mask] = (synthesis + imported - exported + recycled)[mask] - v_charging
		return np.hstack((-dtrna, dtrna, daa))

	# Start from 10% uncharged at the target amino acid concentrations
	total = state['total_trna'][mask]
	c = np.hstack((0.1 * total, 0.9 * total, state['aa_conc'].copy()))
	converged = False
	for _ in range(MAX_STEADY_STATE_ROUNDS):
		sol = solve_ivp(dcdt, [0, STEADY_STATE_TIME], c, method='BDF')
		new_c = np.fmax(sol.y[:, -1], 0)
		change = np.abs(new_c - c) / np.fmax(np.abs(c), 1e-9)
		c = new_c
		if change.max() < STEADY_STATE_TOLERANCE:
			converged = True
			break

	uncharged = np.zeros(n_aas)
	charged = np.zeros(n_aas)
	uncharged[mask] = c[:n_masked]
	charged[mask] = c[n_masked:2 * n_masked]
	aa_conc = c[2 * n_masked:]
	total_all = uncharged + charged
	with np.errstate(divide='ignore', invalid='ignore'):
		fraction_charged = np.where(total_all > 0, charged / total_all, 1.0)
	v_rib = ribosome_elongation_rate(charged[mask], uncharged[mask],
		state['ribosome'], f, p['krta'], p['krtf'], p['max_elong_rate'])
	return {
		'uncharged': uncharged,
		'charged': charged,
		'fraction_charged': fraction_charged,
		'aa_conc': aa_conc,
		'v_rib': v_rib,
		'converged': converged,
		}


def spot_factor_for_target(state, charging):
	# type: (Dict[str, Any], Dict[str, Any]) -> Tuple[float, Dict[str, float]]
	"""
	SpoT synthesis factor that puts ppGpp synthesis and loss (degradation plus
	dilution) in balance at the condition's target ppGpp for this charging.
	Synthesis is linear in the factor, so the factor follows directly.

	Returns:
		factor, and a dict of the rate terms at the target (uM/s)
	"""
	p = state['params']
	mask = state['charging_mask']
	target = state['target_ppgpp']
	v_rela, v_spot_unit, v_deg, _ = ppgpp_rates(
		charging['uncharged'][mask], charging['charged'][mask],
		state['ribosome'], state['f'][mask], charging['v_rib'],
		state['rela'], state['spot'], target,
		p['krta'], p['krtf'], p['max_elong_rate'], p['KD_RelA'][mask],
		p['k_RelA'], p['k_SpoT_syn'], p['k_SpoT_deg'], p['KI_SpoT'][mask])
	v_rela = float(np.sum(v_rela))
	loss = v_deg + state['dilution_rate'] * target
	factor = (loss - v_rela) / v_spot_unit
	return factor, {
		'v_rela_syn': v_rela,
		'v_spot_syn_unscaled': v_spot_unit,
		'v_deg': v_deg,
		'v_dilution': state['dilution_rate'] * target,
		}


def _fit_supply_and_ppgpp_kinetics(sim_data, cell_specs):
	"""
	Amino acid supply constants and ppGpp kinetics spread from the current
	expected bulk counts (same calls and order as final_adjustments).
	"""
	metabolism = sim_data.process.metabolism
	basal_container = create_expected_bulk_container(sim_data)
	with_aa_container = create_expected_bulk_container(sim_data, condition='with_aa')
	metabolism.set_phenomological_supply_constants(sim_data)
	metabolism.set_mechanistic_supply_constants(sim_data, cell_specs,
		basal_container, with_aa_container)
	metabolism.set_mechanistic_export_constants(sim_data, cell_specs, basal_container)
	metabolism.set_mechanistic_uptake_constants(sim_data, cell_specs, with_aa_container)
	sim_data.process.transcription.set_ppgpp_kinetics_parameters(
		basal_container, sim_data.constants)


def set_synthetase_kcats_from_capacity(sim_data,
		conditions=SYNTHETASE_CAPACITY_CONDITIONS):
	# type: (Any, Tuple[str, ...]) -> Optional[Dict[str, Any]]
	"""
	Per-synthetase charging rate constants that give each synthetase at least
	the measured charging capacity (Dykeman 2025 kcat x copies per cell, copies
	interpolated to each condition's growth rate) on the model's own synthetase
	counts:

		k_eff = max over conditions of kcat_Dykeman x copies_Dykeman(mu)/copies_model(mu)

	The maximum is used because the model's synthetase counts rise less with
	growth rate than Dykeman's for some synthetases (e.g. PheRS, on an
	r-protein operon), so a rate set at one growth rate leaves them short of
	capacity at the other. Expression, and therefore synthetase counts, stay as
	fit from RNA-seq. Amino acids without data keep the shared rate constant.
	The headroom floor (MEASURED_HEADROOM_MULTIPLIER x each synthetase's
	measured margin, when the table is present) then raises any rate below
	that margin over the model's own demand. Sets
	sim_data.process.transcription.synthetase_kcats.

	Returns:
		report with the per amino acid values, or None if there is no data
	"""
	transcription = sim_data.process.transcription
	data = getattr(transcription, 'synthetase_kcat_data', None)
	if data is None or np.all(np.isnan(data['kcat'])):
		return None

	shared_rate = sim_data.constants.synthetase_charging_rate.asNumber(1 / units.s)
	per_condition = {}
	for condition in conditions:
		doubling_time = sim_data.condition_to_doubling_time[condition]
		growth_rate = np.log(2) / doubling_time.asNumber(units.h)
		container = create_expected_bulk_container(sim_data, condition=condition)
		model_copies = transcription.aa_from_synthetase @ container.counts(
			transcription.synthetase_names)
		# Copies scale close to exponentially with growth rate: interpolate log copies
		measured_copies = np.array([
			np.exp(np.interp(growth_rate, data['growth_rates_per_h'], np.log(row)))
			if np.all(np.isfinite(row)) else np.nan
			for row in data['copies']
			])
		with np.errstate(divide='ignore', invalid='ignore'):
			rate = data['kcat'] * measured_copies / model_copies
		per_condition[condition] = {
			'growth_rate_per_h': growth_rate,
			'rate': rate,
			'copies_ratio': measured_copies / np.where(model_copies > 0, model_copies, np.nan),
			}

	rates = np.vstack([v['rate'] for v in per_condition.values()])
	has_data = np.all(np.isfinite(rates), axis=0)
	kcats = np.full(rates.shape[1], shared_rate)
	kcats[has_data] = rates[:, has_data].max(axis=0)

	# Headroom floor on this model's own demand in each condition
	floors = {}
	measured_headroom = getattr(transcription, 'synthetase_headroom', None)
	use_measured = measured_headroom is not None
	if use_measured and UNRELIABLE_HEADROOM:
		aa_names = [aa[:-3] for aa in sim_data.molecule_groups.amino_acids]
		replaced = np.array([aa in UNRELIABLE_HEADROOM for aa in aa_names])
		measured_headroom = {
			condition: np.where(replaced,
				np.nanmedian(np.where(replaced, np.nan, values)), values)
			for condition, values in measured_headroom.items()}
	if use_measured:
		for condition in conditions:
			state = condition_charging_state(sim_data, condition)
			elongation = sim_data.process.translation.ribosomeElongationRateDict[
				state['nutrients']].asNumber(units.aa / units.s)
			demand = state['f'] * elongation * state['ribosome']
			if MEASURED_HEADROOM_AT_MODEL_POOLS:
				km_aa = state['params']['KMaa']
				pool_saturation = state['aa_conc'] / (km_aa + state['aa_conc'])
			else:
				pool_saturation = 1.0
			with np.errstate(divide='ignore', invalid='ignore'):
				floor = MEASURED_HEADROOM_MULTIPLIER * measured_headroom[condition] * demand / (
					state['synthetase'] * pool_saturation)
			floor[~np.isfinite(floor) | ~state['charging_mask']] = 0.0
			floors[condition] = floor
			kcats[has_data] = np.fmax(kcats[has_data], floor[has_data])
	transcription.synthetase_kcats = (1 / units.s) * kcats

	aa_ids = [aa[:-3] for aa in sim_data.molecule_groups.amino_acids]
	return {
		'conditions': list(conditions),
		'growth_rate_per_h': {c: v['growth_rate_per_h'] for c, v in per_condition.items()},
		'kcat_effective': dict(zip(aa_ids, kcats)),
		'kcat_by_condition': {c: dict(zip(aa_ids, v['rate'])) for c, v in per_condition.items()},
		'kcat_measured': dict(zip(aa_ids, data['kcat'])),
		'min_headroom': (f'{MEASURED_HEADROOM_MULTIPLIER} x measured, at '
			+ ('model pools' if MEASURED_HEADROOM_AT_MODEL_POOLS else 'saturated pools')
			if use_measured else None),
		'headroom_floor_rate': {c: dict(zip(aa_ids, v)) for c, v in floors.items()},
		'copies_ratio_measured_over_model': {
			c: dict(zip(aa_ids, v['copies_ratio'])) for c, v in per_condition.items()},
		}


def _with_params(state, **changes):
	# type: (Dict[str, Any], Any) -> Dict[str, Any]
	"""Copy of a condition state with some rate parameters replaced."""
	new_state = dict(state)
	new_state['params'] = {**state['params'], **changes}
	return new_state


def mean_charged_fraction(state, charging):
	# type: (Dict[str, Any], Dict[str, Any]) -> float
	"""Unweighted mean charged fraction over the amino acids in charging."""
	return float(charging['fraction_charged'][state['charging_mask']].mean())


def fit_synthetase_rate(sim_data, state, target_charged_fraction, include_recycled=True):
	# type: (Any, Dict[str, Any], float, bool) -> Tuple[float, str]
	"""
	Shared synthetase rate constant kS (1/s) at which this condition's
	steady-state mean charged fraction equals the target. Charging rises
	monotonically with kS.
	"""
	ks_current = state['params']['kS']

	def mismatch(log_scale):
		charging = steady_state_charging(sim_data,
			_with_params(state, kS=ks_current * 10**log_scale),
			include_recycled=include_recycled)
		return mean_charged_fraction(state, charging) - target_charged_fraction

	try:
		log_scale = brentq(mismatch, *np.log10(KS_RANGE), xtol=1e-3)
		return ks_current * 10**log_scale, 'ok'
	except ValueError:
		return ks_current, 'no_solution'


def scaled_rela_null_ppgpp(sim_data, state):
	"""
	relA-null ppGpp (uM) for a condition: DEFAULT_RELA_NULL_PPGPP_UM scaled
	by the condition's target ppGpp over the target at
	RELA_NULL_REFERENCE_DOUBLING_TIME_MIN. Equal to
	DEFAULT_RELA_NULL_PPGPP_UM at that doubling time.
	"""
	reference_target = sim_data.growth_rate_parameters.get_ppGpp_conc(
		RELA_NULL_REFERENCE_DOUBLING_TIME_MIN * units.min).asNumber(UM)
	return DEFAULT_RELA_NULL_PPGPP_UM * (state['target_ppgpp'] / reference_target)


def rela_null_spot_factor(state, charging, rela_null_ppgpp):
	# type: (Dict[str, Any], Dict[str, Any], float) -> Tuple[float, float]
	"""
	SpoT synthesis factor at which SpoT alone (no RelA) holds ppGpp at the
	relA-null level for this charging:
	k_SpoT_syn x factor x SpoT = k_SpoT_deg x SpoT x P0 / (1 + sum(uncharged / KI_SpoT))
	+ dilution x P0.

	Returns:
		factor, and the SpoT degradation inhibition term (1 + sum)
	"""
	p = state['params']
	mask = state['charging_mask']
	inhibition = 1 + float((charging['uncharged'][mask] / p['KI_SpoT'][mask]).sum())
	loss = (p['k_SpoT_deg'] * state['spot'] * rela_null_ppgpp / inhibition
		+ state['dilution_rate'] * rela_null_ppgpp)
	return loss / (p['k_SpoT_syn'] * state['spot']), inhibition


def ki_spot_rescale(sim_data, state, charging, include_recycled=True):
	# type: (Any, Dict[str, Any], Dict[str, Any], bool) -> Tuple[float, float, float]
	"""
	Scale for KI_SpoT (uncharged tRNA inhibition of SpoT ppGpp degradation)
	that keeps its effect at the anchor condition's basal state the same in
	relative terms as when it was set. KI_SpoT (20 uM, "empirically determined
	for stable ppGpp") was set with one shared synthetase rate constant and no
	charged tRNA inhibition, where about 1% of tRNA is uncharged. No measured
	value exists, so it is scaled by the ratio of total uncharged tRNA now to
	total uncharged tRNA under those original charging kinetics.

	Returns:
		scale, reference uncharged total (uM), current uncharged total (uM)
	"""
	mask = state['charging_mask']
	reference_state = _with_params(state,
		kS=sim_data.constants.synthetase_charging_rate.asNumber(1 / units.s),
		KI_charged_trna=np.full(len(mask), np.inf))
	reference = steady_state_charging(sim_data, reference_state,
		include_recycled=include_recycled)
	reference_total = float(reference['uncharged'][mask].sum())
	current_total = float(charging['uncharged'][mask].sum())
	return current_total / reference_total, reference_total, current_total


def fit_kd_rela_scale(state, charging, spot_factor):
	# type: (Dict[str, Any], Dict[str, Any], float) -> Tuple[float, str]
	"""
	Scale on KD_RelA (all amino acids) at which SpoT (with this factor) plus
	RelA holds the condition at its target ppGpp. A larger KD_RelA makes RelA
	less sensitive to ribosomes bound to uncharged tRNA, so the SpoT factor
	needed for the target rises with the scale.
	"""
	kd = state['params']['KD_RelA']

	def mismatch(log_scale):
		needed, _ = spot_factor_for_target(
			_with_params(state, KD_RelA=kd * 10**log_scale), charging)
		return needed - spot_factor

	# Search away from the current value only: at much smaller KD_RelA the
	# product of per amino acid competition terms overflows
	low, high = np.log10(KD_RELA_RANGE)
	at_current = mismatch(0.0)
	if at_current == 0:
		return 1.0, 'ok'
	bracket = (0.0, high) if at_current < 0 else (low, 0.0)
	try:
		log_scale = brentq(mismatch, *bracket, xtol=1e-4)
		return 10**log_scale, 'ok'
	except ValueError:
		return 1.0, 'no_solution'


def fit_charging_attenuation_and_ppgpp(sim_data, cell_specs,
		fit_conditions=DEFAULT_FIT_CONDITIONS,
		anchor_condition=ANCHOR_CONDITION,
		anchor_charged_fraction=DEFAULT_ANCHOR_CHARGED_FRACTION,
		rela_null_ppgpp=None,
		charged_trna_ki_over_km=DEFAULT_CHARGED_TRNA_KI_OVER_KM,
		include_recycled=True,
		verbose=True):
	"""
	Set synthetase rates from measured capacities, make charging, attenuation
	and amino acid supply consistent, then fit KD_RelA and the SpoT synthesis
	factors (see the module docstring). Replaces the supply and ppGpp kinetics
	calls at the end of final_adjustments.

	Args:
		sim_data: SimulationDataEcoli object, modified in place
			(transcription.synthetase_kcats, transcription.charged_trna_kis,
			constants.KD_RelA_ribosome, transcription.KD_RelA,
			constants.KI_SpoT_ppGpp_degradation, transcription.KI_SpoT,
			transcription.spot_synthesis_media_factors, attenuation and amino
			acid supply constants)
		cell_specs: cell specifications from the fitted conditions
		fit_conditions (tuple[str]): conditions to solve charging for. Each
			must have distinct media. Attenuation in other media uses the
			charging of the solved condition with the same media, or all tRNA
			charged if none.
		anchor_condition (str): condition (in fit_conditions) whose relA-null
			level and target set its SpoT factor and KD_RelA (and whose
			measured charging sets kS when anchor_charged_fraction is given)
		anchor_charged_fraction (float | None): measured mean charged fraction
			for the anchor condition. None keeps kS at its current value.
		rela_null_ppgpp (float | None): relA-null ppGpp (uM) in the anchor
			condition. None scales DEFAULT_RELA_NULL_PPGPP_UM with the anchor
			condition's target ppGpp (see scaled_rela_null_ppgpp).
		charged_trna_ki_over_km (float | None): synthetase inhibition constant
			for charged tRNA as a multiple of the uncharged tRNA KM (None: no
			inhibition)
		include_recycled (bool): add amino acids released by protein
			degradation to supply in the charging solve (matches the simulation,
			which adds them to the amino acid targets)
		verbose (bool): print a summary for each iteration and condition

	Returns:
		report stored as sim_data.process.transcription.charging_fit_report
	"""
	transcription = sim_data.process.transcription
	constants = sim_data.constants

	# Synthetase inhibition by charged tRNA
	if charged_trna_ki_over_km is None:
		transcription.charged_trna_kis = None
	else:
		transcription.charged_trna_kis = transcription.trna_kms * charged_trna_ki_over_km

	# Attenuation volume from counts, as in the simulation, before the first fit
	transcription.calculate_attenuation(sim_data, cell_specs,
		charged_fraction=None, volume_from_counts=True)
	_fit_supply_and_ppgpp_kinetics(sim_data, cell_specs)

	# Per-synthetase rate constants from measured capacities (after the first
	# supply fit: the headroom floor needs the supply enzyme lists)
	capacity_report = set_synthetase_kcats_from_capacity(sim_data)
	if verbose and capacity_report is not None:
		print('Synthetase rate constants from Dykeman 2025 capacities (larger of'
			f' {", ".join(capacity_report["conditions"])}) (1/s): '
			+ ', '.join(f'{aa} {k:.1f}' for aa, k in capacity_report['kcat_effective'].items()))

	ks = constants.synthetase_charging_rate.asNumber(1 / units.s)
	ks_initial = ks
	if anchor_charged_fraction is not None and capacity_report is not None:
		raise ValueError('Fitting one shared kS to measured charging would override'
			' the per-synthetase rate constants; use one or the other')
	charged_fraction = None  # type: Optional[Dict[str, np.ndarray]]
	history = []
	converged = False
	ks_status = 'not_fit'
	fallback = None  # type: Optional[Dict[str, Any]]
	for iteration in range(MAX_LOOP_ITERATIONS):
		states = {c: condition_charging_state(sim_data, c) for c in fit_conditions}
		ks_change = 0.0
		if anchor_charged_fraction is not None:
			new_ks, ks_status = fit_synthetase_rate(sim_data, states[anchor_condition],
				anchor_charged_fraction, include_recycled=include_recycled)
			ks_change = abs(new_ks / ks - 1)
			ks = new_ks
			constants.synthetase_charging_rate = ks / units.s
			for state in states.values():
				state['params']['kS'] = ks

		charging = {
			c: steady_state_charging(sim_data, states[c], include_recycled=include_recycled)
			for c in fit_conditions
			}
		new_fraction = {
			states[c]['nutrients']: charging[c]['fraction_charged']
			for c in fit_conditions
			}
		change = (np.inf if charged_fraction is None
			else max(np.abs(new_fraction[m] - charged_fraction[m]).max() for m in new_fraction))
		history.append({
			'iteration': iteration,
			'kS': ks,
			'kS_status': ks_status,
			'max_charged_fraction_change': change,
			**{f'{c}_mean_charged': mean_charged_fraction(states[c], charging[c])
				for c in fit_conditions},
			})
		if verbose:
			print(f'Charging/attenuation loop {iteration}: kS {ks:.2f} /s'
				f' ({ks_status}), max change in charged fraction {change:.2e}; '
				+ ', '.join(f'{c} mean charged {history[-1][f"{c}_mean_charged"]:.3f}'
					for c in fit_conditions))
		if change < CHARGED_FRACTION_TOLERANCE and ks_change < 1e-2:
			converged = True
			break

		# Damped update of the charging that attenuation is calibrated on
		if charged_fraction is None:
			applied = new_fraction
		else:
			applied = {
				m: charged_fraction[m]
					+ ATTENUATION_DAMPING * (new_fraction[m] - charged_fraction[m])
				for m in new_fraction
				}
		try:
			transcription.calculate_attenuation(sim_data, cell_specs,
				charged_fraction=applied, volume_from_counts=True)
			_fit_supply_and_ppgpp_kinetics(sim_data, cell_specs)
		except ValueError as error:
			# Keep the last pass whose supply fit worked: restore its attenuation
			# and supply constants, whose charging was solved above
			transcription.calculate_attenuation(sim_data, cell_specs,
				charged_fraction=charged_fraction, volume_from_counts=True)
			_fit_supply_and_ppgpp_kinetics(sim_data, cell_specs)
			fallback = {'iteration': iteration, 'error': str(error)}
			if verbose:
				print(f'WARNING: charging/attenuation loop did NOT converge. Supply fit'
					f' failed at pass {iteration + 1} ({error}); keeping pass'
					f' {iteration}\'s attenuation, supply constants and charging')
			break
		charged_fraction = applied
	else:
		print(f'WARNING: charging/attenuation loop did NOT converge in'
			f' {MAX_LOOP_ITERATIONS} passes (last max change in charged fraction'
			f' {change:.2e}); the ppGpp fit uses the last pass\'s charging, which'
			' was solved before the last attenuation update')

	# KI_SpoT rescaled to the new charging level when charged tRNA inhibition is on
	ki_spot_scale, ki_reference_uncharged, ki_current_uncharged = 1.0, np.nan, np.nan
	if charged_trna_ki_over_km is not None:
		ki_spot_scale, ki_reference_uncharged, ki_current_uncharged = ki_spot_rescale(
			sim_data, states[anchor_condition], charging[anchor_condition],
			include_recycled=include_recycled)
		constants.KI_SpoT_ppGpp_degradation = constants.KI_SpoT_ppGpp_degradation * ki_spot_scale
		transcription.KI_SpoT = transcription.KI_SpoT * ki_spot_scale
		for state in states.values():
			state['params']['KI_SpoT'] = state['params']['KI_SpoT'] * ki_spot_scale
		if verbose:
			print(f'KI_SpoT x{ki_spot_scale:.2f}: basal uncharged tRNA'
				f' {ki_current_uncharged:.1f} uM vs {ki_reference_uncharged:.1f} uM'
				' under the original charging kinetics')

	# ppGpp: anchor condition sets its SpoT factor (relA-null) and KD_RelA
	factors = dict(getattr(transcription, 'spot_synthesis_media_factors', {}))
	previous_factors = dict(factors)
	fit_media = {states[c]['nutrients'] for c in fit_conditions}
	anchor_state, anchor_charging = states[anchor_condition], charging[anchor_condition]
	if rela_null_ppgpp is None:
		rela_null_ppgpp = scaled_rela_null_ppgpp(sim_data, anchor_state)
	if verbose:
		print(f'relA-null ppGpp anchor: {rela_null_ppgpp:.2f} uM'
			f' ({anchor_condition}, {anchor_state["doubling_time_min"]:.1f} min,'
			f' target {anchor_state["target_ppgpp"]:.2f} uM)')
	anchor_factor, anchor_inhibition = rela_null_spot_factor(
		anchor_state, anchor_charging, rela_null_ppgpp)
	kd_scale, kd_status = fit_kd_rela_scale(anchor_state, anchor_charging, anchor_factor)
	if kd_status != 'ok':
		print(f'WARNING: no KD_RelA scale puts {anchor_condition} ppGpp at its target'
			f' with the relA-null SpoT factor {anchor_factor:.2f} (searched'
			f' x{KD_RELA_RANGE[0]:g} to x{KD_RELA_RANGE[1]:g}); KD_RelA left unscaled')
	constants.KD_RelA_ribosome = constants.KD_RelA_ribosome * kd_scale
	transcription.KD_RelA = transcription.KD_RelA * kd_scale
	for state in states.values():
		state['params']['KD_RelA'] = state['params']['KD_RelA'] * kd_scale

	anchor_media = anchor_state['nutrients']
	previous_anchor_factor = factors.get(anchor_media, 1.0)
	factors[anchor_media] = anchor_factor
	# Media without amino acids that shared the anchor's previous factor follow it
	for media, value in previous_factors.items():
		if value == previous_anchor_factor and media not in fit_media:
			factors[media] = anchor_factor

	condition_reports = {}
	for condition in fit_conditions:
		state, result = states[condition], charging[condition]
		status = 'ok'
		if condition != anchor_condition:
			factor, _ = spot_factor_for_target(state, result)
			if factor <= 0:
				# RelA alone already exceeds the target at this charging
				status = 'target_below_rela_alone'
				factor = 1e-3
			factors[state['nutrients']] = factor
		factor = factors[state['nutrients']]
		needed, terms = spot_factor_for_target(state, result)
		mask = state['charging_mask']
		uncharged_fraction = 1 - result['fraction_charged']
		inhibition = 1 + float(
			(result['uncharged'][mask] / state['params']['KI_SpoT'][mask]).sum())
		condition_reports[condition] = {
			'nutrients': state['nutrients'],
			'doubling_time_min': state['doubling_time_min'],
			'target_ppgpp': state['target_ppgpp'],
			'spot_factor': factor,
			'previous_spot_factor': previous_factors.get(state['nutrients'], 1.0),
			'status': status,
			'steady_state_converged': result['converged'],
			'mean_charged': mean_charged_fraction(state, result),
			'rela_to_spot_synthesis': terms['v_rela_syn'] / (terms['v_spot_syn_unscaled'] * factor),
			'spot_inhibition': inhibition,
			'rate_terms': terms,
			'uncharged_fraction': dict(zip(state['aa_ids'], uncharged_fraction)),
			'aa_conc_over_target': dict(zip(state['aa_ids'], np.where(
				state['aa_conc'] > 0, result['aa_conc'] / state['aa_conc'], np.nan))),
			'v_rib_per_ribosome': result['v_rib'] / state['ribosome'],
			}
		if verbose:
			report = condition_reports[condition]
			top = np.argsort(-uncharged_fraction * mask)[:4]
			print(f'ppGpp fit {condition} ({state["nutrients"]},'
				f' {state["doubling_time_min"]:.0f} min): target'
				f' {state["target_ppgpp"]:.1f} uM, mean charged'
				f' {report["mean_charged"]:.3f}, SpoT factor'
				f' {report["previous_spot_factor"]:.2f} -> {factor:.2f} ({status}),'
				f' RelA/SpoT synthesis {report["rela_to_spot_synthesis"]:.2f},'
				f' SpoT inhibition {inhibition:.2f}, elongation'
				f' {report["v_rib_per_ribosome"]:.1f} aa/s, most uncharged: '
				+ ', '.join(f'{state["aa_ids"][k]} {uncharged_fraction[k]:.2f}' for k in top))
	transcription.spot_synthesis_media_factors = factors

	report = {
		'converged': converged,
		'fallback': fallback,
		'history': history,
		'kS_initial': ks_initial,
		'kS': ks,
		'kS_status': ks_status,
		'anchor_condition': anchor_condition,
		'anchor_charged_fraction': anchor_charged_fraction,
		'synthetase_capacity': capacity_report,
		'charged_trna_ki_over_km': charged_trna_ki_over_km,
		'rela_null_ppgpp': rela_null_ppgpp,
		'anchor_spot_inhibition': anchor_inhibition,
		'KI_SpoT_scale': ki_spot_scale,
		'KI_SpoT_reference_uncharged_uM': ki_reference_uncharged,
		'KI_SpoT_current_uncharged_uM': ki_current_uncharged,
		'KD_RelA_scale': kd_scale,
		'KD_RelA_status': kd_status,
		'KD_RelA_ribosome': constants.KD_RelA_ribosome.asNumber(UM),
		'conditions': condition_reports,
		'include_recycled': include_recycled,
		}
	if verbose:
		print(f'Charging fit: kS {ks_initial:.1f} -> {ks:.2f} /s ({ks_status}),'
			f' KD_RelA x{kd_scale:.3g} -> {report["KD_RelA_ribosome"]:.3g} uM'
			f' ({kd_status}), anchor relA-null SpoT inhibition'
			f' {anchor_inhibition:.2f}, converged {converged}')
	transcription.charging_fit_report = report
	return report
