"""
Total protein count comparison plotly scatter plots

This analysis produces two HTML plots for analyzing differences in monomer
counts between two simulations:
  1. ``..._categorized_sim1_<exp1>_sim2_<exp2>.html``: every monomer colored by
     its functional category (ribosome/RNAP/replisome/transcription factor (TF)
     /complexation/equilibrium/two-component-system (TCS) subunit).
  2. ``..._sim1_<exp1>_sim2_<exp2>.html``: every monomer is plotted in one
     background color except a user-specified highlight list
     (PLOT_PROTEINS_OF_INTEREST), where each highlighted protein gets its own
     color and its own legend entry.

Each plot uses the total monomer count listener to compute the average monomer
counts over all timesteps, then computes log10(avg. count + 1) values for what's
actually plotted. The +1 is added to each raw count to avoid log10(0). Sim 1
(x-axis) is the "reference" sim dir (the second sim dir given in the plot
run command), and Sim 2 (y-axis) is the "input" sim dir.

Each sim's monomer counts are read from their own sim_data/listener in case
the two sims being compared have different numbers of proteins or different
monomer id orderings (i.e. one sim had new proteins added to the model), so the
MonomerCounts listener array (ordered by that sim's own
``sim_data.process.translation.monomer_data["id"]``) is never indexed positionally
by the other sim's id list. Only monomer IDs present in both sims are plotted
(see the diagnostics printed by do_plot() for whether the two reconstructions
differ). A protein that exists in only one sim is dropped from the plot to
avoid assigning an artificial zero for the other sim. do_plot() prints how
many were dropped from each sim. The first subtitle line of each plot,
"<plotted>/<total> proteins plotted", shows the same thing: <total> counts every
protein that appears in either sim, so <plotted> < <total> means some proteins
were dropped, without having to check the printed output.

Note about variants:
This plot is meant for comparing sim dirs that each have a single variant
(wildtype). It reads each sim's root kb/simData.cPickle and averages over
the cells of every variant directory in the sim dir, so for a sim dir with
several variants, cells from different variants would be averaged together,
and the hover half-lives and categories would come from the unmodified sim_data
rather than each variant's simData_Modified.cPickle. Consider making this into
a multivariant plot if one is interested in running this between variants in
the same sim.

Note about hover data construction and category labels:
Category labels, complex membership, and reaction participation info in the hover
lists are built from Sim 1's sim_data only (see get_protein_categories).
Everything is matched by monomer id, never by array position, so this can never
shift points or numbers between proteins if the total number each sim has differs.
The only possible negative side effect of this could be the description labels,
where if Sim 2's reconstruction puts a shared protein in different complexes
or reactions (e.g. a complex was added to the flat files), that protein's
"Category" and hover lists show Sim 1's version only.

Note about replisomes and mechanistic_replisome:
The active_replisome count from UniqueMoleculeCounts is the number of
replisomes (one per replication fork) actually doing replication. Replication
always initiates when cell mass per origin reaches a critical mass (a fixed rule
that is not currently regulated by DnaA in the model). With mechanistic_replisome
set to False (the default), the simulation never checks whether enough replisome
subunits actually exist (it assumes they are not limiting) and never removes any
from the free pool (those replisomes carry no protein mass either). Only when
mechanistic_replisome is set to True does each initiation remove 3 of each
replisome trimer subunit and 1 of each monomer subunit per replisome from the
bulk pool and physically into the unique molecule space (they are then returned
at termination). See chromosome_replication.py for more information on all this.
The MonomerCounts listener only counts subunits as in active replisomes when
mechanistic_replisome is True. This plot reads the option from each sim's
metadata, states it in the plot titles, and prints what it means for each sim.
In the hover breakdown, it describes:
  - with the option on, "In replisomes" is the subunits sequestered in active
    replisomes (part of the marker total)
  - with it off, the same active_replisome * stoich number is instead shown as
    "Assumed counts available to form active replisomes (not actually
    sequestered with MR off)", which is not part of the marker total (those
    subunits were never removed from the free pool)
  - with it off, but in a sim run before the listener fix (when the listener
    added active replisome subunits to the totals regardless of the option),
    the same number is shown as "In replisomes (counted by the pre-fix
    MonomerCounts listener, but not actually sequestered with MR off)", which
    is part of the marker total. Such sims can't be told apart by their
    metadata, so this is detected from the data (see
    relabel_listener_counted_replisomes()), and do_plot() prints a note for
    each sim where it applies.

Hover text overview:
The marker value is the exact MonomerCounts total (which already folds in free
monomer + every complex + ribosomes/RNAPs + DNA-bound TFs, and replisomes when
mechanistic_replisome is on (see models/ecoli/listeners/monomer_counts.py)).
The hover data additionally reports the following per sim for each protein:
  - the total count with its average concentration in µM from the monomerCounts
    listener's ``monomerCounts`` column
  - a per-state breakdown of where that count lives (i.e. free monomer, in
    complexation, equilibrium, TCS complexes, in ribosomes, RNAPs, replisomes,
    DNA-bound as a TF), each with its average µM concentration. Every
    concentration is the mean over timepoints of each timepoint's count times
    that timepoint's countsToMolar (EnzymeKinetics listener), not the mean count
    times the mean countsToMolar, since counts and cell volume change together
    as the cell grows (see _step_means()). The buckets are built from
    separately-averaged free complex bulk counts, active unique-molecule counts,
    and DNA-bound TF counts, so their sum only approximates the exact marker total
    (small differences come from per-timestep integer rounding in the listener)
  - the protein's half life (given per sim since half lives can differ between sims)
  - the specific complex memberships, and the metabolic reactions the monomer
    catalyzes or participates in as a reactant/product (long lists capped at
    MAX_HOVER_IDS). When involvement is mediated by a complex rather than the
    monomer id itself, the entry reads "RXN-ID, via COMPLEX-ID". Each reaction is
    annotated with both sims' average flux, split by direction. Fluxes are
    read from the FBAResults listener's ``reactionFluxes`` column (ids in its
    ``reactionIDs`` attribute), and their units are mM/s (mmol/L/s). A
    direction the reaction does not have in either sim (e.g. "rev" for an
    irreversible, L2R reaction) is left out, and reactions with no FBA flux
    entry at all show no flux. Reactions are listed under their base id (the
    FBA id with any " (reverse)" tag stripped), so a reaction that runs in both
    directions gets one entry with a fwd and a rev line instead of two separate
    entries. Because of this, a reaction marked R2L in metabolic_reactions.tsv
    (written in EcoCyc in the opposite direction to how it runs, so FBA only
    has its "<id> (reverse)" version) is listed under its base id with only
    a rev line. This might look odd, but it's correct: R2L reactions have
    only rev flux.

Note about how metabolic involvement is defined:
  - a monomer "catalyzes" a reaction if its id is a metabolic catalyst (in
    ``sim_data.process.metabolism.catalyst_ids``) or if it is a subunit of a
    complex that is a catalyst (via the complexation/equilibrium complex->monomer
    maps built in get_protein_categories)
  - a monomer is a "reactant/product" if its id (or a complex it is a subunit of)
    appears as a participant key in some reaction's stoichiometry dict, i.e.
    ``sim_data.process.metabolism.reaction_stoich[rxn]`` (this is rare, and
    reactions where this is the case are often not even active in the model)

Stats box overview:
The stats listed in the stats box on the bottom right of plots are all
computed on the log10(count+1) values (see _add_parity_and_stats()):
  - Pearson r: linear correlation of the two sims' log counts (symmetric, so it
    does not depend on which sim is x vs y).
  - Pearson R² = r²: the fraction of variance explained by the best-fit line
    (not the y = x parity line drawn on the plot).
  - Spearman ρ: rank correlation, i.e. how well the two sims agree on the order
    of protein abundances regardless of whether the relationship is linear.
    Ranks are unchanged by the log10(count+1) transform, so unlike Pearson, it
    does not depend on that choice. Tied counts (e.g. the proteins at 0) get
    their average rank.

Note about how counts are averaged:
Monomer counts (and the hover breakdown counts, concentrations, and fluxes) are
averaged as the mean over all timepoints of all cells after the first
SKIP_INITIAL_GENERATIONS generations of each seed (thus giving equal weight per
timepoint). The average doubling time reported in the plot titles is taken over
that same set of cells, but with equal weight per cell, where
one doubling time = finalTime - initialTime per cell, since each cell has a
single doubling time rather than a per-timepoint value. Cells missing any of the
listener tables this plot reads are left out of every average
and of the cell count in the titles, and listed in a printed warning (see
Plot.get_cell_paths). This average could also be computed by getting the
average per gen and then averaging over all those to give equal weight to
each generation, but that is not implemented here so if one is interested in
making that change, they would have to do so manually.
"""

import json
import os
from typing import Tuple

import numpy as np
import plotly.express as px
import plotly.graph_objects as go
from scipy.stats import pearsonr, spearmanr

from models.ecoli.analysis import comparisonAnalysisPlot
from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from reconstruction.ecoli.dataclasses.process.metabolism import REVERSE_TAG
from reconstruction.ecoli.simulation_data import SimulationDataEcoli
from wholecell.analysis.analysis_tools import (
    read_stacked_bulk_molecules,
    read_stacked_columns,
)
from wholecell.io.tablereader import TableReader
from wholecell.utils import constants, units

# USER INPUTS:

# Proteins to highlight (each its own color) in the highlighted plot, as monomer
# ids without the compartment tag. Defaults to every replisome subunit monomer:
PLOT_PROTEINS_OF_INTEREST = [
    "EG10236-MONOMER",  # dnaB (helicase, hexamer)
    "EG10238-MONOMER",  # dnaE (Pol III core)
    "EG10239-MONOMER",  # dnaG (primase)
    "EG10242-MONOMER",  # dnaN (beta clamp, dimer)
    "EG10243-MONOMER",  # dnaQ (Pol III core)
    "EG11412-MONOMER",  # holA (clamp loader)
    "EG11500-MONOMER",  # holB (clamp loader)
    "EG11505-MONOMER",  # holE (Pol III core)
]

# Number of initial generations of each seed to skip when averaging:
SKIP_INITIAL_GENERATIONS = 2

# Cap for long id lists shown in hover text so a given protein's hover data
# doesn't blow up if it catalyzes several reactions (show first N, then "(+M more)"):
MAX_HOVER_IDS = 5

# END USER INPUTS

# USER INPUTS PART 2 (that aren't intended to be edited often but can be
# changed if desired for a different visual style):

# Figure sizing (see _apply_pdp_layout). The width is computed per figure so the
# plot area is square; the legend width is estimated from its longest entry:
FIGURE_HEIGHT = 760
FIGURE_MARGIN = dict(t=160, b=80, l=90, r=20)
LEGEND_PX_PER_CHAR = 7
LEGEND_PADDING_PX = 60

# Fixed colors for the single-category groups. Compound categories (proteins in
# more than one group) get their own color assigned dynamically below, one per
# unique combination actually present:
SINGLE_CATEGORY_COLORS = {
    'Monomer only': 'lightgray',
    'Ribosome subunit': 'darkorange',
    'RNAP subunit': 'magenta',
    'Replisome subunit': 'saddlebrown',
    'Transcription factor subunit': 'burlywood',
    'Complexation complex subunit': 'yellowgreen',
    'Equilibrium complex subunit': 'crimson',
    'Two component system subunit': 'teal',
}

# The two largest groups are drawn as a tiny faint background in the
# categorized plot (see build_category_styles). Everything else is emphasized
# with a larger marker and a thin black outline.
MONOMER_ONLY_LABEL = 'Monomer only'
BACKGROUND_COMPLEX_LABEL = 'Complexation complex subunit'

# Palette compound (multi-category) labels draw their colors from, assigned in
# sorted-label order so re-runs on the same data are deterministic:
COMPOUND_COLOR_PALETTE = (
    px.colors.qualitative.Plotly
    + px.colors.qualitative.Set2
    + px.colors.qualitative.Pastel
)

# Fixed colors for specific compound labels, used instead of their palette color:
COMPOUND_COLOR_OVERRIDES = {
    'Complexation complex subunit + Replisome subunit': 'cornflowerblue',
}

# Palette used to give each highlighted protein of interest its own color/legend
# entry in the highlighted (non-categorized) figure, cycled when there are more
# proteins than colors (the same bright palettes as the compound categories):
HIGHLIGHT_PALETTE = (
    px.colors.qualitative.Plotly
    + px.colors.qualitative.Set2
    + px.colors.qualitative.Pastel
)

# Categories drawn with a square marker in the categorized plot are the core
# molecular-machinery subunits. Everything else (monomer-only and the generic
# complexation/equilibrium/TCS complex subunits) are drawn as a circle. A
# protein whose compound label contains any of these counts as a square (see
# is_square_marker_label).
SQUARE_MARKER_LABELS = {
    'Ribosome subunit',
    'RNAP subunit',
    'Replisome subunit',
    'Transcription factor subunit',
}
# USER INPUTS PART 2 END

def _base_rxn_id(rxn_id):
    """
    Strip the ``REVERSE_TAG`` suffix so forward/reverse ids collapse to one
    base id (the hover then shows both directions' flux under that base).
    """
    if rxn_id.endswith(REVERSE_TAG):
        return rxn_id[: -len(REVERSE_TAG)]
    return rxn_id


def _dedup(seq):
    """
    De-duplicate preserving first-seen order.
    """
    seen = set()
    out = []
    for item in seq:
        if item not in seen:
            seen.add(item)
            out.append(item)
    return out


def get_protein_categories(sim_data, monomer_dict):
    """
    Categorize proteins into functional groups and build the reverse lookup maps
    used by the hover text. monomer_dict maps each plotted monomer id to its index
    in the plotted arrays, and monomers not in it are ignored.

    Returns a dict with:
      - the functional-group membership maps (``ribosome_monomers`` ...
        ``metabolic_enzyme_monomers`` / ``metabolic_reactant_monomers``), each
        ``monomer_id -> monomer_idx``;
      - the reverse hover maps (``monomer_to_*_complexes`` each
        ``monomer_id -> [complex id, ...]``, and ``monomer_to_catalyzed_rxns`` /
        ``monomer_to_reactant_rxns`` each ``monomer_id -> [(base_rxn_id, via), ...]``
        where ``via`` is the mediating complex id or None for direct involvement).

    NOTE: categories are computed from a single sim_data (do_plot passes Sim 1's).
    This assumes the reconstruction-level group membership is identical across the
    two compared sims (only the plotted monomer set is already intersected down to
    ids common to both).
    """

    # COMPLEXATION COMPLEXES
    complexation_molecule_ids = sim_data.process.complexation.molecule_names
    complexation_complex_ids = sim_data.process.complexation.ids_complexes
    complexation_smm = sim_data.process.complexation.stoich_matrix_monomers()
    complexation_monomers = {}
    complexes_to_monomers = {}

    for mol_id in complexation_molecule_ids:
        if mol_id in monomer_dict:
            complexation_monomers[mol_id] = monomer_dict[mol_id]
        elif mol_id in complexation_complex_ids:
            complex_idx = complexation_complex_ids.index(mol_id)
            monomer_idxs = np.where(complexation_smm[:, complex_idx] < 0)[0]
            monomers = []
            for idx in monomer_idxs:
                monomer_id = complexation_molecule_ids[idx]
                if monomer_id in monomer_dict:
                    complexation_monomers[monomer_id] = monomer_dict[monomer_id]
                    monomers.append(monomer_id)
            complexes_to_monomers[mol_id] = monomers

    # EQUILIBRIUM COMPLEXES
    equilibrium_molecule_ids = sim_data.process.equilibrium.molecule_names
    equilibrium_complex_ids = sim_data.process.equilibrium.ids_complexes
    eq_smm = sim_data.process.equilibrium.stoich_matrix_monomers()
    eq_complexes_to_monomers = {}
    equilibrium_monomers = {}

    for mol_id in equilibrium_molecule_ids:
        if mol_id in monomer_dict:
            equilibrium_monomers[mol_id] = monomer_dict[mol_id]
        elif mol_id in complexation_complex_ids:
            if mol_id in complexes_to_monomers:
                for monomer_id in complexes_to_monomers[mol_id]:
                    if monomer_id in monomer_dict:
                        equilibrium_monomers[monomer_id] = monomer_dict[monomer_id]
        elif mol_id in equilibrium_complex_ids:
            eq_complex_idx = equilibrium_complex_ids.index(mol_id)
            molecule_idxs = np.where(eq_smm[:, eq_complex_idx] < 0)[0]
            monomers = []
            for idx in molecule_idxs:
                molecule_name = equilibrium_molecule_ids[idx]
                if molecule_name in monomer_dict:
                    monomers.append(molecule_name)
                    equilibrium_monomers[molecule_name] = monomer_dict[molecule_name]
                elif molecule_name in complexation_complex_ids:
                    if molecule_name in complexes_to_monomers:
                        for monomer_id in complexes_to_monomers[molecule_name]:
                            if monomer_id in monomer_dict:
                                monomers.append(monomer_id)
                                equilibrium_monomers[monomer_id] = monomer_dict[monomer_id]
            eq_complexes_to_monomers[mol_id] = monomers

    # TWO COMPONENT SYSTEM COMPLEXES
    two_component_system_monomers = {}
    tcs_complexes_to_monomers = {}
    two_component_system_complex_ids = []

    try:
        two_component_system_molecule_ids = list(
            sim_data.process.two_component_system.modified_molecules)
        two_component_system_complex_ids = list(
            sim_data.process.two_component_system.complex_to_monomer.keys())
        tcs_complex_to_monomer_dict = sim_data.process.two_component_system.complex_to_monomer

        for mol_id in two_component_system_molecule_ids:
            subunit_ids = []
            if mol_id in monomer_dict:
                two_component_system_monomers[mol_id] = monomer_dict[mol_id]
            elif mol_id in complexation_complex_ids:
                subunit_ids = complexes_to_monomers[mol_id]
                tcs_complexes_to_monomers[mol_id] = subunit_ids
            elif mol_id in equilibrium_complex_ids:
                subunit_ids = eq_complexes_to_monomers[mol_id]
                tcs_complexes_to_monomers[mol_id] = subunit_ids
            elif mol_id in two_component_system_complex_ids:
                subunit_ids = list(tcs_complex_to_monomer_dict[mol_id].keys())
                tcs_complexes_to_monomers[mol_id] = subunit_ids
            for monomer_id in subunit_ids:
                if monomer_id in monomer_dict:
                    two_component_system_monomers[monomer_id] = monomer_dict[monomer_id]
    except AttributeError:
        # Print a warning if the two_component_system.modified_molecules
        # attribute is not found since the listener was added in 03/26
        print("Warning: two_component_system.modified_molecules not found "
              "(likely due to using an older sim version where sim_data version "
              "is out of date with current listeners used in the model).")

    # RIBOSOMES
    ribosome_50s_subunits = sim_data.process.complexation.get_monomers(
        sim_data.molecule_ids.s50_full_complex)
    ribosome_30s_subunits = sim_data.process.complexation.get_monomers(
        sim_data.molecule_ids.s30_full_complex)
    ribosome_subunits = (ribosome_50s_subunits["subunitIds"].tolist() +
                         ribosome_30s_subunits["subunitIds"].tolist())
    ribosome_monomers = {}
    for subunit in ribosome_subunits:
        if subunit in monomer_dict:
            ribosome_monomers[subunit] = monomer_dict[subunit]

    # RNAPs
    rnap_subunits = sim_data.process.complexation.get_monomers(
        sim_data.molecule_ids.full_RNAP)
    rnap_subunit_ids = rnap_subunits["subunitIds"].tolist()
    rnap_monomers = {}
    for subunit in rnap_subunit_ids:
        if subunit in monomer_dict:
            rnap_monomers[subunit] = monomer_dict[subunit]

    # REPLISOMES
    replisome_trimer_subunits = sim_data.molecule_groups.replisome_trimer_subunits
    replisome_monomer_subunits = sim_data.molecule_groups.replisome_monomer_subunits
    replisome_subunit_ids = replisome_trimer_subunits + replisome_monomer_subunits
    replisome_monomers = {}
    for subunit in replisome_subunit_ids:
        if subunit in monomer_dict:
            replisome_monomers[subunit] = monomer_dict[subunit]
        elif subunit in complexation_complex_ids:
            monomer_ids_list = complexes_to_monomers[subunit]
            for mol_id in monomer_ids_list:
                if mol_id in monomer_dict:
                    replisome_monomers[mol_id] = monomer_dict[mol_id]

    # TRANSCRIPTION FACTORS
    tfs = sim_data.process.transcription_regulation.tf_ids
    tf_subunit_ids = [tf_id + f'[{sim_data.getter.get_compartment(tf_id)[0]}]'
                      for tf_id in tfs]
    tf_monomers = {}
    for subunit in tf_subunit_ids:
        if subunit in monomer_dict:
            tf_monomers[subunit] = monomer_dict[subunit]
        elif subunit in complexation_complex_ids:
            monomer_ids_list = complexes_to_monomers[subunit]
            for mol_id in monomer_ids_list:
                if mol_id in monomer_dict:
                    tf_monomers[mol_id] = monomer_dict[mol_id]
        elif subunit in equilibrium_complex_ids:
            monomer_ids_list = eq_complexes_to_monomers[subunit]
            for mol_id in monomer_ids_list:
                if mol_id in monomer_dict:
                    tf_monomers[mol_id] = monomer_dict[mol_id]
        elif subunit in two_component_system_complex_ids:
            monomer_ids_list = tcs_complexes_to_monomers[subunit]
            for mol_id in monomer_ids_list:
                if mol_id in monomer_dict:
                    tf_monomers[mol_id] = monomer_dict[mol_id]

    # METABOLISM
    # A monomer is a "Metabolic enzyme" if its own id is a catalyst, or if it is
    # a subunit of a complex whose id is a catalyst. A monomer is a "Metabolic
    # reactant" if its id appears as a participant key in some reaction's
    # stoichiometry dict (i.e. the monomer is itself consumed/produced).
    metabolism = sim_data.process.metabolism
    catalyst_ids = set(metabolism.catalyst_ids)
    reaction_catalysts = metabolism.reaction_catalysts
    reaction_stoich = metabolism.reaction_stoich

    metabolic_enzyme_monomers = {}
    for cat_id in catalyst_ids:
        if cat_id in monomer_dict:
            metabolic_enzyme_monomers[cat_id] = monomer_dict[cat_id]
    for complex_id, monomers in {
        **complexes_to_monomers,
        **eq_complexes_to_monomers,
    }.items():
        if complex_id in catalyst_ids:
            for monomer_id in monomers:
                if monomer_id in monomer_dict:
                    metabolic_enzyme_monomers[monomer_id] = monomer_dict[monomer_id]

    metabolic_reactant_monomers = {}
    for stoich in reaction_stoich.values():
        for participant_id in stoich:
            if participant_id in monomer_dict:
                metabolic_reactant_monomers[participant_id] = monomer_dict[participant_id]

    # Reverse maps (monomer id -> [membership ids]) for hover detail. Each complex
    # map associates a monomer with every complex it is a subunit of.
    monomer_to_complexation_complexes = {}
    for complex_id, monomers in complexes_to_monomers.items():
        for monomer_id in monomers:
            monomer_to_complexation_complexes.setdefault(monomer_id, []).append(complex_id)
    monomer_to_equilibrium_complexes = {}
    for complex_id, monomers in eq_complexes_to_monomers.items():
        for monomer_id in monomers:
            monomer_to_equilibrium_complexes.setdefault(monomer_id, []).append(complex_id)
    monomer_to_tcs_complexes = {}
    for complex_id, monomers in tcs_complexes_to_monomers.items():
        for monomer_id in monomers:
            monomer_to_tcs_complexes.setdefault(monomer_id, []).append(complex_id)

    # monomer -> reactions it catalyzes/participates in, stored as
    # (base_rxn_id, via) tuples. `via` is the mediating complex id when the
    # involvement is through a complex the monomer is a subunit of, else None.
    # Forward/reverse ids collapse to the base id (the hover shows both
    # directions' flux), and the lists are de-duplicated preserving order.
    monomer_to_catalyzed_rxns = {}
    for rxn_id, catalysts in reaction_catalysts.items():
        base = _base_rxn_id(rxn_id)
        for cat_id in catalysts:
            if cat_id in monomer_dict:
                monomer_to_catalyzed_rxns.setdefault(cat_id, []).append((base, None))
            for monomer_id in (complexes_to_monomers.get(cat_id, [])
                               + eq_complexes_to_monomers.get(cat_id, [])):
                if monomer_id in monomer_dict:
                    monomer_to_catalyzed_rxns.setdefault(monomer_id, []).append(
                        (base, cat_id))
    monomer_to_catalyzed_rxns = {
        k: _dedup(v) for k, v in monomer_to_catalyzed_rxns.items()}

    monomer_to_reactant_rxns = {}
    for rxn_id, stoich in reaction_stoich.items():
        base = _base_rxn_id(rxn_id)
        for participant_id in stoich:
            if participant_id in monomer_dict:
                monomer_to_reactant_rxns.setdefault(participant_id, []).append(
                    (base, None))
            for monomer_id in (complexes_to_monomers.get(participant_id, [])
                               + eq_complexes_to_monomers.get(participant_id, [])):
                if monomer_id in monomer_dict:
                    monomer_to_reactant_rxns.setdefault(monomer_id, []).append(
                        (base, participant_id))
    monomer_to_reactant_rxns = {
        k: _dedup(v) for k, v in monomer_to_reactant_rxns.items()}

    return {
        'ribosome_monomers': ribosome_monomers,
        'rnap_monomers': rnap_monomers,
        'replisome_monomers': replisome_monomers,
        'tf_monomers': tf_monomers,
        'complexation_monomers': complexation_monomers,
        'equilibrium_monomers': equilibrium_monomers,
        'two_component_system_monomers': two_component_system_monomers,
        'metabolic_enzyme_monomers': metabolic_enzyme_monomers,
        'metabolic_reactant_monomers': metabolic_reactant_monomers,
        # Reverse maps for hover detail:
        'monomer_to_complexation_complexes': monomer_to_complexation_complexes,
        'monomer_to_equilibrium_complexes': monomer_to_equilibrium_complexes,
        'monomer_to_tcs_complexes': monomer_to_tcs_complexes,
        'monomer_to_catalyzed_rxns': monomer_to_catalyzed_rxns,
        'monomer_to_reactant_rxns': monomer_to_reactant_rxns,
    }

def classify_molecule_type(monomer_id, monomer_idx, protein_categories):
    """
    Classify a monomer into all functional-group categories it belongs to.

    Returns a single '<group> subunit' label when the monomer falls into exactly
    one group, or a ' + '-joined string when it falls into several (one legend
    entry/color per unique combination, see build_categorized_figure), or
    'Monomer only' when it is in none. Metabolic involvement is not included
    in the category labels to avoid cluttering the legend.
    """
    ribosome_monomers = protein_categories['ribosome_monomers']
    rnap_monomers = protein_categories['rnap_monomers']
    replisome_monomers = protein_categories['replisome_monomers']
    tf_monomers = protein_categories['tf_monomers']
    complexation_monomers = protein_categories['complexation_monomers']
    equilibrium_monomers = protein_categories['equilibrium_monomers']
    two_component_system_monomers = protein_categories['two_component_system_monomers']

    # The `[monomer_id] == monomer_idx` guard confirms the category dict's stored
    # index matches the monomer's plotted-array position:
    subunit_categories = []
    if monomer_id in ribosome_monomers and ribosome_monomers[monomer_id] == monomer_idx:
        subunit_categories.append('Ribosome')
    if monomer_id in rnap_monomers and rnap_monomers[monomer_id] == monomer_idx:
        subunit_categories.append('RNAP')
    if monomer_id in replisome_monomers and replisome_monomers[monomer_id] == monomer_idx:
        subunit_categories.append('Replisome')
    if monomer_id in tf_monomers and tf_monomers[monomer_id] == monomer_idx:
        subunit_categories.append('Transcription factor')
    if monomer_id in complexation_monomers and complexation_monomers[monomer_id] == monomer_idx:
        subunit_categories.append('Complexation complex')
    if monomer_id in equilibrium_monomers and equilibrium_monomers[monomer_id] == monomer_idx:
        subunit_categories.append('Equilibrium complex')
    if monomer_id in two_component_system_monomers and two_component_system_monomers[monomer_id] == monomer_idx:
        subunit_categories.append('Two component system')

    full_labels = [c + ' subunit' for c in subunit_categories]

    if not full_labels:
        return 'Monomer only'
    elif len(full_labels) == 1:
        return full_labels[0]
    else:
        return ' + '.join(sorted(full_labels))


def is_square_marker_label(mol_type: str) -> bool:
    """
    True if any category in mol_type (single or ' + '-joined compound label) is
    a core-machinery subunit (ribosome/RNAP/replisome/transcription factor).
    """
    return any(label in SQUARE_MARKER_LABELS for label in mol_type.split(' + '))

def build_category_styles(molecule_types):
    """
    Assign a color/size/opacity/symbol/line_width to every unique category label
    present in the data, following a visual hierarchy that keeps the big, boring
    groups out of the way and makes the interesting ones pop:

      - 'Monomer only' (by far the biggest group): light grey, tiny, no outline
      - 'Complexation complex subunit' (the next-biggest group): its fixed color
        but tiny and semi-transparent, also effectively background.
      - every other single category: its fixed SINGLE_CATEGORY_COLORS
        color, larger, with a thin black outline so it stands out.
      - compound labels (proteins in more than one category): their own color
        from COMPOUND_COLOR_PALETTE (or COMPOUND_COLOR_OVERRIDES), largest,
        thin black outline.

    Core-machinery subunits (single or compound) are drawn as squares, everything
    else is a circle. Marker size/opacity are set by tier here too, colors for single
    categories come from SINGLE_CATEGORY_COLORS.
    """
    present_labels = sorted(set(molecule_types))

    styles = {}
    compound_i = 0
    for label in present_labels:
        symbol = 'square' if is_square_marker_label(label) else 'circle'
        if label == MONOMER_ONLY_LABEL:
            styles[label] = dict(
                color=SINGLE_CATEGORY_COLORS[label], size=3, opacity=0.4,
                symbol='circle', line_width=0)
        elif label == BACKGROUND_COMPLEX_LABEL:
            styles[label] = dict(
                color=SINGLE_CATEGORY_COLORS[label], size=4, opacity=0.55,
                symbol='circle', line_width=0)
        elif ' + ' in label:
            # The palette index advances even for an overridden label, so an
            # override never shifts the other compound labels' colors:
            color = COMPOUND_COLOR_OVERRIDES.get(
                label,
                COMPOUND_COLOR_PALETTE[compound_i % len(COMPOUND_COLOR_PALETTE)])
            compound_i += 1
            styles[label] = dict(
                color=color, size=11, opacity=0.9, symbol=symbol, line_width=1)
        else:
            styles[label] = dict(
                color=SINGLE_CATEGORY_COLORS[label], size=9, opacity=0.9,
                symbol=symbol, line_width=1)
    return styles


def build_half_life_map(sim_data):
    """
    Map monomer_id -> half life (minutes) from a sim's own degradation rates,
    so hover text can show the value each individual sim actually used (since
    these can differ between sims).
    """
    monomer_data = sim_data.process.translation.monomer_data
    deg_rates = monomer_data["deg_rate"].asNumber(1 / units.s)
    half_life_min = {}
    for monomer_id, deg_rate in zip(monomer_data["id"], deg_rates):
        half_life_min[monomer_id] = (
            np.log(2) / deg_rate / 60.0 if deg_rate > 0 else np.inf
        )
    return half_life_min

def build_monomer_state_breakdown(
    sim_data, monomer_ids, bulk_by_id, active_means, bound_means,
    mechanistic_replisome
):
    """
    Per-monomer average count decomposed by state:
        {monomer_id: dict(free, complexation, equilibrium, tcs,
                          ribosome, rnap, replisome, replisome_assumed, bound)}
    computed from averaged bulk/active-unique/DNA-bound-TF data using each
    registry's stoich_matrix_monomers() (monomers-per-complex). The registry
    buckets use free complex bulk counts, so the buckets are disjoint from each
    other and from the free-monomer/unique-molecule/bound-TF buckets.

    These per-state buckets are rebuilt here because the MonomerCounts listener
    only records each monomer's total, not where that total lives. The total
    shown on the marker and in the hover ("Sim X total count") comes directly from
    MonomerCounts. The buckets below it add up close to that total but not
    exactly, because each one is averaged separately and the listener rounds to
    whole numbers each timestep.

    Active replisome subunits go in 'replisome' only when mechanistic_replisome
    is on (matching the listener). When it is off they go in 'replisome_assumed'
    instead, which is not part of the marker total (see the module docstring
    and relabel_listener_counted_replisomes()).
    """
    monomer_set = set(monomer_ids)
    def registry_data(molecule_names, complex_ids, smm):
        """
        Return (complex_ids, subunit_counts, name_to_row) for one complex
        registry. subunit_counts (monomers per complex) = max(0, -smm), with rows =
        molecule_names and cols = complex_ids; the same matrix feeds BOTH the
        free-bulk buckets and the DNA-bound-TF distribution.
        """
        subunit_counts = np.maximum(0, -np.asarray(smm))
        name_to_row = {m: i for i, m in enumerate(molecule_names)}
        return (list(complex_ids), subunit_counts, name_to_row)

    c_reg = registry_data(
        sim_data.process.complexation.molecule_names,
        sim_data.process.complexation.ids_complexes,
        sim_data.process.complexation.stoich_matrix_monomers(),
    )
    eq_reg = registry_data(
        sim_data.process.equilibrium.molecule_names,
        sim_data.process.equilibrium.ids_complexes,
        sim_data.process.equilibrium.stoich_matrix_monomers(),
    )
    # TCS stoich_matrix_monomers() rows are modified_molecules (not molecule_names)
    # as of the complex-listener merge, so the registry row labels must match:
    tcs_reg = None
    try:
        tcs_reg = registry_data(
            sim_data.process.two_component_system.modified_molecules,
            list(sim_data.process.two_component_system.complex_to_monomer.keys()),
            sim_data.process.two_component_system.stoich_matrix_monomers(),
        )
    except Exception as e:
        print(f"NOTE: TCS state breakdown unavailable: {e}")
        tcs_reg = None

    def free_contrib(reg):
        """
        Monomers sequestered in each FREE complex of this registry, i.e.
        subunit_counts @ (free bulk of each complex). Returns (per-row counts,
        name_to_row).
        """
        complex_ids, subunit_counts, name_to_row = reg
        bulk_vec = np.array(
            [float(bulk_by_id.get(c, 0.0)) for c in complex_ids], dtype=float
        )
        return subunit_counts @ bulk_vec, name_to_row

    c_complex_ids, c_subunit_counts, c_row = c_reg
    c_complex_col = {c: j for j, c in enumerate(c_complex_ids)}

    def unfold_complexation_complexes(contrib, name_to_row):
        """
        The equilibrium (and TCS) stoich_matrix_monomers() keep complexation
        complex subunits (e.g. the CRP dimer PC00004 in cAMP-CRP) as their own
        rows instead of breaking them down into monomers. This passes the
        counts on those rows through the complexation subunit_counts so they
        reach their monomers. Returns the added counts per complexation
        molecule_names row.
        """
        complex_counts = np.zeros(len(c_complex_ids))
        for name, row in name_to_row.items():
            col = c_complex_col.get(name)
            if col is not None:
                complex_counts[col] += contrib[row]
        return c_subunit_counts @ complex_counts

    def per_monomer(contrib, name_to_row, unfold):
        """
        {monomer_id: count} from one registry's per-row counts. With unfold
        set, counts on complexation complex rows are also passed down to
        their monomers (see unfold_complexation_complexes).
        """
        unfolded = (unfold_complexation_complexes(contrib, name_to_row)
                    if unfold else None)
        out = {}
        for m in monomer_ids:
            count = 0.0
            r = name_to_row.get(m)
            if r is not None:
                count += float(contrib[r])
            if unfolded is not None and m in c_row:
                count += float(unfolded[c_row[m]])
            if count:
                out[m] = count
        return out

    c_contrib, _ = free_contrib(c_reg)
    eq_contrib, eq_row = free_contrib(eq_reg)
    complexation_by_monomer = per_monomer(c_contrib, c_row, unfold=False)
    equilibrium_by_monomer = per_monomer(eq_contrib, eq_row, unfold=True)
    tcs_by_monomer = {}
    if tcs_reg is not None:
        tcs_contrib, tcs_row = free_contrib(tcs_reg)
        tcs_by_monomer = per_monomer(tcs_contrib, tcs_row, unfold=True)

    # Ribosome/RNAP/replisome per-monomer stoichiometry (monomers per one
    # assembled unique molecule):
    ribo50 = sim_data.process.complexation.get_monomers(
        sim_data.molecule_ids.s50_full_complex)
    ribo30 = sim_data.process.complexation.get_monomers(
        sim_data.molecule_ids.s30_full_complex)
    ribo_stoich = {}
    for ids, stoichs in (
        (ribo50["subunitIds"], ribo50["subunitStoich"]),
        (ribo30["subunitIds"], ribo30["subunitStoich"]),
    ):
        for s, v in zip(ids, stoichs):
            ribo_stoich[s] = ribo_stoich.get(s, 0.0) + float(v)
    rnap_info = sim_data.process.complexation.get_monomers(
        sim_data.molecule_ids.full_RNAP)
    rnap_stoich = {
        s: float(v) for s, v in zip(rnap_info["subunitIds"], rnap_info["subunitStoich"])
    }
    # Replisome subunits can be complexes themselves, so each subunit's
    # per-replisome stoich (3 per trimer subunit, 1 per monomer subunit) is
    # multiplied through to its monomers:
    repl_stoich = {}
    for subunits, per_replisome in (
        (sim_data.molecule_groups.replisome_trimer_subunits, 3.0),
        (sim_data.molecule_groups.replisome_monomer_subunits, 1.0),
    ):
        for subunit in subunits:
            if subunit in monomer_set:
                repl_stoich[subunit] = repl_stoich.get(subunit, 0.0) + per_replisome
                continue
            subunit_info = sim_data.process.complexation.get_monomers(subunit)
            for s, v in zip(subunit_info["subunitIds"], subunit_info["subunitStoich"]):
                repl_stoich[s] = repl_stoich.get(s, 0.0) + per_replisome * float(v)

    # DNA-bound TF counts, keyed by the TF's active form species id.
    # bound_means is keyed by the bare tf id.
    bound_by_tagged = {}
    for tf_id in sim_data.process.transcription_regulation.tf_ids:
        try:
            tagged = tf_id + f"[{sim_data.getter.get_compartment(tf_id)[0]}]"
        except (KeyError, IndexError):
            continue
        bc = float(bound_means.get(tf_id, 0.0))
        if bc:
            bound_by_tagged[tagged] = bc

    bound_contrib = {}
    # Direct case: the active form is itself a plotted monomer:
    for tagged, bc in bound_by_tagged.items():
        if tagged in monomer_set:
            bound_contrib[tagged] = bound_contrib.get(tagged, 0.0) + bc
    # Complex case: run the bound count through the registry it belongs to:
    for reg in (r for r in (c_reg, eq_reg, tcs_reg) if r is not None):
        complex_ids, subunit_counts, name_to_row = reg
        bvec = np.array(
            [bound_by_tagged.get(c, 0.0) for c in complex_ids], dtype=float
        )
        if not bvec.any():
            continue
        for m, count in per_monomer(
                subunit_counts @ bvec, name_to_row,
                unfold=reg is not c_reg).items():
            bound_contrib[m] = bound_contrib.get(m, 0.0) + count

    ribo_active = float(active_means.get("ribosome", 0.0))
    rnap_active = float(active_means.get("rnap", 0.0))
    repl_active = float(active_means.get("replisome", 0.0))

    replisome_key = "replisome" if mechanistic_replisome else "replisome_assumed"

    breakdown = {}
    for m in monomer_ids:
        breakdown[m] = dict(
            free=float(bulk_by_id.get(m, 0.0)),
            complexation=complexation_by_monomer.get(m, 0.0),
            equilibrium=equilibrium_by_monomer.get(m, 0.0),
            tcs=tcs_by_monomer.get(m, 0.0),
            ribosome=ribo_stoich.get(m, 0.0) * ribo_active,
            rnap=rnap_stoich.get(m, 0.0) * rnap_active,
            bound=float(bound_contrib.get(m, 0.0)),
        )
        breakdown[m][replisome_key] = repl_stoich.get(m, 0.0) * repl_active
    return breakdown


def relabel_listener_counted_replisomes(breakdown, monomer_means):
    """
    Detect whether a sim's MonomerCounts listener added active replisome
    subunits to the totals even with mechanistic_replisome off, as the listener
    did before the fix in models/ecoli/listeners/monomer_counts.py. Such sims
    can't be told apart by their metadata, so this checks the data: for each
    monomer with a 'replisome_assumed' bucket, the gap between its listener
    total (monomer_means) and the sum of its other buckets is compared to its
    assumed count. If the gaps are closer to the assumed counts than to 0, the
    assumed counts are moved to 'replisome_listener' (part of the marker
    total) in breakdown.

    Returns True if the counts were moved.
    """
    mismatch_if_counted = 0.0
    mismatch_if_not_counted = 0.0
    for monomer_id, bd in breakdown.items():
        assumed = bd.get('replisome_assumed', 0.0)
        if assumed <= 0:
            continue
        gap = monomer_means[monomer_id] - sum(
            v for k, v in bd.items() if k != 'replisome_assumed')
        mismatch_if_counted += abs(gap - assumed)
        mismatch_if_not_counted += abs(gap)
    if mismatch_if_counted >= mismatch_if_not_counted:
        return False
    rename_bucket(breakdown, 'replisome_assumed', 'replisome_listener')
    return True


def rename_bucket(breakdown, old_key, new_key):
    """
    Rename the old_key bucket to new_key in every monomer's breakdown.
    """
    for bd in breakdown.values():
        if old_key in bd:
            bd[new_key] = bd.pop(old_key)


def _step_means(arr, ctm_per_step):
    """
    Return (mean counts, mean concentrations in mM) over the timepoints (rows)
    of arr. Each timepoint's counts are multiplied by that timepoint's own
    countsToMolar (ctm_per_step, mM/count) before averaging, since counts and
    cell volume change together as the cell grows, so mean(count) *
    mean(countsToMolar) would be biased. The concentrations are None when
    ctm_per_step is None or does not have one entry per row of arr.
    """
    arr = np.asarray(arr, dtype=float)
    if arr.ndim == 1:
        arr = arr.reshape(-1, 1)
    means = arr.mean(axis=0)
    if ctm_per_step is None:
        return means, None
    if len(ctm_per_step) != arr.shape[0]:
        print(f"NOTE: countsToMolar has {len(ctm_per_step)} timepoints but "
              f"the counts have {arr.shape[0]}, so concentrations are left out.")
        return means, None
    return means, (arr * ctm_per_step[:, None]).mean(axis=0)


def _fmt_id_list(ids):
    """
    Format a list of ids for hover, one per line, capping at MAX_HOVER_IDS.
    """
    ids = list(ids)
    if not ids:
        return ""
    shown = "<br>".join(f"&nbsp;&nbsp;- {x}" for x in ids[:MAX_HOVER_IDS])
    if len(ids) > MAX_HOVER_IDS:
        shown += f"<br>&nbsp;&nbsp;- (+{len(ids) - MAX_HOVER_IDS} more)"
    return shown

def _fmt_flux_val(v):
    """
    Format one flux value (mM/s) for hover; 'n/a' when the reaction has no FBA
    flux entry in that sim/direction.
    """
    if v is None:
        return "n/a"
    return f"{v:.2g}"


def _fmt_rxn_list(entries, flux_map_1, flux_map_2):
    """
    Format reaction entries [(base_rxn_id, via_complex_or_None), ...] for hover.

    Each reaction gets a "<base_rxn>[, via <complex>]" line, followed by one
    indented line per direction, "fwd: S1 .. | S2 .. mM/s" and "rev: S1 .. |
    S2 .. mM/s", with the forward flux looked up at base_rxn and the reverse at
    base_rxn + REVERSE_TAG in each sim's flux map. A direction with no flux entry
    in either sim (e.g. the reverse of an irreversible reaction) is left out, so
    reactions with no flux entry at all get no flux lines. Capped at
    MAX_HOVER_IDS reactions. Fluxes are FBAResults.reactionFluxes, which are in
    mmol/L/s (shown as mM/s).
    """
    entries = list(entries)
    if not entries:
        return ""
    lines = []
    for base_rxn, via in entries[:MAX_HOVER_IDS]:
        head = f"{base_rxn}, via {via}" if via else base_rxn
        lines.append(f"&nbsp;&nbsp;- {head}")
        rev_id = base_rxn + REVERSE_TAG
        fwd1, rev1 = flux_map_1.get(base_rxn), flux_map_1.get(rev_id)
        fwd2, rev2 = flux_map_2.get(base_rxn), flux_map_2.get(rev_id)
        for direction, v1, v2 in (("fwd", fwd1, fwd2), ("rev", rev1, rev2)):
            if v1 is None and v2 is None:
                continue
            lines.append(
                f"&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;{direction}: "
                f"S1 {_fmt_flux_val(v1)} | S2 {_fmt_flux_val(v2)} mM/s")
    if len(entries) > MAX_HOVER_IDS:
        lines.append(f"&nbsp;&nbsp;- (+{len(entries) - MAX_HOVER_IDS} more)")
    return "<br>".join(lines)

def membership_hover_lines(monomer_id, protein_categories, flux_map_1, flux_map_2):
    """
    Build the per-protein hover lines naming the specific complex memberships
    and metabolic reactions (with per-sim, per-direction flux). Returns '' when
    there is nothing to add.
    """
    lines = []
    complexation = protein_categories['monomer_to_complexation_complexes'].get(monomer_id, [])
    equilibrium = protein_categories['monomer_to_equilibrium_complexes'].get(monomer_id, [])
    tcs = protein_categories['monomer_to_tcs_complexes'].get(monomer_id, [])
    catalyzed = protein_categories['monomer_to_catalyzed_rxns'].get(monomer_id, [])
    reactant = protein_categories['monomer_to_reactant_rxns'].get(monomer_id, [])

    if complexation:
        lines.append(f"<b>Complexation complex(es):</b><br>{_fmt_id_list(complexation)}")
    if equilibrium:
        lines.append(f"<b>Equilibrium complex(es):</b><br>{_fmt_id_list(equilibrium)}")
    if tcs:
        lines.append(f"<b>TCS complex(es):</b><br>{_fmt_id_list(tcs)}")
    if catalyzed:
        lines.append(
            f"<b>Catalyzes reaction(s):</b><br>{_fmt_rxn_list(catalyzed, flux_map_1, flux_map_2)}")
    if reactant:
        lines.append(
            f"<b>Reactant/product in reaction(s):</b><br>{_fmt_rxn_list(reactant, flux_map_1, flux_map_2)}")

    return ("<br>" + "<br>".join(lines)) if lines else ""

def _fmt_conc(conc_mM):
    """
    Format an average concentration (mM, see _step_means) for hover as
    ' (<x> µM)'. Returns '' when conc_mM is None (countsToMolar unavailable).
    """
    if conc_mM is None:
        return ""
    return f" ({conc_mM * 1e3:.3g} µM)"

def _fmt_sim_line(sim_label, exp_id, n_cells, dt, mechanistic_replisome):
    """
    Format one sim's title line, as '<sim_label>: <exp_id> | <n> cells |
    <x> min avg. doubling time | mechanistic replisome <on|off>'. The doubling
    time field is left out when dt is None.
    """
    fields = [f"{sim_label}: {exp_id}", f"{n_cells} cells"]
    if dt is not None:
        fields.append(f"{dt:.1f} min avg. doubling time")
    fields.append(
        f"mechanistic replisome {'on' if mechanistic_replisome else 'off'}")
    return " | ".join(fields)

def monomer_breakdown_hover(bd, bd_conc=None):
    """
    Format one monomer's average per-state breakdown for hover (nonzero only),
    each state annotated with its average concentration in µM from bd_conc
    (the same breakdown built from concentrations, see _step_means), when given.
    """
    order = [
        ("Free counts", "free"),
        ("In complexation complexes", "complexation"),
        ("In equilibrium complexes", "equilibrium"),
        ("In TCS complexes", "tcs"),
        ("In ribosomes", "ribosome"),
        ("In RNAPs", "rnap"),
        ("In replisomes", "replisome"),
        ("In replisomes (counted by the pre-fix MonomerCounts listener, but "
         "not actually sequestered with MR off)", "replisome_listener"),
        ("DNA-bound (as TF)", "bound"),
        ("Assumed counts available to form active replisomes (not actually "
         "sequestered with MR off)", "replisome_assumed"),
    ]
    parts = [
        f"&nbsp;&nbsp;{label}: {bd[key]:.1f}"
        f"{_fmt_conc(bd_conc.get(key) if bd_conc else None)}"
        for label, key in order
        if bd.get(key, 0.0) > 0.05
    ]
    if not parts:
        parts = ["&nbsp;&nbsp;(all states ~0)"]
    return "<br>".join(parts)

def _add_parity_and_stats(fig, sim1_log, sim2_log, r_value, pearson_r2, spearman_rho):
    """
    Add the shared y=x line and the Pearson/Spearman stats annotation.
    """
    max_val = max(sim1_log.max(), sim2_log.max())
    fig.add_trace(go.Scatter(
        x=[0, max_val],
        y=[0, max_val],
        mode='lines',
        line=dict(color='black', dash='dash', width=2),
        name='y = x',
        showlegend=True,
        hoverinfo='skip'
    ))

    stats_text = (
        f"Pearson r = {r_value:.3f}<br>"
        f"Pearson R² = {pearson_r2:.3f}<br>"
        f"Spearman ρ = {spearman_rho:.3f}"
    )
    # Placed relative to the axes (not the whole figure) so the box always sits
    # in the bottom-right corner of the plot area:
    fig.add_annotation(
        x=0.97,
        y=0.03,
        xref='x domain',
        yref='y domain',
        xanchor='right',
        yanchor='bottom',
        text=stats_text,
        showarrow=False,
        align='right',
        bgcolor='white',
        bordercolor='gray',
        borderwidth=1,
        borderpad=10,
        font=dict(size=11, family='monospace')
    )

def _apply_pdp_layout(fig, title, xaxis_title, yaxis_title):
    """
    Apply the shared layout: a white plot background with a full black axis
    border (showline + mirror on both axes), a left-aligned (multi-line) title,
    and a square plot area with the legend right beside it.
    """
    plot_side = FIGURE_HEIGHT - FIGURE_MARGIN['t'] - FIGURE_MARGIN['b']
    legend_chars = max(
        (len(trace.name) for trace in fig.data
         if trace.name and trace.showlegend is not False),
        default=0)
    legend_width = LEGEND_PX_PER_CHAR * legend_chars + LEGEND_PADDING_PX
    fig.update_layout(
        title=dict(text=title, x=0.01, xanchor='left', font=dict(size=14)),
        xaxis_title=xaxis_title,
        yaxis_title=yaxis_title,
        autosize=False,
        width=FIGURE_MARGIN['l'] + plot_side + legend_width,
        height=FIGURE_HEIGHT,
        margin=FIGURE_MARGIN,
        template='plotly_white',
        plot_bgcolor='white',
        hovermode='closest',
        showlegend=True,
    )

    # Same padded range on both axes (points at 0 are not cut off by the border):
    axis_max = max(
        float(np.max(np.concatenate([np.asarray(trace.x), np.asarray(trace.y)])))
        for trace in fig.data if trace.x is not None and len(trace.x))
    pad = 0.03 * axis_max
    axis_range = [-pad, axis_max + pad]
    fig.update_xaxes(showline=True, linewidth=1, linecolor='black', mirror=True,
                     range=axis_range, constrain='domain')
    fig.update_yaxes(showline=True, linewidth=1, linecolor='black', mirror=True,
                     range=axis_range, scaleanchor='x', scaleratio=1,
                     constrain='domain')

def build_categorized_figure(sim1_log, sim2_log, molecule_types, hover_texts,
                             r_value, pearson_r2, spearman_rho,
                             title, xaxis_title, yaxis_title):
    """
    Scatter with every protein colored by its functional category. Every unique
    category combination present in the data gets its own trace/legend entry.
    Traces are added largest-category-first, smallest-last. Plotly draws later
    traces on top and lists the legend in trace-addition order, so this puts the
    biggest groups in the background/top of the legend and the smallest on top
    so they do not get buried.
    """
    molecule_types_arr = np.array(molecule_types, dtype=object)
    styles = build_category_styles(molecule_types)

    counts_by_label = {
        label: int((molecule_types_arr == label).sum()) for label in styles
    }
    ordered_labels = sorted(styles, key=lambda label: (-counts_by_label[label], label))

    fig = go.Figure()
    for label in ordered_labels:
        mask = molecule_types_arr == label
        style = styles[label]
        fig.add_trace(go.Scatter(
            x=sim1_log[mask],
            y=sim2_log[mask],
            mode='markers',
            marker=dict(
                color=style['color'],
                size=style['size'],
                opacity=style['opacity'],
                symbol=style['symbol'],
                line=dict(width=style.get('line_width', 0), color='black'),
            ),
            name=f"{label} ({int(mask.sum())})",
            text=[hover_texts[i] for i in np.where(mask)[0]],
            hovertemplate='%{text}<extra></extra>',
            showlegend=True
        ))

    _add_parity_and_stats(fig, sim1_log, sim2_log, r_value, pearson_r2, spearman_rho)
    _apply_pdp_layout(fig, title, xaxis_title, yaxis_title)
    return fig

class Plot(comparisonAnalysisPlot.ComparisonAnalysisPlot):
    """
    Comparison plot of average total monomer counts between two sims.
    """
    def setup(self, inputDir: str) -> Tuple[AnalysisPaths, SimulationDataEcoli]:
        """
        Return the analysis paths and sim_data for one sim.
        """
        ap = AnalysisPaths(inputDir, variant_plot=True)
        sim_data = self.read_sim_data_file(inputDir)
        return ap, sim_data

    def get_cell_paths(self, ap, sim_label):
        """
        Cell paths for all generations after the first SKIP_INITIAL_GENERATIONS
        of each seed, the single cell set every per-sim read below averages over.

        Only cells where every column this plot reads can be fully read
        (i.e. not missing or truncated) are kept, so all reads average over the
        same cells and the cell count in the titles is the number of cells
        actually used. Any cells left out are listed in a printed warning. The
        reads below don't ignore exceptions, so they can never silently drop a
        cell from just one average.
        """
        all_cells = ap.get_cells(
            generation=np.arange(SKIP_INITIAL_GENERATIONS, ap.n_generation)
        )
        # Every listener table this plot reads, and the columns read from it:
        listener_columns = {
            "MonomerCounts": ["monomerCounts"],
            "BulkMolecules": ["counts"],
            "UniqueMoleculeCounts": ["uniqueMoleculeCounts"],
            "RnaSynthProb": ["nActualBound"],
            "EnzymeKinetics": ["countsToMolar"],
            "Main": ["time"],
            "FBAResults": ["reactionFluxes"],
        }
        readable_cells = []
        skipped_cells = []
        for cell_path in all_cells:
            try:
                for table, columns in listener_columns.items():
                    reader = TableReader(
                        os.path.join(cell_path, "simOut", table))
                    for column in columns:
                        reader.readColumn(column, squeeze=False)
            except Exception as e:
                skipped_cells.append((cell_path, e))
                continue
            readable_cells.append(cell_path)
        if skipped_cells:
            print(f"WARNING: {sim_label}: skipped {len(skipped_cells)} of "
                  f"{len(all_cells)} cells with missing or unreadable listener "
                  f"output (not included in the averages or the cell count):")
            for cell_path, e in skipped_cells:
                print(f"  {cell_path}: {e!r}")
        return np.array(readable_cells)

    def read_monomer_means(self, cell_paths, ctm_per_step):
        """
        Return (monomer_ids, mean_count_per_monomer, mean_conc_per_monomer).
        Counts are the mean over all timepoints of all cells, and concentrations
        (mM, or None) are from _step_means. monomer_ids are the MonomerCounts
        listener's own id labels (its emit order).
        """
        avg, conc = _step_means(read_stacked_columns(
            cell_paths, "MonomerCounts", "monomerCounts"
        ), ctm_per_step)
        reader = TableReader(
            os.path.join(cell_paths[0], "simOut", "MonomerCounts")
        )
        monomer_ids = list(reader.readAttribute("monomerIds"))
        return monomer_ids, avg, conc

    def read_bulk_by_id(self, cell_paths, ids, ctm_per_step):
        """
        Average bulk count (over all timepoints) for each requested id that is
        present in this sim's BulkMolecules listener. Returns ({id: mean count},
        {id: mean conc in mM} or None), see _step_means.
        """
        reader = TableReader(os.path.join(cell_paths[0], "simOut", "BulkMolecules"))
        bulk_set = set(reader.readAttribute("objectNames"))
        wanted = [i for i in ids if i in bulk_set]
        if not wanted:
            return {}, None
        (arr,) = read_stacked_bulk_molecules(
            cell_paths, wanted
        )
        means, conc = _step_means(arr, ctm_per_step)
        return (dict(zip(wanted, means.tolist())),
                dict(zip(wanted, conc.tolist())) if conc is not None else None)

    def read_active_means(self, cell_paths, ctm_per_step):
        """
        Average active-unique-molecule counts: ({'ribosome'|'rnap'|'replisome':
        mean count}, the same keyed mean concs in mM or None) from the
        UniqueMoleculeCounts listener, see _step_means.
        """
        try:
            umc = read_stacked_columns(
                cell_paths, "UniqueMoleculeCounts", "uniqueMoleculeCounts"
            )
            reader = TableReader(
                os.path.join(cell_paths[0], "simOut", "UniqueMoleculeCounts"))
            ids = list(reader.readAttribute("uniqueMoleculeIds"))
        except Exception as e:
            print(f"NOTE: could not read UniqueMoleculeCounts: {e!r}")
            return {}, None
        if umc.size == 0:
            return {}, None
        mean, conc = _step_means(umc, ctm_per_step)
        out = {}
        out_conc = {} if conc is not None else None
        for key, name in (
            ("active_ribosome", "ribosome"),
            ("active_RNAP", "rnap"),
            ("active_replisome", "replisome"),
        ):
            if key in ids:
                out[name] = float(mean[ids.index(key)])
                if conc is not None:
                    out_conc[name] = float(conc[ids.index(key)])
        return out, out_conc

    def read_bound_means(self, cell_paths, sim_data, ctm_per_step):
        """
        Average DNA-bound TF counts ({bare_tf_id: mean count}, {bare_tf_id: mean
        conc in mM} or None) from the RnaSynthProb listener's nActualBound column
        (emitted in that sim's tf_ids order), see _step_means.
        """
        tf_ids = list(sim_data.process.transcription_regulation.tf_ids)
        try:
            nab = read_stacked_columns(
                cell_paths, "RnaSynthProb", "nActualBound"
            )
        except Exception as e:
            print(f"NOTE: could not read nActualBound: {e!r}")
            return {}, None
        if nab.size == 0:
            return {}, None
        mean, conc = _step_means(nab, ctm_per_step)
        n = min(len(tf_ids), mean.shape[0])
        return (dict(zip(tf_ids[:n], mean[:n].tolist())),
                dict(zip(tf_ids[:n], conc[:n].tolist())) if conc is not None else None)

    def read_counts_to_molar_per_step(self, cell_paths):
        """
        countsToMolar (mM/count) at every timepoint from the EnzymeKinetics
        listener (one entry per row of the other stacked listener reads), or
        None if unavailable.
        """
        try:
            arr = read_stacked_columns(
                cell_paths, "EnzymeKinetics", "countsToMolar"
            )
        except Exception as e:
            print(f"NOTE: could not read countsToMolar: {e!r}")
            return None
        return np.asarray(arr, dtype=float).reshape(-1) if arr.size else None

    def read_doubling_time(self, cell_paths):
        """
        Average doubling time (minutes) over the same cells the counts are
        averaged over, but giving equal-weight to each cell in the average.
        """
        try:
            per_cell = read_stacked_columns(
                cell_paths, "Main", "time",
                fun=lambda x: np.array([[x[-1, 0] - x[0, 0]]])
            )
        except Exception as e:
            print(f"NOTE: could not read doubling time: {e!r}")
            return None
        return float(per_cell.mean()) / 60.0 if per_cell.size else None

    def read_flux_map(self, cell_paths):
        """
        Average FBA reaction flux (mM/s) {reaction_id: mean} from the FBAResults
        listener's reactionFluxes column (ids in reactionIDs). Reaction ids include
        both forward and ' (reverse)'-tagged entries.
        """
        try:
            flux = read_stacked_columns(
                cell_paths, "FBAResults", "reactionFluxes"
            )
            reader = TableReader(
                os.path.join(cell_paths[0], "simOut", "FBAResults"))
            rxn_ids = list(reader.readAttribute("reactionIDs"))
        except Exception as e:
            print(f"NOTE: could not read FBAResults fluxes: {e!r}")
            return {}
        if flux.size == 0:
            return {}
        mean = flux.mean(axis=0)
        return {rxn_ids[i]: float(mean[i]) for i in range(min(len(rxn_ids), mean.shape[0]))}

    def read_mechanistic_replisome(self, sim_dir):
        """
        Checks whether the sim ran with mechanistic_replisome, from its own metadata
        (see the module docstring). If the option is missing from the metadata
        it means the default was used, which is False
        (as of 10/21/2026).
        """
        metadata_path = os.path.join(
            sim_dir, constants.METADATA_DIR, constants.JSON_METADATA_FILE)
        with open(metadata_path) as f:
            sim_metadata = json.load(f)
        return bool(sim_metadata.get('mechanistic_replisome', False))

    def do_plot(self, reference_sim_dir, plotOutDir, plotOutFileName,
                input_sim_dir, unused, metadata):
        """
        Read both sims, match their monomers by id, build the shared hover
        text, and write the highlighted and categorized HTML plots to
        plotOutDir. Sim 1 = reference (x-axis), Sim 2 = input (y-axis).
        """
        ap1, sim_data1 = self.setup(reference_sim_dir)
        ap2, sim_data2 = self.setup(input_sim_dir)

        if (ap1.n_generation <= SKIP_INITIAL_GENERATIONS
                or ap2.n_generation <= SKIP_INITIAL_GENERATIONS):
            print("Skipping analysis -- not enough sims run.")
            return

        exp_id_1 = reference_sim_dir.split("out/")[-1].rstrip("/")
        exp_id_2 = input_sim_dir.split("out/")[-1].rstrip("/")
        print(f"Comparing {exp_id_1} (Sim 1; x-axis) vs {exp_id_2} (Sim 2; y-axis)")

        # Whether each sim actually sequestered replisome subunits (see the module
        # docstring), stated in the titles and the hover breakdown:
        mechanistic_replisome_1 = self.read_mechanistic_replisome(reference_sim_dir)
        mechanistic_replisome_2 = self.read_mechanistic_replisome(input_sim_dir)
        for sim_label, is_on in (("Sim 1", mechanistic_replisome_1),
                                 ("Sim 2", mechanistic_replisome_2)):
            if is_on:
                print(f"{sim_label}: Mechanistic replisome is on, so each "
                      f"initiation removes replisome subunits from the free pool "
                      f"and their counts include the subunits in active "
                      f"replisomes.")
            else:
                print(f"{sim_label}: Mechanistic replisome is off, so proteins "
                      f"are not actually sequestered to form replisomes and the "
                      f"sim assumes they are in excess.")

        cell_paths_1 = self.get_cell_paths(ap1, "Sim 1")
        cell_paths_2 = self.get_cell_paths(ap2, "Sim 2")
        n_cells_1 = len(cell_paths_1)
        n_cells_2 = len(cell_paths_2)
        if n_cells_1 == 0 or n_cells_2 == 0:
            print("Skipping analysis -- no readable cells after the skipped "
                  "generations in at least one sim.")
            return

        # Per-sim mean monomer counts keyed by that sim's own monomer id. The two
        # sims can differ in their monomer set/order, so counts are never indexed
        # positionally across sims.
        # countsToMolar (mM/count) at every timepoint per sim, to annotate counts
        # with their average µM concentration (see _step_means):
        ctm_1 = self.read_counts_to_molar_per_step(cell_paths_1)
        ctm_2 = self.read_counts_to_molar_per_step(cell_paths_2)

        monomer_ids_1, sim1_all, sim1_all_conc = self.read_monomer_means(
            cell_paths_1, ctm_1)
        monomer_ids_2, sim2_all, sim2_all_conc = self.read_monomer_means(
            cell_paths_2, ctm_2)
        print(f"Sim 1 has {n_cells_1} cells; Sim 2 has {n_cells_2} cells")
        print(f"Sim 1 total monomers: {len(monomer_ids_1)}")
        print(f"Sim 2 total monomers: {len(monomer_ids_2)}")

        means_1 = dict(zip(monomer_ids_1, (float(x) for x in sim1_all)))
        means_2 = dict(zip(monomer_ids_2, (float(x) for x in sim2_all)))
        conc_1 = (dict(zip(monomer_ids_1, sim1_all_conc.tolist()))
                  if sim1_all_conc is not None else {})
        conc_2 = (dict(zip(monomer_ids_2, sim2_all_conc.tolist()))
                  if sim2_all_conc is not None else {})

        # Half-life maps (per sim, since deg rates can differ between sims):
        half_life_map_1 = build_half_life_map(sim_data1)
        half_life_map_2 = build_half_life_map(sim_data2)

        # Plot only monomer ids present in both sims:
        plotted_ids = [mid for mid in means_1 if mid in means_2]
        not_plotted_1 = [mid for mid in means_1 if mid not in means_2]
        not_plotted_2 = [mid for mid in means_2 if mid not in means_1]
        n_total = len(plotted_ids) + len(not_plotted_1) + len(not_plotted_2)

        print(f"Plotted proteins (present in both sims): {len(plotted_ids)}/{n_total}")
        if not_plotted_1:
            print(
                f"NOTE: {len(not_plotted_1)} protein(s) not plotted (retrievable "
                f"only in Sim 1), e.g. {not_plotted_1[:5]}"
            )
        if not_plotted_2:
            print(
                f"NOTE: {len(not_plotted_2)} protein(s) not plotted (retrievable "
                f"only in Sim 2), e.g. {not_plotted_2[:5]}"
            )

        monomer_ids = plotted_ids
        sim1_avg = np.array([means_1[mid] for mid in monomer_ids])
        sim2_avg = np.array([means_2[mid] for mid in monomer_ids])
        monomer_dict = {mol: i for i, mol in enumerate(monomer_ids)}

        # Categorize proteins:
        protein_categories = get_protein_categories(sim_data1, monomer_dict)

        # Per-state hover breakdown inputs (per sim):
        complex_ids_1 = (list(sim_data1.process.complexation.ids_complexes)
                         + list(sim_data1.process.equilibrium.ids_complexes))
        complex_ids_2 = (list(sim_data2.process.complexation.ids_complexes)
                         + list(sim_data2.process.equilibrium.ids_complexes))
        try:
            complex_ids_1 += list(sim_data1.process.two_component_system.complex_to_monomer.keys())
        except AttributeError:
            pass
        try:
            complex_ids_2 += list(sim_data2.process.two_component_system.complex_to_monomer.keys())
        except AttributeError:
            pass
        # Each breakdown is linear in its bulk/active/bound inputs, so building it
        # again from the per-timepoint concentrations gives each bucket's
        # average concentration (None when countsToMolar is unavailable):
        breakdowns = []
        for (sim_data, cell_paths, complex_ids, ctm, mechanistic_replisome,
             means, sim_label) in (
                (sim_data1, cell_paths_1, complex_ids_1, ctm_1,
                 mechanistic_replisome_1, means_1, "Sim 1"),
                (sim_data2, cell_paths_2, complex_ids_2, ctm_2,
                 mechanistic_replisome_2, means_2, "Sim 2")):
            bulk_by_id, bulk_conc = self.read_bulk_by_id(
                cell_paths, monomer_ids + complex_ids, ctm)
            active_means, active_conc = self.read_active_means(cell_paths, ctm)
            bound_means, bound_conc = self.read_bound_means(
                cell_paths, sim_data, ctm)
            breakdown = build_monomer_state_breakdown(
                sim_data, monomer_ids, bulk_by_id, active_means, bound_means,
                mechanistic_replisome)
            breakdown_conc = None
            if bulk_conc is not None:
                breakdown_conc = build_monomer_state_breakdown(
                    sim_data, monomer_ids, bulk_conc, active_conc or {},
                    bound_conc or {}, mechanistic_replisome)
            if relabel_listener_counted_replisomes(breakdown, means):
                if breakdown_conc is not None:
                    rename_bucket(
                        breakdown_conc, 'replisome_assumed', 'replisome_listener')
                print(f"NOTE: {sim_label}'s MonomerCounts totals include "
                      f"active replisome subunits even though mechanistic "
                      f"replisome is off (the sim was run with the pre-fix "
                      f"listener), so the hover shows them as part of the "
                      f"total.")
            breakdowns.append((breakdown, breakdown_conc))
        (breakdown_1, breakdown_conc_1), (breakdown_2, breakdown_conc_2) = breakdowns

        # Average doubling time (minutes) per sim, for the plot titles:
        doubling_1 = self.read_doubling_time(cell_paths_1)
        doubling_2 = self.read_doubling_time(cell_paths_2)

        # Average FBA reaction flux (mM/s) per sim, for the hover reaction lines:
        flux_map_1 = self.read_flux_map(cell_paths_1)
        flux_map_2 = self.read_flux_map(cell_paths_2)

        # Gene symbol per monomer (monomer -> its own cistron -> gene symbol) for
        # hover text (reconstruction-level, so Sim 1's sim_data is fine):
        gene_data = sim_data1.process.replication.gene_data
        cistron_id_to_symbol = dict(zip(gene_data["cistron_id"], gene_data["symbol"]))
        monomer_id_to_cistron_id = dict(zip(
            sim_data1.process.translation.monomer_data["id"],
            sim_data1.process.translation.monomer_data["cistron_id"],
        ))
        monomer_id_to_gene_symbol = {
            mid: cistron_id_to_symbol.get(cid, "")
            for mid, cid in monomer_id_to_cistron_id.items()
        }

        molecule_types = [
            classify_molecule_type(monomer_id, i, protein_categories)
            for i, monomer_id in enumerate(monomer_ids)
        ]

        n_metabolic_enzyme = len(protein_categories['metabolic_enzyme_monomers'])
        n_metabolic_reactant = len(protein_categories['metabolic_reactant_monomers'])
        print(f"Metabolic enzyme monomers (hover-only): {n_metabolic_enzyme}")
        print(f"Metabolic reactant monomers (hover-only): {n_metabolic_reactant}")

        sim1_log = np.log10(sim1_avg + 1)
        sim2_log = np.log10(sim2_avg + 1)

        r_value = pearsonr(sim1_log, sim2_log)[0]
        pearson_r2 = r_value ** 2
        spearman_rho = spearmanr(sim1_log, sim2_log)[0]

        # Highlight list -> compartment-tagged monomer ids (matching plotted ids):
        converted_proteins = [
            protein_id + f'[{sim_data1.getter.get_compartment(protein_id)[0]}]'
            for protein_id in PLOT_PROTEINS_OF_INTEREST
        ]
        highlighted_set = set(converted_proteins)
        background_mask = np.array(
            [mid not in highlighted_set for mid in monomer_ids]
        )

        # Hover text (with per-state breakdown, concentrations, and the specific
        # complex/reaction memberships appended):
        hover_texts = []
        for i, monomer_id in enumerate(monomer_ids):
            half_life_1 = half_life_map_1.get(monomer_id)
            half_life_2 = half_life_map_2.get(monomer_id)
            half_life_1_str = (
                f"{half_life_1:.1f}" if half_life_1 is not None else "N/A"
            )
            half_life_2_str = (
                f"{half_life_2:.1f}" if half_life_2 is not None else "N/A"
            )
            gene_symbol = monomer_id_to_gene_symbol.get(monomer_id, "") or "N/A"
            hover_text = (
                f"<b>{monomer_id}</b><br>"
                f"<b>Gene:</b> {gene_symbol}<br>"
                f"<b>Category:</b> {molecule_types[i]}<br>"
                f"<b>Half-life (min):</b> S1: {half_life_1_str} | S2: {half_life_2_str}<br>"
                f"<b>Log10(counts + 1):</b> S1: {sim1_log[i]:.3f} | S2: {sim2_log[i]:.3f}<br>"
                f"<b>Sim 1 total count:</b> {sim1_avg[i]:.1f}{_fmt_conc(conc_1.get(monomer_id))}<br>"
                f"{monomer_breakdown_hover(breakdown_1[monomer_id], breakdown_conc_1 and breakdown_conc_1[monomer_id])}<br>"
                f"<b>Sim 2 total count:</b> {sim2_avg[i]:.1f}{_fmt_conc(conc_2.get(monomer_id))}<br>"
                f"{monomer_breakdown_hover(breakdown_2[monomer_id], breakdown_conc_2 and breakdown_conc_2[monomer_id])}"
                f"{membership_hover_lines(monomer_id, protein_categories, flux_map_1, flux_map_2)}"
            )
            hover_texts.append(hover_text)

        # Plot 1: highlighted proteins of interest, where each protein gets
        # its own color/legend entry.
        fig = go.Figure()
        fig.add_trace(go.Scatter(
            x=sim1_log[background_mask],
            y=sim2_log[background_mask],
            mode='markers',
            marker=dict(color='lightgray', size=4, opacity=0.4, line=dict(width=0)),
            name='All proteins',
            text=[hover_texts[i] for i in range(len(hover_texts)) if background_mask[i]],
            hovertemplate='%{text}<extra></extra>',
            showlegend=True
        ))

        highlighted_found = []
        highlighted_missing = []
        seen_highlights = set()
        for tagged_id in converted_proteins:
            if tagged_id in seen_highlights:
                continue
            seen_highlights.add(tagged_id)
            pos = monomer_dict.get(tagged_id)
            if pos is None:
                highlighted_missing.append(tagged_id)
                continue
            color = HIGHLIGHT_PALETTE[len(highlighted_found) % len(HIGHLIGHT_PALETTE)]
            highlighted_found.append(tagged_id)
            symbol = monomer_id_to_gene_symbol.get(tagged_id, "")
            legend_name = f"{tagged_id} ({symbol})" if symbol else tagged_id
            fig.add_trace(go.Scatter(
                x=[sim1_log[pos]],
                y=[sim2_log[pos]],
                mode='markers',
                marker=dict(
                    color=color, size=11, opacity=0.95,
                    line=dict(width=1, color='black')
                ),
                name=legend_name,
                text=[hover_texts[pos]],
                hovertemplate='%{text}<extra></extra>',
                showlegend=True
            ))
        _add_parity_and_stats(fig, sim1_log, sim2_log, r_value, pearson_r2, spearman_rho)

        # Title shared by both plots:
        plot_title = (
            f'Total Protein Count Comparison<br>'
            f'<sub>{len(monomer_ids)}/{n_total} proteins plotted | '
            f'skipped first {SKIP_INITIAL_GENERATIONS} gens from each seed when '
            f'calculating averages<br>'
            f'{_fmt_sim_line("Sim 1 (x)", exp_id_1, n_cells_1, doubling_1, mechanistic_replisome_1)}<br>'
            f'{_fmt_sim_line("Sim 2 (y)", exp_id_2, n_cells_2, doubling_2, mechanistic_replisome_2)}</sub>'
        )
        _apply_pdp_layout(
            fig, plot_title,
            'log10(Sim 1 Counts + 1)', 'log10(Sim 2 Counts + 1)'
        )
        output_filename = os.path.join(
            plotOutDir,
            f"{plotOutFileName}_sim1_{exp_id_1}_sim2_{exp_id_2}.html",
        )
        fig.write_html(output_filename)
        print(f"Saved plot to {output_filename}")
        print(
            f"Highlighted proteins: {len(highlighted_found)} of "
            f"{len(seen_highlights)} requested (each drawn as its own legend entry)"
        )
        if highlighted_missing:
            print(
                f"NOTE: {len(highlighted_missing)} highlighted id(s) not among the "
                f"plotted (shared) proteins, so not drawn: {highlighted_missing[:10]}"
                + (" ..." if len(highlighted_missing) > 10 else "")
            )

        # Plot 2: colored by protein category (no protein of interest highlighting)
        fig_categorized = build_categorized_figure(
            sim1_log, sim2_log, molecule_types, hover_texts,
            r_value, pearson_r2, spearman_rho,
            plot_title,
            'log10(Sim 1 Counts + 1)',
            'log10(Sim 2 Counts + 1)'
        )
        categorized_filename = os.path.join(
            plotOutDir,
            f"{plotOutFileName}_categorized_sim1_{exp_id_1}_sim2_{exp_id_2}.html",
        )
        fig_categorized.write_html(categorized_filename)
        print(f"Saved categorized plot to {categorized_filename}")

if __name__ == "__main__":
    Plot().cli()
