"""
Total protein (monomer) count comparison plotly scatter plots!

Produces two HTML plots:
  1. ``..._categorized_sim1_<exp1>_sim2_<exp2>.html``: every monomer colored by
     its functional category (ribosome / RNAP / replisome / transcription factor
     / complexation / equilibrium / two-component-system subunit), with the number
     of each in the legend. Metabolic involvement is shown only in the hover text
     (see below), NOT as a legend category, to keep the legend from being
     overloaded. Proteins that belong to more than one functional group get their
     own compound-category legend entry (e.g. "Ribosome subunit + Complexation
     complex subunit") rather than being folded into a catch-all.
  2. ``..._sim1_<exp1>_sim2_<exp2>.html``: every monomer is plotted in one
     background color (lightseagreen) except a user-specified highlight list
     (PLOT_PROTEINS_OF_INTEREST), where EACH highlighted protein gets its own
     color and its own legend entry (labeled "<gene symbol> (<monomer id>)"),
     cycling HIGHLIGHT_PALETTE -- the same per-item highlighting scheme as the
     flux comparison plots, so several proteins of interest can be told apart at
     a glance instead of all sharing one color. The legend is placed to the right
     since the entry count can be large.

Plots the average monomer counts from two simulations against each other.
Averages are log10-transformed (+1 pseudocount) and plotted on linear axes with a
square layout, a y = x parity line, and a Pearson r / Pearson R² / COD R² stats
box computed over all plotted proteins. A +1 is added to each raw count to avoid
log10(0).

Each sim's monomer counts are resolved from its own sim_data/listener in case
the two sims being compared have different numbers of proteins or different
monomer id orderings (e.g. one sim had new proteins added to the model), so the
MonomerCounts listener array (ordered by that sim's own
``sim_data.process.translation.monomer_data["id"]``) is never indexed positionally
by the other sim's id list. Only monomer ids present in both sims are plotted (see
the diagnostics printed by do_plot() for whether the two reconstructions differ).

Sim 1 (x-axis) is the reference sim dir; Sim 2 (y-axis) is the input sim dir.

HOVER TEXT:
The marker value is the exact MonomerCounts total (which already folds in free
monomer + every complex + ribosomes/RNAPs/replisomes + DNA-bound TFs -- see
models/ecoli/listeners/monomer_counts.py). The hover additionally reports, per sim:
  - the total count with its average concentration in µM;
  - a per-state breakdown of WHERE that count lives (free monomer / in
    complexation, equilibrium, TCS complexes / in ribosomes, RNAPs, replisomes /
    DNA-bound as a TF), each with its µM concentration. The buckets are built from
    separately-averaged free complex bulk counts, active unique-molecule counts,
    and DNA-bound TF counts, so their sum only APPROXIMATES the exact marker total
    (small differences come from per-timestep integer rounding in the listener);
  - the protein's half life (per sim -- deg rates can differ between sims);
  - the specific complex memberships, and the metabolic reactions the monomer
    catalyzes or participates in as a reactant/product (long lists capped at
    MAX_HOVER_IDS). When involvement is mediated by a complex rather than the
    monomer id itself, the entry reads "RXN-ID, via COMPLEX-ID". Each reaction is
    annotated with both sims' average flux, split by direction -- forward and
    reverse are kept as separate FBA reaction ids in the model, so both are shown:
    "(fwd S1:.. S2:.. | rev S1:.. S2:.. mM/s)". Fluxes are read from the FBAResults
    listener's ``reactionFluxes`` column (ids in its ``reactionIDs`` attribute);
    their unit is mM/s (mmol/L/s). Reactions with no FBA flux entry show no flux.

NOTE: here is how we define if/how a protein has metabolic involvement:
  - a monomer "catalyzes" a reaction if its id is a metabolic catalyst (in
    ``sim_data.process.metabolism.catalyst_ids``) or if it is a subunit of a
    complex that is a catalyst (via the complexation/equilibrium complex->monomer
    maps built in get_protein_categories);
  - a monomer is a "reactant/product" if its id (or a complex it is a subunit of)
    appears as a participant key in some
    ``sim_data.process.metabolism.reaction_stoich[rxn]`` (this is rare, and
    reactions where this is the case may not even be active in the model).


STATS BOX (all three computed on the log10(count+1) values, see
_add_parity_and_stats):
  - Pearson r: linear correlation of the two sims' log counts (symmetric, so it
    does not depend on which sim is x vs y).
  - Pearson R² = r**2: the fraction of variance explained by the best-fit line,
    which is NOT the y = x parity line drawn on the plot.
  - COD R² (coefficient of determination via sklearn r2_score, measured against
    the y = x line): how close the points actually are to y=x. It is computed
    as r2_score(y_true=sim2_log, y_pred=sim1_log), i.e. the fraction of Sim 2's
    variance the y = x line accounts for. Because the denominator is the y-axis
    sim's variance it is asymmetric, swapping x and y changes it, and it can
    go negative when the sims disagree worse than predicting the mean. A high
    Pearson R² with a low/negative COD R² means the sims are well correlated but
    systematically off the y=x line.

AVERAGING NOTE:
Monomer counts (and the hover breakdown counts, fluxes, and counts_to_molar) are
averaged as the mean over ALL timepoints of ALL cells after the first
SKIP_INITIAL_GENERATIONS generations of each seed (thus giving equal weight per
timepoint). The average doubling time reported in the plot titles is taken over
that same set of cells, but with equal weight per cell
(one doubling time = finalTime - initialTime per cell), since each cell has a
single doubling time rather than a per-timepoint value.

TF note: while there are technically many TFs in the model, the only TFs that get
a "Transcription factor subunit" category here are those actively modeled (the
list in reconstruction/ecoli/flat/condition/tf_condition.tsv).
"""

import os
from typing import Tuple

import numpy as np
import plotly.express as px
import plotly.graph_objects as go
from scipy.stats import pearsonr
from sklearn.metrics import r2_score

from models.ecoli.analysis import comparisonAnalysisPlot
from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from reconstruction.ecoli.dataclasses.process.metabolism import REVERSE_TAG
from reconstruction.ecoli.simulation_data import SimulationDataEcoli
from validation.ecoli.validation_data import ValidationDataEcoli
from wholecell.analysis.analysis_tools import (
    read_stacked_bulk_molecules,
    read_stacked_columns,
)
from wholecell.io.tablereader import TableReader
from wholecell.utils import units


# Cap for long id lists shown in hover text (show first N, then "(+M more)"):
MAX_HOVER_IDS = 5

# Number of initial generations of each seed to skip when averaging.
SKIP_INITIAL_GENERATIONS = 2

# Proteins to highlight (each its own color) in the highlighted plot:
PLOT_PROTEINS_OF_INTEREST = ["RPOD-MONOMER", "EG10239-MONOMER"]

# Flux unit string shown in hover text (FBAResults.reactionFluxes are mmol/L/s):
FLUX_UNIT = "mM/s"

# Fixed colors for the single-category groups. Compound categories (proteins in
# more than one group) get their own color assigned dynamically below, one per
# unique combination actually present -- no catch-all "Multiple categories".
# NOTE: metabolic involvement is intentionally NOT a category here (it would swamp
# the legend); it is surfaced only in the hover text. See the module docstring.
SINGLE_CATEGORY_SPECS = [
    ('Monomer only', 'lightgray', 3, 0.4),
    ('Ribosome subunit', 'darkorange', 8, 0.85),
    ('RNAP subunit', 'magenta', 8, 0.85),
    ('Replisome subunit', 'royalblue', 8, 0.85),
    ('Transcription factor subunit', 'green', 8, 0.85),
    ('Complexation complex subunit', 'yellowgreen', 8, 0.85),
    ('Equilibrium complex subunit', 'crimson', 8, 0.85),
    ('Two component system subunit', 'mediumpurple', 8, 0.85),
]

# The two largest, least interesting groups, drawn as a tiny faint background in
# the categorized plot (see build_category_styles). Everything else is emphasized
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

# Palette used to give each highlighted protein of interest its own color/legend
# entry in the highlighted (non-categorized) figure, cycled when there are more
# proteins than colors:
HIGHLIGHT_PALETTE = (
    px.colors.qualitative.Dark24 + px.colors.qualitative.Light24
)

# Categories drawn with a square marker in the categorized plot: the core
# molecular-machinery subunits. Everything else (monomer-only AND the generic
# complexation / equilibrium / TCS complex subunits) is drawn as a circle. A
# protein whose compound label contains ANY of these counts as a square (see
# is_square_marker_label).
SQUARE_MARKER_LABELS = {
    'Ribosome subunit',
    'RNAP subunit',
    'Replisome subunit',
    'Transcription factor subunit',
}


def _base_rxn_id(rxn_id):
    """Strip the ``REVERSE_TAG`` suffix so forward/reverse ids collapse to one
    base id (the hover then shows both directions' flux under that base)."""
    if rxn_id.endswith(REVERSE_TAG):
        return rxn_id[: -len(REVERSE_TAG)]
    return rxn_id


def _dedup(seq):
    """De-duplicate preserving first-seen order."""
    seen = set()
    out = []
    for item in seq:
        if item not in seen:
            seen.add(item)
            out.append(item)
    return out


def get_protein_categories(sim_data, monomer_ids, monomer_dict):
    """
    Categorize proteins into functional groups and build the reverse lookup maps
    used by the hover text.

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

    def get_monomer_indexes(keys):
        indices = []
        for x in keys:
            if x in monomer_dict:
                indices.append(monomer_dict[x])
        return np.array(indices)

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
            if mol_id in monomer_ids:
                two_component_system_monomers[mol_id] = get_monomer_indexes([mol_id])[0]
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
                if monomer_id in monomer_ids:
                    two_component_system_monomers[monomer_id] = (
                        get_monomer_indexes([monomer_id]))[0]
    except AttributeError:
        print("Warning: two_component_system.modified_molecules not found "
              "(likely due to using an older sim where sim_data version is out "
              "of date with current listeners used in the model)")

    # RIBOSOMES
    ribosome_50s_subunits = sim_data.process.complexation.get_monomers(
        sim_data.molecule_ids.s50_full_complex)
    ribosome_30s_subunits = sim_data.process.complexation.get_monomers(
        sim_data.molecule_ids.s30_full_complex)
    ribosome_subunits = (ribosome_50s_subunits["subunitIds"].tolist() +
                         ribosome_30s_subunits["subunitIds"].tolist())
    ribosome_monomers = {}
    for subunit in ribosome_subunits:
        if subunit in monomer_ids:
            ribosome_monomers[subunit] = get_monomer_indexes([subunit])[0]

    # RNAPs
    rnap_subunits = sim_data.process.complexation.get_monomers(
        sim_data.molecule_ids.full_RNAP)
    rnap_subunit_ids = rnap_subunits["subunitIds"].tolist()
    rnap_monomers = {}
    for subunit in rnap_subunit_ids:
        if subunit in monomer_ids:
            rnap_monomers[subunit] = get_monomer_indexes([subunit])[0]
        if subunit in complexation_complex_ids:
            monomer_ids_list = complexes_to_monomers[subunit]
            for mol_id in monomer_ids_list:
                if mol_id in monomer_ids:
                    rnap_monomers[mol_id] = get_monomer_indexes([mol_id])[0]

    # REPLISOMES
    replisome_trimer_subunits = sim_data.molecule_groups.replisome_trimer_subunits
    replisome_monomer_subunits = sim_data.molecule_groups.replisome_monomer_subunits
    replisome_subunit_ids = replisome_trimer_subunits + replisome_monomer_subunits
    replisome_monomers = {}
    for subunit in replisome_subunit_ids:
        if subunit in monomer_ids:
            replisome_monomers[subunit] = get_monomer_indexes([subunit])[0]
        elif subunit in complexation_complex_ids:
            monomer_ids_list = complexes_to_monomers[subunit]
            for mol_id in monomer_ids_list:
                if mol_id in monomer_ids:
                    replisome_monomers[mol_id] = get_monomer_indexes([mol_id])[0]

    # TRANSCRIPTION FACTORS
    tfs = sim_data.process.transcription_regulation.tf_ids
    tf_subunit_ids = [tf_id + f'[{sim_data.getter.get_compartment(tf_id)[0]}]'
                      for tf_id in tfs]
    tf_monomers = {}
    for subunit in tf_subunit_ids:
        if subunit in monomer_ids:
            tf_monomers[subunit] = get_monomer_indexes([subunit])[0]
        elif subunit in complexation_complex_ids:
            monomer_ids_list = complexes_to_monomers[subunit]
            for mol_id in monomer_ids_list:
                if mol_id in monomer_ids:
                    tf_monomers[mol_id] = get_monomer_indexes([mol_id])[0]
        elif subunit in equilibrium_complex_ids:
            monomer_ids_list = eq_complexes_to_monomers[subunit]
            for mol_id in monomer_ids_list:
                if mol_id in monomer_ids:
                    tf_monomers[mol_id] = get_monomer_indexes([mol_id])[0]
        elif subunit in two_component_system_complex_ids:
            monomer_ids_list = tcs_complexes_to_monomers[subunit]
            for mol_id in monomer_ids_list:
                if mol_id in monomer_ids:
                    tf_monomers[mol_id] = get_monomer_indexes([mol_id])[0]

    # METABOLISM (hover-only; not a plotted category)
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

    # monomer -> reactions it catalyzes / participates in, stored as
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
    """Classify a monomer into ALL functional-group categories it belongs to.

    Returns a single '<group> subunit' label when the monomer falls into exactly
    one group, or a ' + '-joined string when it falls into several (one legend
    entry/color per unique combination -- see build_categorized_figure), or
    'Monomer only' when it is in none. Metabolic involvement is intentionally NOT
    part of the label -- it is hover-only (see the module docstring)."""
    ribosome_monomers = protein_categories['ribosome_monomers']
    rnap_monomers = protein_categories['rnap_monomers']
    replisome_monomers = protein_categories['replisome_monomers']
    tf_monomers = protein_categories['tf_monomers']
    complexation_monomers = protein_categories['complexation_monomers']
    equilibrium_monomers = protein_categories['equilibrium_monomers']
    two_component_system_monomers = protein_categories['two_component_system_monomers']

    # The `[monomer_id] == monomer_idx` guard confirms the category dict's stored
    # index matches THIS monomer's plotted-array position -- a cheap consistency
    # check that the membership map and the plotted arrays are in the same id
    # space (they always should be, since both are keyed off monomer_dict).
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
    """True if any category in mol_type (single or ' + '-joined compound label) is
    a core-machinery subunit (ribosome / RNAP / replisome / transcription factor)
    -- these get the square marker; everything else stays a circle."""
    return any(label in SQUARE_MARKER_LABELS for label in mol_type.split(' + '))


def build_category_styles(molecule_types):
    """Assign a color/size/opacity/symbol/line_width to every unique category label
    present in the data, following a visual hierarchy that keeps the big, boring
    groups out of the way and makes the interesting ones pop:

      - 'Monomer only' (by far the biggest group): light grey, tiny, no outline --
        a faint background cloud.
      - 'Complexation complex subunit' (the next-biggest group): its fixed color
        but tiny and semi-transparent, also effectively background.
      - every other single category (ribosome / RNAP / replisome / transcription
        factor / equilibrium / TCS subunit): its fixed SINGLE_CATEGORY_SPECS color,
        larger, with a thin black outline so it stands out.
      - compound labels (proteins in more than one category -- the most
        interesting): their own color from COMPOUND_COLOR_PALETTE, largest, thin
        black outline.

    Core-machinery subunits (single or compound) are drawn as squares; everything
    else is a circle. Marker size/opacity are set by tier here; colors for single
    categories still come from SINGLE_CATEGORY_SPECS."""
    single_colors = {label: color for label, color, *_ in SINGLE_CATEGORY_SPECS}
    present_labels = sorted(set(molecule_types))

    styles = {}
    compound_i = 0
    for label in present_labels:
        symbol = 'square' if is_square_marker_label(label) else 'circle'
        if label == MONOMER_ONLY_LABEL:
            styles[label] = dict(
                color='lightgray', size=3, opacity=0.4,
                symbol='circle', line_width=0)
        elif label == BACKGROUND_COMPLEX_LABEL:
            styles[label] = dict(
                color=single_colors.get(label, 'goldenrod'), size=4, opacity=0.55,
                symbol='circle', line_width=0)
        elif ' + ' in label:
            color = COMPOUND_COLOR_PALETTE[compound_i % len(COMPOUND_COLOR_PALETTE)]
            compound_i += 1
            styles[label] = dict(
                color=color, size=11, opacity=0.9, symbol=symbol, line_width=1)
        else:
            styles[label] = dict(
                color=single_colors.get(label, 'black'), size=9, opacity=0.9,
                symbol=symbol, line_width=1)
    return styles


def build_half_life_map(sim_data):
    """Map monomer_id -> half life (minutes) from a sim's own degradation rates,
    so hover text can show the value each individual sim actually used (these can
    differ between sims, e.g. via protein half-life adjustments)."""
    monomer_data = sim_data.process.translation.monomer_data
    deg_rates = monomer_data["deg_rate"].asNumber(1 / units.s)
    half_life_min = {}
    for monomer_id, deg_rate in zip(monomer_data["id"], deg_rates):
        half_life_min[monomer_id] = (
            np.log(2) / deg_rate / 60.0 if deg_rate > 0 else np.inf
        )
    return half_life_min


def build_monomer_state_breakdown(
    sim_data, monomer_ids, bulk_by_id, active_means, bound_means
):
    """Per-monomer average count decomposed by state:
        {monomer_id: dict(free, complexation, equilibrium, tcs,
                          ribosome, rnap, replisome, bound)}
    computed from averaged bulk / active-unique / DNA-bound-TF data using each
    registry's stoich_matrix_monomers() (monomers-per-complex). The registry
    buckets use FREE complex bulk counts, so the buckets are disjoint from each
    other and from the free-monomer / unique-molecule / bound-TF buckets; their
    sum APPROXIMATES the exact monomer_counts total shown on the marker (small
    differences come from per-timestep integer rounding inside the listener)."""

    monomer_set = set(monomer_ids)

    def registry_data(molecule_names, complex_ids, smm):
        # subunit counts (monomers per complex) = max(0, -smm); rows =
        # molecule_names, cols = complex_ids. Returned so the same matrix feeds
        # BOTH the free-bulk buckets and the DNA-bound-TF distribution.
        molecule_names = list(molecule_names)
        subunit_counts = np.maximum(0, -np.asarray(smm))
        name_to_row = {m: i for i, m in enumerate(molecule_names)}
        return (molecule_names, list(complex_ids), subunit_counts, name_to_row)

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
    # TCS stoich_matrix_monomers() rows are modified_molecules (NOT molecule_names)
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
        # Monomers sequestered in each FREE complex of this registry:
        # subunit_counts @ (free bulk of each complex).
        _, complex_ids, subunit_counts, name_to_row = reg
        bulk_vec = np.array(
            [float(bulk_by_id.get(c, 0.0)) for c in complex_ids], dtype=float
        )
        return subunit_counts @ bulk_vec, name_to_row

    c_contrib, c_row = free_contrib(c_reg)
    eq_contrib, eq_row = free_contrib(eq_reg)
    if tcs_reg is not None:
        tcs_contrib, tcs_row = free_contrib(tcs_reg)
    else:
        tcs_contrib, tcs_row = None, {}

    # Ribosome / RNAP / replisome per-monomer stoichiometry (monomers per one
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
    trimer = list(sim_data.molecule_groups.replisome_trimer_subunits)
    monos = list(sim_data.molecule_groups.replisome_monomer_subunits)
    repl_stoich = {**{s: 3.0 for s in trimer}, **{s: 1.0 for s in monos}}

    # DNA-bound TF counts, keyed by the TF ACTIVE-FORM species id (compartment-
    # tagged). A TF's active form is usually a COMPLEX (e.g. a dimer), so the
    # bound count is distributed to that complex's monomers via stoichiometry --
    # exactly what the monomer_counts listener does. bound_means is keyed by the
    # bare tf id.
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
    # Direct case: the active form is itself a plotted monomer.
    for tagged, bc in bound_by_tagged.items():
        if tagged in monomer_set:
            bound_contrib[tagged] = bound_contrib.get(tagged, 0.0) + bc
    # Complex case: run the bound count through the registry it belongs to.
    for reg in (r for r in (c_reg, eq_reg, tcs_reg) if r is not None):
        _, complex_ids, subunit_counts, name_to_row = reg
        bvec = np.array(
            [bound_by_tagged.get(c, 0.0) for c in complex_ids], dtype=float
        )
        if not bvec.any():
            continue
        contrib_vec = subunit_counts @ bvec
        for m in monomer_ids:
            r = name_to_row.get(m)
            if r is not None and contrib_vec[r] > 0:
                bound_contrib[m] = bound_contrib.get(m, 0.0) + float(contrib_vec[r])

    ribo_active = float(active_means.get("ribosome", 0.0))
    rnap_active = float(active_means.get("rnap", 0.0))
    repl_active = float(active_means.get("replisome", 0.0))

    breakdown = {}
    for m in monomer_ids:
        breakdown[m] = dict(
            free=float(bulk_by_id.get(m, 0.0)),
            complexation=float(c_contrib[c_row[m]]) if m in c_row else 0.0,
            equilibrium=float(eq_contrib[eq_row[m]]) if m in eq_row else 0.0,
            tcs=(
                float(tcs_contrib[tcs_row[m]])
                if (tcs_contrib is not None and m in tcs_row)
                else 0.0
            ),
            ribosome=ribo_stoich.get(m, 0.0) * ribo_active,
            rnap=rnap_stoich.get(m, 0.0) * rnap_active,
            replisome=repl_stoich.get(m, 0.0) * repl_active,
            bound=float(bound_contrib.get(m, 0.0)),
        )
    return breakdown


def _fmt_id_list(ids):
    """Format a list of ids for hover, one per line, capping at MAX_HOVER_IDS."""
    ids = list(ids)
    if not ids:
        return ""
    shown = "<br>".join(f"&nbsp;&nbsp;- {x}" for x in ids[:MAX_HOVER_IDS])
    if len(ids) > MAX_HOVER_IDS:
        shown += f"<br>&nbsp;&nbsp;- (+{len(ids) - MAX_HOVER_IDS} more)"
    return shown


def _fmt_flux_val(v):
    """Format one flux value (mM/s) for hover; 'n/a' when the reaction has no FBA
    flux entry in that sim/direction."""
    if v is None:
        return "n/a"
    return f"{v:.2g}"


def _fmt_rxn_list(entries, flux_map_1, flux_map_2):
    """Format reaction entries [(base_rxn_id, via_complex_or_None), ...] for hover.

    Each line is "<base_rxn>[, via <complex>] (fwd S1:.. S2:.. | rev S1:.. S2:..
    mM/s)", with the forward flux looked up at base_rxn and the reverse at
    base_rxn + REVERSE_TAG in each sim's flux map. Reactions with no flux entry in
    either sim/direction show no flux parenthetical. Capped at MAX_HOVER_IDS."""
    entries = list(entries)
    if not entries:
        return ""
    lines = []
    for base_rxn, via in entries[:MAX_HOVER_IDS]:
        head = f"{base_rxn}, via {via}" if via else base_rxn
        rev_id = base_rxn + REVERSE_TAG
        fwd1, rev1 = flux_map_1.get(base_rxn), flux_map_1.get(rev_id)
        fwd2, rev2 = flux_map_2.get(base_rxn), flux_map_2.get(rev_id)
        if any(v is not None for v in (fwd1, rev1, fwd2, rev2)):
            flux_str = (
                f" (fwd S1:{_fmt_flux_val(fwd1)} S2:{_fmt_flux_val(fwd2)}"
                f" | rev S1:{_fmt_flux_val(rev1)} S2:{_fmt_flux_val(rev2)} {FLUX_UNIT})"
            )
        else:
            flux_str = ""
        lines.append(f"&nbsp;&nbsp;- {head}{flux_str}")
    if len(entries) > MAX_HOVER_IDS:
        lines.append(f"&nbsp;&nbsp;- (+{len(entries) - MAX_HOVER_IDS} more)")
    return "<br>".join(lines)


def membership_hover_lines(monomer_id, protein_categories, flux_map_1, flux_map_2):
    """Build the per-protein hover lines naming the specific complex memberships
    and metabolic reactions (with per-sim, per-direction flux). Returns '' when
    there is nothing to add."""
    lines = []
    complexation = protein_categories['monomer_to_complexation_complexes'].get(monomer_id, [])
    equilibrium = protein_categories['monomer_to_equilibrium_complexes'].get(monomer_id, [])
    tcs = protein_categories['monomer_to_tcs_complexes'].get(monomer_id, [])
    catalyzed = protein_categories['monomer_to_catalyzed_rxns'].get(monomer_id, [])
    reactant = protein_categories['monomer_to_reactant_rxns'].get(monomer_id, [])

    if complexation:
        lines.append(f"Complexation complex(es):<br>{_fmt_id_list(complexation)}")
    if equilibrium:
        lines.append(f"Equilibrium complex(es):<br>{_fmt_id_list(equilibrium)}")
    if tcs:
        lines.append(f"TCS complex(es):<br>{_fmt_id_list(tcs)}")
    if catalyzed:
        lines.append(
            f"Catalyzes reaction(s):<br>{_fmt_rxn_list(catalyzed, flux_map_1, flux_map_2)}")
    if reactant:
        lines.append(
            f"Reactant/product in reaction(s):<br>{_fmt_rxn_list(reactant, flux_map_1, flux_map_2)}")

    return ("<br>" + "<br>".join(lines)) if lines else ""


def _fmt_conc(count, ctm):
    """Format the average concentration of a count for hover, as ' (<x> µM)'.

    ctm is that sim's average countsToMolar in mM/count (mmol/L per molecule);
    count * ctm is a concentration in mM, and *1e3 converts to µM. Returns '' when
    ctm is None (listener unavailable)."""
    if ctm is None:
        return ""
    return f" ({count * ctm * 1e3:.3g} µM)"


def _fmt_doubling(dt):
    """Format an average doubling time (minutes) for the title parenthetical, as
    ', <x> min avg doubling time'. Returns '' when dt is None."""
    if dt is None:
        return ""
    return f", {dt:.1f} min avg doubling time"


def monomer_breakdown_hover(bd, ctm=None):
    """Format one monomer's average per-state breakdown for hover (nonzero only),
    each state annotated with its average concentration in µM (see _fmt_conc)."""
    order = [
        ("Free monomer", "free"),
        ("In complexation complexes", "complexation"),
        ("In equilibrium complexes", "equilibrium"),
        ("In TCS complexes", "tcs"),
        ("In ribosomes", "ribosome"),
        ("In RNAPs", "rnap"),
        ("In replisomes", "replisome"),
        ("DNA-bound (as TF)", "bound"),
    ]
    parts = [
        f"&nbsp;&nbsp;{label}: {bd[key]:.1f}{_fmt_conc(bd[key], ctm)}"
        for label, key in order
        if bd.get(key, 0.0) > 0.05
    ]
    if not parts:
        parts = ["&nbsp;&nbsp;(all states ~0)"]
    return "<br>".join(parts)


def _add_parity_and_stats(fig, sim1_log, sim2_log, r_value, pearson_r2, cod_r2):
    """Add the shared y=x parity line and the Pearson/COD stats annotation."""
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
        f"<b>Statistics:</b><br>"
        f"Pearson r = {r_value:.3f}<br>"
        f"Pearson R² = {pearson_r2:.3f}<br>"
        f"COD R² = {cod_r2:.3f}"
    )
    fig.add_annotation(
        x=0.95,
        y=0.05,
        xref='paper',
        yref='paper',
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
    """Apply the shared bounded-box layout: a fixed 1100x600 rectangle with a
    white plot background and a full black axis border (showline + mirror on both
    axes), left-aligned (multi-line) title, and the default legend placement
    (plotly auto-reserves right-margin space for it, so a long legend reads
    cleanly and scrolls rather than spilling over the plot). The axes are NOT
    forced square -- the box stays compact and the legend stays readable."""
    fig.update_layout(
        title=dict(text=title, x=0.01, xanchor='left', font=dict(size=14)),
        xaxis_title=xaxis_title,
        yaxis_title=yaxis_title,
        autosize=False,
        width=1200,
        height=720,
        margin=dict(t=160, b=80, l=90),
        template='plotly_white',
        plot_bgcolor='white',
        hovermode='closest',
        showlegend=True,
    )
    fig.update_xaxes(showline=True, linewidth=1, linecolor='black', mirror=True)
    fig.update_yaxes(showline=True, linewidth=1, linecolor='black', mirror=True)


def build_categorized_figure(sim1_log, sim2_log, molecule_types, hover_texts,
                             r_value, pearson_r2, cod_r2,
                             title, xaxis_title, yaxis_title):
    """Scatter with every protein colored by its functional category (no
    highlighting). Every unique category combination present in the data gets its
    own trace/legend entry/color (see build_category_styles); core-machinery
    subunits are drawn as squares, everything else as circles.

    Traces are added largest-category-first, smallest-last: Plotly draws later
    traces on top and lists the legend in trace-addition order, so this puts the
    biggest groups in the background/top of the legend and the smallest (easiest
    to lose under a big point cloud) on top/bottom of the legend, visible rather
    than buried."""
    molecule_types_arr = np.array(molecule_types, dtype=object)
    styles = build_category_styles(molecule_types)

    counts_by_label = {
        label: int((molecule_types_arr == label).sum()) for label in styles
    }
    ordered_labels = sorted(styles, key=lambda label: (-counts_by_label[label], label))

    fig = go.Figure()
    for label in ordered_labels:
        mask = molecule_types_arr == label
        if mask.sum() == 0:
            continue
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

    _add_parity_and_stats(fig, sim1_log, sim2_log, r_value, pearson_r2, cod_r2)
    _apply_pdp_layout(fig, title, xaxis_title, yaxis_title)
    return fig


class Plot(comparisonAnalysisPlot.ComparisonAnalysisPlot):

    def setup(self, inputDir: str) -> Tuple[
        AnalysisPaths, SimulationDataEcoli, ValidationDataEcoli]:
        """Return objects used for analyzing a single sim."""
        ap = AnalysisPaths(inputDir, variant_plot=True)
        sim_data = self.read_sim_data_file(inputDir)
        validation_data = self.read_validation_data_file(inputDir)
        return ap, sim_data, validation_data

    def get_cell_paths(self, ap):
        """Cell paths for all generations after the first SKIP_INITIAL_GENERATIONS
        of each seed -- the single cell set every per-sim read below averages over."""
        return ap.get_cells(
            generation=np.arange(SKIP_INITIAL_GENERATIONS, ap.n_generation)
        )

    def read_monomer_means(self, cell_paths):
        """Return (monomer_ids, mean_count_per_monomer). Counts are the mean over
        all timepoints of all cells; monomer_ids are the MonomerCounts listener's
        own id labels (its emit order)."""
        avg = read_stacked_columns(
            cell_paths, "MonomerCounts", "monomerCounts", ignore_exception=True
        ).mean(axis=0)
        reader = TableReader(
            os.path.join(cell_paths[0], "simOut", "MonomerCounts")
        )
        monomer_ids = list(reader.readAttribute("monomerIds"))
        return monomer_ids, np.asarray(avg, dtype=float)

    def read_bulk_by_id(self, cell_paths, ids):
        """Average bulk count (over all timepoints) for each requested id that is
        present in this sim's BulkMolecules listener. Returns {id: mean}."""
        reader = TableReader(os.path.join(cell_paths[0], "simOut", "BulkMolecules"))
        bulk_set = set(reader.readAttribute("objectNames"))
        wanted = [i for i in ids if i in bulk_set]
        if not wanted:
            return {}
        (arr,) = read_stacked_bulk_molecules(
            cell_paths, wanted, ignore_exception=True
        )
        arr = np.asarray(arr, dtype=float)
        if arr.ndim == 1:
            arr = arr.reshape(-1, 1)
        means = arr.mean(axis=0)
        return {wanted[i]: float(means[i]) for i in range(len(wanted))}

    def read_active_means(self, cell_paths):
        """Average active-unique-molecule counts: {'ribosome'|'rnap'|'replisome':
        mean} from the UniqueMoleculeCounts listener."""
        try:
            umc = read_stacked_columns(
                cell_paths, "UniqueMoleculeCounts", "uniqueMoleculeCounts",
                ignore_exception=True,
            )
            reader = TableReader(
                os.path.join(cell_paths[0], "simOut", "UniqueMoleculeCounts"))
            ids = list(reader.readAttribute("uniqueMoleculeIds"))
        except Exception as e:
            print(f"NOTE: could not read UniqueMoleculeCounts: {e!r}")
            return {}
        if umc.size == 0:
            return {}
        mean = umc.mean(axis=0)
        out = {}
        for key, name in (
            ("active_ribosome", "ribosome"),
            ("active_RNAP", "rnap"),
            ("active_replisome", "replisome"),
        ):
            if key in ids:
                out[name] = float(mean[ids.index(key)])
        return out

    def read_bound_means(self, cell_paths, sim_data):
        """Average DNA-bound TF counts {bare_tf_id: mean} from the RnaSynthProb
        listener's nActualBound column (emitted in that sim's tf_ids order)."""
        tf_ids = list(sim_data.process.transcription_regulation.tf_ids)
        try:
            nab = read_stacked_columns(
                cell_paths, "RnaSynthProb", "nActualBound", ignore_exception=True
            )
        except Exception as e:
            print(f"NOTE: could not read nActualBound: {e!r}")
            return {}
        if nab.size == 0:
            return {}
        mean = nab.mean(axis=0)
        return {tf_ids[i]: float(mean[i]) for i in range(min(len(tf_ids), mean.shape[0]))}

    def read_counts_to_molar(self, cell_paths):
        """Average countsToMolar (mM/count) from the EnzymeKinetics listener, or
        None if unavailable."""
        try:
            arr = read_stacked_columns(
                cell_paths, "EnzymeKinetics", "countsToMolar", ignore_exception=True
            )
        except Exception as e:
            print(f"NOTE: could not read countsToMolar: {e!r}")
            return None
        return float(arr.mean()) if arr.size else None

    def read_doubling_time(self, cell_paths):
        """Average doubling time (minutes) over the same cells the counts are
        averaged over, but equal-weight PER CELL: each cell's doubling time =
        finalTime - initialTime (Main/time), then averaged across cells."""
        try:
            per_cell = read_stacked_columns(
                cell_paths, "Main", "time",
                fun=lambda x: np.array([[x[-1, 0] - x[0, 0]]]),
                ignore_exception=True,
            )
        except Exception as e:
            print(f"NOTE: could not read doubling time: {e!r}")
            return None
        return float(per_cell.mean()) / 60.0 if per_cell.size else None

    def read_flux_map(self, cell_paths):
        """Average FBA reaction flux (mM/s) {reaction_id: mean} from the FBAResults
        listener's reactionFluxes column (ids in reactionIDs). Reaction ids include
        both forward and ' (reverse)'-tagged entries."""
        try:
            flux = read_stacked_columns(
                cell_paths, "FBAResults", "reactionFluxes", ignore_exception=True
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

    def do_plot(self, reference_sim_dir, plotOutDir, plotOutFileName,
                input_sim_dir, unused, metadata):
        # Sim 1 = reference (x-axis), Sim 2 = input (y-axis).
        ap1, sim_data1, _ = self.setup(reference_sim_dir)
        ap2, sim_data2, _ = self.setup(input_sim_dir)

        if ap1.n_generation <= 2 or ap2.n_generation <= 2:
            print("Skipping analysis -- not enough sims run.")
            return

        exp_id_1 = reference_sim_dir.split("out/")[-1].rstrip("/")
        exp_id_2 = input_sim_dir.split("out/")[-1].rstrip("/")
        print(f"Comparing {exp_id_1} (Sim 1; x-axis) vs {exp_id_2} (Sim 2; y-axis)")

        cell_paths_1 = self.get_cell_paths(ap1)
        cell_paths_2 = self.get_cell_paths(ap2)
        n_cells_1 = len(cell_paths_1)
        n_cells_2 = len(cell_paths_2)

        # Per-sim mean monomer counts keyed by that sim's own monomer id. The two
        # sims can differ in their monomer set/order, so counts are never indexed
        # positionally across sims; only ids present in both are plotted.
        monomer_ids_1, sim1_all = self.read_monomer_means(cell_paths_1)
        monomer_ids_2, sim2_all = self.read_monomer_means(cell_paths_2)
        print(f"Sim 1 has {n_cells_1} cells; Sim 2 has {n_cells_2} cells")
        print(f"Sim 1 total monomers: {len(monomer_ids_1)}")
        print(f"Sim 2 total monomers: {len(monomer_ids_2)}")

        means_1 = dict(zip(monomer_ids_1, (float(x) for x in sim1_all)))
        means_2 = dict(zip(monomer_ids_2, (float(x) for x in sim2_all)))

        # Half-life maps (per sim; deg rates can differ between sims):
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

        # Categorize proteins (from Sim 1's sim_data; both sims share the same
        # reconstruction-level category structure, and only ids present in both
        # are plotted above):
        protein_categories = get_protein_categories(
            sim_data1, monomer_ids, monomer_dict
        )

        # Per-state hover breakdown inputs (per sim). The plotted marker value is
        # UNCHANGED (the exact monomer_counts total); these reads only power the
        # hover decomposition of WHERE each monomer's average counts live.
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
        bulk_by_id_1 = self.read_bulk_by_id(cell_paths_1, monomer_ids + complex_ids_1)
        bulk_by_id_2 = self.read_bulk_by_id(cell_paths_2, monomer_ids + complex_ids_2)
        active_means_1 = self.read_active_means(cell_paths_1)
        active_means_2 = self.read_active_means(cell_paths_2)
        bound_means_1 = self.read_bound_means(cell_paths_1, sim_data1)
        bound_means_2 = self.read_bound_means(cell_paths_2, sim_data2)
        breakdown_1 = build_monomer_state_breakdown(
            sim_data1, monomer_ids, bulk_by_id_1, active_means_1, bound_means_1)
        breakdown_2 = build_monomer_state_breakdown(
            sim_data2, monomer_ids, bulk_by_id_2, active_means_2, bound_means_2)

        # Average counts_to_molar (mM/count) per sim, to annotate counts with µM:
        ctm_1 = self.read_counts_to_molar(cell_paths_1)
        ctm_2 = self.read_counts_to_molar(cell_paths_2)

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
        cod_r2 = r2_score(sim2_log, sim1_log)

        # Highlight list -> compartment-tagged monomer ids (matching plotted ids):
        converted_proteins = [
            tf_id + f'[{sim_data1.getter.get_compartment(tf_id)[0]}]'
            for tf_id in PLOT_PROTEINS_OF_INTEREST
        ]
        highlighted_set = set(converted_proteins)
        background_mask = np.array(
            [mid not in highlighted_set for mid in monomer_ids]
        )

        # Hover text (with per-state breakdown, concentrations, and the specific
        # complex / reaction memberships appended):
        hover_texts = []
        for i, monomer_id in enumerate(monomer_ids):
            half_life_1 = half_life_map_1.get(monomer_id)
            half_life_2 = half_life_map_2.get(monomer_id)
            half_life_1_str = (
                f"{half_life_1:.1f} min" if half_life_1 is not None else "N/A"
            )
            half_life_2_str = (
                f"{half_life_2:.1f} min" if half_life_2 is not None else "N/A"
            )
            gene_symbol = monomer_id_to_gene_symbol.get(monomer_id, "") or "N/A"
            hover_text = (
                f"<b>{monomer_id}</b><br>"
                f"Gene: {gene_symbol}<br>"
                f"Category: {molecule_types[i]}<br>"
                f"Sim 1 count: {sim1_avg[i]:.1f}{_fmt_conc(sim1_avg[i], ctm_1)}<br>"
                f"{monomer_breakdown_hover(breakdown_1[monomer_id], ctm_1)}<br>"
                f"Sim 2 count: {sim2_avg[i]:.1f}{_fmt_conc(sim2_avg[i], ctm_2)}<br>"
                f"{monomer_breakdown_hover(breakdown_2[monomer_id], ctm_2)}<br>"
                f"Sim 1 log: {sim1_log[i]:.3f}<br>"
                f"Sim 2 log: {sim2_log[i]:.3f}<br>"
                f"Sim 1 half life: {half_life_1_str}<br>"
                f"Sim 2 half life: {half_life_2_str}"
                f"{membership_hover_lines(monomer_id, protein_categories, flux_map_1, flux_map_2)}"
            )
            hover_texts.append(hover_text)

        # Plot 1: highlighted proteins of interest, each its OWN color/legend entry.
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
            legend_name = f"{symbol} ({tagged_id})" if symbol else tagged_id
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
        _add_parity_and_stats(fig, sim1_log, sim2_log, r_value, pearson_r2, cod_r2)

        highlighted_title = (
            f'Total Protein Count Comparison<br>'
            f'<sub>Sim 1 (x): {exp_id_1} (avg over {n_cells_1} cells'
            f'{_fmt_doubling(doubling_1)})<br>'
            f'Sim 2 (y): {exp_id_2} (avg over {n_cells_2} cells'
            f'{_fmt_doubling(doubling_2)})<br>'
            f'{len(monomer_ids)}/{n_total} proteins plotted | '
            f'skipped {SKIP_INITIAL_GENERATIONS} gens from each seed when '
            f'calculating the averages</sub>'
        )
        _apply_pdp_layout(
            fig, highlighted_title,
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

        # Plot 2: colored by protein category (no highlighting).
        categorized_title = (
            f'Total Protein Count Comparison<br>'
            f'<sub>Sim 1 (x): {exp_id_1} ({n_cells_1} cells'
            f'{_fmt_doubling(doubling_1)})<br>'
            f'Sim 2 (y): {exp_id_2} ({n_cells_2} cells'
            f'{_fmt_doubling(doubling_2)})<br>'
            f'{len(monomer_ids)}/{n_total} proteins plotted | '
            f'skipped {SKIP_INITIAL_GENERATIONS} gens from each seed when '
            f'calculating the averages</sub>'
        )
        fig_categorized = build_categorized_figure(
            sim1_log, sim2_log, molecule_types, hover_texts,
            r_value, pearson_r2, cod_r2,
            categorized_title,
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
