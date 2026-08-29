"""
Comparison RNA synthesis probability scatter plotly

For every transcription unit (TU) present in both sims, the average *actual* RNA
synthesis probability (`RnaSynthProb/actual_rna_synth_prob`, indexed by the
listener attribute `rnaIds`) of Sim 1 is plotted on x against Sim 2 on y. The
"target" synthesis probability is shown in hover text only.

Each sim resolves its TU ids from its own listener `rnaIds` attribute, and values
are aligned across sims by id (never positionally) in case the two sims carry
different TU sets/orders. Only ids present in both are plotted; the diagnostics
printed by do_plot() report the counts and any TUs unique to one sim.

Every per-TU quantity (actual/target synth prob, and the per-(TU, TF) bound-TF
counts) is the mean over all timepoints of all cells after the first
SKIP_INITIAL_GENERATIONS generations of each seed. `nActualBound` (per TF) is
summed over that same window to decide which TFs are active (see below).

Synthesis probabilities are small non-negative fractions that sum to ~1 across
all ~n_TU TUs, so individual values are tiny (~1e-4 typical) and span several
orders of magnitude. Plotting them linearly would crush almost every point into
the bottom-left corner. So instead, this plots log10 of the raw probability with a
small floor (SYNTH_PROB_FLOOR) applied first so that TUs that are effectively
never transcribed (avg prob 0) land at a finite left/bottom edge rather than
-inf.

The categorized plot:
Each TU is grouped by the set of transcription factors that regulate it,
considering only TFs that are active in the (reference) sim:

  * "Active TF" is defined data-drivenly as any TF that is bound to at least one
    promoter at some point in Sim 1 (its summed `RnaSynthProb/nActualBound` over
    the averaging window is > 0). This directly reflects the sim's condition
    (binding probabilities are condition-dependent) and sidesteps the ambiguity
    of the nominal `condition_active_tfs` list, which in e.g. the basal condition
    is empty even though many TFs still bind. Categories are defined from the
    reference sim (Sim 1), and both sims are assumed to be run in the same condition.
    # TODO (mia): check if the other conditions actually use the different
    calculated synth probs

  * Which TFs regulate a TU comes from the parCa delta-prob structure
    (`transcription_regulation.delta_prob`): its (deltaI = TU index, deltaJ = TF
    index) nonzero pattern lists, per TU, every regulating TF. TU rows are in
    `transcription.rna_data["id"]` order (== the listener's `rnaIds`). TF columns
    are in `transcription_regulation.tf_ids` order (== the listener's `tf_ids`).

  * Groups: TUs regulated by none of the active TFs -> "No TF" (grey).
    TUs regulated by exactly one active TF -> that TF's common name. TUs
    regulated by several active TFs -> one unique entry per unique combination,
    e.g. "ArgR + Lrp", "ArgR + LexA + Lrp".

  * TF common names are resolved by mapping each internal TF id (the "active TF"
    column, e.g. "CPLX0-228") to its abbreviation via
    `transcription_regulation.tf_to_gene_id` (e.g. "argR") and capitalising the
    first letter ("ArgR"). Falls back to the raw id if unmapped. See
    tf_condition.tsv for the list of active TF forms/model IDs.

Highlight plot:
Like the protein plot's proteins-of-interest plot: every TU is gray,
and any TU id listed in the module-level PLOT_TUS_OF_INTEREST is drawn in red on
top. Leave the list empty to colour everything uniformly.

TF focused plot:
Unlike the categorized plot (where each TU appears exactly once, in its unique
active-TF-combination group), this plot draws one trace per active TF containing
every TU that TF regulates. A TU regulated by several active TFs is drawn once per
regulating TF, so points intentionally overlap at the same (x, y). TUs
regulated by no active TF are shown as a faint gray "No active TF" background
group. See build_per_tf_figure.

Hover text information (same in all three graphs):
  1. TU name (the rnaId).
  2. The gene symbol(s) it encodes (a TU can be polycistronic): TU -> constituent
     cistrons via `transcription.cistron_tu_mapping_matrix` -> gene symbol via
     `replication.gene_data` (cistron_id -> symbol).
  3. Average actual and target RNA synth prob for Sim 1 and Sim 2 (the axes
     compare the actual value only though).
  4. The parCa `basal_prob` for the TU (reference sim).
  5. The parCa default `delta_prob` for each active regulating TF, by common name.
  6. The "% promoters bound" for each active regulating TF, per sim: the
     window-averaged number of that TF bound to the TU's promoters (time-averaged
     `n_bound_TF_per_TU`) over the window-averaged number of promoter copies of
     that TU (time-averaged `promoter_copy_number`), shown as
     "bound/available promoter copies (%)". The denominator is the TU's
     instantaneous promoter/gene copy number, which already accounts for DNA
     replication (a replicated locus has >1 promoter copy), so the ratio is the
     fraction of available promoter-copy opportunities that this TF actually
     occupied (the ratio of the two means, not the mean of per-timestep ratios).

Important note:
TU <-> gene, TU <-> TF and basal/delta values are all taken from the reference
sim's reconstruction (both sims share the same reconstruction structure and only
ids present in both are plotted).
"""

import os
from collections import defaultdict
from typing import Tuple

import numpy as np
import plotly.express as px
import plotly.graph_objects as go
from scipy.stats import pearsonr
from sklearn.metrics import r2_score

from models.ecoli.analysis import comparisonAnalysisPlot
from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from reconstruction.ecoli.simulation_data import SimulationDataEcoli
from validation.ecoli.validation_data import ValidationDataEcoli
from wholecell.analysis.analysis_tools import read_stacked_columns
from wholecell.io.tablereader import TableReader


# Cap for long id lists shown in hover text (show first N, then "(+M more)"):
MAX_HOVER_IDS = 8

# Number of initial generations of each seed to skip when averaging over cells:
SKIP_INITIAL_GENERATIONS = 2

# delta that gets added to each TU's average synth prob before taking log10:
SYNTH_PROB_FLOOR = 1e-8

# TUs (rnaIds) to highlight in red in the highlighted plot:
plot_set = lexA_TUs = ["TU0-12809[c]", "TU0-12849[c]", "TU0-12853[c]", "TU0-12854[c]", "TU0-13128[c]", "TU0-13438[c]", "TU0-13487[c]", "TU0-13744[c]", "TU0-14002[c]", "TU0-14047[c]", "TU0-1422[c]", "TU0-1423[c]", "TU0-14368[c]", "TU0-14447[c]", "TU0-14460[c]", "TU0-14679[c]", "TU0-14730[c]", "TU0-2[c]", "TU0-3941[c]", "TU0-4525[c]", "TU0-45271[c]", "TU0-6301[c]", "TU0-6559[c]", "TU0-6563[c]", "TU0-6661[c]", "TU0-6683[c]", "TU0-6686[c]", "TU0-7043[c]", "TU0-7581[c]", "TU0-8261[c]", "TU0-8264[c]", "TU0-8265[c]", "TU0-8267[c]", "TU00037[c]", "TU00062[c]", "TU00066[c]", "TU00069[c]", "TU00180[c]", "TU00352[c]", "TU365[c]", "EG10214_RNA[c]", "EG10341_RNA[c]", "EG10344_RNA[c]", "EG10604_RNA[c]", "EG10619_RNA[c]", "EG10620_RNA[c]", "EG10621_RNA[c]", "EG10622_RNA[c]", "EG10623_RNA[c]", "EG11086_RNA[c]"]
PLOT_TUS_OF_INTEREST = plot_set

# Fixed styling for the "No TF" group in the categorized plot.
NO_TF_LABEL = 'No TF'
NO_TF_COLOR = 'grey'

# TODO (mia): maybe make these ascend?
# Colors for single-active-TF groups:
SINGLE_TF_COLOR_PALETTE = (
    px.colors.qualitative.Dark24
    + px.colors.qualitative.Light24
)

# Colors for compound (multi-TF) groups, assigned in sorted-label order:
COMPOUND_COLOR_PALETTE = (
    px.colors.qualitative.Plotly
    + px.colors.qualitative.Set2
    + px.colors.qualitative.Pastel
)


def _tf_common_name(tf_id, tf_to_gene_id):
    """
    Map an internal TF id to a human-readable common name via the abbreviation
    in tf_to_gene_id. Falls back to the raw id when no abbreviation is available.
    """
    abbr = tf_to_gene_id.get(tf_id, '')
    if not abbr:
        return tf_id
    return abbr[0].upper() + abbr[1:]


def _fmt_id_list(ids):
    """
    Format a list of ids/strings for hover, one per line, capping at
    MAX_HOVER_IDS with a '(+M more)' overflow line.
    """
    ids = list(ids)
    if not ids:
        return ""
    shown = "<br>".join(f"&nbsp;&nbsp;- {x}" for x in ids[:MAX_HOVER_IDS])
    if len(ids) > MAX_HOVER_IDS:
        shown += f"<br>&nbsp;&nbsp;- (+{len(ids) - MAX_HOVER_IDS} more)"
    return shown


def _fmt_pct_bound(n_bound, n_copies):
    """
    Format one sim's per-(TU, TF) occupancy as 'bound/copies (pct%)' for hover.

    n_bound is the window-averaged n_bound_TF_per_TU for this (TU, TF); n_copies is
    the window-averaged promoter_copy_number availability for this TU. The
    percentage is the ratio of the two means (not the mean of per-timestep ratios)
    over the whole averaging window, of all the (promoter copy x timestep)
    opportunities for this TU, i.e. what fraction had this TF bound. Defaults to
    'N/A' (missing numerator) or 'bound/N/A' (missing/zero denominator) so the
    hover line still renders when a column was not emitted.
    """
    if n_bound is None:
        return "N/A"
    if not n_copies:  # None or 0.0 -> no usable denominator
        return f"{n_bound:.2f}/N/A"
    return f"{n_bound:.2f}/{n_copies:.2f} ({100.0 * n_bound / n_copies:.0f}%)"


def _add_parity_and_stats(fig, x_vals, y_vals, r_value, pearson_r2, cod_r2):
    """
    Add the shared y=x line and the Pearson/COD stats annotation.

    The three  stats are computed by the caller on the log values, not
    recomputed from raw counts here (since we care more about the magnitude
    correlation than the raw-prob correlation):
      - Pearson r: linear correlation of the two log-prob vectors (symmetric).
      - Pearson R² = r**2: variance explained by the best-fit line (NOT y = x).
      - COD R² via sklearn r2_score(y_true=sim2, y_pred=sim1): agreement with the
        y = x line. This is asymmetric (swapping the sims changes it), and can go
        negative when the sims disagree worse than predicting the mean would
        (so be careful interpreting it). NTS (mia): consider removing COD.
    """
    lo = float(min(x_vals.min(), y_vals.min()))
    hi = float(max(x_vals.max(), y_vals.max()))
    pad = 0.02 * (hi - lo) if hi > lo else 1.0
    fig.add_trace(go.Scatter(
        x=[lo - pad, hi + pad],
        y=[lo - pad, hi + pad],
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


def _apply_square_layout(fig, title, xaxis_title, yaxis_title):
    """
    Apply the shared square, plotly_white layout used by the highlighted plot,
    with the legend pinned inside the plot area (top-left).
    """
    fig.update_layout(
        title=dict(text=title, x=0.5, xanchor='center'),
        xaxis_title=xaxis_title,
        yaxis_title=yaxis_title,
        width=900,
        height=900,
        template='plotly_white',
        hovermode='closest',
        showlegend=True,
        legend=dict(
            x=0.02,
            y=0.98,
            bgcolor='rgba(255, 255, 255, 0.8)',
            bordercolor='gray',
            borderwidth=1
        )
    )
    fig.update_xaxes(scaleanchor="y", scaleratio=1, constrain='domain')
    fig.update_yaxes(constrain='domain')


def _apply_square_layout_legend_right(fig, title, xaxis_title, yaxis_title):
    """
    Like _apply_square_layout, but places the legend outside/right of the plot
    area. One legend entry per unique TF combination can produce many entries, so
    the figure is widened with a right margin to keep the legend from clipping.
    """
    fig.update_layout(
        title=dict(text=title, x=0.4, xanchor='center'),
        xaxis_title=xaxis_title,
        yaxis_title=yaxis_title,
        width=1200,
        height=900,
        margin=dict(r=300),
        template='plotly_white',
        hovermode='closest',
        showlegend=True,
        legend=dict(
            x=1.02,
            y=1,
            xanchor='left',
            yanchor='top',
            bgcolor='rgba(255, 255, 255, 0.8)',
            bordercolor='gray',
            borderwidth=1
        )
    )
    fig.update_xaxes(scaleanchor="y", scaleratio=1, constrain='domain')
    fig.update_yaxes(constrain='domain')


def build_category_styles(molecule_labels):
    """
    Assign a color/size/opacity/symbol to every unique category label present.

    "No TF" gets a fixed grey circle; each single-active-TF label gets a
    color from SINGLE_TF_COLOR_PALETTE (assigned in sorted-name order); each
    compound ('A + B') label gets its own color from COMPOUND_COLOR_PALETTE and a
    diamond marker (to visually separate multi-TF TUs). sorted() fixes both the
    palette assignment and the legend order so re-runs on the same active-TF
    set are deterministic.
    """
    present_labels = sorted(set(molecule_labels))
    styles = {}
    single_i = 0
    compound_i = 0
    for label in present_labels:
        if label == NO_TF_LABEL:
            styles[label] = dict(
                color=NO_TF_COLOR, size=6, opacity=0.5, symbol='circle')
        elif ' + ' in label:
            color = COMPOUND_COLOR_PALETTE[compound_i % len(COMPOUND_COLOR_PALETTE)]
            compound_i += 1
            styles[label] = dict(
                color=color, size=9, opacity=0.85, symbol='diamond')
        else:
            color = SINGLE_TF_COLOR_PALETTE[single_i % len(SINGLE_TF_COLOR_PALETTE)]
            single_i += 1
            styles[label] = dict(
                color=color, size=8, opacity=0.85, symbol='circle')
    return styles


def build_categorized_figure(x_vals, y_vals, molecule_labels, hover_texts,
                             r_value, pearson_r2, cod_r2,
                             title, xaxis_title, yaxis_title):
    """
    Scatter with every TU colored by its active-TF category.

    Traces are added largest-category-first, smallest-last: Plotly draws later
    traces on top and lists the legend in trace-addition order, so the biggest
    groups sit in the background/top of the legend and the smallest (easiest to
    lose under a big point cloud) sit on top/bottom of the legend, visible.
    """
    labels_arr = np.array(molecule_labels, dtype=object)
    styles = build_category_styles(molecule_labels)

    counts_by_label = {
        label: int((labels_arr == label).sum()) for label in styles
    }
    ordered_labels = sorted(styles, key=lambda label: (-counts_by_label[label], label))

    fig = go.Figure()
    for label in ordered_labels:
        mask = labels_arr == label
        if mask.sum() == 0:
            continue
        style = styles[label]

        # Get the indices where mask is True
        mask_indices = np.where(mask)[0]

        fig.add_trace(go.Scatter(
            x=x_vals[mask],
            y=y_vals[mask],
            mode='markers',
            marker=dict(
                color=style['color'],
                size=style['size'],
                opacity=style['opacity'],
                symbol=style['symbol'],
                line=dict(width=0),
            ),
            name=f"{label} (n={int(mask.sum())})",
            text=[hover_texts[i] for i in mask_indices],
            hovertemplate='%{text}<extra></extra>',
            showlegend=True
        ))

    _add_parity_and_stats(fig, x_vals, y_vals, r_value, pearson_r2, cod_r2)
    _apply_square_layout_legend_right(fig, title, xaxis_title, yaxis_title)
    return fig


def build_per_tf_figure(x_vals, y_vals, tu_ids, active_reg_by_tu, tf_names,
                        hover_texts, r_value, pearson_r2, cod_r2,
                        title, xaxis_title, yaxis_title,
                        include_no_tf_background=True):
    """
    Scatter with one trace per active TF holding every TU that TF regulates.

    This differs from build_categorized_figure, where each TU appears exactly once
    (in the trace for its unique combination of active TFs). Here a TU appears once
    per active TF that regulates it, so a TU regulated by several active TFs is
    drawn once in each of those TFs' traces (points intentionally overlap at the
    same (x, y)). The legend count for a TF is therefore the total number of
    plotted TUs that TF regulates (i.e. the exact number of TUs each TF is
    responsible for) regardless of what else co-regulates them (but co-regulation
    is still visible in the hover data info).

    active_reg_by_tu maps tu_id -> [active regulator tf indices] (the same per-TU
    active-regulator lists used for categorization/hover). tf_names[j] is the
    common name of TF index j. hover_texts[i] is the shared per-TU hover (it
    already lists every active TF regulating that TU, so the hover on a point is
    identical no matter which TF's trace it sits in).

    A TU regulated by no active TF is optionally shown as a faint gray "No active
    TF" background group (include_no_tf_background) so the full point cloud stays
    visible. Those TUs are never counted toward any TF's n=[].

    Colors are assigned to TFs in sorted-name order (deterministic across runs on
    the same active-TF set); traces are added largest-group-first so small groups
    draw on top and stay visible under the overlap.
    """
    # Per-active-TF membership among plotted TUs (a TU joins every regulator's
    # list), plus the "no active regulator" set for the optional background:
    tf_to_point_idxs = defaultdict(list)
    no_tf_idxs = []
    for i, tu_id in enumerate(tu_ids):
        regs = active_reg_by_tu.get(tu_id, [])
        if not regs:
            no_tf_idxs.append(i)
            continue
        for j in regs:
            tf_to_point_idxs[j].append(i)

    # Assign color per TF (sorted by common name), independent of trace draw
    # order below:
    active_js = sorted(tf_to_point_idxs.keys(), key=lambda j: tf_names[j])
    color_for_j = {
        j: SINGLE_TF_COLOR_PALETTE[k % len(SINGLE_TF_COLOR_PALETTE)]
        for k, j in enumerate(active_js)
    }

    fig = go.Figure()

    if include_no_tf_background and no_tf_idxs:
        fig.add_trace(go.Scatter(
            x=x_vals[no_tf_idxs],
            y=y_vals[no_tf_idxs],
            mode='markers',
            marker=dict(color='lightgray', size=5, opacity=0.3, line=dict(width=0)),
            name=f"No active TF (n={len(no_tf_idxs)})",
            text=[hover_texts[i] for i in no_tf_idxs],
            hovertemplate='%{text}<extra></extra>',
            showlegend=True
        ))

    # One trace per active TF, largest group first (so small groups sit on top):
    ordered_js = sorted(
        tf_to_point_idxs, key=lambda j: (-len(tf_to_point_idxs[j]), tf_names[j]))
    for j in ordered_js:
        idxs = tf_to_point_idxs[j]
        fig.add_trace(go.Scatter(
            x=x_vals[idxs],
            y=y_vals[idxs],
            mode='markers',
            marker=dict(color=color_for_j[j], size=7, opacity=0.6, line=dict(width=0)),
            name=f"{tf_names[j]} (n={len(idxs)})",
            text=[hover_texts[i] for i in idxs],
            hovertemplate='%{text}<extra></extra>',
            showlegend=True
        ))

    _add_parity_and_stats(fig, x_vals, y_vals, r_value, pearson_r2, cod_r2)
    _apply_square_layout_legend_right(fig, title, xaxis_title, yaxis_title)
    return fig


class Plot(comparisonAnalysisPlot.ComparisonAnalysisPlot):

    def setup(self, inputDir: str) -> Tuple[
        AnalysisPaths, SimulationDataEcoli, ValidationDataEcoli]:
        """
        Return objects used for analyzing a single sim.
        """
        ap = AnalysisPaths(inputDir, variant_plot=True)
        sim_data = self.read_sim_data_file(inputDir)
        validation_data = self.read_validation_data_file(inputDir)
        return ap, sim_data, validation_data

    def read_synth_prob_means(self, ap):
        """
        Return time/cell-averaged RnaSynthProb quantities for one sim.

        All averages are the mean over all timepoints of all cells after the
        first SKIP_INITIAL_GENERATIONS generations of each seed; nActualBound is
        summed over the same window. Ids come from the listener's own attributes.
        NTS (mia): this differs from how this plot averages in vEcoli.

        Returns a dict with:
          rna_ids           : list of TU ids (listener 'rnaIds' order)
          tf_ids            : list of TF ids (listener 'tf_ids' order)
          actual_avg        : (n_TU,) mean actual synth prob
          target_avg        : (n_TU,) mean target synth prob
          n_bound_avg       : (n_TU, n_TF) mean # of each TF bound per TU
          promoter_copy_avg : (n_TU,) mean promoter copy number per TU
          nbound_tf_total   : (n_TF,) summed nActualBound per TF (activity test)
          n_cells           : number of cells averaged
        """
        cell_paths = ap.get_cells(
            generation=np.arange(SKIP_INITIAL_GENERATIONS, ap.n_generation)
        )
        n_cells = len(cell_paths)

        reader = TableReader(
            os.path.join(cell_paths[0], "simOut", "RnaSynthProb")
        )
        rna_ids = list(reader.readAttribute("rnaIds"))
        tf_ids = list(reader.readAttribute("tf_ids"))
        n_TU = len(rna_ids)
        n_TF = len(tf_ids)

        actual_avg = read_stacked_columns(
            cell_paths, "RnaSynthProb", "actual_rna_synth_prob",
            ignore_exception=True).mean(axis=0)
        target_avg = read_stacked_columns(
            cell_paths, "RnaSynthProb", "target_rna_synth_prob",
            ignore_exception=True).mean(axis=0)

        # n_bound_TF_per_TU is stored flattened (n_TU * n_TF) per timestep; the
        # time-mean is reshaped back to (n_TU, n_TF):
        n_bound_flat = read_stacked_columns(
            cell_paths, "RnaSynthProb", "n_bound_TF_per_TU",
            ignore_exception=True).mean(axis=0)
        n_bound_avg = np.asarray(n_bound_flat, dtype=float).reshape(n_TU, n_TF)

        # promoter_copy_number is a 1-D per-TU column (the TU's instantaneous
        # promoter/gene copy number, which rises with DNA replication); it shares
        # rnaIds order with actual_rna_synth_prob. It is the denominator of the
        # per-(TU, TF) "% bound" hover detail:
        promoter_copy_avg = read_stacked_columns(
            cell_paths, "RnaSynthProb", "promoter_copy_number",
            ignore_exception=True).mean(axis=0)

        nbound_tf_total = read_stacked_columns(
            cell_paths, "RnaSynthProb", "nActualBound",
            ignore_exception=True).sum(axis=0)

        return dict(
            rna_ids=rna_ids,
            tf_ids=tf_ids,
            actual_avg=np.asarray(actual_avg, dtype=float),
            target_avg=np.asarray(target_avg, dtype=float),
            n_bound_avg=n_bound_avg,
            promoter_copy_avg=np.asarray(promoter_copy_avg, dtype=float),
            nbound_tf_total=np.asarray(nbound_tf_total, dtype=float),
            n_cells=n_cells,
        )

    def do_plot(self, reference_sim_dir, plotOutDir, plotOutFileName,
                input_sim_dir, unused, metadata):
        """
        Make the plots!!!
        """
        # Sim 1 = reference (x-axis), Sim 2 = input (y-axis).
        ap1, sim_data1, _ = self.setup(reference_sim_dir)
        ap2, sim_data2, _ = self.setup(input_sim_dir)

        if ap1.n_generation <= 2 or ap2.n_generation <= 2:
            print("Skipping analysis -- not enough sims run.")
            return

        # Obtain the simulation names:
        exp_id_1 = reference_sim_dir.split("out/")[-1].rstrip("/")
        exp_id_2 = input_sim_dir.split("out/")[-1].rstrip("/")
        print(f"Comparing {exp_id_1} (Sim 1; x-axis) vs {exp_id_2} (Sim 2; y-axis)")

        # Obtain time/cell-averaged RnaSynthProb data, keyed by each sim's own TU
        # id:
        data1 = self.read_synth_prob_means(ap1)
        data2 = self.read_synth_prob_means(ap2)
        print(f"Sim 1 has {data1['n_cells']} cells; Sim 2 has {data2['n_cells']} cells")
        print(f"Sim 1 total TUs: {len(data1['rna_ids'])}")
        print(f"Sim 2 total TUs: {len(data2['rna_ids'])}")

        actual_1 = dict(zip(data1['rna_ids'], data1['actual_avg']))
        actual_2 = dict(zip(data2['rna_ids'], data2['actual_avg']))
        target_1 = dict(zip(data1['rna_ids'], data1['target_avg']))
        target_2 = dict(zip(data2['rna_ids'], data2['target_avg']))
        # n_bound rows aligned to Sim 1's TU ids:
        n_bound_by_tu = dict(zip(data1['rna_ids'], data1['n_bound_avg']))
        # Calculate promoter copy number (the "% bound" denominator):
        promoter_copy_1 = dict(zip(data1['rna_ids'], data1['promoter_copy_avg']))
        promoter_copy_2 = dict(zip(data2['rna_ids'], data2['promoter_copy_avg']))

        plotted_ids = [tid for tid in data1['rna_ids'] if tid in actual_2]
        not_plotted_1 = [tid for tid in data1['rna_ids'] if tid not in actual_2]
        not_plotted_2 = [tid for tid in data2['rna_ids'] if tid not in actual_1]
        n_total = len(plotted_ids) + len(not_plotted_1) + len(not_plotted_2)
        print(f"Plotted TUs (present in both sims): {len(plotted_ids)}/{n_total}")
        if not_plotted_1:
            print(f"NOTE: {len(not_plotted_1)} TU(s) only in Sim 1, e.g. "
                  f"{not_plotted_1[:5]}")
        if not_plotted_2:
            print(f"NOTE: {len(not_plotted_2)} TU(s) only in Sim 2, e.g. "
                  f"{not_plotted_2[:5]}")

        tu_ids = plotted_ids
        sim1_actual = np.array([actual_1[t] for t in tu_ids])
        sim2_actual = np.array([actual_2[t] for t in tu_ids])
        sim1_target = np.array([target_1[t] for t in tu_ids])
        sim2_target = np.array([target_2[t] for t in tu_ids])

        # Assuming both sims have the same reconstruction structure, use the
        # reference sim's reconstruction structure (i.e. variable names used
        # to extract basal/delta, TF names, gene mapping) to obtain values:
        # TODO (mia): consider making each sim use its own reconstruction
        #  structure and just keying purely by TU ID?
        tr = sim_data1.process.transcription_regulation
        transcription = sim_data1.process.transcription

        tf_ids_reg = list(tr.tf_ids)  # == data1['tf_ids'] order
        tf_names = [_tf_common_name(tf, tr.tf_to_gene_id) for tf in tf_ids_reg]

        # Determine active TFs (defined as >=1 in Sim 1, should be 23 for the
        # default condition):
        nbound_by_tf = dict(zip(data1['tf_ids'], data1['nbound_tf_total']))
        active_tf = np.array(
            [nbound_by_tf.get(tf, 0.0) > 0 for tf in tf_ids_reg], dtype=bool)
        active_tf_names = [tf_names[j] for j in range(len(tf_ids_reg)) if active_tf[j]]
        print(f"Active TFs in Sim 1 (bound >=1x): {len(active_tf_names)} of "
              f"{len(tf_ids_reg)}: {sorted(active_tf_names)}")

        # TU index (rna_data order) for basal/delta alignment:
        rna_data_ids = list(transcription.rna_data["id"])
        tu_id_to_index = {tid: i for i, tid in enumerate(rna_data_ids)}
        basal_prob = np.asarray(tr.basal_prob, dtype=float)

        # Determine which regulating TFs are assigned to each TU and their
        # default delta_prob:
        dp = tr.delta_prob
        tu_to_tf_delta = defaultdict(dict)
        for ti, tj, tv in zip(dp['deltaI'], dp['deltaJ'], dp['deltaV']):
            tu_to_tf_delta[int(ti)][int(tj)] = float(tv)

        # Gather TU -> constituent cistrons -> gene symbols:
        cistron_tu_mapping = transcription.cistron_tu_mapping_matrix
        cistron_ids = list(transcription.cistron_data["id"])
        gene_data = sim_data1.process.replication.gene_data
        cistron_id_to_symbol = dict(zip(gene_data["cistron_id"], gene_data["symbol"]))

        # Obtain Sim 2 reconstruction structures, so basal_prob / delta_prob can be
        # shown per sim in hover (the two sims can still have different parca
        # output values but again this assumes the same structure/variable
        # naming scheme):
        tr2 = sim_data2.process.transcription_regulation
        transcription2 = sim_data2.process.transcription
        tf_id_to_idx2 = {tf: i for i, tf in enumerate(tr2.tf_ids)}
        rna_data_ids2 = list(transcription2.rna_data["id"])
        tu_id_to_index2 = {tid: i for i, tid in enumerate(rna_data_ids2)}
        basal_prob2 = np.asarray(tr2.basal_prob, dtype=float)
        dp2 = tr2.delta_prob
        tu_to_tf_delta2 = defaultdict(dict)
        for ti, tj, tv in zip(dp2['deltaI'], dp2['deltaJ'], dp2['deltaV']):
            tu_to_tf_delta2[int(ti)][int(tj)] = float(tv)

        # Compute avg-bound rows + TF-id -> listener column maps:
        n_bound_by_tu_2 = dict(zip(data2['rna_ids'], data2['n_bound_avg']))
        tf_id_to_col1 = {tf: i for i, tf in enumerate(data1['tf_ids'])}
        tf_id_to_col2 = {tf: i for i, tf in enumerate(data2['tf_ids'])}

        # Build each TU category label + hover text:
        molecule_labels = []
        hover_texts = []
        active_reg_by_tu = {}  # tu_id -> [active regulator tf idxs]
        for i, tu_id in enumerate(tu_ids):
            tu_idx = tu_id_to_index.get(tu_id)

            # Active regulating TFs of this TU (sorted by common name):
            reg = tu_to_tf_delta.get(tu_idx, {}) if tu_idx is not None else {}
            active_reg = sorted(
                (j for j in reg if active_tf[j]), key=lambda j: tf_names[j])
            active_reg_by_tu[tu_id] = active_reg

            if not active_reg:
                label = NO_TF_LABEL
            else:
                label = ' + '.join(tf_names[j] for j in active_reg)
            molecule_labels.append(label)

            # Gene symbols encoded by the TU (polycistronic -> several):
            if tu_idx is not None:
                cistron_idxs = cistron_tu_mapping.getcol(tu_idx).nonzero()[0]
                genes = [cistron_id_to_symbol.get(
                    cistron_ids[ci], cistron_ids[ci]) for ci in cistron_idxs]
            else:
                genes = []

            # Extract basal_prob values:
            tu_idx2 = tu_id_to_index2.get(tu_id)
            basal1_str = (f"{basal_prob[tu_idx]:.3e}"
                          if tu_idx is not None else "N/A")
            basal2_str = (f"{basal_prob2[tu_idx2]:.3e}"
                          if tu_idx2 is not None else "N/A")

            # For each active TF, obtain delta_prob and avg # bound:
            reg2 = tu_to_tf_delta2.get(tu_idx2, {}) if tu_idx2 is not None else {}
            n_bound_row1 = n_bound_by_tu.get(tu_id)
            n_bound_row2 = n_bound_by_tu_2.get(tu_id)
            delta_lines = []
            bound_lines = []
            for j in active_reg:
                tf_id = tf_ids_reg[j]
                # delta_prob
                d1 = reg.get(j)
                jj2 = tf_id_to_idx2.get(tf_id)
                d2 = reg2.get(jj2) if jj2 is not None else None
                d1s = f"{d1:+.3e}" if d1 is not None else "N/A"
                d2s = f"{d2:+.3e}" if d2 is not None else "N/A"
                delta_lines.append(
                    f"&nbsp;&nbsp;- {tf_names[j]}: Sim 1 {d1s}, Sim 2 {d2s}")
                # % promoters bound: avg n_bound_TF_per_TU for this
                # (TU, TF) over the avg promoter_copy_number for this TU, shown as
                # "bound/copies (%)":
                col1 = tf_id_to_col1.get(tf_id)
                col2 = tf_id_to_col2.get(tf_id)
                b1 = (n_bound_row1[col1]
                      if n_bound_row1 is not None and col1 is not None else None)
                b2 = (n_bound_row2[col2]
                      if n_bound_row2 is not None and col2 is not None else None)
                copies1 = promoter_copy_1.get(tu_id)
                copies2 = promoter_copy_2.get(tu_id)
                bound_lines.append(
                    f"&nbsp;&nbsp;- {tf_names[j]}: "
                    f"Sim 1 {_fmt_pct_bound(b1, copies1)}, "
                    f"Sim 2 {_fmt_pct_bound(b2, copies2)}")

            genes_str = _fmt_id_list(genes) if genes else "&nbsp;&nbsp;- N/A"
            delta_block = ("<br>" + "<br>".join(delta_lines)) if delta_lines else " N/A"
            bound_block = ("<br>" + "<br>".join(bound_lines)) if bound_lines else " N/A"

            hover_text = (
                f"<b>{tu_id}</b><br>"
                f"Gene(s):<br>{genes_str}<br>"
                f"Category: {label}<br>"
                f"<br>"
                f"<b>Sim 1:</b><br>"
                f"&nbsp;&nbsp;actual prob: {sim1_actual[i]:.3e}<br>"
                f"&nbsp;&nbsp;target prob: {sim1_target[i]:.3e}<br>"
                f"&nbsp;&nbsp;basal_prob: {basal1_str}<br>"
                f"<br>"
                f"<b>Sim 2:</b><br>"
                f"&nbsp;&nbsp;actual prob: {sim2_actual[i]:.3e}<br>"
                f"&nbsp;&nbsp;target prob: {sim2_target[i]:.3e}<br>"
                f"&nbsp;&nbsp;basal_prob: {basal2_str}<br>"
                f"<br>"
                f"<b>Delta prob per active TF:</b>{delta_block}<br>"
                f"<br>"
                f"<b>% promoters bound per active TF "
                f"(avg bound / avg copies):</b>{bound_block}"
            )
            hover_texts.append(hover_text)

        # Diagnostic: single-active-TF group sizes (sanity check vs expectations
        # like ArgR/Lrp/LexA):
        labels_arr = np.array(molecule_labels, dtype=object)
        single_counts = {
            name: int((labels_arr == name).sum())
            for name in active_tf_names
            if int((labels_arr == name).sum()) > 0
        }
        if single_counts:
            print("Single-active-TF group sizes (TUs regulated by only that TF):")
            for name in sorted(single_counts, key=lambda n: (-single_counts[n], n)):
                print(f"  {name}: n={single_counts[name]}")

        # Compute log10 of the actual synth prob for each sim (with the added
        # delta floor value):
        sim1_log = np.log10(np.maximum(sim1_actual, SYNTH_PROB_FLOOR))
        sim2_log = np.log10(np.maximum(sim2_actual, SYNTH_PROB_FLOOR))

        r_value = pearsonr(sim1_log, sim2_log)[0]
        pearson_r2 = r_value ** 2
        cod_r2 = r2_score(sim2_log, sim1_log)

        xaxis_title = 'log10(avg actual RNA synth prob) — Sim 1'
        yaxis_title = 'log10(avg actual RNA synth prob) — Sim 2'

        highlighted_mask = np.array([t in PLOT_TUS_OF_INTEREST for t in tu_ids])
        background_mask = ~highlighted_mask

        # Plot 1: highlighted TUs of interest
        fig = go.Figure()
        fig.add_trace(go.Scatter(
            x=sim1_log[background_mask],
            y=sim2_log[background_mask],
            mode='markers',
            marker=dict(color='lightseagreen', size=6, opacity=0.6, line=dict(width=0)),
            name='All TUs',
            text=[hover_texts[i] for i in range(len(hover_texts)) if background_mask[i]],
            hovertemplate='%{text}<extra></extra>',
            showlegend=True
        ))
        if highlighted_mask.sum() > 0:
            fig.add_trace(go.Scatter(
                x=sim1_log[highlighted_mask],
                y=sim2_log[highlighted_mask],
                mode='markers',
                marker=dict(
                    color='red', size=10, opacity=0.9, line=dict(width=1, color='darkred')
                ),
                name='TUs of interest',
                text=[hover_texts[i] for i in range(len(hover_texts)) if highlighted_mask[i]],
                hovertemplate='%{text}<extra></extra>',
                showlegend=True
            ))
        _add_parity_and_stats(fig, sim1_log, sim2_log, r_value, pearson_r2, cod_r2)

        highlighted_title = (
            f'Average Actual RNA Synthesis Probability Comparison<br>'
            f'<sub>Sim 1 (x): {exp_id_1} (avg over {data1["n_cells"]} cells) vs.'
            f'<br>Sim 2 (y): {exp_id_2} (avg over {data2["n_cells"]} cells) '
            f'<br>{len(tu_ids)}/{n_total} TUs plotted</sub>'
        )
        _apply_square_layout(fig, highlighted_title, xaxis_title, yaxis_title)
        output_filename = os.path.join(
            plotOutDir,
            f"{plotOutFileName}_sim1_{exp_id_1}_sim2_{exp_id_2}_additive_only.html",
        )
        fig.write_html(output_filename)
        print(f"Saved plot to {output_filename}")
        print(f"Highlighted TUs: {int(highlighted_mask.sum())} of "
              f"{len(PLOT_TUS_OF_INTEREST)} requested")

        # Plot 2: colored by active-TF category
        categorized_title = (
            f'Average Actual RNA Synthesis Probability Comparison<br>'
            f'<sub>Colored by active-TF binding'
            f'<br>Sim 1 (x): {exp_id_1} ({data1["n_cells"]} cells) vs.'
            f'<br>Sim 2 (y): {exp_id_2} ({data2["n_cells"]} cells) '
            f'<br>{len(tu_ids)}/{n_total} TUs plotted</sub>'
        )
        fig_categorized = build_categorized_figure(
            sim1_log, sim2_log, molecule_labels, hover_texts,
            r_value, pearson_r2, cod_r2,
            categorized_title, xaxis_title, yaxis_title
        )
        categorized_filename = os.path.join(
            plotOutDir,
            f"{plotOutFileName}_categorized_sim1_{exp_id_1}_sim2_{exp_id_2}.html",
        )
        fig_categorized.write_html(categorized_filename)
        print(f"Saved categorized plot to {categorized_filename}")

        # Plot 3: TUs regulated by each active TF (one trace per active TF):
        per_tf_title = (
            f'Average Actual RNA Synthesis Probability Comparison<br>'
            f'<sub>Sorted by active TF regulation'
            f'<br>Sim 1 (x): {exp_id_1} ({data1["n_cells"]} cells) vs.'
            f'<br>Sim 2 (y): {exp_id_2} ({data2["n_cells"]} cells)'
            f'<br>{len(active_tf_names)} active TFs | {len(tu_ids)}/{n_total} TUs plotted'
            f'</sub>'
        )
        fig_per_tf = build_per_tf_figure(
            sim1_log, sim2_log, tu_ids, active_reg_by_tu, tf_names, hover_texts,
            r_value, pearson_r2, cod_r2,
            per_tf_title, xaxis_title, yaxis_title
        )
        per_tf_filename = os.path.join(
            plotOutDir,
            f"{plotOutFileName}_per_tf_sim1_{exp_id_1}_sim2_{exp_id_2}.html",
        )
        fig_per_tf.write_html(per_tf_filename)
        print(f"Saved per-active-TF plot to {per_tf_filename}")


if __name__ == "__main__":
    Plot().cli()
