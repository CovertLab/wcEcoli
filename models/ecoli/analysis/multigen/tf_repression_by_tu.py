"""
tf_repression_by_tu matplotlib plot

Produces one figure per transcription unit (TU) listed in ``TU_IDS`` below, a
tall 5-panel stack sharing a time axis over the duration of a seed.

Background on how RNA synth prob is modeled:
``target_rna_synth_prob`` is a normalized multinomial distribution over all
TUs in the model that sums to 1.0 every timestep. It is each TU's share of
RNAP initiation events, not an absolute transcription rate. Thus, a TU's synth
prob can fall because other TUs rose, with no change to its own regulation.
Note that in the ParCa, a basal_prob is fit for each TU, and each TF that can
bind to that TU (if applicable) has a default delta_prob that is also fit.

Panels plotted:
  1. RNA synthesis probability: ``target_rna_synth_prob`` vs ``actual_rna_synth_prob``.
  2. Promoters bound per TF, as a stacked bar per timestep. Legend labels show
     each TF's default delta_prob calculated in the ParCa. The black outline
     marks the promoter-copy number available for this TU at each timestep
     (typically ranging from 1-2). The stack total can exceed that line since
     one promoter copy can be bound by several TFs at once (distinct sites),
     so note that Σ_TF n_bound is not capped by the copy number.
  3. Transcript counts: total/full/partial mRNA existance. The title reports
     the TU's mRNA half-life as well.
  4. ``max_p``: the per-timestep cap on any single promoter's share of the
     multinomial (the most initiations one promoter can physically fit in a
     timestep, expressed as a fraction of the RNAPs being activated). On the
     same normalized scale as the synth prob.
  5. Indication of RNAP overcrowding (as a binary present/absent strip). Plots
     both the real per-TU flag ``tu_is_overcrowded`` (red = this TU hit the cap)
     and the cell-wide ``actual !=target`` (gray). Note that the capping
     renormalization applies to the whole cell, so ``actual != target`` fires
     for a TU being plotted that may never be crowded itself here merely because
     something else in the cell did. Only the red flag means "this TU hit the
     cap". The panel title reports the occurance of both over the sim as %s.

When ppgpp_regulation is on, the basal_prob is recomputed from the current ppGpp
concentration every timestep, and the TF delta is scaled by that same new basal,
so per promoter: init_prob = B_ppgpp · (1 + Σ δ_norm·bound). Both terms scale
with B_ppgpp together — that is why fold repression is invariant to ppGpp. The
ParCa-fit reference numbers (basal_prob for this TU and each TF's default
delta_prob) are still shown in the figure title for reference.

# TODOs:
- add the availablity to plot tRNA and rRNA (have the mature vs immature tRNA counts shown instead
of full/partial/etc.) in panel 3
"""

import os
import pickle

from matplotlib import pyplot as plt
import numpy as np

from models.ecoli.analysis import multigenAnalysisPlot
from wholecell.analysis.analysis_tools import exportFigure, read_stacked_columns
from wholecell.io.tablereader import TableReader
from wholecell.utils import units

# USER INPUTS
# List the TU(s) to be plotted (with or without compartment tag):
TU_IDS = ["TU103[c]", "TU355[c]"]

# Short names for TFs to be used in the figure legend (otherwise the raw TF id
# is used, see tf_condition.tsv for the full list of TFs in the model):
# TODO (mia): put all names in
TF_SHORT_NAMES = {
    "CPLX0-228": "ArgR",
    "PC00010": "lexA",
    "MONOMER0-155": "Lrp",
    "CPLX0-226": "crp",
    "CPLX0-7669": "ArgP",
"PHOSPHO-ARCA": "ArcA",
"MONOMER0-160": "dnaA",
"CPLX0-7705":"fis",
    "PUTA-CPLX":"putA"
}

# Colors for each TF in the stacked bar plot (panel 2):
# TODO (mia): figure out the max # of TFs that regulate a single TU and make the options that long
TF_PALETTE = ["#d62728", "#1f77b4", "#2ca02c", "#9467bd", "#ff7f0e", "#8c564b"]

STEP = "steps-post"

# END USER INPUTS

class Plot(multigenAnalysisPlot.MultigenAnalysisPlot):
    def do_plot(self, seedOutDir, plotOutDir, plotOutFileName, simDataFile,
                validationDataFile, metadata):
        with open(simDataFile, 'rb') as f:
            sim_data = pickle.load(f)

        cell_paths = self.ap.get_cells()
        sim_out_dir = os.path.join(cell_paths[0], 'simOut')

        treg = sim_data.process.transcription_regulation
        rel = sim_data.relation
        trans = sim_data.process.transcription
        # Obtain delta_probs and basal_prob:
        dpm = treg.get_delta_prob_matrix(dense=True)  # (n_TU, n_TF)
        basal_prob = treg.basal_prob

        # TU mRNA half-life (min) = ln2 / deg_rate, for the transcript-count title:
        deg = trans.rna_data['deg_rate'].asNumber(1 / units.s)
        with np.errstate(divide='ignore', invalid='ignore'):
            hl_min = np.where(deg > 0, np.log(2) / deg / 60.0, np.nan)
        tu_half_life_min = {
            rid: float(hl) for rid, hl in zip(trans.rna_data['id'], hl_min)}

        # Extract listener data:
        rsp = TableReader(os.path.join(sim_out_dir, 'RnaSynthProb'))
        rna_ids = list(rsp.readAttribute('rnaIds'))
        listener_tf_ids = list(rsp.readAttribute('tf_ids'))
        n_tf = len(listener_tf_ids)
        rc = TableReader(os.path.join(sim_out_dir, 'RNACounts'))
        mrna_ids = list(rc.readAttribute('mRNA_ids'))

        # TF index maps: the listener flattens n_bound_TF_per_TU as (n_TU, n_TF)
        # row-major, so the flat column for (tu, tf) is tu_idx * n_TF + tf_idx
        # using the LISTENER tf order. The delta_prob matrix columns use the
        # sim_data treg tf order, so keep a separate map for that.
        tf_listener_idx = {tf: i for i, tf in enumerate(listener_tf_ids)}
        tf_treg_idx = {tf: i for i, tf in enumerate(treg.tf_ids)}
        mrna_idx = {m: i for i, m in enumerate(mrna_ids)}

        # Resolve user TU ids (accept with/without compartment tag):
        idx_by_full = {r: i for i, r in enumerate(rna_ids)}
        idx_by_bare = {}
        for i, r in enumerate(rna_ids):
            bare = r[:-3] if r.endswith(']') and '[' in r else r
            idx_by_bare.setdefault(bare, i)

        resolved_tus = []  # list of (tu_full_id, tu_index)
        for tu in TU_IDS:
            if tu in idx_by_full:
                resolved_tus.append((tu, idx_by_full[tu]))
            else:
                bare = tu[:-3] if tu.endswith(']') and '[' in tu else tu
                if bare in idx_by_bare:
                    i = idx_by_bare[bare]
                    resolved_tus.append((rna_ids[i], i))
                else:
                    print(f"NOTE: {tu} not found in RnaSynthProb "
                          f"rnaIds; skipping.")
        if not resolved_tus:
            print("WARNING: none of TU_IDS were valid, thus, there is nothing "
                  "to plot. Please check the list of TU_IDS and try again.")
            return

        def tf_label(tf_id):
            return TF_SHORT_NAMES.get(tf_id, tf_id)

        def fmt(x):
            return "0" if abs(x) < 1e-30 else f"{x:.2e}"

        # Build the time axis + generation boundaries:
        time = read_stacked_columns(cell_paths, 'Main', 'time').squeeze()
        t_min = time / 60.0
        # Determine the final time of each generation to create the internal
        # division boundaries:
        gen_end_times = np.atleast_1d(read_stacked_columns(
            cell_paths, 'Main', 'time', fun=lambda x: x[-1]).squeeze())
        gen_boundaries_min = (gen_end_times / 60.0)[:-1]  # drop the sim end

        # Deterine ppGpp regulation state
        # In wcEcoli ``ppgpp_regulation`` is always written to metadata by
        # runSim/fw_queue (it is in scriptBase.METADATA_KEYS):
        ppgpp_on = bool(metadata.get('ppgpp_regulation', True))
        ppgpp_title = "ppgpp reglation ON" if ppgpp_on else "ppgpp regulation OFF"

        # These listener columns are very wide (n_bound_TF_per_TU is n_TU*n_TF ≈
        # 75k columns), so subselect the needed TU (and TU×TF) columns per
        # generation via the `fun` hook before stacking (reading them in full
        # across all generations would be many GB):
        tu_indices = [ti for _, ti in resolved_tus]
        tu_pos = {tu_full: p for p, (tu_full, _) in enumerate(resolved_tus)}

        target_all = read_stacked_columns(
            cell_paths, 'RnaSynthProb', 'target_rna_synth_prob',
            fun=lambda x: x[:, tu_indices])
        actual_all = read_stacked_columns(
            cell_paths, 'RnaSynthProb', 'actual_rna_synth_prob',
            fun=lambda x: x[:, tu_indices])
        max_p = read_stacked_columns(
            cell_paths, 'RnaSynthProb', 'max_p').squeeze()
        pcn_all = read_stacked_columns(
            cell_paths, 'RnaSynthProb', 'promoter_copy_number',
            fun=lambda x: x[:, tu_indices])
        toc_all = read_stacked_columns(
            cell_paths, 'RnaSynthProb', 'tu_is_overcrowded',
            fun=lambda x: x[:, tu_indices])
        # n_bound flat columns for every resolved TU × every listener TF:
        nbound_cols = [ti * n_tf + j for ti in tu_indices for j in range(n_tf)]
        nbound_all = read_stacked_columns(
            cell_paths, 'RnaSynthProb', 'n_bound_TF_per_TU',
            fun=lambda x: x[:, nbound_cols])

        # RNA counts: only the resolved TUs that are mRNAs (have a count column):
        # TODO (mia): for tRNA and rRNA, make this plot the mature vs immature
        #  counts, or the # in ribosomes etc.
        count_cols = []
        count_pos = {}
        for tu_full, _ in resolved_tus:
            if tu_full in mrna_idx:
                count_pos[tu_full] = len(count_cols)
                count_cols.append(mrna_idx[tu_full])
        if count_cols:
            total_all = read_stacked_columns(
                cell_paths, 'RNACounts', 'mRNA_counts',
                fun=lambda x: x[:, count_cols])
            full_all = read_stacked_columns(
                cell_paths, 'RNACounts', 'full_mRNA_counts',
                fun=lambda x: x[:, count_cols])
            partial_all = read_stacked_columns(
                cell_paths, 'RNACounts', 'partial_mRNA_counts',
                fun=lambda x: x[:, count_cols])

        exp_id = metadata.get('description', '?')

        # Create one figure per TU inputted:
        for tu_full, tu_idx in resolved_tus:
            label = tu_full[:-3] if tu_full.endswith(']') else tu_full
            p = tu_pos[tu_full]  # position in the subselected arrays

            target = target_all[:, p]
            actual = actual_all[:, p]
            copies = pcn_all[:, p].astype(float)
            b0 = float(basal_prob[tu_idx])

            hl = tu_half_life_min.get(tu_full, float("nan"))
            hl_txt = f"{hl:.1f} min" if np.isfinite(hl) else "n/a"

            # TFs regulating this TU that are present in the listener TF order:
            regulating_tfs = [
                tf for tf in rel.rna_id_to_regulating_tfs.get(tu_full, [])
                if tf in tf_listener_idx]
            tf_color = {
                tf: TF_PALETTE[i % len(TF_PALETTE)]
                for i, tf in enumerate(regulating_tfs)}

            def delta(tf):
                if tf in tf_treg_idx:
                    return float(dpm[tu_idx, tf_treg_idx[tf]])
                return 0.0

            # Compute n_bound per TF for this TU. In the subselected array the
            # block for this TU starts at column p * n_tf (the TF offset is its
            # listener idx):
            nbound = {
                tf: nbound_all[:, p * n_tf + tf_listener_idx[tf]].astype(float)
                for tf in regulating_tfs}
            have_binding = len(regulating_tfs) > 0

            # obtain transcript counts (only if this TU is an mRNA):
            # TODO (mia): expand to tRNA and rRNA here too
            have_counts = tu_full in count_pos
            if have_counts:
                ci = count_pos[tu_full]
                total_c = total_all[:, ci]
                full_c = full_all[:, ci]
                partial_c = partial_all[:, ci]

            fig, axes = plt.subplots(
                5, 1, figsize=(12, 17), sharex=True,
                gridspec_kw={"height_ratios": [0.25, 0.20, 0.20, 0.14, 0.09]})
            ax1, ax2, ax3, ax4, ax5 = axes

            # Create generation (division) boundary lines on all panels:
            for gb in gen_boundaries_min:
                for ax in axes:
                    ax.axvline(gb, color="0.6", lw=1, ls=":", alpha=0.6)

            # panel 1: synth prob (target / actual)
            ax1.plot(t_min, target, color="0.2", lw=.4,
                     label="target synth prob")
            ax1.plot(t_min, actual, color="#1f77b4", lw=1, alpha=0.2,
                     label="actual synth prob")
            ax1.set_ylabel("synth prob")
            ax1.set_title(
                f"RNA synthesis probability  "
                f"(mean target={np.nanmean(target):.2e}, "
                f"mean actual={np.nanmean(actual):.2e})",
                loc="left")
            ax1.set_ylim(bottom=0)
            ax1.legend(loc="upper right", fontsize=8, framealpha=0.9)

            # panel 2: promoters bound per TF (stacked bars)
            if have_binding:
                dt = float(np.median(np.diff(t_min))) if t_min.size > 1 else 1.0
                bottom = np.zeros_like(t_min, dtype=float)
                for tf in regulating_tfs:
                    vals = nbound[tf]
                    ax2.bar(
                        t_min, vals, bottom=bottom, width=dt, align="edge",
                        color=tf_color[tf], linewidth=0,
                        label=f"{tf_label(tf)} (default delta_prob={fmt(delta(tf))})")
                    bottom += vals
                # Max promoter regions available for this TU at each timestep:
                ax2.plot(t_min, copies, color="black", lw=1.2, drawstyle=STEP,
                         label="promoter copies available (max)")
                ax2.legend(loc="upper right", fontsize=8, framealpha=0.9)
            ax2.set_ylabel("# promoters bound")
            ax2.set_title("Promoters bound per TF (stacked)", loc="left")
            ax2.set_ylim(bottom=0)

            # panel 3: transcript counts
            if have_counts:
                for series, color, ls, tag in [
                    (total_c, "#2c782c", "-", "total"),
                    (full_c, "#1f77b4", "--", "full"),
                    (partial_c, "#d62728", ":", "partial"),
                ]:
                    ax3.plot(t_min, series, color=color, lw=1.8, ls=ls,
                             label=f"{tag} mRNA")
                ax3.legend(loc="upper right", fontsize=8, framealpha=0.9)
                counts_txt = (
                    f", mean total={np.nanmean(total_c):.2f}"
                    f", full={np.nanmean(full_c):.2f}")
            else:
                counts_txt = " — not an mRNA (no counts)"
            ax3.set_ylabel("mRNA count")
            ax3.set_title(
                f"Transcript counts (mRNA half-life ≈ {hl_txt}{counts_txt})",
                loc="left")
            ax3.set_ylim(bottom=0)

            # panel 4: max_p
            ax4.plot(t_min, max_p, color="0.3", lw=1.5, label="max_p")
            ax4.set_ylabel("max_p")
            ax4.set_title("max_p (RNAP initiation-probability cap)", loc="left")
            ax4.set_ylim(bottom=0)
            ax4.legend(loc="upper right", fontsize=8, framealpha=0.9)

            # panel 5: RNAP overcrowding
            # tu_is_overcrowded = the TU plotted actually hit the max_p cap
            # actual != target only tells you the cap fired somewhere in
            # the cell, because the rescale multiplies every non-capped promoter.
            cellwide = ~np.isclose(actual, target, rtol=1e-6, atol=0.0)
            ax5.fill_between(
                t_min, 0, cellwide.astype(int), step="post", color="0.6",
                alpha=0.55, linewidth=0,
                label=f"cell-wide renorm fired ({100 * cellwide.mean():.1f}%)")
            this_tu = toc_all[:, p].astype(bool)
            ax5.fill_between(
                t_min, 0, this_tu.astype(int), step="post", color="red",
                alpha=0.85, linewidth=0,
                label=f"THIS TU capped ({100 * this_tu.mean():.1f}%)")
            ax5.set_ylabel("overcrowding")
            ax5.set_ylim(-0.05, 1.05)
            ax5.set_yticks([0, 1])
            ax5.set_yticklabels(["no", "yes"])
            ax5.set_title(
                f"RNAP overcrowding | this TU capped {this_tu.sum()}/{this_tu.size} "
                f"timesteps ({100 * this_tu.mean():.1f}%) | "
                f"cell-wide renorm {100 * cellwide.mean():.1f}%",
                loc="left")
            ax5.legend(loc="upper right", fontsize=8, framealpha=0.9)
            ax5.set_xlabel("Time (min)")

            # ParCa-fit reference numbers (basal_prob + each TF's default delta):
            tf_delta_bits = ", ".join(
                f"{tf_label(tf)} (default delta_prob={fmt(delta(tf))})"
                for tf in regulating_tfs) or "none"
            ppgpp_note = (
                "\nppGpp active → basal_prob & TF deltas rescaled each timestep"
                if ppgpp_on else "ppgpp regulation inactive")
            fig.suptitle(
                f"{label}\n{exp_id} · {len(gen_end_times)} generation(s) · "
                f"{ppgpp_title}\n"
                f"ParCa-fit basal_prob={fmt(b0)}; {tf_delta_bits}{ppgpp_note}",
                fontsize=11)
            fig.tight_layout(rect=(0, 0, 1, 0.94))

            exportFigure(plt, plotOutDir, f"{plotOutFileName}_{label}", metadata)
            plt.close(fig)
            print(f"Plot successful: saved figure for {tu_full}")


if __name__ == '__main__':
    Plot().cli()
