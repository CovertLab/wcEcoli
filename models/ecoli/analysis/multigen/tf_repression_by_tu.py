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
  3. Transcript/molecule counts, branched by RNA type:
       - mRNA: total/full/partial transcript counts (per-TU). The title
         reports the TU's mRNA half-life as well.
       - rRNA: nascent partial (this TU), free mature rRNA (all rRNA), rRNA
         held in ribosomes (all rRNA). The mature and in-ribosome series are
         cell-wide pools because RnaMaturation consolidates all rRNA variants
         into one shared 16S/23S/5S species, so the counts aren't attributable
         to a single operon/TU (see more info on this below).
       - tRNA: unprocessed precursor (this TU), mature uncharged, mature
         charged, summed over the tRNA cistrons this TU encodes. The mature
         uncharged/charged pools are cell-wide because each matured tRNA is a
         single bulk-molecule pool shared across the whole cell (and, for some
         tRNAs, produced by more than one TU), so the count isn't attributable
         to this operon alone (see more info on this below).
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

Note on stable-RNA (tRNA/rRNA) counts in panel 3:
Unlike mRNAs (unique molecules with genuine per-TU transcript counts) matured
tRNAs/rRNAs are bulk molecules held in shared, cell-wide pools. rRNA variants
are consolidated by the RnaMaturation process into a single 16S/23S/5S species
used by ribosomal subunits, so an rRNA TU's "mature"/"in-ribosome" counts are
whole-cell pools, not attributable to that operon. For tRNA, the mature molecule
keeps its own cistron id (no consolidation), so most tRNA species trace to a
single TU, but the count is still the cell-wide pool of that species. Only the
nascent (partial) and unprocessed-precursor stages are truly per-TU. Panel 3
labels each series accordingly; the mature/in-ribosome means are on a much
larger scale than the nascent/precursor series, so a twin axis is used.

# TODOs:
- random, but figure out if it is ok that partial tRNA do not appear to be
counted anywhere while rRNA do?
- consider changing the axis sharing methodology in panel 3 for rRNAs and tRNAs
"""

import os
import pickle

from matplotlib import pyplot as plt
import numpy as np

from models.ecoli.analysis import multigenAnalysisPlot
from wholecell.analysis.analysis_tools import (
    exportFigure, read_stacked_columns, read_stacked_bulk_molecules)
from wholecell.io.tablereader import TableReader
from wholecell.utils import units

# USER INPUTS
# List the TU(s) to be plotted (with or without compartment tag):
TU_IDS = ["TU103[c]", "TU00415[c]", "TU0-1181[c]", "TU0-13035[c]", "TU00507",
          "TU0-13001[c]", "TU0-1182[c]"]

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

        def favg(series):
            # mean of a time series, formatted with magnitude-appropriate
            # precision, for the panel-3 legend labels
            m = float(np.nanmean(series))
            if abs(m) >= 100:
                return f"{m:.0f}"
            if abs(m) >= 1:
                return f"{m:.1f}"
            return f"{m:.2f}"

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

        # mRNA transcript counts: only the resolved TUs that are mRNAs (have a
        # count column). tRNA/rRNA counts are handled separately below.
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

        # Stable-RNA (tRNA/rRNA) count sources for panel 3
        # RnaSynthProb rnaIds are aligned with transcription.rna_data, so a TU's
        # index (tu_idx) indexes the type flags directly. Unlike mRNAs (unique
        # molecules with per-TU transcript counts), matured tRNA/rRNA live as
        # bulk molecules in shared, cell-wide pools (rRNA variants are even
        # consolidated into a single 16S/23S/5S species), so their "mature" and
        # "in-ribosome" counts are not attributable to a single TU. Only the
        # nascent/unprocessed-precursor stages remain per-TU.
        rna_data = trans.rna_data
        cistron_ids = list(trans.cistron_data['id'])
        cistron_is_tRNA = np.asarray(trans.cistron_data['is_tRNA'])
        ctu = trans.cistron_tu_mapping_matrix  # (n_cistron, n_TU)
        uncharged_names = list(trans.uncharged_trna_names)
        charged_names = list(trans.charged_trna_names)
        uncharged_pos = {n: i for i, n in enumerate(uncharged_names)}

        # rRNA mature/ribosome species (cell-wide, consolidated pools):
        mg = sim_data.molecule_groups
        mi = sim_data.molecule_ids
        rrna_16s_ids = list(mg.s30_16s_rRNA)
        rrna_23s_ids = list(mg.s50_23s_rRNA)
        rrna_5s_ids = list(mg.s50_5s_rRNA)
        s30_id = mi.s30_full_complex
        s50_id = mi.s50_full_complex

        def trna_species_for_tu(ti):
            # mature uncharged/charged bulk ids for the tRNA cistrons in this TU
            col = np.asarray(ctu[:, ti].todense()).ravel()
            unids, chids = [], []
            for c in np.nonzero(col)[0]:
                if not cistron_is_tRNA[c]:
                    continue
                key = f"{cistron_ids[c]}[c]"
                if key in uncharged_pos:
                    unids.append(key)
                    chids.append(charged_names[uncharged_pos[key]])
            return unids, chids

        # Gather every bulk molecule id needed across the resolved TUs, so they
        # can be read in a single pass and indexed per TU below:
        trna_species = {}   # tu_full -> (uncharged_ids, charged_ids)
        precursor_id = {}   # tu_full -> unprocessed-precursor bulk id
        any_rRNA = False
        bulk_ids_needed = set()
        for tu_full, ti in resolved_tus:
            if rna_data['is_unprocessed'][ti]:
                precursor_id[tu_full] = tu_full
                bulk_ids_needed.add(tu_full)
            if rna_data['is_tRNA'][ti]:
                unids, chids = trna_species_for_tu(ti)
                trna_species[tu_full] = (unids, chids)
                bulk_ids_needed.update(unids)
                bulk_ids_needed.update(chids)
            if rna_data['is_rRNA'][ti]:
                any_rRNA = True
        if any_rRNA:
            bulk_ids_needed.update(rrna_16s_ids)
            bulk_ids_needed.update(rrna_23s_ids)
            bulk_ids_needed.update(rrna_5s_ids)
            bulk_ids_needed.update([s30_id, s50_id])

        bulk_ids_needed = sorted(bulk_ids_needed)
        bulk_col = {name: i for i, name in enumerate(bulk_ids_needed)}
        if bulk_ids_needed:
            bulk_counts = read_stacked_bulk_molecules(
                cell_paths, (bulk_ids_needed,))[0].astype(float)
            if bulk_counts.ndim == 1:
                bulk_counts = bulk_counts[:, np.newaxis]
        else:
            bulk_counts = np.zeros((t_min.size, 0))

        def bulk_sum(ids):
            if not ids:
                return np.zeros_like(t_min, dtype=float)
            return bulk_counts[:, [bulk_col[i] for i in ids]].sum(axis=1)

        # rRNA nascent (partial) per-TU counts + active-ribosome counts:
        if any_rRNA:
            rrna_ids_attr = list(rc.readAttribute('rRNA_ids'))
            rrna_partial_all = read_stacked_columns(
                cell_paths, 'RNACounts', 'partial_rRNA_counts').astype(float)
            rrna_partial_col = {r: i for i, r in enumerate(rrna_ids_attr)}
            umc = TableReader(os.path.join(sim_out_dir, 'UniqueMoleculeCounts'))
            unique_ids = list(umc.readAttribute('uniqueMoleculeIds'))
            active_rib_idx = unique_ids.index('active_ribosome')
            active_ribosome = read_stacked_columns(
                cell_paths, 'UniqueMoleculeCounts', 'uniqueMoleculeCounts',
                fun=lambda x: x[:, [active_rib_idx]]).squeeze().astype(float)

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

            # obtain mRNA transcript counts (only if this TU is an mRNA):
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

            # Calculate avg. absolute deviation of actual from target over the sim
            # (how much the multinomial capping/renorm pushed this TU off target):
            mean_abs_dev = np.nanmean(np.abs(actual - target))
            ax1.set_title(
                f"RNA synthesis probability  "
                f"(mean target={np.nanmean(target):.2e}, "
                f"mean actual={np.nanmean(actual):.2e}, "
                f"mean |actual−target|={mean_abs_dev:.2e})",
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
            # Per-TF occupancy over the whole sim, matching the scatter hover
            # stat (see
            # models/ecoli/analysis/comparison/rna_synth_prob_comparison_plotly.py):
            # the ratio of the two time-means (mean n_bound over mean
            # promoter copies), i.e. the fraction of available promoter-copy
            # opportunities this TF occupied, shown as "bound/copies (pct%)".
            #
            if have_binding:
                copies_mean = float(np.mean(copies))
                occ_bits = []
                for tf in regulating_tfs:
                    bound_mean = float(np.mean(nbound[tf]))
                    if copies_mean > 0:
                        occ = (f"{tf_label(tf)} {bound_mean:.2f}/{copies_mean:.2f} "
                               f"({100.0 * bound_mean / copies_mean:.0f}%)")
                    else:
                        occ = f"{tf_label(tf)} {bound_mean:.2f}/N/A"
                    occ_bits.append(occ)
                occ_txt = "; ".join(occ_bits)
            else:
                occ_txt = "no regulating TFs"
            ax2.set_title(
                f"Promoters bound per TF | occupancy "
                f"(avg. bound/avg. available promoter copies): {occ_txt}",
                loc="left")
            ax2.set_ylim(bottom=0)

            # panel 3: transcript/molecule counts (structure based on RNA type)
            # mRNA is per-TU, tRNA/rRNA mature & in-ribosome/charged tRNA
            # series are shared cell-wide pools (see docstring), with only the
            # nascent/unprocessed-precursor stage attributable to this TU.
            is_rrna = bool(rna_data['is_rRNA'][tu_idx])
            is_trna = bool(rna_data['is_tRNA'][tu_idx])
            # For tRNA/rRNA the counts span ~3 orders of magnitude, so the
            # panel uses a right-hand twin axis (ax3_twin). Series are grouped by
            # magnitude, not scope, so every line is visible: the dominant
            # reservoir (rRNA-in-ribosomes / charged tRNA) sits alone on the
            # right axis, and the smaller pools share the left axis.
            ax3_twin = None
            if have_counts:
                # mRNA: per-TU total / full / partial transcript counts
                for series, color, ls, tag in [
                    (total_c, "#2c782c", "-", "total"),
                    (full_c, "#1f77b4", "--", "full"),
                    (partial_c, "#d62728", ":", "partial"),
                ]:
                    ax3.plot(t_min, series, color=color, lw=1.8, ls=ls,
                             label=f"{tag} mRNA")
                counts_txt = (
                    f", mean total={np.nanmean(total_c):.2f}"
                    f", full={np.nanmean(full_c):.2f}")
                ax3.set_ylabel("mRNA count")
                ax3.set_title(
                    f"Transcript counts (mRNA half-life ≈ {hl_txt}{counts_txt})",
                    loc="left")
            elif is_rrna:
                # rRNA: nascent (this TU) + free mature (cell-wide) share the
                # left axis, and the dominant in-ribosome pool (cell-wide) is
                # alone on the right twin axis:
                nascent = (rrna_partial_all[:, rrna_partial_col[tu_full]]
                           if tu_full in rrna_partial_col
                           else np.zeros_like(t_min, dtype=float))
                free_mature = (bulk_sum(rrna_16s_ids) + bulk_sum(rrna_23s_ids)
                               + bulk_sum(rrna_5s_ids))
                # rRNA molecules held in ribosomal machinery: a 30S carries one
                # rRNA (16S), a 50S two (23S+5S), and an active 70S all three:
                in_ribosomes = (bulk_sum([s30_id]) + 2.0 * bulk_sum([s50_id])
                                + 3.0 * active_ribosome)
                for series, color, ls, tag in [
                    (nascent, "#d62728", ":",
                     f"nascent partial (this TU, avg. = {favg(nascent)})"),
                    (free_mature, "#1f77b4", "--",
                     f"free mature rRNA (cell-wide, avg. = {favg(free_mature)})"),
                ]:
                    ax3.plot(t_min, series, color=color, lw=1.8, ls=ls, label=tag)
                ax3.set_ylabel("count: nascent + free mature")
                ax3_twin = ax3.twinx()
                ax3_twin.plot(
                    t_min, in_ribosomes, color="#2c782c", lw=1.8, ls="-",
                    label=f"rRNA in ribosomes (cell-wide, avg. = {favg(in_ribosomes)})")
                ax3_twin.set_ylabel("count: in ribosomes\n(cell-wide)")
                ax3.set_title("rRNA counts", loc="left")
            elif is_trna:
                # tRNA: unprocessed precursor (this TU) + mature uncharged
                # (cell-wide) share the left axis; the dominant charged pool
                # (cell-wide) is alone on the right twin axis. Both mature pools
                # are summed over this TU's tRNA cistrons:
                unids, chids = trna_species.get(tu_full, ([], []))
                uncharged = bulk_sum(unids)
                charged = bulk_sum(chids)
                precursor = (bulk_sum([precursor_id[tu_full]])
                             if tu_full in precursor_id else None)
                if precursor is not None:
                    ax3.plot(
                        t_min, precursor, color="#d62728", lw=1.8, ls=":",
                        label=f"unprocessed precursor (this TU, avg. = {favg(precursor)})")
                ax3.plot(
                    t_min, uncharged, color="#1f77b4", lw=1.8, ls="--",
                    label=f"mature uncharged (cell-wide, avg. = {favg(uncharged)})")
                ax3.set_ylabel(
                    "count: uncharged"
                    + ("" if precursor is None else " + precursor"))
                ax3_twin = ax3.twinx()
                ax3_twin.plot(
                    t_min, charged, color="#ff7f0e", lw=1.8, ls="-",
                    label=f"mature charged (cell-wide, avg. = {favg(charged)})")
                ax3_twin.set_ylabel("count: charged\n(cell-wide)")
                ax3.set_title("tRNA counts", loc="left")
            else:
                ax3.set_ylabel("count")
                ax3.set_title(
                    "Transcript counts — not an mRNA/tRNA/rRNA (no counts)",
                    loc="left")
            # Combined legend (fold in the twin axis handles when it is in use):
            if have_counts or is_rrna or is_trna:
                handles, labels = ax3.get_legend_handles_labels()
                if ax3_twin is not None:
                    h2, l2 = ax3_twin.get_legend_handles_labels()
                    handles += h2
                    labels += l2
                ax3.legend(handles, labels, loc="upper right", fontsize=8,
                           framealpha=0.9)
            ax3.set_ylim(bottom=0)
            if ax3_twin is not None:
                ax3_twin.set_ylim(bottom=0)

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
                label=f"this TU capped ({100 * this_tu.mean():.1f}%)")
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
                for tf in regulating_tfs) or "No regulating TFs"
            ppgpp_note = (
                "\nppGpp active → basal_prob & TF deltas rescaled each timestep"
                if ppgpp_on else "ppgpp regulation inactive")
            fig.suptitle(
                f"{label}\n{exp_id} · {len(gen_end_times)} generation(s) · "
                f"{ppgpp_title}\n"
                f"ParCa-fit basal_prob={fmt(b0)}; {tf_delta_bits}{ppgpp_note}",
                fontsize=11, y=0.995)
            fig.tight_layout(rect=(0, 0, 1, 0.985))

            exportFigure(plt, plotOutDir, f"{plotOutFileName}_{label}", metadata)
            plt.close(fig)
            print(f"Plot successful: saved figure for {tu_full}")


if __name__ == '__main__':
    Plot().cli()
