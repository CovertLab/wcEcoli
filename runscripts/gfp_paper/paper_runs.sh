#!/usr/bin/env bash
# Simulation batches for the GFP burden paper.
#
# ONLY ONE BATCH RUNS AT A TIME. Run these in order, waiting for each to
# finish before starting the next, and record every batch in batches.tsv as it
# is submitted and again when it completes. See README.md for why the order is
# what it is.
#
# All batches share one sim_data: EGFP inserted at the lambda attB site,
# 24 generations, induction at generation 8, on CPU generations with AVX-512.
#
# STATUS: skeleton. The commands are filled in once the ParCa (2.1) and the
# smoke test (2.2) have passed, so that nothing here encodes a setup that has
# not been checked.

set -euo pipefail

NEW_GENES="egfp"
N_GENS=24
# Reproducibility: these CPU generations support AVX-512, so NumPy and
# OpenBLAS take the same code paths on every node. See vEcoli PR #444.
# Verified bitwise across all four, five CPU models from AMD and Intel
# (Sherlock job 1, 2026-11-19). Skylake-SP also has AVX-512 but was not tested,
# so it stays out.
#
# The constraint lives in the per-user ~/wcEcoli/my_qadapter.yaml, not in the
# lab-wide template (wholecell/fireworks/templates/my_qadapter.yaml), so it
# does not change anyone else's runs. A regenerated qadapter would drop it,
# so every batch checks for it before submitting anything.
CPU_CONSTRAINT="CPU_GEN:BGM|CPU_GEN:SIE|CPU_GEN:GEN|CPU_GEN:SPR"
QADAPTER="${QADAPTER:-${HOME}/wcEcoli/my_qadapter.yaml}"

require_cpu_pin() {
	if ! grep -Eq "^constraint:[[:space:]]*\"?${CPU_CONSTRAINT//|/\\|}\"?[[:space:]]*$" "${QADAPTER}"; then
		echo "STOP: ${QADAPTER} does not pin the CPU generations." >&2
		echo "Expected the line: constraint: \"${CPU_CONSTRAINT}\"" >&2
		exit 1
	fi
}
require_cpu_pin

echo "Not yet filled in. Blocked on tasks 2.1 (ParCa) and 2.2 (smoke test)."
echo "New genes: ${NEW_GENES}, generations: ${N_GENS}"
echo "CPU constraint: ${CPU_CONSTRAINT}"
exit 1

## 1. Smoke test -- 108 sims
#    Checks every listener the panels need is recorded, that seed-matched
#    pre-induction generations are identical across variants, that a rerun
#    reproduces bit for bit, and how far translation sits below the crowding
#    cap after the footprint fix. DONE 2026-11-19 (20261119.154823__gfp_smoke).

## 2. B1 scan -- 46 variants x 5 seeds x 24 gens = 5,520. Expression factors
#    [0, 7, 8, 8.33, 8.66, 9] and translation efficiencies
#    [20, 10, 7, 5, 3, 2, 1, 0.1, 0] in
#    models/ecoli/sim/variants/new_gene_internal_shift.py.

## 3. W wildtype -- match 20260302.164308__wildtype_minimal_media exactly
#    (~773 cells). Do NOT enlarge: the kcat estimator is a maximum, so more
#    cells raise the estimates without improving them. See the plan.

## 4. B2 focus -- 4 designs x 60 seeds x 24 gens = 5,760, using
#    models/ecoli/sim/variants/new_gene_design_list.py. The control must use
#    the same seed numbers as the induced designs.

## 5. B1deep -- extra seeds on the subset of B1 chosen from the scan.

## 6. M1 -- wildtype with forced perceived-enzyme reduction.

## 7. M2 -- GFP with and without kcat constraints, and with forced reduction.
