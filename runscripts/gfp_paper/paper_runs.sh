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
# OpenBLAS take the same code paths on every node. See vEcoli PR #444, and
# the matching constraint in wholecell/fireworks/templates/my_qadapter.yaml.
CPU_CONSTRAINT="CPU_GEN:BGM|CPU_GEN:SIE|CPU_GEN:GEN|CPU_GEN:SPR"

echo "Not yet filled in. Blocked on tasks 2.1 (ParCa) and 2.2 (smoke test)."
echo "New genes: ${NEW_GENES}, generations: ${N_GENS}"
echo "CPU constraint: ${CPU_CONSTRAINT}"
exit 1

## 1. Smoke test -- 108 sims
#    Checks every listener the panels need is recorded, that seed-matched
#    pre-induction generations are identical across variants, that a rerun
#    reproduces bit for bit, and that crowding flags fall after the footprint
#    fix.

## 2. B1 scan -- 49 variants x 5 seeds x 24 gens = 5,880, plus 5 extra control
#    seeds = 6,000. Expression factors [0, 7, 8, 8.33, 8.66, 9, 10] and
#    translation efficiencies [10, 7, 5, 3, 2, 1, 0.1, 0] in
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
