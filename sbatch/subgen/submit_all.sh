#!/usr/bin/bash
# Submit the full subgenerational pipeline with the correct dependencies.
#
#   extract      (heavy simOut pass; the classification analyses depend on it)
#     -> cached  (afterok:extract) -- canonical tables + figures
#          -> cli (afterany:cached) -- the standalone argparse CLIs
#   independent  (no dependency; reads simOut directly, runs alongside)
#
# Two notes on that DAG, both changed from the old four-file layout:
#
#   * `cli` depends only on `cached`, not on `independent`. Every one of its
#     inputs comes from `extract` and `cached`; nothing in `independent` feeds
#     it, so making it wait on `independent` only delayed it.
#   * `cli` uses afterany, not afterok. subgen_summary returns 1 if ANY single
#     script in a stage failed, so afterok would cancel `cli` outright over one
#     unrelated failure -- and `cli` already skips gracefully when its per-gene
#     table is missing.
#
# Usage:
#   bash submit_all.sh                          # Set 1 profile (see subgen_env.sh)
#   bash submit_all.sh set2                     # Set 2, all three variants
#   SUBGEN_OUT=... SUBGEN_VARIANTS="-v 0" bash submit_all.sh
#
# Configuration lives in subgen_env.sh and is passed to every job through the
# environment, so the sim directory and variant range are set ONCE.
#
# To re-run a single stage instead of the whole pipeline:
#   SUBGEN_STAGE=cli sbatch subgen_stage.sbatch

set -euo pipefail
cd "$(dirname "$0")"

case "${1:-}" in
	set1)
		export SUBGEN_OUT="$SCRATCH/wcEcoli_out/20260706.223256__subgen_paper_sims_1"
		export SUBGEN_VARIANTS="-v 0"
		;;
	set2)
		export SUBGEN_OUT="$SCRATCH/wcEcoli_out/20260804.115726__subgen_paper_sims_2"
		# 0 2 INCLUSIVE = basal, with_aa, acetate. The original Set 2 run used
		# `0 1`, silently excluding condition_000002 (acetate). Acetate's sims are
		# incomplete (1-13 of 32 generations for sampled seeds), so its outputs will
		# be thin -- include it anyway so the gap is recorded in run_metadata.json
		# instead of being invisible.
		export SUBGEN_VARIANTS="--variant-range 0 2"
		;;
	'')
		;;  # use whatever subgen_env.sh defaults to / the caller exported
	*)
		echo "Unknown profile '$1'. Use: set1 | set2 | (nothing)" >&2
		exit 1
		;;
esac

echo "Submitting with:"
echo "  SUBGEN_OUT      = ${SUBGEN_OUT:-<subgen_env.sh default>}"
echo "  SUBGEN_VARIANTS = ${SUBGEN_VARIANTS:-<subgen_env.sh default>}"

# Per-stage resources. One job script serves all four stages, so the request is
# made here on the sbatch command line, which overrides subgen_stage.sbatch's
# #SBATCH directives.
submit() {  # submit <stage> <time> <mem> [extra sbatch args...]
	local stage="$1" time="$2" mem="$3"; shift 3
	SUBGEN_STAGE="$stage" sbatch --parsable --export=ALL \
		--job-name="subgen_$stage" --time="$time" --mem="$mem" \
		"$@" subgen_stage.sbatch
}

jid_extract=$(submit extract 12:00:00 128GB)
echo "extract     -> job $jid_extract (no dependency)"

jid_independent=$(submit independent 8:00:00 128GB)
echo "independent -> job $jid_independent (no dependency)"

jid_cached=$(submit cached 4:00:00 32GB \
	--dependency=afterok:"$jid_extract")
echo "cached      -> job $jid_cached (afterok:$jid_extract)"

jid_cli=$(submit cli 8:00:00 64GB \
	--dependency=afterany:"$jid_cached")
echo "cli         -> job $jid_cli (afterany:$jid_cached)"
