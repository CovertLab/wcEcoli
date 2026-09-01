#!/usr/bin/bash
# Shared configuration for the subgenerational analysis pipeline.
#
# Sourced by subgen_stage.sbatch, so the sim directory and variant selection are
# set in ONE place instead of being hand-edited per stage (which is how Set 1 and
# Set 2 ended up analyzed with different variant ranges).
#
# Override per submission by exporting before you submit, or by using
# submit_all.sh, which does it for you:
#
#   SUBGEN_OUT=$SCRATCH/wcEcoli_out/20260804.115726__subgen_paper_sims_2 \
#   SUBGEN_VARIANTS="--variant-range 0 2" \
#   bash submit_all.sh
#
# Variant selection notes:
#   -v 0                -> single variant, index 0 (e.g. wildtype_000000). Use for
#                          sims with no variant applied.
#   --variant-range 0 2 -> variants 0,1,2 INCLUSIVE
#                          (wholecell/utils/scriptBase.py: range(start, end + 1)).
#                          Set 2 has three condition variants, so it needs `0 2`;
#                          it was originally run with `0 1`, which silently
#                          excluded condition_000002 (acetate).

set -euo pipefail

REPO="${SUBGEN_REPO:-$HOME/repos/wcEcoli}"

# Default profile: Set 1 (single wildtype variant).
OUT="${SUBGEN_OUT:-$SCRATCH/wcEcoli_out/20260706.223256__subgen_paper_sims_1}"
VARIANTS="${SUBGEN_VARIANTS:--v 0}"

# The subgen suite is cohort-only: both multigen Figure-5 scripts were retired to
# models/ecoli/analysis/multigen/old_subgen_scripts/, so there is no $MULTIGEN and
# no -s/--seed output-directory selection to configure any more.
COH="$REPO/models/ecoli/analysis/cohort"

cd "$REPO"
export PYTHONPATH="$PWD"
export OPENBLAS_NUM_THREADS=1

echo "=== subgen pipeline configuration ==="
echo "  REPO     = $REPO"
echo "  OUT      = $OUT"
echo "  VARIANTS = $VARIANTS"
echo "  STAGE    = ${SUBGEN_STAGE:-<unset>}"
echo "  git      = $(git -C "$REPO" rev-parse --short HEAD) \
$(git -C "$REPO" rev-parse --abbrev-ref HEAD)\
$([ -n "$(git -C "$REPO" status --porcelain)" ] && echo ' (DIRTY)' || echo ' (clean)')"
echo "====================================="

# Every run_metadata.json records git_dirty. A dirty tree means the outputs are not
# reproducible, which is exactly why Set 1 and Set 2 cannot be compared today.
if [ -n "$(git -C "$REPO" status --porcelain)" ]; then
	echo "WARNING: the working tree is DIRTY. Outputs will be tagged git_dirty=true"
	echo "         and will not be reproducible. Commit before a publication run."
fi

# Run one analysis script, reporting its exit status without aborting the job.
# `set -e` would stop the whole loop on the first failure; these scripts are
# independent, so a single failure should not cost the rest of the job. The summary
# at the end makes failures impossible to miss (the old loops printed a traceback
# and carried on silently).
SUBGEN_FAILED=()
run_analysis() {
	local script="$1"; shift
	local status=0
	echo "=== Running $(basename "$script") ==="
	# Capture the status BEFORE echoing: inside `if python ...; then`, `$?` is the
	# status of the `if` construct, not of python, so the old code always
	# reported "exit 0" on a failure.
	python "$script" "$@" || status=$?
	if [ "$status" -eq 0 ]; then
		echo "--- OK: $(basename "$script")"
	else
		echo "--- FAILED (exit $status): $(basename "$script")"
		SUBGEN_FAILED+=("$(basename "$script")")
	fi
}

subgen_summary() {
	if [ ${#SUBGEN_FAILED[@]} -eq 0 ]; then
		echo "=== All scripts completed successfully. ==="
	else
		echo "=== ${#SUBGEN_FAILED[@]} script(s) FAILED: ${SUBGEN_FAILED[*]} ==="
		return 1
	fi
}
