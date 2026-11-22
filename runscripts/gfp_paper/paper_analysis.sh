#!/usr/bin/env bash
# Paper analyses for the GFP burden paper.
#
# Runs on the mcovert partition, which has 32 cores with 8 held by Jenkins.
# Ask before using more than the free 24. Analyses may run while a simulation
# batch is going, because the sims hold owners and these hold mcovert.
#
# Every panel lands in <batch>/plotOut/paper_figures/ as PNG, PDF and SVG with
# a JSON provenance sidecar, so one download gets the whole set.
#
# Panels are listed in panels.tsv; the single-batch ones run through the
# GFP_PAPER analysis tags, the multi-batch ones (4F, 4H) through their own
# scripts in this directory.

set -euo pipefail

# Usage:
#   paper_analysis.sh batch <batch_dir> [cpus]
#       every single-batch panel (variant and multigen GFP_PAPER tags);
#       panels that do not apply to the batch (for example 4E on a GFP
#       batch) print why and skip
#   paper_analysis.sh fig4f <batch_dir> [<batch_dir> ...]
#       4F, GFP with and without a kcat bound, across one or more batches
#   paper_analysis.sh fig4h <kcat_batch_dir> <gfp_batch_dir>
#       4H, ppGpp vs growth, kcat tightening (M1) vs GFP (B1)
#
# Run from the repository root with PYTHONPATH set to it. The first panel on
# a batch builds the shared cache in <batch>/plotOut/paper_figures/_cache/;
# later panels reuse it.

usage() {
	sed -n '/^# Usage:/,/^# later panels/p' "$0"
	exit 1
}

[ $# -ge 2 ] || usage
MODE="$1"
shift
export OPENBLAS_NUM_THREADS=1

case "${MODE}" in
	batch)
		BATCH="$1"
		CPUS="${2:-24}"
		python runscripts/manual/analysisVariant.py "${BATCH}" \
			-p GFP_PAPER --cpus "${CPUS}"
		# Multigen panels (5A, 6A) are one file per lineage, and the runner
		# analyses one variant and seed per call, so run every lineage, in
		# parallel. Limit with MULTIGEN_VARIANTS="4 5" and/or
		# MULTIGEN_SEEDS="0 1 2"; by default every variant and seed found.
		VARIANTS="${MULTIGEN_VARIANTS:-$(ls -d "${BATCH}"/*_[0-9][0-9][0-9][0-9][0-9][0-9] \
			| sed -E 's/.*_0*([0-9]+)$/\1/' | sort -n | uniq | tr '\n' ' ')}"
		FIRST_VARIANT_DIR="$(ls -d "${BATCH}"/*_[0-9][0-9][0-9][0-9][0-9][0-9] | head -1)"
		SEEDS="${MULTIGEN_SEEDS:-$(ls -d "${FIRST_VARIANT_DIR}"/[0-9][0-9][0-9][0-9][0-9][0-9] \
			| sed -E 's/.*\/0*([0-9]+)$/\1/' | sort -n | tr '\n' ' ')}"
		echo "Multigen panels for variants: ${VARIANTS}; seeds: ${SEEDS}"
		for v in ${VARIANTS}; do for s in ${SEEDS}; do echo "$v $s"; done; done \
			| xargs -P "${CPUS}" -n 2 sh -c 'python runscripts/manual/analysisMultigen.py \
				"$0" -p GFP_PAPER -v "$1" -s "$2" > /dev/null 2>&1 \
				|| echo "multigen failed: variant $1 seed $2"' "${BATCH}"
		;;
	fig4f)
		python runscripts/gfp_paper/gfp_paper_gfp_kcat_constraints.py "$@" \
			--cpus "${CPUS:-24}"
		;;
	fig4h)
		[ $# -eq 2 ] || usage
		python runscripts/gfp_paper/gfp_paper_ppgpp_vs_growth_kcat_tightening_vs_gfp.py \
			--kcat-batch "$1" --gfp-batch "$2" --cpus "${CPUS:-24}"
		;;
	*)
		usage
		;;
esac
