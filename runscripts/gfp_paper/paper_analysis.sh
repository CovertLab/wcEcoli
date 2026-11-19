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
# STATUS: skeleton. Filled in as the per-panel scripts listed in panels.tsv
# are written.

set -euo pipefail

if [ $# -lt 1 ]; then
	echo "usage: $0 <batch_dir>    # e.g. \$SCRATCH/wcEcoli_out/<batch>"
	exit 1
fi
BATCH="$1"

echo "Not yet filled in. See panels.tsv for the scripts still to be written."
echo "Batch: ${BATCH}"
exit 1

## Variant analyses:
#    python runscripts/manual/analysisVariant.py "${BATCH}" -p <script>
## Multigen analyses:
#    python runscripts/manual/analysisMultiGen.py "${BATCH}" -p <script>
