set -e

source runscripts/jenkins/setup-environment.sh
sh runscripts/jenkins/fireworks-config.sh optional

echo y | lpad reset

DESC="No growth rate control" TRNA_ATTENUATION=0 PPGPP_REGULATION=0 \
  MECHANISTIC_TRANSLATION_SUPPLY=0 MECHANISTIC_AA_TRANSPORT=0 \
  AA_SUPPLY_IN_CHARGING=0 D_PERIOD_DIVISION=0 MECHANISTIC_REPLISOME=1 \
  N_GENS=4 PLOTS=ACTIVE WC_ANALYZE_FAST=1 \
  python runscripts/fireworks/fw_queue.py
DESC="No operons" OPERONS="off" N_GENS=4 PLOTS=ACTIVE WC_ANALYZE_FAST=1 \
  python runscripts/fireworks/fw_queue.py
DESC="Superhelical Densities" SUPERHELICAL_DENSITIES=1 N_GENS=4 \
  PARALLEL_PARCA=1 SINGLE_DAUGHTERS=1 COMPRESS_OUTPUT=1 RAISE_ON_TIME_LIMIT=1 \
  PLOTS=ACTIVE WC_ANALYZE_FAST=1 \
  python runscripts/fireworks/fw_queue.py
DESC="Causality Network" BUILD_CAUSALITY_NETWORK=1 N_GENS=2 SEED=$RANDOM \
  PARALLEL_PARCA=1 SINGLE_DAUGHTERS=1 COMPRESS_OUTPUT=1 RAISE_ON_TIME_LIMIT=1 \
  WC_ANALYZE_FAST=1 \
  python runscripts/fireworks/fw_queue.py
DESC="New Gene Existence" N_GENS=4 NEW_GENES="gfp" PLOTS=ACTIVE \
  COMPRESS_OUTPUT=1 RAISE_ON_TIME_LIMIT=1 WC_ANALYZE_FAST=1 \
  python runscripts/fireworks/fw_queue.py
# NOTE: The version of global variables from new_gene_internal_shift.py on the
# master branch will be used for this Jenkins test.
# Variant 2 is intended to correspond to expression factor 7, translation
# efficiency value 1.
DESC="New Gene Expression and Translation Efficiency Internal Shift" \
  VARIANT="new_gene_internal_shift" \
  FIRST_VARIANT_INDEX=2 LAST_VARIANT_INDEX=2 N_GENS=4 NEW_GENES="gfp" \
  PLOTS=ACTIVE COMPRESS_OUTPUT=1 RAISE_ON_TIME_LIMIT=1 WC_ANALYZE_FAST=1 \
  python runscripts/fireworks/fw_queue.py

bash runscripts/jenkins/run-fireworks.sh

runscripts/jenkins/save_output.sh out/ /scratch/groups/mcovert/wc_ecoli/optional_features/
