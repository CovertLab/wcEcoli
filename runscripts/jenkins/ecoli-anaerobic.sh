set -e

source runscripts/jenkins/setup-environment.sh
sh runscripts/jenkins/fireworks-config.sh anaerobic

echo y | lpad reset

DESC="Anaerobic." VARIANT="condition" FIRST_VARIANT_INDEX=4 \
  LAST_VARIANT_INDEX=4 SINGLE_DAUGHTERS=1 N_GENS=8 MASS_DISTRIBUTION=0 \
  TRNA_ATTENUATION=0 PPGPP_REGULATION=0 MECHANISTIC_TRANSLATION_SUPPLY=0 \
  MECHANISTIC_AA_TRANSPORT=0 AA_SUPPLY_IN_CHARGING=0 D_PERIOD_DIVISION=0 \
  MECHANISTIC_REPLISOME=1 TIMESTEP_MAX=2 COMPRESS_OUTPUT=1 PLOTS=ACTIVE \
  RAISE_ON_TIME_LIMIT=1 python runscripts/fireworks/fw_queue.py

bash runscripts/jenkins/run-fireworks.sh

runscripts/jenkins/save_output.sh out/ /scratch/groups/mcovert/wc_ecoli/anaerobic/
