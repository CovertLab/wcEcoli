#! /usr/bin/env bash
set -e

source runscripts/jenkins/setup-environment.sh
echo y | lpad reset

DESC="Modeling tab files run" VARIANT="condition" FIRST_VARIANT_INDEX=0 \
LAST_VARIANT_INDEX=3 SINGLE_DAUGHTERS=1 N_GENS=4 N_INIT_SIMS=64 \
RUN_AGGREGATE_ANALYSIS=0 EXPORT_ECOCYC_FILES=1 \
python runscripts/fireworks/fw_queue.py
