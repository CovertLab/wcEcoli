# Script to run simulations used to generate data for figures in the paper
# Data will be uncompressed to allow running of new analysis scripts
# Requires setup outlined in requirements.txt and in README.md
# Suggestion: run one set at a time by commenting out all other fw_queue.py lines

lpad reset
make clean compile

## Set A - basal condition
# Used for figure 4
DESC="SET A 32 gens 8 seeds basal with growth noise and D period" VARIANT="wildtype" FIRST_VARIANT_INDEX=0 LAST_VARIANT_INDEX=0 SINGLE_DAUGHTERS=1 N_GENS=32 N_INIT_SIMS=8 MASS_DISTRIBUTION=1 GROWTH_RATE_NOISE=1 D_PERIOD_DIVISION=1 python runscripts/fw_queue.py

## Set B - nutrient shift from minimal to minimal + AA
# Used for figure 1
DESC="SET B 9 gens 8 seeds shift to plus AA without growth noise with D period" VARIANT="nutrientTimeSeries" FIRST_VARIANT_INDEX=2 LAST_VARIANT_INDEX=2 SINGLE_DAUGHTERS=1 N_GENS=9 N_INIT_SIMS=8 MASS_DISTRIBUTION=1 GROWTH_RATE_NOISE=0 D_PERIOD_DIVISION=1 python runscripts/fw_queue.py

## Set C - 3 growth rates from different conditions
# Used for figure 2
DESC="SET C 4 gens 256 seeds 3 conditions with growth noise and D period" VARIANT="condition" FIRST_VARIANT_INDEX=0 LAST_VARIANT_INDEX=2 SINGLE_DAUGHTERS=1 N_GENS=4 N_INIT_SIMS=256 MASS_DISTRIBUTION=1 GROWTH_RATE_NOISE=1 D_PERIOD_DIVISION=1 python runscripts/fw_queue.py

## Set D - changes to RNAP and ribosome expression
# Used for figure 2
# TODO - John fill this out

## Set E - metabolism objective weighting
# Used for figure 3
DESC="SET E 8 gens 8 seeds 10 metabolism weighting values" VARIANT="metabolism_kinetic_objective_weight" FIRST_VARIANT_INDEX=0 LAST_VARIANT_INDEX=9 SINGLE_DAUGHTERS=1 N_GENS=8 N_INIT_SIMS=8 MASS_DISTRIBUTION=1 GROWTH_RATE_NOISE=0 D_PERIOD_DIVISION=0 python runscripts/fw_queue.py

## Set F - MenE expression
# Used in supplemental figure
DESC="SET_F 8 gens 8 seeds 9 menE expression values" N_INIT_SIMS=8 SINGLE_DAUGHTERS=1 N_GENS=8 RUN_VARIANT=meneParams FIRST_VARIANT_INDEX=0 LAST_VARIANT_INDEX=8 MASS_DISTRIBUTION=1 GROWTH_RATE_NOISE=1 D_PERIOD_DIVISION=1 python runscripts/fw_queue.py


## Launch the fireworks created with fw_queue.py
# rlaunch is interactive, qlaunch is distributed
# rlaunch rapidfire
qlaunch -r rapidfire --nlaunches infinite --sleep 5 --maxjobs_queue 100