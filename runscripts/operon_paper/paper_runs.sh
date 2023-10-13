# Script to run simulations used to generate data for figures in the operon
# paper.

lpad reset
make clean compile

## Set A - no operons, glucose minimal media
# Used for comparisons in Figures 1, 2, 3, 4, 5
DESC="SET A 8 gens 128 seeds operons off with glucose minimal media" \
OPERONS="off" \
VARIANT="wildtype" FIRST_VARIANT_INDEX=0 LAST_VARIANT_INDEX=0 \
TRNA_ATTENUATION=0 PPGPP_REGULATION=0 \
MECHANISTIC_TRANSLATION_SUPPLY=0 MECHANISTIC_AA_TRANSPORT=0 \
AA_SUPPLY_IN_CHARGING=0 D_PERIOD_DIVISION=0 \
SINGLE_DAUGHTERS=1 N_GENS=8 N_INIT_SIMS=128 \
RUN_AGGREGATE_ANALYSIS=0 \
python runscripts/fireworks/fw_queue.py

## Set B - no operons, rich media
# Used for comparisons in Figures 1, 2, 3
DESC="SET B 8 gens 128 seeds operons off with rich media" \
OPERONS="off" \
VARIANT="condition" FIRST_VARIANT_INDEX=1 LAST_VARIANT_INDEX=1 \
TRNA_ATTENUATION=0 PPGPP_REGULATION=0 \
MECHANISTIC_TRANSLATION_SUPPLY=0 MECHANISTIC_AA_TRANSPORT=0 \
AA_SUPPLY_IN_CHARGING=0 D_PERIOD_DIVISION=0 \
SINGLE_DAUGHTERS=1 N_GENS=8 N_INIT_SIMS=128 \
RUN_AGGREGATE_ANALYSIS=0 \
python runscripts/fireworks/fw_queue.py

## Set C - operons v1, glucose minimal media
# Used for comparisons in Figures 2
DESC="SET C 8 gens 128 seeds operons v1 with glucose minimal media" \
OPERONS="v1" \
VARIANT="wildtype" FIRST_VARIANT_INDEX=0 LAST_VARIANT_INDEX=0 \
TRNA_ATTENUATION=0 PPGPP_REGULATION=0 \
MECHANISTIC_TRANSLATION_SUPPLY=0 MECHANISTIC_AA_TRANSPORT=0 \
AA_SUPPLY_IN_CHARGING=0 D_PERIOD_DIVISION=0 \
SINGLE_DAUGHTERS=1 N_GENS=8 N_INIT_SIMS=128 \
RUN_AGGREGATE_ANALYSIS=0 \
python runscripts/fireworks/fw_queue.py

## Set D - operons v1, rich media
# Used for comparisons in Figures 2
DESC="SET D 8 gens 128 seeds operons v1 with rich media" \
OPERONS="v1" \
VARIANT="condition" FIRST_VARIANT_INDEX=1 LAST_VARIANT_INDEX=1 \
TRNA_ATTENUATION=0 PPGPP_REGULATION=0 \
MECHANISTIC_TRANSLATION_SUPPLY=0 MECHANISTIC_AA_TRANSPORT=0 \
AA_SUPPLY_IN_CHARGING=0 D_PERIOD_DIVISION=0 \
SINGLE_DAUGHTERS=1 N_GENS=8 N_INIT_SIMS=128 \
RUN_AGGREGATE_ANALYSIS=0 \
python runscripts/fireworks/fw_queue.py

## Set E - operons v2, glucose minimal media
# Used for comparisons in Figures 2, 3
DESC="SET E 8 gens 128 seeds operons v2 with glucose minimal media" \
OPERONS="v2" \
VARIANT="wildtype" FIRST_VARIANT_INDEX=0 LAST_VARIANT_INDEX=0 \
TRNA_ATTENUATION=0 PPGPP_REGULATION=0 \
MECHANISTIC_TRANSLATION_SUPPLY=0 MECHANISTIC_AA_TRANSPORT=0 \
AA_SUPPLY_IN_CHARGING=0 D_PERIOD_DIVISION=0 \
SINGLE_DAUGHTERS=1 N_GENS=8 N_INIT_SIMS=128 \
RUN_AGGREGATE_ANALYSIS=0 \
python runscripts/fireworks/fw_queue.py

## Set F - operons v2, rich media
# Used for comparisons in Figures 2, 3
DESC="SET F 8 gens 128 seeds operons v2 with rich media" \
OPERONS="v2" \
VARIANT="condition" FIRST_VARIANT_INDEX=1 LAST_VARIANT_INDEX=1 \
TRNA_ATTENUATION=0 PPGPP_REGULATION=0 \
MECHANISTIC_TRANSLATION_SUPPLY=0 MECHANISTIC_AA_TRANSPORT=0 \
AA_SUPPLY_IN_CHARGING=0 D_PERIOD_DIVISION=0 \
SINGLE_DAUGHTERS=1 N_GENS=8 N_INIT_SIMS=128 \
RUN_AGGREGATE_ANALYSIS=0 \
python runscripts/fireworks/fw_queue.py

## Set G - operons v3, glucose minimal media
# Used for comparisons in Figure 3
DESC="SET G 8 gens 128 seeds operons v3 with glucose minimal media" \
OPERONS="v3" \
VARIANT="wildtype" FIRST_VARIANT_INDEX=0 LAST_VARIANT_INDEX=0 \
TRNA_ATTENUATION=0 PPGPP_REGULATION=0 \
MECHANISTIC_TRANSLATION_SUPPLY=0 MECHANISTIC_AA_TRANSPORT=0 \
AA_SUPPLY_IN_CHARGING=0 D_PERIOD_DIVISION=0 \
SINGLE_DAUGHTERS=1 N_GENS=8 N_INIT_SIMS=128 \
RUN_AGGREGATE_ANALYSIS=0 \
python runscripts/fireworks/fw_queue.py

## Set H - operons v3, rich media
# Used for comparisons in Figure 3
DESC="SET H 8 gens 128 seeds operons v3 with rich media" \
OPERONS="v3" \
VARIANT="condition" FIRST_VARIANT_INDEX=1 LAST_VARIANT_INDEX=1 \
TRNA_ATTENUATION=0 PPGPP_REGULATION=0 \
MECHANISTIC_TRANSLATION_SUPPLY=0 MECHANISTIC_AA_TRANSPORT=0 \
AA_SUPPLY_IN_CHARGING=0 D_PERIOD_DIVISION=0 \
SINGLE_DAUGHTERS=1 N_GENS=8 N_INIT_SIMS=128 \
RUN_AGGREGATE_ANALYSIS=0 \
python runscripts/fireworks/fw_queue.py

## Set I - operons on (final version), glucose minimal media
# Used for comparisons in Figures 1, 4, 5
DESC="SET I 8 gens 128 seeds operons on with glucose minimal media" \
OPERONS="on" \
VARIANT="wildtype" FIRST_VARIANT_INDEX=0 LAST_VARIANT_INDEX=0 \
TRNA_ATTENUATION=0 PPGPP_REGULATION=0 \
MECHANISTIC_TRANSLATION_SUPPLY=0 MECHANISTIC_AA_TRANSPORT=0 \
AA_SUPPLY_IN_CHARGING=0 D_PERIOD_DIVISION=0 \
SINGLE_DAUGHTERS=1 N_GENS=8 N_INIT_SIMS=128 \
RUN_AGGREGATE_ANALYSIS=0 \
python runscripts/fireworks/fw_queue.py

## Set J - operons on (final version), rich media
# Used for comparisons in Figure 4
DESC="SET J 8 gens 128 seeds operons on with rich media" \
OPERONS="on" \
VARIANT="condition" FIRST_VARIANT_INDEX=1 LAST_VARIANT_INDEX=1 \
TRNA_ATTENUATION=0 PPGPP_REGULATION=0 \
MECHANISTIC_TRANSLATION_SUPPLY=0 MECHANISTIC_AA_TRANSPORT=0 \
AA_SUPPLY_IN_CHARGING=0 D_PERIOD_DIVISION=0 \
SINGLE_DAUGHTERS=1 N_GENS=8 N_INIT_SIMS=128 \
RUN_AGGREGATE_ANALYSIS=0 \
python runscripts/fireworks/fw_queue.py

## Set K - operons off, glucose minimal media, longer sims
# Used for comparisons in Figures 5
DESC="SET K 32 gens 8 seeds operons off with glucose minimal media" \
OPERONS="off" \
VARIANT="wildtype" FIRST_VARIANT_INDEX=0 LAST_VARIANT_INDEX=0 \
TRNA_ATTENUATION=0 PPGPP_REGULATION=0 \
MECHANISTIC_TRANSLATION_SUPPLY=0 MECHANISTIC_AA_TRANSPORT=0 \
AA_SUPPLY_IN_CHARGING=0 D_PERIOD_DIVISION=0 \
SINGLE_DAUGHTERS=1 N_GENS=32 N_INIT_SIMS=8 \
RUN_AGGREGATE_ANALYSIS=0 \
python runscripts/fireworks/fw_queue.py

## Set L - operons on, glucose minimal media, longer sims
# Used for comparisons in Figure 5
DESC="SET L 32 gens 8 seeds operons on with glucose minimal media" \
OPERONS="on" \
VARIANT="wildtype" FIRST_VARIANT_INDEX=0 LAST_VARIANT_INDEX=0 \
TRNA_ATTENUATION=0 PPGPP_REGULATION=0 \
MECHANISTIC_TRANSLATION_SUPPLY=0 MECHANISTIC_AA_TRANSPORT=0 \
AA_SUPPLY_IN_CHARGING=0 D_PERIOD_DIVISION=0 \
SINGLE_DAUGHTERS=1 N_GENS=32 N_INIT_SIMS=8 \
RUN_AGGREGATE_ANALYSIS=0 \
python runscripts/fireworks/fw_queue.py


## Launch the fireworks created with fw_queue.py
# Uncomment one method - rlaunch is interactive, qlaunch is distributed
# rlaunch rapidfire
# qlaunch -r rapidfire --nlaunches infinite --sleep 5 --maxjobs_queue 100
