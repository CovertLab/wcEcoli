# Script to run simulations used to generate data for figures in the rRNA paper.

lpad reset
make clean compile

## Set A - WT, no changes to rRNA operons
# Set A1 - glucose minimal media
DESC="SET A1 16 gens 32 seeds WT with glucose minimal media" \
VARIANT="wildtype" FIRST_VARIANT_INDEX=0 LAST_VARIANT_INDEX=0 \
SINGLE_DAUGHTERS=1 N_GENS=16 N_INIT_SIMS=32 \
RUN_AGGREGATE_ANALYSIS=0 \
python runscripts/fireworks/fw_queue.py

# Set A2 - rich media
DESC="SET A2 16 gens 32 seeds WT with rich media" \
VARIANT="condition" FIRST_VARIANT_INDEX=1 LAST_VARIANT_INDEX=1 \
SINGLE_DAUGHTERS=1 N_GENS=16 N_INIT_SIMS=32 \
RUN_AGGREGATE_ANALYSIS=0 \
python runscripts/fireworks/fw_queue.py

# Set A3 - shift from minimal to rich media
DESC="SET A3 24 gens 32 seeds WT with shift from minimal to rich media" \
VARIANT="timelines" FIRST_VARIANT_INDEX=28 LAST_VARIANT_INDEX=28 \
SINGLE_DAUGHTERS=1 N_GENS=24 N_INIT_SIMS=32 \
RUN_AGGREGATE_ANALYSIS=0 \
python runscripts/fireworks/fw_queue.py


## Set B - All rRNA genes are transcribed as monocistronic transcription units
# Set B1 - glucose minimal media
DESC="SET B1 16 gens 32 seeds monocistronic rRNA TUs with glucose minimal media" \
VARIANT="wildtype" FIRST_VARIANT_INDEX=0 LAST_VARIANT_INDEX=0 \
REMOVE_RRNA_OPERONS=1 \
SINGLE_DAUGHTERS=1 N_GENS=16 N_INIT_SIMS=32 \
RUN_AGGREGATE_ANALYSIS=0 \
python runscripts/fireworks/fw_queue.py

# Set B2 - rich media
DESC="SET B2 16 gens 32 seeds monocistronic rRNA TUs with rich media" \
VARIANT="condition" FIRST_VARIANT_INDEX=1 LAST_VARIANT_INDEX=1 \
REMOVE_RRNA_OPERONS=1 \
SINGLE_DAUGHTERS=1 N_GENS=16 N_INIT_SIMS=32 \
RUN_AGGREGATE_ANALYSIS=0 \
python runscripts/fireworks/fw_queue.py

# Set B3 - shift from minimal to rich media
DESC="SET B3 24 gens 32 seeds monocistronic rRNA TUs with shift from minimal to rich media" \
VARIANT="timelines" FIRST_VARIANT_INDEX=28 LAST_VARIANT_INDEX=28 \
REMOVE_RRNA_OPERONS=1 \
SINGLE_DAUGHTERS=1 N_GENS=24 N_INIT_SIMS=32 \
RUN_AGGREGATE_ANALYSIS=0 \
python runscripts/fireworks/fw_queue.py


## Set C - rRNA operon knockouts
# Set C1 - glucose minimal media
DESC="SET C1 16 gens 32 seeds rRNA operon knockouts with glucose minimal media" \
VARIANT="rrna_operon_knockout" FIRST_VARIANT_INDEX=1 LAST_VARIANT_INDEX=6 \
SINGLE_DAUGHTERS=1 N_GENS=16 N_INIT_SIMS=32 \
RUN_AGGREGATE_ANALYSIS=0 \
python runscripts/fireworks/fw_queue.py

# Set C2 - rich media
DESC="SET C2 16 gens 32 seeds rRNA operon knockouts with rich media" \
VARIANT="rrna_operon_knockout" FIRST_VARIANT_INDEX=7 LAST_VARIANT_INDEX=12 \
SINGLE_DAUGHTERS=1 N_GENS=16 N_INIT_SIMS=32 \
RUN_AGGREGATE_ANALYSIS=0 \
python runscripts/fireworks/fw_queue.py

# Set C3 - shift from minimal to rich media
DESC="SET C3 24 gens 32 seeds rRNA operon knockouts with shift from minimal to rich media" \
VARIANT="rrna_operon_knockout" FIRST_VARIANT_INDEX=13 LAST_VARIANT_INDEX=18 \
SINGLE_DAUGHTERS=1 N_GENS=24 N_INIT_SIMS=32 \
RUN_AGGREGATE_ANALYSIS=0 \
python runscripts/fireworks/fw_queue.py


## Set D - Flip chromosomal locations of rRNA operons
# Set D1 - glucose minimal media
DESC="SET D1 16 gens 32 seeds flipped rRNA locations with glucose minimal media" \
VARIANT="rrna_location" FIRST_VARIANT_INDEX=1 LAST_VARIANT_INDEX=1 \
SINGLE_DAUGHTERS=1 N_GENS=16 N_INIT_SIMS=32 \
RUN_AGGREGATE_ANALYSIS=0 \
python runscripts/fireworks/fw_queue.py

# Set D2 - rich media
DESC="SET D2 16 gens 32 seeds flipped rRNA locations with rich media" \
VARIANT="rrna_location" FIRST_VARIANT_INDEX=2 LAST_VARIANT_INDEX=2 \
SINGLE_DAUGHTERS=1 N_GENS=16 N_INIT_SIMS=32 \
RUN_AGGREGATE_ANALYSIS=0 \
python runscripts/fireworks/fw_queue.py

# Set D3 - shift from minimal to rich media
DESC="SET D3 24 gens 32 seeds flipped rRNA locations with shift from minimal to rich media" \
VARIANT="rrna_location" FIRST_VARIANT_INDEX=3 LAST_VARIANT_INDEX=3 \
SINGLE_DAUGHTERS=1 N_GENS=24 N_INIT_SIMS=32 \
RUN_AGGREGATE_ANALYSIS=0 \
python runscripts/fireworks/fw_queue.py


## Set E - Flip the orientation of rRNA operons
# Set E1 - glucose minimal media
DESC="SET E1 16 gens 32 seeds flipped rRNA orientations with glucose minimal media" \
VARIANT="rrna_orientation" FIRST_VARIANT_INDEX=1 LAST_VARIANT_INDEX=1 \
SINGLE_DAUGHTERS=1 N_GENS=16 N_INIT_SIMS=32 \
RUN_AGGREGATE_ANALYSIS=0 \
python runscripts/fireworks/fw_queue.py

# Set E2 - rich media
DESC="SET E2 16 gens 32 seeds flipped rRNA orientations with rich media" \
VARIANT="rrna_orientation" FIRST_VARIANT_INDEX=2 LAST_VARIANT_INDEX=2 \
SINGLE_DAUGHTERS=1 N_GENS=16 N_INIT_SIMS=32 \
RUN_AGGREGATE_ANALYSIS=0 \
python runscripts/fireworks/fw_queue.py

# Set E3 - shift from minimal to rich media
DESC="SET E3 24 gens 32 seeds flipped rRNA orientations with shift from minimal to rich media" \
VARIANT="rrna_orientation" FIRST_VARIANT_INDEX=3 LAST_VARIANT_INDEX=3 \
SINGLE_DAUGHTERS=1 N_GENS=24 N_INIT_SIMS=32 \
RUN_AGGREGATE_ANALYSIS=0 \
python runscripts/fireworks/fw_queue.py


## Set F - Remove the extra 5S gene (rrfF)
# Set F1 - glucose minimal media
DESC="SET F1 16 gens 32 seed remove extra 5S with glucose minimal media" \
VARIANT="wildtype" FIRST_VARIANT_INDEX=0 LAST_VARIANT_INDEX=0 \
REMOVE_RRFF=1 \
SINGLE_DAUGHTERS=1 N_GENS=16 N_INIT_SIMS=32 \
RUN_AGGREGATE_ANALYSIS=0 \
python runscripts/fireworks/fw_queue.py

# Set F2 - rich media
DESC="SET F2 16 gens 32 seed remove extra 5S with rich media" \
VARIANT="condition" FIRST_VARIANT_INDEX=1 LAST_VARIANT_INDEX=1 \
REMOVE_RRFF=1 \
SINGLE_DAUGHTERS=1 N_GENS=16 N_INIT_SIMS=32 \
RUN_AGGREGATE_ANALYSIS=0 \
python runscripts/fireworks/fw_queue.py

# Set F3 - shift from minimal to rich media
DESC="SET F3 24 gens 32 seed remove extra 5S with shift from minimal to rich media" \
VARIANT="timelines" FIRST_VARIANT_INDEX=28 LAST_VARIANT_INDEX=28 \
REMOVE_RRFF=1 \
SINGLE_DAUGHTERS=1 N_GENS=24 N_INIT_SIMS=32 \
RUN_AGGREGATE_ANALYSIS=0 \
python runscripts/fireworks/fw_queue.py

## Launch the fireworks created with fw_queue.py
# Uncomment one method - rlaunch is interactive, qlaunch is distributed
# rlaunch rapidfire
# qlaunch -r rapidfire --nlaunches infinite --sleep 5 --maxjobs_queue 100
