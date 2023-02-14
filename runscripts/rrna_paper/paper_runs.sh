# Script to run simulations used to generate data for figures in the rRNA paper.

lpad reset
make clean compile

## Set A - WT, glucose minimal media
# No changes to rRNA operons
DESC="SET A 8 gens 1 seed WT with glucose minimal media" \
VARIANT="wildtype" FIRST_VARIANT_INDEX=0 LAST_VARIANT_INDEX=0 \
SINGLE_DAUGHTERS=1 N_GENS=8 N_INIT_SIMS=1 \
RUN_AGGREGATE_ANALYSIS=0 \
python runscripts/fireworks/fw_queue.py

## Set B - no rRNA operons, glucose minimal media
# All rRNA genes are transcribed as single-gene transcription units
DESC="SET B 8 gens 1 seed single gene rRNA tus with glucose minimal media" \
VARIANT="wildtype" FIRST_VARIANT_INDEX=0 LAST_VARIANT_INDEX=0 \
REMOVE_RRNA_OPERONS=1 \
SINGLE_DAUGHTERS=1 N_GENS=8 N_INIT_SIMS=1 \
RUN_AGGREGATE_ANALYSIS=0 \
python runscripts/fireworks/fw_queue.py

## Set C - rRNA operon knockouts, glucose minimal media
# Remove one rRNA operon
DESC="SET C 8 gens 1 seed 6 rRNA operons with glucose minimal media" \
VARIANT="rrna_operon_knockout" FIRST_VARIANT_INDEX=1 LAST_VARIANT_INDEX=1 \
SINGLE_DAUGHTERS=1 N_GENS=8 N_INIT_SIMS=1 \
RUN_AGGREGATE_ANALYSIS=0 \
python runscripts/fireworks/fw_queue.py

# Remove two rRNA operons
DESC="SET C 8 gens 1 seed 5 rRNA operons with glucose minimal media" \
VARIANT="rrna_operon_knockout" FIRST_VARIANT_INDEX=8 LAST_VARIANT_INDEX=8 \
SINGLE_DAUGHTERS=1 N_GENS=8 N_INIT_SIMS=1 \
RUN_AGGREGATE_ANALYSIS=0 \
python runscripts/fireworks/fw_queue.py

# Remove three rRNA operons
DESC="SET C 8 gens 1 seed 4 rRNA operons with glucose minimal media" \
VARIANT="rrna_operon_knockout" FIRST_VARIANT_INDEX=29 LAST_VARIANT_INDEX=29 \
SINGLE_DAUGHTERS=1 N_GENS=8 N_INIT_SIMS=1 \
RUN_AGGREGATE_ANALYSIS=0 \
python runscripts/fireworks/fw_queue.py

# Remove four rRNA operons
DESC="SET C 8 gens 1 seed 3 rRNA operons with glucose minimal media" \
VARIANT="rrna_operon_knockout" FIRST_VARIANT_INDEX=64 LAST_VARIANT_INDEX=64 \
SINGLE_DAUGHTERS=1 N_GENS=8 N_INIT_SIMS=1 \
RUN_AGGREGATE_ANALYSIS=0 \
python runscripts/fireworks/fw_queue.py

# Remove five rRNA operons
DESC="SET C 8 gens 1 seed 2 rRNA operons with glucose minimal media" \
VARIANT="rrna_operon_knockout" FIRST_VARIANT_INDEX=99 LAST_VARIANT_INDEX=99 \
SINGLE_DAUGHTERS=1 N_GENS=8 N_INIT_SIMS=1 \
RUN_AGGREGATE_ANALYSIS=0 \
python runscripts/fireworks/fw_queue.py

# Remove six rRNA operons
DESC="SET C 8 gens 1 seed 1 rRNA operon with glucose minimal media" \
VARIANT="rrna_operon_knockout" FIRST_VARIANT_INDEX=120 LAST_VARIANT_INDEX=120 \
SINGLE_DAUGHTERS=1 N_GENS=8 N_INIT_SIMS=1 \
RUN_AGGREGATE_ANALYSIS=0 \
python runscripts/fireworks/fw_queue.py

## Set D - flipped locations of rRNA operons, glucose minimal media
# Change the chromosomal location of rRNA operons
DESC="SET D 8 gens 1 seed flipped rRNA locations with glucose minimal media" \
VARIANT="rrna_location" FIRST_VARIANT_INDEX=1 LAST_VARIANT_INDEX=1 \
SINGLE_DAUGHTERS=1 N_GENS=8 N_INIT_SIMS=1 \
RUN_AGGREGATE_ANALYSIS=0 \
python runscripts/fireworks/fw_queue.py

## Set E - flipped orientations of rRNA operons, glucose minimal media
# Change the orientation of rRNA operons
DESC="SET E 8 gens 1 seed flipped rRNA orientations with glucose minimal media" \
VARIANT="rrna_orientation" FIRST_VARIANT_INDEX=1 LAST_VARIANT_INDEX=1 \
SINGLE_DAUGHTERS=1 N_GENS=8 N_INIT_SIMS=1 \
RUN_AGGREGATE_ANALYSIS=0 \
python runscripts/fireworks/fw_queue.py

## Set F - remove extra 5S gene, glucose minimal media
# Remove the extra 5S gene (rrfF)
DESC="SET F 8 gens 1 seed remove extra 5S with glucose minimal media" \
VARIANT="wildtype" FIRST_VARIANT_INDEX=0 LAST_VARIANT_INDEX=0 \
REMOVE_RRFF=1 \
SINGLE_DAUGHTERS=1 N_GENS=8 N_INIT_SIMS=1 \
RUN_AGGREGATE_ANALYSIS=0 \
python runscripts/fireworks/fw_queue.py

## Launch the fireworks created with fw_queue.py
# Uncomment one method - rlaunch is interactive, qlaunch is distributed
# rlaunch rapidfire
# qlaunch -r rapidfire --nlaunches infinite --sleep 5 --maxjobs_queue 100
