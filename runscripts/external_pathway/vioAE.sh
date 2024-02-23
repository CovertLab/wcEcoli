## To Run
# chmod u+x FireWorksBox_1.sh
# ./FireWorksBox_1.sh

## Initialise the Database (is this really the only way to do this?)
cd /user/home/ig13470/wholecell3/wcm2024/wcEcoli/wholecell/fireworks/
lpad -l my_launchpad.yaml reset

## Pre Simulation Requirements
cd /user/home/ig13470/wholecell3/wcm2024/wcEcoli/

##Compile Cython code
make clean compile

# Set the $PYTHONPATH:
export PYTHONPATH="/user/home/ig13470/wholecell3/wcm2024/wcEcoli/"

## Create and Add FireWorks to the Launchpad
cd /user/home/ig13470/wholecell3/wcm2024/wcEcoli/wholecell/fireworks/

#DESC="new_gene_expression_pred_ks_338_031" PLOTS="NEW_GENE" \
#VARIANT="new_gene_expression" FIRST_VARIANT_INDEX=4 LAST_VARIANT_INDEX=7 \
#NEW_GENES="vioAE_meta" \
#CACHED_SIM_DATA=0 PARALLEL_PARCA=1 \
#SINGLE_DAUGHTERS=1 N_GENS=8 N_INIT_SIMS=1 \
#LAUNCHPAD_FILE='/user/home/ig13470/wholecell3/wcm2024/wcEcoli/wholecell/fireworks/my_launchpad.yaml' \
#python /user/home/ig13470/wholecell3/wcm2024/wcEcoli/runscripts/fireworks/fw_queue.py
#
#DESC="new_gene_expression_no_react" PLOTS="NEW_GENE" \
#VARIANT="new_gene_expression" FIRST_VARIANT_INDEX=4 LAST_VARIANT_INDEX=7 \
#NEW_GENES="vioAE" \
#CACHED_SIM_DATA=0 PARALLEL_PARCA=1 \
#SINGLE_DAUGHTERS=1 N_GENS=8 N_INIT_SIMS=1 \
#LAUNCHPAD_FILE='/user/home/ig13470/wholecell3/wcm2024/wcEcoli/wholecell/fireworks/my_launchpad.yaml' \
#python /user/home/ig13470/wholecell3/wcm2024/wcEcoli/runscripts/fireworks/fw_queue.py

DESC="wildtype 3 seeds" \
FIRST_VARIANT_INDEX=0 LAST_VARIANT_INDEX=0 \
CACHED_SIM_DATA=0 PARALLEL_PARCA=1 \
SINGLE_DAUGHTERS=1 N_GENS=8 N_INIT_SIMS=3 \
LAUNCHPAD_FILE='/user/home/ig13470/wholecell3/wcm2024/wcEcoli/wholecell/fireworks/my_launchpad.yaml' \
python /user/home/ig13470/wholecell3/wcm2024/wcEcoli/runscripts/fireworks/fw_queue.py

## Launch Rockets to collect the FireWorks (created by fw_queue.py) from the Launchpad (MongoDB hosted on MLab via Heroku) and run them on the Queue (BlueCrystal)
nohup qlaunch -r -l my_launchpad.yaml -w my_fworker.yaml -q my_qadapter.yaml rapidfire --nlaunches infinite --sleep 30 --maxjobs_queue 100
