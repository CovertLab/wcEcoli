## To Run
# chmod u+x FireWorksBox_1.sh
# ./FireWorksBox_1.sh


## Create and Add FireWorks to the Launchpad
cd wholecell/fireworks/

DESC="new_gene_expression_pred_ks_338_031" PLOTS="NEW_GENE" \
VARIANT="new_gene_internal_shift" FIRST_VARIANT_INDEX=1 LAST_VARIANT_INDEX=3 \
NEW_GENES="vioAE_meta" \
CACHED_SIM_DATA=0 PARALLEL_PARCA=1 \
SINGLE_DAUGHTERS=1 N_GENS=8 N_INIT_SIMS=1 \
LAUNCHPAD_FILE='/wholecell/fireworks/my_launchpad.yaml' \
python runscripts/fireworks/fw_queue.py

DESC="new_gene_expression_no_react" PLOTS="NEW_GENE" \
VARIANT="new_gene_internal_shift" FIRST_VARIANT_INDEX=1 LAST_VARIANT_INDEX=3 \
NEW_GENES="vioAE" \
CACHED_SIM_DATA=0 PARALLEL_PARCA=1 \
SINGLE_DAUGHTERS=1 N_GENS=8 N_INIT_SIMS=1 \
LAUNCHPAD_FILE='/wholecell/fireworks/my_launchpad.yaml' \
python runscripts/fireworks/fw_queue.py


## Launch Rockets to collect the FireWorks (created by fw_queue.py) from the Launchpad (MongoDB hosted on MLab via Heroku) and run them on the Queue 
nohup qlaunch -r -l my_launchpad.yaml -w my_fworker.yaml -q my_qadapter.yaml rapidfire --nlaunches infinite --sleep 30 --maxjobs_queue 100
