## To Run
# chmod u+x FireWorksBox_1.sh
# ./FireWorksBox_1.sh


## Create and Add FireWorks to the Launchpad
cd $HOME/wcEcoli/wholecell/fireworks

# Reset the database
yes | lpad -l my_launchpad.yaml reset

export DESC="bpsA_KI_20gen" \
PLOTS="NEW_GENE" \
VARIANT="new_gene_internal_shift" \
FIRST_VARIANT_INDEX=1 \
LAST_VARIANT_INDEX=1 \
NEW_GENES="bpsA" \
CACHED_SIM_DATA=0 \
PARALLEL_PARCA=1 \
SINGLE_DAUGHTERS=1 \
N_GENS=20 \
N_INIT_SIMS=1 \
LAUNCHPAD_FILE='/user/home/co18263/wcEcoli/wholecell/fireworks/my_launchpad.yaml' \

python $HOME/wcEcoli/runscripts/fireworks/fw_queue.py

## Launch Rockets to collect the FireWorks (created by fw_queue.py) from the Launchpad (MongoDB hosted on MLab via Heroku) and run them on the Queue 
nohup qlaunch -r -l my_launchpad.yaml -w my_fworker.yaml -q my_qadapter.yaml rapidfire --nlaunches infinite --sleep 30 --maxjobs_queue 100
