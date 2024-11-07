## Initialise the Database (is this really the only way to do this?)
cd /user/home/ig13470/wholecell3/wcm_master/wcEcoli/wholecell/fireworks/
lpad -l my_launchpad.yaml reset

## Pre Simulation Requirements
cd /user/home/ig13470/wholecell3/wcm_master/wcEcoli/

##Compile Cython code
make clean compile

# Set the $PYTHONPATH:
export PYTHONPATH="/user/home/ig13470/wholecell3/wcm_master/wcEcoli/"

## Create and Add FireWorks to the Launchpad
cd /user/home/ig13470/wholecell3/wcm_master/wcEcoli/wholecell/fireworks/

DESC="new_gene_internal_shift_vioA_opoff" PLOTS="NEW_GENE" \
VARIANT="new_gene_internal_shift" FIRST_VARIANT_INDEX=1 LAST_VARIANT_INDEX=1 \
NEW_GENES="vioA" OPERONS="off" \
SINGLE_DAUGHTERS=1 N_GENS=1 SEED=0 N_INIT_SIMS=1 \
LAUNCHPAD_FILE='/user/home/ig13470/wholecell3/wcm_master/wcEcoli/wholecell/fireworks/my_launchpad.yaml' \
python /user/home/ig13470/wholecell3/wcm_master/wcEcoli/runscripts/fireworks/fw_queue.py

DESC="new_gene_internal_shift_vioA_opon" PLOTS="NEW_GENE" \
VARIANT="new_gene_internal_shift" FIRST_VARIANT_INDEX=1 LAST_VARIANT_INDEX=1 \
NEW_GENES="vioA" OPERONS="on" \
SINGLE_DAUGHTERS=1 N_GENS=1 SEED=0 N_INIT_SIMS=1 \
LAUNCHPAD_FILE='/user/home/ig13470/wholecell3/wcm_master/wcEcoli/wholecell/fireworks/my_launchpad.yaml' \
python /user/home/ig13470/wholecell3/wcm_master/wcEcoli/runscripts/fireworks/fw_queue.py

DESC="new_gene_internal_shift_vioAB_opoff" PLOTS="NEW_GENE" \
VARIANT="new_gene_internal_shift" FIRST_VARIANT_INDEX=1 LAST_VARIANT_INDEX=1 \
NEW_GENES="vioAB" OPERONS="off" \
SINGLE_DAUGHTERS=1 N_GENS=1 SEED=0 N_INIT_SIMS=1 \
LAUNCHPAD_FILE='/user/home/ig13470/wholecell3/wcm_master/wcEcoli/wholecell/fireworks/my_launchpad.yaml' \
python /user/home/ig13470/wholecell3/wcm_master/wcEcoli/runscripts/fireworks/fw_queue.py

DESC="new_gene_internal_shift_vioAB_opon" PLOTS="NEW_GENE" \
VARIANT="new_gene_internal_shift" FIRST_VARIANT_INDEX=1 LAST_VARIANT_INDEX=1 \
NEW_GENES="vioAB" OPERONS="on" \
SINGLE_DAUGHTERS=1 N_GENS=1 SEED=0 N_INIT_SIMS=1 \
LAUNCHPAD_FILE='/user/home/ig13470/wholecell3/wcm_master/wcEcoli/wholecell/fireworks/my_launchpad.yaml' \
python /user/home/ig13470/wholecell3/wcm_master/wcEcoli/runscripts/fireworks/fw_queue.py

DESC="new_gene_internal_shift_vioABC_opoff" PLOTS="NEW_GENE" \
VARIANT="new_gene_internal_shift" FIRST_VARIANT_INDEX=1 LAST_VARIANT_INDEX=1 \
NEW_GENES="vioABC" OPERONS="off" \
SINGLE_DAUGHTERS=1 N_GENS=1 SEED=0 N_INIT_SIMS=1 \
LAUNCHPAD_FILE='/user/home/ig13470/wholecell3/wcm_master/wcEcoli/wholecell/fireworks/my_launchpad.yaml' \
python /user/home/ig13470/wholecell3/wcm_master/wcEcoli/runscripts/fireworks/fw_queue.py

DESC="new_gene_internal_shift_vioABC_opon" PLOTS="NEW_GENE" \
VARIANT="new_gene_internal_shift" FIRST_VARIANT_INDEX=1 LAST_VARIANT_INDEX=1 \
NEW_GENES="vioABC" OPERONS="on" \
SINGLE_DAUGHTERS=1 N_GENS=1 SEED=0 N_INIT_SIMS=1 \
LAUNCHPAD_FILE='/user/home/ig13470/wholecell3/wcm_master/wcEcoli/wholecell/fireworks/my_launchpad.yaml' \
python /user/home/ig13470/wholecell3/wcm_master/wcEcoli/runscripts/fireworks/fw_queue.py

DESC="new_gene_internal_shift_vioABCDE_opoff" PLOTS="NEW_GENE" \
VARIANT="new_gene_internal_shift" FIRST_VARIANT_INDEX=1 LAST_VARIANT_INDEX=1 \
NEW_GENES="vioABCDE" OPERONS="off" \
SINGLE_DAUGHTERS=1 N_GENS=1 SEED=0 N_INIT_SIMS=1 \
LAUNCHPAD_FILE='/user/home/ig13470/wholecell3/wcm_master/wcEcoli/wholecell/fireworks/my_launchpad.yaml' \
python /user/home/ig13470/wholecell3/wcm_master/wcEcoli/runscripts/fireworks/fw_queue.py

DESC="new_gene_internal_shift_vioABCDE_opon" PLOTS="NEW_GENE" \
VARIANT="new_gene_internal_shift" FIRST_VARIANT_INDEX=1 LAST_VARIANT_INDEX=1 \
NEW_GENES="vioABCDE" OPERONS="on" \
SINGLE_DAUGHTERS=1 N_GENS=1 SEED=0 N_INIT_SIMS=1 \
LAUNCHPAD_FILE='/user/home/ig13470/wholecell3/wcm_master/wcEcoli/wholecell/fireworks/my_launchpad.yaml' \
python /user/home/ig13470/wholecell3/wcm_master/wcEcoli/runscripts/fireworks/fw_queue.py


## Launch Rockets to collect the FireWorks (created by fw_queue.py) from the Launchpad (MongoDB hosted on MLab via Heroku) and run them on the Queue (BlueCrystal)
nohup qlaunch -r -l my_launchpad.yaml -w my_fworker.yaml -q my_qadapter.yaml rapidfire --nlaunches infinite --sleep 30 --maxjobs_queue 100