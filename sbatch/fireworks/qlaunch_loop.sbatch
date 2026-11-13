#!/usr/bin/bash
#SBATCH --job-name=qlaunch_loop
#SBATCH --time=2-00:00:00
#SBATCH -p normal
#SBATCH -c 1
#SBATCH --mem=4GB
#SBATCH -o /scratch/users/aj0204/fw_logs/qlaunch_loop-%j.out
source /home/groups/mcovert/pyenv/versions/wcEcoli3/bin/activate
cd /home/users/aj0204/repos/wcEcoli

python runscripts/fireworks/run_qlaunch.py 64
