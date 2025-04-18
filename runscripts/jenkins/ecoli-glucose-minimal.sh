set -e

source runscripts/jenkins/setup-environment.sh
sh runscripts/jenkins/fireworks-config.sh glucose

echo y | lpad reset

DESC="Daily build." SINGLE_DAUGHTERS=1 N_GENS=16 MASS_DISTRIBUTION=0 COMPRESS_OUTPUT=1 PLOTS=ACTIVE RAISE_ON_TIME_LIMIT=1 python runscripts/fireworks/fw_queue.py

bash runscripts/jenkins/run-fireworks.sh

cp out/2*/kb/rawData.cPickle.bz2 /scratch/groups/mcovert/wc_ecoli/cached/
bunzip2 -f /scratch/groups/mcovert/wc_ecoli/cached/rawData.cPickle.bz2
chmod 444 /scratch/groups/mcovert/wc_ecoli/cached/rawData.cPickle
cp out/2*/kb/simData.cPickle.bz2 /scratch/groups/mcovert/wc_ecoli/cached/
bunzip2 -f /scratch/groups/mcovert/wc_ecoli/cached/simData.cPickle.bz2
chmod 444 /scratch/groups/mcovert/wc_ecoli/cached/simData.cPickle

runscripts/jenkins/save_output.sh out/ /scratch/groups/mcovert/wc_ecoli/daily_build/
