#! /bin/bash
# Script to run many variants in parallel for parameter sensitivity analysis.
# Fireworks is too slow and checks too many task dependencies for this scale.

# Usage (from wcEcoli home directory):
#   runscripts/paper/sensitivity.sh [output dir] [start variant] [final variant]

set -eu

# Required arguments
out_dir=$1
start_var=$2
end_var=$3

date

# Run simulation and remove sim_data to save space
function simulation {
	variant="param_sensitivity"
	python runscripts/manual/runSim.py $1 --variant $variant $2 $2 > out/$1/metadata/$2.log
	rm out/$1/${variant}_$(printf "%06d" $2)/kb/simData_Modified.cPickle
}
export -f simulation

# Create sim_data
python runscripts/manual/runFitter.py --disable-ribosome-fitting --disable-rnapoly-fitting $out_dir

# Run simulation variants in parallel
seq $start_var $end_var | parallel simulation $out_dir

date
