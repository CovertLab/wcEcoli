#! /usr/bin/env bash

set -eu

script_dir=$(dirname $(dirname $(realpath $0)))
out_dir=$(dirname $script_dir)/out
dir1="check-reproducibility-1"
dir2="check-reproducibility-2"
sim_out_dir="wildtype_000000/000000/generation_000001/000000/simOut"

# Run parca twice and check that output is consistent from run to run
python $script_dir/manual/runParca.py -c4 $dir1 --save
python $script_dir/manual/runParca.py -c4 $dir2 --save
$script_dir/debug/comparePickles.py out/$dir1/kb out/$dir2/kb

# Run entire simulation for each parca output
python $script_dir/manual/runSim.py $dir1
python $script_dir/manual/runSim.py $dir2

# Run short daughters to check that everything including division is consistent
python $script_dir/manual/runDaughter.py $dir1 --length-sec 10
python $script_dir/manual/runDaughter.py $dir2 --length-sec 10
$script_dir/debug/diff_simouts.py $out_dir/$dir1/$sim_out_dir $out_dir/$dir2/$sim_out_dir
