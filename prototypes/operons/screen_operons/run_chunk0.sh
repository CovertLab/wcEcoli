#!/bin/bash

num_sim=0
sim_base_name="screen_5_tus_chunk0"

for entry in prototypes/operons/polycistron_chunks_shuffled/chunk_0/*
do
  echo "$entry"
  sim_name="${sim_base_name}_${num_sim}"
  echo $sim_name

  python runscripts/reconstruction/polycistronic_rnas/create_tu_files.py $entry
  python runscripts/manual/runParca.py $sim_name
  python runscripts/manual/runSim.py -g 10 $sim_name


  num_sim=$((num_sim+1))


done
