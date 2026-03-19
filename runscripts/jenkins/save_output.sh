#! /bin/bash

# Script to move completed Jenkins runs to an output directory.

# Usage: ./save_output.sh src_dir dest_dir
#   src_dir: will save all directories inside
#   dest_dir: directory to move sims to

set -e

src_dir=$1
dest_dir=$2

# Check correct number of inputs
if [ $# -ne 2 ]; then
    echo "Must supply src and dest dir"
    exit
fi

# Make dest_dir an absolute path for when directory is changed
dest_dir=$(realpath $dest_dir)

# Ensure dest_dir exists
mkdir -p $dest_dir

cd $src_dir
mv * $dest_dir
