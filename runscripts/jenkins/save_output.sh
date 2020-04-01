#! /bin/bash

# Script to move completed Jenkins runs to an output directory.
# Also adds the output to a tar archive in the output directory 3x per month
# to prevent purging of the archive and loss of data.

# Usage: ./save_output.sh src_dir dest_dir
#   src_dir: can be relative and will save all directories inside
#   dest_dir: must be absolute path and end in /

set -eu

src_dir=$1
dest_dir=$2

# Ensure dest_dir exists
mkdir -p $dest_dir

# Get date string to only save in tar every 11 days (3x per month)
# (eg 20041 for all dates between 4/11/20 - 4/21/20)
date_str="$(date +%y%m)$((10#$(date +%d) / 11))"

# Files for tar backup
date_file="${dest_dir}tar_dates.txt"
tar_file="${dest_dir}sims.tar"

# Ensure date_file exists if a new directory
touch $date_file

# If current date string does not have a sim saved in tar,
# then add the output from this run to the archive
cd $src_dir
if [ -z "$(grep $date_str $date_file)" ]; then
    echo $date_str >> $date_file
    tar rf $tar_file *
fi

mv * $dest_dir
