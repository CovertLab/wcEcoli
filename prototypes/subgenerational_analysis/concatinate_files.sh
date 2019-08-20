#!/bin/bash


for i in `ls /Users/mialydefelice/Documents/code_repositories/wcEcoli/out/counts/wildtype_000000/count_out/000000/*_gen_data_rna.tsv | sort -V`; 
do arr+=($i)
#echo $i; 
done;
#declare -p arr

for ((n=0; n<=32; n++))
do
	echo "${arr[$n]}"
done;

cat $(ls /Users/mialydefelice/Documents/code_repositories/wcEcoli/out/counts/wildtype_000000/count_out/000000/*_gen_data_rna.tsv | sort -V) > /Users/mialydefelice/Documents/code_repositories/wcEcoli/out/counts/wildtype_000000/count_out/000000/merged_rna.tsv