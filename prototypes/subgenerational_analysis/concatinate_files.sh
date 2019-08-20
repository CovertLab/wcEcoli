#!/bin/bash

for i in `ls /Users/mialydefelice/Documents/code_repositories/wcEcoli/out/counts/wildtype_000000/count_out/000000/*_gen_data_rna.tsv | sort -V`; 
do arr+=($i)
echo $(cat $i | wc -l); 
done;
#declare -p arr


cat $(ls /Users/mialydefelice/Documents/code_repositories/wcEcoli/out/counts/wildtype_000000/count_out/000000/*_gen_data_rna.tsv | sort -V) > /Users/mialydefelice/Documents/code_repositories/wcEcoli/out/counts/wildtype_000000/count_out/000000/merged_rna.tsv

echo $(cat /Users/mialydefelice/Documents/code_repositories/wcEcoli/out/counts/wildtype_000000/count_out/000000/merged_rna.tsv | wc -l)