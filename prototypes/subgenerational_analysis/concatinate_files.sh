#!/bin/bash
base_output_dir="/Users/mialydefelice/Documents/code_repositories/wcEcoli/out/counts/wildtype_000000/count_out/"
for i in $(find "${base_output_dir}"* -maxdepth 1 -type d); do
	for j in $(find "${i}"* -maxdepth 1 -type d); do
		output_file_name_rna="$j""/merged_rna.tsv";
		output_file_name_protein="$j""/merged_protein.tsv";
		#create multigen merge of data for protein and rna
		cat $(ls $(find "$i" -iname '*_gen_data_rna.tsv' | sort -V)) > $output_file_name_rna;
		cat $(ls $(find "$i" -iname '*_gen_data_protein.tsv' | sort -V)) > $output_file_name_protein;
	done;
done;




#for i in `ls /*_gen_data_rna.tsv | sort -V`; 
#do arr+=($i)
#echo $(cat $i | wc -l);
#echo $(head -1 $i | sed 's/,/\t/g' | wc -w)
#done;
#declare -p arr


#cat $(ls /Users/mialydefelice/Documents/code_repositories/wcEcoli/out/counts/wildtype_000000/count_out/000000/*_gen_data_rna.tsv | sort -V) > /Users/mialydefelice/Documents/code_repositories/wcEcoli/out/counts/wildtype_000000/count_out/000000/merged_rna.tsv



#echo $(cat /Users/mialydefelice/Documents/code_repositories/wcEcoli/out/counts/wildtype_000000/count_out/000000/merged_rna.tsv | wc -l)
#echo $(head -1 /Users/mialydefelice/Documents/code_repositories/wcEcoli/out/counts/wildtype_000000/count_out/000000/merged_rna.tsv | sed 's/,/\t/g' | wc -w)