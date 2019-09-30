#!/bin/bash

: '
Author: Mialy DeFelice
Date: 082119
Function:
Basically will take outputs downloaded from the cloud and merge them all together.
There are no checks to make sure the file exists due to a problem with incomplete files not being rewritten.

Add in additional checks to make sure the vectors and matrices are of the proper length before proceeding.
'
base_output_dir="/Users/mialydefelice/Documents/code_repositories/wcEcoli/out/counts/wildtype_000000/count_out_2"
save_dir="/Users/mialydefelice/Documents/code_repositories/wcEcoli/out/counts/wildtype_000000/count_out_concatenated/"

concatenate_data_files () {
	echo "Making $2"
	cat $(ls $(find "$i" -iname $1 | sort -V)) > $2
}
	
for i in $(find "${base_output_dir}"* -depth 1 -type d); do
	base=`basename $i`
	num_timepoints=$(cat $"$i/time_info.tsv" | wc -l)
	# Merge counts data from multiple generations into single files
	# For RNA, protein and complexes
	
	#num_files_to_merge=$(find "$i" -maxdepth 1 -name '*_gen_data_rna.tsv' | wc -l)
	#echo $num_files_to_merge

	for j in $(find "${i}"* -maxdepth 1 -type d); do
		#create multigen merge of data for protein and rna
		#check if the merged file is the correct length, if its not run till its made.
		#Note this can be dangerous, so run while watching after it.

		# --- RNA --- 
		concatenate_data_files '*_gen_data_rna.tsv' "$j/merged_rna_counts.tsv"
		# The following is a check to ensure the concatenation file made is the correct lenght.
		while [ $(cat $"$j/merged_rna_counts.tsv" | wc -l) != $num_timepoints ]; do
			echo "WARNING: There was an issue making $j/merged_rna_counts.tsv, reconcatenating this file now!"
			concatenate_data_files '*_gen_data_rna.tsv' "$j/merged_rna_counts.tsv"
		done;

		echo "$j/merged_rna_counts.tsv file successfully made"

		# --- Protein --- 
		concatenate_data_files '*_gen_data_protein.tsv' "$j/merged_protein_counts.tsv"
		while [ $(cat $"$j/merged_protein_counts.tsv" | wc -l) != $num_timepoints ]; do
			echo "WARNING: There was an issue making $j/merged_protein_counts.tsv, reconcatenating this file now!"
			concatenate_data_files '*_gen_data_protein.tsv' "$j/merged_protein_counts.tsv"
		done;

		echo "$j/merged_protein_counts.tsv file successfully made"
		
		# --- Complexes --- 
		concatenate_data_files '*_gen_data_complex.tsv' "$j/merged_complex_counts.tsv"
		while [ $(cat $"$j/merged_complex_counts.tsv" | wc -l) != $num_timepoints ]; do
			echo "WARNING: There was an issue making $j/merged_complex_counts.tsv, reconcatenating this file now!"
			concatenate_data_files '*_gen_data_complex.tsv' "$j/merged_complex_counts.tsv"
		done;

		echo "$j/merged_complex_counts.tsv file successfully made"
		
	done;

#----- Add a header to time_info file (assumes one does not exist)
	echo "Making $i/gen_info_w_header.tsv"
	gen_info_file=$(find "$i" -iname 'gen_info.tsv')
	echo 'gen' > "$i/gen_info_w_header.tsv"
	cat $gen_info_file >> "$i/gen_info_w_header.tsv"


#----- Add a header to time_info file (assumes one does not exist):
	echo "Making $i/time_info_w_header.tsv"
	time_info_file=$(find "$i" -iname 'time_info.tsv')
	echo 'time' > "$i/time_info_w_header.tsv" 
	cat $time_info_file >> "$i/time_info_w_header.tsv"


#----- Transpose RNA ID File so it can be merged with the counts later:
	echo "Making $i/transposed_ids_rnas.tsv"
	cut -f1 "$i/ids_rnas.tsv" | paste -s - > "$i/transposed_ids_rnas.tsv"


#----- Transpose Protein ID File so it can be merged with the counts later:
	echo "Making $i/transposed_ids_proteins.tsv"
	cut -f1 "$i/ids_proteins.tsv" | paste -s - > "$i/transposed_ids_proteins.tsv"



#----- Transpose Complex ID File so it can be merged with the counts later:
	echo "Making $i/transposed_ids_complex.tsv"
	cut -f1 "$i/ids_complex.tsv" | paste -s - > "$i/transposed_ids_complex.tsv"


#----- Merge RNA IDs to Counts
output_file_rna_id="$i/rna_ids_and_counts.tsv"
	echo "Making $output_file_rna_id"
	rna_id_file=$(find "$i" -iname 'transposed_ids_rnas.tsv')
	merged_rna_counts_file="$i/merged_rna_counts.tsv"
	cat $rna_id_file $merged_rna_counts_file >> $output_file_rna_id


#----- Merge Protein IDs to Counts
output_file_protein_id="$i/protein_ids_and_counts.tsv"
	echo "Making $output_file_protein_id"
	protein_id_file=$(find "$i" -iname 'transposed_ids_proteins.tsv')
	merged_protein_counts_file="$i/merged_protein_counts.tsv"
	cat $protein_id_file $merged_protein_counts_file >> $output_file_protein_id


#----- Merge Complex IDs to Counts
output_file_complex_id="$i/complex_ids_and_counts.tsv"
	echo "Making $output_file_complex_id"
	complex_id_file=$(find "$i" -iname 'transposed_ids_complex.tsv')
	merged_complex_counts_file="$i/merged_complex_counts.tsv"
	cat $complex_id_file $merged_complex_counts_file >> $output_file_complex_id


#----- Add Time and Gen info to each merged file.
output_all_rna_data=""$save_dir${base}_all_rna_data.tsv""
	echo "Making $output_all_rna_data"
	paste "$i/time_info_w_header.tsv" "$i/gen_info_w_header.tsv" "$i/rna_ids_and_counts.tsv" >> $output_all_rna_data 


#----- Add Time and Gen info to each merged file.
output_all_protein_data=""$save_dir${base}_all_protein_data.tsv""
	echo "Making $output_all_protein_data"
	paste "$i/time_info_w_header.tsv" "$i/gen_info_w_header.tsv" "$i/protein_ids_and_counts.tsv" >> $output_all_protein_data 


#----- Add Time and Gen info to each merged file.
output_all_complex_data=""$save_dir${base}_all_complex_data.tsv""
	echo "Making $output_all_complex_data"
	paste "$i/time_info_w_header.tsv" "$i/gen_info_w_header.tsv" "$i/complex_ids_and_counts.tsv" >> $output_all_complex_data 


#---- To save space, remove some precursor files.
if [ -f $output_file_rna_id ]; then
	rm $output_file_rna_id
else
	echo "$output_file_rna_id does not exist"
fi

#---- To save space, remove some precursor files.
if [ -f $output_file_protein_id ]; then
	rm $output_file_protein_id
else
	echo "$output_file_protein_id does not exist"
fi

#---- To save space, remove some precursor files.
if [ -f $output_file_complex_id ]; then
	rm $output_file_complex_id
else
	echo "$output_file_complex_id does not exist"
fi

	#---- Use for checking file dimensions

	: '
	
	echo "rows and columns of merged rna counts id file"
	echo $(cat $"$i/merged_rna_ids.tsv" | wc -l)
	echo $(head -1 $"$i/merged_rna_ids.tsv" | sed ''s/,/\t/g'' | wc -w)
	'

done;