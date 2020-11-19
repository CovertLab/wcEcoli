#!/bin/bash

: '
Author: Mialy DeFelice
Date: 082119
Function:
Yay this is my first bash script!
Basically will take outputs downloaded from the cloud and merge them all together.
The way I did this was to save as much data as possible and to avoid overwriting
anything, to be safe. But would be more space efficient to overwrite many of the files
or delete them in the end.
Its also longer than it needs to be, so I can put in some print statements, and do 
	a lot of gratitious checks.
TODO:
-Take base_output_dir as an argument?
'
base_output_dir="/Users/taryn/operons/2019-08-22_100seed_subgen/simdata/"
save_dir="/Users/taryn/operons/2019-08-22_100seed_subgen/concatenated_files/"

for i in $(find "${base_output_dir}"* -maxdepth 1 -type d); do

	# Merge counts data from multiple generations into single files
	# For RNA and protein

	for j in $(find "${i}"* -maxdepth 1 -type d); do
		#create multigen merge of data for protein and rna
		#only if it doesnt already exist.
		if [ ! -f "$j/merged_rna_counts.tsv" ]; then
			echo "Making $j/merged_rna_counts.tsv"
			cat $(ls $(find "$i" -iname '*_gen_data_rna.tsv' | sort -V)) > "$j/merged_rna_counts.tsv";
		else
			echo "$j/merged_rna_counts.tsv file already exists"
		fi
		if [ ! -f "$j/merged_protein_counts.tsv" ]; then
			echo "Making $j/merged_protein_counts.tsv"
			cat $(ls $(find "$i" -iname '*_gen_data_protein.tsv' | sort -V)) > "$j/merged_protein_counts.tsv";
		else
			echo "$j/merged_protein_counts.tsv.tsv file already exists"
		fi
	done;

	#----- Add a header to time_info file (assumes one does not exist)
	if [ ! -f "$i/gen_info_w_header.tsv" ]; then
		echo "Making $i/gen_info_w_header.tsv"
    	gen_info_file=$(find "$i" -iname 'gen_info.tsv')
		echo 'gen' > "$i/gen_info_w_header.tsv"
		cat $gen_info_file >> "$i/gen_info_w_header.tsv"
    else
    	echo "$i/gen_info_w_header.tsv file already exists"
	fi

	#----- Add a header to time_info file (assumes one does not exist):
	if [ ! -f "$i/time_info_w_header.tsv" ]; then
		echo "Making $i/time_info_w_header.tsv"
    	time_info_file=$(find "$i" -iname 'time_info.tsv')
		echo 'time' > "$i/time_info_w_header.tsv" 
		cat $time_info_file >> "$i/time_info_w_header.tsv"
    else
    	echo "$i/time_info_w_header.tsv file already exists"
	fi

	#----- Transpose RNA ID File so it can be merged with the counts later:
	if [ ! -f "$i/transposed_ids_rnas.tsv" ]; then
		echo "Making $i/transposed_ids_rnas.tsv"
    	cut -f1 "$i/ids_rnas.tsv" | paste -s - > "$i/transposed_ids_rnas.tsv"
    else
    	echo "$i/transposed_ids_rnas.tsv already exists"
	fi

	#----- Transpose Protein ID File so it can be merged with the counts later:
	if [ ! -f "$i/transposed_ids_proteins.tsv" ]; then
		echo "Making $i/transposed_ids_proteins.tsv"
    	cut -f1 "$i/ids_proteins.tsv" | paste -s - > "$i/transposed_ids_proteins.tsv"
    else
    	echo "transposed_ids_proteins.tsv already exists"
	fi

	#----- Merge RNA IDs to RNA Counts
	merged_rna_ids_output="$i/rna_ids_and_counts.tsv"
	if [ ! -f $merged_rna_ids_output ]; then
		echo "Making $merged_rna_ids_output"
    	rna_id_file=$(find "$i" -iname 'transposed_ids_rnas.tsv')
		merged_rna_counts_file="$i/merged_rna_counts.tsv"
		cat $rna_id_file $merged_rna_counts_file >> $merged_rna_ids_output
    else
    	echo "$merged_rna_ids_output file already exists"
	fi

	#----- Merge Protein IDs to Protein Counts
	output_file_protein_id="$i/protein_ids_and_counts.tsv"
	if [ ! -f $output_file_protein_id ]; then
		echo "Making $output_file_protein_id"
    	protein_id_file=$(find "$i" -iname 'transposed_ids_proteins.tsv')
		merged_protein_counts_file="$i/merged_protein_counts.tsv"
		cat $protein_id_file $merged_protein_counts_file >> $output_file_protein_id
    else
    	echo "$output_file_protein_id file already exists"
	fi

	base=`basename $i`
	output_multi_gen_rna_protein_counts_ids="$save_dir${base}_multi_gen_rna_protein_counts_ids.tsv"
	if [ ! -f $output_multi_gen_rna_protein_counts_ids ]; then
		echo "Making $output_multi_gen_rna_protein_counts_ids"
		paste "$i/time_info_w_header.tsv" "$i/gen_info_w_header.tsv" "$i/protein_ids_and_counts.tsv" "$i/rna_ids_and_counts.tsv" >> $output_multi_gen_rna_protein_counts_ids 
    else
    	echo "$output_multi_gen_rna_protein_counts_ids file already exists"
	fi

	#---- Use for checking file dimensions
	
	: '
	echo "rows and columns of merged rna counts id file"
	echo $(cat $"$i/merged_rna_ids.tsv" | wc -l)
	echo $(head -1 $"$i/merged_rna_ids.tsv" | sed 's/,/\t/g' | wc -w)
	'

done;
