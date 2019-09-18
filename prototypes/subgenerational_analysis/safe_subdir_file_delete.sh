: '
Author: Mialy DeFelice
Date: 091819
Function:
A safer way to delete unwanted files.
'
base_output_dir="/Users/mialydefelice/Documents/code_repositories/wcEcoli/out/counts/wildtype_000000/count_out_2"

for i in $(find "${base_output_dir}"* -depth 1 -type d); do

	if [ -f "$i/merged_rna_counts.tsv" ]; then
		echo "Removing $i/merged_rna_counts.tsv"
		rm "$i/merged_rna_counts.tsv"
	else
		echo "$i/merged_rna_counts.tsv does not exist"
	fi

	#---- To save space, remove some precursor files.
	if [ -f "$i/merged_protein_counts.tsv" ]; then
		echo "Removing $i/merged_protein_counts.tsv"
		rm "$i/merged_protein_counts.tsv"
	else
		echo "$i/merged_protein_counts.tsv does not exist"
	fi

	#---- To save space, remove some precursor files.
	if [ -f "$i/merged_complex_counts.tsv" ]; then
		echo "Removing $i/merged_complex_counts.tsv"
		rm "$i/merged_complex_counts.tsv"
	else
		echo "$i/merged_complex_counts.tsv does not exist"
	fi

done;