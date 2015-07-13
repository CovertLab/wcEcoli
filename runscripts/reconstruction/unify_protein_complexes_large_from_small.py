"""
Util function to copy any rows from proteinComplexes.tsv into
proteinComplexes_large.tsv where they differ.

Description:

Copies any rows in proteinComplexes.tsv into proteinComplexes_large.tsv if
they differ, joining on the "name" column. Only touches these two files.

Usage: 

python runscripts/reconstruction/unify_protein_complexes_large_from_small.py


@author: Morgan Paull
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 07/13/2015
"""


with open('/home/mpaull/wcEcoli/reconstruction/ecoli/flat/proteinComplexes.tsv', 'r') as infile_small:

	# Read the proteinComplexes.tsv file (the smaller file)	
	num_modified = 0
	nameMap = {}

	# Make a dict mapping reaction name to the row in proteinComplexes.tsv
	for row in infile_small:
		columns = row.split('\"')

		name = columns[1]

		nameMap[name] = row


with open('/home/mpaull/wcEcoli/reconstruction/ecoli/flat/proteinComplexes_large_new.tsv', 'w') as output:
	with open('/home/mpaull/wcEcoli/reconstruction/ecoli/flat/proteinComplexes_large.tsv', 'r') as infile_large:
		num_modified = 0

		# Read proteinComplexes_large
		for row in infile_large:
			columns = row.split('\"')
			name = columns[1]
			
			# For every line in the file, check if it is the same as the line in proteinComplexes.tsv
			if name in nameMap:
				# If not, change it to that row.
				output_row = nameMap[name]
				num_modified =+ 1
			else:
				output_row = row

			# Write the resulting row, modified or unmodified, to the output file
			output.write(output_row)

		print("Modified %i rows. " % num_modified)