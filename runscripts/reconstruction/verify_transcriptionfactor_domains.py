"""
Util function checks TF subunits in tfSingleComponent.tsv are in the appropriate DNA-binding/non-binding column 


Description:

Checks that transcription factor subunits listed in tfSingleComponent.tsv are
in the appropriate column (first column should be DNA binding, second column
should be non-DNA-binding).

Assesses correctness by examining the transcription factor name in
proteinComplexes.tsv or proteins.tsv (checks for the tf ID in both), and seeing
if that transcription factor name contains "DNA binding" or a similar string.



Usage: 

python runscripts/reconstruction/verify_transcriptionfactor_domains.py


@author: Morgan Paull
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 07/20/2015
"""

import numpy as np
import re

import reconstruction.ecoli.knowledge_base_raw
kbr = reconstruction.ecoli.knowledge_base_raw.KnowledgeBaseEcoli()
from wholecell.utils import units
import reconstruction.ecoli.simulation_data



# DNA-binding IDs 
DNA_binding_list = []
DNA_binding_dict =  {}

# Non-DNA-binding IDs 
non_DNA_binding_list = []
non_DNA_binding_dict = {}


# Loop through the transciption factors and make lists of all members of the
# DNA and non-DNA binding columns. Also make a dict from ID to the line in
# tfSingleComponent.tsv
for i, transciptionFactorInfo in enumerate(kbr.tfSingleComponent):

	# Loop through the DNA-binding frame ids in this transcription factor's info
	for DNA_binding_ID in transciptionFactorInfo["DNA-binding frame ids"]:
		# Make a map from DNA binding domains to their lines in tfSingleComponent.tsv
		DNA_binding_dict[DNA_binding_ID] = i
		# Add all of the DNA binding domains to the list
		DNA_binding_list.append(DNA_binding_ID)

	# Loop through the DNA-binding frame ids in this transcription factor's info
	for DNA_binding_ID in transciptionFactorInfo["Non-DNA-binding frame ids"]:
		# Make a map from DNA binding domains to their lines in tfSingleComponent.tsv
		non_DNA_binding_dict[DNA_binding_ID] = i
		# Add all of the DNA binding domains to the list
		non_DNA_binding_list.append(DNA_binding_ID)


# List of proteins in the DNA-binding column without "DNA-binding" in the name
DNA_column_not_named = []
# List of proteins in the non-DNA-binding column WITH "DNA-binding" in the name
non_DNA_column_named = []


# Loop over all proteinComplexes.tsv
for i, proteinComplex in enumerate(kbr.proteinComplexes):
	# When a protein complex with an ID in the DNA_binding_list appears, check
	# the name of that protein complex for 'DNA-binding'
	if proteinComplex["id"] in DNA_binding_list:
		name = proteinComplex["name"]
		m = re.search('DNA-binding|DNA binding', name)
		if m == None:
			DNA_column_not_named.append(('proteinComplexes.tsv line ' + str(i+2), 'tfSingleComponent.tsv line ' + str(DNA_binding_dict[proteinComplex["id"]] + 2), proteinComplex["id"],name))
	
	# When a protein complex with an ID in the non_DNA_binding_list appears,
	# check that the name of that protein complex doesn't contain 'DNA-binding'
	if proteinComplex["id"] in non_DNA_binding_list:
		name = proteinComplex["name"]
		m = re.search('DNA-binding|DNA binding', name)
		if m != None:
			non_DNA_column_named.append(('proteinComplexes.tsv line ' + str(i+2), 'tfSingleComponent.tsv line ' + str(non_DNA_binding_dict[proteinComplex["id"]] + 2),proteinComplex["id"],name))


# Loop over all proteins.tsv when a protein with an ID in the DNA_binding_list
# appears, check the name of that protein
for i, protein in enumerate(kbr.proteins):
	# When a protein with an ID in the DNA_binding_list appears, check the name
	# of that protein for 'DNA-binding'
	if protein["id"] in DNA_binding_list:
		name = protein["name"]
		m = re.search('DNA-binding|DNA binding', name)
		if m == None:
			DNA_column_not_named.append(('proteins.tsv line ' + str(i+2), 'tfSingleComponent.tsv line ' + str(DNA_binding_dict[protein["id"]] + 2),protein["id"],name))
	
	# When a protein with an ID in the DNA_binding_list appears, check the name
	# of that protein for 'DNA-binding'
	if protein["id"] in non_DNA_binding_list:
		name = protein["name"]
		m = re.search('DNA-binding|DNA binding', name)
		if m != None:
			non_DNA_column_named.append(('proteins.tsv line ' + str(i+2), 'tfSingleComponent.tsv line ' + str(non_DNA_binding_dict[protein["id"]] + 2),protein["id"],name))


print "\n In DNA-binding column and maybe shouldn't be: "
for entry in DNA_column_not_named:
	print entry

print "\n In non-DNA-binding column and maybe shouldn't be: "
for entry in non_DNA_column_named:
	print entry