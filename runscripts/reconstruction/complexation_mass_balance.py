"""
Util function to check that all complexation reactions are mass balanced.


Description:

Prints a vector of reactionIDs where the mass of a reaction is more imbalanced
than the 'numerical_zero' constant defined in this program.

Attempts to fix any entries in proteinComplexes.tsv which are not balanced
by ensuring their metabolite mass (the eigth entry in the 'mw' vector) matches
the mass of the corresponding metabolite in metabolites.tsv. This will often
not work if there should have been multiple metabolites combining into a
complex. Run the script TWICE and manually check any reactions not working
the second time.

The attempt to fix metabolite masses will output into a new file,
proteinComplexes_new.tsv. You must manually RENAME THIS FILE to
proteinComplexes.tsv (overwriting the original) if you want to accept the changes.

Use the unify_protein_complexes_large_from_small.py script (currently in 
wcEcoli/runscripts/reconstruction/unify_protein_complexes_large_from_small.py)
to copy any rows in proteinComplexes.tsv into proteinComplexes_large.tsv if
they differ, joining on the "name" column.


Usage: 

python runscripts/reconstruction/complexation_mass_balance.py

(rename proteinComplexes_new.tsv to proteinComplexes.tsv)

python runscripts/reconstruction/complexation_mass_balance.py

(manually change any reactions/rows in proteinComplexes.tsv) which still error

(optional: to update proteinComplexes_large.tsv):
python runscripts/reconstruction/unify_protein_complexes_large_from_small.py



@author: Morgan Paull
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 07/13/2015
"""

import numpy as np
import ast

import reconstruction.ecoli.knowledge_base_raw
kbr = reconstruction.ecoli.knowledge_base_raw.KnowledgeBaseEcoli()
from wholecell.utils import units
import reconstruction.ecoli.simulation_data

# Constants
numerical_zero = 1e-9

# Load the knowledge base
kb = reconstruction.ecoli.simulation_data.SimulationDataEcoli()
kb.initialize(doubling_time = 60. * units.min, raw_data = kbr, media_conditions = "M9 Glucose minus AAs")

# Find the reactions which are not mass balanced
mba = np.array(kb.process.complexation.massBalance())
unbalancedIdxs = np.where(np.abs(mba) > numerical_zero)[0]

metaboliteWeightMap = {}

# Build dictionary of metabolite ID to molecular weight
for metabolite in kbr.metabolites:
	metaboliteId = metabolite["id"]
	metaboliteWeightMap[metaboliteId] = metabolite["mw7.2"]

# Find the proper weight of the metabolites found to be unbalanced
smallWeightMap = {}
reactionNames = []

# Build a dictionary from the name of a reaction which needs rebalancing
# to the new metabolite weight
for idx in unbalancedIdxs:
	reaction = kbr.complexationReactions[idx]
	reactionName = reaction["id"]
	reactionNames.append(reactionName)
	for molecule in reaction["stoichiometry"]:
		if molecule["type"] == 'metabolite':
			moleculeName = molecule["molecule"]
			smallWeightMap[reactionName] = metaboliteWeightMap[moleculeName]

# Print the indexes and names of reactions which need changing
print unbalancedIdxs
print reactionNames


# If no reactions were unbalanced, report this and stop
if len(smallWeightMap) < 1:
	print("No modifications needed.")
else:
	update_rows = []
	# Loop over the protein complexes file and find the rows in which a metabolite
	#  weight needs to be changed.
	for proteinComplex in kbr.proteinComplexes:
		row = []
		name = proteinComplex["reactionId"]
		molecularWeight = proteinComplex["mw"][7]
		if name in reactionNames:
			if molecularWeight != smallWeightMap[name]:
				# Print out a table of complex name and new metabolite weight,
				# allowing manual changes if desired
				print(name + "\t\t\t" + str(smallWeightMap[name]))
				row = [name,smallWeightMap[name]]
				update_rows.append(row)



	# Attempt to correct the metabolite massses in proteinComplexes.tsv, 
	#  writing into proteinComplexes_new.tsv
	with open('/home/mpaull/wcEcoli/reconstruction/ecoli/flat/proteinComplexes_new.tsv', 'w') as output:

		with open('/home/mpaull/wcEcoli/reconstruction/ecoli/flat/proteinComplexes.tsv', 'rw') as f:
				num_modified = 0

				# Loop through every row in the file
				for row in f:
					columns = row.split('\"')

					# Write the header row completely unmodified
					if columns[1] == 'name':
						output.write(row)
						continue

					# Read the reaction ID
					reactionId = columns[7]

					# Extract the string list of masses, trim junk, and convert it to an actual list
					masses = ast.literal_eval(columns[4][:-2][1:])

					if reactionId in smallWeightMap:
						masses[7] = smallWeightMap[reactionId]

					# Put back the junk around the vector that was originally there
					columns[4] = '\t' + str(masses) + '\t['

					new_row = '\"'.join(columns)

					# Note that a row was modified
					if row != new_row: 
						num_modified += 1

					# Write the row to output whether it was changed or not
	 				output.write(new_row)

	# Print the number of rows modified
	print("Modified %i rows" % num_modified)


