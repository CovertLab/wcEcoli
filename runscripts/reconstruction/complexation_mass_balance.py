import numpy as np
import ast

import reconstruction.ecoli.knowledge_base_raw
kbr = reconstruction.ecoli.knowledge_base_raw.KnowledgeBaseEcoli()
from wholecell.utils import units
import reconstruction.ecoli.simulation_data
kb = reconstruction.ecoli.simulation_data.SimulationDataEcoli()
kb.initialize(doubling_time = 60. * units.min, raw_data = kbr, media_conditions = "M9 Glucose minus AAs")

mba = np.array(kb.process.complexation.massBalance())
unbalancedIdxs = np.where(np.abs(mba) > 1e-9)[0]

metaboliteWeightMap = {}

# Build dictionary of metabolite ID to molecular weight
for metabolite in kbr.metabolites:
	metaboliteId = metabolite["id"]
	metaboliteWeightMap[metaboliteId] = metabolite["mw7.2"]


smallWeightMap = {}
reactionNames = []
# Find the proper weight of the metabolites found to be unbalanced

print unbalancedIdxs

for idx in unbalancedIdxs:
	reaction = kbr.complexationReactions[idx]
	reactionName = reaction["id"]
	reactionNames.append(reactionName)
	for molecule in reaction["stoichiometry"]:
		if molecule["type"] == 'metabolite':
			moleculeName = molecule["molecule"]
			smallWeightMap[reactionName] = metaboliteWeightMap[moleculeName]


if len(smallWeightMap) < 1:
	print("No modifications needed.")
else:
	update_rows = []
	# Loop over the protein complexes file and find the rows in which a metabolite
	#  weight needs to be changed
	for proteinComplex in kbr.proteinComplexes:
		row = []
		name = proteinComplex["reactionId"]
		molecularWeight = proteinComplex["mw"][7]
		if name in reactionNames:
			if molecularWeight != smallWeightMap[name]:
				print(name + "\t\t\t" + str(smallWeightMap[name]))
				row = [name,smallWeightMap[name]]
				update_rows.append(row)



	with open('/home/mpaull/wcEcoli/reconstruction/ecoli/flat/proteinComplexes_new.tsv', 'w') as output:

		with open('/home/mpaull/wcEcoli/reconstruction/ecoli/flat/proteinComplexes.tsv', 'rw') as f:
				num_modified = 0

				for row in f:
					columns = row.split('\"')

					# Write the header row completely unmodified
					if columns[1] == 'name':
						output.write(row)
						continue

					# Read the reaction ID
					reactionId = columns[7]

					# import ipdb; ipdb.set_trace()

					# Extract the string list of masses and convert it to an actual list
					masses = ast.literal_eval(columns[4][:-2][1:])

					if reactionId in smallWeightMap:
						# print('Original is: %s' % str(masses[7]))
						# print('New is: %s' % str(smallWeightMap[reactionId]))
						masses[7] = smallWeightMap[reactionId]


					columns[4] = '\t' + str(masses) + '\t['

					new_row = '\"'.join(columns)

					if row != new_row: 
						num_modified += 1

	 				output.write(new_row)

	print("Modified %i rows" % num_modified)
				# print(num_modified)


