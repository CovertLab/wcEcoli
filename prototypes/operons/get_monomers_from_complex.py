from reconstruction.spreadsheets import JsonReader, JsonWriter
import numpy as np
import copy

def monomers_to_complexes():
	# extract each line of the complexes tsv
	tsv_file = '/Users/taryn/GoogleDrive/code/wcEcoli/reconstruction/ecoli/flat/complexationReactions.tsv'
	tsv_list = []
	with open(tsv_file) as tsvfile:
		reader = JsonReader((row for row in tsvfile if not row.startswith('#')), dialect="excel-tab")
		fieldnames = reader.fieldnames
		for row in reader:
			tsv_list.append(row)
			# import ipdb; ipdb.set_trace()

	# get all complex ids
	# '[' + str(x['stoichiometry'][0]['location']) + ']'
	complex_dict = {x['stoichiometry'][0]['molecule'] + '[' + str(x['stoichiometry'][0]['location']) + ']'\
				 : x['stoichiometry'][1:] for x in tsv_list}

	# CREATE MAPPING OF COMPLEXES TO MONOMERS
	complex_to_monomers = {}
	for c in complex_dict.keys():
		complex_to_monomers[c] = {}
		subunit_list = copy.deepcopy(complex_dict[c])
		subunit_multipliers = {}
		i = 0

		# make a list of the components of the complex, check if they are monomers
		# or complexes. If complex, add the components of that complex to the list
		# and continue to iterate until you have fully decomposed to monomers
		while len(subunit_list) > 0:
			for subunit in subunit_list:
				# if this subunit is a monomer, add it to dict for the complex name

				if subunit['molecule'] not in complex_dict.keys():
					#check if monomer is already in complex dict
					if subunit['molecule'] not in complex_to_monomers[c]:
						complex_to_monomers[c][subunit['molecule']] = 0

					# check if multiplier exits
					if subunit['molecule'] in subunit_multipliers:
						su_ratio = subunit_multipliers[subunit['molecule']]
					else:
						su_ratio = 1.0

					# increase count by the number used in the complex 
					complex_to_monomers[c][subunit['molecule']] += -subunit['coeff'] * su_ratio
					subunit_list.remove(subunit)

				# if this subunit is a complex, add its components to the list
				elif subunit['molecule'] in complex_dict.keys():
					# append subunits of the subcomplex and delete the subcomplex from the list
					subunit_list.extend(complex_dict[subunit['molecule']])
					# keep stoichiometry for the subcomplex to multiply the monomers by
					for su in complex_dict[subunit['molecule']]:
						subunit_multipliers[su['molecule']] = -subunit['coeff']

					subunit_list.remove(subunit)


	# CREATE MAPPING OF MONOMER TO COMPLEX
	monomer_to_complex = {}
	for c in complex_to_monomers:
		for m in complex_to_monomers[c]:
			# add key to dict if it doesn't exist
			if m not in monomer_to_complex:
				monomer_to_complex[m] = {}
			# store the amount monomer used in each complex
			monomer_to_complex[m][c] = complex_to_monomers[c][m]


	return monomer_to_complex








