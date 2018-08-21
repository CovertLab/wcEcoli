'''
Print id, mw7.2, location for charged tRNAs to add to modifiedForms.tsv
Output file is charged_data.tsv in the same directory as the script

Notes:
- Masses of small molecules are added to the tRNA mass so that the charging
reactions are mass balanced
- If sim_data.cp file is in the same directory, the script will load sim_data
from there instead of running fitSimData_1
'''

import cPickle
import csv
import os

import numpy as np

from reconstruction.ecoli.fit_sim_data_1 import fitSimData_1
from reconstruction.ecoli.knowledge_base_raw import KnowledgeBaseEcoli
from reconstruction.spreadsheets import JsonWriter

# file paths
file_loc = os.path.dirname(__file__)
sim_data_filename = os.path.join(file_loc, "sim_data.cp")
output_filename = os.path.join(file_loc, 'charged_data.tsv')

# suppress scientific notation output
np.set_printoptions(suppress = True)

# get raw and sim data
raw_data = KnowledgeBaseEcoli()

if os.path.exists(sim_data_filename):
	with open(sim_data_filename, 'rb') as f:
		sim_data = cPickle.load(f)
else:
	sim_data = fitSimData_1(raw_data)

# determine masses and write to output file
with open(output_filename, 'w') as out:
	writer = JsonWriter(out, ["id", "mw7.2", "location"], dialect = "excel-tab")

	trnas = sim_data.process.transcription.rnaData['id'][sim_data.process.transcription.rnaData['isTRna']]
	charged = [x['modifiedForms'] for x in raw_data.rnas if x['id']+'[c]' in trnas]
	filtered_charged = []
	for c1 in charged:
		for c2 in c1:
			if 'FMET' in c2 or 'modified' in c2:
				continue
			filtered_charged += [c2 + '[c]']

	mol_names = sim_data.internal_state.bulkMolecules.bulkData['id']
	mws = sim_data.internal_state.bulkMolecules.bulkData['mass']
	for rxn in raw_data.modificationReactions:
		reactants = []
		products = []
		for mol in rxn['stoichiometry']:
			if mol['coeff'] == -1:
				reactants += ['%s[%s]' % (mol['molecule'], mol['location'])]
			else:
				products += ['%s[%s]' % (mol['molecule'], mol['location'])]

		for trna, ctrna in zip(trnas, filtered_charged):
			if trna in reactants and ctrna in products:
				mass = 0
				for reactant in reactants:
					if reactant in mol_names:
						mass += mws[np.where(mol_names == reactant)[0][0]].asNumber()
					else:
						print 'could not get mass for %s' % (reactant)
				for product in products:
					if product == ctrna:
						continue

					if product in mol_names:
						mass -= mws[np.where(mol_names == product)[0][0]].asNumber()
					else:
						print 'could not get mass for %s' % (product)

				writer.writerow({
					'id': ctrna[:-3],
					'mw7.2': mass,
					'location': [ctrna[-2]],
					})
