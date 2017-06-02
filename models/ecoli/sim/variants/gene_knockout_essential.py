import cPickle
import numpy as np
import pandas as pd
import wholecell.states.bulk_molecules
from wholecell.utils.fitting import normalize

CONTROL_OUTPUT = dict(
	shortName = "control",
	desc = "Control simulation"
	)

def geneKnockoutEssential(sim_data, validation_data, index):
	
	#Read in all essential genes as defined in validation data
	esGenes = validation_data.essentialGenes.essentialGenes
	proteinIDs = validation_data.essentialGenes.essentialProteins

	nConditions = len(esGenes)

	if index % nConditions == 0:
		return CONTROL_OUTPUT, sim_data

	esGene = esGenes[index - 1] + '_RNA[c]'
	protID = proteinIDs[index-1]

	scaling_factor = 0.5
	KO_dict = {esGene: scaling_factor}
	sim_data.genetic_perturbations = KO_dict
	ind = np.where(sim_data.process.transcription.rnaData["id"] == esGene)[0][0]

	sim_data.process.transcription.rnaExpression[sim_data.condition][ind] = 0.0
	normalize(sim_data.process.transcription.rnaExpression[sim_data.condition])

	try:
		rxns = [x for x in sim_data.process.metabolism.reactionCatalysts.keys() if protID in sim_data.process.metabolism.reactionCatalysts[x]]
		RXN_dict = dict()
		for rxn in rxns:
			RXN_dict[rxn] = scaling_factor
		sim_data.rxn_perturbations = RXN_dict
	
	except:
		return dict(
			shortName = "{}_KO".format(esGene),
			desc = "Complete knockout of {}.".format(esGene)
			), sim_data


	return dict(
			shortName = "{}_KO".format(esGene),
			desc = "Complete knockout of {}.".format(esGene)
			), sim_data

