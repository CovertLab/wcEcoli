
CONTROL_OUTPUT = dict(
	shortName = "control",
	desc = "Control simulation"
	)

def geneKnockoutTotalIndices(sim_data):
	nGenes = sim_data.process.transcription.rnaData.fullArray().size
	nConditions = nGenes + 1
	return nConditions


def geneKnockout(sim_data, index):
	# Knocks-out genes in order
	nConditions = geneKnockoutTotalIndices(sim_data)

	if index % nConditions == 0:
		return CONTROL_OUTPUT, sim_data

	geneIndex = (index - 1) % nConditions
	geneID = sim_data.process.transcription.rnaData["id"][geneIndex]
	KO_dict = {geneID: 0.0}
	sim_data.genetic_perturbations = KO_dict

    #Test in conditions with amino acids added, per experimental data
	sim_data.nutrientsTimeSeriesLabel = "000003_aa"

	sim_data.process.transcription.rnaExpression[sim_data.condition][geneIndex] = 0.0

	return dict(
		shortName = "{}_KO".format(geneID),
		desc = "Complete knockout of {}.".format(geneID)
		), sim_data

