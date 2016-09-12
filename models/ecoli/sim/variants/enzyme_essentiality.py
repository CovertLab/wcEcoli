
CONTROL_OUTPUT = dict(
	shortName = "control",
	desc = "Control simulation"
	)

def enzymeEssentialityTotalIndices(sim_data):
	nEnzymes = len(sim_data.process.metabolism.rescueEnzymes)
	nConditions = nEnzymes + 1
	return nConditions


def enzymeEssentiality(sim_data, index):
	# Knocks-out genes in order

	nConditions = enzymeEssentialityTotalIndices(sim_data)

	if index % nConditions == 0:
		return CONTROL_OUTPUT, sim_data

	enzymeIndex = (index - 1) % nConditions

	enzymeID = sorted(sim_data.process.metabolism.rescueEnzymes)[enzymeIndex]

	sim_data.process.metabolism.rescueEnzymes.remove(enzymeID)

	sim_data.process.metabolism.knockoutEnzymes.add(enzymeID)

	return dict(
		shortName = "{}_not_rescued".format(geneID),
		desc = "No enzyme rescue of {}.".format(geneID)
		), sim_data
