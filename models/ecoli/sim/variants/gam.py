from reconstruction.ecoli.knowledge_base_raw import KnowledgeBaseEcoli
from reconstruction.ecoli.fit_sim_data_1 import fitSimData_1

# fraction of original GAM to set
GAMS = [0.5, 0.75, 0.9, 1.1, 1.25, 1.5, 2]

CONTROL_OUTPUT = dict(
	shortName = "control",
	desc = "Control simulation"
	)

def gamTotalIndices(sim_data):
	nEnvironments = len(GAMS) + 1
	return nEnvironments


def gam(sim_data, index):
	nEnvironments = gamTotalIndices(sim_data)

	if index % nEnvironments == 0:
		return CONTROL_OUTPUT, sim_data

	raw_data = KnowledgeBaseEcoli()
	sim_data = fitSimData_1(raw_data, GAM = GAMS[index - 1])

	return dict(
		shortName = "GAM {}".format(sim_data.constants.growthAssociatedMaintenance),
		desc = "Setting GAM to {}.".format(sim_data.constants.growthAssociatedMaintenance)
		), sim_data
