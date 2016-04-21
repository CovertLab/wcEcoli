from wholecell.utils import units
from reconstruction.ecoli.knowledge_base_raw import KnowledgeBaseEcoli
from reconstruction.ecoli.fit_sim_data_1 import fitSimData_1

CONTROL_OUTPUT = dict(
	shortName = "control",
	desc = "Control simulation"
	)

def uptakeRateTotalIndices(sim_data):
	nEnvironments = len(sim_data.envDict)
	return nEnvironments


def uptakeRate(sim_data, index):
	# if index == 0:
	# 	return CONTROL_OUTPUT, sim_data

	environments = {
		1: "000004_unbound_glucose", 
		2: "000005_unbound_glucose_aa",
		3: "000006_unbound_glycerol",
		4: "000007_unbound_glycerol_aa",
		5: "000008_unbound_succinate",
		6: "000009_unbound_succinate_aa",
	}

	growthRates = units.min * [44., 25., 59., 33., 80., 37.]
	
	raw_data = KnowledgeBaseEcoli()
	sim_data = fitSimData_1(raw_data, doubling_time = growthRates[index-1], env = environments[index])

	return dict(
		shortName = "{}_uptake".format(environments[index]),
		desc = "Simulation of uptake rate {}.".format(environments[index])
		), sim_data
