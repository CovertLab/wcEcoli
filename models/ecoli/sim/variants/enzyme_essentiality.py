knockoutEnzymes = [
"FDNG-MONOMER[p]",
"FORMATEDEHYDROGH-MONOMER[c]",
"G7106-MONOMER[c]",
"MONOMER0-12[i]",
"XANTHOSINEPHOSPHORY-MONOMER[c]",
"ACPSUB-MONOMER[c]",
"ATOA-MONOMER[c]",
"CARBODEHYDRAT-MONOMER[c]",
"DMSA-MONOMER[i]",
"EG11234-MONOMER[c]",
"EG11245-MONOMER[c]",
"EG12279-MONOMER[c]",
"FUCPALDOL-MONOMER[c]",
"G6140-MONOMER[c]",
"G6710-MONOMER[c]",
"G6716-MONOMER[m]",
"G7145-MONOMER[c]",
"GLYOCARBOLIG-MONOMER[c]",
"LYXK-MONOMER[c]",
"MONOMER-162[c]",
"NARK-MONOMER[i]",
"NIKA-MONOMER[i]",
"NIKB-MONOMER[i]",
"NIKC-MONOMER[i]",
"RIBULPEPIM-MONOMER[c]",
"TAGAALDOL1-MONOMER[c]",
"XYLH-MONOMER[i]",
"YCBM-MONOMER[i]",
"YJCW-MONOMER[i]",
"BETAGALACTOSID-MONOMER[c]",
"G6277-MONOMER[c]",
"G6711-MONOMER[c]",
"G7234-MONOMER[c]",
"G7921-MONOMER[o]",
"GALACTOACETYLTRAN-MONOMER[c]",
"KDPB-MONOMER[i]",
"TNAB-MONOMER[i]",
]

CONTROL_OUTPUT = dict(
	shortName = "control",
	desc = "Control simulation"
	)

def enzymeEssentialityTotalIndices(sim_data):
	if len(knockoutEnzymes) > 0:
		nEnzymes = len(knockoutEnzymes)
	else:
		nEnzymes = len(sim_data.process.metabolism.rescueEnzymes)

	nConditions = nEnzymes + 1
	return nConditions


def enzymeEssentiality(sim_data, index):
	# Knocks-out genes in order

	nConditions = enzymeEssentialityTotalIndices(sim_data)

	if index % nConditions == 0:
		return CONTROL_OUTPUT, sim_data

	enzymeIndex = (index - 1) % nConditions

	if len(knockoutEnzymes) > 0:
		enzymeID = knockoutEnzymes[enzymeIndex]
	else:
		enzymeID = sorted(sim_data.process.metabolism.rescueEnzymes)[enzymeIndex]

	sim_data.process.metabolism.rescueEnzymes.remove(enzymeID)

	sim_data.process.metabolism.knockoutEnzymes.add(enzymeID)

	return dict(
		shortName = "{}_not_rescued".format(geneID),
		desc = "No enzyme rescue of {}.".format(geneID)
		), sim_data
