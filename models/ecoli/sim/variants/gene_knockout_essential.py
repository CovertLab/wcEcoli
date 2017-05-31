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
	#esGenes = validation_data.essentialGenes.essentialGenes
	
	#Reading in list of genes that were not impacted in first generation
	#gene_file = '~/wcEcoli/function_no_impact_geneIDs.txt'
	gene_file = '~/wcEcoli/basal_genes_large_impact.txt'
	prot_file = '~/wcEcoli/validation/ecoli/flat/essentialGenes.tsv'
	df_esGenes = pd.read_table(gene_file, sep='\t')
	df_prot = pd.read_table(prot_file, sep='\t')
	esGenes = df_esGenes['Genes'].tolist()
	proteinIDs = df_prot['proteinID'].tolist()

	nConditions = len(esGenes)

	if index % nConditions == 0:
		#sim_data.nutrientsTimeSeriesLabel = "000003_aa"
		#sim_data.condition = "with_aa"
		return CONTROL_OUTPUT, sim_data

	esGene = esGenes[index - 1] + '_RNA[c]'
	protID = proteinIDs[np.where(df_prot['FrameID'] == esGenes[index - 1])[0]]

	scaling_factor = 0.5
	KO_dict = {esGene: scaling_factor}
	sim_data.genetic_perturbations = KO_dict
	ind = np.where(sim_data.process.transcription.rnaData["id"] == esGene)[0][0]

	sim_data.process.transcription.rnaExpression[sim_data.condition][ind] = 0.0
	normalize(sim_data.process.transcription.rnaExpression[sim_data.condition])

	try:
		enzyme = protID + '[c]'
		rxns = [x for x in sim_data.process.metabolism.reactionCatalysts.keys() if enzyme in sim_data.process.metabolism.reactionCatalysts[x]]
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

