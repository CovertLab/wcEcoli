"""
Compare the impacts of increasing expression level and
increasing or decreasing the translational efficiency of new genes.

Modifies:
	sim_data.process.transcription.rna_synth_prob
	sim_data.process.transcription.rna_expression
	sim_data.process.transcription.exp_free
	sim_data.process.transcription.exp_ppgpp
	sim_data.process.transcription.attenuation_basal_prob_adjustments
	sim_data.process.transcription_regulation.basal_prob
	sim_data.process.transcription_regulation.delta_prob
	sim_data.process.translation.translation_efficiencies_by_monomer

Expected variant indices (int, positive):
	0: control (knockout new gene expression)
	z > 0: converted to an index for a list of
		new gene expression variant factor and an index for a list of new gene
		translation efficiency values
		separator = number of translation efficiency values to try
		expression index = index div separator + 1
		translation efficiency index = index mod separator

New gene expression factor:
	x > 0: multiply new gene expression by a factor of 10^(x-1)

New gene translation efficiency value (normalized):
	y > 0: set new gene translation efficiency to y
"""

CONTROL_OUTPUT = dict(
	shortName = "control",
	desc = "Control simulation"
	)

# The variant index will be split into an index for each of these lists
# which are written to simulation metadata for later use in analysis scripts

NEW_GENE_EXPRESSION_FACTORS = [0, 7, 8, 9, 10, 11, 12]
NEW_GENE_TRANSLATION_EFFICIENCY_VALUES = [25, 10, 5, 2.5, 1, 0.5, 0.1, 0.05, 0.01, 0]

SEPARATOR = len(NEW_GENE_TRANSLATION_EFFICIENCY_VALUES)

assert NEW_GENE_EXPRESSION_FACTORS[0] == 0, \
	"The first new gene expression factor should always be the control sim"

def new_gene_expression_and_translation_efficiency(sim_data, index):
	# Determine ids and indices of new genes
	mRNA_sim_data = sim_data.process.transcription.cistron_data.struct_array
	monomer_sim_data = sim_data.process.translation.monomer_data.struct_array
	new_gene_mRNA_ids = mRNA_sim_data[mRNA_sim_data['is_new_gene']]['id'].tolist()
	mRNA_monomer_id_dict = dict(
		zip(monomer_sim_data['cistron_id'], monomer_sim_data['id']))
	new_gene_monomer_ids = [
		mRNA_monomer_id_dict.get(mRNA_id) for mRNA_id in new_gene_mRNA_ids]
	if len(new_gene_mRNA_ids) == 0:
		raise Exception("This variant  is intended to be run on simulations "
			"where the new gene option was enabled, but no new gene mRNAs were "
			"found.")
	if len(new_gene_monomer_ids) == 0:
		raise Exception("This variant is intended to be run on simulations where"
			" the new gene option was enabled, but no new gene proteins "
			"were found.")
	assert len(new_gene_monomer_ids) == len(new_gene_mRNA_ids), \
		'number of new gene monomers and mRNAs should be equal'
	rna_data = sim_data.process.transcription.rna_data
	mRNA_idx_dict = {rna[:-3]: i for i, rna in enumerate(rna_data['id'])}
	new_gene_indices = [
		mRNA_idx_dict.get(mRNA_id) for mRNA_id in new_gene_mRNA_ids]
	monomer_idx_dict = {
		monomer: i for i, monomer in enumerate(monomer_sim_data['id'])}
	new_gene_monomer_indices = [
		monomer_idx_dict.get(monomer_id) for monomer_id in new_gene_monomer_ids]
	trl_eff_data = sim_data.process.translation.translation_efficiencies_by_monomer

	# Determine factor for new gene expression and value for new
	# gene translation efficiency
	if index == 0:
		expression_factor = NEW_GENE_EXPRESSION_FACTORS[0]
		# Note: this value should not matter since gene is knocked out
		trl_eff_value = min(trl_eff_data)
	else:
		trl_eff_list_index = index % SEPARATOR
		if trl_eff_list_index == 0:
			expression_list_index = index // SEPARATOR
		else:
			expression_list_index = index // SEPARATOR + 1

		expression_factor = (
			10**(NEW_GENE_EXPRESSION_FACTORS[expression_list_index] - 1))
		trl_eff_value = NEW_GENE_TRANSLATION_EFFICIENCY_VALUES[trl_eff_list_index]

	# Modify expression and translation efficiency for new genes
	for i in range(len(new_gene_indices)):
		gene_index = new_gene_indices[i]
		monomer_index = new_gene_monomer_indices[i]

		sim_data.adjust_final_expression([gene_index], [expression_factor])
		sim_data.process.translation.translation_efficiencies_by_monomer[
			monomer_index] = trl_eff_value

	if index == 0:
		return CONTROL_OUTPUT, sim_data

	return dict(
		shortName = "{}_NGEXP_TRLEFF",
		desc = "Changed expression of new genes by variant index {}.".format(
			expression_factor) + ". Set translation efficiency of new genes "
								 "to {}".format(trl_eff_value)
		), sim_data
