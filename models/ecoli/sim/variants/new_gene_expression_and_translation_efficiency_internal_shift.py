"""
Compare the impacts of increasing expression level and
increasing or decreasing the translational efficiency of new genes, when
induced after a particular generation. This variant assumes that new genes are
very minimally expressed in the wildtype simulation (i.e. does not explicitly
knockout new genes).

Modifies (after shift):
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

# Generation to induce new gene expression
NEW_GENE_INDUCTION_GEN = 8
NEW_GENE_KNOCKOUT_GEN = 128

def get_new_gene_expression_factor_and_translation_efficiency(sim_data, index):
	"""
	Maps variant index to new gene expression factor and translation effieincy

	Returns:
		expression_factor: will multiply new gene expression by
			10^(expression_factor-1)
		trl_eff_value: translation efficiency to use for the new genes
	"""
	# Determine factor for new gene expression and value for new
	# gene translation efficiency
	trl_eff_data = sim_data.process.translation.translation_efficiencies_by_monomer
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

	return expression_factor, trl_eff_value

def induce_new_genes(sim_data, index):
	"""
	'Induces' new genes by setting their expression levels and translation
	effiencies to the values specified using the variant index.
	"""
	# Map variant index to expression factor and tranlsation efficiency value
	expression_factor, trl_eff_value = get_new_gene_expression_factor_and_translation_efficiency(
		sim_data, index)

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

	# Modify expression and translation efficiency for new genes
	for i in range(len(new_gene_indices)):
		gene_index = new_gene_indices[i]
		monomer_index = new_gene_monomer_indices[i]

		sim_data.adjust_final_expression([gene_index], [expression_factor])
		sim_data.process.translation.translation_efficiencies_by_monomer[
			monomer_index] = trl_eff_value

def new_gene_expression_and_translation_efficiency_internal_shift(sim_data, index):
	"""
	Apply variant. Specifies that at the beginning of NEW_GENE_INDUCTION_GEN,
	the new gene expression and translation efficiency values will be set to
	a specific parameter combination determined by the variant index.
	"""
	# Map variant index to expression factor and tranlsation efficiency value
	expression_factor, trl_eff_value = get_new_gene_expression_factor_and_translation_efficiency(
		sim_data, index)

	# Initialize internal shift dictionary
	setattr(sim_data, 'internal_shift_dict', {})

	# Add the new gene induction to the internal_shift instructions
	for gen in range(NEW_GENE_INDUCTION_GEN, NEW_GENE_KNOCKOUT_GEN):
		sim_data.internal_shift_dict[gen] = [(induce_new_genes, index)]

	# Variant descriptions to save to metadata
	if index == 0:
		return CONTROL_OUTPUT, sim_data

	return dict(
		shortName = "{}_NGEXP_TRLEFF_INTERNAL_SHIFT",
		desc = "Changed expression of new genes by variant index {}.".format(
			expression_factor) + ". Set translation efficiency of new genes "
								 "to {}".format(trl_eff_value)
		), sim_data