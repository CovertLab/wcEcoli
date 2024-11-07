"""
Compare the impacts of increasing expression level and
increasing or decreasing the translational efficiency of new genes, with the
option to induce and/or knock out new gene expression after a particular
generation. This variant assumes that new genes are present but not expressed in
the wildtype simulation (new genes are knocked out by default).

- From generations [0, `NEW_GENE_INDUCTION_GEN`), the new genes will be
	transcribed and translated based upon the default parameters in wildtype
	simulation (i.e. from the flat files).
- From generations [`NEW_GENE_INDUCTION_GEN`, `NEW_GENE_KNOCKOUT_GEN`),
	the new genes will be transcribed and translated using the expression
	factor and translation efficiency value from the variant index.
- From generations [`NEW_GENE_KNOCKOUT_GEN`, `FINAL_SHIFT_GEN`), new
	gene expression probabilities will be set to 0 corresponding to new
	gene knockout.

Modifies:
	sim_data.condition
	sim_data.external_state.current_timeline_id

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
	z > 0: converted to an index for a media condition, an index for a list of
		new gene expression variant factors, and an index for a list of new gene
		translation efficiency values

		media condition index = index div 1000
		index remainder = index - media condition index

		separator = number of translation efficiency values to try
		expression index = remainder index div separator + 1
		translation efficiency index = remainder index mod separator

New gene expression factor:
	x > 0: multiply new gene expression by a factor of 10^(x-1)

New gene translation efficiency value (normalized):
	y > 0: set new gene translation efficiency to y
"""

CONTROL_OUTPUT = dict(
	shortName = "control",
	desc = "Control simulation"
	)

# NOTE: The version of these global variables on the master branch will
# be used for Jenkins testing

# NOTE: If these values are greater than the number of generations you are
# running, you will not see their effects.
# OPTION: Set = -1 to skip induction or knockout shifts.
NEW_GENE_INDUCTION_GEN = 2 # Generation index to induce new gene expression
NEW_GENE_KNOCKOUT_GEN = -1 # Generation index to knock out new gene expression
assert NEW_GENE_INDUCTION_GEN != 0, (
	"New genes must be induced after the first generation.")
if NEW_GENE_KNOCKOUT_GEN != -1:
	assert NEW_GENE_KNOCKOUT_GEN > NEW_GENE_INDUCTION_GEN, (
		"New genes are knocked out by default, so induction should happen"
		" before knockout.")

# The variant index will be split into an index for each of these lists
# which are written to simulation metadata for later use in analysis scripts

NEW_GENE_EXPRESSION_FACTORS = [0, 7]
NEW_GENE_TRANSLATION_EFFICIENCY_VALUES = [1]

SEPARATOR = len(NEW_GENE_TRANSLATION_EFFICIENCY_VALUES)
assert NEW_GENE_EXPRESSION_FACTORS[0] == 0, (
	"The first new gene expression factor should always be the control sim")


def condition(sim_data, condition_index):
	"""
	Condition variant for simulations in different environmental conditions

	Modifies:
		sim_data.condition
		sim_data.external_state.current_timeline_id

	Expected variant indices (dependent on sim_data.ordered_conditions and should
	be the same order as rows in condition_defs.tsv):
		0: control (minimal media)
		1: with amino acids
		2: acetate
		3: succinate
		4: minimal media (anaerobic)
	"""
	condition_labels = sim_data.ordered_conditions
	condition_label = condition_labels[condition_index]
	sim_data.condition = condition_label
	sim_data.external_state.current_timeline_id = condition_label
	sim_data.external_state.saved_timelines[condition_label] = [(
		0, sim_data.conditions[condition_label]["nutrients"])]


def determine_new_gene_ids_and_indices(sim_data):
	"""
	Determines the ids and indices of new gene mRNAs and proteins using the new
	gene flag in sim_data.

	Returns:
		new_gene_mRNA_ids: ids corresponding to new gene mRNAs
		new_gene_mRNA_indices: indexes in rna_data table corresponding to new
			gene mRNAs
		new_gene_monomer_ids: ids corresponding to new gene monomers
		new_gene_monomer_indices: indexes in monomer_data table corresponding to
			new gene monomers
	"""
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

	return new_gene_mRNA_ids, new_gene_indices, new_gene_monomer_ids, new_gene_monomer_indices


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
	# Map variant index to expression factor and translation efficiency value
	expression_factor, trl_eff_value = get_new_gene_expression_factor_and_translation_efficiency(
		sim_data, index)

	# Determine ids and indices of new genes
	new_gene_mRNA_ids, new_gene_indices, new_gene_monomer_ids, \
		new_gene_monomer_indices = determine_new_gene_ids_and_indices(
		sim_data)

	# Modify expression and translation efficiency for new genes
	for i in range(len(new_gene_indices)):
		gene_index = new_gene_indices[i]
		monomer_index = new_gene_monomer_indices[i]

		sim_data.adjust_new_gene_final_expression([gene_index], [expression_factor])
		sim_data.process.translation.translation_efficiencies_by_monomer[
			monomer_index] = trl_eff_value


def knockout_induced_new_gene_expression(sim_data, index):
	"""
	Intended to be used after induce_new_gene_expression has been used on
	previous generations. Knocks out new genes by setting their expression
	levels to 0. Maintains translation effiencies at the values specified using
	the variant index so that way any pre-existing new gene mRNAs have the same
	translation efficiency as when they were created.
	"""
	# Map variant index to expression factor and tranlsation efficiency value
	expression_factor, trl_eff_value = get_new_gene_expression_factor_and_translation_efficiency(
		sim_data, index)

	# Determine ids and indices of new genes
	new_gene_mRNA_ids, new_gene_indices, new_gene_monomer_ids, \
		new_gene_monomer_indices = determine_new_gene_ids_and_indices(
		sim_data)

	# Knockout expression and maintain translation efficiency for new genes
	for i in range(len(new_gene_indices)):
		gene_index = new_gene_indices[i]
		monomer_index = new_gene_monomer_indices[i]

		sim_data.adjust_final_expression([gene_index], [0])
		sim_data.process.translation.translation_efficiencies_by_monomer[
			monomer_index] = trl_eff_value


def new_gene_internal_shift(sim_data, index):
	"""
	Apply variant. Specifies that from NEW_GENE_INDUCTION_GEN to
	NEW_GENE_KNOCKOUT_GEN, the new gene expression and translation efficiency
	values will be set to a specific parameter combination determined by the
	variant index. Then, from NEW_GENE_KNOCKOUT_GEN to FINAL_SHIFT_GEN, the new
	gene expression values will be set to 0, while translation efficiency will
	maintain the values from the variant index, so that already created new gene
	mRNAs have the same translation efficiency as they did when they were
	transcribed.
	"""
	# Set media condition
	condition_index = index // 1000
	condition(sim_data, condition_index)

	# Map variant index to expression factor and tranlsation efficiency value
	index_remainder = index - condition_index * 1000
	expression_factor, trl_eff_value = get_new_gene_expression_factor_and_translation_efficiency(
		sim_data, index_remainder)

	# Initialize internal shift dictionary
	setattr(sim_data, 'internal_shift_dict', {})

	# Add the new gene induction to the internal_shift instructions
	# Note: Must add an entry for each non wildtype gen because sim_data is
	# reloaded from the file between generations
	if NEW_GENE_INDUCTION_GEN != -1:
		sim_data.internal_shift_dict[NEW_GENE_INDUCTION_GEN] = [
			(induce_new_genes, index_remainder)]
	if NEW_GENE_KNOCKOUT_GEN != -1:
		sim_data.internal_shift_dict[NEW_GENE_KNOCKOUT_GEN] = [
			(knockout_induced_new_gene_expression, index_remainder)]

	# Variant descriptions to save to metadata
	if index == 0:
		return CONTROL_OUTPUT, sim_data

	return dict(
		shortName = "{}_NGEXP_TRLEFF_INTERNAL_SHIFT."
			+ "{}_env".format(sim_data.ordered_conditions[condition_index]),
		desc = "Changed expression of new genes by variant index {}.".format(
			expression_factor) + ". Set translation efficiency of new genes "
			"to {}".format(trl_eff_value)
			+ "Simulation of condition {}.".format(sim_data.ordered_conditions[condition_index])
		), sim_data
