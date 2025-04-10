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
	sim_data.process.transcription.RNA_synth_prob
	sim_data.process.transcription.RNA_expression
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

import numpy as np

CONTROL_OUTPUT = dict(
	shortName = "control",
	desc = "Control simulation"
	)

# NOTE: The version of these global variables on the master branch will
# be used for Jenkins testing

# Renormalization options
RNAP_RENORMALIZATION = False # Preserve expression of RNAP subunits
RIBOSOME_RENORMALIZATION = False # Preserve expression of ribosomal components
NEW_GENE_RENORMALIZATION = True

# NOTE: If these values are greater than the number of generations you are
# running, you will not see their effects.
# OPTION: Set = -1 to skip induction or knockout shifts.
NEW_GENE_INDUCTION_GEN = 8 # Generation index to induce new gene expression
NEW_GENE_KNOCKOUT_GEN = -1 # Generation index to knock out new gene expression
assert NEW_GENE_INDUCTION_GEN != 0, (
	"New genes must be induced after the first generation.")
if NEW_GENE_KNOCKOUT_GEN != -1:
	assert NEW_GENE_KNOCKOUT_GEN > NEW_GENE_INDUCTION_GEN, (
		"New genes are knocked out by default, so induction should happen"
		" before knockout.")

# The variant index will be split into an index for each of these lists
# which are written to simulation metadata for later use in analysis scripts
NEW_GENE_EXPRESSION_FACTORS = [0, 5, 6, 7, 8, 9, 10, 11, 12]
NEW_GENE_TRANSLATION_EFFICIENCY_VALUES = [10, 5, 1, 0.1, 0] # Should be descending

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
		new_gene_mRNA_indices: indexes in RNA_data table corresponding to new
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
		raise Exception("This is intended to be run on simulations "
			"where the new gene option was enabled, but no new gene mRNAs were "
			"found.")
	if len(new_gene_monomer_ids) == 0:
		raise Exception("This is intended to be run on simulations where"
			" the new gene option was enabled, but no new gene proteins "
			"were found.")
	assert len(new_gene_monomer_ids) == len(new_gene_mRNA_ids), \
		'number of new gene monomers and mRNAs should be equal'
	RNA_data = sim_data.process.transcription.rna_data
	mRNA_idx_dict = {RNA[:-3]: i for i, RNA in enumerate(RNA_data['id'])}
	new_gene_indices = [
		mRNA_idx_dict.get(mRNA_id) for mRNA_id in new_gene_mRNA_ids]
	monomer_idx_dict = {
		monomer: i for i, monomer in enumerate(monomer_sim_data['id'])}
	new_gene_monomer_indices = [
		monomer_idx_dict.get(monomer_id) for monomer_id in new_gene_monomer_ids]

	# TODO: fix the above - only works if the new gene is not in a TU with anything else

	return new_gene_mRNA_ids, new_gene_indices, new_gene_monomer_ids, new_gene_monomer_indices


def determine_RNAP_component_ids_and_indices(sim_data):
	"""
	Determines the ids and indices of RNAP components using the RNAP flag in
	sim_data.

	Returns:
		RNAP_component_ids: ids corresponding to RNAP components
		RNAP_component_indices: indexes in monomer_data table corresponding to
			RNAP components
	"""

	monomer_sim_data = sim_data.process.translation.monomer_data.struct_array
	monomer_id_to_cistron_id_dict = dict(
		zip(monomer_sim_data['id'], monomer_sim_data['cistron_id']))
	RNA_data = sim_data.process.transcription.rna_data
	RNA_idx_to_id_dict = {i: RNA[:-3] for i, RNA in enumerate(RNA_data['id'])}
	monomer_idx_dict = {
		monomer: i for i, monomer in enumerate(monomer_sim_data['id'])}

	RNAP_component_monomer_ids = sim_data.molecule_groups.RNAP_subunits
	RNAP_component_cistron_ids = [
		monomer_id_to_cistron_id_dict.get(monomer_id) for monomer_id in RNAP_component_monomer_ids]
	RNAP_component_RNA_indices = np.array([], dtype=int)
	for cistron_id in RNAP_component_cistron_ids:
		RNAP_component_RNA_indices = np.append(
			RNAP_component_RNA_indices,
			sim_data.process.transcription.cistron_id_to_rna_indexes(cistron_id))
	RNAP_component_RNA_indices = np.unique(RNAP_component_RNA_indices)
	RNAP_component_RNA_ids = np.array([
		RNA_idx_to_id_dict[RNA_index] for RNA_index in RNAP_component_RNA_indices])
	RNAP_component_monomer_indices = np.array([
		monomer_idx_dict.get(monomer_id) for monomer_id in RNAP_component_monomer_ids])

	assert len(RNAP_component_RNA_ids) == len(RNAP_component_RNA_indices)
	assert len(RNAP_component_monomer_ids) == len(RNAP_component_monomer_indices)

	return RNAP_component_RNA_ids, RNAP_component_RNA_indices, RNAP_component_monomer_ids, RNAP_component_monomer_indices


def determine_ribosome_component_ids_and_indices(sim_data):
	"""
	Determines the ids and indices of ribosome components using the ribosome
	flag in sim_data.

	Returns:
		ribosome_component_ids: ids corresponding to ribosome components
		ribosome_component_indices: indexes in monomer_data table corresponding
			to ribosome components
	"""

	monomer_sim_data = sim_data.process.translation.monomer_data.struct_array
	monomer_id_to_cistron_id_dict = dict(
		zip(monomer_sim_data['id'], monomer_sim_data['cistron_id']))
	RNA_data = sim_data.process.transcription.rna_data
	RNA_idx_to_id_dict = {i: RNA[:-3] for i, RNA in enumerate(RNA_data['id'])}
	monomer_idx_dict = {
		monomer: i for i, monomer in enumerate(monomer_sim_data['id'])}

	ribosome_component_monomer_ids = sim_data.molecule_groups.ribosomal_proteins
	ribosome_component_cistron_ids = [
		monomer_id_to_cistron_id_dict.get(monomer_id) for monomer_id in ribosome_component_monomer_ids]
	ribosome_component_RNA_indices = np.array([], dtype=int)
	for cistron_id in ribosome_component_cistron_ids:
		ribosome_component_RNA_indices = np.append(
			ribosome_component_RNA_indices,
			sim_data.process.transcription.cistron_id_to_rna_indexes(cistron_id))
	ribosome_component_RNA_ids = np.array([
		RNA_idx_to_id_dict[RNA_index] for RNA_index in ribosome_component_RNA_indices])
	ribosome_component_monomer_indices = np.array([
		monomer_idx_dict.get(monomer_id) for monomer_id in ribosome_component_monomer_ids])

	gene_sim_data = sim_data.process.replication.gene_data
	gene_id_to_cistron_id_dict = dict(
		zip(gene_sim_data['name'], gene_sim_data['cistron_id']))

	rRNA_operons = sim_data.molecule_groups.rRNA_operons  # this is a list of operon common names
	rRNA_operons_genes = []
	for rRNA_operon in rRNA_operons:
		rRNA_operons_genes.extend(sim_data.molecule_groups.__dict__[rRNA_operon])
	rRNA_operons_cistron_ids = np.array([gene_id_to_cistron_id_dict[gene_id] for gene_id in rRNA_operons_genes])
	rRNA_operons_RNA_indices = np.array([], dtype=int)
	for cistron_id in rRNA_operons_cistron_ids:
		rRNA_operons_RNA_indices = np.append(
			rRNA_operons_RNA_indices,
			sim_data.process.transcription.cistron_id_to_rna_indexes(cistron_id))
	rRNA_operons_RNA_indices = np.unique(rRNA_operons_RNA_indices)
	rRNA_operons_RNA_ids = np.array([
		RNA_idx_to_id_dict[RNA_index] for RNA_index in rRNA_operons_RNA_indices])

	ribosome_component_RNA_indices = np.append(ribosome_component_RNA_indices, rRNA_operons_RNA_indices)
	ribosome_component_RNA_ids = np.append(ribosome_component_RNA_ids, rRNA_operons_RNA_ids)

	assert len(ribosome_component_RNA_ids) == len(ribosome_component_RNA_indices)
	assert len(ribosome_component_monomer_ids) == len(ribosome_component_monomer_indices)

	return ribosome_component_RNA_ids, ribosome_component_RNA_indices, ribosome_component_monomer_ids, ribosome_component_monomer_indices


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
		expression_list_index = 0
		expression_factor = NEW_GENE_EXPRESSION_FACTORS[expression_list_index]
		# Note: Translation efficiency should not matter since gene is knocked out,
		# so we will set it to the smallest translation efficiency value, which
		# comes last in NEW_GENE_TRANSLATION_EFFICIENCY_VALUES.
		trl_eff_list_index = -1
		trl_eff_value = NEW_GENE_TRANSLATION_EFFICIENCY_VALUES[trl_eff_list_index]
	else:
		trl_eff_list_index = index % SEPARATOR
		if trl_eff_list_index == 0:
			expression_list_index = index // SEPARATOR
		else:
			expression_list_index = index // SEPARATOR + 1

		expression_factor = (
			10**(NEW_GENE_EXPRESSION_FACTORS[expression_list_index] - 1))
		trl_eff_value = NEW_GENE_TRANSLATION_EFFICIENCY_VALUES[trl_eff_list_index]

	return expression_list_index, trl_eff_list_index, expression_factor, trl_eff_value


def induce_new_genes(sim_data, index):
	"""
	'Induces' new genes by setting their expression levels and translation
	effiencies to the values specified using the variant index.
	"""
	# Map variant index to expression factor and translation efficiency value
	(expression_list_index, trl_eff_list_index, expression_factor,
	trl_eff_value) = get_new_gene_expression_factor_and_translation_efficiency(
		sim_data, index)

	# Determine ids and indices of new genes
	new_gene_mRNA_ids, new_gene_indices, new_gene_monomer_ids, \
		new_gene_monomer_indices = determine_new_gene_ids_and_indices(
		sim_data)

	print("New gene indices: ", new_gene_indices)

	# Determine ids and indices of RNAP and ribosome components
	RNAP_component_RNA_ids, RNAP_component_RNA_indices, \
		RNAP_component_monomer_ids, RNAP_component_monomer_indices = determine_RNAP_component_ids_and_indices(sim_data)
	ribosome_component_RNA_ids, ribosome_component_RNA_indices, \
		ribosome_component_monomer_ids, ribosome_component_monomer_indices = determine_ribosome_component_ids_and_indices(sim_data)

	indices_to_not_adjust = np.append(RNAP_component_RNA_indices, ribosome_component_RNA_indices)
	indices_to_not_adjust = np.unique(indices_to_not_adjust)
	indices_to_not_adjust = indices_to_not_adjust.tolist()

	# Modify expression and translation efficiency for new genes
	for i in range(len(new_gene_indices)):
		gene_index = new_gene_indices[i]
		monomer_index = new_gene_monomer_indices[i]

		if (RNAP_RENORMALIZATION and RIBOSOME_RENORMALIZATION and NEW_GENE_RENORMALIZATION) or (index ==0):
			sim_data.adjust_new_gene_final_expression([gene_index], [expression_factor])
			sim_data.process.translation.translation_efficiencies_by_monomer[
				monomer_index] = trl_eff_value
		elif not RNAP_RENORMALIZATION and not RIBOSOME_RENORMALIZATION:
			sim_data.adjust_new_gene_final_expression_modified_renormalization(
				[gene_index], [expression_factor], indices_to_not_adjust, NEW_GENE_RENORMALIZATION)
			sim_data.process.translation.translation_efficiencies_by_monomer[ # TODO: change this to prevent renormalization of new gene, ribosome, and RNAP components
				monomer_index] = trl_eff_value
		else:
			raise Exception("RNAP_RENORMALIZATION, RIBOSOME_RENORMALIZATION, and NEW_GENE_RENORMALIZATION must be set to all True or all False")


def knockout_induced_new_gene_expression(sim_data, index):
	"""
	Intended to be used after induce_new_gene_expression has been used on
	previous generations. Knocks out new genes by setting their expression
	levels to 0. Maintains translation effiencies at the values specified using
	the variant index so that way any pre-existing new gene mRNAs have the same
	translation efficiency as when they were created.
	"""
	# Map variant index to expression factor and tranlsation efficiency value
	(expression_list_index, trl_eff_list_index, expression_factor,
	trl_eff_value) = get_new_gene_expression_factor_and_translation_efficiency(
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
	(expression_list_index, trl_eff_list_index, expression_factor,
	trl_eff_value) = get_new_gene_expression_factor_and_translation_efficiency(
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
			+ " Renormalization of RNAP: " + str(RNAP_RENORMALIZATION)
			+ " Renormalization of ribosome: " + str(RIBOSOME_RENORMALIZATION)
			+ " Renormalization of new genes: " + str(NEW_GENE_RENORMALIZATION),
		), sim_data
