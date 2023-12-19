import importlib


variants = [
	'aa_synthesis_ko',
	'aa_synthesis_ko_shift',
	'aa_synthesis_sensitivity',
	'aa_uptake_sensitivity',
	'add_one_aa',
	'add_one_aa_shift',
	'condition',
	'gene_knockout',
	'mene_params',
	'metabolism_kinetic_objective_weight',
	'metabolism_secretion_penalty',
	'new_gene_expression',
	'new_gene_expression_and_translation_efficiency',
	'new_gene_expression_and_translation_efficiency_internal_shift',
	'param_sensitivity',
	'ppgpp_conc',
	'ppgpp_limitations',
	'ppgpp_limitations_ribosome',
	'remove_aa_inhibition',
	'remove_aas_shift',
	'remove_one_aa',
	'remove_one_aa_shift',
	'rrna_operon_knockout',
	'rrna_location',
	'rrna_orientation',
	'tf_activity',
	'time_step',
	'timelines',
	'wildtype',
	]

def get_function(variant):
	module = importlib.import_module(f'models.ecoli.sim.variants.{variant}')
	return getattr(module, variant)

nameToFunctionMapping = {v: get_function(v) for v in variants}

# Support the old names for compatibility with existing shell scripts.
nameToFunctionMapping.update({
	'geneKnockout': get_function('gene_knockout'),
	'meneParams': get_function('mene_params'),
	'nutrientTimeSeries': get_function('timelines'),
	'tfActivity': get_function('tf_activity'),
	})
