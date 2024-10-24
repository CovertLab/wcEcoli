# Active analysis modules to run, in this order.
# Tip: Edit this during development to run ones you're working on first.
ACTIVE = [
	'aa_synthesis_enzymes.py',
	'aa_synthesis_pathways.py',
	'all_operon_tu_structures.py',
	'amino_acid_uptake_rates.py',
	'corrected_rnaseq_read_counts.py',
	'expression_probabilities.py',
	'fit_cistron_degradation_ls_residuals.py',
	'fit_cistron_expression_ls_residuals.py',
	'fold_changes.py',
	'growth_expression_comparison.py',
	'interpolation.py',
	'metabolite_concentrations.py',
	'mRNA_cistron_expression.py',
	'mRNA_transcript_table.py',
	'ppgpp_expression.py',
	'rRNA_operon_structures.py',
	'start_codon_distribution.py',
	'target_vs_actual_ppgpp_fold_change.py',
	'tf_target.py',
	'tRNA_cistron_expression.py',
]

TAGS = {
	'ACTIVE': ACTIVE,   # all active analyses in this category
	'CORE': [           # the default list to run in development
		'metabolite_concentrations.py',
		'interpolation.py',
		],
	'METABOLISM': [
		'aa_synthesis_enzymes.py',
		'aa_synthesis_pathways.py',
		'amino_acid_uptake_rates.py',
		'metabolite_concentrations.py',
		],
	'REGULATION': [
		'expression_probabilities.py',
		'fold_changes.py',
		'growth_expression_comparison.py',
		'ppgpp_expression.py',
		'tf_target.py',
		],
	'VALIDATION': [
		'aa_synthesis_enzymes.py',
		'amino_acid_uptake_rates.py',
		'fit_cistron_expression_ls_residuals.py',
		'growth_expression_comparison.py',
		],
	}
