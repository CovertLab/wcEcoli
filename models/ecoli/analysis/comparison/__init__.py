# Active analysis modules to run, in this order.
# Tip: Edit this during development to run ones you're working on first.
ACTIVE = [
	"active_ribosome_counts_histogram.py",
	"coexpression_probabilities.py",
	"doubling_time_histogram.py",
	"excess_protein_monomers.py",
	"gene_position_vs_expression_change.py",
	"growth_histograms.py",
	"massFractionComparison.py",
	"mRNA_copy_numbers.py",
	"mRNA_copy_numbers_growth_genes.py",
	"mRNA_copy_numbers_short_genes.py",
	"mRNA_counts_histogram.py",
	"mRNA_length_histogram.py",
	"mRNA_mass_histogram.py",
	"polycistronic_transcription.py",
	"polycistronic_transcription_extended.py",
	"ppgpp_concentrations.py",
	"ppgpp_expression_growth_genes.py",
	"protein_mass_histogram.py",
	"protein_stoichiometry.py",
	"proteomics_fluxomics_comparison.py",
	"proteomics_fluxomics_validation.py",
	"rrna_to_ribosome_yield.py",
	"tRNA_cistron_expression.py",
]

TAGS = {
	'ACTIVE': ACTIVE,   # all active analyses in this category
	'CORE': [           # the default list to run in development
		"growth_histograms.py",
		"massFractionComparison.py",
		],
	}
