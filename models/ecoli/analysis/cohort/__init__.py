# Active analysis modules to run, in this order.
# Tip: Edit this during development to run ones you're working on first.
ACTIVE = [
	"aa_conc.py",
	"centralCarbonMetabolismCorrelationTimeCourse.py",
	"centralCarbonMetabolismScatter.py",
	"doubling_times_histogram_all.py",
	"ecocyc_table.py",
	"expression_dynamics.py",
	"growthDynamics.py",
	"growth_time_series.py",
	"histogramDoublingTime.py",
	"histogramFinalMass.py",
	"histogramGrowthRate.py",
	"initialVsFinalMass.py",
	"kinetics_flux_comparison.py",
	"mass_fraction_instantaneous_growth_rates.py",
	"promoter_probabilities.py",
	"proteinCopyNumberDistribution.py",
	"replisome_rnap_collision_per_gene.py",
	"rnaCopyNumberDistribution.py",
	"subgen_extract.py",
	"subgen_seed_ci.py",
	"subgen_expression_table.py",
	"subgen_gen_expression.py",
	"subgen_definitions.py",
	"subgen_peak_counts.py",
	"subgen_protein_distribution.py",
	"subgen_transcription_prob.py",
	"subgen_monomer_dynamics.py",
	"subgen_subsample_mrna_protein_10k.py",
	"subgen_cell_health_metrics.py",
	]

TAGS = {
	'ACTIVE': ACTIVE,   # all active analyses in this category
	'CORE': [           # the default list to run in development
		"proteinCopyNumberDistribution.py",  # TODO(jerry): an empty CORE list could be annoying, so include this?
		],
	'DIVISION': [
		"initialVsFinalMass.py",
		],
	'ECOCYC': [
		"ecocyc_table.py",
		],
	'GROWTH': [
		"aa_conc.py",
		"growth_time_series.py",
		],
	'HETEROGENEITY': [
		"proteinCopyNumberDistribution.py",
		"rnaCopyNumberDistribution.py",
		],
	'METABOLISM': [
		"aa_conc.py",
		"centralCarbonMetabolismCorrelationTimeCourse.py",
		"centralCarbonMetabolismScatter.py",
		"kinetics_flux_comparison.py",
		],
	
	'SUBGEN': [   
		"subgen_extract.py",
		"subgen_seed_ci.py",
		"subgen_expression_table.py",
		"subgen_gen_expression.py",
		"subgen_definitions.py",
		"subgen_cell_health_metrics.py",
		"subgen_peak_counts.py",
		"subgen_protein_distribution.py",
		"subgen_transcription_prob.py",
		"subgen_monomer_dynamics.py",
		"subgen_subsample_mrna_protein_10k.py",
		],
	'PAPER': [
		"centralCarbonMetabolismScatter.py",
		"doubling_times_histogram_all.py",
		"expression_dynamics.py",
		"kinetics_flux_comparison.py",
		"histogramDoublingTime.py",
		"histogramFinalMass.py",
		"histogramGrowthRate.py",
		"mass_fraction_instantaneous_growth_rates.py",
		],
	'TRANSCRIPTION': [
		# "transcriptionGenomeCoverage.py",
		# "transcriptionGenomeCoverageSecondHalf.py",
        # "transcriptFrequency.py",
		],
	}
