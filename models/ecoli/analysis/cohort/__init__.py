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
	# "proteinFoldChangeVsTranscriptionFrequency.py",
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
	"subgenerational_expression_table.py",
	# "transcriptFrequency.py",
	# "transcriptionGenomeCoverage.py",
	# "transcriptionGenomeCoverageSecondHalf.py",
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
		# "transcriptFrequency.py",
		# "transcriptionGenomeCoverage.py",
		# "transcriptionGenomeCoverageSecondHalf.py",
		],
	}
